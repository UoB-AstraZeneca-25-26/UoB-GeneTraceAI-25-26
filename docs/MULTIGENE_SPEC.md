# MULTIGENE_SPEC — the Kleene multi-gene set-membership layer

**Status: built and shipping.** Implementation `src/pipeline/multi_gene_kleene.py`.
Measurements `src/pipeline/measure_multigene_coverage.py` →
`src/pipeline/outputs/multigene_coverage_results.json`. Coverage detail in
[MULTIGENE_COVERAGE.md](MULTIGENE_COVERAGE.md).

```bash
python src/pipeline/multi_gene_kleene.py --expressed EGFR --not-altered KRAS
```

This layer makes **no prediction**. It reports set membership over measured
facts, with a third truth value for "not measured". There is nothing in it a
later audit can retract.

---

## A1. Why three-valued logic, not boolean

`NOT altered(B)` carries two different meanings that boolean logic collapses:

| situation | should mean | boolean gives |
|---|---|---|
| B sequenced, no driver alteration | real exclusion → **passes** | passes |
| B never sequenced | nothing known → **must not pass** | **passes** |

Under `AND NOT` the second case is worse than in scoring: an unmeasured line
**silently enters the confident match set**. This is the absence-as-negative
fault, already found three times in this project — `gate_audit/03` (never-screened
lines counted as confirmed non-dependencies, which manufactured 98% of a
headline), `gate_audit/07A`, and the coverage handling in the similarity query
path. This layer exists so it does not happen a fourth time, in the UI.

**Kleene strong three-valued logic**, TRUE / FALSE / UNKNOWN:

```
AND  | T  F  U          NOT
-----|---------         ---------
 T   | T  F  U          T -> F
 F   | F  F  F          F -> T
 U   | U  F  U          U -> U
```

`TRUE AND UNKNOWN = UNKNOWN` is the line that matters: an unresolved line can
never reach Matches. `FALSE` dominates, because one hard exclusion is enough
regardless of what else is unknown.

Implemented as `k_and`, `k_not`, `k_and_all`.

---

## A2. What resolves each term

### `expressed(G)` on line L

| state | condition |
|---|---|
| TRUE | measured, value **> EXPRESSED_MIN + BAND** |
| FALSE | measured, value **< EXPRESSED_MIN − BAND** |
| UNKNOWN | L absent from the expression matrix |
| UNKNOWN | G fails the validity guard (`frac_expressed < 0.20`) |
| UNKNOWN | \|value − EXPRESSED_MIN\| ≤ BAND — inside measurement uncertainty |

`EXPRESSED_MIN = 1.0` in log2(TPM+1), the TPM ≥ 1 detection convention
(`docs/TRANSCRIPTOMICS_SPEC.md` E2).

### BAND is measured, not asserted

16 DepMap models carry two independent RNA profiles
(`reference/depmap_profiles.parquet`). Over **64,000 paired gene observations**
from `data/parquet/data_clean/depmap_expr_clean.parquet`:

| quantity | value |
|---|---|
| median \|paired difference\| | 0.2036 log2 |
| SD of paired difference | 0.5042 log2 |
| **implied single-measurement SD** | **0.3565** = 0.5042/√2 |
| **BAND = 1.96 × SD** | **0.699 log2 units** |

So a call is unresolvable when the value lies in **[0.301, 1.699]**.

**Caveat, stated because it matters:** n = 16 models is thin, and both profiles of
a model share a library-prep batch, so 0.3565 is a **lower bound** on true
measurement error. A wider band would move more calls to UNKNOWN, never fewer —
the layer errs toward caution in the direction the uncertainty runs.

### How many calls does the band move?

Measured on 2,000 valid genes × 1,673 lines = 3,346,000 measured calls:

| | value |
|---|---|
| calls inside the band | **449,542 (13.44%)** |
| median **per gene** | **3.44%** |

**The answer is not "very few", so the rule stays.** The mean/median gap is the
finding: most genes have almost nothing near the threshold (median 3.4%), but a
minority sit astride it and contribute the bulk of the 13.4%. Dropping the band
would silently harden ~450,000 coin-flips into confident TRUE/FALSE calls, and
they would concentrate in exactly the genes where the call is least safe.

### `altered(G)` on line L

| state | condition |
|---|---|
| TRUE | L sequenced **and** a driver alteration recorded for G |
| FALSE | L sequenced **and** no driver alteration recorded, **and** the MODERATE rule below does not fire |
| UNKNOWN | L never sequenced |
| **UNKNOWN** | **L sequenced, no HIGH-impact row, ≥1 MODERATE-impact row, and G is on the uncertainty list** |

#### The MODERATE rule — why the fourth row exists

`any_driver` is effectively `max_vep_rank == 3`: **100%** of HIGH-impact rows
carry it against **0.5%** of MODERATE rows. VEP scores an in-frame deletion as
MODERATE, so activating in-frame indels returned a confident FALSE — EGFR
exon-19 deletion in PC9 and HCC827 both did.

That is worse than an UNKNOWN. `NOT altered(EGFR)` would have admitted PC9 to
**Matches**, and a false Match is the one thing this layer exists to prevent. So
where there is *evidence of something the annotation could not rank*, the layer
abstains.

**The rule fires only when all four hold:** line sequenced · no HIGH-impact row ·
≥1 MODERATE-impact row · gene on the uncertainty list. A sequenced line with no
rows at all still returns FALSE.

**The uncertainty list is pre-committed**, generated by
`build_altered_uncertainty_list.py` *before* any call was inspected:

> COSMIC CGC **Tier 1** AND role contains **oncogene** AND `MUTATION_TYPES`
> contains **"O"** → **22 genes**

Source is `Cosmic_CancerGeneCensus_v104_GRCh37.tsv` (768 genes, Tier 1 = 592).
*Note: the GRCh38 copy in this repo holds only 99 genes and is missing EGFR,
TP53, KRAS, PTEN and CTNNB1 — a fragment, not used.*

The list: ATP1A1, BRAF, CBL, CREBBP, CTNNB1, CUX1, EGFR, ERBB2, FLT3, IL6ST,
IL7R, JAK2, KIT, KMT2A, MAP3K1, NOTCH1, PDGFRA, RHOA, STAT3, STAT5B, TBX3, UBTF.

It contains the genes where in-frame indels are a known activating mechanism —
EGFR exon 19, KIT exon 11, FLT3-ITD, NOTCH1 PEST — and **excludes KRAS**
(`MUTATION_TYPES = "Mis"` only), which is correct: a MODERATE KRAS row is not a
missed activating indel.

**Two honest limits on the instrument:**

1. **Gene-level, not variant-level.** It flags genes where MODERATE calls are
   untrustworthy. It does not identify activating variants and makes no claim
   about any specific mutation.
2. **"O" = in-frame indel is an interpretation.** CGC codes A/D/F/Mis/N/O/S/T
   with no dedicated in-frame-indel class, so such events fall under "Other".
   That reading is not a documented CGC statement.

#### What it costs — measured before adoption

`measure_fix1_abstention_cost.py`, 100 sampled two-term queries per arm:

| arm | queries affected | Matches lost, median | max |
|---|---|---|---|
| random modifier B | **0%** | 0.00% | 0.00% |
| B drawn from the 22-gene list | 100% | **3.05%** | 8.33% |

The rule can touch a median of 46 lines per uncertainty gene (3.2% of the 1,424
sequenced). On ground truth exactly two rows moved — PC9 and HCC827, FALSE →
UNKNOWN — and the 27 HIGH-impact rows that passed still pass.

**Adopted.** ~3% of Matches in the worst case, nothing at all in the common case,
in exchange for the guarantee that a Match is real.

#### Still not covered — declared, not fixed

| | |
|---|---|
| **covered** | driver-flagged point mutations and small indels; in-frame gene fusions; MODERATE-impact rows on the 22 uncertainty genes (as UNKNOWN) |
| **not covered** | copy-number amplification; homozygous deletion; enhancer-juxtaposing translocations that make no protein fusion (IGH-MYC in Raji); exon-level structural deletions absent from the SNV table (CTNNB1 exon 3 in HepG2) |

`altered() == FALSE` means **"no SNV, indel or in-frame fusion found"**, not
"not altered". This legend is printed in the query output whenever an `altered()`
term is used, because the two readings differ and three of twelve ground-truth
failures live in the gap. Wiring the copy-number layer would make amplification
and exon-level deletion visible; not done here.

Sources: `cleaned_track_data/mutations_collapsed.parquet` (`any_driver`,
`max_vep_rank`) and `cleaned_track_data/fusions_gene_level.parquet`
(`any_in_frame`). "Sequenced" = the line appears anywhere in the mutation table.

**Boolean only.** `mutations_collapsed` also carries `variant_count` and
`max_pathogenicity`, which aggregate by max and therefore drift upward with
variant count. A continuous alteration score would import that bias; this uses a
flag and nothing else. (Task B trap table, "Alteration aggregation bias".)

---

## A3. Output buckets — three, never two

| Bucket | Meaning |
|---|---|
| **Matches** | every clause resolved on measured data, conjunction TRUE |
| **Possible** | at least one clause UNKNOWN — the offending gene and reason are named per line |
| **Excluded** | some clause definitively FALSE |

An UNKNOWN line never enters Matches. Worked example
(`expressed(EGFR) AND NOT altered(KRAS)`):

```
MATCHES 792    POSSIBLE 318    EXCLUDED 563
```

**Those 318 lines are the point of the layer.** Under boolean `AND NOT` they
would have been distributed into Matches and Excluded with no signal that
anything was missing.

Per-line blame is reported:

```
ach-001339   expressed(EGFR)     within measurement uncertainty of the threshold (+/-0.699 log2)
ach-001979   NOT altered(KRAS)   line was never sequenced
```

---

## A4. Coverage reporting

Printed at the top of every query, never buried:

```
COVERAGE  1,355 of 1,673 lines resolved (81.0%)
  expressed(EGFR)      1,512 resolved (90.4%)   161 unknown -- within measurement uncertainty
  NOT altered(KRAS)    1,424 resolved (85.1%)   249 unknown -- line was never sequenced
```

Three things per query: the joint figure, the per-term breakdown (so a user sees
**which gene is costing them coverage**), and per-line reasons for the
unresolved.

### Measured joint coverage — and independence was tested, not assumed

Genes that pass the validity guard, line-level resolution, 200 random draws per
arity:

| terms | observed | independence predicts | ratio |
|---|---|---|---|
| 1 | 85.2% | 85.2% | 1.00× |
| 2 | 73.7% | 73.5% | 1.00× |
| 3 | 62.4% | 62.2% | 1.00× |
| 4 | 59.7% | 59.3% | 1.01× |

Mixed `expressed(A) AND NOT altered(B)`: observed **74.2%** against 74.1%
predicted.

**Finding, contrary to the brief's expectation:** for expression terms, coverage
is *effectively independent* — the observed and predicted curves agree to within
1%. The co-variation in which lines were profiled does not measurably help. The
independence approximation happens to be safe here, but it was checked rather
than assumed, and the check is what licenses using it.

**Two different denominators, do not confuse them:**

- **70.1%** (13,047 of 18,623 genes) is the *gene-level* validity rate — how
  often a randomly chosen gene is usable at all.
- **85.2%** is the *line-level* resolution rate for a gene that already passes.

A query on a randomly chosen gene faces both: ~30% chance the term is
UNKNOWN everywhere, and ~15% of lines unresolved if it is not.

---

## A5. Abstention

| condition | trigger | measured rate over 3,000 random genes |
|---|---|---|
| `SILENT_GENE` | `frac_expressed < 0.20`; term is UNKNOWN everywhere and constrains nothing | **29.6%** |
| `NEAR_UNIVERSAL_GENE` | `NOT expressed(G)` where G is expressed in > 95% of lines; excludes almost everything | **43.8%** |
| `TRIVIALLY_SATISFIED` | `NOT altered(G)` where fewer than 10 lines carry a driver alteration | **78.2%** |
| `LOSS_IS_PROTEIN_LEVEL` | `NOT expressed(G)` where G is a CGC tumour suppressor — truncating mutations leave transcript intact, so the term cannot find G-null lines | 324 of 18,623 genes |

`LOSS_IS_PROTEIN_LEVEL` suggests the alternative rather than only flagging the
problem:

```
!! ABSTENTION [LOSS_IS_PROTEIN_LEVEL] on NOT expressed(TP53)
   TP53 loss is typically protein-level. Truncating mutations leave transcript
   intact, so this term will not find TP53-null lines. Consider altered(TP53)
   instead.
```

It is measured, not assumed: the `deletion_null` ground-truth category scores
**25%** agreement — TP53/SAOS-2 at 5.09, PTEN/U87MG at 4.93, RB1/SAOS-2 at 3.40
are all abundantly transcribed. Genuine transcript-level losses do pass
(CDKN2A/U2OS 0.30, RB1/WERI-Rb-1 0.00), so the boundary is homozygous deletion
versus truncating mutation.

Two-gene `expressed(A) AND expressed(B)`: at least one abstention fires in
**53.2%** of random gene pairs.

The 29.6% reconciles with the 30.5% in `docs/TRANSCRIPTOMICS_SPEC.md` F2 (that
figure is on a 19,173-gene universe, this on 18,623). Full breakdown in
[MULTIGENE_COVERAGE.md](MULTIGENE_COVERAGE.md).

---

## A6. Clustered output, not a lineage filter

The user supplies no lineage. Results are grouped by lineage **after the fact**,
from `cleaned_track_data/sample_info.parquet` (30-value controlled vocabulary).

- **Ordering: by match count, descending.** Recommended because it puts the
  lineages a user would act on first. The alternative is alphabetical, which is
  stable across queries and easier to scan when comparing two result sets;
  available via `cluster_order="alpha"`.
- Every cluster reports **matches / possible / excluded / n_lines**, so a
  heavily-unresolved lineage is visible as information rather than as absence.

```
lineage                     match   poss   excl  lines
lung                          118     31     64    213
central_nervous_system         74      4      8     86
upper_aerodigestive            55      2      1     58
```

A cluster with 40 matches (`lung`, 118) reads as a place to look. A cluster with
one match reads as a single candidate, and is shown with the same columns so the
difference is in the numbers rather than in the presentation. A cluster like
`unknown` (29 matches, **156 possible**, 37.3% resolved) is a warning that its
lines are mostly unsequenced, not that they failed.

**This clustering is descriptive.** It is a grouping by a metadata label, not a
claim about biology.

---

## A7. Ordering within a cluster

Each gene's value appears in **its own column**. Nothing is combined.

```
TOP MATCHES (each term in its own column; nothing combined)
             expressed(EGFR)   NOT altered(KRAS)          lineage
ach-001113             6.180               0.000             lung
ach-001289             1.816               0.000      soft_tissue
```

The user sorts by whichever term they care about, and the ordering is
transparently "sorted by a number you can see" rather than by a judgement.
Nothing combined means nothing to validate and nothing to retract.

`query()` returns `combined_score: None` explicitly, with a note pointing at
[COMBINED_SCORE_PRECONDITION.md](COMBINED_SCORE_PRECONDITION.md). **No combined
score reaches the interface until Task B Stage 2 clears.**

---

## A8. Alternatives per line

Each matched line carries the validated similarity claim from
`src/pipeline/cell_similarity_bridge.py`: **+0.186 [+0.149, +0.227]** drug-response
concordance over a same-tissue baseline (GDSC2, 649 lines × 150 compounds,
cluster bootstrap over cell lines).

**The synergy, stated:** the similarity graph is strongest **within lineage**
(**+0.2519** against +0.1856 overall — `cell_similarity/08_drug_arm.py` §8C). The
clustered output of A6 is therefore exactly the setting where the validated
component performs best, and the within-lineage figure is reported alongside
clustered results rather than only the pooled one.

---

## A9. Performance — confirmed, not assumed

Expression is 1,673 × 18,623. Measured on this machine, `benchmark()`:

| | before | after |
|---|---|---|
| cold load | 13,835 ms | **3,144 ms** |
| warm query (median of 5) | — | **27–107 ms** |

The boolean work was never the cost — the 507 MB source CSV parse was. The matrix
is cached to `cell_similarity/outputs/expression_24q4.parquet` on first use.
Real-time confirmed at ~0.03–0.1 s per query once loaded.

---

## What this layer does not do

- No combined score, no ranking by a derived number (A7).
- No dependency prediction. That claim died in `gate_audit/` and does not
  reappear here.
- No lineage filtering — grouping only (A6).
- No claim that a Match will respond to anything. Membership in a measured set is
  all that is asserted.
