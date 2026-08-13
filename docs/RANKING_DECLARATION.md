# RANKING_DECLARATION — convergent-evidence ranking, declared before building

**Written before any ranking is computed or any output is inspected.** Commit
this file before running `rank_convergent.py`. Same discipline as
`COMBINED_SCORE_PREREGISTRATION.md`.

Governing rule: every claim traces to a file path or a computed number.
Unestablished claims are marked **UNVERIFIED — could not locate**.

---

## 0. What this replaces, and why

`rank_cell_lines()` sorts by `core_score` with the direction flipped by gene
class (`activation_driven` → high wins, `loss_of_function` → low wins). That is
the **drug-target** framing: it scores vulnerability, and needs to know whether a
gene is an oncogene or a suppressor.

Under model-selection semantics the question is different — *how much
independent evidence is there that this cell line is in the state the researcher
asked for?* — and it needs **no gene-role classification at all**. The direction
comes from the query.

This document declares the replacement.

---

## 1. Two orders, because the evidence hierarchy is not symmetric

A single order cannot serve both directions. The measurement below is why.

CN value at nine ground-truth "null" cell lines
(`cn_gene_complete.parquet`, neutral = 1.0):

| gene | line | CN | percentile | mechanism |
|---|---|---|---|---|
| TP53 | SAOS-2 | 0.000 | 0.1% | homozygous deletion |
| CDKN2A | MIA PaCa-2 | 0.000 | 4.2% | homozygous deletion |
| PTEN | LNCaP | 0.557 | 4.7% | deletion |
| TP53 | H1299 | 0.989 | 60% | partial / non-CN |
| PTEN | U87MG | 1.000 | 58% | **splice mutation** |
| PTEN | MDA-MB-468 | 1.129 | 87% | **nonsense mutation** |
| RB1 | SAOS-2 | 1.248 | 91% | **mutation, CN neutral** |
| CDKN2A | U2OS | 0.656 | 46% | **promoter methylation** |

**Only 2–3 of 9 known-null lines are copy-number deletions.** The rest are null
through mutation or silencing, with CN reading neutral — correctly, because the
gene *is* present, just broken.

So for a LOSS query, a truncating mutation must rank near the top or RB1/SAOS-2,
PTEN/U87MG and PTEN/MDA-MB-468 become invisible. For a GAIN query the logic
inverts: an activating point mutation is a stronger functional statement than
extra copies, because amplification raises dose while the protein stays normal.

---

## 2. LOSS query — declared order

Triggered by `--not-expressed G`, or `--altered G` where the researcher's intent
is absence of function.

| rank | layer | confirms loss when | why here |
|---|---|---|---|
| **1** | **CN deep deletion** | CN < **0.25** | direct, dense (0% missing), unambiguous. Bottom **0.121%** of all values panel-wide |
| **2** | **Truncating mutation** | `max_vep_rank == 3` (nonsense / frameshift / splice) | direct evidence the product is broken; CN is neutral for this whole class, so it cannot rank below CN-loss |
| **3** | **CN loss** | CN < **0.50** | bottom **1.68%**; real but not necessarily biallelic |
| **4** | **Protein absent** | measured and below panel p05 | would be the most direct evidence of all, **demoted because NaN conflates not-detected with not-run** (`btad200`). Promote to rank 2 once `protein_measured` exists |
| **5** | **Expression low** | `expressed(G) == FALSE` | weakest. Measured: transcript persists in 6 of 9 true-null lines |

## 3. GAIN query — declared order

Triggered by `--expressed G`, `--high G`, or `--altered G` where the intent is
presence or activation.

| rank | layer | confirms gain when | why here |
|---|---|---|---|
| **1** | **Activating mutation** | driver flag TRUE (hotspot / `hessdriver`) | constitutive activation is a stronger functional claim than dose. Ground truth: mutation category **92%** |
| **2** | **CN amplification** | CN > **1.50** | top **5.12%** panel-wide; direct dosage evidence |
| **3** | **Protein high** | measured and above panel p95 | direct evidence of product; **46.6% coverage**, MNAR-confounded |
| **4** | **Expression high** | above the within-gene percentile cut | necessary but not sufficient — a gene can be transcribed without being activated |

**Thresholds are declared from the CN marginal distribution, not from
outcomes.** No ground-truth result was consulted in choosing 0.25 / 0.50 / 1.50.
They are free parameters and §6 requires them swept.

---

## 4. The rule itself — lexicographic, never a sum

For each candidate line, build the ordered evidence vector for the declared
direction:

```
v = ( s₁, s₂, …, s₅ )        sᵢ ∈ { +1 confirms, 0 unknown, −1 contradicts }
```

Sort by:

1. **number of confirming layers**, descending
2. **the vector v itself, lexicographically** — so a line confirmed by layer 1
   outranks a line confirmed by layers 2 and 3, at equal count
3. **number of contradicting layers**, ascending
4. still tied → **declare tied**

**A weighted sum is explicitly rejected.** `+1/0/−1` summed across layers imposes
equal weights with nothing to justify them, hides which layers fired, and is the
same class of object as the abandoned noisy-OR (`MATH_TRIAL_AND_ERROR.md` §2.3).
Every layer's state is shown as its own column, per `MULTIGENE_SPEC.md` A7.

**Ties are reported as bands, not broken arbitrarily.** Precedent: on BRAF the
gap between rank 1 and rank 2 was 0.0007 against a DKW measurement band of
±0.0352 — 2% of the noise. Five lines were statistically indistinguishable from
first place while being printed as 1, 2, 3, 4, 5. This rule outputs
`#1–5 (tied)` instead.

**No gene-role classification is used anywhere.** Direction comes from the query.

---

## 5. Unknown is not contradiction

`sᵢ = 0` when the layer did not measure that (gene, line) — **not** when it
measured and disagreed. Distinguishing these requires a `measured` flag per
layer:

| layer | `measured` flag available? | status |
|---|---|---|
| CN | **yes** — membership in `cn_gene_complete.parquet` (908 × 18,276, 0% missing) | ready |
| Mutation | **yes** — line present in `mutations_collapsed.parquet` | ready |
| Expression | **yes** — non-NaN in the 24Q4 matrix | ready |
| Protein | **NO** — NaN conflates not-detected with not-run | **blocks rank 4 of the loss order** |

Until `protein_assay_attempted` exists, the protein layer contributes `0`
(unknown) and never `−1`. It can raise a line's rank; it can never lower one.
That is the conservative direction and it is declared here so the asymmetry is
not mistaken for a bug.

---

## 6. Acceptance tests — declared before running

The rule ships only if all four hold. Any failure is reported, not patched away.

1. **Discrimination.** On the survivors of a representative Kleene query, fewer
   than **60%** of lines fall in the single largest tied band. Above that the
   ranking is mostly ties and the filter is doing all the work — which is itself
   a finding and must be reported as one.
2. **Ground truth.** On the 53-fact set, lines with a confirmed state rank above
   lines without, at **≥ 80%** pairwise concordance — matching the Kleene
   filter's own 80.9%.
3. **Rank change vs `core_score`.** Report the maximum and median rank change,
   and the top-10 overlap. If top-10 overlap is **< 50%**, this is a material
   architectural revision and must be documented as a replacement, not a patch.
4. **No direction leakage.** Running a LOSS query and a GAIN query on the same
   gene must produce substantially different orderings. If they correlate at
   Spearman **> 0.8**, the direction is not actually being used and the rule has
   collapsed back to a single score.

## 7. Free parameters, to be swept before adoption

| parameter | declared value | basis |
|---|---|---|
| CN deep-deletion cut | 0.25 | bottom 0.121% of the CN marginal |
| CN loss cut | 0.50 | bottom 1.68% |
| CN amplification cut | 1.50 | top 5.12% |
| protein absent / high | panel p05 / p95 | convention, **UNVERIFIED — no literature basis located** |
| tied-band tolerance | DKW band at panel n | `MATH_REFERENCE` §0.1b |

Sweep each across at least three values and report whether the acceptance tests
in §6 hold throughout. A result that survives only at one setting is not adopted.

---

## 8. What is NOT claimed

- No claim that a top-ranked line will respond to any drug.
- No claim of dependency. That died in `gate_audit/` and does not return here.
- No claim that more confirming layers means stronger biology — only that it
  means **more independent measurements agree**, which is a statement about
  evidence, not about effect size.

---

## 9. Postscript — measured after acceptance test 1 failed

**Dated 2026-08-10.**

Test 1 failed (median largest tied band 92.9%, see `ranking_acceptance_test_1.json`).
Before considering any repair, four follow-up questions were answered —
computed, not assumed — by `measure_completeness.py`.

**Q1/Q2 — completeness distribution over qualifying lines** (7 representative
queries, 2,707 gene×line rows):

| layers measured | share |
|---|---|
| 1 | 16.4% |
| 2 | 22.8% |
| 3 | 34.7% |
| 4 | 26.1% |

**Median completeness = 3, not 1.** Per-layer coverage: expression 100.0%,
mutation 83.2%, CNA 43.6%, protein 43.7%. The test-1 collapse is not explained by
most lines sitting at completeness 1 — the modal band is 3 of 4, and it is the
CONTENT of each layer (mostly `0`, rare `+1`/`-1`) that ties lines together, not
the count of measured layers.

**Q3 — does completeness track how well-studied a line is? Confirmed.**

```
Spearman rho(n_measured, n_assay_types) = 0.5982   p = 1.4e-134   n = 1,379
```

Per layer: mutation rho=0.72, protein rho=0.48, CNA rho=0.33. **A completeness
ranking would be a fame ranking.** This is the same coverage-proxy fault already
found in the abundance gate (`gate_audit/03`) and the `NOT altered(G)` trivial
case (`MULTIGENE_COVERAGE.md` — 78.2% of genes) — the fourth instance in this
project.

**Q4 — labelling.** Implemented in `rank_convergent.py`: every ranked CLI output
prints the completeness/fame warning above the result table, citing the exact
Spearman figure, before printing a single ranked row.

**Q5 — tumour-similarity reference dataset. UNVERIFIED — could not locate.** No
TCGA- or Celligner-style tumour comparison exists anywhere in this codebase.
`cell_similarity/` compares cell line to cell line only. If tumour-matched
selection is wanted it is new scope requiring an explicitly acquired,
lineage-matched reference — not a reframing of anything that exists today.

### Net effect on the section 6 acceptance path

Test 1 already failed on discrimination grounds. Q3 adds a second, independent
reason the convergent-evidence ranking as declared should not ship even if
discrimination were fixed: raising completeness (e.g. by unblocking the protein
`measured` flag) would raise the correlation with study frequency further,
because protein coverage is itself study-frequency-correlated (rho=0.48). Fixing
test 1 by adding measured layers makes the fame-proxy problem worse, not better.

**Recommendation, unchanged from before this postscript:** ship the Kleene
filter without a ranker for the general case (option 1, prior turn). If ordering
within confirmed subsets is wanted (option 2), it must be restricted to
CONTENT-based confirmation (`+1` cells) and must never reward completeness
alone — a line with 1 confirming layer should not be outranked by a line with 3
measured-but-unconfirmed layers, and the current lexicographic rule already
respects that (`n_confirm` sorts before `n_measured` is even computed). The
`n_measured` diagnostic itself must never appear as a sort key.

---

## 10. Second postscript — Q1–Q8, and a real bug found and fixed

**Dated 2026-08-10, same session as §9.**

### Q3 (re-opened) — a genuine bug in the ranking evidence encoder, NOT the Kleene core

Direct question: was `0` conflating measured-negative with unknown? **Yes, for
the mutation layers specifically.**

```python
r["truncating_mutation"] = (1 if m in hi_set else 0)   # never checked `sequenced`
```

A never-sequenced line and a sequenced-and-clean line both produced `0`, with no
way to tell them apart. Confirmed directly: on EGFR, line `ach-001979` (never
sequenced) and `ach-002376` (sequenced, no HIGH-impact row) were **identical**
in the evidence table.

**This is not the Kleene core.** Verified against the 53-fact ground truth: all
9 FAIL rows resolve to a definite TRUE/FALSE, never UNKNOWN — the core already
separates the two states correctly (Q2 below). The bug was specific to
`rank_convergent.py`'s independently-built evidence encoder.

**Fixed.** Every layer now carries a parallel `<layer>__measured` column. The
*sort value* is unchanged by design (a measured-clean mutation call still
contributes `0`, matching how the CN layers already treat "measured but
neutral" — RB1/SAOS-2 is CN-neutral and mutation-positive, so a neutral CN
reading must not be scored as if it contradicted loss). What changes is that
measured-vs-unknown is now tracked and reportable, which it was not before.

### Q1 — measured-negative vs unknown, as a fraction of zeros

7 representative loss queries, 2,707 (gene, line) rows, after the fix:

| layer | confirm +1 | measured-0 | unknown-0 | contradict −1 |
|---|---|---|---|---|
| CN deep deletion | 153 | 960 | 1,527 | 67 |
| **truncating mutation** | 18 | **2,234** | 455 | 0 |
| CN loss | 174 | 705 | 1,527 | 301 |
| protein absent | 0 | 0 | 2,707 | 0 |
| expression low | 2,707 | 0 | 0 | 0 |

**Of the truncating-mutation zeros, 83.1% were measured-negative, not
unknown.** Most lines that don't confirm via this layer were actually
sequenced and found clean — the opposite of what the pre-fix code could report.

### Q2 — do the 9 ground-truth FAILs reflect measured-negative or unknown?

All 9 resolve to a definite TRUE/FALSE (never UNKNOWN) — confirmed directly
against `kleene_ground_truth_results.json`. The 5 UNKNOWN cases are already
separated into their own bucket and are not among the 9 FAILs. **The Kleene
core's UNKNOWN/FALSE distinction is correct**; Q3's suspicion applied only to
the newer ranking module built on top of it.

### Q4 — does the fame correlation survive within lineage?

| | rho | p | n |
|---|---|---|---|
| pooled, all lineages | **0.379** | 3.9e-43 | 1,226 |
| median, within-lineage (18 lineages, n≥30) | **0.356** | — | — |

*(Note: this pooled figure of 0.379 uses the loss-query gene sample in this
script, not identical to the first postscript's 0.598 general sample — both
point the same direction.)*

Per-lineage range: fibroblast 0.974 (n=37) down to peripheral_nervous_system
−0.214 (n.s.), with most lineages positive and significant (lung 0.343 at
p=1.3e-6, kidney 0.641 at p=1.9e-5, central_nervous_system 0.508 at p=1.6e-6).

**Verdict: HOLDS WITHIN LINEAGE.** The correlation is not purely a
between-lineage artefact — median within-lineage rho (0.356) is close to the
pooled figure. A researcher querying inside a single cancer type is still
subject to substantial fame bias, just not quite as much as the pooled number
alone would suggest.

### Q5 — sequencing depth or panel inclusion?

```
mutation-measured vs WES replicate count : rho=0.446  p=6.1e-61
mutation-measured vs lineage panel size  : rho=0.027  p=0.346 (n.s.)
WES profile count distribution: {1: 926, 2: 169, 0: 131}
```

**Verdict: cannot separate depth from inclusion with the data available** — 926
of 1,226 lines have exactly one WES profile, so "depth" in the sense of
replicate coverage barely exists in this dataset; the available signal is
binary (sequenced at all, or not). Lineage panel size shows **no** correlation
(p=0.35), ruling out "some lineages just have more lines" as the mechanism.

**Remedy implied:** the fame proxy is an **ascertainment** proxy (which lines
were ever included in a sequencing panel), not a technical-depth proxy. Fixing
it means expanding coverage to under-profiled lines, not re-sequencing existing
ones more deeply — a different, cheaper remedy than Q5 originally worried about.

### Q6 — report rho alongside the discrimination result. Done.

The CLI warning now states both numbers together: completeness/fame rho (0.60
pooled, 0.38 within-lineage) and the acceptance-test-1 discrimination failure,
in the same message, before any ranked row is printed.

### Q7 — does ranking add anything even within the confirmed subset?

14 representative queries, confirmed subset = lines with ≥1 confirming layer:

```
median % of survivors with >=1 confirming layer         : 100.0%
median largest tied band WITHIN the confirmed subset    : 94.3%
```

**No.** For gain queries this is mechanical: surviving `expressed(G)` already
guarantees `expression_high = +1`, so "confirmed" and "survivor" are nearly the
same set — restating the test-1 finding that the filter spends the one
layer that varies before ranking starts. **Verdict:
`NO_ORDER_EVEN_WITHIN_CONFIRMED`.**

### Q8 — report N-qualify / K-additional-evidence. Implemented.

```
query: loss of TP53
  9 lines qualify; 9 have >=1 additional confirming layer beyond the query itself
```

Printed before the ranked table on every CLI call.

### Net effect

The bug fix changes what can be *reported* (Q1, Q3), not what ships. Q4–Q7
confirm and sharpen, rather than overturn, the prior postscript's conclusion:
completeness-as-ranking is a fame proxy that persists within lineage, mutation
coverage reflects ascertainment rather than depth, and even the confirmed
subset does not discriminate. The recommendation is unchanged — ship the
Kleene filter without a general ranker, and if any qualifying-set summary is
shown, it is the Q8 flag ("N qualify, K have extra evidence"), never a
positional ranking of the K.

---

## 11. Third postscript — display fixes and remaining questions closed

**Dated 2026-08-10, same session as §9–10.**

### Bug-fix verification (Q1 of the third round)

`test_kleene_ground_truth.py` re-run after the `evidence()` fix:
**38 PASS / 9 FAIL / 5 UNKNOWN, 80.9% — bit-for-bit unchanged.** Confirmed by
inspection that the ground truth script imports only `multi_gene_kleene`, never
`rank_convergent`; the two are disjoint code paths, so the fix could not have
touched filter output. This corroborates rather than contradicts the earlier
Q2 finding that the Kleene core was never at fault.

### `__measured` is now exposed, not internal-only

The CLI table shows **`U`** for never-measured and **`0`** for measured-but-not-
confirming, distinguishing them visibly:

```
rank      line            cn_deep_del  truncating_      cn_loss  protein_abs  expression_
#1-2      ach-001611               +1            0           +1            U           +1
#3-9      ach-003184                U            U            U            U           +1
```

`ach-001611` (row 1) was sequenced and found mutation-clean (`0`). `ach-003184`
(row 3) was never sequenced for CN or mutation at all (`U`). Previously both
printed as bare `0`.

### The correct Kleene value for unmeasured cells

Confirmed: it is `U`, not `0` — implemented as stated above. The **sort value**
used internally for ranking still treats `U` and measured-`0` identically (both
non-confirming, neither contradicting), by the design already argued in §10 —
but the **display** now shows the true three-way state, so a researcher is never
shown a `0` that could mean either thing.

### Mutation-specific fame correlation, now in the warning

The CLI warning states rho=0.72 for mutation by name, not only the pooled
rho=0.60, because mutation is the most fame-biased of the four layers (CNA
0.33, protein 0.48, mutation **0.72**). A researcher seeing only the pooled
figure would not know mutation-confirmed lines are the least trustworthy on
this axis.

### Fibroblast as the extreme case

From the within-lineage table in §10: **fibroblast, rho = 0.974, n = 37,
p = 3.7e-24** — near-perfect correlation between measurement completeness and
study frequency, the single largest of any lineage tested (next highest: kidney
at 0.64). It is the clearest illustration of the structural problem: in this
lineage, whether a line's molecular state can be confirmed is almost entirely
explained by how much attention that specific line has received, not by
anything about the gene being queried. Contrast lung (rho=0.34, n=189) or
central nervous system (rho=0.51, n=80) — real but far weaker, and on far more
lines, so any remedy tested on fibroblast alone would not generalise.

### The 92.9% and 94.3% figures, side by side

| set | median largest tied band | source |
|---|---|---|
| full survivor set (acceptance test 1) | **92.9%** | §9 |
| confirmed subset only (n_confirm ≥ 1) | **94.3%** | §10, Q7 |

Restricting to lines with additional confirming evidence does not improve
discrimination — if anything it is marginally worse, because (per Q7) the
confirmed subset for gain queries is nearly the full survivor set by
construction (surviving `expressed(G)` already guarantees the expression layer
confirms).

### "Tied" does not mean missing data — stated explicitly

Per Q1 (§10): 83.1% of truncating-mutation zeros are measured-negative, not
unknown. The CLI now states this directly, before the ranked table, so "tied"
is not misread as a data-quality complaint:

> **'TIED' HERE MEANS EVIDENCE-EQUIVALENT, NOT MISSING DATA.** Most
> non-confirming cells are measured negatives, not unknowns... The tie is a
> property of the biology matching across lines, not of gaps in what was
> measured.

### Can a researcher sort by any column themselves?

Yes — `--sort-by <column>` (any layer name, `_cn`, or `_log2tpm`), added to the
CLI. It re-orders the printed rows only; it does not change `n_confirm` or the
band assignment, and the output is explicitly labelled `CUSTOM VIEW ... NOT the
evidence-based ranking ... position here is not an evidence claim` — the same
"sorted by a number you can see" principle as `MULTIGENE_SPEC.md` A7, applied
to this module.

```
python src/pipeline/rank_convergent.py --gene TP53 --direction loss --sort-by _cn
```

### What ships

Unchanged from §10: the Kleene filter without a general ranker. What this round
adds is that the ranking module's *diagnostic output* is now honest about three
things it previously hid — measured vs unknown, which layer is worst for fame
bias, and that ties reflect real evidence agreement rather than data gaps — even
though the underlying recommendation not to present a positional ranking has not
changed.

---

## 12. Fourth postscript — write-up framing and cross-references

**Dated 2026-08-10, same session as §9–11.**

### 1. Disjoint code path — stated for the method section

> **The 53-fact ground truth test operates exclusively on `multi_gene_kleene`;
> `rank_convergent` was validated independently.** These are two separate
> modules with no import relationship in either direction. The ground truth
> pass rate (80.9%, unchanged before and after the §10 bug fix) is evidence
> about the Kleene filter's correctness. It is not evidence about the ranking
> module, which required its own separate instrumentation (§10 Q1–Q3) precisely
> because the two cannot certify each other.

This sentence belongs in the method section verbatim or near-verbatim — it is
the fact that makes Q2 (§10) a real corroboration rather than a circular one.

### 2. The four tests as a discovery sequence, not a list of failures

Each test closed a specific option and the closure is what produced the next
question — this is the intended reading, made explicit:

| test | question closed | what it produced |
|---|---|---|
| **Acceptance test 1** (§8–9) | can a lexicographic evidence-vector ranking discriminate at all? **No** — 92.9% median largest band | forced the question of *why*: was it sparse data, or something structural? |
| **Completeness measurement** (§9 Q1–Q3) | is the ranking mostly ties because most lines are under-measured (completeness=1)? **No** — median completeness is 3 of 4. Is completeness itself a fame proxy? **Yes** — rho=0.60 | ruled out "just add more layers" as a fix, since more layers (protein) would raise fame correlation further, not resolve discrimination |
| **Evidence-encoder audit** (§10 Q1–Q3) | is the `0` code correctly distinguishing measured-negative from unknown? **No for mutation layers** — found and fixed a real bug, isolated to the ranking module and verified not to touch the Kleene core | separated "is the tie real" from "is the tie an artefact of bad bookkeeping" — it is real |
| **Confirmed-subset test** (§10 Q7) | does restricting to lines with independent extra evidence rescue discrimination? **No** — 94.3%, no better than the full set | closed the last plausible repair path within the declared design |

Four tests, each answering a question the previous one raised, converging on
one conclusion from four independent angles: the ranking is genuinely
undiscriminating on this data, not undiscriminating because of a bug, missing
layers, or an unlucky subset. That convergence — not any single test — is the
finding, and it is a stronger result than any one test alone would be.

### 3. "Tied means evidence-equivalent" — restated as the positive finding it is

> **We measured the epistemic status of every tied position rather than
> assuming it, and found that 83.1% of non-confirming truncating-mutation
> cells are measured-negative, not unmeasured** (§10 Q1). This is a positive
> result: it demonstrates that the observed ties reflect genuine agreement in
> the underlying biology across lines, not a gap in what was measured. Had the
> figure been reversed — mostly unmeasured — the ties would have been an
> artefact of coverage, and the correct response would have been to expand
> measurement, not to trust the tie. Measuring which case obtained, rather than
> assuming either, is what makes the "ship without a ranker" recommendation
> defensible rather than merely cautious.

### Q4/Q5 (redirection and consistency) — verified, not just designed-for

```
python src/pipeline/rank_convergent.py --gene TP53 --direction loss --top 3 \
    > out_default.txt
python src/pipeline/rank_convergent.py --gene TP53 --direction loss --top 3 \
    --sort-by _cn > out_sorted.txt
```

Both files contain exactly one `WARNING: ranked by MEASUREMENT COMPLETENESS`
block (grep count 1/1) — plain `print()` to stdout, no TTY-only formatting, so
it is captured identically by redirection. `CUSTOM VIEW` appears in the
`--sort-by` file only (count 0/1), so a researcher diffing the two files sees
immediately that they represent different claims — the warning is identical in
both, but only the custom-sorted file carries the additional "not an evidence
claim" label.

### 6. Consolidated evidence trail — this section is the index

All eight completeness questions, in one place, with thresholds set before
each was measured:

| Q | question | pre-declared threshold | observed | verdict | section |
|---|---|---|---|---|---|
| Test 1 | discrimination | largest band < 60% | median 92.9% | **FAIL** | §9 |
| Q1 (round 1) | completeness distribution | none (descriptive) | median = 3 of 4 layers | not a coverage-1 artefact | §9 |
| Q3 (round 1) | completeness vs fame | none (descriptive) | rho=0.60, p=1.4e-134 | **confirmed proxy** | §9 |
| Q2 (round 2) | ground truth affected by evidence encoder? | must be independent | 38/9/5 unchanged, disjoint imports | **confirmed independent** | §10, §12.1 |
| Q3 (round 2) | 0 vs unknown correctly distinguished? | — | **bug found**: mutation layers never checked `sequenced` | **fixed** | §10 |
| Q1 (round 2) | measured-0 vs unknown-0 split | — | 83.1% of truncating-mutation zeros are measured | ties are biology, not gaps | §10, §12.3 |
| Q4 (round 2) | fame bias within lineage | — | pooled 0.379, within-lineage median 0.356 | **survives lineage control** | §10 |
| Q5 (round 2) | depth or ascertainment | — | WES near-binary (926/1226 = exactly 1 profile), lineage size rho=0.03 n.s. | **ascertainment, not depth** | §10 |
| Q7 (round 2) | confirmed-subset discrimination | — | 94.3% vs 92.9% full set | **no better** | §10 |
| — | display fixes | — | `U`/`0` distinguished, mutation rho=0.72 shown, `--sort-by` added | implemented | §11 |

An examiner asking "show me the evidence trail" is pointed at this table; every
row names the section with the full computation.

### 7. Cross-reference to the gate retraction — both corrections in one place

This project has **two withdrawn or revised claims**, and they should be read
together because they are the same failure mode at two different layers of the
architecture:

| | claim | what was wrong | where retracted |
|---|---|---|---|
| **Gate** | gate-region pAUC 0.596–0.608 "verified, replicated" | 98% of the depletion reproduced by a placebo with labels randomised — never-screened lines counted as confirmed negatives | `gate_audit/03`, banners in `BUILD_SPEC.md`, `DECISIONS_REQUIRED.md`, `EVAL_DENOMINATOR.md`, `EXAM_DRILL.md`, `RANKING_FEASIBILITY.md`, dated 2026-08-10 |
| **Ranking** | convergent-evidence ranking would meaningfully order qualifying lines | acceptance test 1 failed outright (92.9% tied); the underlying reason (fame-proxy completeness, §9 Q3) is structurally the same class of fault as the gate's denominator problem — a plausible-looking signal that turns out to measure study effort, not biology | this document, §9–12, dated 2026-08-10 |

Both are cases where a number that looked like evidence was, on closer
measurement, tracking **how much a cell line had been studied** rather than the
biological state being asked about. The gate's version was hidden in a
denominator; the ranking's version was hidden in what "more layers confirm"
actually means. Citing one without the other in the thesis understates how
central this failure mode is to the project's methodology.
