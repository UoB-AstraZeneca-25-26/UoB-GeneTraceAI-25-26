# Query layer — abstention states and joint (multi-gene) selection

Two additions to Stage 7, both built against
*CellLineSelector — Confidence Engine: Design Record (v1)*:

1. the architecture now **says outright when it has nothing to say**, instead of
   returning a confident-looking ranking with a caveat underneath
2. a **joint query** — cell lines that suit *several* genes at once, bounded by the
   weakest link

Code: `src/pipeline/evidence_state.py`, `src/pipeline/multi_gene.py`, wired into
`src/pipeline/cli.py`. Nothing here writes to or alters `core_score` — role
separation (C7) holds.

---

## 1. Why abstention was needed

`core_score` is a **within-gene percentile**. Every gene is rescaled onto the same
0–1 spread, so a gene that barely varies across cell lines produces a top-10 that
reads exactly as confidently as a sharply tissue-restricted one:

| gene | top line's score | top line vs a typical line |
|---|---:|---:|
| ERBB2 | 0.9975 | **76×** |
| GAPDH | 0.9989 | **2.8×** |

Identical numbers, completely different biology. GAPDH is a housekeeping gene —
near-uniformly expressed by design, which is why it is used as a loading control.
Its "top" cell line is not meaningfully better than its median one.

The pipeline already *measured* this (`gene_dispersion.signal_spread`) and the CLI
already printed a warning. That was not enough: a warning sits under a ranked table
and a 0.9989, and the table wins the reader's attention. Per **C1**, the fix is a
distinct **state**, not a lower score on the same axis.

## 2. The states

| state | meaning | ranking shown? |
|---|---|---|
| `NO_EVIDENCE` | no layer measures this gene (or this gene/line pair) | **no** — suppressed |
| `NO_EVIDENCE` (too few) | scored in under 30 lines; a percentile is not stable | **no** — suppressed |
| `UNINFORMATIVE` | measured, but does not discriminate between lines | behind an explicit refusal |
| `CONFLICTING` | abundance ranking runs opposite to measured CRISPR dependency | behind an explicit refusal |
| `RANKED` | the ordering is worth acting on | yes, clean |

Design-record provenance:

- **C1** — ordinal bands *plus* distinct non-ordinal states. The record already uses
  `Conflicting` for source disagreement; `NO_EVIDENCE` and `UNINFORMATIVE` are the
  two "cannot answer" cases, which are categorically different from "the answer is
  Low".
- **§5 — "missing modalities widen uncertainty, never lower the band (missing ≠
  against)."** A gene with no protein layer is not *low confidence* for that reason;
  it is *less corroborated*. Absence is reported as absence. The CLI now says
  "Not a low score — no score."
- **L4** — the record flags that the reliability proxy diverges "in a patterned way
  at ubiquitously-high and flat-profile genes". `UNINFORMATIVE` is exactly that
  population, named and handled rather than left as a known distortion.
- **C3** — evidence is tallied, never averaged. Each state is triggered by a named,
  independently checkable fact.

### Measured behaviour

```
> gene GAPDH
==========================================================================
NO DISCRIMINATING SIGNAL: gapdh does not separate cell lines.
==========================================================================
  The best cell line expresses gapdh only 2.8x above a typical one...
  core_score is a within-gene percentile, so it still spans 0-1 and the top row
  still reads ~0.99. That number reflects the rescaling, not a real gap.
  Any ordering below is dominated by measurement noise. Do not select cell lines on it.
  Median expression is high (log2TPM 12.0) and near-uniform -- the signature of a
  housekeeping gene.

  The table below is shown for completeness only. It is NOT a recommendation.

> gene ERBB2
Gene: erbb2 (ensg00000141736) | class: abundance_tracking | 1632 candidates
      [no banner — the ranking stands on its own]
```

**Thresholds** (`evidence_state.py`) are documented design choices, not literature
values — the same stance §6 of the design record takes.

The flat test **defers to `gene_dispersion`'s own `narrow` band** rather than
applying an absolute fold cutoff. That distinction matters and was found by testing:
BRAF sits at **2.9×** top-vs-median and GAPDH at **2.8×**, but the banding (derived
in Cell 7c against genes that are actually expressed) calls BRAF `moderate` and
GAPDH `narrow`. An absolute ceiling near 3 refuses to rank BRAF — the project's own
anchor gene. `FLAT_FOLD_CEILING = 3.0` is therefore only a **fallback for when the
band is missing**, never a second gate. `MIN_LINES_FOR_RANKING = 30`.

**Activation-driven genes are exempt.** Stage 4 ranks them on the driver flag with
`core_score` only as support, so a flat *abundance* profile does not invalidate a
ranking that never rested on abundance. Those genes return `RANKED` with an
explanatory note instead of a refusal.

---

## 3. Joint (multi-gene) selection

```
genes <A> <B> [C ...]
```

Returns a **short list** of cell lines that suit every gene asked for.

### The combination rule — weakest link, not average

Design record **§5**: *"Conjunctive claims (needs several facts jointly) are bounded
by the weakest necessary link (a min-like rule / Fréchet–Hoeffding lower bound)."*

```
joint(line) = min over genes of core_score(gene, line)
```

A line scoring 0.99 for gene A and 0.10 for gene B gets **0.10**, not 0.55. A mean
would let a strong gene carry a weak one and surface lines that suit neither claim
jointly — precisely the failure a conjunctive query must avoid. The minimum is the
Fréchet–Hoeffding lower bound on P(A ∧ B): the most the conjunction can be
guaranteed when the dependence between the two is unknown.

`core_score` is used rather than `stratum_rank` because it is already a within-gene
percentile in [0,1] and so is comparable across genes; `stratum_rank` is computed
within (gene, `n_layers`) strata and is not.

Each row reports its **limiting gene** — which of the requested genes is holding it
back — so the constraint is visible rather than buried in one number.

### Confidence is also weakest-link

The joint verdict is the **worst** per-gene verdict. If any requested gene is
`NO_EVIDENCE` or `UNINFORMATIVE`, the conjunction inherits it and says so, rather
than quietly dropping that gene from the query:

```
> genes ERBB2 GAPDH
  erbb2      ok                    1,632 lines scored
  gapdh      NOT DISCRIMINATING    1,746 lines scored
==========================================================================
NO DISCRIMINATING SIGNAL: gapdh does not separate cell lines.
==========================================================================
  This is the limiting gene, so the combined answer inherits it.
```

### Redundancy is measured, not corrected (L6)

The record notes that correlated sources must not be double-counted and leaves
`n_eff` as future work. The same applies to genes: asking for two co-expressed genes
is close to asking for one. `multi_gene.py` **measures and reports** the pairwise
score correlation and flags |ρ| > 0.7 as near-redundant. It does not reweight —
reweighting without a validated `n_eff` would be a fabricated correction.

```
> genes ERBB2 MUC1
  Gene-gene score correlation: rho = +0.412   (largely independent constraints)
  1,530 cell lines are scored for every gene; 473 clear 0.50 on all of them.

  model_id      joint    erbb2       muc1        limiting
  ach-000679    0.9820   0.9959      0.9820      muc1
  ach-000568    0.9803   0.9803      0.9939      erbb2
  ...
```

### The floor

A line must clear **0.50 on every gene** to be offered as a joint hit. Without a
floor, `min()` still returns an ordering when every line is mediocre for one gene,
and the caller cannot tell the difference. When nothing clears it, the list is still
shown but explicitly labelled *"NOT a set of good joint matches"*.

---

## 4. What this does not claim

The returned lines are **not** experimentally validated for the combination. This
reports which cell lines the existing evidence ranks highly for every gene asked
for, with the weakest link named. Consistent with **C5** (documented-not-validated),
there is no calibration behind these numbers and none is implied.

---

## 5. Metadata in the output

`build_warehouse_views.py` now builds two annotation views so every consumer
describes an entity the same way instead of re-joining four tables:

- **`v_cell_line`** — one row per line: identity axes (RRID/CVCL, SIDM, COSMIC,
  CCLE, profiles, GEO), DepMap clinical annotation (lineage, disease, subtype,
  sex, age, primary/metastasis, collection site), the independent Sanger/GDSC
  annotation (tissue, cancer type, MSI, ploidy, mutational burden, growth), and
  which of the 11 measured layers reach it.
- **`v_gene`** — one row per gene: HGNC identity, locus type, COSMIC role and CGC
  tier, UniProt accessions, aliases, plus `gene_dispersion` (`signal_spread`,
  `top_vs_median_fold`) because that is what decides whether a ranking means
  anything.

`src/pipeline/metadata.py` reads them and degrades gracefully — with no warehouse
present the CLI falls back to ids rather than failing.

## 6. Web export — `export_web.py`

```bash
python src/pipeline/export_web.py --genes BRAF ERBB2 EGFR --top 25
```

Writes one JSON per gene (~48 KB) plus `index.json` into the React app's
`public/data/`. Static files, no backend: `core_score.parquet` is 29.8M rows and
cannot ship to a browser, but one gene's top-25 with full metadata is a few tens
of KB.

Each payload carries `gene`, `class`, `n_candidates`, the **`verdict`** block
(state / headline / detail / show_ranking / qualify_ranking) and `lines[]` with
metadata and per-layer `tracks`.

**One deliberate omission.** The payload sends `confidence_tier` (ordinal) and
`core_score` (a real within-gene percentile) but **no 0–1 confidence float**. The
design record is explicit in C1 and C5 that the reported unit is an ordinal band
with no calibration behind it. The React prototype currently renders a fabricated
`confidence: 0.96`; reproducing that from real data would invent a number the
pipeline does not have.
