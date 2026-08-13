# Rebuild runbook — ranking-layer fixes

Every change described here is **source-only**. Until this runbook is run, the
parquet artifacts in `src/pipeline/outputs/` are the *old* ones and the CLI still
serves the old numbers. The one exception is the sort fix, which is re-derived at
query time (`explain_pair.sort_gene_rows`) and is therefore already live.

Run from the repo root:

```bash
cd C:\Disertation\UoB-GeneTraceAI-25-26
```

## Step 0 — harmonisation (run first, or the rest rebuilds on stale identity)

Stage 0 is now two notebooks writing one DuckDB warehouse. Architecture and
measured outputs: [HARMONISATION_V2.md](HARMONISATION_V2.md).

| # | notebook | writes |
|---|---|---|
| 0 | `src/pipeline/00_harmonisation.ipynb` | `celllineselector.db` (18 core tables) + `harmonised.parquet`, `gene.parquet`, `coverage_matrix.parquet` |
| 0b | `src/pipeline/00b_enriched_harmonisation.ipynb` | 10 enriched tables into the same DB + the `*_enriched.parquet` exports |

Run both top to bottom, **0 before 0b** — 0b opens the database 0 creates and
fails fast with a clear message if it is missing. Both are idempotent
(`CREATE OR REPLACE`), so a partial run is safe to repeat.

Expected console figures, as measured on 2026-08-03:

```
00   hub 2,127 ACHs | gene 20,163 | lines with all 9 modalities: 237
00b  hub 2,145 ACHs (+18 from the newer sources)
     procan  948/948 lines resolved via the Sanger axis (100%)
     proteomics reach  375 -> 1,036 lines (+176%)
     v1 protein-coding inference vs HGNC: precision 95.3%, recall 99.9%
```

Wall time: ~4 min for step 0 (the 24 M-row `hpa_rna` load dominates), ~1 min for 0b.
The database lands at ~1.5 GB and is gitignored.

## Step 0c — `01_lookup_metadata_join.ipynb`

Runs unmodified. Verified against the v2 harmonisation: (2,127 × 26) written to
`harmonised_enriched.parquet`.

Its `<-- INVESTIGATE KEY MISMATCH` warning on `canonical_name` / `gap_class`
(85.5%) is **expected and not a fault**. The absolute match count is unchanged at
1,818; the rate fell only because the hub widened from 1,840 to 2,127 lines and
`union_lookup_entity.parquet` has no row for the 309 outside the DepMap roster.

## Step 0c-bis — `build_gene_regime.py` (run before Stage 4)

```bash
python src/pipeline/build_gene_regime.py --compare
```

Writes `outputs/gene_regime.parquet`, which until 2026-08-04 had **no producer** and
was the single blocker stopping the chain at Stage 5. Read the module docstring
before trusting the output: the classification rule is **reconstructed**, not
recovered, and the script documents both what is certain (the gene set and
`n_sensitive_gdsc` reproduce exactly) and what is inferred (the threshold).

**Ordering matters.** It reads `core_score.parquet`, and `build_cna_layer.py` reads
*it*, and Stage 4 reads `cna_flags.parquet`. The correct order is:

```
02_core_score  ->  build_gene_regime  ->  build_cna_layer  ->  04  ->  05  ->  ...
```

Running Stage 4 before rebuilding `cna_flags` silently uses the previous CNA layer.

## Step 0d — expose the outputs as SQL

```bash
python src/pipeline/build_warehouse_views.py --verify
```

Run **after** the stage that produced the parquet you want to query; it registers
views, so re-run it only when a file is added or removed, not when one is rebuilt.
Mind the casing trap documented in
[HARMONISATION_V2.md](HARMONISATION_V2.md#51-querying-it--build_warehouse_viewspy) —
use `v_score_annotated` rather than joining `out.*` to the warehouse by hand.

## Why there are two passes

`build_chronos_validation.py` reads `predictions_with_confidence.parquet` (it
needs `class == 'unknown'` and `core_score`), and `build_full_predictions.py`
reads `chronos_validation.parquet`. That is a genuine cycle in the existing
design. Because `core_score` itself changes in step 1, the CHRONOS correlations
have to be recomputed against the *new* scores, which means predictions are built
twice: once to carry the new `core_score` into a form CHRONOS can read, then
again to pick up the corrected `chronos_check` labels.

## Step 1 — `02_core_score.ipynb`

Run the notebook top to bottom. Four cells were added; the originals they
supersede were left in place and still execute first, which is intentional — each
patch cell is auditable against the version it replaces.

| cell | does |
|---|---|
| `Cell 3b` | UniProt→ENSG map: isoform suffixes (`p15941-2`) + `\|`-separated multi-accession genes |
| `Cell 7b` | `percentile_rank` ties → `method="min"` (zero-expression blocks to the floor) |
| `Cell 7c` | per-gene dynamic range → **new** `gene_dispersion.parquet` |
| `Cell 7d` | ProCan proteomics merged at the percentile level |

Outputs: `core_score.parquet` (rewritten), `gene_dispersion.parquet` (new).

Expected console figures, as verified against the real files in isolation:

```
genes with a protein layer     : 8,564 -> 9,698
protein-layer cell lines       :   375 -> 1,025   (CCLE 375 + ProCan 941, 291 shared)
ProCan lines mapped to ACH-    :   941 / 948
```

**Memory**: this is the heavy step. It holds the full expression matrix plus both
proteomics matrices. Close other Python processes first — the machine has 15.6 GB
and this session saw it drop below 1 GB free with one stray job running.

## Step 2 — predictions, pass 1

```bash
python src/pipeline/build_full_predictions.py
```

Fails fast with a clear message if `gene_dispersion.parquet` is missing, i.e. if
step 1 was skipped.

## Step 3 — CHRONOS, with the corrected effect-size floor

```bash
python src/pipeline/build_chronos_validation.py
```

Expect a large, intended reduction. The old floor (|ρ| > 0.10, ~1% of variance)
admitted almost anything, because at the median overlap of n=823 lines p<0.05 is
reached by |ρ| ≥ 0.068 — below the ρ floor, so significance was doing no work:

```
validated : 769  ->   29     (floor now |rho| >= 0.30, Cohen medium)
inverted  : 1258 ->   14
weak_positive / weak_negative : the 0.10-0.30 band, recorded but not acted on
```

## Step 4 — predictions, pass 2 (final)

```bash
python src/pipeline/build_full_predictions.py
```

Now picks up the corrected `chronos_check`. This is the artifact the CLI serves.

## Step 5 — evidence ledger

```bash
python src/pipeline/build_evidence_ledger.py
```

Requires the new `signal_spread` and `gene_top_vs_median_fold` columns from step
4, and gains `signal_spread` as a RANKING-tier ledger field.

## Step 6 — regression tests

```bash
python src/pipeline/rank_cell_lines.py
```

Six checks. Tests 1–4 are pre-existing (BRAF/A375 driver detection and top-5%
position, abundance_tracking pure-score sort, `rank_cell_lines`/`explain_pair`
agreement, `contradicted` agreement). Tests 5–6 are new and lock in the sort fix:

- **Test 5** — for `class='unknown'` genes, `core_score` must be non-increasing
  down the ranking. Catches any return of the driver hard-gate.
- **Test 6** — for `activation_driven` / `loss_of_function`, driver-positive rows
  must still sort above driver-negative ones. Catches over-application of the fix.

Both are structural over a sample of genes, not assertions about named genes, so
they keep their meaning as the gene set changes.

## What is NOT rebuilt by this runbook

`flags_with_driver.parquet`, `cna_flags.parquet` and the mutation/fusion layers
are untouched — no change in this batch affects them.

---

# Pending flow changes mandated by the §C.15–C.18 validation round

**Status: none of this is implemented.** The validation round wrote only
`test_run_*.py` diagnostics; no build script, notebook or export was modified.
The runbook above still rebuilds the pipeline as it currently stands. This section
records what the measurements say must change, so the two do not drift apart.

The single sentence version: **`core_score` moves from the sort step to the
exclusion step.** It gates (pAUC 0.593–0.601, replicated across two screens); it
does not rank (pAUC 0.500–0.518, at chance). See `MATH_REFERENCE.md` §C.15.

## Stage 2 — meaning changes, formula does not

`core_score.parquet` stays as it is. The percentile survived its challenge: it is
not an encoded detection boolean (§C.16, incremental depletion 0.810–0.827 across
six calibration rules), and the absolute per-line alternative that would have
replaced it lost under half-sample cross-validation. **No change to
`02_core_score.ipynb` is required.**

What changes is the label: the column is an *exclusion* statistic, not a ranking
statistic, and the docstrings and `provenance` notes should say so.

Protein's verdict depends on the metric (§C.18; `PROTEOMICS_DEFENCE.md` §6.2).
RNA alone beats every two-layer variant on **global AUROC** (−0.0086 GDSC /
−0.0042 Chronos). But on **pAUC** — the metric that matches a top-N recommender,
which is what Stage 7/8 actually expose to a user — the production two-layer
score *beats* RNA alone (+0.0234 GDSC / +0.0171 Chronos). Demoting proteomics to
display-only is therefore **not** indicated by the pAUC result; it would need to
be re-justified on a metric other than the one this product serves before that
change is made.

## Stage 3b — a new gate layer (does not exist)

Two independent exclusion rules, each with a measured scope:

| gate | basis | scope | caveat |
|---|---|---|---|
| A — abundance floor | bottom decile of `core_score` | all genes | floor effect not separable from a profiling confound, `p = 0.246` (§C.15) |
| B — homozygous deletion | `cna_flags` | **69 genes** with ≥10 deleted lines | residual 3.04% above label FPR 0.93%, `check_calls` (§C.18) |

Output shape: one row per `(gene, line)` carrying retained/excluded and which gate
fired. Gate B is not general — at 0.49 lines excluded per gene it is only usable
restricted to those 69 genes (CDKN2A 224, CDKN2B 202, MTAP 116).

## Stages 4–5 — two rules are contradicted by measurement

- `(regime_source=="measured") & (class=="abundance_tracking") & (n_layers==2) →
  "high"` at [build_full_predictions.py:151](../src/pipeline/build_full_predictions.py)
  and [05_confidence_tiers.ipynb:117](../src/pipeline/05_confidence_tiers.ipynb).
  §C.18 shows two-layer coverage is worse evidence on global AUROC but *better*
  evidence on pAUC (`PROTEOMICS_DEFENCE.md` §6.2) — so this promotes rows on a
  directionally justified basis, not a rejected one. It still has not been
  calibrated: no check exists that "high"-tier rows actually hit more often than
  "moderate"-tier ones.
- The tier ladder should carry the gate outcome, not only regime and CNA.

## Stage 7 — the sort key has no measured support

`explain_pair.sort_gene_rows` ([explain_pair.py:99-122](../src/pipeline/explain_pair.py))
uses `core_score` as primary key for `abundance_tracking` and as a tiebreak
elsewhere. With ranker-region pAUC at 0.500 there is no evidence for any ordering
*within* the retained set. The defensible presentation is the retained set as a
set, ordered by mechanistic flags where those exist.

`evidence_state.py` needs a state beyond
`NO_EVIDENCE / UNINFORMATIVE / CONFLICTING / RANKED` — something like `SCREENED`:
"these lines are excluded, for this reason; the rest are retained". That is what
the evidence supports. Note also that `MIN_LINES_FOR_RANKING = 30` fires for
exactly **one** gene in 19,177, so it is a correct safety rail doing no real work;
the gate that bites is `signal_spread == "narrow"`.

## Stage 8 — export contract

[export_web.py:233](../src/pipeline/export_web.py) does
`sort_values("core_score").head(top_n)`. The top-25 *ordering* is no longer
defensible as a ranking. The payload should carry the retained set plus the
excluded set with reasons, and the `provenance` note should state that the score
supports exclusion and not ordering.

Separately, and independent of all the above: `export_web.layer_magnitudes` already
reads `geo_expr` and `hpa_rna` and averages them into the **display** expression
track, while `core_score` uses neither. That gap between what the user sees and
what they are ranked by is currently undocumented in the payload.

## Rebuild order once these are written

```
02_core_score.ipynb          (only if protein is demoted)
  → build_gate_layer.py      (new — Stage 3b)
  → build_full_predictions.py (tier rules + gate outcome)
  → build_evidence_ledger.py
  → export_web.py
```

Nothing upstream of Stage 2 is affected; harmonisation and the flags do not move.
