# Confidence intervals & rank-stability for cell-line recommendations

This module adds an uncertainty layer on top of the `core_score` stage of
CellLineSelector. For any queried gene it produces, per cell line:

- a **95% confidence interval** on the `core_score`, and
- two **rank-stability** measures — how often a line lands in the top-10, and
  how often it is the single best line ("win-rate").

It does this by perturbing the scoring inputs (RNA and protein z-scores) by an
amount equal to their *measured* uncertainty, re-scoring and re-ranking under
that noise many times (Monte-Carlo), and summarising the spread. The intervals
are a statement of **reproducibility** — how much a recommendation moves under
realistic perturbation of its own inputs — not a calibrated probability that a
line is biologically the best model.

## Files

| File | What it is |
|---|---|
| `confidence_intervals_final.ipynb` | The method itself: CI/stability machinery, real-warehouse wiring (`run_gene_ci`), and the figures. |
| `ci_checks.py` | The verification & robustness suite (nine checks) that validates the method against the real warehouse. It does **not** define the method; it imports the notebook's names and tests them. |

## Requirements

- Python 3.10+
- `numpy`, `pandas`, `scipy`, `matplotlib`, `pyarrow` (for parquet)

```bash
pip install numpy pandas scipy matplotlib pyarrow
```

The **real-data path** additionally needs the CellLineSelector pipeline
environment on the import path:

- Modules: `config` (supplies the `RNA_Z`, `PROT_Z`, `PROT_TIER`, `CORE_SCORE`
  paths), `common` (`load_lineage`), and `protein_scorer`.
- An open warehouse connection, referred to as `con` (DuckDB).
- Pipeline outputs on disk: the core-score parquet plus the per-source exports
  `bulk_rna_z_by_source.parquet` and `bulk_prot_z_by_source.parquet` (read from
  next to `RNA_Z` / `PROT_Z`).
- **`ci_globals.json`** sitting next to `CORE_SCORE`. This holds three panel-wide
  scalars (global RNA SD, protein SD, protein-residual SD) that cannot be
  recovered from a single gene. It must be persisted once by `core_score.py`'s
  `run()`; see §2.3 of the notebook for the exact block to add. **The real-data
  path will not run until this file exists.**

## How to run

### 1. Real gene against the warehouse

First run the notebook's definition cells top-to-bottom (the constants/knobs cell
and the cells defining `_stratum_constants`, `build_anchor`, `core_score_chain`,
`monte_carlo_ci`, and the real-warehouse wiring). Then, with the pipeline
environment importable, a live `con`, and `ci_globals.json` in place:

```python
result, scores, ids = run_gene_ci(con, "EGFR", N=1000, k=10, seed=42)
result.head(10)   # core_score, ci_lo, ci_hi, ci_width, top1_winrate, ...
```

`run_gene_ci` reads the gene's leaves and noise, rebuilds the per-gene anchor,
attaches lineage, and runs the Monte-Carlo propagation. Increase `N` for
smoother win-rates (see the convergence check); `N=1000` is the working default.

### 2. Verification suite (`ci_checks.py`)

The checks reuse the notebook's functions, so make those names available first —
either run the notebook's definition cells in the same kernel, or import them:

```python
from confidence_intervals_final import (
    run_gene_ci, _gene_leaves, build_anchor, core_score_chain,
    monte_carlo_ci, _lookup, _resolve_gid, CORE_SCORE, RNA_Z,
    RNA_WINSOR_BOUND, PROT_STD_WINSOR_BOUND, PROT_RESID_WINSOR_BOUND,
    REGIME3_SIG_THRESH, RNA_BASE_FLOOR, PROT_BASE, RNA_TAU_PRIOR,
)
import ci_checks as chk
```

Most checks take either explicit callables/paths or a small `ctx` / `ns` dict of
those notebook names (each function's docstring lists exactly what it expects).
Typical calls, over a gene sample:

```python
genes = chk._sample_genes(CORE_SCORE, n=30, seed=0)

# 1. noise-free recompute must equal the stored core_score (the key check)
chk.centre_draw_identity(con, run_gene_ci, _resolve_gid, CORE_SCORE)

# 2-4. interval/ranking soundness and Monte-Carlo convergence
chk.order_consistency(con, genes, run_gene_ci)
chk.self_consistency(con, genes, run_gene_ci)
chk.convergence(con, "EGFR", run_gene_ci)

# 5. robustness to the noise knobs (needs the knob namespace)
chk.knob_sweep(con, genes, run_gene_ci, ns=globals())

# 6-8. variance split, regime flips, correlated-error sensitivity (need a ctx dict)
ctx = dict(_gene_leaves=_gene_leaves, build_anchor=build_anchor,
           core_score_chain=core_score_chain, _lookup=_lookup,
           CORE_SCORE=CORE_SCORE, RNA_WINSOR_BOUND=RNA_WINSOR_BOUND,
           PROT_STD_WINSOR_BOUND=PROT_STD_WINSOR_BOUND,
           REGIME3_SIG_THRESH=REGIME3_SIG_THRESH)
chk.variance_decomposition(con, "EGFR", ctx)
chk.regime3_flip(con, genes, ctx)
chk.correlated_error_sensitivity(con, genes, ctx)

# 9. is the score saturation the Φ transform or the winsor clip? (needs ns)
chk.winsor_test(con, genes, run_gene_ci, RNA_Z, _resolve_gid, ns=globals())
```

Each check prints its headline number; several also return a DataFrame for
inspection. The reference results obtained on the real pipeline are recorded in
the module docstring at the top of `ci_checks.py`.

## Reading the output

- **`core_score`** — the pipeline fit score (0–1); identical to production. It
  saturates toward 1.0 near the top (it is Φ of a combined z), so the interval
  and win-rate carry information the point score cannot.
- **`ci_lo`, `ci_hi`, `ci_width`** — the 95% interval and its width. Wide means
  the score would move a lot under realistic measurement noise.
- **`top1_winrate`** — how often (%) this line is the single best across all
  lines: the honest "is it *the* best?" number.
- **`topk_global_pct` vs `topk_lineage_pct`** — top-10 frequency against all
  lines vs. against the line's own tissue.

## Notes & gotchas

- The constants block in the notebook mirrors `core_score.py` exactly and must
  **not** be edited independently — the centre-draw identity check exists to
  catch any drift (target max|Δ| ≈ 1e-8).
- The only genuinely new parameters are the noise knobs; they are deliberately
  conservative and are stress-tested by `knob_sweep` and `winsor_test`.
- Single-source / single-platform lines have no observable disagreement, so
  their perturbation falls back to a corpus prior (RNA) or the gene-level
  platform tier (protein). Interval width therefore tracks *measured*
  disagreement, not raw source count.
