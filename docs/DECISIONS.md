# DECISIONS

Logged decisions about statistical and architectural choices, with rationale and
what was considered and rejected. Chronological within each section.

---

## Proteomics platform merge (2026-08-10)

### Context

Two proteomics datasets are used: ProCan (DIA-MS, 948 lines, absolute intensities)
and CCLE/Gygi (TMT, 375 lines, log-ratios to pooled reference). They use different
technologies and different measurement scales. `test_run_proteomics_platform_overlap.py`
measured per-protein Spearman rho across 291 shared lines: median rho = 0.373.

### Decision 1 — Apply protein-platform conflict tier filter (IMPLEMENTED)

**What**: For proteins where the two platforms actively disagree (rho < 0.3, the
"conflicting" tier: 2,246 proteins, 36.7%), do not blend the two percentile scores.
Use ProCan alone for those proteins on shared lines. ProCan is preferred because it
has more lines (948 vs 375) and uses DIA-MS (absolute, more sensitive).

**Why**: Cell 7e (before this fix) computed the mean of the two percentile scores
for all shared lines × shared proteins, including proteins where the platforms
disagree on the cell-line ordering. Averaging two disagreeing signals does not
produce a better signal — it produces a noisier one. The RNA pipeline already
applies an equivalent guard via Cochran's Q / I²; this brings proteomics into
parity.

**Implemented in**: `src/pipeline/02_core_score.ipynb` Cell 7f.
**Sidecar**: `src/pipeline/outputs/protein_platform_tier.parquet`

---

### Decision 2 — Restrict mean-blend to Q3/Q4 abundance proteins (IMPLEMENTED)

**What**: Even among non-conflicting proteins, restrict the mean-blend to the
top half of ProCan abundance (Q3/Q4, median abundance >= 50th percentile). Bottom-
half proteins (Q1/Q2) are reverted to single-platform (ProCan preferred) even if
their tier is "cautious" (rho 0.3–0.5).

**Why**: Step 6 of `test_run_proteomics_platform_overlap.py` showed a clear
abundance gradient in per-protein rho:
  - Q1 (low): median rho = 0.235 — below the "do not merge" threshold
  - Q2:        median rho = 0.322
  - Q3:        median rho = 0.430
  - Q4 (high): median rho = 0.511 — above the safe-merge threshold

The 0.373 overall median is pulled down by low-abundance proteins where DIA-MS
detects signal that TMT cannot reliably quantify relative to its pooled reference.
For those proteins, the two platforms are not measuring the same quantity; blending
them is not justified by the data.

**Combined rule after Decisions 1 + 2**:
  - Mean-blend kept: consistent/cautious tier AND Q3/Q4 abundance
  - Reverted to ProCan: conflicting tier OR Q1/Q2 abundance

**Implemented in**: `src/pipeline/02_core_score.ipynb` Cell 7f (same cell as Decision 1).

---

### Decision 3 — Do NOT gate on per-cell-line cross-protein agreement (NOT ACTIONED)

**What was found**: Step 5 of `test_run_proteomics_platform_overlap.py` measured
Spearman rho across proteins within each shared cell line (i.e., do the platforms
agree on which proteins are high vs low in a given cell line?). Median rho = 0.097;
all 291 shared lines fell below the 0.3 threshold.

**Why this was NOT actioned**: The pipeline's percentile ranking is gene-wise across
cell lines, not protein-wise within a cell line. The question the scoring asks is:
"Is this cell line in the top N% of lines for gene X?" — not "Is protein X more
abundant than protein Y in this line?". The Step 5 result answers the second
question, which is irrelevant to the scoring design.

Excluding lines based on the within-line cross-protein ρ would remove data for a
reason that does not apply to the pipeline's actual inference. The low within-line
rho is a known property of DIA-MS vs TMT: the two technologies have systematically
different relative protein-abundance orderings within a sample, driven by differences
in ionisation efficiency, fractionation, and reference-pool normalisation. This does
not affect gene-wise cross-line percentile rankings.

**What this finding IS useful for**: Any downstream analysis that uses the
proteomics layer to ask about relative protein abundance within a cell line (e.g.,
pathway enrichment using a single line's profile) should note this limitation. The
core scoring pipeline is not such an analysis.

---

## Tie handling in percentile normalisation (earlier, logged here for completeness)

**Decision**: Use `rank(pct=True, method='min')` throughout. Ties take the lowest
rank in the group (floor) rather than the midpoint (method='average').

**Why**: Tissue-restricted genes can have 30–60% of lines with exactly zero
expression. Under method='average', those lines score at the midpoint of the tie
block (e.g., 56.7% zeros → 28th percentile for KLK4). Under method='min', they
score at the floor (0.07th percentile). Absence of expression is the informative
signal for tissue-restricted genes; the floor is the correct placement.

**Implemented in**: `02_core_score.ipynb` Cell 7b (expr_pct) and Cell 9b (stratum_rank).

---

## S5 — "NAN" cell-line name bug: negative finding (2026-08-16)

### Investigation

The CoWork spec (v7.0 task S5) asked to confirm and fix a bug where a cell line
literally named "NAN" causes `pandas.read_csv` to coerce it to `NaN` (float),
losing the row. The fix was specified as `keep_default_na=False` at the point of
ingestion.

**Verdict: bug cannot be confirmed as described.**

Searches performed:

1. `sample_info` table in `celllineselector.db`: SQL query
   `WHERE upper(trim(cell_line_name)) = 'NAN'` — **0 rows**.
2. `cell_line_lookup.parquet`: pandas filter `cell_line_name.str.upper() == 'NAN'`
   — **0 rows**.
3. `data/DepMap_24Q4/Model.csv` (raw source): read with `keep_default_na=False`,
   searched all name/model columns for literal "NAN" — **0 rows**.

No cell line named "NAN" (or "nan", "Na", "NA") exists in any data source in
this project.

### Related finding: 12,751 NaN coercions in auxiliary fields

`data/DepMap_24Q4/Model.csv` read with default `keep_default_na=True` produces
**12,751 NaN cells** that are not NaN under `keep_default_na=False`. These come
from **auxiliary string columns** (ModelDerivationMaterial, ModelTreatment,
PatientCountry, etc.) where empty strings are coerced. **The model ID columns
and cell-line name columns are unaffected** — no model_id or cell_line_name
becomes NaN due to this coercion.

The harmonisation pipeline reads `Model.csv` via `data_utils.py`, which uses
`pd.read_parquet` (not `pd.read_csv`) on pre-cleaned parquet files. The coercion
therefore has no effect on any scored or ranked output.

### Decision

No fix is implemented: the bug as specified does not exist in this project's
data. If a future data update introduces a cell line with a pandas sentinel name
(NAN, NA, N/A, NULL, None, etc.), the correct fix is `keep_default_na=False` in
`src/scripts/data_utils.py` at the `pd.read_csv` call for Model.csv. The
relevant sentinel list is `pd.io.parsers.readers.STR_NA_VALUES`.
