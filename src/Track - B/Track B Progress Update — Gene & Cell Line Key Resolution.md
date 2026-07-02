---
title: Track B Progress Update — Gene & Cell Line Key Resolution

---

# Track B Progress Update — Gene & Cell Line Key Resolution
**Owner:** Musa
**Track:** B — Cross-source transcriptomic expression
**Files owned:** `depmap_expr`, `hpa_rna`, `hpa_desc`, `geo_expr`
**Status:** In progress — identity resolution underway, GEO transform pipeline not yet built

---

## 1. What this work actually is

Before any of the four files I own can be merged into the shared dataset, every row needs two things resolved: a canonical **gene ID** (bare ENSG) and a canonical **cell line ID** (ACH- `model_id`, cross-referenced to CVCL). The cleaning functions already written by the team (`cleaning_functions.py`) only did string-level tidying — lowercasing, whitespace stripping, version-suffix removal. They did **not** verify that the resulting keys actually join correctly against the shared lookup tables (`sample_info`, `depmap_profiles`, `geo_info`, `cellosaurus`).

My job on this track is that verification and resolution layer: for each of my four files, identify the real join key, test it against the shared lookups, measure how much resolves cleanly, and surface every mismatch, duplicate, or structural gap I find — fixing what's local to my own files, and flagging anything that touches a shared lookup table back to whoever owns it.

This work lives in a new notebook: `src/notebooks/02b_key_resolution_trackB.ipynb`, sitting between the team's `02_eda.ipynb` and `03_data_integration.ipynb` in the existing pipeline sequence.

---

## 2. Files involved and why each one matters

### Files I own directly

| File | Role | Why it's on my track |
|---|---|---|
| `depmap_expr` | DepMap RNA-seq expression matrix (TPM log2+1) | Primary RNA expression backbone — largest gene coverage of the three expression sources |
| `hpa_rna` | Human Protein Atlas RNA expression (TPM/nTPM) | Second RNA source — independent cross-check against DepMap |
| `hpa_desc` | HPA cell line metadata (disease, tissue, CVCL accession) | Exists purely to help resolve `hpa_rna`'s cell line identities when direct name matching fails |
| `geo_expr` | GEO microarray expression across 3,267 samples | Third, supplementary expression source — different platform (microarray, not RNA-seq) |

### Shared lookup files I read but do not own (read-only reference)

| File | Why I need it | Who owns it |
|---|---|---|
| `depmap_profiles` | Bridges `depmap_expr`'s PR- ProfileIDs to ACH- ModelIDs — without this, DepMap's row index is meaningless | Track A |
| `sample_info` | The canonical cell line registry — every resolved ID ultimately gets checked against this | Track A |
| `geo_info` | Bridges `geo_expr`'s GSM column headers to cell line identity (CVCL) — `geo_expr` has no cell line column of its own at all | Track D (per current diagram) |
| `cellosaurus` | Name-authority fallback for cell line synonyms not found anywhere else | Track A |

I read these tables to verify my own joins. I do not edit their cleaning logic — anything I find broken inside them gets written up and flagged back to the owner, not silently patched.

---

## 3. Gene key resolution — what I checked

**Canonical key:** bare ENSG (no version suffix), e.g. `ensg00000146648`

All three of my expression files carry gene identity differently in their raw form:

- `depmap_expr` — gene ID is embedded in **column headers** (`"TSPAN6 (ENSG00000000003)"`), already partially extracted by the team's `clean_depmap_expr()` regex
- `hpa_rna` — gene ID is a **row value** in a `gene` column, already in clean ENSG format from source
- `geo_expr` — gene ID is also a **row value**, clean ENSG format, no version suffixes

What I did: confirmed each file's gene column/headers actually match the expected `ensg\d+` pattern after cleaning, checked for any columns/rows that the existing regex-based cleaning silently failed to catch, and measured pairwise and three-way gene overlap across `depmap_expr`, `hpa_rna`, and `geo_expr` to establish the cross-source "high-confidence gene set" — the genes present in all three sources, which is the most reliable subset for downstream confidence scoring.

---

## 4. Cell line key resolution — the four different chains

This was the bulk of the work, because each file arrives with a completely different native cell line identifier.

### depmap_expr — two-hop bridge
```
PR- ProfileID (row index)
   → depmap_profiles (filter Datatype == "rna")
   → ACH- ModelID
   → confirm against sample_info.depmap_id
```
Checked match rate at each hop, and specifically checked for the known issue of cell lines with two RNA profiles (multiple PR- IDs pointing at the same ACH-), which need a tie-break decision before collapsing to one row per cell line.

### hpa_rna — three-tier name resolution
```
Tier 1: normalised name → exact match vs sample_info.cell_line_name
Tier 2: normalised name → match vs sample_info.stripped_cell_line_name
Tier 3: normalised name → look up CVCL via hpa_desc → match vs sample_info.rrid
Unresolved → logged, not dropped
```
This is where most of the debugging time went (see Section 5). The existing `clean_hpa_rna()` function only lowercased cell line names — it did not strip hyphens or slashes, so names like `"23132/87"` survived cleaning with punctuation intact and would fail an exact-match join against `sample_info` unless properly normalised first.

### geo_expr — header-based resolution (no cell line column exists)
```
GSM column header → geo_info.geo_accession → cellosaurus_id → CVCL
```
`geo_expr` is structurally different from the other two: there is no cell line column anywhere in the file. The identity lives entirely in the 3,267 GSM-prefixed column headers, and is completely unresolvable without joining `geo_info` first.

---

## 5. Challenges faced

### Challenge 1 — `geo_expr`'s cleaning function does far less than it appears to

The README describes `geo_expr` cleaning as "lowercase, strip Ensembl version suffix from gene." On inspection, that is **literally all it does**. The cleaned output is still in raw wide format (genes as rows, 3,267 GSM columns), with completely untransformed values (raw microarray intensities, no log transform, no z-scoring). This file needs a near-complete rebuild of its cleaning step — transpose/melt to long format, `log2(x+1)`, and critically, **z-scoring within each GSE study separately**, because this data is RMA-normalised microarray data from 18 different studies with real batch effects between them, not RNA-seq. None of this exists yet. This is the single largest piece of remaining work on my track.

### Challenge 2 — relative path errors from notebook location

Initial attempts to load parquet files failed with `FileNotFoundError` because the notebook's working directory (`src/notebooks/`) is two levels below the `data/` folder at the repo root. Relative paths needed `../../` prefixes, or better, routing through the shared `load_clean_parquets()` utility in `data_utils.py` instead of hardcoding paths — which is more robust if the data folder structure ever changes.

### Challenge 3 — inconsistent normalised-column naming across dataframes caused a silent join failure

While building the three-tier `hpa_rna` resolution chain, I named the normalised cell line column `cell_line_norm` on the `hpa` dataframe, but `cln_norm` on `sample_info` and `hpa_desc`. The merge step (`on="cln_norm"`) failed with `KeyError: 'cln_norm'` because that exact column name never existed on the `hpa`-derived dataframe — it only existed on the other two. This passed an earlier isolated test (which only checked `sample_info` against `hpa_desc` directly) because that test never touched the inconsistently-named column. The bug only surfaced once the full three-tier chain ran end to end. **Lesson:** testing pipeline pieces in isolation can pass while the full chain still has a bug — the connections between pieces need their own explicit check.

### Challenge 4 — Jupyter cell execution order caused a `KeyError` that looked like a missing column but wasn't

A separate `KeyError: 'cln_norm'` occurred even after the column genuinely existed in the dataframe, traced to having reloaded `hpa_desc` fresh via `pd.read_parquet()` in a later cell without re-running the cell that added the normalised column — silently reverting to a version of the dataframe without it. Resolved by restarting the kernel and running all cells top to bottom, and by auditing the notebook for duplicate `pd.read_parquet()` calls on the same variable name.

### Challenge 5 — placeholder column names in early draft code needed verification against real data

Draft resolution code initially assumed column names like `"cellosaurus id"` and `"cell line"` in `hpa_desc` without confirming them against the actual cleaned parquet. This is a general risk when writing cross-file join logic from memory of EDA done earlier in the project — every column reference needs to be confirmed with `.columns.tolist()` against the live file before being trusted, not assumed from documentation or prior exploration.

### Challenge 6 — `geo_info` has three candidate fallback name columns with very different match rates

Cross-checking `geo_info` (not my file, but required for `geo_expr` resolution) revealed it has three different cell-line-name columns (`cellline`, `cell_line`, `cell_line_trimmed`) with substantially different match rates against `sample_info` (60.5%, 22.8%, and 14.9% respectively). The resolution chain doesn't currently specify which one to use at the fallback step — this is now a flagged question back to whoever owns `geo_info`, since the choice materially affects how many GEO samples resolve.

### Challenge 7 — a likely sample_info version mismatch affecting orphan rates

Comparing match rates across multiple datasets against `sample_info` suggests `sample_info` may be a slightly different DepMap release than `depmap_profiles` and other DepMap-origin files — IDs that fail to match `sample_info` often do exist in `depmap_profiles` or `signatures`. This isn't something I can fix from my track alone since I don't own `sample_info`, but it directly affects my `depmap_expr` resolution numbers, so it's flagged to Track A.

---

## 6. Current status

| File | Gene key | Cell line key | Status |
|---|---|---|---|
| `depmap_expr` | Verified clean | Two-hop chain tested, dual-profile cell lines identified, awaiting Track A clarification on sample_info version mismatch | Resolution logic done, blocked on shared-table question |
| `hpa_rna` | Verified clean | Three-tier chain built and debugged, match rates measured | Resolution logic done, needs cleaning function extended (hyphen/slash stripping, log2 transform, drop tpm/ptpm) |
| `hpa_desc` | N/A (metadata only) | Confirmed as working fallback source for hpa_rna's Tier 3 | Verified, no further work needed on this file itself |
| `geo_expr` | Verified clean | Resolution chain identified (GSM → geo_info → CVCL) | **Not yet built** — needs transpose, melt, log2 transform, and z-score-per-GSE-study, none of which exist in current cleaning function |

---

## 7. Open questions to raise with the team

1. Is `geo_info` actually owned by Track D as shown in the latest diagram, or still Track B? Need to confirm before I either fix it myself or wait on someone else.
2. Which of `geo_info`'s three candidate name columns should be the standard Tier-3 fallback for any chain using it? My data points to `cellline` as the strongest option.
3. Can Track A confirm whether `sample_info` is the same DepMap release as `depmap_profiles` and `signatures`? The orphan-rate pattern across multiple datasets suggests a version mismatch.
4. Should I extend `clean_hpa_rna()` and rebuild `clean_geo_expr()` directly in the shared `cleaning_functions.py`, or keep my fixes local to my own notebook until reviewed?

---

## 8. Next steps

1. Rebuild `geo_expr` cleaning: transpose/melt to long format, drop excluded GSM columns (null CVCL, and the unresolved `cvcl_7082` anomaly block), apply `log2(x+1)`, z-score within each GSE study, collapse multi-study cell lines to one row each.
2. Extend `clean_hpa_rna()`: drop `tpm`/`ptpm`, keep only `ntpm`, apply `log2(ntpm+1)`, and fix cell line name normalisation to strip hyphens/slashes, not just lowercase.
3. Finalise and document match-rate numbers for all four files once the above fixes are in place.
4. Write up findings for the shared Integration Contract document, including the `geo_info` fallback-column recommendation and the `sample_info` version-mismatch flag for Track A.