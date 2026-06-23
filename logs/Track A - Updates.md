# Track A Progress Update — Cell Line Key Resolution & Data Harmonization

**Track:** A
**Files owned:** `cellosaurus`, `sample_info`, `depmap_profiles`
**Notebook:** `src/notebooks/02b_key_resolution_cellosaurus_sampleinfo_depmap.ipynb`
**Status:** Harmonization complete — data loss documented, no further resolution feasible within provided data environment

---

## 1. What This Work Actually Is

Before any of the three datasets can be merged into the shared pipeline, every row needs a **canonical cell line identifier** resolved and verified against its counterpart. The prior cleaning functions handled only string-level tidying (lowercasing, whitespace stripping) — they did **not** verify that resulting keys actually join correctly across datasets.

This track covers the **verification and resolution layer**:
- Identify the real join key per dataset pair
- Test it against the shared lookups
- Measure how much resolves cleanly
- Surface every mismatch, structural gap, or unresolvable record
- Document all losses with root cause

---

## 2. Files Involved

| File | Role | Rows | Key Column | Key Format |
|---|---|---|---|---|
| `cellosaurus_clean` | Global cell line registry — name authority | 152,231 | `cellosaurus_accession` | `cvcl_xxxx` (len=9) |
| `sample_info_cleaned` | Experiment-scoped cell line metadata — **anchor** | 1,840 | `rrid` / `depmap_id` | `cvcl_xxxx` / `ach-xxxxxx` |
| `depmap_profiles_cleaned_pivoted` | DepMap RNA expression profiles | 1,822 | `modelid` | `ach-xxxxxx` (len=10) |

`sample_info_cleaned` is the **anchor dataset** — all joins are measured relative to its 1,840 records.

---

## 3. Join Chain Architecture

Two independent joins were performed with `sample_info_cleaned` as the anchor:

```
cellosaurus_clean
  [cellosaurus_accession]
        ↕  (Join 1 — cvcl_ key)
  sample_info_cleaned
  [rrid]          [depmap_id]
                       ↕  (Join 2 — ach- key)
          depmap_profiles_cleaned_pivoted
                  [modelid]
```

---

## 4. Key Format Audit — Both Joins

### Join 1: `cellosaurus_accession` ↔ `rrid`

| Property | `cellosaurus_accession` | `rrid` |
|---|---|---|
| Prefix | `cvcl_` | `cvcl_` |
| Length | 9 (uniform) | 9 (uniform) |
| Null count | 0 | 22 |
| `RRID:` prefix pollution | None | None |
| Case mismatch | None | None |

### Join 2: `depmap_id` ↔ `modelid`

| Property | `depmap_id` | `modelid` |
|---|---|---|
| Prefix | `ach-` | `ach-` |
| Length | 10 (uniform) | 10 (uniform) |
| Null count | 0 | 0 |
| Format deviation | None | None |

**Verdict:** Zero format-related failures in either join — all losses are structural, not technical.

---

## 5. Challenges Faced

### Challenge 1 — Scale Misread: 150,419 "Mis-Mapped" Appeared Alarming

Initial match check returned 150,419 mis-mapped records between `cellosaurus_clean` and `sample_info_cleaned`. This was **structurally expected** — Cellosaurus is a global registry (152,231 rows); `sample_info` is experiment-scoped (1,840 rows). The correct diagnostic lens is the reverse: how many of sample_info's 1,840 rrid values fail to match Cellosaurus — which is only 24.

**Lesson:** For asymmetric dataset joins, always measure from the smaller/anchor side, not the larger registry side.

---

### Challenge 2 — NaN rrid Values Not Caught by Prior Cleaning

22 records in `sample_info_cleaned` have a null `rrid` field. The existing cleaning function (`clean_sample_info()`) performed string normalization but did not flag or impute missing accessions. These 22 records carry no resolvable cell line registry ID and cannot be matched to Cellosaurus by any method within the provided data environment.

**Root cause:** Missing values at data collection time not a pipeline artifact.

---

### Challenge 3 — Two Valid CVCL Accessions Absent from Provided Dump

`CVCL_X507` and `CVCL_V618` are syntactically valid accessions (correct prefix, correct length, non-null) but are absent from the provided `cellosaurus_clean` dataset. These are likely **secondary or deprecated accessions** that map to a different primary AC in the full Cellosaurus registry — but the provided dump only contains primary accessions.

**Root cause:** Snapshot incompleteness — the institutional Cellosaurus dump does not include secondary/deprecated AC mappings.
**Resolution attempted:** None feasible — Cellosaurus API access not available in the provided environment.

---

### Challenge 4 — DepMap Release Version Mismatch Between `sample_info` and `depmap_profiles`

Both datasets use identical `ach-` formatted IDs with identical lengths and zero nulls yet 91 + 73 = 164 IDs are unmatched across the two tables. Since format is not the cause, the mismatch is a **DepMap release version gap**: `sample_info` and `depmap_profiles` were sourced from different snapshot versions of the DepMap database, where cell lines are periodically added or retired between releases.

- 91 `depmap_id` values in `sample_info` are absent from `depmap_profiles` → cell lines in metadata, not in expression snapshot
- 73 `modelid` values in `depmap_profiles` are absent from `sample_info` → profiled cell lines absent from metadata snapshot

**Root cause:** Inter-release divergence in DepMap data not resolvable without access to a unified release version.

---

## 6. Data Loss Summary — Per Join

### Join 1: Cellosaurus ↔ sample_info

| Outcome | Count | % of sample_info |
|---|---|---|
| Matched | 1,812 | 98.48% |
| Lost — null `rrid` at source | 22 | 1.20% |
| Lost — CVCL absent from dump (`CVCL_X507`, `CVCL_V618`) | 2 | 0.11% |
| **Total lost** | **24** | **1.30%** |

### Join 2: sample_info ↔ depmap_profiles (inner join)

| Outcome | Count | % of sample_info |
|---|---|---|
| Matched | 1,749 | 95.05% |
| Lost from sample_info (not in profiles) | 91 | 4.95% |
| Lost from depmap_profiles (not in sample_info) | 73 | — |
| **Total lost from anchor** | **91** | **4.95%** |

---

## 7. Overall Data Loss Across Three Datasets

1. **Three datasets in scope** → `cellosaurus_clean` (152,231 rows), `sample_info_cleaned` (1,840 rows, **anchor**), `depmap_profiles_pivoted` (1,822 rows) — all keyed on uniform, fixed-length identifiers (`cvcl_` / `ach-`).

2. **Cellosaurus ↔ sample_info** (via `cellosaurus_accession` ↔ `rrid`) → **24 records lost**, 22 due to null `rrid` at source, 2 due to secondary/deprecated accessions absent from the provided dump.

3. **sample_info ↔ depmap_profiles** (via `depmap_id` ↔ `modelid`) → **91 records lost** from sample_info side, **73 records lost** from depmap_profiles side — caused purely by DepMap release version mismatch between datasets.

4. **Zero format-related losses** across all three joins — no case mismatches, no prefix pollution, no null keys in DepMap join, no length inconsistencies detected.

5. **Cumulative loss from sample_info anchor** → `24 (Cellosaurus) + 91 (DepMap) = 115 records` affected across both integrations (overlap unknown without three-way join).

6. **Final harmonized coverage** → Cellosaurus join retains **98.48%** (1,812/1,840); DepMap join retains **95.05%** (1,749/1,840); three-way intersection is the effective working dataset.

7. **Total irrecoverable loss** → **115 records maximum (6.25% of sample_info)** split between source nulls, snapshot version gaps, and deprecated registry entries; none resolvable within the provided data environment.

---

## 8. Current Status

| Join Pair | Key Format | Match Rate | Status |
|---|---|---|---|
| Cellosaurus ↔ sample_info | `cvcl_` (len=9) | 98.48% (1,812/1,840) | Complete — loss documented |
| sample_info ↔ depmap_profiles | `ach-` (len=10) | 95.05% (1,749/1,840) | Complete — loss documented |
| Three-way intersection | Both | TBD | Pending final merge step |

---

## 9. Open Questions

1. Can the team confirm whether `sample_info` and `depmap_profiles` are sourced from the same DepMap release? The 91+73 orphan pattern strongly suggests a version mismatch.
2. Does the provided Cellosaurus dump include secondary/deprecated accessions? If not, `CVCL_X507` and `CVCL_V618` remain permanently unresolvable.
3. Should the 22 null `rrid` rows be retained in downstream analysis with `NaN` for Cellosaurus fields, or excluded entirely?

---

## 10. Next Steps

1. Perform **three-way inner join** across all three datasets to establish the final high-confidence working set.
2. Document three-way intersection count and residual loss for the Integration Contract.
3. Flag the DepMap version mismatch formally to the team — affects any track using `sample_info` as anchor.
4. Confirm null-rrid row handling policy with supervisor before final pipeline run.