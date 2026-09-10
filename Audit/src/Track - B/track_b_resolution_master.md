# Track B — ID Resolution Master Document

**Owner:** Musa
**Track:** B — Cross-source transcriptomic expression
**Files:** `depmap_expr`, `hpa_rna`, `hpa_desc`, `geo_expr`
**Pipeline:** `src/notebooks/02b_key_resolution_trackB.ipynb`
**Final script:** `src/scripts/trackb_final_resolution_pipeline.py`
**Status:** Resolution complete for `hpa_rna`, `hpa_desc`, `geo_expr`. `depmap_expr` gene resolution complete; cell line resolution validated and documented, held as staging pending team alignment on PR- bridge approach.

---

## 1. What this work is and why it was needed

Before any of the four files on this track can feed into the shared feature store or the Neo4j knowledge graph, every row needs two things resolved: a canonical **gene key** (`ensg_id`, bare ENSG) and a canonical **cell line key** (`model_id`, ACH- format). These two keys are what allows DepMap expression data, HPA expression data, GEO expression data, mutation data, proteomics, and everything else on this project to be compared against each other as referring to the same genes in the same cell lines, regardless of how each original file happened to label them.

The four files on this track each arrive with their identifiers in completely different formats — embedded in column headers, as free-text names, as GSM sample accession codes — and none of them agree with each other's naming conventions. This document records what was found, what was fixed, and what the final outputs contain.

---

## 2. Shared lookup tables — validation findings

Before resolving any file, both shared lookup tables were validated. Several issues were found that affect every track, not just Track B.

### gene_lookup (19,446 rows)

| Check | Result |
|---|---|
| Primary key `ensg_id` unique | ✓ confirmed |
| Case convention | **UPPERCASE on delivery — must be lowercased on load** to match cleaned parquets |
| Coverage vs known protein-coding universe | **19,446 vs ~20,151–20,162 confirmed protein-coding genes in HPA** |
| Real genes missing from gene_lookup | **4,404 genes present in both HPA and GEO (confirmed protein-coding) are absent from gene_lookup** |

The 4,404 gap means any file filtered through `gene_lookup["ensg_id"]` silently loses real, legitimate protein-coding genes. This is not a per-file problem — it is a shared table quality issue. Logged to `logs/lookup_issue_genes_missing_from_lookup_depmap_expr.csv` and flagged to Track A.

### cell_line_lookup (1,840 rows)

| Check | Result |
|---|---|
| Primary key `model_id` unique | ✓ confirmed |
| Case convention | **UPPERCASE on delivery — must be lowercased on load** |
| Scope vs sample_info | Exact 1:1 copy of sample_info's model_id set (confirmed True, 0 asymmetric differences) — inherits all sample_info staleness issues |
| Duplicate CVCL accessions | **4 CVCLs claimed by 2 model_ids each** |
| Deprecated model_ids present | **11 rows** — see deprecation classification below |

**Duplicate CVCLs in cell_line_lookup:**

| CVCL | model_id 1 | name 1 | model_id 2 | name 2 |
|---|---|---|---|---|
| cvcl_0041 | ach-000833 | rh-30 | ach-001189 | none |
| cvcl_1122 | ach-001024 | chl-1 | ach-001041 | chl-1-dm |
| cvcl_1150 | ach-001737 | ctv-1-dm | ach-002222 | ctv-1 |
| cvcl_1337 | ach-001543 | kosc-2 | ach-002260 | none |

Logged to `logs/lookup_issue_duplicate_cvcl.csv`. Flagged to Track A.

### Deprecation classification (built once, applied to every file)

`sample_info.cellosaurus_issues` flags 11 model_ids as "removed from depmap." Cross-checking each against `signatures` (an independent DepMap-derived file) splits them into two categories with fundamentally different implications:

```
confirmed_gone  (3 IDs): flagged removed AND absent from signatures
                          → safe to exclude

disputed        (8 IDs): flagged removed BUT still active in signatures
                          → sample_info and signatures disagree
                          → do NOT silently drop
                          → hold out, flag to Track A
```

Disputed IDs logged to `logs/lookup_issue_disputed_deprecation.csv`. Full evidence in `flag_to_track_a_deprecation_disagreement.md`.

---

## 3. depmap_expr

### File structure
- **Raw shape:** 1,495 rows × 53,961 columns
- **Row identity:** PR- ProfileID (sequencing run — one level more specific than a cell line)
- **Column identity:** Gene symbol + ENSG embedded in header, e.g. `"TSPAN6 (ENSG00000000003)"`
- **Values:** Already log2(TPM+1) — no transform needed

### Gene resolution

Gene identity lives in column headers, not row values. Requires regex extraction:

```python
col_to_ensg = {col: re.search(r"ensg\d+", str(col)).group(0)
               for col in depmap_raw.columns
               if re.search(r"ensg\d+", str(col))}
```

| Result | Count |
|---|---|
| Total columns | 53,961 |
| Successfully extracted | 53,961 (0 unparseable) |
| Present in gene_lookup | 19,383 |
| Missing from gene_lookup | 34,578 — of which **4,404 are confirmed real per HPA/GEO cross-check** |

**Status:** gene_lookup gap flagged to Track A. Current output uses gene_lookup as-is; the 4,404 real genes being dropped are a known, documented limitation pending a corrected gene_lookup.

### Cell line resolution — two-hop chain

PR- IDs cannot be resolved directly to a cell line — they require a bridge through `depmap_profiles`:

```
PR- ProfileID (row index)
   → depmap_profiles (filter datatype == "rna") → modelid (ACH-)
   → cell_line_lookup.model_id (confirmed, + deprecation check)
```

**Why filter to `datatype == "rna"`:** `depmap_expr` is an RNA expression file; all its PR- IDs are exclusively RNA-datatype profiles. Confirmed empirically: 100% of `depmap_expr`'s PR- IDs appear only under `datatype == "rna"` in `depmap_profiles`, never under `wes` or `wgs`.

**Problem 1 — dual-profile cell lines:** 16 cell lines had two separate RNA ProfileIDs both mapping to the same `model_id` (two sequencing runs of the same biological sample). Without handling this, downstream output would contain duplicate (ensg_id, model_id) pairs. Fixed by keeping the profile with the highest total expression coverage (proxy for sequencing depth) and discarding the other.

**Problem 2 — 114 cell lines orphaned against cell_line_lookup:** After bridging PR- → ACH-, 114 model_ids found no match in `cell_line_lookup`. This is higher than the original EDA estimate of 67 because `cell_line_lookup` (being a 1:1 copy of sample_info) inherits sample_info's staleness — confirmed by comparing max ACH- IDs:

```
sample_info / cell_line_lookup max ID:  ACH-002926
depmap_profiles max ID:                 ACH-003161
signatures max ID:                      ACH-003480
```

**Recovery via Cellosaurus cross-reference:** Cellosaurus's `cross-references` field (parsed for `depmap; ach-XXXXXX` pattern) recovered **73 of 114** with full identity (name + CVCL). The remaining 41 split into:
- 25 tight block (ACH-003132–ACH-003161): exceed sample_info's ceiling entirely — too new for any current file. Absent from signatures AND Cellosaurus.
- 16 scattered gap: absent from sample_info but confirmed real via signatures; Cellosaurus cross-ref repeated same deprecated ID (circular, not new information).

**Problem 3 — deprecated model_ids in resolved rows:** Post-resolution audit found 7 resolved rows resting on deprecated model_ids. These split into:
- 1 `confirmed_gone` (no active alternative anywhere)
- 6 `disputed` (flagged removed in sample_info, but still active in signatures)

All 7 moved out of clean output and logged separately.

### Final cell line resolution counts

| Category | Count | Reason code |
|---|---|---|
| Clean resolved | 1,431 | model_id confirmed, not deprecated |
| Disputed (held out) | 7 | deprecated_disputed_vs_signatures |
| Confirmed gone | 0 | — |
| Orphaned — too new | 25 | unresolvable_too_new_confirmed_3way |
| Orphaned — stale snapshot | 16 | verified_real_missing_from_sample_info |
| Cellosaurus-recovered | 26 (of original 114) | cellosaurus_xref_fallback |
| **Total cell lines in depmap_expr** | **1,479** | after dual-profile dedup |

### Output status

```
depmap_expr_gene_resolved_STAGING.parquet
  → ensg_id resolved, profile_id retained as native key
  → STATUS: STAGING — cell line resolution pending team alignment
  → Shape: ~26M rows × 3 cols (ensg_id, profile_id, tpm_log2)

depmap_expr_cell_line_crosswalk.csv
  → 1,431 clean resolved model_ids, ready to merge
  → One .merge(on="profile_id") from final output once team aligned

logs/unmapped_cells_depmap_expr_orphans.csv
logs/disputed_deprecation_depmap_expr.csv
logs/confirmed_deprecated_depmap_expr.csv
```

---

## 4. hpa_rna

### File structure
- **Raw shape:** 24,315,372 rows × 6 columns
- **Format:** Already long — one row per (gene, cell line) pair
- **Row identity (gene):** `gene` column, already bare ENSG
- **Row identity (cell line):** `cell line` column, free-text hyphenated name
- **Values:** `tpm`, `ptpm`, `ntpm` (all three kept in output)

### Gene resolution

`gene` column is already clean bare ENSG — no extraction or translation needed. Validated against `gene_lookup`:

| Result | Count |
|---|---|
| Total unique genes | 20,162 |
| Present in gene_lookup | 19,366 |
| Missing from gene_lookup | 796 (confirmed real — same gene_lookup gap issue) |

### Cell line resolution — four-tier chain

**Key finding: hpa_desc promoted to Tier 1.** Confirmed 100% same-source name overlap between `hpa_rna["cell line"]` and `hpa_desc["cell line"]` (1,206/1,206 exact match). Since `hpa_desc` already carries Cellosaurus IDs (CVCL) for 99% of HPA's cell lines, going via CVCL directly is faster, more reliable, and more precise than name-matching against `sample_info` first.

**Original tier order (Integration Contract):**
1. Exact name → sample_info
2. Stripped name → sample_info
3. hpa_desc CVCL → sample_info (as afterthought)

**Corrected tier order:**
1. hpa_desc CVCL → cell_line_lookup (same-source, 100% name overlap confirmed)
2. Exact name → cell_line_lookup
3. Stripped name → cell_line_lookup
4. Cellosaurus name → cell_line_lookup (disambiguation-aware — see below)

**Critical bug found and fixed — bracket-collision in Cellosaurus name matching:**
Stripping bracket annotations for normalisation (e.g. `"BJ [human fibroblast]"` → `"bj"`) silently collapses multiple biologically distinct cell lines that share a base name. Found when three completely different cell lines (`BJ [human fibroblast]`, `BJ [human B-cell IHW]`, `BJ [human pancreatic adenocarcinoma]`) all collapsed to the same key. Audit of Cellosaurus confirmed **730 base names have this same collision risk**.

Fix: retain bracket content as a disambiguating signal. When multiple Cellosaurus entries share a base name, match on bracket content to select the correct one. If disambiguation still can't isolate a single entry, log as unresolved rather than guessing.

**Additional normalisation fixes:**
- Bracket content stripping added to normalise_base() — `[...]` removed before punctuation stripping
- `+` character added to stripped set (found in engineered subline names like `"BJ hTERT+ SV40 Large T+"`)

**Deprecation check on resolved rows:** All resolved model_ids checked against the three-way deprecation classification. 7 rows affected, same cell lines as depmap_expr (shared CVCL issue).

### Final cell line resolution counts

| Category | Count | % | Resolution method |
|---|---|---|---|
| Resolved — usable | 1,096 | 90.9% | hpa_desc_tier1 (1,103 before depr. check) |
| Identified, no DepMap equivalent | 97 | 8.0% | hpa_only_no_depmap_equivalent |
| Deprecated — disputed | 7 | 0.6% | deprecated_disputed_vs_signatures |
| Genuinely unresolved | 5 | 0.4% | unresolved |
| **Total** | **1,206** | **100%** | |

**The 5 genuinely unresolved, individually verified:**

| Name | Explanation |
|---|---|
| asc2telo differentiated | Engineered/differentiated derivative, no Cellosaurus entry |
| bj htert+ sv40 large t+ | Engineered subline, confirmed no match with `+` stripped |
| bj htert+ sv40 large t+ rasg12v | Further-engineered subline, confirmed no match |
| hhstec | Likely primary cell type, not an established catalogued line |
| hskmc | Likely primary cell type, not an established catalogued line |

Each checked against Cellosaurus and cell_line_lookup with full normalisation before being classified as unresolved.

### Output files

```
resolved/hpa_rna_resolved.parquet
  → 21,244,502 rows
  → columns: ensg_id, model_id, gene, gene name, cell line, tpm, ptpm, ntpm
  → (ensg_id and model_id added as new columns; all originals kept)

logs/disputed_deprecation_hpa_rna.csv
logs/no_depmap_equivalent_hpa_rna.csv
logs/unmapped_cells_hpa_rna.csv
```

---

## 5. hpa_desc

### File structure
- **Shape:** 1,206 rows × 7 columns (one row per HPA cell line)
- **No gene dimension** — output only needs `model_id`

### Resolution

`hpa_desc["cell line"]` confirmed identical to `hpa_rna["cell line"]` at 100% overlap. No separate resolution chain needed — the crosswalk built for `hpa_rna` is reused directly via a left merge.

```python
hpa_desc_resolved = hpa_desc.merge(
    hpa_crosswalk[["hpa_name", "model_id"]],
    left_on="cell line", right_on="hpa_name", how="left"
)
```

| Result | Count |
|---|---|
| Resolved with model_id | 1,104 |
| No model_id (no DepMap equivalent or deprecated) | 102 |

### Output files

```
resolved/hpa_desc_resolved.parquet
  → 1,104 rows
  → columns: model_id, cell line, disease, disease subtype,
             patient, primary/metastasis, sample collection site
  → (model_id added as first column; all originals kept)
```

---

## 6. geo_expr

### File structure
- **Raw shape:** 19,914 rows × 3,268 columns
- **Row identity (gene):** `gene` column, already bare ENSG
- **Column identity (cell line):** 3,267 GSM-prefixed column headers — no row-level cell line information exists

### Gene resolution

| Result | Count |
|---|---|
| Total unique genes | 19,914 |
| Present in gene_lookup | 15,988 |
| Missing from gene_lookup | 3,926 (same gene_lookup gap) |

### Cell line resolution — pre-exclusion mandatory

Two classes of GSMs must be excluded BEFORE any resolution attempt:

**Class 1 — null cellosaurus_id (108 GSMs):**
These have no cell line identity at all. Initially invisible to a standard `.isna()` check because they were stored as the literal text string `"none"` rather than a true pandas null — a gap in `_clean_str_cols()`, which only replaces `"nan"` with `pd.NA`, not `"none"`. Required explicit string check: `geo_info["cellosaurus_id"].astype(str).str.lower() == "none"`.

**Class 2 — cvcl_7082 (478 GSMs):**
These carry a CVCL value but every supporting field (`cellline`, `gse_id`, `matching_type`) is null — no evidence for how the CVCL was determined. An unverifiable identity assertion is treated the same as a missing one: excluded, not trusted.

Total excluded before resolution: **586 of 3,267 (17.9%)**.

**Three candidate fallback name columns in geo_info — measured, not assumed:**
`geo_info` has three different cell-line-name columns with very different match rates:

```
cellline           60.5% vs sample_info.cell_line_name
cell_line          22.8%
cell_line_trimmed  14.9%
```

`cellline` confirmed as the correct fallback. In practice this fallback was never needed — every usable GSM resolved via `cellosaurus_id` directly once the 586 excluded rows were removed.

**Resolution chain:**
```
GSM column header → geo_info.geo_accession → cellosaurus_id
                  → cell_line_lookup.rrid → model_id
                  (+ deprecation check)
                  (+ title-field disambiguation for duplicate CVCLs)
```

**Ambiguous duplicate CVCL cases found:**

**CVCL_1150 (CTV-1 / CTV-1-DM):** Two real, biologically distinct cell lines sharing one CVCL. Resolved using `geo_info.title` field — GEO's own submitter text confirmed this specific GSM is `"ctv-1"` (not `-dm`) → resolved to `ach-002222`.

**CVCL_0041 (RH-30):** Two candidate model_ids both flagged "removed from depmap" in `sample_info.cellosaurus_issues`. Checked against `signatures` — both absent, confirming genuine retirement. Classified `cvcl_exists_but_only_deprecated_ach_ids`.

**Discovery: sample_info deprecation flag unreliable (significant finding):**
Checking the deprecation flag across all resolved rows found 14 rows across geo_expr and the other two crosswalks where model_ids were flagged "removed from depmap" — but 6 of those 8 unique cell lines were still present and active in `signatures`. This is a direct contradiction between two DepMap-derived files. Full investigation documented in `flag_to_track_a_deprecation_disagreement.md`.

**Multi-study cell lines:** 276 of 586 resolved unique cell lines appear in more than one GSE study. These will need within-study z-score then cross-study averaging at the transform stage. `gse_id` is carried through the crosswalk specifically for this reason.

### Final resolution counts

| Category | Count | Resolution method |
|---|---|---|
| Clean resolved | 2,294 | resolved / resolved_via_title |
| Disputed (held out) | 13 | deprecated_disputed_vs_signatures |
| No DepMap equivalent / deprecated | 373 | no_depmap_equivalent / cvcl_exists_but_only_deprecated_ach_ids |
| Ambiguous — genuinely unresolvable | 1 | ambiguous_unresolved |
| Pre-excluded | 586 | null_cvcl / cvcl_7082_unverified |
| **Total** | **3,267** | **reconciliation confirmed** |

### Output files

```
resolved/geo_expr_resolved.parquet
  → 36,676,472 rows
  → columns: ensg_id, model_id, gsm, rma_intensity
  → (wide-to-long melt required; gsm kept as native audit trail)

logs/excluded_cells_geo_expr.csv
logs/disputed_deprecation_geo_expr.csv
logs/no_depmap_equivalent_geo_expr.csv
logs/ambiguous_unresolved_geo_expr.csv
```

---

## 7. Cross-file issues flagged to Track A

All flags documented in full in their respective files. Summary:

| Issue | Files affected | Document |
|---|---|---|
| sample_info is a stale snapshot — missing 41+ valid model_ids present in depmap_profiles/signatures | depmap_expr | flag_to_track_a_sample_info.md |
| cell_line_lookup has 4 duplicate CVCL accessions (CTV-1/CTV-1-DM, RH-30/None, CHL-1/CHL-1-DM, KOSC-2/None) | all files | logs/lookup_issue_duplicate_cvcl.csv |
| sample_info.cellosaurus_issues deprecation flag disagrees with signatures for 8 of 11 flagged model_ids | depmap_expr, hpa_rna, geo_expr | flag_to_track_a_deprecation_disagreement.md |
| gene_lookup has 19,446 genes vs ~20,151–20,162 confirmed protein-coding — 4,404 real genes silently dropped | depmap_expr, hpa_rna, geo_expr | logs/lookup_issue_genes_missing_from_lookup_*.csv |
| _clean_str_cols() null handling only replaces "nan" not "none" — affected geo_info.cellosaurus_id | geo_expr | documented in geo_expr writeup |
| Both lookup tables delivered UPPERCASE — must be lowercased on load to match cleaned parquets | all files | fixed in pipeline script |

---

## 8. Key methodological findings (reusable by other tracks)

**Bracket-collision in Cellosaurus name matching (730 affected base names):** Naive name-normalisation that strips bracket annotations can silently merge multiple biologically distinct cell lines. Any Cellosaurus-based fallback matching must retain bracket content for disambiguation before resolving.

**"none"-as-string null (affects _clean_str_cols()):** The shared cleaning function only converts `"nan"` to `pd.NA`. The literal text `"none"` (and potentially `"na"`, `""`, `"null"`) survives unchecked. Any `.isna()` relying on this function may undercount nulls silently.

**Three-way deprecation classification:** The `cellosaurus_issues` flag in `sample_info` should never be applied as a binary drop. It must be cross-referenced against `signatures` and split into `confirmed_gone` (safe to exclude) vs `disputed` (flagged removed but still active elsewhere — hold out and escalate). Treating it as a reliable binary flag resulted in removing cell lines that are provably still active in the DepMap ecosystem.

**Lookup table case convention:** The shared `gene_lookup` and `cell_line_lookup` tables are UPPERCASE; all cleaned parquets are lowercase. Every join between them must include a lowercasing step or it will silently return zero matches with no error.

**hpa_desc as Tier 1 for hpa_rna:** The Integration Contract originally proposed this as a late fallback. Empirical verification (100% same-source name overlap, 99% CVCL coverage) confirms it should be the primary resolution path, not a fallback — it's faster, more direct, and avoids the name-normalisation ambiguity of sample_info matching entirely.

---

## 9. Overall resolution summary

| File | Resolved rows | Primary output | Status |
|---|---|---|---|
| depmap_expr | 1,431 clean cell lines (gene: 19,383 genes) | depmap_expr_gene_resolved_STAGING.parquet | Cell line pending team alignment |
| hpa_rna | 1,096 cell lines × 19,366 genes = 21.2M rows | hpa_rna_resolved.parquet | Complete |
| hpa_desc | 1,104 cell lines | hpa_desc_resolved.parquet | Complete |
| geo_expr | 2,294 GSMs × 15,988 genes = 36.7M rows | geo_expr_resolved.parquet | Complete |

---

## 10. Pending before Track B is fully closed

1. **Team alignment on depmap_expr cell line resolution** — specifically: confirm PR- → ACH- bridge via depmap_profiles is the agreed approach, and confirm the highest-coverage tiebreak rule for the 16 dual-profile cell lines. Once confirmed, the final output is one `.merge()` call using the already-built `depmap_expr_cell_line_crosswalk.csv`.

2. **Updated gene_lookup from Track A** — once a corrected gene_lookup is published (covering the full ~20,151 protein-coding genes rather than 19,446), re-run the final assembly for all four files to recover the 4,404 currently-dropped real genes.

3. **Track A decision on disputed deprecation IDs (8 model_ids)** — once Track A confirms whether `sample_info.cellosaurus_issues` or `signatures` is the authoritative deprecation signal, the 7 disputed rows in `depmap_expr` and the 7 in `hpa_rna` can be either restored to the clean crosswalk or permanently excluded.
