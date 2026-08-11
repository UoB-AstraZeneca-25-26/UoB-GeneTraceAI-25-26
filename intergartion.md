# Integration Contract
**Project:** CellLineFinder  
**Status:** Draft — awaiting Track A sign-off  
**Date:** [today's date]

---

## 1. Canonical Column Schemas

### gene_lookup
| Column         | Type   | Format rule                                      |
|----------------|--------|--------------------------------------------------|
| ensg_id        | string | Bare ENSG, no version suffix — PRIMARY KEY       |
| hgnc_symbol    | string | Current approved symbol                          |
| prev_symbols   | string | Pipe-delimited list                              |
<<<<<<< HEAD
<<<<<<< HEAD
| alias_symbols  | string | Pipe-delimited list                              |   
=======
| alias_symbols  | string | Pipe-delimited list                              |
>>>>>>> a6ea596 (New Folder structure with requirements and readme files)
=======
| alias_symbols  | string | Pipe-delimited list                              |   
>>>>>>> c2ee8d3 (updated the track A harmonization according to the integration file design system.)
| hgnc_id        | string |                                                  |
| biotype        | string | Always protein_coding                            |

### cell_line_lookup
| Column                  | Type   | Format rule                           |
|-------------------------|--------|---------------------------------------|
| model_id                | string | ACH-format — PRIMARY KEY              |
| cell_line_name          | string |                                       |
| stripped_cell_line_name | string |                                       |
| cvcl_accession          | string | CVCL-format                           |
| synonyms                | string | Pipe-delimited list                   |
| rrid                    | string |                                       |

---

## 2. Non-Negotiable Rules

- Primary keys (`ensg_id`, `model_id`) are **never null, never duplicated**
- All list-type columns use pipe `|` as delimiter — no commas, no semicolons
- Nothing is silently dropped — every unresolved row is written to an unmapped log

---

## 3. Gene Key — Resolution Chain

**Dependency:** `depmap_profiles` table (bridge table linking PR- → ACH-)  
must be loaded before this chain can run.
**Canonical key:** `ensg_id` (bare ENSG, no version suffix e.g. ENSG00000141510)

All ENSG version suffixes (e.g. `.15`, `.21`) must be stripped at cleaning time.

**Proven gene coverage from EDA:**
- DepMap: 53,961 genes (includes non-coding — filtered to ~20,000 after biotype check)
- HPA: 20,162 genes
- GEO: 19,914 genes
- High-confidence set (all three sources): **16,060 genes**
- Fusions covered by DepMap: 25,222/25,424 (99.2%)
- Mutations covered by DepMap: 19,635/19,798 (99.2%)

**Format by dataset:**

| Dataset    | Raw format                  | Cleaning action                          |
|------------|-----------------------------|------------------------------------------|
| depmap_expr | `TSPAN6 (ENSG00000000003)` | Extract ENSG from brackets, strip suffix |
| hpa_rna    | Bare ENSG                   | Strip version suffix only                |
| geo_expr   | Bare ENSG                   | Strip version suffix only                |
| fusions    | `ENSG00000075711.21`        | Strip version suffix                     |
| mutations  | Bare ENSG in EnsemblGeneID  | Strip version suffix                     |

---

## 4. Cell Line Key — Resolution Chains

**Canonical key:** `model_id` (ACH-format, sourced from sample_info.DepMap_ID)

### Chain 1 — Direct ACH- match
**Datasets:** proteomics, fusions, signatures, metabolomics

These datasets already carry ACH-format IDs. Match directly to `sample_info.DepMap_ID`.

| Dataset      | Match rate          | Orphaned |
|--------------|---------------------|----------|
| Proteomics   | 375/375 (100%)      | 0        |
| Metabolomics | 927/927 (100%)      | 0        |
| Fusions      | 1,428/1,699 (84%)   | 271      |
| Signatures   | 1,704/1,955 (87%)   | 251      |

Known gaps: Fusions (271 orphaned) and Signatures (251 orphaned) 
 have elevated orphan rates. These IDs exist in the dataset but have no match in sample_info. Flagged for Track A to investigate whether these are deprecated cell line IDs or a sample_info version mismatch.
No fallback needed. Unmatched rows (if any) go to unmapped log.

---

### Chain 2 — PR- profile ID → ACH- via depmap_profiles
**Datasets:** depmap_expr (RNA), mutations (WES/WGS)

PR- IDs are profile-level, not cell-line-level. Must be resolved first.

**Resolution order:**
1. Look up ProfileID in `depmap_profiles`
   - For depmap_expr: filter `Datatype == "rna"`
   - For mutations: filter `Datatype` in `["wes", "wgs"]`
2. Get ModelID (ACH-) from result
3. Match ModelID to `sample_info.DepMap_ID`
4. If no match at any step → write to unmapped log

| Dataset              | Match rate          | Orphaned |
|----------------------|---------------------|----------|
| depmap_expr (RNA)    | 1,412/1,479 (95%)   | 67       |
| Mutations (WES/WGS)  | 1,697/1,744 (97%)   | 47       |

---

### Chain 3 — Cell line name → ACH- via name matching
**Datasets:** hpa_rna

**Resolution order:**
1. Exact match vs `sample_info.cell_line_name` → resolves **862/1,206**
2. Match vs `sample_info.stripped_cell_line_name` → resolves additional **63**
3. Look up name in `hpa_desc` to get CVCL accession → match vs `sample_info.RRID` → resolves additional **179**
4. If still unmatched → write to unmapped log

| Step                   | Resolved | Cumulative |
|------------------------|----------|------------|
| Exact name match       | 862      | 862        |
| Stripped name match    | 63       | 925        |
| CVCL → RRID fallback   | 179      | 1,104      |
| Permanently unresolved | 102      | logged     |

**Total recovery: 1,104/1,206 (91.5%)**

---

### Chain 4 — GSM accession → cell line → ACH-
**Datasets:** geo_expr

GSM IDs are column headers in geo_expr. Must be mapped via geo_info first.

**Resolution order:**
1. Match GSM column name to `geo_info.geo_accession` → **3,267/3,267 (100%)**
2. If row has `Cellosaurus_ID` (3,159 GSMs) → match vs `sample_info.RRID`
3. If row has `cell_line` name only (2,329 GSMs) → match vs `sample_info.cell_line_name`, then `stripped_cell_line_name`
4. If unresolved → write to unmapped log

Note: 0 unmapped GSMs at step 1 — all GSM columns are accounted for in geo_info.

---

### Chain 5 — CCLE-style column headers → ACH-
**Datasets:** mirna

miRNA has ~954 cell-line columns in format `CELLNAME_TISSUE` (e.g. `DMS53_LUNG`).

**Resolution order:**
1. Strip everything after first underscore → get base name (e.g. `DMS53`)
2. Match base name vs `sample_info.stripped_cell_line_name`
3. If unmatched → write to unmapped log

---

## 5. Unmapped Log Specification

Every dataset's cleaning function must write unresolved rows to:
`/data/logs/unmapped_{dataset_name}.csv`

Columns: `original_value | dataset | resolution_step_reached | reason`

 ⚠️ Match rate for miRNA CCLE columns not yet quantified.
---
