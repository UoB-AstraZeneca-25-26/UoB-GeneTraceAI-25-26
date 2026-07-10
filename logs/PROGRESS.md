# GeneTraceAI — Project Progress Log

Newest entries at the top. Detailed per-track documents live in their track folders.

---

## Quick links to detailed track docs

| Track | Detailed document |
|---|---|
| Track A — Cell line identity spine | [`src/Track - A/Track A - Updates.md`](../src/Track%20-%20A/Track%20A%20-%20Updates.md) |
| Track B — Transcriptomic expression | [`src/Track - B/track_b_resolution_master.md`](../src/Track%20-%20B/track_b_resolution_master.md) |
| Track B → Track A flag (sample_info staleness) | [`src/Track - B/flag_to_track_a_sample_info.md`](../src/Track%20-%20B/flag_to_track_a_sample_info.md) |
| Track C — Mutations, fusions, signatures, miRNA | [`src/Track - C/track-C.md`](../src/Track%20-%20C/track-C.md) |
| Integration contract | [`integration_contract.md`](../integration_contract.md) |
| Lookup table how-to | [`reference/work_with_lookup.md`](../reference/work_with_lookup.md) |

---

## Track status summary (as of 2026-07-09)

| Track | Owner | Files | Status |
|---|---|---|---|
| A | Chaithali | cellosaurus, sample_info, depmap_profiles | Complete — cell_line_lookup v2 shipped (union_lookup) |
| B | Musa | depmap_expr, hpa_rna, hpa_desc, geo_expr | hpa_rna, hpa_desc, geo_expr complete; depmap_expr staging pending team alignment |
| C | — | mutations, fusions, signatures, miRNA | In progress |
| D | — | proteomics, metabolomics, geo_info | In progress |

---

## Cross-track open issues

| Issue | Raised by | Status |
|---|---|---|
| `sample_info` is a stale DepMap snapshot — 41 valid model_ids missing | Track B | Partially resolved via Cellosaurus xref (+26); 41 remain |
| `gene_lookup` has 19,446 genes vs ~20,151 confirmed protein-coding — 4,404 real genes silently dropped | Track B | Open — awaiting updated gene_lookup from Track A |
| `sample_info.cellosaurus_issues` deprecation flag disagrees with `signatures` for 8 of 11 flagged model_ids | Track B | Open — Track A needs to determine authoritative deprecation signal |
| `cell_line_lookup` has 4 duplicate CVCL accessions | Track B | Logged to `src/Initial Data Audit and Mapping/logs/lookup_issue_duplicate_cvcl.csv` |
| `_clean_str_cols()` does not replace literal `"none"` strings with pd.NA | Track B | Open — affects geo_info.cellosaurus_id and potentially others |

---

## Work log — Triveni (project setup & EDA)

### 2026-06-12 — Triveni
**Focus:** D16 — Developing Data Cleaning Pipeline
**Done:**
- Developed a data cleaning pipeline for stripping whitespace, double spaces, and lowercasing across all 14 tables.

---

### 2026-06-11 — Triveni
**Focus:** D15 — Team Meet
**Done:**
- Completed Level 0 EDA; compiled findings into a shared [Data Architecture Sheet](https://docs.google.com/spreadsheets/d/1YAw23fSNvva7y_1YMlRYrkLlrCURqWsdQak9FKL5bPs/edit?usp=sharing).

---

### 2026-06-10 — Triveni
**Focus:** D14 — Full EDA (3/3) — cell line universe & audit report
**Done:**
- Full cell-line-ID overlap analysis across all 6 data-bearing sources vs Sample Info master registry.
- Characterised HPA cell-line name matching: 862 exact, 63 via stripped name, 281 unmatched.
- Characterised GEO sample coverage: 3,267/3,267 in GEO Info, 3,159/3,267 with direct Cellosaurus ID.
- Identified the PR- vs ACH- ID resolution dependency on OmicsProfiles for DepMap expression and mutations.
- Confirmed proteomics and metabolomics at 100% cell-line coverage with zero orphans.
- Identified 271 orphaned fusion ModelIDs and 251 orphaned signature ModelIDs not in Sample Info.
- Produced the full data audit report (12 sections).

---

### 2026-06-09 — Triveni
**Focus:** D13 — Full EDA (2/3) — gene universe
**Done:**
- Full gene-ID overlap analysis across DepMap, HPA RNA, and GEO expression.
- Classified all 30,231 DepMap-only genes by symbol pattern — almost entirely non-protein-coding.
- Confirmed version suffixes are not the cause of gene-overlap discrepancy.
- Identified expression-scale inconsistency: DepMap log2(TPM+1), HPA nTPM linear, GEO linear TPM/FPKM.

---

### 2026-06-08 — Triveni
**Focus:** D12 — Full EDA (1/3) — column decisions & structural fixes
**Done:**
- Documented all column-level decisions across all 14 files — 459 total columns reduced to ~90 retained.
- Confirmed GEO Info separator issue (incorrect delimiter on read).
- Flagged column-format differences: UniProt(SYMBOL) in proteomics vs SYMBOL(ENSG) in DepMap.

---

### 2026-06-07 — Triveni
**Focus:** D11 — Per-column audit
**Done:**
- Ran per-column audit: dtype, missing %, sample values per column.
- Produced full missing-value summary table with worst-column callout per dataset.

---

### 2026-06-06 — Triveni
**Focus:** D10 — Data booklet (remaining 7 files)
**Done:** Built and shared the data booklet for the remaining 7 of 14 tables/CSVs.

---

### 2026-06-05 — Triveni
**Focus:** D9 — Data booklet (first 7 files)
**Done:** Built the data booklet for the first 7 of 14 tables/CSVs.

---

### 2026-06-04 — Triveni
**Focus:** D8 — Git organisation live
**Done:** Git org up and running; project board initialised with Gantt chart, to-do list, and progress board.

---

### 2026-06-03 — Triveni
**Focus:** D7 — Initial data mining & understanding
**Done:** Downloaded all data files; created local environment; previewed all 14 files; checked missing-value % across all 14.

---

### 2026-06-02 — Triveni
**Focus:** D6 — Literature review
**Done:** Shortlisted 18 papers for the literature review.

---

### 2026-06-01 — Triveni
**Focus:** D5 — Git repo setup
**Done:** Wrote SOP for setting up the UoB/AZ GitHub organisation.

---

### 2026-05-31 — Triveni
**Focus:** D4 — Project plan
**Done:** Drafted the project plan document for the team.

---

### 2026-05-29/30 — Triveni
**Focus:** D2/D3 — Topic exploration
**Done:** Studied cell biology and cell-line formation; built up domain knowledge background.

---

### 2026-05-28 — Triveni
**Focus:** D1 — Kickoff meeting with Daniel, Jens and AZ
**Done:** Attended kickoff session; debriefed with team and scoped initial setup.
