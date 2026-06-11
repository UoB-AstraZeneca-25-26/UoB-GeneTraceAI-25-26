# CellLineSelector — Work Log
Newest entries at the top. One entry per work session. Keep it factual; link commits/PRs where relevant.

<!-- Copy this block for each new entry:
## YYYY-MM-DD — <initials>
**Focus:** <task ID + short description, e.g. D4 — HPA RNA EDA>
**Done:**
-
**Decisions:**
-->

---

## 2026-05-28 — Triveni
**Focus:** D1 — First meeting with Daniel, Jens and AZ
**Done:**
- Attended the kickoff session with AZ; understood the problem statement and project needs.
- Debriefed with the team afterwards and scoped the initial setup.

---

## 2026-05-29 — Triveni
**Focus:** D2 — Topic exploration
**Done:**
- Studied cell biology and cell-line formation; built up background domain knowledge.

---

## 2026-05-30 — Triveni
**Focus:** D3 — Topic exploration (cont.)
**Done:**
- Continued background reading on cell biology and the problem domain.

---

## 2026-05-31 — Triveni
**Focus:** D4 — Project plan
**Done:**
- Drafted the project plan document for the team.

---

## 2026-06-01 — Triveni
**Focus:** D5 — Git repo setup
**Done:**
- Wrote an SOP for setting up the UoB/AZ GitHub organisation for the team.

---

## 2026-06-02 — Triveni
**Focus:** D6— Literature review
**Done:**
- Shortlisted 18 papers for the literature review.

---

## 2026-06-03 — Triveni
**Focus:** D7 — Initial data mining & understanding
**Done:**
- Downloaded all data files.
- Created the local environment for the project.
- Completed data preview: converted all text files to CSV, previewed all 14 files, and checked missing-value % across all 14.
- Compiled the findings into a document to share with the team.

---

## 2026-06-04 — Triveni
**Focus:** D8 — Git organisation live
**Done:**
- Git org up and running.
- Initialised the project board with Gantt chart, to-do list, and a progress board for the team to add their tasks.

---

## 2026-06-05 — Triveni
**Focus:** D9 — Data booklet
**Done:**
- Built the data booklet for the first 7 of the 14 tables/CSVs (column name, meaning, and an example/interpretation per column).

---

## 2026-06-06 — Triveni
**Focus:** D10 — Data booklet (cont.)
**Done:**
- Built the data booklet for the remaining 7 of the 14 tables/CSVs.
- Shared the completed booklet with the team for better understanding of the data.

---

## 2026-06-07 — Triveni
**Focus:** D11 — Deeper data understanding (01) — per-column audit
**Done:**
- Ran a per-column audit using an `audit()` helper (dtype, missing %, sample values per column).
- Produced the full missing-value summary table with a worst-column callout per dataset.

---

## 2026-06-08 — Triveni
**Focus:** 12 — Full EDA (1/3) — column decisions & structural fixes
**Done:**
- Documented all column-level decisions (keep / drop / rename / transform) across all 14 files — 459 total columns reduced to ~90 retained.
- Confirmed the GEO Info separator issue (file loaded as a single column); root cause identified as an incorrect delimiter on read.
- Flagged column-format differences: UniProt(SYMBOL) in proteomics vs SYMBOL(ENSG) in DepMap, and Ensembl version suffixes in the fusions file (e.g. ENSG00000075711.21).

---

## 2026-06-09 — Triveni
**Focus:** 13 — Full EDA (2/3) — gene universe
**Done:**
- Ran full gene-ID overlap analysis across DepMap, HPA RNA, and GEO expression.
- Classified all 30,231 DepMap-only genes by symbol pattern (uncharacterised loci, pseudogenes, lncRNA, snoRNA, snRNA, antisense RNA, ncRNA, miRNA, uncharacterised family, bare ENSG) — almost entirely non-protein-coding.
- Confirmed version suffixes are not the cause of the gene-overlap discrepancy — stripping them made zero difference to counts.
- Identified expression-scale inconsistency across sources: DepMap log2(TPM+1), HPA nTPM linear, GEO linear TPM/FPKM.

---

## 2026-06-10 — Triveni
**Focus:** D14 — Full EDA (3/3) — cell line universe & audit report
**Done:**
- Ran full cell-line-ID overlap analysis across all 6 data-bearing sources against the Sample Info master registry.
- Characterised HPA cell-line name matching: 862 exact, 63 via stripped name, 281 unmatched.
- Characterised GEO sample coverage: 3,267/3,267 in GEO Info, 3,159/3,267 with a direct Cellosaurus ID.
- Identified the PR- vs ACH- ID resolution dependency on OmicsProfiles for DepMap expression and mutations.
- Confirmed proteomics and metabolomics at 100% cell-line coverage with zero orphans.
- Identified 271 orphaned fusion ModelIDs and 251 orphaned signature ModelIDs not in Sample Info.
- Produced the full data audit report (12 sections): provenance, structure, identifiers, missing data, gene universe, cell line universe, expression scales, column reduction, tier classification, and open questions.