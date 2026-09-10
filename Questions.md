# questions.md — Team Q&A Log (GeneTraceAI / CellLineFinder)

> Cross-track questions raised during development, with owner and resolution.
> Team: Thiruvel · Musa ·
>       Chaithali · Triveni ·Daniel (supervisor)
>
> STATUS: [RESOLVED] [OPEN] [FOR-DANIEL] [FOR-TEAM]
> RAISED-BY / OWNER attributed per question.

---

## 1.Achitecture questions.

### Should miRNA contribute to gene-level scoring?
- **Raised by:**Thiruvel (regrading the signifance)· **Owner:** Chaithali → Daniel
- **Status:** [FOR-DANIEL]
- `mirna_model_level.parquet` has no `ensg_id` axis. Needs miRTarBase target
  mapping + gene-class gating.
  [Resolved] Trivani's new harmonisation method attacheds both Gene Axis and Model ID's

### Flat merge vs. hub-and-spoke?
- **Raised by:** Chaithali · **Owner:** Triveni
- **Status:** [RESOLVED] By Chaithali
- Flat merge (~24M rows, mostly-null master cols) corrupts grain. Hub-and-spoke
  with canonical keys locked. Datasets joined through lookup tables, not merged.

### Adding new data to the new Architecture
- **Raised by:** Chaithali . **Owner:** Team
- **Status** [OPEN] By Chaithali
- Adding MiRNA and metabolomics is not possible as the datasets have no gene axis information and cannot be pulled into the Layer 1 architecture.I tried 2 different methods to make it work, one where i tried connecting by "model_ID", but that just results in equal value to all the genes.
---

## 2. Identity & Coverage



### Is gene_lookup missing ~4,400 genes?
- **Raised by:** Musa · **Owner:** Thiruvel 
- **Status:** [OPEN]
- Direct challenge to completeness of the table the universe filter depends on.

### Four CVCL accessions each claimed by two model_ids?
- **Raised by:** Triveni, confirmed w/ Musa · **Owner:** Chaithali
- **Status:** [RESOLVED]
- Duplicate accessions in cell_line_lookup.They contain many to many relations.

### Why is Track A's harmonised output unused across all tracks?
- **Raised by:** Chaithali (cross-track review) · **Owner:** Team
- **Status:** [OPEN]
- `trackA_harmonization.parquet` (1,840×13, with coverage flags) sitting unused.
  Raised as live architectural question.
- The issue is to be considered null as the harmonisation superseds the question
---

## 3. Validation & Calibration

### Which datasets validate the scoring stack, and how?
- **Raised by:** Chaithali · **Owner:** Chaithali → Daniel (scope sign-off)
- **Status:** [RESOLVED — pending Daniel confirm as primary metric]
- Chronos = oncogene arm (PR-AUC). GDSC = TSG/fusion arm (osimertinib/EGFR,
  olaparib/BRCA1/2). Anchor pairs (BRAF/A375, KRAS/HCT116) = calibration, NOT
  validation.
-*Architecture in Pause

---

## 4. Team / Process

---

## 5. Brief-Requirement Gaps

