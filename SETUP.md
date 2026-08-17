# GeneTraceAI — Teammate Setup

## What you need
- Python 3.11+ — https://www.python.org/downloads/
- Git — https://git-scm.com/downloads

---

## Step 1 — Clone the repo

```bash
git clone git@github.com:UoB-AstraZeneca-25-26/UoB-GeneTraceAI-25-26.git
cd UoB-GeneTraceAI-25-26
```

> Ask Chaithali to add your GitHub account to the organisation if clone fails.

---

## Step 2 — Run the setup script

**Windows (PowerShell):**
```powershell
.\setup.ps1
```

**Mac / Linux:**
```bash
bash setup.sh
```

The script installs Python packages and downloads the data files (~2.3 GB total).  
Grab a coffee — the database file is ~2 GB.

---

## Step 3 — Run the CLI

```bash
# Rank all cell lines for a gene
python final_pipeline/Ranking/cli.py gene BRAF

# Detail view — one gene + one cell line
python final_pipeline/Ranking/cli.py gene BRAF ACH-000219

# Multi-gene co-selection
python final_pipeline/Ranking/cli.py genes BRAF KRAS

# Exclusion query — high gene A, low gene B
python final_pipeline/Ranking/cli.py exclude BRAF EGFR
```

---

## Quick reference — useful genes to try

| Gene | Type | Try because |
|------|------|-------------|
| BRAF | Oncogene | Classic melanoma driver |
| KRAS | Oncogene | Most common cancer mutation |
| TP53 | TSG | Most frequently mutated gene |
| ERBB2 | Oncogene | Breast cancer amplification |
| EGFR | Oncogene | Lung cancer target |
| MYC  | Oncogene | Broad amplification across cancers |
