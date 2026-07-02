# GeneTraceAI — Data Preparation Pipeline

This folder contains the notebooks and scripts for loading, cleaning, exploring, and integrating the 14 multi-omics datasets used in the GeneTraceAI project.

---

## Setup

```bash
pip install -r requirements.txt
```

---

## Folder Structure

```
src/
├── notebooks/
│   ├── 00_data_loading.ipynb
│   ├── 01_data_cleaning.ipynb
│   ├── 02_eda.ipynb
│   └── 03_data_integration.ipynb
├── scripts/
│   ├── data_utils.py
│   ├── cleaning_functions.py
│   └── export_data.py
└── requirements.txt
```

---

## Data Flow

```
UoB-GeneTraceAI-25-26/data/              (raw source files — gitignored)
        |
        v  00_data_loading.ipynb
UoB-GeneTraceAI-25-26/data/parquet/raw_data/      (14 raw parquets — gitignored)
        |
        v  01_data_cleaning.ipynb
UoB-GeneTraceAI-25-26/data/parquet/data_clean/    (15 cleaned parquets — gitignored)
        |
        v  03_data_integration.ipynb
UoB-GeneTraceAI-25-26/outputs/processed/          (15 CSVs — gitignored)
UoB-GeneTraceAI-25-26/outputs/clean_data_report.xlsx  (tracked — small)
```

## What is tracked in git

```
UoB-GeneTraceAI-25-26/
├── src/              ← all notebooks and scripts (tracked)
├── outputs/
│   └── clean_data_report.xlsx  (tracked)
├── README.md         (tracked)
└── .gitignore        (tracked)
```

Everything under `data/`, all `.parquet` files, and `outputs/processed/` are gitignored because they are too large for git.

---

## Notebooks

### `00_data_loading.ipynb`
Reads all 14 raw source files from `C:/Disertation/Data/`, converts `.txt` files to `.csv`, performs initial structure and format checks (shape, dtypes, missing %), and saves every dataset as a Parquet file under `data/parquet/raw_data/`. Run this first.

### `01_data_cleaning.ipynb`
Loads the raw parquets and applies a cleaning function to each dataset (imported from `scripts/cleaning_functions.py`). Cleaned tables are saved as Parquet files under `data/parquet/data_clean/`. Includes a post-cleaning shape and missing value summary.

### `02_eda.ipynb`
Loads the cleaned parquets and performs exploratory data analysis: per-column quality audits (missing %, unique counts, sample values), gene ID overlap across transcriptomics datasets (HPA, GEO, DepMap), and cell line ID overlap across all datasets (ACH-, CVCL, cell line names).

### `03_data_integration.ipynb`
Loads the cleaned parquets and discovers join connections between tables by comparing actual values in candidate key columns. Confirms pairwise linkages (ACH- IDs, PR- IDs, ENSG IDs). Exports the final processed CSVs to `data/processed/` and generates a per-column quality Excel report.

---

## Scripts

### `scripts/data_utils.py`
Utility functions shared across all notebooks.
- `preview(df, name)` — prints shape, dtypes, missing % and first 5 rows.
- `load_raw_parquets()` — loads all 14 raw parquet files into a dict.
- `load_clean_parquets()` — loads all 15 cleaned parquet files into a dict.
- `save_dataframes_to_csv(tables, out_dir)` — saves a dict of DataFrames to CSV.

### `scripts/cleaning_functions.py`
One cleaning function per dataset, plus a `clean_all()` convenience wrapper. Every function applies the same baseline steps and any dataset-specific transformations:

| Dataset | Key cleaning steps |
|---|---|
| `hpa_rna` | Lowercase columns/values, strip Ensembl version suffix from `gene` |
| `depmap_expr` | Parse `GENE (ENSGXXX)` column headers → bare ENSG IDs |
| `geo_expr` | Lowercase, strip Ensembl version suffix from `gene` |
| `proteomics` | Parse `UNIPROT (GENE)` headers → bare UniProt IDs; returns a `protein_map` lookup table |
| `fusions` | Split `gene1(ENS ID)` / `gene2(ENS ID)` into separate gene name and ENSG ID columns |
| `mutations` | Lowercase, strip Ensembl version suffix from `ensemblgeneid` / `ensemblfeatureid` |
| `cellosaurus` | Rename verbose key columns to `cellosaurus_cell_line_name` / `cellosaurus_accession` |
| `depmap_profiles` | Lowercase columns and values |
| `sample_info` | Lowercase columns and values |
| `geo_info` | Lowercase, strip whitespace (fixes silent join failures on `cell_line` field) |
| `hpa_desc` | Lowercase columns and values |
| `metabolomics` | Lowercase columns and values |
| `mirna` | Lowercase, strip trailing decimal suffix from `name` column |
| `signatures` | Rename `unnamed: 0` → `signature_index`, lowercase |

All functions also: lowercase all column names, lowercase all string values, strip/collapse extra whitespace, and replace `"nan"` strings with `pd.NA`.

### `scripts/export_data.py`
Functions for saving and reporting processed data.
- `save_cleaned_parquets(cleaned)` — saves all cleaned DataFrames to `data/parquet/data_clean/`.
- `build_table_summary(df)` — returns a per-column quality DataFrame (missing %, unique count, sample values, primary key flag).
- `export_clean_report(tables)` — writes an Excel file with one quality summary sheet per dataset.
- `export_processed_csvs(tables)` — exports all cleaned DataFrames as CSVs to `data/processed/`.

---

## Datasets

| # | File | Content |
|---|---|---|
| 1 | `1_4_hpa_rna_celline.tsv` | HPA RNA expression by cell line (TPM) |
| 2 | `2_DepMap_OmicsExpressionAllGenesTPMLogp1Profile.csv` | DepMap RNA expression matrix (TPM log+1) |
| 3 | `3_GEOexpression.csv` | GEO microarray expression across 3,267 samples |
| 4 | `4_Harmonized_MS_CCLE_Gygi_subsetted.csv` | CCLE/Gygi mass spectrometry proteomics |
| 5 | `5_OmicsFusionFilteredSupplementary.csv` | DepMap gene fusions |
| 6 | `6_OmicsSomaticMutationsProfile.csv` | DepMap somatic mutations (~1M rows) |
| 7 | `7_cellosaurus.csv` | Cellosaurus cell line dictionary (152k entries) |
| 8 | `8_DepMap_OmicsProfiles.csv` | DepMap profile-to-model bridge table (PR- → ACH-) |
| 9 | `9_DepMap_sample_info.csv` | DepMap sample metadata and lineage annotations |
| 10 | `10_GEOInfo.txt` | GEO series/sample metadata |
| 11 | `11_hpa_rna_celline_description.tsv` | HPA cell line disease and origin annotations |
| 12 | `12_CCLE_metabolomics_20190502.csv` | CCLE metabolomics profiling (225 metabolites) |
| 13 | `13_CCLE_miRNA_20181103.gct` | CCLE miRNA expression (GCT format) |
| 14 | `14_OmicsGlobalSignatures.csv` | DepMap global omics signatures (MSI, CIN, ploidy) |
