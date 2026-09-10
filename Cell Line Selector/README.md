# Cell Line Selector

Cell Line Selector is a multi-omics data engineering and analysis project designed to transform heterogeneous biomedical sources into a consistent, auditable, human cell-line-focused resource for downstream cell-line selection and comparison.

The project integrates transcriptomic, proteomic, metabolomic, mutation, fusion, copy-number, and cell-line metadata into a common analytical framework. It includes dataset-specific cleaning, identifier harmonisation, species filtering, genomic annotation, quality audits, DuckDB-based storage, and statistical integration steps.

---

## Project overview

The repository is structured around a staged pipeline:

1. Raw data discovery and loading
2. Dataset-specific cleaning
3. Cell-line and gene harmonisation
4. Database loading and validation
5. Exploratory data analysis and figure generation
6. Transcriptomic and protein scoring integration

The end goal is to create analysis-ready data tables and summary resources that can be queried consistently across multiple molecular layers.

---

## Repository structure

### Root-level folders and files

- `data/`  
  Project data storage. Typically contains raw and cleaned data outputs, as well as intermediate tabular data products.

  - `data/raw data/`  
    Original or minimally processed source files.

  - `data/clean data/`  
    Cleaned datasets after source-specific preprocessing.

  - `data/clean_csv/`  
    Cleaned CSV exports for inspection and downstream reuse.

- `db/`  
  Persistent analytical database storage. The main project database is `db/celllineselector.duckdb`.

- `figs/`  
  Figure outputs and numerical audit summaries. Includes publication-ready plots and audit tables such as `numerical_audit.csv`.

- `geo_meta/`  
  GEO sample metadata files (for example `GSM*.txt`) used for parsing sample information and processing metadata.

- `logs/`  
  Pipeline logs and audit artifacts generated during execution.

- `pipelines/`  
  Executable pipeline scripts that define the analysis workflow across data loading, cleaning, harmonisation, EDA, transcriptomics, and protein scoring.

- `reports/`  
  Reports, figures, and summary outputs generated during exploratory analysis and validation.

- `results/`  
  Results from statistical and integration stages, including calibration outputs, model summaries, and selected parameter files.

- `src/`  
  Shared source code for data utilities, cleaning, harmonisation, transcriptomics, proteomics, EDA, and DuckDB access.

- `tests/`  
  Validation tests for key logic, especially transcriptomic gene filtering and Ensembl identifier filtering.

- `content.txt`  
  Repository transcript-style summary of key database tables and records.

- `project_contributions.txt`  
  Detailed project contribution narrative and technical summary.

- `protein-coding_gene.txt`  
  Protein-coding gene reference information used in validation and filtering steps.

- `requirements.txt`  
  Core requirement list for the project environment.

- `setup_env.sh`  
  Environment bootstrap script that creates a local virtual environment and installs the dependencies.

- `pytest.ini`  
  Pytest configuration for reliable repository-root import resolution.

---

## Source code organisation

### `src/scripts/`
This directory contains the reusable processing code used across the project.

- `cleaning_functions.py`  
  Core cleaning layer for dataset-specific transformations. Handles GEO, HPA, Cellosaurus, fusions, model metadata, and non-human filtering logic.

- `data_utils.py`  
  Shared paths and dataset utilities for consistent file access across the project.

- `duckdb_utils.py`  
  Functions for connecting to and populating the DuckDB database.

- `eda_functions.py`  
  Utility functions for profiling, missingness analysis, numeric summaries, and exploratory visualisations.

- `export_data.py`  
  Export routines for writing cleaned or harmonised tables to disk.

- `harmonisation_functions.py`  
  Identifier harmonisation logic for cell lines, genes, Ensembl IDs, UniProt mappings, and cross-source bridges.

- `logging_utils.py`  
  Logging helpers used to record pipeline execution and anomalies.

- `proteonomics_m1.py`  
  Protein-scoring method 1 implementation.

- `proteonomics_m2.py`  
  Protein-scoring method 2 implementation.

- `transcriptomics.py`  
  Canonical transcriptomic processing logic, including scale harmonisation, detection-threshold handling, and source integration.

---

## Pipeline scripts

The `pipelines/` directory contains the workflow steps that orchestrate the project.

- `00_data_loading.py`  
  Finds raw input files, loads them into a common representation, and prepares intermediate Parquet outputs.

- `01_data_cleaning.py`  
  Executes dataset-specific cleaning routines and produces cleaned datasets.

- `02_data_harmonisation.py`  
  Builds gene rosters, cell-line rosters, identifier bridges, and harmonised tables.

- `03_load_duckdb.py`  
  Loads cleaned outputs into DuckDB for queryable analysis.

- `04_eda.py`  
  Produces missingness reports, statistical summaries, and exploratory plots.

- `05_transcriptomics_stale.py`  
  Older or stale transcriptomic implementation retained for reference and comparison.

- `05_transcriptomics_stats_layer.py`  
  Main transcriptomic statistical integration workflow.

- `06_protein_score M1.py`  
  Protein scoring pipeline using Method 1.

- `06_protein_score M2.py`  
  Protein scoring pipeline using Method 2.

---

## Data products and core outputs

The project creates a suite of tables and intermediate products that support analysis and validation.

### Core harmonised tables

- `cell_line_roster`  
  One-row-per-model representation of cell-line identifiers and associated metadata.

- `gene_roster`  
  One-row-per-gene reference table used for identifier harmonisation and filtering.

- `geo_info` and `geo_expr`  
  GEO sample metadata and expression matrices.

- `hpa_rna` and `hpa_desc`  
  Human Protein Atlas RNA expression and descriptive metadata.

- `depmap_expr` and `depmap_profiles`  
  DepMap expression matrix and model metadata.

- `cellosaurus`  
  Cellosaurus annotation data for cell-line identity and species filtering.

- `fusions`  
  Fusion events and gene-level partner information.

- `metabolomics`  
  Metabolite abundance tables.

- `cosmic_cancergenecensus_v104_grch37`  
  Cancer-gene census data used for gene bridge analysis.

### Generated analytical outputs

- `db/celllineselector.duckdb`  
  Persistent database used for analysis.

- `results/`  
  Statistical results, calibration outputs, chosen parameters, and summary metrics.

- `reports/`  
  EDA and validation outputs such as missingness and distribution summaries.

- `figs/`  
  Publication-style figures and audits.

---

## Data validation and auditing

The project includes a strong validation layer to assess identifier quality and data consistency.

- `tests/test_transcriptomics_gene_filter.py`  
  Tests for ENSG filtering, version normalisation, and protein-coding gene filtering.

- `audit.py`  
  Independent audit for non-human leaks and identifier integrity.

- `harmonisation_validation.py`  
  Checks across gene roster, model coverage, and internal consistency.

- `expression_audit.py`  
  Audits ENSG coverage and expression-table validation.

- `cosmic_cancergene_census_audit.py`  
  Checks for nulls, duplicates, and mapping drift in COSMIC gene records.

- `transcriptomics_z_audit.py`  
  Verifies generated transcriptomics z-score outputs.

- `audit_protein_z.py`  
  Checks protein z-score table structure and source coverage.

These scripts are designed to make the workflow inspectable rather than relying on opaque batch transformations.

---

## Environment and setup

To create the project environment:

```bash
./setup_env.sh
```

This script creates a local virtual environment named `.venv` and installs the dependencies listed in `requirements.txt`.

To activate it manually:

```bash
source .venv/bin/activate
```

---

## Core dependencies

The declared requirements are minimal but central to the project:

- pandas
- duckdb
- pyarrow

The codebase also relies on scientific and plotting libraries during analysis and figure generation, and these may need to be added to the environment explicitly depending on how the project is executed in a fresh installation.

---

## Typical workflow

A standard analytical path is as follows:

1. Run raw data loading
2. Run cleaning
3. Run harmonisation
4. Load results into DuckDB
5. Run EDA and validation
6. Execute transcriptomic or protein scoring stages
7. Review logs, figures, and outputs

This staged design keeps each transformation inspectable and makes debugging easier when a stage fails or produces unexpected data quality patterns.

---

## Notes

- The project is built around reproducibility, data provenance, and auditability.
- Wide expression and omics tables are handled with a combination of Parquet storage, DuckDB access, and chunked analysis to keep processing tractable.
- Several scripts and outputs in the repository show both current and legacy research stages; the canonical workflow should be interpreted in context and validated against the latest pipeline execution.

---

## Summary

Cell Line Selector is a comprehensive data-processing and multi-omics integration project for cell-line selection. It combines quality-controlled raw data, harmonised gene and model identifiers, cross-source analytical scoring, and reproducible validation to prepare a robust resource for biomedical cell-line analysis.
