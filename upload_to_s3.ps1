# GeneTraceAI — S3 Upload Script
# Run from the project root: .\upload_to_s3.ps1
# Prerequisites: aws configure (with your IAM key + secret)

$BUCKET  = "genetraceai-data"
$REGION  = "eu-north-1"
$ROOT    = "C:\Disertation\UoB-GeneTraceAI-25-26"

function S3Sync($local, $s3key) {
    Write-Host "`n>>> Syncing: $local  -->  s3://$BUCKET/$s3key" -ForegroundColor Cyan
    aws s3 sync $local "s3://$BUCKET/$s3key" --region $REGION
}

# ─── RAW DATA ──────────────────────────────────────────────────────────────────

# DepMap 24Q4  (~900 MB: CRISPRGeneEffect, OmicsExpression, Model, Achilles controls)
S3Sync "$ROOT\data\DepMap_24Q4"         "raw/DepMap_24Q4"

# DEMETER2 RNAi  (~490 MB: D2_combined + Achilles + DRIVE scores)
S3Sync "$ROOT\data\DEMETER2"            "raw/DEMETER2"

# COSMIC v104  (~16 MB active files; ~55 READMEs)
S3Sync "$ROOT\data\COSMIC"             "raw/COSMIC"

# GDSC  (~330 MB: Project Score TSVs, dose-response XLSX, model list)
S3Sync "$ROOT\data\GDSC"               "raw/GDSC"

# ProCan proteomics  (~49 MB: Protein_matrix_averaged)
S3Sync "$ROOT\data\proteomics_procan"  "raw/proteomics_procan"

# DepMap Chronos  (~340 MB: HDF5 copy-number + fitness effects)
S3Sync "$ROOT\data\DepMap_Chronos"     "raw/DepMap_Chronos"

# Ensembl paralogy  (~21 MB parquet + PROVENANCE.json)
S3Sync "$ROOT\data\ensembl_paralogy"   "raw/ensembl_paralogy"

# Gene expression raw CSVs  (~3 GB — largest chunk, takes a while)
# Includes: HPA 1 GB, DepMap expression 866 MB, GEO 1 GB, mutations 479 MB
S3Sync "$ROOT\data\nomenclature"       "raw/nomenclature"

# Parquet numbered inputs (raw_data/)  — the originals the pipeline was built on
S3Sync "$ROOT\data\parquet\raw_data"   "raw/parquet/raw"

# Parquet cleaned versions (data_clean/)
S3Sync "$ROOT\data\parquet\data_clean" "raw/parquet/clean"

# ─── REFERENCE LOOKUP TABLES ───────────────────────────────────────────────────
# cell_line_lookup, gene_lookup, union_lookup_entity, gene_symbol_resolver.json
# depmap_profiles, protein_map
S3Sync "$ROOT\reference"               "reference"

# ─── GATE AUDIT SCRIPTS ────────────────────────────────────────────────────────
# 00_build_panel.py through 09_*, common.py, run_all.py, README, outputs/
S3Sync "$ROOT\gate_audit"              "gate_audit"

# ─── PIPELINE OUTPUTS (current on-disk state) ──────────────────────────────────
# Uploads whatever is already in outputs/ so teammates start with a baseline.
# The pipeline will overwrite these when it runs on EC2.
S3Sync "$ROOT\src\pipeline\outputs"    "pipeline-outputs"

Write-Host "`n=== Upload complete ===" -ForegroundColor Green
Write-Host "Bucket: s3://$BUCKET  (region: $REGION)" -ForegroundColor Green
