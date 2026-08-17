# GeneTraceAI — Windows Setup Script
# Run from the repo root: .\setup.ps1

$ErrorActionPreference = "Stop"
$BUCKET = "https://genetraceai-data.s3.eu-north-1.amazonaws.com"

function Download($url, $dest) {
    $dir = Split-Path $dest
    if ($dir) { New-Item -ItemType Directory -Force $dir | Out-Null }
    if (Test-Path $dest) {
        Write-Host "  already exists, skipping: $dest" -ForegroundColor DarkGray
        return
    }
    Write-Host "  Downloading $dest ..." -ForegroundColor Cyan
    Invoke-WebRequest -Uri $url -OutFile $dest -UseBasicParsing
}

# ── 1. Python packages ────────────────────────────────────────────────────────
Write-Host "`n[1/3] Installing Python packages..." -ForegroundColor Yellow
pip install pandas pyarrow duckdb --quiet
Write-Host "      Done." -ForegroundColor Green

# ── 2. Data files from S3 ────────────────────────────────────────────────────
Write-Host "`n[2/3] Downloading data files (~2.3 GB)..." -ForegroundColor Yellow

Download "$BUCKET/pipeline-outputs/predictions_with_confidence.parquet" `
         "final_pipeline\outputs\predictions_with_confidence.parquet"

Download "$BUCKET/pipeline-outputs/celllineselector.db" `
         "final_pipeline\outputs\celllineselector.db"

Download "$BUCKET/final-reference/gene_lookup.parquet" `
         "final_pipeline\reference\gene_lookup.parquet"

Download "$BUCKET/final-reference/cell_line_lookup.parquet" `
         "final_pipeline\reference\cell_line_lookup.parquet"

Download "$BUCKET/cell-similarity/rna_neighbours.parquet" `
         "cell_similarity\outputs\rna_neighbours.parquet"

# ── 3. Smoke test ─────────────────────────────────────────────────────────────
Write-Host "`n[3/3] Running smoke test (gene BRAF)..." -ForegroundColor Yellow
python final_pipeline/Ranking/cli.py gene BRAF

Write-Host "`n=== Setup complete. See SETUP.md for usage. ===" -ForegroundColor Green
