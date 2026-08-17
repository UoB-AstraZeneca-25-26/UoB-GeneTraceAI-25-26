#!/usr/bin/env bash
# GeneTraceAI — Mac / Linux Setup Script
# Run from the repo root: bash setup.sh
set -e

BUCKET="https://genetraceai-data.s3.eu-north-1.amazonaws.com"

download() {
    local url=$1 dest=$2
    mkdir -p "$(dirname "$dest")"
    if [ -f "$dest" ]; then
        echo "  already exists, skipping: $dest"
        return
    fi
    echo "  Downloading $dest ..."
    curl -L --progress-bar "$url" -o "$dest"
}

# ── 1. Python packages ────────────────────────────────────────────────────────
echo ""
echo "[1/3] Installing Python packages..."
pip install pandas pyarrow duckdb --quiet
echo "      Done."

# ── 2. Data files from S3 ────────────────────────────────────────────────────
echo ""
echo "[2/3] Downloading data files (~2.3 GB)..."

download "$BUCKET/pipeline-outputs/predictions_with_confidence.parquet" \
         "final_pipeline/outputs/predictions_with_confidence.parquet"

download "$BUCKET/pipeline-outputs/celllineselector.db" \
         "final_pipeline/outputs/celllineselector.db"

download "$BUCKET/final-reference/gene_lookup.parquet" \
         "final_pipeline/reference/gene_lookup.parquet"

download "$BUCKET/final-reference/cell_line_lookup.parquet" \
         "final_pipeline/reference/cell_line_lookup.parquet"

download "$BUCKET/cell-similarity/rna_neighbours.parquet" \
         "cell_similarity/outputs/rna_neighbours.parquet"

# ── 3. Smoke test ─────────────────────────────────────────────────────────────
echo ""
echo "[3/3] Running smoke test (gene BRAF)..."
python final_pipeline/Ranking/cli.py gene BRAF

echo ""
echo "=== Setup complete. See SETUP.md for usage. ==="
