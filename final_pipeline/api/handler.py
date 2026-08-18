"""
api/handler.py
--------------
AWS Lambda entry point. Routes API Gateway requests to query functions.

Endpoints:
  GET /gene?query=BRAF
  GET /gene?query=BRAF&cell_line=ACH-000219
  GET /genes?query=BRAF,KRAS
  GET /exclude?gene_a=BRAF&gene_b=EGFR
"""
from __future__ import annotations
import json
import logging
import os
from pathlib import Path

import boto3

logger = logging.getLogger()
logger.setLevel(logging.INFO)

# ── S3 → /tmp file cache (runs once per cold start) ───────────────────────────
BUCKET   = os.environ.get("DATA_BUCKET", "genetraceai-data")
TMP      = Path("/tmp/genetraceai")

_S3_FILES = {
    "predictions_with_confidence.parquet": "pipeline-outputs/predictions_with_confidence.parquet",
    "gene_lookup.parquet":                 "final-reference/gene_lookup.parquet",
    "cell_line_lookup.parquet":            "final-reference/cell_line_lookup.parquet",
    "rna_neighbours.parquet":              "cell-similarity/rna_neighbours.parquet",
    # celllineselector.db (~2 GB) omitted — too large for Lambda memory.
    # Lineage/sex metadata fields return empty gracefully when DB is absent.
}

_files_ready = False


def _ensure_files():
    global _files_ready
    if _files_ready:
        return
    TMP.mkdir(parents=True, exist_ok=True)
    s3 = boto3.client("s3")
    for filename, s3_key in _S3_FILES.items():
        dest = TMP / filename
        if dest.exists():
            continue
        logger.info(f"Downloading s3://{BUCKET}/{s3_key} → {dest}")
        s3.download_file(BUCKET, s3_key, str(dest))
        logger.info(f"Downloaded {filename} ({dest.stat().st_size / 1e6:.1f} MB)")

    # Point queries.py at /tmp
    os.environ["DATA_DIR"] = str(TMP)
    os.environ["REF_DIR"]  = str(TMP)
    os.environ["SIM_DIR"]  = str(TMP)
    _files_ready = True


def _response(status: int, body: dict) -> dict:
    return {
        "statusCode": status,
        "headers": {
            "Content-Type": "application/json",
            "Access-Control-Allow-Origin": "*",
        },
        "body": json.dumps(body, default=str),
    }


def _params(event: dict) -> dict:
    return event.get("queryStringParameters") or {}


def handler(event, context):
    try:
        _ensure_files()
        from queries import query_gene, query_genes, query_exclude
    except Exception as e:
        logger.exception("Startup failed")
        return _response(500, {"error": f"Startup error: {e}"})

    path   = event.get("rawPath") or event.get("path", "/")
    params = _params(event)
    route  = path.rstrip("/").split("/")[-1]   # last segment: gene | genes | exclude

    try:
        if route == "gene":
            gene = params.get("query") or params.get("gene")
            if not gene:
                return _response(400, {"error": "Missing ?query=<GENE>"})
            cell_line = params.get("cell_line")
            result = query_gene(gene, cell_line)

        elif route == "genes":
            raw = params.get("query") or params.get("genes", "")
            gene_list = [g.strip() for g in raw.split(",") if g.strip()]
            if len(gene_list) < 2:
                return _response(400, {"error": "Provide at least two genes: ?query=BRAF,KRAS"})
            result = query_genes(gene_list)

        elif route == "exclude":
            gene_a = params.get("gene_a")
            gene_b = params.get("gene_b")
            if not gene_a or not gene_b:
                return _response(400, {"error": "Missing ?gene_a=<GENE>&gene_b=<GENE>"})
            result = query_exclude(gene_a, gene_b)

        else:
            return _response(404, {
                "error": f"Unknown route: {path}",
                "valid_routes": ["/gene", "/genes", "/exclude"],
            })

        if "error" in result:
            return _response(404, result)
        return _response(200, result)

    except ValueError as e:
        return _response(400, {"error": str(e)})
    except Exception as e:
        logger.exception("Query failed")
        return _response(500, {"error": str(e)})
