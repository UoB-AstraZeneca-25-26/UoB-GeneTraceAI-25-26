"""
final_pipeline/api/app.py
--------------------------
FastAPI application for the GeneTraceAI ranking endpoint.

Deployed as a Lambda container image via handler.py (Mangum adapter).
Run locally with:  uvicorn final_pipeline.api.app:app --reload
"""

from __future__ import annotations

import logging
import os
from contextlib import asynccontextmanager
from pathlib import Path

from fastapi import FastAPI

from .ranking import load_ranking_data
from .routes import router

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger("genetraceai.api")

_PIPELINE = Path(__file__).resolve().parent.parent
_REPO = _PIPELINE.parent


def _path(env_var: str, default: Path) -> Path:
    """Return env var as Path if set, else the default. Allows Lambda overrides."""
    val = os.getenv(env_var)
    return Path(val) if val else default


@asynccontextmanager
async def lifespan(app: FastAPI):
    load_ranking_data(
        predictions_path=_path("PREDICTIONS_PATH",
                               _PIPELINE / "outputs" / "predictions_with_confidence.parquet"),
        gene_lookup_path=_path("GENE_LOOKUP_PATH",
                               _PIPELINE / "reference" / "gene_lookup.parquet"),
        cell_lookup_path=_path("CELL_LOOKUP_PATH",
                               _PIPELINE / "reference" / "cell_line_lookup.parquet"),
        db_path=_path("RANKING_DB_PATH",
                      _PIPELINE / "outputs" / "celllineselector.db"),
    )
    yield


def create_app() -> FastAPI:
    app = FastAPI(
        title="GeneTraceAI Ranking API",
        version="1.0.0",
        lifespan=lifespan,
    )
    app.include_router(router)
    return app


app = create_app()
