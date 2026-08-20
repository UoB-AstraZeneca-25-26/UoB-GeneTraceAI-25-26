"""
final_pipeline/api/app.py
--------------------------
FastAPI application for the GeneTraceAI ranking endpoint.

Deployed as a Lambda container image via handler.py (Mangum adapter).
Run locally with:  uvicorn final_pipeline.api.app:app --reload
"""

from __future__ import annotations

import logging
from contextlib import asynccontextmanager

from fastapi import FastAPI

from .ranking import ensure_loaded
from .routes import router

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger("genetraceai.api")


@asynccontextmanager
async def lifespan(app: FastAPI):
    ensure_loaded()
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
