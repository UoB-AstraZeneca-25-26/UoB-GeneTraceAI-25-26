"""
architecture/api/app.py
-------------------------
FastAPI application for the GeneTraceAI ranking endpoint.

Deployed as a Lambda container image via handler.py (Mangum adapter).
Run locally with:  uvicorn architecture.api.app:app --reload
"""

from __future__ import annotations

import logging
from contextlib import asynccontextmanager

from fastapi import FastAPI
from fastapi.middleware.cors import CORSMiddleware

from .ranking import ensure_loaded
from .routes import router

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger("genetraceai.api")


@asynccontextmanager
async def lifespan(app: FastAPI):
    """
    Application lifespan: warm the ranking data once per cold start.

    Calls :func:`ranking.ensure_loaded` before the app starts accepting
    requests, so the httpx client, gene index, and any lazily-fetched
    S3 data are already loaded by the time the first request arrives.
    On Lambda this cost is paid once per container (cold start), not
    once per request — every subsequent warm invocation reuses the
    same loaded state.

    Parameters
    ----------
    app : fastapi.FastAPI
        The application instance being started. Unused directly, but
        required by FastAPI's lifespan protocol.

    Returns
    -------
    None
        This is an async generator used as a context manager: code
        before ``yield`` runs at startup, code after would run at
        shutdown (none here).
    """
    ensure_loaded()
    yield


def create_app() -> FastAPI:
    """
    Build and configure the FastAPI application.

    Wires up the lifespan hook (see :func:`lifespan`), adds a
    permissive CORS policy (GET/OPTIONS only, no cookies, so an
    allow-all origin list is safe), and mounts the ranking router.

    Parameters
    ----------
    None

    Returns
    -------
    fastapi.FastAPI
        Configured but not yet running application instance. The
        module-level ``app`` object below is this function's result —
        import that, not this function, when wiring up a server.
    """
    app = FastAPI(
        title="GeneTraceAI Ranking API",
        version="1.0.0",
        lifespan=lifespan,
    )
    # The hosted UI runs on a different origin than this API, so browsers require
    # CORS headers on every response or they block the call. GET-only, no cookies,
    # so a permissive allow-list is fine.
    app.add_middleware(
        CORSMiddleware,
        allow_origins=["*"],
        allow_methods=["GET", "OPTIONS"],
        allow_headers=["*"],
    )
    app.include_router(router)
    return app


app = create_app()
