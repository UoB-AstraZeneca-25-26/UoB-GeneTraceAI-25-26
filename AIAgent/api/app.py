"""FastAPI app factory + lifespan.

Builds everything that must exist exactly once per process — the shared
httpx client, the local gene index, a warmed LLM connection — and tears it
down cleanly on shutdown. The service is stateless otherwise (H2): any
replica can serve any request.
"""

from __future__ import annotations

import logging
from contextlib import asynccontextmanager

import httpx
from fastapi import FastAPI
from slowapi.errors import RateLimitExceeded
from slowapi.middleware import SlowAPIMiddleware

from . import _pathsetup  # noqa: F401  side effect: sys.path for the flat modules below

from gene_index import load_gene_index  # noqa: E402
import FunctionCalling  # noqa: E402
import http_utils  # noqa: E402

from .errors import error_response, install_exception_handlers
from .middleware import BodySizeLimitMiddleware, RequestIDMiddleware, TimingMiddleware, install_cors
from .routes import router
from .security import limiter
from .settings import get_settings

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger("celllinefinder.startup")


@asynccontextmanager
async def lifespan(app: FastAPI):
    settings = get_settings()

    client = httpx.AsyncClient(
        timeout=httpx.Timeout(settings.http_timeout, connect=settings.http_connect_timeout),
        limits=httpx.Limits(max_keepalive_connections=20, max_connections=100),
        headers={"User-Agent": "CellLineFinder/1.0 (dissertation; UoB)"},
    )
    http_utils.set_http_client(client)
    app.state.http = client

    index = load_gene_index(settings.gene_lookup_path)
    FunctionCalling.set_gene_index(index)
    app.state.gene_index = index
    logger.info("gene index loaded: %d genes from %s", len(index), settings.gene_lookup_path)

    # Prewarm: pay the Groq TLS handshake now instead of on the first real
    # request. Best-effort — an upstream blip at startup must not crash the
    # container (the health check must stay reachable regardless).
    try:
        from BusinessFlow import llm as _llm

        await _llm.ainvoke("ping")
    except Exception as exc:  # noqa: BLE001
        logger.warning("LLM prewarm failed (continuing): %s", exc)

    yield

    await http_utils.close_http_client()


def create_app() -> FastAPI:
    settings = get_settings()

    app = FastAPI(title="CellLineFinder Agent", version="1.0.0", lifespan=lifespan)

    app.state.limiter = limiter
    app.add_middleware(SlowAPIMiddleware)

    async def _rate_limit_handler(request, exc: RateLimitExceeded):
        response = error_response(
            request,
            type_="client_rate_limited",
            title="Too Many Requests",
            status=429,
            detail="Rate limit exceeded. Slow down and retry shortly.",
        )
        response.headers["Retry-After"] = "60"
        return response

    app.add_exception_handler(RateLimitExceeded, _rate_limit_handler)
    install_exception_handlers(app)

    install_cors(app, settings.allowed_origins)
    app.add_middleware(BodySizeLimitMiddleware, max_bytes=settings.max_body_bytes)
    app.add_middleware(TimingMiddleware)
    app.add_middleware(RequestIDMiddleware)

    app.include_router(router)

    return app


app = create_app()
