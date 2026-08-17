"""Cross-cutting middleware: request ID propagation, timing/structured logs,
CORS, and a body-size cap.
"""

from __future__ import annotations

import time
import uuid

import structlog
from fastapi import FastAPI
from starlette.middleware.base import BaseHTTPMiddleware
from starlette.middleware.cors import CORSMiddleware
from starlette.requests import Request
from starlette.responses import JSONResponse, Response

logger = structlog.get_logger("celllinefinder.api")


class RequestIDMiddleware(BaseHTTPMiddleware):
    """Accepts an inbound X-Request-ID or generates a UUID4; echoes it back
    and attaches it to `request.state` so routes and exception handlers can
    tag every log line and error body with it."""

    async def dispatch(self, request: Request, call_next):
        request_id = request.headers.get("X-Request-ID") or str(uuid.uuid4())
        request.state.request_id = request_id
        response: Response = await call_next(request)
        response.headers["X-Request-ID"] = request_id
        return response


class TimingMiddleware(BaseHTTPMiddleware):
    """Records duration and emits one structured JSON log line per request."""

    async def dispatch(self, request: Request, call_next):
        start = time.perf_counter()
        response: Response = await call_next(request)
        duration_ms = int((time.perf_counter() - start) * 1000)
        response.headers["X-Response-Time-Ms"] = str(duration_ms)
        logger.info(
            "request",
            request_id=getattr(request.state, "request_id", None),
            method=request.method,
            path=request.url.path,
            status_code=response.status_code,
            latency_ms=duration_ms,
        )
        return response


class BodySizeLimitMiddleware(BaseHTTPMiddleware):
    """Rejects payloads over `max_bytes` before they reach a route (G-level
    hardening: an unbounded body is a free way to burn LLM context/cost)."""

    def __init__(self, app, max_bytes: int):
        super().__init__(app)
        self.max_bytes = max_bytes

    async def dispatch(self, request: Request, call_next):
        content_length = request.headers.get("content-length")
        if content_length is not None and int(content_length) > self.max_bytes:
            return JSONResponse(
                status_code=413,
                content={
                    "type": "payload_too_large",
                    "title": "Payload Too Large",
                    "status": 413,
                    "detail": f"Body exceeds the {self.max_bytes}-byte limit.",
                    "request_id": getattr(request.state, "request_id", "unknown"),
                },
            )
        return await call_next(request)


def install_cors(app: FastAPI, allowed_origins: list[str]) -> None:
    # Never pair allow_origins=["*"] with allow_credentials=True — this
    # stays an explicit allowlist and credential-free by design.
    app.add_middleware(
        CORSMiddleware,
        allow_origins=allowed_origins,
        allow_credentials=False,
        allow_methods=["GET", "POST", "OPTIONS"],
        allow_headers=["*"],
        expose_headers=["X-Request-ID", "X-Response-Time-Ms"],
    )
