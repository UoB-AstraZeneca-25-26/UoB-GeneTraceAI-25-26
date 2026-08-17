"""Typed exceptions + a single error envelope (RFC 9457-flavoured).

Every error response has the same shape, so a caller only ever needs one
parser:

    {type, title, status, detail, request_id}

Stack traces and secrets never leave the process — they're logged
server-side against `request_id`; only the id crosses the wire.
"""

from __future__ import annotations

import logging

from fastapi import FastAPI, Request
from fastapi.exceptions import RequestValidationError
from fastapi.responses import JSONResponse
from pydantic import BaseModel
from starlette.exceptions import HTTPException as StarletteHTTPException

logger = logging.getLogger("celllinefinder.errors")

_STATUS_TYPES = {
    401: "unauthorized",
    403: "forbidden",
    404: "not_found",
    413: "payload_too_large",
    429: "client_rate_limited",
}


class ErrorResponse(BaseModel):
    type: str
    title: str
    status: int
    detail: str
    request_id: str


class UpstreamRateLimitedError(Exception):
    """Raised when the LLM provider (Groq) 429s or its quota is exhausted."""

    def __init__(
        self,
        detail: str = "The LLM provider's rate limit was hit.",
        retry_after: int = 5,
    ):
        super().__init__(detail)
        self.detail = detail
        self.retry_after = retry_after


class AgentTimeoutError(Exception):
    """Raised when the agent exceeds max_iterations / max_execution_time."""

    def __init__(self, detail: str = "The agent exceeded its execution budget."):
        super().__init__(detail)
        self.detail = detail


def request_id_of(request: Request) -> str:
    return getattr(request.state, "request_id", "unknown")


def error_response(
    request: Request, *, type_: str, title: str, status: int, detail: str
) -> JSONResponse:
    body = ErrorResponse(
        type=type_, title=title, status=status, detail=detail, request_id=request_id_of(request)
    )
    return JSONResponse(status_code=status, content=body.model_dump())


def install_exception_handlers(app: FastAPI) -> None:
    @app.exception_handler(RequestValidationError)
    async def _validation_error(request: Request, exc: RequestValidationError) -> JSONResponse:
        return error_response(
            request,
            type_="validation_error",
            title="Validation Error",
            status=422,
            detail=str(exc.errors()),
        )

    @app.exception_handler(UpstreamRateLimitedError)
    async def _rate_limited(request: Request, exc: UpstreamRateLimitedError) -> JSONResponse:
        response = error_response(
            request, type_="rate_limited", title="Upstream Rate Limited", status=503, detail=exc.detail
        )
        response.headers["Retry-After"] = str(exc.retry_after)
        return response

    @app.exception_handler(AgentTimeoutError)
    async def _agent_timeout(request: Request, exc: AgentTimeoutError) -> JSONResponse:
        return error_response(
            request, type_="agent_timeout", title="Agent Timeout", status=504, detail=exc.detail
        )

    @app.exception_handler(StarletteHTTPException)
    async def _http_exception(request: Request, exc: StarletteHTTPException) -> JSONResponse:
        return error_response(
            request,
            type_=_STATUS_TYPES.get(exc.status_code, "http_error"),
            title=exc.detail if isinstance(exc.detail, str) else "HTTP Error",
            status=exc.status_code,
            detail=str(exc.detail),
        )

    @app.exception_handler(Exception)
    async def _unhandled(request: Request, exc: Exception) -> JSONResponse:
        logger.exception("unhandled_error request_id=%s", request_id_of(request))
        return error_response(
            request,
            type_="internal_error",
            title="Internal Server Error",
            status=500,
            detail="An unexpected error occurred. Reference the request_id when reporting this.",
        )
