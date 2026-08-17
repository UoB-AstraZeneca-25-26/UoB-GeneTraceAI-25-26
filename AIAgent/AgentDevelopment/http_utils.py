"""Shared outbound HTTP plumbing for the CellLineFinder tools.

Provides:
- a process-wide pooled ``httpx.AsyncClient`` (fixes G14 — no more TLS
  handshake per call),
- ``safe_get``: retry with exponential backoff + jitter, never raises
  (fixes G2),
- ``CircuitBreaker``: stops hammering a dead upstream (R3).
"""

from __future__ import annotations

import asyncio
import random
import time

import httpx

_RETRYABLE = (httpx.TimeoutException, httpx.ConnectError, httpx.ReadError)
_RETRYABLE_STATUS = {429, 500, 502, 503, 504}

_client: httpx.AsyncClient | None = None


def get_http_client() -> httpx.AsyncClient:
    """Return the process-wide pooled client, creating it on first use."""
    global _client
    if _client is None or _client.is_closed:
        _client = httpx.AsyncClient(
            timeout=httpx.Timeout(10.0, connect=3.0),
            limits=httpx.Limits(max_keepalive_connections=20, max_connections=100),
            headers={"User-Agent": "CellLineFinder/1.0 (dissertation; UoB)"},
        )
    return _client


def set_http_client(client: httpx.AsyncClient) -> None:
    """Install an externally-owned client (used by the FastAPI lifespan)."""
    global _client
    _client = client


async def close_http_client() -> None:
    global _client
    if _client is not None and not _client.is_closed:
        await _client.aclose()
    _client = None


async def safe_get(
    url: str,
    *,
    headers: dict[str, str] | None = None,
    client: httpx.AsyncClient | None = None,
    attempts: int = 3,
) -> dict | None:
    """GET with exponential backoff + jitter. Returns ``None`` on any
    failure instead of raising, so a dead upstream degrades gracefully.
    """
    http = client or get_http_client()
    for attempt in range(attempts):
        try:
            response = await http.get(url, headers=headers)
        except _RETRYABLE:
            if attempt < attempts - 1:
                await asyncio.sleep(0.25 * 2**attempt + random.uniform(0, 0.1))
                continue
            return None

        if response.status_code == 200:
            return response.json()
        if response.status_code in _RETRYABLE_STATUS and attempt < attempts - 1:
            await asyncio.sleep(0.25 * 2**attempt + random.uniform(0, 0.1))
            continue
        return None
    return None


class CircuitBreaker:
    """Trips after N consecutive failures; skips the upstream for a cooldown
    window; allows one probe request when the window elapses (half-open).
    """

    def __init__(self, failure_threshold: int = 5, reset_after: float = 60.0):
        self.failure_threshold = failure_threshold
        self.reset_after = reset_after
        self._failures = 0
        self._opened_at: float | None = None

    def is_open(self) -> bool:
        if self._opened_at is None:
            return False
        if time.monotonic() - self._opened_at >= self.reset_after:
            # half-open: allow exactly one probe through
            self._opened_at = None
            self._failures = self.failure_threshold - 1
            return False
        return True

    def record_success(self) -> None:
        self._failures = 0
        self._opened_at = None

    def record_failure(self) -> None:
        self._failures += 1
        if self._failures >= self.failure_threshold and self._opened_at is None:
            self._opened_at = time.monotonic()
