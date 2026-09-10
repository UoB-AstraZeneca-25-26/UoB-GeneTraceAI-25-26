"""API-key auth and per-IP rate limiting (closes gap G12).

Both are opt-in via Settings. An empty `api_keys` list disables auth
entirely — acceptable for a single-replica dissertation demo behind a
private URL — while rate limiting stays on with a generous default so a
runaway client can't silently burn the whole Bedrock quota.
"""

from __future__ import annotations

from fastapi import Header, HTTPException
from slowapi import Limiter
from slowapi.util import get_remote_address

from .settings import get_settings

limiter = Limiter(key_func=get_remote_address)


async def require_api_key(x_api_key: str | None = Header(default=None, alias="X-API-Key")) -> None:
    settings = get_settings()
    if not settings.api_keys:
        return
    if x_api_key not in settings.api_keys:
        raise HTTPException(status_code=401, detail="Missing or invalid API key.")
