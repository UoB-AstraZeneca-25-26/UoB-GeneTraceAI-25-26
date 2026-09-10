"""Enforces the latency budget from plan.md's Latency Architecture section.

Every test stubs the LLM/executor, so these measure this codebase's own
overhead — never Bedrock's real generation time — and must always pass
offline. A regression here silently destroys the "responds in
milliseconds" guarantee even though the demo still "looks like it works".

The three SSE-timing tests below iterate `_stream_events` directly instead
of going through TestClient/httpx: httpx's ASGITransport collects the
*entire* ASGI response before returning anything (see
httpx._transports.asgi.ASGITransport.handle_async_request), so it cannot
measure real time-to-first-byte — it would report ~5000ms for every SSE
response regardless of how fast this codebase actually flushes events.
Iterating the generator directly tests exactly the logic this design is
about, without that transport-layer artifact.
"""

import asyncio
import json
import time
from types import SimpleNamespace

import pytest
from fastapi.testclient import TestClient

import api.routes as routes_module
import http_utils
from api.app import app
from api.schemas import QueryRequest


def _tool_end_event(output: dict) -> dict:
    return {
        "event": "on_tool_end",
        "name": "gene_alias_lookup",
        "data": {"input": {"query": "TP53"}, "output": json.dumps(output)},
    }


def _token_event(text: str) -> dict:
    return {"event": "on_chat_model_stream", "data": {"chunk": SimpleNamespace(content=text)}}


class SlowFakeExecutor:
    """Emits the tool result immediately, then blocks 5s before any token —
    models a slow Bedrock call. Tests below only read the events they need and
    close the generator before the sleep would ever resolve."""

    async def astream_events(self, payload, version="v2"):
        yield _tool_end_event({"found": True, "symbol": "TP53"})
        await asyncio.sleep(5.0)
        yield _token_event("TP53 is p53.")


class HappyExecutor:
    async def astream_events(self, payload, version="v2"):
        yield _tool_end_event({"found": True, "symbol": "TP53"})
        yield _token_event("TP53 is p53.")


class FailingAfterToolExecutor:
    async def astream_events(self, payload, version="v2"):
        yield _tool_end_event({"found": True, "symbol": "TP53"})
        raise RuntimeError("bedrock is down")


def _fake_request() -> SimpleNamespace:
    return SimpleNamespace(headers={}, state=SimpleNamespace(request_id="test-request-id"))


async def _collect(gen, n):
    events = []
    async for evt in gen:
        events.append(evt)
        if len(events) >= n:
            break
    return events


@pytest.mark.asyncio
async def test_time_to_first_byte_under_50ms(monkeypatch):
    monkeypatch.setattr(routes_module, "executor", SlowFakeExecutor())
    gen = routes_module._stream_events(QueryRequest(query="Tell me about TP53"), _fake_request())

    start = time.perf_counter()
    first = await gen.__anext__()
    elapsed_ms = (time.perf_counter() - start) * 1000
    await gen.aclose()

    assert first["event"] == "accepted"
    assert elapsed_ms < 50, f"first byte took {elapsed_ms:.1f}ms"


@pytest.mark.asyncio
async def test_data_arrives_before_any_token(monkeypatch):
    monkeypatch.setattr(routes_module, "executor", SlowFakeExecutor())
    gen = routes_module._stream_events(QueryRequest(query="Tell me about TP53"), _fake_request())

    events = await _collect(gen, 3)
    await gen.aclose()

    # The narration stub sleeps 5s before any token, so if a "token" ever
    # made it into this 3-event window, ordering would already be broken.
    assert [e["event"] for e in events] == ["accepted", "data", "status"]


@pytest.mark.asyncio
async def test_time_to_data_event_under_100ms(monkeypatch):
    monkeypatch.setattr(routes_module, "executor", SlowFakeExecutor())
    gen = routes_module._stream_events(QueryRequest(query="Tell me about TP53"), _fake_request())

    start = time.perf_counter()
    events = await _collect(gen, 2)  # accepted, data
    elapsed_ms = (time.perf_counter() - start) * 1000
    await gen.aclose()

    assert events[1]["event"] == "data"
    assert elapsed_ms < 100, f"data event took {elapsed_ms:.1f}ms"


def test_sse_ping_interval_matches_heartbeat_budget():
    # A live 25s-gap heartbeat test would make the suite slow on every CI
    # run; what actually regresses is the configured ping interval, so
    # assert that directly instead of waiting out a real 10s heartbeat.
    result = routes_module._stream(QueryRequest(query="Tell me about TP53"), _fake_request())

    assert result.ping_interval == 10  # must stay under nginx's 60s proxy_read_timeout


def test_llm_failure_still_returns_data(monkeypatch):
    monkeypatch.setattr(routes_module, "executor", FailingAfterToolExecutor())
    routes_module._response_cache.clear()

    with TestClient(app) as c:
        response = c.post("/v1/agent/query", json={"query": "Tell me about TP53"})

    assert response.status_code == 200
    body = response.json()
    assert body["answer"] is None
    assert body["data"] == {"found": True, "symbol": "TP53"}
    assert body["degraded"] == ["llm_unavailable"]


def test_cache_hit_under_10ms_and_marked_cached(monkeypatch):
    monkeypatch.setattr(routes_module, "executor", HappyExecutor())
    routes_module._response_cache.clear()

    with TestClient(app) as c:
        c.post("/v1/agent/query", json={"query": "Tell me about TP53 cache warm"})

        start = time.perf_counter()
        response = c.post("/v1/agent/query", json={"query": "Tell me about TP53 cache warm"})
        elapsed_ms = (time.perf_counter() - start) * 1000

    assert response.json()["cached"] is True
    assert elapsed_ms < 10, f"cache hit took {elapsed_ms:.1f}ms"


def test_health_is_fast_and_makes_zero_outbound_calls(monkeypatch):
    def _boom(*args, **kwargs):
        raise AssertionError("no HTTP client should be touched for /health")

    monkeypatch.setattr(http_utils, "get_http_client", _boom)

    with TestClient(app) as c:
        start = time.perf_counter()
        response = c.get("/health")
        elapsed_ms = (time.perf_counter() - start) * 1000

    assert response.status_code == 200
    assert elapsed_ms < 50, f"/health took {elapsed_ms:.1f}ms"


@pytest.mark.asyncio
async def test_fast_path_triggers_exactly_one_llm_call(monkeypatch):
    call_count = {"n": 0}

    class FakeLLM:
        async def astream(self, messages):
            call_count["n"] += 1
            for text in ("TP53 ", "is p53."):
                yield SimpleNamespace(content=text)

    def _executor_must_not_be_called(*args, **kwargs):
        raise AssertionError("fast path must not invoke the full AgentExecutor")

    monkeypatch.setattr(routes_module, "_llm", FakeLLM())
    monkeypatch.setattr(routes_module, "executor", SimpleNamespace(astream_events=_executor_must_not_be_called))

    events = [evt async for evt in routes_module._agent_events("ENSG00000141510")]

    assert call_count["n"] == 1
    data_events = [e for e in events if e.kind == "data"]
    assert data_events and data_events[0].payload["tool"] == "gene_alias_lookup"
