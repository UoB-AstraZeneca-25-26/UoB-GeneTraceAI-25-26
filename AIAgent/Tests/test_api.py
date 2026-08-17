"""HTTP-layer tests (closes gap G13). The agent executor is mocked, so
these run offline and never touch Groq/Ensembl/HGNC.
"""

import json
from types import SimpleNamespace

import httpx
import pytest
from fastapi.testclient import TestClient
from groq import RateLimitError

import api.routes as routes_module
import http_utils
from api.app import app
from api.schemas import QueryResponse


class FakeExecutor:
    """Stands in for BusinessFlow.executor: replays a canned event sequence
    from astream_events, optionally raising at the end."""

    def __init__(self, events=None, exc=None):
        self._events = events or []
        self._exc = exc

    async def astream_events(self, payload, version="v2"):
        for event in self._events:
            yield event
        if self._exc is not None:
            raise self._exc


def _tool_end_event(name: str, tool_input: dict, output: dict) -> dict:
    return {"event": "on_tool_end", "name": name, "data": {"input": tool_input, "output": json.dumps(output)}}


def _token_event(text: str) -> dict:
    return {"event": "on_chat_model_stream", "data": {"chunk": SimpleNamespace(content=text)}}


GENE_OUTPUT = {
    "found": True, "query": "TP53", "symbol": "TP53", "ensembl_id": "ENSG00000141510",
    "full_name": "tumor protein p53", "previous_symbols": [], "synonyms": ["p53"],
    "cross_validated": True, "sources": ["local_index"], "degraded": [],
}

HAPPY_EVENTS = [
    _tool_end_event("gene_alias_lookup", {"query": "TP53"}, GENE_OUTPUT),
    _token_event("TP53 "),
    _token_event("is also known as p53."),
]


@pytest.fixture()
def client(monkeypatch):
    monkeypatch.setattr(routes_module, "executor", FakeExecutor(events=HAPPY_EVENTS))
    routes_module._response_cache.clear()
    with TestClient(app) as c:
        yield c


def test_valid_query_returns_200_matching_schema(client):
    response = client.post("/v1/agent/query", json={"query": "Tell me about TP53"})

    assert response.status_code == 200
    body = QueryResponse.model_validate(response.json())
    assert body.data == GENE_OUTPUT
    assert body.answer == "TP53 is also known as p53."
    assert body.tool_calls[0].tool == "gene_alias_lookup"


def test_empty_body_returns_422(client):
    response = client.post("/v1/agent/query", json={})

    assert response.status_code == 422
    assert response.json()["type"] == "validation_error"


def test_oversized_query_returns_422(client):
    response = client.post("/v1/agent/query", json={"query": "x" * 5000})

    assert response.status_code == 422


def test_health_returns_200_without_touching_network(monkeypatch):
    def _boom(*args, **kwargs):
        raise AssertionError("health must not open an HTTP client")

    monkeypatch.setattr(http_utils, "get_http_client", _boom)
    with TestClient(app) as c:
        response = c.get("/health")

    assert response.status_code == 200
    assert response.json() == {"status": "ok"}


def test_groq_rate_limit_returns_503_with_retry_after(monkeypatch):
    exc = RateLimitError(
        "rate limited",
        response=httpx.Response(429, request=httpx.Request("POST", "http://x")),
        body=None,
    )
    monkeypatch.setattr(routes_module, "executor", FakeExecutor(events=[], exc=exc))
    routes_module._response_cache.clear()

    with TestClient(app) as c:
        response = c.post("/v1/agent/query", json={"query": "Tell me about some rare gene"})

    assert response.status_code == 503
    assert response.json()["type"] == "rate_limited"
    assert "Retry-After" in response.headers


def test_tool_exception_returns_500_with_no_stack_trace(monkeypatch):
    monkeypatch.setattr(
        routes_module, "executor", FakeExecutor(events=[], exc=ValueError("boom: secret internal detail"))
    )
    routes_module._response_cache.clear()

    # raise_server_exceptions=False: Starlette's TestClient re-raises the
    # original exception to the test process for debugging by default, even
    # though the app already sent a proper 500 response over the wire. We
    # want to assert on that response, not on the raw exception.
    with TestClient(app, raise_server_exceptions=False) as c:
        response = c.post("/v1/agent/query", json={"query": "Tell me about some rare gene"})

    assert response.status_code == 500
    body = response.json()
    assert body["type"] == "internal_error"
    assert "boom" not in body["detail"]
    assert "Traceback" not in json.dumps(body)


def test_request_id_echoed_back_unchanged(client):
    response = client.post(
        "/v1/agent/query",
        json={"query": "Tell me about TP53"},
        headers={"X-Request-ID": "test-fixed-id-123"},
    )

    assert response.headers["X-Request-ID"] == "test-fixed-id-123"
    assert response.json()["request_id"] == "test-fixed-id-123"


def test_cors_preflight_allows_configured_origin(client):
    response = client.options(
        "/v1/agent/query",
        headers={
            "Origin": "http://localhost:5173",
            "Access-Control-Request-Method": "POST",
        },
    )

    assert response.status_code in (200, 204)
    assert response.headers["access-control-allow-origin"] == "http://localhost:5173"
