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


# ---------------------------------------------------------------------------
# score_explainer -> QueryResponse.methodology (the UI's step chart)
# ---------------------------------------------------------------------------

SCORE_OUTPUT = {
    "found": True, "ensg_id": "ENSG00000157764", "model_id": "ACH-000219",
    "core_score": None, "confidence_tier": None, "is_tsg": False,
    "driver_gated": False, "message": "Methodology explanation.",
}


def _methodology_json(status: str = "skipped") -> str:
    return json.dumps({
        "steps": [
            {"key": f"Step {i}", "value": f"Explanation number {i}.", "formula": "f(x)",
             "status": status if i == 5 else "active"}
            for i in range(1, 8)
        ]
    })


def _score_client(monkeypatch, *answer_chunks):
    events = [_tool_end_event("score_explainer", {"ensg_id": "ENSG00000157764", "model_id": "ACH-000219"}, SCORE_OUTPUT)]
    events += [_token_event(chunk) for chunk in answer_chunks]
    monkeypatch.setattr(routes_module, "executor", FakeExecutor(events=events))
    routes_module._response_cache.clear()
    return TestClient(app)


def test_score_explainer_json_answer_becomes_methodology(monkeypatch):
    with _score_client(monkeypatch, _methodology_json()) as c:
        response = c.post("/v1/agent/query", json={"query": "Explain the score for ENSG00000157764 in ACH-000219"})

    body = QueryResponse.model_validate(response.json())
    assert body.answer is None                      # prose suppressed: the UI draws the chart
    assert len(body.methodology.steps) == 7
    assert body.methodology.steps[4].status == "skipped"
    assert body.data == SCORE_OUTPUT                # deterministic payload still returned


def test_score_explainer_json_survives_fences_and_reasoning(monkeypatch):
    wrapped = "<think>weighing the layers</think>\n```json\n" + _methodology_json("inverted") + "\n```"
    with _score_client(monkeypatch, wrapped) as c:
        response = c.post("/v1/agent/query", json={"query": "Explain the score for ENSG00000141510 in ACH-000219"})

    body = QueryResponse.model_validate(response.json())
    assert body.methodology is not None
    assert body.methodology.steps[4].status == "inverted"


def test_score_explainer_prose_falls_back_to_answer(monkeypatch):
    with _score_client(monkeypatch, "Step 1 - raw measurements are ", "converted to percentiles.") as c:
        response = c.post("/v1/agent/query", json={"query": "Explain the score for ENSG00000157764 in ACH-000000"})

    body = QueryResponse.model_validate(response.json())
    assert body.methodology is None
    assert body.answer == "Step 1 - raw measurements are converted to percentiles."


def test_non_score_tool_never_populates_methodology(client):
    response = client.post("/v1/agent/query", json={"query": "Tell me about TP53"})

    body = QueryResponse.model_validate(response.json())
    assert body.methodology is None
    assert body.answer == "TP53 is also known as p53."


@pytest.mark.asyncio
async def test_sse_withholds_json_tokens_and_ships_methodology_on_done(monkeypatch):
    """Half a JSON object is not renderable text, so the SSE path must not
    stream score_explainer tokens — the parsed steps ride the `done` event."""
    from api.schemas import QueryRequest

    chunks = _methodology_json()
    events = [_tool_end_event("score_explainer", {"ensg_id": "ENSG00000157764", "model_id": "ACH-000219"}, SCORE_OUTPUT)]
    events += [_token_event(chunks[i:i + 60]) for i in range(0, len(chunks), 60)]
    monkeypatch.setattr(routes_module, "executor", FakeExecutor(events=events))
    routes_module._response_cache.clear()

    request = SimpleNamespace(state=SimpleNamespace(request_id="sse-test-id"))
    emitted = [
        evt async for evt in routes_module._stream_events(
            QueryRequest(query="Explain the score for ENSG00000157764 in ACH-000219"), request
        )
    ]

    # FakeExecutor emits no on_tool_start, so there is no leading "resolving"
    # status: the tool result lands first, then "explaining", then done.
    assert [e["event"] for e in emitted] == ["accepted", "data", "status", "done"]
    body = QueryResponse.model_validate_json(emitted[-1]["data"])
    assert body.answer is None
    assert len(body.methodology.steps) == 7
