"""Tests for the score_explainer tool.

The tool has two modes, both covered here offline:
- methodology mode (SCORING_API_URL unset) — the pair's identity with every
  numeric field null, so the LLM explains the method rather than a result;
- API mode (SCORING_API_URL set) — the fields come from the scoring API,
  mocked here so nothing hits the network.

There are no fixture scores anywhere: a number in a ScoreResult only ever
comes from the API.

score_explainer does not gate on the local gene index -- gene validation is
gene_alias_lookup's job. It always returns found=True with either real
scoring-API values or a methodology explanation, regardless of whether the
gene is in the local 19,213-gene panel.
"""

import json

import pytest

import FunctionCalling as fc
from FunctionCalling import ScoreResult, score_explainer

_NUMERIC_FIELDS = (
    "raw_expression",
    "raw_proteomics",
    "pit_expression",
    "pit_proteomics",
    "rho_ep",
    "weight_expression",
    "weight_proteomics",
    "core_score",
)

# TP53 — happens to be in the reference panel, but that no longer matters
# to score_explainer: it does not consult the local index at all.
_KNOWN_GENE = "ENSG00000141510"
_KNOWN_MODEL = "ACH-000001"

# Deliberately not a real Ensembl ID / not in the local 19,213-gene panel.
_UNKNOWN_GENE = "ENSG_NOT_IN_PANEL"
_UNKNOWN_MODEL = "ACH-999999"


@pytest.fixture(autouse=True)
def _no_scoring_api(monkeypatch):
    """Default every test to methodology mode; the API-mode tests opt in.

    Without this, a SCORING_API_URL in the developer's environment would
    silently turn the offline tests into live ones.
    """
    monkeypatch.delenv("SCORING_API_URL", raising=False)
    monkeypatch.setattr(fc, "_scoring_api_url", lambda: "")


# ---------------------------------------------------------------------------
# Mode B — no scoring API configured
# ---------------------------------------------------------------------------

@pytest.mark.asyncio
async def test_methodology_mode_returns_null_numerics():
    raw = await score_explainer.ainvoke({"ensg_id": _KNOWN_GENE, "model_id": _KNOWN_MODEL})
    data = json.loads(raw)
    print(f"\n[test_methodology_mode_returns_null_numerics] {json.dumps(data, indent=2)}")

    assert data["found"] is True
    assert data["ensg_id"] == _KNOWN_GENE
    assert data["model_id"] == _KNOWN_MODEL
    for field in (*_NUMERIC_FIELDS, "confidence_tier"):
        assert data[field] is None, f"{field} must be null without a scoring API"
    assert "Methodology" in data["message"]


@pytest.mark.asyncio
async def test_methodology_mode_leaves_is_tsg_null():
    """is_tsg can no longer be determined without a scoring API, since
    score_explainer does not consult the local gene index."""
    raw = await score_explainer.ainvoke({"ensg_id": _KNOWN_GENE, "model_id": _KNOWN_MODEL})
    data = json.loads(raw)

    assert data["is_tsg"] is None
    assert data["driver_gated"] is False


@pytest.mark.asyncio
async def test_gene_outside_local_panel_still_returns_found_true():
    """Gene validation belongs to gene_alias_lookup, not score_explainer.
    A gene missing from the local 19,213-gene panel must still get a
    methodology explanation, not found=False."""
    raw = await score_explainer.ainvoke({"ensg_id": _UNKNOWN_GENE, "model_id": _UNKNOWN_MODEL})
    data = json.loads(raw)
    print(f"\n[test_gene_outside_local_panel_still_returns_found_true] {json.dumps(data, indent=2)}")

    assert data["found"] is True
    assert data["ensg_id"] == _UNKNOWN_GENE
    assert data["model_id"] == _UNKNOWN_MODEL
    assert "Methodology" in data["message"]
    for field in _NUMERIC_FIELDS:
        assert data[field] is None


@pytest.mark.asyncio
async def test_pydantic_validation_field_types():
    raw = await score_explainer.ainvoke({"ensg_id": _KNOWN_GENE, "model_id": _KNOWN_MODEL})
    result = ScoreResult.model_validate_json(raw)
    print(f"\n[test_pydantic_validation_field_types] {result.model_dump_json(indent=2)}")

    assert isinstance(result.found, bool)
    assert result.is_tsg is None
    assert isinstance(result.driver_gated, bool)
    assert isinstance(result.message, str)
    assert result.core_score is None


# ---------------------------------------------------------------------------
# Mode A — scoring API configured
# ---------------------------------------------------------------------------

_API_PAYLOAD = {
    "raw_expression": 1.5,
    "raw_proteomics": -0.25,
    "pit_expression": 0.5,
    "pit_proteomics": 0.25,
    "rho_ep": 0.2,
    "weight_expression": 0.5,
    "weight_proteomics": 0.5,
    "core_score": 0.375,
    "confidence_tier": "moderate",
    "is_tsg": True,
    "driver_gated": True,
}


@pytest.fixture()
def scoring_api(monkeypatch):
    """Point the tool at a scoring API and capture the URL it calls."""
    monkeypatch.setattr(fc, "_scoring_api_url", lambda: "https://scoring.internal/v1/score")
    calls: list[str] = []

    def _install(payload):
        async def fake_safe_get(url, *, headers=None, client=None, attempts=3):
            calls.append(url)
            return payload

        monkeypatch.setattr(fc, "safe_get", fake_safe_get)
        return calls

    return _install


@pytest.mark.asyncio
async def test_api_mode_populates_fields_from_response(scoring_api):
    calls = scoring_api(_API_PAYLOAD)

    raw = await score_explainer.ainvoke({"ensg_id": _KNOWN_GENE, "model_id": _KNOWN_MODEL})
    data = json.loads(raw)
    print(f"\n[test_api_mode_populates_fields_from_response] {json.dumps(data, indent=2)}")

    assert calls == [f"https://scoring.internal/v1/score?gene={_KNOWN_GENE}&model={_KNOWN_MODEL}"]
    assert data["found"] is True
    for field, expected in _API_PAYLOAD.items():
        assert data[field] == expected
    assert "Methodology" not in data["message"]


@pytest.mark.asyncio
async def test_api_mode_leaves_omitted_fields_null(scoring_api):
    """A partial payload must not be back-filled with anything invented."""
    scoring_api({"pit_expression": 0.5, "core_score": 0.5})

    raw = await score_explainer.ainvoke({"ensg_id": _KNOWN_GENE, "model_id": _KNOWN_MODEL})
    data = json.loads(raw)

    assert data["pit_expression"] == 0.5
    assert data["core_score"] == 0.5
    assert data["raw_expression"] is None
    assert data["rho_ep"] is None
    assert data["confidence_tier"] is None


@pytest.mark.asyncio
async def test_api_unreachable_falls_back_to_methodology(scoring_api):
    """safe_get returns None on a dead upstream — degrade to methodology
    mode rather than surfacing an error or a stale number."""
    scoring_api(None)

    raw = await score_explainer.ainvoke({"ensg_id": _KNOWN_GENE, "model_id": _KNOWN_MODEL})
    data = json.loads(raw)
    print(f"\n[test_api_unreachable_falls_back_to_methodology] {json.dumps(data, indent=2)}")

    assert data["found"] is True
    assert "Methodology" in data["message"]
    for field in _NUMERIC_FIELDS:
        assert data[field] is None


@pytest.mark.asyncio
async def test_api_mode_still_queries_gene_outside_local_panel(scoring_api):
    """No local-panel gate: the scoring API is queried regardless of
    whether the gene is in the local index, and its answer wins."""
    calls = scoring_api(_API_PAYLOAD)

    raw = await score_explainer.ainvoke({"ensg_id": _UNKNOWN_GENE, "model_id": _UNKNOWN_MODEL})
    data = json.loads(raw)

    assert data["found"] is True
    assert calls == [f"https://scoring.internal/v1/score?gene={_UNKNOWN_GENE}&model={_UNKNOWN_MODEL}"]


@pytest.mark.asyncio
async def test_api_mode_through_real_safe_get(monkeypatch):
    """End-to-end through http_utils.safe_get with an httpx MockTransport,
    so the URL construction and JSON handling are exercised for real."""
    import httpx

    import http_utils

    monkeypatch.setattr(fc, "_scoring_api_url", lambda: "https://scoring.internal/v1/score")
    seen: list[httpx.URL] = []

    def handler(request: httpx.Request) -> httpx.Response:
        seen.append(request.url)
        return httpx.Response(200, json=_API_PAYLOAD)

    http_utils.set_http_client(httpx.AsyncClient(transport=httpx.MockTransport(handler)))

    raw = await score_explainer.ainvoke({"ensg_id": _KNOWN_GENE, "model_id": _KNOWN_MODEL})
    data = json.loads(raw)

    assert seen[0].params["gene"] == _KNOWN_GENE
    assert seen[0].params["model"] == _KNOWN_MODEL
    assert data["core_score"] == _API_PAYLOAD["core_score"]
    assert data["confidence_tier"] == "moderate"


# ---------------------------------------------------------------------------
# Regression guard
# ---------------------------------------------------------------------------

def test_no_hardcoded_scores_in_module():
    """The mock score table is gone for good: no scoring numbers may be
    baked into the tool module."""
    source = fc.__file__
    with open(source, encoding="utf-8") as handle:
        text = handle.read()

    assert "_MOCK_SCORES" not in text
    for value in ("4.32", "0.85", "0.702", "0.72", "0.68", "0.41"):
        assert value not in text, f"hardcoded scoring value {value} in {source}"
