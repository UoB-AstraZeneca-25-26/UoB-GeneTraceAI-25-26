"""Live tests against the deployed Lambda Function URL.

Unlike test_agent.py (which invokes the AgentExecutor in-process and
requires local AWS Bedrock credentials), these tests exercise the
already-deployed Lambda over HTTP. The Lambda has its own execution
role with Bedrock permissions, so no local credentials are needed —
only network access to the Function URL.

Run with: pytest Tests/test_agent_lambda.py -v
"""

import json
import os

import httpx
import pytest

LAMBDA_URL = os.getenv(
    "AGENT_LAMBDA_URL",
    "https://mhyvpwmjma3gqy44dk4wskkxje0gqfud.lambda-url.eu-west-2.on.aws",
).rstrip("/")

# Cold starts have measured ~25s init + a few seconds for the LLM call;
# give real headroom rather than flaking on a slow start.
_TIMEOUT = httpx.Timeout(60.0, connect=10.0)

pytestmark = pytest.mark.live


def _post(query: str) -> dict:
    """POST a query to the deployed agent and return the parsed JSON body."""
    resp = httpx.post(
        f"{LAMBDA_URL}/v1/agent/query",
        json={"query": query},
        headers={"Content-Type": "application/json"},
        timeout=_TIMEOUT,
    )
    resp.raise_for_status()
    return resp.json()


def _print_result(name: str, query: str, result: dict) -> None:
    print(f"\n[{name}] input={query!r}")
    print(f"  tool_calls: {json.dumps(result.get('tool_calls', []), indent=2)}")
    print(f"  answer: {result.get('answer')}")
    print(f"  data: {json.dumps(result.get('data'), indent=2)}")
    print(f"  degraded: {result.get('degraded')}")
    print(f"  latency_ms: {result.get('latency_ms')}")


def test_lambda_health_check():
    """Sanity check the Function URL is reachable before running the
    business-case tests below — fail fast with a clear message if the
    endpoint itself is down, rather than four confusing timeouts."""
    resp = httpx.get(f"{LAMBDA_URL}/health", timeout=_TIMEOUT)
    assert resp.status_code == 200


def test_lambda_calls_gene_alias_lookup():
    """BC1 — gene alias resolution, against the deployed Lambda."""
    query = "What are the aliases for TP53"
    result = _post(query)
    _print_result("test_lambda_calls_gene_alias_lookup", query, result)

    tool_names = {tc["tool"] for tc in result.get("tool_calls", [])}
    assert "gene_alias_lookup" in tool_names

    data = result.get("data", {})
    assert data.get("found") is True
    assert data.get("symbol") == "TP53"
    assert data.get("ensembl_id") == "ENSG00000141510"

    answer = (result.get("answer") or "").lower()
    assert "p53" in answer or "tp53" in answer


def test_lambda_calls_dataset_info():
    """BC3 — dataset provenance, against the deployed Lambda."""
    query = "What is DepMap?"
    result = _post(query)
    _print_result("test_lambda_calls_dataset_info", query, result)

    tool_names = {tc["tool"] for tc in result.get("tool_calls", [])}
    assert "dataset_info" in tool_names

    data = result.get("data", {})
    assert data.get("name") == "depmap"

    answer = (result.get("answer") or "").lower()
    assert "depmap" in answer or "dependency" in answer


def test_lambda_calls_score_explainer():
    """BC2 — methodology explanation, against the deployed Lambda."""
    query = "Explain the score for ENSG00000141510 in ACH-000001"
    result = _post(query)
    _print_result("test_lambda_calls_score_explainer", query, result)

    tool_names = {tc["tool"] for tc in result.get("tool_calls", [])}
    assert "score_explainer" in tool_names

    # Structured methodology output, if the agent returned it, otherwise
    # fall back to checking the prose answer references the pipeline.
    methodology = result.get("methodology")
    if methodology and methodology.get("steps"):
        steps = methodology["steps"]
        assert len(steps) == 7
        blob = json.dumps(steps).lower()
        # The corrected pipeline description: robust z, within-lineage
        # transform, and a step 5 that combines the residualised protein
        # signal rather than raw prot_z.
        assert "robust_z" in blob, "step 1 is a robust z-score, not log2(TPM+1)"
        assert "lineage" in blob, "step 2 ranks within lineage, not the whole panel"
        assert "resid" in blob or "e*" in blob, "step 5 must combine e*, not raw prot_z"
    else:
        answer = (result.get("answer") or "").lower()
        assert "core" in answer or "lineage" in answer or "robust" in answer


def test_lambda_does_not_hallucinate_on_unknown_pair():
    """BC4 — hallucination guard, against the deployed Lambda."""
    query = "Explain the score for ENSG_FAKE in ACH-FAKE"
    result = _post(query)
    _print_result("test_lambda_does_not_hallucinate_on_unknown_pair", query, result)

    answer = (result.get("answer") or "").lower()
    # Must not fabricate a numeric core score for a nonexistent pair.
    assert "0." not in answer or "no scoring data" in answer or "not found" in answer


def test_lambda_response_shape():
    """Confirm the deployed response matches the documented contract —
    guards against a future change silently breaking the frontend."""
    result = _post("What is DepMap?")
    for field in ("answer", "data", "tool_calls", "degraded", "cached", "latency_ms", "request_id"):
        assert field in result, f"missing expected field: {field}"
    assert isinstance(result["tool_calls"], list)
    assert isinstance(result["degraded"], list)
    assert isinstance(result["latency_ms"], int)
