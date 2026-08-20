"""Full agent flow tests. Require valid AWS credentials + Bedrock access (live LLM calls)."""

import asyncio
import json
import os

import boto3
import pytest

def _has_bedrock_access() -> bool:
    try:
        boto3.client("bedrock", region_name=os.getenv("AWS_REGION", "us-east-1")).list_foundation_models()
        return True
    except Exception:
        return False

_HAS_BEDROCK = _has_bedrock_access()

pytestmark = [
    pytest.mark.live,
    pytest.mark.skipif(not _HAS_BEDROCK, reason="AWS credentials / Bedrock access not available — skipping live agent tests"),
]


@pytest.fixture(scope="module")
def executor():
    from BusinessFlow import executor as agent_executor

    return agent_executor


def _print_result(label: str, result: dict) -> None:
    steps = result.get("intermediate_steps", [])
    print(f"\n[{label}] input={result.get('input')!r}")
    for action, observation in steps:
        print(f"  tool call: {action.tool}({json.dumps(action.tool_input)})")
        try:
            print(f"  tool output: {json.dumps(json.loads(observation), indent=2)}")
        except (TypeError, json.JSONDecodeError):
            print(f"  tool output: {observation}")
    print(f"  final output: {result['output']}")


_STATUSES = {"active", "skipped", "inverted"}


def _parse_steps(output: str):
    """Pull the step contract out of a score_explainer answer.

    Mirrors api.routes._parse_methodology: tolerate a reasoning preamble and
    code fences, and return None for anything that is not the contract.
    """
    if "</think>" in output:
        output = output.rsplit("</think>", 1)[1]
    start = output.find("{")
    if start < 0:
        return None
    try:
        parsed, _ = json.JSONDecoder().raw_decode(output[start:])
    except json.JSONDecodeError:
        return None
    steps = parsed.get("steps") if isinstance(parsed, dict) else None
    return steps if isinstance(steps, list) and steps else None


@pytest.mark.asyncio
async def test_agent_calls_gene_alias_lookup(executor):
    result = await executor.ainvoke({"input": "What are the aliases for TP53?"})
    _print_result("test_agent_calls_gene_alias_lookup", result)
    steps = result.get("intermediate_steps", [])
    tool_names = {step[0].tool for step in steps}

    assert "gene_alias_lookup" in tool_names
    assert result["output"]


@pytest.mark.asyncio
async def test_agent_calls_dataset_info(executor):
    await asyncio.sleep(20)
    result = await executor.ainvoke({"input": "What is DepMap?"})
    _print_result("test_agent_calls_dataset_info", result)
    steps = result.get("intermediate_steps", [])
    tool_names = {step[0].tool for step in steps}

    assert "dataset_info" in tool_names or "depmap" in result["output"].lower()


@pytest.mark.asyncio
async def test_agent_calls_score_explainer(executor):
    """BC2. With no scoring API connected the tool returns null numerics, so
    the answer has to be the methodology — and for score_explainer that answer
    is the structured 7-step contract the UI's chart renders, not prose."""
    await asyncio.sleep(20)
    result = await executor.ainvoke(
        {"input": "Explain the score for ENSG00000141510 in ACH-000001"}
    )
    _print_result("test_agent_calls_score_explainer", result)
    steps = result.get("intermediate_steps", [])
    tool_names = {step[0].tool for step in steps}

    assert "score_explainer" in tool_names

    parsed = _parse_steps(result["output"])
    assert parsed is not None, "score_explainer answer is not the JSON step contract"
    assert len(parsed) == 7, f"expected 7 steps, got {len(parsed)}"
    assert all({"key", "value"} <= set(step) for step in parsed)
    assert all(step.get("status", "active") in _STATUSES for step in parsed)

    # TP53 is a tumour suppressor, so the inversion step applies to it.
    inversion = next((s for s in parsed if "inver" in s["key"].lower()), None)
    assert inversion is not None, "no TSG inversion step in the answer"
    assert inversion.get("status") == "inverted"

    blob = json.dumps(parsed).lower()
    for keyword in ("pit", "percentile", "correlation", "weight"):
        assert keyword in blob, f"methodology term {keyword!r} missing from the answer"


@pytest.mark.asyncio
async def test_agent_does_not_hallucinate_on_empty_result(executor):
    await asyncio.sleep(20)
    result = await executor.ainvoke(
        {"input": "Explain the score for ENSG_FAKE in ACH-FAKE"}
    )
    _print_result("test_agent_does_not_hallucinate_on_empty_result", result)
    output = result["output"].lower()

    assert "0.7" not in output and "0.702" not in output
