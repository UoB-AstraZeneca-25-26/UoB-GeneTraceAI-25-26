"""Full agent flow tests. Require a real GROQ_API_KEY (live LLM calls)."""

import asyncio
import json
import os

import pytest
from dotenv import load_dotenv

load_dotenv()

_HAS_REAL_KEY = bool(os.getenv("GROQ_API_KEY")) and os.getenv("GROQ_API_KEY") != "your_groq_api_key_here"

pytestmark = pytest.mark.skipif(
    not _HAS_REAL_KEY,
    reason="GROQ_API_KEY not set — skipping live agent tests",
)


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
    await asyncio.sleep(20)
    result = await executor.ainvoke(
        {"input": "Explain the score for ENSG00000141510 in ACH-000001"}
    )
    _print_result("test_agent_calls_score_explainer", result)
    steps = result.get("intermediate_steps", [])
    tool_names = {step[0].tool for step in steps}

    assert "score_explainer" in tool_names


@pytest.mark.asyncio
async def test_agent_does_not_hallucinate_on_empty_result(executor):
    await asyncio.sleep(20)
    result = await executor.ainvoke(
        {"input": "Explain the score for ENSG_FAKE in ACH-FAKE"}
    )
    _print_result("test_agent_does_not_hallucinate_on_empty_result", result)
    output = result["output"].lower()

    assert "0.7" not in output and "0.702" not in output
