import pytest
from langchain_aws import ChatBedrock
from langchain_core.tools import tool


@tool
def echo_tool(query: str) -> str:
    """A test tool that echoes input."""
    return f"echo: {query}"


@pytest.mark.live
def test_llm_resolves_tool_call():
    llm = ChatBedrock(
        model_id="anthropic.claude-haiku-4-5",
        model_kwargs={"temperature": 0, "max_tokens": 512},
        region_name="us-east-1",
    )
    llm_with_tools = llm.bind_tools([echo_tool])
    response = llm_with_tools.invoke("Call the echo_tool tool with 'hello'")

    assert response.tool_calls
    assert response.tool_calls[0]["name"] == "echo_tool"
    assert response.tool_calls[0]["args"] == {"query": "hello"}
