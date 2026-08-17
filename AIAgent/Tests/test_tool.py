from dotenv import load_dotenv
from langchain_groq import ChatGroq
from langchain_core.tools import tool

load_dotenv()

@tool
def echo_tool(query: str) -> str:
    """A test tool that echoes input."""
    return f"echo: {query}"


def test_llm_resolves_tool_call():
    llm = ChatGroq(model="openai/gpt-oss-120b", temperature=0)
    llm_with_tools = llm.bind_tools([echo_tool])
    response = llm_with_tools.invoke("Call the echo_tool tool with 'hello'")

    assert response.tool_calls
    assert response.tool_calls[0]["name"] == "echo_tool"
    assert response.tool_calls[0]["args"] == {"query": "hello"}