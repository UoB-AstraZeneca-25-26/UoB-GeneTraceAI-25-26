"""Wires the CellLineFinder agent together and exposes it via FastAPI.

The agent (LLM + tools + prompt) is built ONCE at module level. The
AgentExecutor handles the full cycle: LLM decides tool -> tool executes ->
LLM formats response. No manual routing logic.
"""

import sys
from pathlib import Path

_AGENT_DEV_DIR = Path(__file__).resolve().parent.parent / "AgentDevelopment"
sys.path.insert(0, str(_AGENT_DEV_DIR))

from langchain_classic.agents import AgentExecutor, create_tool_calling_agent
from langchain_core.messages import SystemMessage
from langchain_core.prompts import ChatPromptTemplate, MessagesPlaceholder

from FunctionCalling import dataset_info, gene_alias_lookup, score_explainer
from llm import get_llm
from Prompts import build_system_prompt

tools = [gene_alias_lookup, score_explainer, dataset_info]
llm = get_llm()
system_prompt = build_system_prompt()

# System prompt is passed as a SystemMessage, not a ("system", text) template
# tuple, so the embedded JSON's literal { } braces are never run through the
# f-string template parser (which chokes on nested braces).
prompt = ChatPromptTemplate.from_messages([
    SystemMessage(content=system_prompt),
    ("human", "{input}"),
    MessagesPlaceholder("agent_scratchpad"),
])

agent = create_tool_calling_agent(llm, tools, prompt)
executor = AgentExecutor(
    agent=agent, tools=tools, verbose=True, return_intermediate_steps=True
)


from fastapi import FastAPI
from pydantic import BaseModel

app = FastAPI(title="CellLineFinder Agent")


class QueryRequest(BaseModel):
    query: str


class QueryResponse(BaseModel):
    answer: str


@app.post("/agent/query", response_model=QueryResponse)
async def agent_query(request: QueryRequest) -> QueryResponse:
    result = await executor.ainvoke({"input": request.query})
    return QueryResponse(answer=result["output"])


if __name__ == "__main__":
    import uvicorn

    uvicorn.run(app, host="0.0.0.0", port=8000)
