"""Wires the CellLineFinder agent together: LLM + tools + prompt.

The agent is built ONCE at module level. The AgentExecutor handles the full
cycle: LLM decides tool -> tool executes -> LLM formats response. No manual
routing logic.

The FastAPI surface lives in `api/app.py` / `api/routes.py`, which import
`executor` from this module. This file only builds the agent.
"""

import os
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

# R5: bound the agent loop. Without limits, a model that loops on tool
# calls burns quota and holds the connection until the client gives up.
# verbose defaults to False (H6) — it prints raw LLM output (and query
# content) straight to stdout on every request; gate it behind DEBUG.
executor = AgentExecutor(
    agent=agent,
    tools=tools,
    verbose=os.getenv("DEBUG", "false").lower() in ("1", "true", "yes"),
    return_intermediate_steps=True,
    max_iterations=4,
    max_execution_time=45.0,
)
