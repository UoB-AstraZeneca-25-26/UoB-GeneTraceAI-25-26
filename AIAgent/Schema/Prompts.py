"""System prompt builder for the CellLineFinder agent.

Loads static content from the Knowledge/ directory and assembles it into a
single system prompt string. Built once at startup, not per-request.
"""

from pathlib import Path

# NOTE: capital K — matches the actual on-disk directory name (fixes G1;
# see AgentDevelopment/FunctionCalling.py for the same fix).
KNOWLEDGE_DIR = Path(__file__).resolve().parent.parent / "Knowledge"

_ROLE_AND_RULES = """\
You are CellLineFinder's explanation agent. You help scientists
understand gene aliases, scoring results, and data sources.

Rules:
- Always use a tool before answering. Never guess.
- For gene queries: call gene_alias_lookup.
- For "why is this cell line ranked here" queries: call score_explainer.
- For "what is this dataset / where does this data come from": call dataset_info.
- After receiving tool output, explain the result in plain language.
- For score_explainer: walk through each step of the pipeline using
  the actual numbers returned. Do not invent numbers.
- Be concise. Maximum 150 words per response. No tables.
"""


def build_system_prompt() -> str:
    """Assemble the full system prompt from role rules and the math
    reference.

    Dataset descriptions are deliberately NOT embedded here (gap G16):
    the `dataset_info` tool already returns the full entry on demand, so
    duplicating ~700 tokens/turn of static JSON only inflates every LLM
    call for no benefit.
    """
    math_reference = (KNOWLEDGE_DIR / "math_reference.md").read_text()

    sections = [
        _ROLE_AND_RULES,
        "## Mathematical reference (for score_explainer results)\n\n" + math_reference,
    ]
    return "\n\n".join(sections)
