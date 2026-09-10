"""Response cache and fast-path routing helpers (Latency Architecture).

`detect_fast_path` matches unambiguous inputs (a bare ENSG id, an ENSG+ACH
pair, a known dataset name) so the route can call the tool directly and
skip the agent's tool-selection LLM round trip — only narration still needs
the LLM (Step 15: halves latency on the most common query shapes).
"""

from __future__ import annotations

import re

from cachetools import TTLCache

_ENSG_RE = re.compile(r"^ENSG\d{11}(?:\.\d+)?$", re.IGNORECASE)
_ENSG_ACH_RE = re.compile(
    r"^(?P<ensg>ENSG\d{11})(?:\.\d+)?\s+(?:in|for)?\s*(?P<ach>ACH-\d+)$", re.IGNORECASE
)


def normalise_query(query: str) -> str:
    return " ".join(query.strip().lower().split())


def make_response_cache(maxsize: int, ttl: int) -> TTLCache:
    return TTLCache(maxsize=maxsize, ttl=ttl)


def detect_fast_path(query: str, dataset_names: set[str]) -> tuple[str, dict] | None:
    """Return (tool_name, tool_args) for an unmistakable input, else None."""
    text = query.strip()

    combined = _ENSG_ACH_RE.match(text)
    if combined:
        return "score_explainer", {
            "ensg_id": combined.group("ensg").upper(),
            "model_id": combined.group("ach").upper(),
        }

    if _ENSG_RE.match(text):
        return "gene_alias_lookup", {"query": text.split(".")[0].upper()}

    lowered = text.lower()
    if lowered in dataset_names:
        return "dataset_info", {"dataset_name": lowered}

    return None
