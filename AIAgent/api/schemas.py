"""Request/response models for the /v1/agent/query surface."""

from typing import Literal

from pydantic import BaseModel, Field


class QueryRequest(BaseModel):
    query: str = Field(..., min_length=1, max_length=2000)
    session_id: str | None = Field(None, max_length=64)


class ToolCall(BaseModel):
    tool: str
    input: dict
    output: dict


class MethodologyStep(BaseModel):
    """One node of the scoring-methodology chart the UI renders.

    `status` drives the node's visual state: a step that did not apply to
    this gene (TSG inversion on an oncogene) is drawn greyed and dashed
    rather than silently dropped, so the pipeline stays auditable.
    """

    key: str                                    # chart node title (3-5 words)
    value: str                                  # explanation shown on click
    formula: str | None = None                  # short formula, always visible
    status: Literal["active", "skipped", "inverted"] = "active"


class Methodology(BaseModel):
    steps: list[MethodologyStep]


class QueryResponse(BaseModel):
    answer: str | None = None       # LLM narration; null when degraded or structured
    data: dict | None = None        # deterministic tool payload — renderable alone
    methodology: Methodology | None = None  # score_explainer only; chart-renderable
    request_id: str
    model: str
    tool_calls: list[ToolCall] = []
    degraded: list[str] = []        # e.g. ["ensembl", "llm_unavailable"]
    cached: bool = False
    latency_ms: int


class HealthResponse(BaseModel):
    status: str = "ok"
