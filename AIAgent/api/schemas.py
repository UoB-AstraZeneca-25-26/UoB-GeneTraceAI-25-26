"""Request/response models for the /v1/agent/query surface."""

from pydantic import BaseModel, Field


class QueryRequest(BaseModel):
    query: str = Field(..., min_length=1, max_length=2000)
    session_id: str | None = Field(None, max_length=64)


class ToolCall(BaseModel):
    tool: str
    input: dict
    output: dict


class QueryResponse(BaseModel):
    answer: str | None = None       # LLM narration; null when degraded
    data: dict | None = None        # deterministic tool payload — renderable alone
    request_id: str
    model: str
    tool_calls: list[ToolCall] = []
    degraded: list[str] = []        # e.g. ["ensembl", "llm_unavailable"]
    cached: bool = False
    latency_ms: int


class HealthResponse(BaseModel):
    status: str = "ok"
