"""The API surface: one business endpoint + /health.

POST /v1/agent/query negotiates on the Accept header — text/event-stream
gets an SSE stream (data before token, per the Latency Architecture in
plan.md), anything else gets a single buffered JSON body. Both paths share
one internal event generator (`_agent_events`) so the tool-selection /
fast-path / degradation logic is written exactly once.
"""

from __future__ import annotations

import json
import time
from dataclasses import dataclass, field
from typing import Any, AsyncIterator

from fastapi import APIRouter, Depends, Request
from fastapi.encoders import jsonable_encoder
from fastapi.responses import JSONResponse
from botocore.exceptions import ClientError as BotoClientError
from groq import RateLimitError as GroqRateLimitError
from langchain_core.messages import HumanMessage, SystemMessage
from pydantic import ValidationError
from sse_starlette.sse import EventSourceResponse

from . import _pathsetup  # noqa: F401  side effect: sys.path for the flat modules below

from BusinessFlow import executor, llm as _llm, system_prompt, tools  # noqa: E402
from FunctionCalling import _KB  # noqa: E402

from .cache import detect_fast_path, make_response_cache, normalise_query
from .errors import UpstreamRateLimitedError
from .schemas import HealthResponse, Methodology, QueryRequest, QueryResponse, ToolCall
from .security import limiter, require_api_key
from .settings import get_settings

router = APIRouter()

_settings = get_settings()
_response_cache = make_response_cache(_settings.response_cache_maxsize, _settings.response_cache_ttl)
_TOOLS_BY_NAME = {t.name: t for t in tools}
_DATASET_NAMES: set[str] = set(_KB.keys())

# The one tool whose answer is a structured contract rather than prose.
SCORE_EXPLAINER = "score_explainer"


# ---------------------------------------------------------------------------
# Shared event model
# ---------------------------------------------------------------------------

@dataclass
class _Event:
    kind: str  # "status" | "data" | "token" | "error"
    payload: dict[str, Any] = field(default_factory=dict)


async def _invoke_tool(tool_name: str, args: dict) -> dict:
    raw = await _TOOLS_BY_NAME[tool_name].ainvoke(args)
    return json.loads(raw) if isinstance(raw, str) else raw


def _chunk_to_text(content: Any) -> str:
    """Normalise a streamed chat-model chunk's ``.content`` to plain text.

    Groq/OpenAI-compatible backends yield ``content`` as a plain str.
    Bedrock Converse (via langchain_aws) yields a list of content blocks
    like ``[{"type": "text", "text": "...", "index": 0}]``, which
    ``"".join(...)`` chokes on downstream if left unnormalised.
    """
    if isinstance(content, str):
        return content
    if isinstance(content, list):
        return "".join(
            block.get("text", "") for block in content
            if isinstance(block, dict) and block.get("type") == "text"
        )
    if isinstance(content, dict):
        return content.get("text", "")
    return ""


# score_explainer answers are consumed by the UI's methodology chart, so the
# narration instruction has to ask for the JSON contract rather than prose --
# otherwise the fast path contradicts the system prompt and the model splits
# the difference.
_NARRATION_INSTRUCTION = (
    "Explain this result to the scientist in plain language, "
    "following your instructions."
)
_STRUCTURED_INSTRUCTION = (
    "Return the 7-step methodology as JSON, following the output contract "
    "in your instructions. Output the JSON object only."
)


async def _narrate_tool_output(tool_name: str, output: dict) -> AsyncIterator[str]:
    instruction = (
        _STRUCTURED_INSTRUCTION if tool_name == SCORE_EXPLAINER else _NARRATION_INSTRUCTION
    )
    messages = [
        SystemMessage(content=system_prompt),
        HumanMessage(
            content=(
                f"The `{tool_name}` tool returned this JSON result:\n"
                f"{json.dumps(output)}\n\n" + instruction
            )
        ),
    ]
    async for chunk in _llm.astream(messages):
        text = _chunk_to_text(getattr(chunk, "content", ""))
        if text:
            yield text


_VALID_STATUS = {"active", "skipped", "inverted"}


def _strip_reasoning(text: str) -> str:
    """Drop a reasoning model's <think> preamble.

    The default model (qwen-qwq-32b) emits its chain of thought before the
    answer; without this the JSON object never starts at a parseable offset.
    """
    marker = "</think>"
    if marker in text:
        text = text.rsplit(marker, 1)[1]
    return text.strip()


def _parse_methodology(answer_text: str | None) -> Methodology | None:
    """Read a score_explainer answer as the structured step contract.

    Returns None for anything that is not the contract -- prose, a truncated
    object, an unexpected shape -- so the caller can keep the raw text and the
    UI falls back to rendering it as a paragraph. Deliberately tolerant of the
    wrappers models add (code fences, a trailing sentence): the JSON object is
    located by hand and decoded with raw_decode so trailing text is ignored.
    """
    if not answer_text:
        return None

    text = _strip_reasoning(answer_text)
    start = text.find("{")
    if start < 0:
        return None
    try:
        parsed, _ = json.JSONDecoder().raw_decode(text[start:])
    except json.JSONDecodeError:
        return None

    if not isinstance(parsed, dict):
        return None
    steps = parsed.get("steps")
    if not isinstance(steps, list) or not steps:
        return None

    # An out-of-contract status is a cosmetic error, not a reason to throw the
    # whole chart away -- fall back to the neutral node state.
    normalised = []
    for step in steps:
        if not isinstance(step, dict):
            return None
        step = dict(step)
        if step.get("status") not in _VALID_STATUS:
            step["status"] = "active"
        normalised.append(step)

    try:
        return Methodology.model_validate({"steps": normalised})
    except ValidationError:
        return None


def _split_answer(answer: str | None, tool_calls: list[ToolCall]) -> tuple[str | None, Methodology | None]:
    """Route a score_explainer answer into `methodology`, everything else
    into `answer`. A score_explainer answer that will not parse stays prose."""
    if not any(tc.tool == SCORE_EXPLAINER for tc in tool_calls):
        return answer, None
    methodology = _parse_methodology(answer)
    return (None, methodology) if methodology is not None else (answer, None)


def _parse_tool_output(raw_output: Any) -> dict:
    text = getattr(raw_output, "content", raw_output)
    try:
        return json.loads(text) if isinstance(text, str) else (text or {})
    except (TypeError, json.JSONDecodeError):
        return {"raw": text}


async def _agent_events(query: str) -> AsyncIterator[_Event]:
    """Drives either the fast path (direct tool call + narration-only LLM
    call) or the full AgentExecutor, yielding a common event shape that both
    the SSE and buffered handlers below consume.
    """
    fast = detect_fast_path(query, _DATASET_NAMES)

    if fast is not None:
        tool_name, tool_args = fast
        yield _Event("status", {"stage": "resolving", "label": f"Looking up {tool_name.replace('_', ' ')}..."})
        try:
            output = await _invoke_tool(tool_name, tool_args)
        except Exception as exc:  # tool itself failed -> no data at all
            yield _Event("error", {"exc": exc})
            return

        yield _Event("data", {"tool": tool_name, "input": tool_args, "output": output})
        yield _Event("status", {"stage": "explaining", "label": "Writing explanation..."})
        try:
            async for text in _narrate_tool_output(tool_name, output):
                yield _Event("token", {"text": text})
        except GroqRateLimitError as exc:
            yield _Event("error", {"exc": UpstreamRateLimitedError(detail=str(exc))})
        except BotoClientError as exc:
            if exc.response.get("Error", {}).get("Code") == "ThrottlingException":
                yield _Event("error", {"exc": UpstreamRateLimitedError(detail=str(exc))})
            else:
                yield _Event("error", {"exc": exc})
        except Exception as exc:
            yield _Event("error", {"exc": exc})
        return

    try:
        async for event in executor.astream_events({"input": query}, version="v2"):
            kind = event["event"]
            if kind == "on_tool_start":
                yield _Event("status", {"stage": "resolving", "label": f"Calling {event.get('name')}..."})
            elif kind == "on_tool_end":
                output = _parse_tool_output(event["data"].get("output"))
                yield _Event("data", {
                    "tool": event.get("name", "unknown"),
                    "input": event["data"].get("input", {}) or {},
                    "output": output,
                })
                yield _Event("status", {"stage": "explaining", "label": "Writing explanation..."})
            elif kind == "on_chat_model_stream":
                chunk = event["data"].get("chunk")
                text = _chunk_to_text(getattr(chunk, "content", "")) if chunk is not None else ""
                if text:
                    yield _Event("token", {"text": text})
    except GroqRateLimitError as exc:
        yield _Event("error", {"exc": UpstreamRateLimitedError(detail=str(exc))})
    except BotoClientError as exc:
        if exc.response.get("Error", {}).get("Code") == "ThrottlingException":
            yield _Event("error", {"exc": UpstreamRateLimitedError(detail=str(exc))})
        else:
            yield _Event("error", {"exc": exc})
    except Exception as exc:
        yield _Event("error", {"exc": exc})


def _collect_degraded(tool_calls: list[ToolCall]) -> list[str]:
    degraded: list[str] = []
    for tc in tool_calls:
        for flag in tc.output.get("degraded", []) if isinstance(tc.output, dict) else []:
            if flag not in degraded:
                degraded.append(flag)
    return degraded


# ---------------------------------------------------------------------------
# Buffered (Postman / curl / tests)
# ---------------------------------------------------------------------------

async def _buffered(req: QueryRequest, request: Request) -> JSONResponse:
    settings = get_settings()
    request_id = request.state.request_id
    start = time.perf_counter()
    normalised = normalise_query(req.query)

    cached = _response_cache.get(normalised)
    if cached is not None:
        payload = cached.model_copy(update={
            "request_id": request_id,
            "cached": True,
            "latency_ms": int((time.perf_counter() - start) * 1000),
        })
        return JSONResponse(status_code=200, content=jsonable_encoder(payload))

    tool_calls: list[ToolCall] = []
    answer_chunks: list[str] = []
    terminal_error: Exception | None = None

    async for evt in _agent_events(req.query):
        if evt.kind == "data":
            tool_calls.append(ToolCall(tool=evt.payload["tool"], input=evt.payload["input"], output=evt.payload["output"]))
        elif evt.kind == "token":
            answer_chunks.append(evt.payload["text"])
        elif evt.kind == "error":
            terminal_error = evt.payload["exc"]

    if terminal_error is not None and not tool_calls:
        raise terminal_error  # no data at all -> a genuine error (503/500 via global handlers)

    degraded = _collect_degraded(tool_calls)
    if terminal_error is not None:
        degraded.append("llm_unavailable")

    answer = "".join(answer_chunks).strip() or None
    answer, methodology = _split_answer(answer, tool_calls)
    data = tool_calls[-1].output if tool_calls else None

    response = QueryResponse(
        answer=answer,
        data=data,
        methodology=methodology,
        request_id=request_id,
        model=settings.llm_model,
        tool_calls=tool_calls,
        degraded=degraded,
        cached=False,
        latency_ms=int((time.perf_counter() - start) * 1000),
    )
    if terminal_error is None:
        _response_cache[normalised] = response
    return JSONResponse(status_code=200, content=jsonable_encoder(response))


# ---------------------------------------------------------------------------
# SSE (React UI)
# ---------------------------------------------------------------------------

async def _stream_events(req: QueryRequest, request: Request) -> AsyncIterator[dict]:
    """The SSE body, factored out of `_stream` so tests can iterate it
    directly (bypassing httpx's ASGITransport, which buffers a streaming
    response until the ASGI call fully completes and so cannot measure
    real time-to-first-byte — see Tests/test_latency.py)."""
    settings = get_settings()
    request_id = request.state.request_id
    start = time.perf_counter()
    yield {"event": "accepted", "data": json.dumps({"request_id": request_id})}

    tool_calls: list[ToolCall] = []
    answer_chunks: list[str] = []
    terminal_error: Exception | None = None

    async for evt in _agent_events(req.query):
        if evt.kind == "status":
            yield {"event": "status", "data": json.dumps(evt.payload)}
        elif evt.kind == "data":
            tool_calls.append(ToolCall(tool=evt.payload["tool"], input=evt.payload["input"], output=evt.payload["output"]))
            yield {"event": "data", "data": json.dumps(evt.payload["output"])}
        elif evt.kind == "token":
            answer_chunks.append(evt.payload["text"])
            # A score_explainer answer is a JSON object, not prose: half a
            # JSON object is nothing a client can render, so those tokens are
            # withheld and the parsed steps ship in the `done` event instead.
            if not any(tc.tool == SCORE_EXPLAINER for tc in tool_calls):
                yield {"event": "token", "data": json.dumps(evt.payload)}
        elif evt.kind == "error":
            terminal_error = evt.payload["exc"]

    if terminal_error is not None and not tool_calls:
        status = 503 if isinstance(terminal_error, UpstreamRateLimitedError) else 500
        yield {"event": "error", "data": json.dumps({
            "type": "rate_limited" if status == 503 else "internal_error",
            "status": status,
            "detail": str(terminal_error),
            "request_id": request_id,
        })}
        return

    degraded = _collect_degraded(tool_calls)
    if terminal_error is not None:
        degraded.append("llm_unavailable")

    answer = "".join(answer_chunks).strip() or None
    answer, methodology = _split_answer(answer, tool_calls)
    data = tool_calls[-1].output if tool_calls else None
    response = QueryResponse(
        answer=answer,
        data=data,
        methodology=methodology,
        request_id=request_id,
        model=settings.llm_model,
        tool_calls=tool_calls,
        degraded=degraded,
        cached=False,
        latency_ms=int((time.perf_counter() - start) * 1000),
    )
    if terminal_error is None:
        _response_cache[normalise_query(req.query)] = response
    yield {"event": "done", "data": response.model_dump_json()}


def _stream(req: QueryRequest, request: Request) -> EventSourceResponse:
    # ping=10: nginx's default proxy_read_timeout is 60s and most managed
    # platforms kill idle connections around the same mark; the heartbeat
    # keeps a slow LLM response from getting the connection severed.
    return EventSourceResponse(
        _stream_events(req, request), ping=10, headers={"X-Accel-Buffering": "no"}
    )


# ---------------------------------------------------------------------------
# Routes
# ---------------------------------------------------------------------------

async def _agent_query(req: QueryRequest, request: Request):
    accept = request.headers.get("accept", "")
    if "text/event-stream" in accept:
        return _stream(req, request)
    return await _buffered(req, request)


@router.post("/v1/agent/query", response_model=None, dependencies=[Depends(require_api_key)])
@limiter.limit(_settings.rate_limit)
async def agent_query_v1(req: QueryRequest, request: Request):
    return await _agent_query(req, request)


@router.post("/agent/query", response_model=None, dependencies=[Depends(require_api_key)], deprecated=True)
@limiter.limit(_settings.rate_limit)
async def agent_query_deprecated(req: QueryRequest, request: Request):
    response = await _agent_query(req, request)
    response.headers["Deprecation"] = "true"
    response.headers["Link"] = '</v1/agent/query>; rel="successor-version"'
    return response


@router.get("/health", response_model=HealthResponse)
async def health() -> HealthResponse:
    # Dependency-free by design (H6/G9): pinging Groq or Ensembl here would
    # let an upstream blip make the orchestrator kill a healthy container.
    return HealthResponse(status="ok")
