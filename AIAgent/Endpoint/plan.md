# CellLineFinder AI Agent — Implementation Plan

## Overview
**Stack:** LangChain + ChatGroq (Qwen 2.5) + FastAPI + Pydantic v2 + httpx
**File structure** (actual on-disk layout — note the capitalised directory names):

```
AIAgent/
├── AgentDevelopment/
│   ├── llm.py               # LLM backend config + initialisation
│   └── FunctionCalling.py   # 3 tool definitions + Pydantic I/O schemas
├── Schema/
│   ├── Prompts.py           # System prompt builder
│   └── BusinessFlow.py      # AgentExecutor wiring + FastAPI app
├── Endpoint/
│   └── plan.md              # This file
├── Knowledge/               # NOTE: capital K — see "Hosting" section
│   ├── datasets.json        # Pre-extracted dataset descriptions (static)
│   └── math_reference.md    # Scoring pipeline equations (static)
├── Tests/
│   ├── test_gene_alias.py
│   ├── test_score_explainer.py
│   ├── test_dataset_info.py
│   ├── test_tool.py
│   └── test_agent.py
├── main.py                  # Entrypoint (currently EMPTY — see Gap G7)
├── pytest.ini
└── requirement_agents.txt
```

**Planned additions for production hosting** (see the Production Readiness sections):
```
AIAgent/
├── api/
│   ├── app.py               # FastAPI app factory + lifespan
│   ├── routes.py            # Versioned /v1 routes
│   ├── schemas.py           # Request/response models
│   ├── errors.py            # Exception handlers + error envelope
│   ├── middleware.py        # Request ID, timing, CORS
│   └── settings.py          # pydantic-settings config
├── Dockerfile
└── requirements.lock.txt    # Pinned, hash-checked
```

---

## File 1: `llm.py`

### Purpose
Single function that returns a configured LangChain ChatModel.
Swap provider by changing env vars — no code change.

### Implementation steps

1. Load env vars using `python-dotenv`:
   - `LLM_PROVIDER` — one of: `groq`, `huggingface`, `openai`
   - `LLM_MODEL` — model string (default: `qwen-qwq-32b`)
   - `LLM_TEMPERATURE` — float (default: `0.0`)
   - `GROQ_API_KEY`
   - `HUGGINGFACE_API_KEY` (optional, for future swap)
   - `OPENAI_API_KEY` (optional, for future swap)

2. Define function `get_llm()` that returns a `BaseChatModel`:
   ```python
   def get_llm() -> BaseChatModel:
   ```
   - If provider is `groq`: return `ChatGroq(model=..., temperature=..., api_key=...)`
   - If provider is `huggingface`: return `ChatHuggingFace(...)` (future)
   - If provider is `openai`: return `ChatOpenAI(...)` (future)
   - Else: raise `ValueError`

3. Create `.env` file in project root with:
   ```
   LLM_PROVIDER=groq
   LLM_MODEL=qwen-qwq-32b
   LLM_TEMPERATURE=0.0
   GROQ_API_KEY=<key>
   ```

### Constraints
- No global state. `get_llm()` reads env fresh each call.
- Only `langchain-groq` is required now. Other providers are if/elif stubs.
- Function must be importable by `BusinessFlow.py`.

---

## File 2: `FunctionCalling.py`

### Purpose
Defines 3 tool functions with `@tool` decorator + Pydantic I/O schemas.
Each tool is self-contained — no cross-tool dependencies.

### Tool 1: `gene_alias_lookup`

**What it does:** Takes a gene name/symbol/Ensembl ID, queries Ensembl and HGNC REST APIs, cross-validates, returns structured aliases.

**Pydantic output schema — `GeneAliasResult`:**
```python
class GeneAliasResult(BaseModel):
    found: bool
    query: str                          # What the user typed
    symbol: str | None                  # HGNC-approved symbol
    ensembl_id: str | None              # Bare ENSG ID (no dot version)
    full_name: str | None               # Full gene name
    previous_symbols: list[str]         # Retired symbols
    synonyms: list[str]                 # Known aliases
    cross_validated: bool               # True if Ensembl + HGNC agree
    sources: list[str]                  # ["ensembl", "hgnc"]
```

**Implementation steps:**

1. Define async helper `_query_ensembl(query: str) -> dict | None`:
   - URL: `https://rest.ensembl.org/lookup/symbol/homo_sapiens/{query}`
   - Headers: `{"Content-Type": "application/json"}`
   - Use `httpx.AsyncClient`, timeout 10s
   - If 200: return JSON. Else: return None.
   - Extract: `id`, `display_name`, `description`

2. Define async helper `_query_ensembl_by_id(ensg_id: str) -> dict | None`:
   - URL: `https://rest.ensembl.org/lookup/id/{ensg_id}`
   - Same pattern as above.
   - This handles the case where user provides an Ensembl ID directly.

3. Define async helper `_query_hgnc(query: str) -> dict | None`:
   - URL: `https://rest.genenames.org/search/symbol/{query}`
   - Headers: `{"Accept": "application/json"}`
   - Parse: `response.docs[0]` if exists
   - Extract: `symbol`, `name`, `ensembl_gene_id`, `prev_symbol`, `alias_symbol`

4. Define the tool:
   ```python
   @tool
   async def gene_alias_lookup(query: str) -> str:
       """Look up gene aliases and synonyms. Takes a gene name, symbol,
       or Ensembl ID. Returns validated aliases from Ensembl and HGNC."""
   ```
   - Detect if query looks like an Ensembl ID (starts with "ENSG") → use `_query_ensembl_by_id`
   - Otherwise → use `_query_ensembl` by symbol
   - Always call `_query_hgnc` in parallel
   - Cross-validate: if both return an ensembl_id, check they match
   - Build `GeneAliasResult`, return `.model_dump_json()`

### Tool 2: `score_explainer`

**What it does:** Takes an ensg_id + model_id pair, reads pre-computed scores from pipeline parquet output, returns intermediate values for LLM to explain.

**Pydantic output schema — `ScoreResult`:**
```python
class ScoreResult(BaseModel):
    found: bool
    ensg_id: str
    model_id: str
    raw_expression: float | None        # log2(TPM+1)
    raw_proteomics: float | None        # log2 ratio
    pit_expression: float | None        # PIT percentile rank
    pit_proteomics: float | None        # PIT percentile rank
    rho_ep: float | None                # Pearson correlation E vs P
    weight_expression: float | None     # Correlation-penalized weight
    weight_proteomics: float | None     # Correlation-penalized weight
    core_score: float | None            # Weighted sum
    confidence_tier: str | None         # high | medium | low | insufficient
    is_tsg: bool                        # Tumour suppressor flag
    driver_gated: bool                  # Driver gating applied
    message: str
```

**Implementation steps:**

1. Define helper `_load_scores(ensg_id: str, model_id: str) -> ScoreResult`:
   - Read from pipeline output parquet file (path from env or config)
   - Filter by ensg_id + model_id
   - If found: populate all fields from row
   - If not found: return `ScoreResult(found=False, message="...")`
   - NOTE: for initial scaffold, use a hardcoded mock dict for testing.
     Replace with actual parquet read when pipeline output is ready.

2. Define the tool:
   ```python
   @tool
   def score_explainer(ensg_id: str, model_id: str) -> str:
       """Retrieve scoring intermediates for a gene–cell line pair.
       Returns raw values, PIT ranks, weights, core score, and tier."""
   ```
   - Call `_load_scores()`
   - Return `.model_dump_json()`

### Tool 3: `dataset_info`

**What it does:** Takes a dataset name, looks it up in `knowledge/datasets.json`, returns the description.

**Implementation steps:**

1. Load `knowledge/datasets.json` at module level (loaded once, static):
   ```python
   _KB = json.loads(Path("knowledge/datasets.json").read_text())
   ```

2. Define the tool:
   ```python
   @tool
   def dataset_info(dataset_name: str) -> str:
       """Get a description of a CellLineFinder data source.
       Valid: depmap, hpa, geo, cellosaurus, ensembl, cosmic, hgnc"""
   ```
   - Normalise key: `dataset_name.lower().strip()`
   - Exact match → return JSON
   - Substring match → return JSON
   - No match → return available keys

### Constraints for all tools
- Every `@tool` function returns a `str` (JSON-serialised Pydantic model).
- Docstrings are critical — the LLM reads them to decide which tool to call.
- No tool calls another tool. No shared state.
- `gene_alias_lookup` is async. `score_explainer` and `dataset_info` are sync.

---

## File 3: `Prompts.py`

### Purpose
Builds the system prompt string. Loads static content from `knowledge/` directory.

### Implementation steps

1. Define function `build_system_prompt() -> str` that constructs the full prompt.

2. The system prompt has 3 sections:

   **Section A — Role and behaviour:**
   ```
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
   ```

   **Section B — Mathematical reference (for score explanation):**
   Load from `knowledge/math_reference.md` and embed verbatim.
   This file contains the 7-step scoring pipeline equations:
   1. Raw inputs: log2(TPM+1) for expression, log2 ratio for proteomics
   2. PIT normalisation: F̂(x) = rank(x) / (n+1), per gene across cell lines
   3. Correlation-penalized weights: w_raw = 1/(1+mean|ρ|), normalised
   4. Core score: weighted sum of PIT ranks
   5. TSG inversion: if gene is TSG, invert rank (low expression = high score)
   6. Driver gating: filter by known cancer driver gene lists
   7. Confidence tier assignment: based on data coverage and concordance

   NOTE: Create `knowledge/math_reference.md` with these equations written out.
   The LLM reads this to explain what each number means — it never computes.

   **Section C — Dataset descriptions (for dataset_info context):**
   Load `knowledge/datasets.json` and embed as a JSON block.
   This allows the LLM to answer dataset questions even without calling the
   tool, since the full knowledge base is small enough (~1500 tokens).

3. Return the concatenated string.

### Constraints
- System prompt is built once at startup, not per-request.
- Total prompt size target: under 3000 tokens.
- No f-strings with user input in the system prompt — user input goes in the message, never the system prompt.

---

## File 4: `BusinessFlow.py`

### Purpose
Wires everything together: creates the agent, exposes it via FastAPI.

### Implementation steps

1. **Imports:**
   ```python
   from llm import get_llm
   from FunctionCalling import gene_alias_lookup, score_explainer, dataset_info
   from Prompts import build_system_prompt
   ```

2. **Build the agent** (module-level, created once):
   ```python
   from langchain_core.prompts import ChatPromptTemplate, MessagesPlaceholder
   from langchain.agents import create_tool_calling_agent, AgentExecutor

   tools = [gene_alias_lookup, score_explainer, dataset_info]
   llm = get_llm()
   system_prompt = build_system_prompt()

   prompt = ChatPromptTemplate.from_messages([
       ("system", system_prompt),
       ("human", "{input}"),
       MessagesPlaceholder("agent_scratchpad"),
   ])

   agent = create_tool_calling_agent(llm, tools, prompt)
   executor = AgentExecutor(agent=agent, tools=tools, verbose=True)
   ```

3. **FastAPI app — current scaffold (WORKS, but not production-ready):**
   ```python
   from fastapi import FastAPI
   from pydantic import BaseModel

   app = FastAPI(title="CellLineFinder Agent")

   class QueryRequest(BaseModel):
       query: str

   class QueryResponse(BaseModel):
       answer: str

   @app.post("/agent/query", response_model=QueryResponse)
   async def agent_query(request: QueryRequest):
       result = await executor.ainvoke({"input": request.query})
       return QueryResponse(answer=result["output"])
   ```

4. **Run with uvicorn:**
   ```python
   if __name__ == "__main__":
       import uvicorn
       uvicorn.run(app, host="0.0.0.0", port=8000)
   ```

### Constraints
- Agent is created ONCE at module level, not per request.
- `verbose=True` during development for debugging. Set to `False` for demo.
- The `AgentExecutor` handles the full cycle: LLM decides tool → tool executes → LLM formats response.
- No manual routing logic. The LLM picks the tool based on docstrings.

> **Status:** the scaffold above is implemented and verified — the app imports,
> `POST /agent/query` is registered, and an invalid body correctly returns 422.
> Everything below upgrades it to a hostable, scalable service. Sections
> "API Design Standards" through "Hosting & Scalability" supersede this
> minimal version.

---

# PRODUCTION READINESS

## Verified gap analysis

Findings below were confirmed by running the code, not by reading it.

| # | Gap | Severity | Evidence |
|---|-----|----------|----------|
| G1 | `KNOWLEDGE_DIR` points at `knowledge/` (lowercase) but the directory is `Knowledge/` | **Blocker for hosting** | `FunctionCalling.py:15`, `Prompts.py:10` vs actual `AIAgent/Knowledge/`. Works on Windows (case-insensitive NTFS), raises `FileNotFoundError` at import on Linux/Docker |
| G2 | Network exceptions from Ensembl/HGNC propagate uncaught | **Critical** | `httpx.ReadTimeout` raised out of `_query_ensembl` during testing; status-code checks at `FunctionCalling.py:39,49,59` do not catch exceptions |
| G3 | Ensembl REST is an unreliable hard dependency | **Critical** | Sampled 8 live calls: 7× HTTP 500, 1× ReadTimeout. 0% success rate at time of review |
| G4 | 4 of 7 gene-alias tests currently FAIL because of G2/G3 | **Critical** | `pytest Tests/test_gene_alias.py` → `4 failed, 3 passed`, despite `AIAgent-testResults.md` recording all-pass |
| G5 | ENSG-ID lookups have no fallback when Ensembl is down | High | `_query_hgnc` is called with the raw ENSG string against `fetch/symbol/`, which can never match. `test_gene_alias_lookup_by_ensembl_id` returns `found=False` |
| G6 | Local `reference/gene_lookup.parquet` (19,213 genes) is unused | High | Contains `ensg_id`, `hgnc_symbol`, `approved_name`, `prev_symbols`, `alias_symbols` — the exact fields the tool returns. Lookup measured at 10 ms vs a 500 ms–10 s network round trip |
| G7 | `main.py` is empty (0 bytes) | Medium | No canonical entrypoint; app only reachable via `Schema/BusinessFlow.py` |
| G8 | No CORS middleware | **Blocker for React UI** | Verified absent. Browser will block every call |
| G9 | No `/health` endpoint | High | `GET /health` → 404. Load balancers and Docker `HEALTHCHECK` require one |
| G10 | No error handling on the route | High | Any tool or LLM failure returns a bare 500 with a stack trace |
| G11 | No API versioning | Medium | `/agent/query` is unversioned; breaking changes will break the UI |
| G12 | No auth, no rate limiting | High for public hosting | Endpoint is open; each call costs an LLM quota unit |
| G13 | No HTTP-layer tests | Medium | No `TestClient` anywhere in the repo; endpoint is the only untested component |
| G14 | New `httpx.AsyncClient` per call | Medium (perf) | `FunctionCalling.py:37,47,57` — full TLS handshake per request, no keep-alive |
| G15 | Dependencies unpinned | Medium | `requirement_agents.txt` uses `>=` ranges — builds are not reproducible |
| G16 | System prompt ~1,870 tokens resent on every LLM turn | Medium (perf/cost) | Measured 7,473 chars. An agent turn is ≥2 LLM calls → ~3.7k prompt tokens per query |
| G17 | `.env` loaded via bare `load_dotenv()` | Medium | `llm.py:11` — CWD-dependent; silently yields no API key if launched from repo root |
| G18 | No request timeout budget | Medium | A slow LLM call can hold a connection indefinitely |

Not an issue: `AIAgent/.env` **is** correctly gitignored (`.gitignore:16`), so the
Groq key is not at risk of being committed.

---

## API Design Standards

### Versioning and routes — ONE business endpoint

The API surface is deliberately minimal. All four business cases route through
a single endpoint; the LLM selects the tool, so there is no reason to expose
one route per business case. Diagnostic routes (`/ready`, `/v1/meta`) are
**removed** — Postman against the business endpoint covers that need.

| Method | Path | Purpose |
|--------|------|---------|
| POST | `/v1/agent/query` | The business endpoint. All 4 cases. |
| GET | `/health` | Liveness probe — infrastructure, not a feature |

`/health` is retained because it is not a developer convenience: Docker
`HEALTHCHECK`, and every hosting platform's load balancer, require a cheap
liveness URL to decide whether a container is alive. Without it the platform
cannot route traffic or restart a hung process. It must stay dependency-free —
if it pings Groq, an upstream blip makes the orchestrator kill healthy
containers.

The unversioned `/agent/query` stays one release as a deprecated alias with a
`Deprecation` header, then goes.

### One endpoint, two representations

Rather than adding a second `/stream` route, the same endpoint returns either
a stream or a single JSON body, chosen by the `Accept` header. This is
standard HTTP content negotiation and keeps the surface at one route.

| `Accept` header | Behaviour | Use |
|-----------------|-----------|-----|
| `text/event-stream` | SSE stream, first byte in ms | React UI |
| `application/json` (or absent) | Buffered single JSON body | Postman, curl, tests |

Postman therefore works exactly as you expect — send the POST with no special
header and get one complete JSON response.

### Request schema

```python
class QueryRequest(BaseModel):
    query: str = Field(..., min_length=1, max_length=2000)
    session_id: str | None = Field(None, max_length=64)
```

`max_length` is mandatory — without it a caller can paste megabytes of text
straight into a paid LLM context.

### Response schema

The current `{answer: str}` is too thin for a scientific UI. Scientists need
provenance, and the frontend needs to render degradation and progress states.

```python
class ToolCall(BaseModel):
    tool: str
    input: dict
    output: dict              # the FULL deterministic payload, not a summary

class QueryResponse(BaseModel):
    answer: str               # LLM narration
    data: dict | None         # deterministic tool payload — renderable alone
    request_id: str
    model: str
    tool_calls: list[ToolCall] = []
    degraded: list[str] = []  # e.g. ["ensembl_unavailable"]
    cached: bool = False
    latency_ms: int
```

Two additions carry the design:

- **`data`** — the deterministic tool output, separated from the prose. This is
  what makes millisecond responses possible (see Latency Architecture below).
- **`degraded`** — today, when Ensembl is down (as it is now), the agent
  silently answers from HGNC alone with `cross_validated=false` and nothing
  tells the user a source was missing. For a tool making scientific claims,
  silent degradation is the most dangerous failure mode.

### Error envelope

One consistent shape for every error, modelled on RFC 9457:

```python
class ErrorResponse(BaseModel):
    type: str        # "upstream_unavailable" | "rate_limited" | ...
    title: str
    status: int
    detail: str
    request_id: str
```

| Condition | Status | `type` |
|-----------|--------|--------|
| Body fails validation | 422 | `validation_error` (FastAPI default) |
| Groq 429 / quota exhausted | 503 + `Retry-After` | `rate_limited` |
| Ensembl + HGNC both unreachable | 200 with `degraded` | — (not an error; degrade) |
| Agent exceeded `max_iterations` | 504 | `agent_timeout` |
| Anything unhandled | 500, generic message | `internal_error` |

Never leak a stack trace or the Groq key in an error body. Log the traceback
server-side against the `request_id`; return only the ID to the caller.

### Middleware

1. **Request ID** — accept inbound `X-Request-ID` or generate a UUID4; attach
   to every log line and echo in the response header.
2. **CORS** — explicit origin allowlist from `ALLOWED_ORIGINS` env var. Never
   `allow_origins=["*"]` together with credentials.
3. **Timing** — record duration, emit a structured log line per request.
4. **Body size limit** — reject payloads over 64 KB at the edge.

---

## Latency Architecture

This section addresses the core requirement: **the endpoint must respond in
milliseconds, and the user must never think the system is down.**

### The constraint, stated honestly

An LLM cannot produce a complete answer in milliseconds. Measured from the
existing suite: `test_agent.py` ran 4 queries in 93.76 s, of which 60 s was
deliberate `asyncio.sleep(20)` rate-limit padding — so **~8.4 s of real
generation per query**, matching the 10–20 s noted in the test results. That
time is spent inside Groq's servers and is not something this codebase can
optimise away.

So the design does not try to make the LLM fast. It makes **the endpoint**
fast, and it keeps the user continuously informed. Those are different
problems, and conflating them is what produces a UI that looks hung.

### The key insight: every business case is deterministic data + LLM prose

All four business cases have the same shape. The *facts* come from a local
lookup; only the *explanation* needs the LLM.

| Business case | Deterministic part | Measured | LLM part |
|---------------|--------------------|----------|----------|
| BC1 Alias resolution | `gene_lookup.parquet` | **10 ms** | Phrasing the aliases |
| BC2 Score explanation | Score parquet row | **~10 ms** | Walking the 7 steps |
| BC3 Dataset info | `datasets.json` in memory | **<1 ms** | Summarising the entry |
| BC4 Hallucination guard | `found=false` | **<1 ms** | Explaining the miss |

The scientifically meaningful content — the ENSG ID, the aliases, the PIT
ranks, the core score, the confidence tier — is available in **~10 ms**. Only
the narration takes 8 s.

Therefore: **send the data as soon as it exists; stream the prose after it.**
The user sees real, correct, renderable content almost instantly, and the
explanation fills in underneath. Nothing ever looks hung, because nothing ever
is.

### Response timeline

```
t + 0 ms      event: accepted   {request_id}
                → HTTP 200 headers flushed. The endpoint has responded.
                  This is the "milliseconds" requirement, satisfied.

t + ~5 ms     event: status     {stage: "resolving", label: "Looking up TP53…"}

t + ~15 ms    event: data       {full GeneAliasResult / ScoreResult JSON}
                → UI renders the actual answer: symbol, ENSG ID, aliases,
                  scores, tier. The user already has what they came for.

t + ~20 ms    event: status     {stage: "explaining", label: "Writing explanation…"}

t + 0.5–8 s   event: token      {text: "…"}    ← streamed, word by word
                → visible motion the entire time

t + ~8 s      event: done       {latency_ms, degraded, cached}
```

Time to first byte: **milliseconds**. Time to usable content: **~15 ms**. Time
to complete prose: unchanged at ~8 s, but now it is *visibly progressing*
rather than a blank screen.

### Implementation

```python
@router.post("/v1/agent/query")
async def agent_query(req: QueryRequest, request: Request):
    if "text/event-stream" in request.headers.get("accept", ""):
        return EventSourceResponse(_stream(req), ping=10)
    return await _buffered(req)          # Postman / curl / tests
```

The generator drives `executor.astream_events(..., version="v2")` and maps
LangChain events onto the wire protocol:

| LangChain event | Emitted as |
|-----------------|------------|
| `on_tool_start` | `status` — "Looking up TP53…" |
| `on_tool_end` | **`data`** — the full deterministic payload |
| `on_chat_model_stream` | `token` |
| `on_chain_end` | `done` |

`ping=10` sends an SSE comment heartbeat every 10 s. This is not optional:
nginx's default `proxy_read_timeout` is 60 s and most managed platforms kill
idle connections around the same mark. Without heartbeats a slow LLM response
gets the connection severed mid-answer — which the user would correctly read
as "the system went down".

### Guaranteed-fast paths

Three mechanisms make repeat and simple queries genuinely millisecond-scale
end to end, not just at first byte.

**1. Response cache.** Key on the normalised query string; store the complete
`QueryResponse`. A demo asks about TP53 repeatedly, so hit rates are high.

```python
_CACHE = TTLCache(maxsize=1024, ttl=3600)
```

A cache hit returns the full answer, LLM prose included, in **<5 ms** with
`cached: true`. This is the only way to get a *complete* answer in
milliseconds, and it is worth doing before the demo — prewarm the cache with
the queries you intend to show.

**2. Prewarm at startup.** The first request otherwise pays TLS handshake to
Groq plus parquet load. Do both in the lifespan hook: load the gene index and
fire one throwaway LLM call. Startup gets ~2 s slower; the first real user
stops seeing a 12 s outlier.

**3. Skip the agent loop when intent is unambiguous.** The agent currently
makes **two** LLM round trips — one to choose the tool, one to narrate. The
tool-selection call adds seconds while producing nothing the user sees. For
inputs that are unmistakable (a bare `ENSG…` ID, an `ACH-…` model ID, a known
dataset name), match with a regex, call the tool directly, and use the LLM
only for narration. That halves latency on the most common queries. Ambiguous
input still goes through the full agent.

### Trimming the LLM time itself

| Change | Effect |
|--------|--------|
| Drop Section C from the system prompt | ~700 tokens/turn saved. It duplicates what `dataset_info` already returns |
| `max_tokens=512` (already set) with "max 150 words" in the prompt | Generation time is roughly linear in output length |
| `max_iterations=2` | Bounds the worst case to one tool call + one narration |
| Fast-path routing (above) | Removes one full round trip |

Prompt caching is not available on Groq's free tier; do not design around it.

### When the LLM fails or times out — still answer

The most important failure rule in this design. If Groq 429s, times out, or
errors **after** the deterministic data was already sent, the request does
**not** become an error. It completes with the data and an explicit note:

```
event: data       {…full ScoreResult…}
event: degraded   {reason: "llm_unavailable"}
event: done       {answer: null, degraded: ["llm_unavailable"]}
```

The UI renders the scores with a banner: *"Explanation unavailable — showing
raw results."* The scientist still gets every number they asked for. A quota
outage degrades the experience instead of destroying it.

For the buffered (Postman) path the same rule applies: HTTP 200 with
`answer: null`, `data` populated, `degraded: ["llm_unavailable"]`. Only a
failure that prevents even the deterministic lookup returns 5xx.

### What the frontend must do

The endpoint's speed is wasted if the UI waits for `done` before rendering.
The contract:

1. Render on `data`, not on `done`.
2. Show the `status` label as a live progress line — naming the actual step
   ("Looking up TP53…") is what tells the user the system is working. A
   generic spinner does not.
3. Append `token` events incrementally.
4. Treat `degraded` as a banner, never as an error page.
5. Client timeout 60 s, but never show a "failed" state while events are still
   arriving.

### Latency budget

| Metric | Target | Basis |
|--------|--------|-------|
| Time to first byte | **< 50 ms** | Hard requirement |
| Time to `data` event | **< 100 ms** | Local parquet + in-memory dicts |
| Cached full response | **< 10 ms** | TTLCache hit |
| Gap between any two events | **< 10 s** | SSE heartbeat interval |
| Full narration (p95) | ~10 s | Groq-bound; not optimisable here |
| `/health` | **< 5 ms** | No dependency checks |

Log and alert on time-to-first-byte and time-to-`data`. Those are the numbers
this design controls. Total latency is a property of Groq and should be
tracked but not treated as a regression when it moves.

---

## Reliability

### R1 — Make the local parquet the primary gene source

This is the single highest-value change, addressing G3, G4, G5 and G6 at once.

`reference/gene_lookup.parquet` already holds 19,213 genes with exactly the
fields `gene_alias_lookup` returns. Measured: 719 ms to load once at startup,
10 ms per lookup (sub-millisecond once indexed into a dict).

Revised resolution order:

```
1. Local gene_lookup index      → ~0 ms, always available, offline
2. If miss → Ensembl + HGNC     → enrichment for genes outside the panel
3. If both fail → found=false, degraded=["ensembl","hgnc"]
```

Load into two dicts at startup, keyed by `ensg_id` and by upper-cased
`hgnc_symbol`, with alias and previous symbols folded in as secondary keys.
Memory cost is a few MB.

This also fixes G5: ENSG lookups resolve locally instead of depending on a
single flaky endpoint with no fallback.

### R2 — Wrap every outbound call

```python
_RETRYABLE = (httpx.TimeoutException, httpx.ConnectError, httpx.ReadError)

async def _safe_get(client, url, headers, *, attempts=3) -> dict | None:
    for i in range(attempts):
        try:
            r = await client.get(url, headers=headers)
            if r.status_code == 200:
                return r.json()
            if r.status_code in (429, 500, 502, 503, 504):
                await asyncio.sleep(0.25 * 2**i + random.uniform(0, 0.1))
                continue
            return None                      # 404 etc. — do not retry
        except _RETRYABLE:
            await asyncio.sleep(0.25 * 2**i + random.uniform(0, 0.1))
    return None
```

Exponential backoff with jitter, capped at 3 attempts. Returns `None` rather
than raising, so a dead upstream degrades instead of 500-ing (fixes G2).

### R3 — Circuit breaker on Ensembl

Given the measured 0% success rate, retrying every request wastes ~3 s per
call. After 5 consecutive failures, open the circuit for 60 s and skip Ensembl
entirely, recording `degraded=["ensembl"]`. Half-open with a single probe.

### R4 — Correct the Ensembl request header

`FunctionCalling.py:36,46` send `Content-Type: application/json` on a GET with
no body. The correct header for a GET is `Accept: application/json`. Harmless
today, but it is not what the Ensembl REST docs specify.

### R5 — Bound the agent loop

`AgentExecutor(max_iterations=4, max_execution_time=45)`. Without a bound, a
model that loops on tool calls burns quota and holds the connection until the
client gives up.

### R6 — Timeout budget

```
Client (React)        60 s
  └─ API request      50 s   (middleware-enforced)
      └─ Agent        45 s   (max_execution_time)
          └─ LLM call 30 s
          └─ HTTP     10 s per attempt, 3 attempts
```

Each layer's timeout must be strictly shorter than its caller's, or the outer
layer gives up while inner work continues and leaks resources.

### R7 — Restore test integrity

G4 means the recorded results in `AIAgent-testResults.md` no longer reflect
reality. Split the suite:

- **Unit tests** — mock `httpx`, run offline, must always pass. These gate CI.
- **Integration tests** — hit live Ensembl/HGNC, marked `@pytest.mark.live`,
  excluded from CI, allowed to fail when an upstream is down.

A dissertation cannot claim "all tests passed" from a suite whose result
depends on a third party's uptime on the day.

---

## Performance

Measured baseline: **10–20 s per query** (from `AIAgent-testResults.md`),
dominated by LLM round trips.

| Fix | Saving | Notes |
|-----|--------|-------|
| Local parquet lookup (R1) | 0.5–10 s | Removes network entirely from the common path |
| Shared `httpx.AsyncClient` | 200–400 ms/call | Connection pooling + keep-alive; fixes G14 |
| TTL cache on gene lookups | ~100% on repeats | Gene aliases change on release cycles, not hourly. `cachetools.TTLCache(maxsize=4096, ttl=86400)` |
| Trim system prompt | ~700 tokens/turn | Section C duplicates what `dataset_info` returns. Drop it and let the tool do its job (G16) |
| Fast-path routing | One full LLM round trip | Skips tool-selection call for unambiguous input — see Latency Architecture |
| Response cache | Full answer in <5 ms | Repeat queries bypass the LLM entirely |
| SSE + `data`-first ordering | TTFB → <50 ms | The endpoint-latency requirement; see Latency Architecture |

### Shared HTTP client via lifespan

```python
@asynccontextmanager
async def lifespan(app: FastAPI):
    app.state.http = httpx.AsyncClient(
        timeout=httpx.Timeout(10.0, connect=3.0),
        limits=httpx.Limits(max_keepalive_connections=20, max_connections=100),
        headers={"User-Agent": "CellLineFinder/1.0 (dissertation; UoB)"},
    )
    app.state.genes = load_gene_index()   # parquet → dicts, once
    yield
    await app.state.http.aclose()
```

A `User-Agent` is required courtesy for Ensembl and HGNC — anonymous clients
are the first to be throttled.

### Streaming

Streaming is **not** a separate endpoint. `POST /v1/agent/query` negotiates on
the `Accept` header and streams the same work when the client asks for
`text/event-stream`. The event protocol, ordering guarantees, heartbeat
interval, and failure semantics are specified in **Latency Architecture**
above, which is the authoritative section for anything response-time related.

---

## Hosting & Scalability

### H1 — Fix the case-sensitivity blocker first (G1)

`AIAgent/Knowledge/` is capitalised on disk; the code asks for `knowledge/`.
Linux will not forgive this. Either rename the directory to `knowledge/` or
correct both constants. **Nothing else in this section matters until this is
fixed — the container will not start.**

### H2 — Stateless by construction

The service holds no per-user state: the agent, prompt, and gene index are
immutable and built at startup. Any replica can serve any request, so
horizontal scaling is just "run more containers". Preserve this — if
conversation memory is added later, put it in Redis, never in process memory.

### H3 — Container

```dockerfile
FROM python:3.12-slim
WORKDIR /app
ENV PYTHONUNBUFFERED=1 PIP_NO_CACHE_DIR=1
COPY requirements.lock.txt .
RUN pip install --no-cache-dir -r requirements.lock.txt
COPY AIAgent/ ./AIAgent/
COPY reference/gene_lookup.parquet ./reference/
RUN useradd -m app && chown -R app /app
USER app
EXPOSE 8000
HEALTHCHECK --interval=30s --timeout=3s CMD python -c "import urllib.request;urllib.request.urlopen('http://localhost:8000/health')"
CMD ["uvicorn", "AIAgent.api.app:app", "--host", "0.0.0.0", "--port", "8000", "--workers", "2"]
```

Pin to Python 3.12, not 3.14 — `langchain-core` emits a Pydantic v1
incompatibility warning on 3.14 in the current environment.

Never bake `.env` into the image. Inject `GROQ_API_KEY` as a runtime secret.

### H4 — The real scaling ceiling is the LLM, not the app

This is the most important scalability fact in the design. FastAPI on one
container handles thousands of concurrent requests; **Groq's free tier allows
~30 requests/minute**. Adding replicas multiplies pressure on a shared quota
and simply converts 200s into 429s.

Consequences:

1. **The quota is a global resource, so guard it globally.** A per-process
   semaphore does not work across replicas. Use a Redis token bucket, or run
   exactly one replica until a paid tier is in place.
2. **Fail fast and honestly.** When the bucket is empty, return 503 with
   `Retry-After` immediately rather than queueing — a queued request will
   exceed the client timeout anyway.
3. **Scale the LLM before scaling the API.** Ordered by effort: paid Groq tier
   → self-hosted Ollama/vLLM → SageMaker endpoint. The `get_llm()` provider
   switch already makes this a config change (`llm.py:19-39`), which is the
   design's strongest scalability property. Keep it that way.

### H5 — Concurrency model

Every hot path must be genuinely async — one blocking call stalls the whole
event loop and destroys throughput. `score_explainer` and `dataset_info` are
sync `@tool`s; that is fine while they read in-memory dicts, but the moment
`score_explainer` reads parquet from disk or S3 it must move to `async def`
with the read in a thread (`asyncio.to_thread`).

Workers: `2 × CPU cores`, capped at 4. This workload is I/O-bound (waiting on
Groq), so more workers mainly add memory, and each one holds its own copy of
the gene index.

### H6 — Observability

- **Structured JSON logs**, one line per request: `request_id`, route, status,
  `latency_ms`, `tool_calls`, `degraded`, token counts.
- Turn `verbose=True` **off** (`BusinessFlow.py:37`). It prints raw LLM output
  to stdout on every request — noisy, slow, and it leaks query content into
  container logs. Gate it behind `DEBUG`.
- Track four metrics: request rate, p95 latency, error rate, and 429 rate from
  Groq. The last one is the leading indicator of the H4 ceiling.

### H7 — Configuration

Replace bare `load_dotenv()` (G17) with `pydantic-settings`:

```python
class Settings(BaseSettings):
    model_config = SettingsConfigDict(env_file=".env", extra="ignore")
    groq_api_key: SecretStr
    llm_provider: str = "groq"
    llm_model: str = "openai/gpt-oss-120b"
    allowed_origins: list[str] = ["http://localhost:5173"]
    gene_lookup_path: Path = Path("reference/gene_lookup.parquet")
    debug: bool = False
```

Settings are validated at startup, so a missing key fails the container
immediately and visibly instead of surfacing as a confusing 500 on the first
user query.

### H8 — Deployment topology

```
React UI (static, CDN)
      │ HTTPS
      ▼
Reverse proxy / TLS + edge rate limit   (nginx, Caddy, or platform LB)
      │
      ▼
FastAPI containers ×N  (stateless, /health liveness probe)
      │
      ├── in-process: gene index, datasets.json, math_reference.md
      ├── Redis: global LLM token bucket + response cache   [when N > 1]
      └── Groq API  ← the actual bottleneck (H4)
```

For a dissertation demo, N=1 on Render/Railway/Fly.io with a managed TLS
certificate is sufficient and keeps the Groq quota coherent without Redis.

### H9 — Proxy configuration for SSE (easy to get wrong)

A reverse proxy will silently undo the entire latency design unless told not
to. nginx buffers responses by default, so it holds every SSE event until the
stream closes — the client then receives everything at once after 8 s, exactly
the behaviour the design exists to prevent, and it will look like a code bug.

```nginx
location /v1/agent/query {
    proxy_pass http://app:8000;
    proxy_buffering off;          # REQUIRED — else events arrive in one batch
    proxy_cache off;
    proxy_read_timeout 120s;      # must exceed the longest LLM response
    proxy_set_header Connection '';
    proxy_http_version 1.1;
    chunked_transfer_encoding on;
}
```

Also set `X-Accel-Buffering: no` on the SSE response header — managed
platforms that front you with nginx honour it without needing config access.

**Verify after every deployment**, because this is invisible locally:

```bash
curl -N -H "Accept: text/event-stream" \
     -H "Content-Type: application/json" \
     -d '{"query":"What are the aliases for TP53?"}' \
     https://<host>/v1/agent/query
```

Events must appear progressively. If they all land at once, buffering is on.

---

## Revised implementation order

Sequenced so that each step leaves the service in a working state.

| Step | Change | Addresses | Effort |
|------|--------|-----------|--------|
| 1 | Fix `Knowledge/` casing | G1 | 5 min |
| 2 | Add CORS + `/health` | G8, G9 | 20 min |
| 3 | Wrap outbound calls in `_safe_get`, add retries | G2, G3 | 1 h |
| 4 | Load `gene_lookup.parquet` as primary source | G3–G6 | 2 h |
| 5 | Split unit vs live tests; re-record results | G4, R7 | 1 h |
| **6** | **SSE via `Accept` negotiation, `data` before `token`** | **Latency** | **2 h** |
| **7** | **LLM-failure degradation (200 + `data` + `degraded`)** | **Latency** | **45 min** |
| **8** | **Response cache + startup prewarm** | **Latency** | **45 min** |
| 9 | Error envelope + global exception handlers | G10 | 1 h |
| 10 | Move to `/v1`, enrich response schema with `data` | G11 | 1 h |
| 11 | Shared `httpx` client + gene index via lifespan | G14 | 1 h |
| 12 | `test_api.py` + `test_latency.py` | G13 | 1.5 h |
| 13 | `pydantic-settings`, populate `main.py` | G7, G17 | 45 min |
| 14 | Pin dependencies, write Dockerfile | G15 | 1 h |
| 15 | Fast-path routing for unambiguous input | Latency | 1.5 h |
| 16 | API key auth + rate limiting | G12 | 1.5 h |

**Steps 1–5** make the service start on Linux and survive a dead upstream.
**Steps 6–8** deliver the millisecond-endpoint requirement and are no longer
optional polish — they are the difference between a demo that looks broken and
one that looks fast. Everything from step 9 is hardening.

Step 15 is deliberately late: it is a real latency win, but it adds a routing
path that must be kept in sync with the tools, so it should land only once the
streaming behaviour is stable and measured.

---

## `knowledge/datasets.json`

Already created. Contains 7 entries: depmap, hpa, geo, cellosaurus, ensembl, cosmic, hgnc.
Each entry has: name, full_name, url, what_it_measures, data_types, cell_line_count, version_used, caveats, reference.

## `knowledge/math_reference.md`

Create this file with the scoring pipeline equations from MATH_REFERENCE.md.
The agent's system prompt embeds this verbatim so the LLM can reference exact equations when explaining scores.

Minimum content needed:
- PIT formula: F̂(x) = rank(x) / (n+1)
- Weight formula: w_raw = 1 / (1 + mean|ρ|), then normalise
- Core score: Σ(w_i × PIT_i)
- TSG inversion rule
- Driver gating rule
- Confidence tier thresholds

---

## Testing plan

### `tests/test_gene_alias.py`
- Test `_query_ensembl("TP53")` returns valid JSON with id field
- Test `_query_hgnc("TP53")` returns valid JSON with symbol field
- Test `_query_ensembl("NONEXISTENT_GENE")` returns None
- Test `gene_alias_lookup("TP53")` returns cross-validated result
- Test `gene_alias_lookup("ENSG00000141510")` works with Ensembl ID input
- Test `gene_alias_lookup("ERBB2")` returns HER2 as synonym

### `tests/test_score_explainer.py`
- Test with mock parquet data: known gene-cell line pair returns all fields
- Test with unknown pair: returns `found=False`
- Test Pydantic validation: all fields have correct types

### `tests/test_dataset_info.py`
- Test `dataset_info("depmap")` returns valid description
- Test `dataset_info("DepMap")` case-insensitive
- Test `dataset_info("unknown")` returns available keys
- Test all 7 datasets return valid JSON

### `tests/test_agent.py`
- Test full agent flow: "What are the aliases for TP53?" → calls gene_alias_lookup
- Test full agent flow: "What is DepMap?" → calls dataset_info
- Test full agent flow: "Explain the score for ENSG00000141510 in ACH-000001" → calls score_explainer
- Test agent does NOT hallucinate when tool returns no results

### `tests/test_api.py` (NEW — closes gap G13)

The HTTP layer is currently the only untested component. Using
`fastapi.testclient.TestClient` with the agent executor mocked:

- `POST /v1/agent/query` with a valid body → 200, response matches `QueryResponse`
- `POST /v1/agent/query` with `{}` → 422 (already verified against the scaffold)
- `POST /v1/agent/query` with a 5,000-character query → 422 (`max_length`)
- `GET /health` → 200 without touching Groq or Ensembl
- Simulated Groq 429 → 503 with a `Retry-After` header
- Simulated tool exception → 500 with the error envelope, **no stack trace in the body**
- `X-Request-ID` sent by the client is echoed back unchanged
- CORS preflight `OPTIONS` from an allowed origin → correct `Access-Control-Allow-Origin`

### `tests/test_latency.py` (NEW — enforces the latency budget)

These are the tests that keep the millisecond guarantee from silently
regressing. All run offline with the LLM stubbed to a slow fake, so they
measure *this codebase's* overhead, never Groq's.

- **Time to first byte < 50 ms** — stub the LLM with a 5 s delay; assert the
  first SSE byte still arrives in under 50 ms. This is the headline assertion:
  it proves endpoint latency is decoupled from LLM latency.
- **`data` event arrives before any `token` event** — ordering is the whole
  design; a regression here silently destroys the UX.
- **Time to `data` event < 100 ms** with the same slow LLM stub.
- **Heartbeat present** — with a 25 s LLM stub, assert no gap between events
  exceeds 10 s.
- **LLM failure still yields data** — stub Groq to raise; assert 200,
  `data` populated, `answer` null, `degraded == ["llm_unavailable"]`.
- **Cache hit < 10 ms** and returns `cached: true` on the second identical query.
- **`GET /health` < 5 ms** and makes zero outbound calls (assert with a
  patched `httpx` client that records requests).
- **Fast-path routing** — a bare `ENSG00000141510` triggers exactly one LLM
  call, not two.

### Test split (closes gap G4/R7)

```ini
# pytest.ini
markers =
    live: hits third-party APIs; excluded from CI
```

```bash
pytest Tests/ -m "not live"    # offline, deterministic — must always pass
pytest Tests/ -m live          # live upstreams — may fail when Ensembl is down
```

Gene-alias unit tests mock `httpx` so they exercise the parsing and
cross-validation logic without a network dependency. The live variants stay
for manual verification but never gate a result claim in the write-up.

### Run tests
```bash
pytest Tests/ -v -m "not live"
```

---

## Implementation order (original scaffold — Steps 1–8 COMPLETE)

```
Step 1: llm.py + .env                         → DONE  verify Groq connection
Step 2: FunctionCalling.py (dataset_info only) → DONE  simplest tool, test standalone
Step 3: Prompts.py                             → DONE  build system prompt
Step 4: BusinessFlow.py                        → DONE  wire agent, test with dataset_info
Step 5: FunctionCalling.py (gene_alias_lookup) → DONE  add API calls, test standalone
Step 6: FunctionCalling.py (score_explainer)   → PARTIAL: still mock data, no parquet read
Step 7: Test full agent with all 3 tools       → DONE  (4/4 agent tests passed)
Step 8: FastAPI endpoint                       → DONE  POST /agent/query verified
Step 9: Connect to React UI                    → BLOCKED on CORS (gap G8)
```

Step 9 cannot proceed until the CORS middleware is added. Continue from the
**Revised implementation order** above, which sequences the production work.

---

## Environment setup

```bash
pip install langchain langchain-groq langchain-core httpx pydantic pydantic-settings python-dotenv fastapi uvicorn pytest
```

Additional packages for the production API:

```bash
pip install "uvicorn[standard]" sse-starlette cachetools tenacity slowapi structlog pyarrow pandas
```

| Package | Purpose |
|---------|---------|
| `uvicorn[standard]` | `uvloop` + `httptools` — measurably faster event loop |
| `sse-starlette` | `EventSourceResponse` with built-in heartbeat (`ping=`) |
| `cachetools` | `TTLCache` for gene lookups |
| `tenacity` | Declarative retry/backoff (alternative to hand-rolled `_safe_get`) |
| `slowapi` | Per-IP rate limiting at the route level |
| `structlog` | Structured JSON logging with request-ID binding |
| `pyarrow`, `pandas` | Reading `gene_lookup.parquet` and the score parquet |

**Pin everything before hosting (G15).** Generate a lockfile so the deployed
image matches what was tested:

```bash
pip freeze > requirements.lock.txt
```

Target **Python 3.12** in the container. The current dev environment is 3.14,
where `langchain-core` raises a Pydantic v1 incompatibility warning.

## Key decisions documented

| Decision | Choice | Reason |
|----------|--------|--------|
| Framework | LangChain single agent | Production-grade, LLM-agnostic, `bind_tools()` native |
| LLM provider | Groq (free tier) | Free, fast, reliable tool calling, 30 req/min |
| Model | Qwen 2.5 32B | Best tool-calling accuracy at this size |
| Gene lookup | REST API (not web search) | Deterministic, sub-second, zero hallucination |
| Score explanation | Pre-computed values + system prompt | LLM explains, never computes |
| Dataset descriptions | Static JSON in system prompt | Corpus is <2000 tokens, no retrieval needed |
| API framework | FastAPI | Async, Pydantic-native, production-standard |
| Gene source priority | Local parquet first, REST as enrichment | Ensembl measured at 0% success during review; local table has the same 19,213 genes at 10 ms |
| Upstream failure mode | Degrade with a `degraded[]` flag, never 500 | A missing source must be visible to the scientist, not silent |
| API versioning | `/v1` prefix from the start | Frontend contract stability; retrofitting versioning is painful |
| API surface | One business endpoint + `/health` | 4 business cases share one route; LLM does the routing. No diagnostic endpoints — Postman covers that |
| Streaming | Content negotiation on one route, not a `/stream` route | Keeps the surface minimal; standard HTTP |
| Health check | `/health` only, dependency-free | Required by Docker/LB. Checking Groq would let an upstream blip kill healthy containers |
| Response ordering | Deterministic `data` before LLM `token`s | Makes the answer usable in ~15 ms while prose streams for ~8 s |
| LLM failure | 200 with `data` + `degraded`, not 5xx | Scientist still gets every number; quota outage degrades rather than destroys |
| Scaling limit | LLM quota, not app throughput | Groq free tier ~30 req/min is the binding constraint; replicas share it |
| State | Stateless; Redis only if N > 1 | Keeps horizontal scaling trivial |
| Config | `pydantic-settings`, secrets injected at runtime | Fails fast at startup; `.env` never in the image |