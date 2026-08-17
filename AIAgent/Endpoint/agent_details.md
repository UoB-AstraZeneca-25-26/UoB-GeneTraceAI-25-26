# CellLineFinder AI Agent — Business & Technical Flows

**Scope:** the `AIAgent/` subtree of the GeneTraceAI project.
**Status:** implemented and verified as of 2026-08-17. Offline suite `pytest Tests/ -m "not live"` → **35 passed, 9 deselected**.

This document describes what the agent does for the scientist (business flows),
how a request actually travels through the code (technical flows), and how the
whole thing is wired into runnable infrastructure. `Endpoint/plan.md` is the
design rationale and gap analysis; this document is the description of the
system as built.

---

## 1. What the agent is for

CellLineFinder's scoring pipeline ranks cancer cell lines for a given gene. It
emits numbers — PIT ranks, correlation-penalised weights, a core score, a
confidence tier. Those numbers are correct but not self-explanatory, and the
identifiers involved (gene symbols, ENSG IDs, `ACH-…` model IDs) are a
well-known source of confusion in this domain.

The agent is a **narration and lookup layer over deterministic data**. It never
computes a score, never resolves a gene by guessing, and never answers from the
model's own parametric knowledge. Every user-visible fact originates in a tool
call against a local table or a curated knowledge file; the LLM's only job is to
turn that JSON into prose a scientist can read.

That single constraint — *the LLM explains, it never computes* — is what the
rest of the architecture is built around, and it is what makes the latency
design (§4) possible at all.

### The four business cases

| # | Business case | User question, typically | Tool | Ground truth |
|---|---------------|--------------------------|------|--------------|
| BC1 | Alias resolution | "What are the aliases for TP53?" / "What is ENSG00000141510?" | `gene_alias_lookup` | `reference/gene_lookup.parquet` (19,213 genes), Ensembl + HGNC REST as enrichment |
| BC2 | Score explanation | "Why is ACH-000001 ranked here for TP53?" | `score_explainer` | Pipeline scoring output (currently mock — see §9) |
| BC3 | Data provenance | "What is DepMap?" / "Where does the expression data come from?" | `dataset_info` | `Knowledge/datasets.json` (7 curated sources) |
| BC4 | Hallucination guard | Any of the above with an identifier that does not exist | whichever tool matched | `found: false` + an explicit message |

BC4 is not a separate feature — it is the *required behaviour* of the other
three when the lookup misses. It is called out as a business case because for a
scientific tool, "says it doesn't know" is a deliverable, not an error path.

---

## 2. Business flow, end to end

```
Scientist (React UI / Postman / curl)
   │  "What are the aliases for TP53?"
   ▼
POST /v1/agent/query
   │
   ├─ 1. Is this an exact repeat of a recent query?  ──── yes ──► return cached
   │                                                              answer (<10 ms)
   ├─ 2. Is the intent unmistakable?                 ──── yes ──► call the tool
   │     (bare ENSG id / ENSG+ACH pair /                          directly, skip
   │      a known dataset name)                                   tool-selection
   │                                                              LLM round trip
   └─ 3. Otherwise ─────────────────────────────────► full AgentExecutor:
                                                        LLM picks the tool from
                                                        the tool docstrings
   ▼
Tool executes against local data  ─────────────► DETERMINISTIC RESULT (~ms)
   │                                              emitted to the client first
   ▼
LLM narrates that JSON in plain language ──────► PROSE (~8 s, streamed)
   │
   ├─ LLM unavailable / rate-limited?  ─────────► still HTTP 200:
   │                                              data present, answer null,
   │                                              degraded: ["llm_unavailable"]
   └─ upstream source missing?         ─────────► degraded: ["ensembl"] etc.
   ▼
Response: {answer, data, tool_calls, degraded, cached, latency_ms, request_id}
```

Two properties of this flow are the design, not incidental:

**The facts arrive before the prose.** The scientifically meaningful content —
symbol, ENSG ID, aliases, PIT ranks, core score, confidence tier — is available
in milliseconds because it comes from an in-memory dict. The narration takes
seconds because it comes from an LLM. The response is ordered accordingly, so
the UI can render a real answer almost immediately and fill in the explanation
underneath.

**Missing sources are visible, never silent.** If Ensembl is unreachable, the
answer is still returned, but `degraded` names the source that was missing. A
tool that makes scientific claims must not quietly narrow its evidence base
without saying so.

### BC1 — Alias resolution, in detail

Resolution order is local-first:

```
1. Local gene index      dict lookup, ~3 µs, offline, always available
2. Ensembl + HGNC        only on a local miss; concurrent, retried, circuit-broken
3. found = false         with degraded: ["ensembl", "hgnc"] if both were down
```

The local index (`AgentDevelopment/gene_index.py`) loads
`reference/gene_lookup.parquet` once at startup — measured **222 ms to load,
19,213 ENSG keys, 61,894 symbol keys** (approved symbols plus previous symbols
and aliases folded in as secondary keys), then **~3 µs per lookup**. A local hit
is marked `cross_validated: true, sources: ["local_index"]`, because the
reference table was itself built by reconciling Ensembl and HGNC upstream in the
pipeline's lookup-builder stage.

Ensembl and HGNC are therefore *enrichment for genes outside the panel*, not the
critical path. This matters: during the plan.md review Ensembl REST was measured
at 0% success over 8 sampled calls. The service is unaffected by that for any
gene in the panel.

### BC2 — Score explanation, in detail

`score_explainer(ensg_id, model_id)` returns every intermediate the pipeline
produced: `raw_expression`, `raw_proteomics`, `pit_expression`,
`pit_proteomics`, `rho_ep`, `weight_expression`, `weight_proteomics`,
`core_score`, `confidence_tier`, `is_tsg`, `driver_gated`.

The system prompt embeds `Knowledge/math_reference.md` verbatim, so the LLM has
the exact definitions to hand — PIT normalisation, the correlation-penalised
weight formula, the core-score sum, TSG inversion, driver gating, and the four
confidence tiers — and is instructed to walk the steps *using the numbers the
tool returned*. Verified in the recorded run: the model reproduced
`0.55 × 0.72 + 0.45 × 0.68 = 0.702` and the TSG inversion `1 − 0.702 = 0.298`
correctly, with no invented values.

### BC3 — Data provenance, in detail

`Knowledge/datasets.json` holds seven curated entries — `depmap`, `hpa`, `geo`,
`cellosaurus`, `ensembl`, `cosmic`, `hgnc` — each with `full_name`, `url`,
`what_it_measures`, `data_types`, `cell_line_count`, `version_used`, `caveats`,
`reference`. Matching is exact key first, then substring, then a listing of the
available keys.

Notably, these descriptions are **not** embedded in the system prompt. An
earlier design did embed them (~700 tokens resent on every LLM turn); that was
removed because `dataset_info` already returns the full entry on demand. The
knowledge base is queried, not memorised.

### BC4 — Hallucination guard, in detail

Every tool returns a well-formed Pydantic model with `found: false` and a
human-readable `message` rather than raising or returning nothing. The system
prompt's rules ("Always use a tool before answering. Never guess.") plus a
structured miss give the model something concrete to narrate. The recorded
behaviour for `ENSG_FAKE` / `ACH-FAKE`: acknowledged the miss, listed which
inputs were unavailable, suggested verifying the identifier, and did **not**
fall back to a different tool or invent a score.

---

## 3. Technical flow — the request lifecycle

### 3.1 Middleware chain

Middleware is added in `api/app.py:create_app()`. Starlette applies middleware
in reverse registration order, so the effective inbound order is:

```
inbound  ──► RequestIDMiddleware      accept X-Request-ID or mint a UUID4;
             │                        attach to request.state
             ▼
             TimingMiddleware         perf counter start
             │
             ▼
             BodySizeLimitMiddleware  reject > 64 KB with a 413 envelope
             │
             ▼
             CORSMiddleware           explicit origin allowlist, credential-free
             │
             ▼
             SlowAPIMiddleware        per-IP rate limit (default 30/minute)
             │
             ▼
             route: require_api_key ──► handler
outbound ◄── X-Request-ID, X-Response-Time-Ms headers + one structured log line
```

The request ID is the thread that ties everything together: it is on every log
line, echoed in the response header, present in the response body, and is the
*only* internal detail an error response ever exposes.

### 3.2 Content negotiation — one route, two representations

`POST /v1/agent/query` inspects the `Accept` header:

| `Accept` | Path | Consumer |
|----------|------|----------|
| contains `text/event-stream` | `_stream()` → `EventSourceResponse` | React UI |
| anything else (or absent) | `_buffered()` → single `JSONResponse` | Postman, curl, tests |

There is deliberately no separate `/stream` route. Both representations are
driven by **one** internal async generator, `_agent_events()`, so the fast-path
logic, degradation rules, and tool handling exist in exactly one place and
cannot drift between the two surfaces.

Route inventory is intentionally minimal:

| Method | Path | Notes |
|--------|------|-------|
| POST | `/v1/agent/query` | the business endpoint; all four business cases |
| POST | `/agent/query` | deprecated alias; sets `Deprecation: true` and a `Link: rel="successor-version"` header |
| GET | `/health` | liveness only; dependency-free by design |

`/health` performs zero outbound calls and is asserted to do so in
`Tests/test_api.py` and `Tests/test_latency.py`. If it pinged Groq, an upstream
blip would make the orchestrator kill healthy containers.

### 3.3 `_agent_events()` — the shared engine

```
_agent_events(query)
   │
   ├── detect_fast_path(query, dataset_names)
   │      ^ENSG\d{11}(\.\d+)?$              ──► gene_alias_lookup
   │      ^ENSG\d{11} (in|for)? ACH-\d+$    ──► score_explainer
   │      query is a known dataset key      ──► dataset_info
   │
   ├── FAST PATH (one LLM call — narration only)
   │      status  → tool.ainvoke(args) → data → status → llm.astream(...) → token*
   │      tool itself raised?  ──► error   (no data at all)
   │      LLM raised?          ──► error   (data already sent → degrades, not fails)
   │
   └── AGENT PATH (two LLM calls — selection, then narration)
          executor.astream_events(version="v2")
             on_tool_start        ──► status  "Calling gene_alias_lookup…"
             on_tool_end          ──► data    (full deterministic payload) + status
             on_chat_model_stream ──► token
          GroqRateLimitError      ──► error(UpstreamRateLimitedError)
          any other exception     ──► error
```

The fast path exists because the agent path spends a full LLM round trip
deciding which tool to call — seconds of latency producing nothing the user
sees. For inputs where intent is unmistakable, that round trip is pure waste.
Ambiguous natural-language input still goes through the full `AgentExecutor`; no
routing heuristic is applied to anything it cannot match with certainty.

`Tests/test_latency.py::test_fast_path_triggers_exactly_one_llm_call` asserts
both halves of this: exactly one LLM call, and the `AgentExecutor` is never
invoked (the stubbed executor raises `AssertionError` if touched).

### 3.4 SSE event protocol

```
t + 0 ms      event: accepted   {request_id}
                → 200 headers flushed; the endpoint has responded

t + ~5 ms     event: status     {stage: "resolving", label: "Looking up gene alias lookup…"}

t + ~15 ms    event: data       {full GeneAliasResult / ScoreResult JSON}
                → the UI renders the actual answer here

t + ~20 ms    event: status     {stage: "explaining", label: "Writing explanation…"}

t + 0.5–8 s   event: token      {text: "…"}     ← streamed incrementally

t + ~8 s      event: done       {full QueryResponse}
```

`ping=10` sends an SSE heartbeat every 10 s, and `X-Accel-Buffering: no` is set
on the response. Both are load-bearing: nginx's default `proxy_read_timeout` is
60 s and buffers responses by default, so without these a slow LLM answer either
gets its connection severed mid-stream or is delivered as one batch at the end —
which is exactly the "looks hung" behaviour the design exists to prevent.

The frontend contract that makes this worth anything: **render on `data`, not on
`done`**; show the `status` label as a live progress line; append `token` events
incrementally; treat `degraded` as a banner, never an error page.

### 3.5 Buffered path

Same generator, accumulated instead of streamed. `_buffered()`:

1. Normalises the query (`" ".join(query.strip().lower().split())`) and checks
   the TTL response cache. A hit returns the complete answer — prose included —
   with `cached: true`, in **under 10 ms** (asserted in `test_latency.py`).
2. Otherwise drains `_agent_events()`, collecting `data` events into
   `tool_calls` and `token` events into the answer.
3. Resolves the terminal state (§3.6).
4. Caches the response only if no error occurred, and returns HTTP 200.

### 3.6 Failure semantics — the most important rule

```
terminal_error AND no tool_calls   ──► raise ──► 503 (Groq 429) or 500, error envelope
terminal_error AND tool_calls      ──► 200: data populated,
                                            answer = null,
                                            degraded += "llm_unavailable"
no terminal_error                  ──► 200: answer + data, cached for next time
```

An LLM outage **after** the deterministic data was obtained is not an error. The
scientist still receives every number they asked for; only the prose is missing,
and its absence is declared explicitly. A quota outage degrades the experience
instead of destroying it. Only a failure that prevents even the local lookup
produces a 5xx.

`degraded` is assembled from two sources: flags the tool itself reported (e.g.
`["ensembl"]` when that upstream was skipped or failed), plus `llm_unavailable`
appended by the route.

### 3.7 Error envelope

Every error — validation, auth, rate limit, timeout, unhandled — has one shape
(RFC 9457-flavoured), so a caller needs exactly one parser:

```json
{"type": "...", "title": "...", "status": 0, "detail": "...", "request_id": "..."}
```

| Condition | Status | `type` |
|-----------|--------|--------|
| Body fails validation | 422 | `validation_error` |
| Missing/invalid `X-API-Key` (when auth enabled) | 401 | `unauthorized` |
| Body > 64 KB | 413 | `payload_too_large` |
| Per-IP rate limit exceeded | 429 + `Retry-After` | `client_rate_limited` |
| Groq 429 / quota exhausted, no data obtained | 503 + `Retry-After` | `rate_limited` |
| Agent exceeded its execution budget | 504 | `agent_timeout` |
| Anything unhandled | 500, generic detail | `internal_error` |

Stack traces and secrets never cross the wire. The unhandled handler calls
`logger.exception(...)` with the request ID and returns a fixed generic string.
`Tests/test_api.py::test_tool_exception_returns_500_with_no_stack_trace` asserts
that an exception message containing `"boom: secret internal detail"` does not
appear anywhere in the response body.

---

## 4. Latency architecture

The honest constraint: an LLM cannot produce a complete answer in milliseconds.
Measured from the original agent suite — 4 queries in 93.76 s, of which 60 s was
deliberate rate-limit padding — real generation is **~8.4 s per query**, spent
inside Groq's servers and not optimisable from this codebase.

So the design does not try to make the LLM fast. It makes **the endpoint** fast
and keeps the user continuously informed. Four mechanisms:

| Mechanism | Effect | Where |
|-----------|--------|-------|
| Local-first gene resolution | removes 0.5–10 s of network from the common path | `gene_index.py`, `FunctionCalling.gene_alias_lookup` |
| `data` before `token` ordering | usable content at ~15 ms while prose streams for ~8 s | `routes._stream_events` |
| Response cache (TTLCache, 1024 entries / 1 h) | full repeat answers in <10 ms, `cached: true` | `api/cache.py`, `routes._buffered` |
| Fast-path routing | removes one full LLM round trip on unambiguous input | `api/cache.detect_fast_path` |
| Startup prewarm | first real user stops seeing a 12 s TLS-handshake outlier | `api/app.lifespan` |
| System-prompt trim | ~700 tokens/turn saved by not embedding `datasets.json` | `Schema/Prompts.build_system_prompt` |

### Budget, and how it is enforced

| Metric | Target | Enforcing test |
|--------|--------|----------------|
| Time to first byte | < 50 ms | `test_time_to_first_byte_under_50ms` (LLM stubbed with a 5 s delay) |
| Time to `data` event | < 100 ms | `test_time_to_data_event_under_100ms` |
| `data` strictly before any `token` | ordering | `test_data_arrives_before_any_token` |
| Cached full response | < 10 ms, `cached: true` | `test_cache_hit_under_10ms_and_marked_cached` |
| SSE heartbeat interval | 10 s (< nginx's 60 s) | `test_sse_ping_interval_matches_heartbeat_budget` |
| `GET /health` | fast, zero outbound calls | `test_health_is_fast_and_makes_zero_outbound_calls` |
| LLM failure still yields data | 200 + `degraded` | `test_llm_failure_still_returns_data` |
| Full narration (p95) | ~10 s | Groq-bound; tracked, not gated |

Every latency test stubs the LLM, so the suite measures *this codebase's*
overhead and never Groq's. The three SSE-timing tests iterate `_stream_events`
directly rather than going through `TestClient`: httpx's `ASGITransport` collects
the entire ASGI response before returning, so it would report ~5000 ms for every
SSE response regardless of how fast events are actually flushed. `_stream_events`
was factored out of `_stream` specifically to make that measurable.

---

## 5. Reliability

| Mechanism | Implementation | Rationale |
|-----------|----------------|-----------|
| Local primary source | `gene_index.load_gene_index` | Ensembl measured at 0% success during review; the local table has the same fields for all 19,213 panel genes |
| Never-raising HTTP | `http_utils.safe_get` — 3 attempts, exponential backoff `0.25 × 2ⁿ` + jitter, retries on `{429,500,502,503,504}` and on timeout/connect/read errors, returns `None` otherwise | A dead upstream must degrade, not 500 |
| Circuit breaker | `http_utils.CircuitBreaker(failure_threshold=5, reset_after=60)` on Ensembl, half-open with a single probe | Retrying a confirmed-dead upstream wastes ~3 s per call |
| Correct REST headers | `Accept: application/json` on GETs (not `Content-Type`) | What the Ensembl REST docs specify |
| Bounded agent loop | `AgentExecutor(max_iterations=4, max_execution_time=45)` | A model that loops on tool calls otherwise burns quota and holds the connection |
| Graceful missing data | `load_gene_index` returns an *empty* index if the parquet is absent | A checkout without reference data still starts; lookups fall through to network |
| Missing-source visibility | `degraded[]` on both the tool result and the response | Silent degradation is the dangerous failure mode for a scientific tool |

### Timeout budget

Each layer is strictly shorter than its caller, or the outer layer gives up
while inner work continues and leaks resources:

```
Client (React)        60 s
  └─ API request      50 s    api_request_timeout
      └─ Agent        45 s    AgentExecutor.max_execution_time
          └─ LLM call 30 s    llm_call_timeout → ChatGroq(timeout=…)
          └─ HTTP     10 s per attempt (3 s connect), 3 attempts
```

---

## 6. Component reference

```
AIAgent/
├── main.py                       canonical entrypoint: python main.py
├── Dockerfile                    python:3.12-slim, non-root, HEALTHCHECK
├── pytest.ini                    asyncio_mode=auto; `live` marker
├── requirement_agents.txt        loose ranges, for dev installs
├── requirements.lock.txt         pinned versions, for the image
│
├── api/                          ── the production HTTP surface ──
│   ├── app.py                    app factory + lifespan (shared client, gene index, prewarm)
│   ├── routes.py                 /v1/agent/query, deprecated alias, /health, _agent_events
│   ├── schemas.py                QueryRequest / QueryResponse / ToolCall / HealthResponse
│   ├── errors.py                 typed exceptions + the single error envelope
│   ├── middleware.py             request ID, timing/structlog, body cap, CORS
│   ├── security.py               X-API-Key dependency + slowapi per-IP limiter
│   ├── settings.py               pydantic-settings; validated once at import
│   ├── cache.py                  TTLCache factory, query normalisation, fast-path regexes
│   └── _pathsetup.py             puts AgentDevelopment/ and Schema/ on sys.path
│
├── AgentDevelopment/             ── tools and their plumbing ──
│   ├── FunctionCalling.py        the 3 @tool definitions + Pydantic I/O schemas
│   ├── gene_index.py             parquet → in-memory GeneIndex (by ENSG, by symbol)
│   ├── http_utils.py             pooled AsyncClient, safe_get, CircuitBreaker
│   └── llm.py                    get_llm() — provider switch via env vars
│
├── Schema/                       ── agent assembly ──
│   ├── BusinessFlow.py           tools + llm + prompt → AgentExecutor (module-level, once)
│   └── Prompts.py                build_system_prompt(): role rules + math_reference.md
│
├── Knowledge/                    ── static, curated ground truth ──
│   ├── datasets.json             7 data sources, extracted once from AZ documentation
│   └── math_reference.md         the 7 scoring steps, embedded verbatim in the prompt
│
├── Tests/                        44 tests: 35 offline + 9 live
└── Endpoint/
    ├── plan.md                   design rationale, gap analysis, implementation order
    └── agent_details.md          this document
```

### Notes on a few of these

**`api/_pathsetup.py`** — `AgentDevelopment/` and `Schema/` are plain script
folders, not installable packages, so `import FunctionCalling` needs them on
`sys.path` first. Importing `_pathsetup` *is* the side effect; it mirrors what
`Tests/conftest.py` does for the test runtime. The `# noqa: F401` / `# noqa: E402`
comments in `app.py` and `routes.py` mark this deliberate import ordering.

**`Schema/BusinessFlow.py`** — the system prompt is passed as a
`SystemMessage(content=…)` rather than a `("system", text)` template tuple,
because the embedded content contains literal `{ }` braces that
`ChatPromptTemplate`'s f-string parser would try to interpret. `verbose` is
gated behind the `DEBUG` env var: left on, `AgentExecutor` prints raw LLM output
— including user query content — straight to stdout on every request.

**`AgentDevelopment/llm.py`** — `get_llm()` reads env fresh on each call and
returns a `BaseChatModel`. Only the `groq` branch is implemented;
`huggingface` and `openai` raise `NotImplementedError` as declared stubs. This
one function is the single point of change for swapping the LLM backend, which
is the design's strongest scalability property (§7.4).

**`Tests/conftest.py`** — an autouse fixture drops the module-level
`http_utils._client` after every test. The process-wide singleton client is
correct in production (one event loop for the app's lifetime) but unsafe across
pytest-asyncio's per-test loops, where a client holding keep-alive connections
from a closed loop raises `RuntimeError: Event loop is closed`. The fixture also
reconfigures stdout/stderr to UTF-8, because Windows' cp1252 default cannot
encode the Unicode punctuation Groq models emit.

---

## 7. Infrastructure

### 7.1 Configuration

All config flows through `api/settings.py` (`pydantic-settings`), cached with
`@lru_cache` so env is parsed exactly once per process. `.env` is located
absolutely, relative to the package — not via a bare `load_dotenv()`, which is
CWD-dependent and silently yields no API key when launched from the repo root.

Settings are validated at import, so a missing `GROQ_API_KEY` fails the
container immediately and visibly rather than surfacing as a confusing 500 on
the first user query.

| Group | Keys | Defaults |
|-------|------|----------|
| LLM | `groq_api_key` (`SecretStr`, required), `llm_provider`, `llm_model`, `llm_temperature`, `llm_max_tokens` | groq / `openai/gpt-oss-120b` / 0.0 / 512 |
| CORS | `allowed_origins` (CSV-coercing validator) | `["http://localhost:5173"]` |
| Data | `gene_lookup_path`, `knowledge_dir` | `../reference/gene_lookup.parquet`, `AIAgent/Knowledge` |
| Budget | `agent_max_iterations`, `agent_max_execution_time`, `api_request_timeout`, `llm_call_timeout`, `http_timeout`, `http_connect_timeout` | 4 / 45 / 50 / 30 / 10 / 3 s |
| Cache | `response_cache_ttl`, `response_cache_maxsize`, `gene_cache_ttl`, `gene_cache_maxsize` | 3600 s / 1024, 86400 s / 4096 |
| Auth | `api_keys` (empty ⇒ auth disabled), `rate_limit` | `[]`, `30/minute` |
| Misc | `max_body_bytes`, `debug`, `host`, `port` | 64 KB, false, `0.0.0.0`, 8000 |

`groq_api_key` is a `SecretStr`, so it cannot be accidentally serialised into a
log line or an error body. `AIAgent/.env` is gitignored.

### 7.2 Process lifecycle

```python
@asynccontextmanager
async def lifespan(app):
    client = httpx.AsyncClient(timeout=…, limits=Limits(20 keepalive / 100 max),
                               headers={"User-Agent": "CellLineFinder/1.0 (dissertation; UoB)"})
    http_utils.set_http_client(client)          # one pooled client, process-wide
    index = load_gene_index(settings.gene_lookup_path)
    FunctionCalling.set_gene_index(index)       # 19,213 genes, ~222 ms, once
    try:    await llm.ainvoke("ping")           # prewarm: pay the TLS handshake now
    except: log warning and continue            # best-effort — never crash startup
    yield
    await http_utils.close_http_client()
```

Three things happen exactly once per process: the pooled HTTP client (a new
`AsyncClient` per call meant a full TLS handshake every time), the gene index,
and a throwaway LLM call to warm the Groq connection. The prewarm is wrapped in
a bare `except` on purpose — an upstream blip at boot must not stop the
container coming up, because `/health` has to stay reachable for the
orchestrator regardless.

A `User-Agent` is set as required courtesy for Ensembl and HGNC; anonymous
clients are the first to be throttled.

### 7.3 Container

```dockerfile
FROM python:3.12-slim
COPY AIAgent/requirements.lock.txt .
RUN pip install --no-cache-dir -r requirements.lock.txt
COPY AIAgent/ ./AIAgent/
COPY reference/gene_lookup.parquet ./reference/gene_lookup.parquet
RUN useradd -m app && chown -R app /app
USER app
WORKDIR /app/AIAgent
HEALTHCHECK --interval=30s --timeout=3s CMD python -c "…urlopen('http://localhost:8000/health')"
CMD ["uvicorn", "api.app:app", "--host", "0.0.0.0", "--port", "8000", "--workers", "2"]
```

Points worth noting:

- **Python 3.12, not 3.14.** `langchain-core` emits a Pydantic v1 incompatibility
  warning on 3.14 (visible in the local test run, which uses 3.14).
- **The sibling layout is preserved.** `reference/gene_lookup.parquet` is copied
  to `/app/reference/`, one level above `/app/AIAgent/`, matching the path the
  code computes.
- **Non-root user**, lockfile copied before the source so the pip layer caches.
- **`.env` is never baked in.** `GROQ_API_KEY` is injected as a runtime secret.
- **Directory casing matters.** `Knowledge/` is capitalised on disk, and both
  `FunctionCalling.py` and `Prompts.py` now reference it with the capital K.
  Windows' case-insensitive NTFS hides this class of bug; Linux does not, and it
  is a `FileNotFoundError` at import — the container would not start at all.

### 7.4 Deployment topology and the real scaling ceiling

```
React UI (static, CDN)
      │ HTTPS
      ▼
Reverse proxy / TLS + edge rate limit    (nginx, Caddy, or a platform LB)
      │  proxy_buffering off  ← REQUIRED for SSE
      ▼
FastAPI containers ×N   (stateless; /health liveness probe)
      │
      ├── in-process: gene index, datasets.json, math_reference.md, response cache
      ├── Redis: global LLM token bucket + shared response cache   [only when N > 1]
      └── Groq API  ← the actual bottleneck
```

**The service is stateless by construction.** The agent, prompt, and gene index
are immutable and built at startup; nothing per-user is retained. Any replica
can serve any request, so horizontal scaling is just "run more containers". The
one caveat: the response cache is currently in-process, so with N > 1 each
replica warms its own. That is acceptable at N=1 and is the first thing to move
to Redis if replicas are added — along with conversation memory, if that is ever
introduced.

**The binding constraint is not the app.** FastAPI on one container handles
thousands of concurrent requests; Groq's free tier allows roughly 30
requests/minute. Adding replicas multiplies pressure on a shared quota and
simply converts 200s into 429s. Consequences:

1. The quota is a global resource, so guard it globally. A per-process semaphore
   does not work across replicas — use a Redis token bucket, or run exactly one
   replica until a paid tier is in place.
2. Fail fast and honestly. An empty bucket returns 503 with `Retry-After`
   immediately rather than queueing; a queued request would exceed the client
   timeout anyway.
3. Scale the LLM before scaling the API. Ordered by effort: paid Groq tier →
   self-hosted Ollama/vLLM → SageMaker endpoint. `get_llm()` already makes this
   a config change rather than a code change.

Workers are set to 2. This workload is I/O-bound (waiting on Groq), so more
workers mainly add memory — each holds its own copy of the gene index. `2 × CPU
cores`, capped at 4, is the ceiling worth using.

**Concurrency caution:** every hot path must be genuinely async — one blocking
call stalls the event loop and destroys throughput. `score_explainer` and
`dataset_info` are sync `@tool`s, which is fine while they read in-memory dicts,
but the moment `score_explainer` reads parquet from disk or S3 it must become
`async def` with the read in `asyncio.to_thread`.

### 7.5 Proxy configuration for SSE

A reverse proxy silently undoes the entire latency design unless told not to.
nginx buffers responses by default, holding every SSE event until the stream
closes — the client then receives everything at once after 8 s, which looks
exactly like a code bug.

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

The app also sets `X-Accel-Buffering: no` on the SSE response, which managed
platforms fronting you with nginx honour without needing config access.

This is invisible locally, so verify after every deployment:

```bash
curl -N -H "Accept: text/event-stream" -H "Content-Type: application/json" \
     -d '{"query":"What are the aliases for TP53?"}' https://<host>/v1/agent/query
```

Events must appear progressively. If they all land at once, buffering is on.

### 7.6 Security posture

| Control | Implementation | Default |
|---------|----------------|---------|
| API key | `X-API-Key` header checked against `settings.api_keys` | **disabled** (empty list) — acceptable for a single-replica demo behind a private URL |
| Per-IP rate limit | `slowapi`, keyed on remote address | `30/minute`, 429 + `Retry-After` |
| Body size cap | `BodySizeLimitMiddleware` | 64 KB → 413 |
| Input length cap | `QueryRequest.query` `max_length=2000` | 422 |
| CORS | explicit origin allowlist, `allow_credentials=False` | never `["*"]` with credentials |
| Secret handling | `SecretStr`; `.env` gitignored and never in the image | — |
| Error leakage | generic 500 detail; traceback logged against `request_id` only | asserted by test |

The input cap is not cosmetic: without it a caller can paste megabytes of text
straight into a paid LLM context.

### 7.7 Observability

`structlog` emits one structured line per request from `TimingMiddleware` —
`request_id`, method, path, `status_code`, `latency_ms` — and the response
carries `X-Request-ID` and `X-Response-Time-Ms` headers. `QueryResponse` itself
reports `latency_ms`, `tool_calls`, `degraded`, and `cached`, so client-side
telemetry needs no separate channel.

Four metrics are worth alerting on: request rate, p95 latency, error rate, and
the 429 rate from Groq — the last being the leading indicator of the quota
ceiling in §7.4. Log and alert specifically on time-to-first-byte and
time-to-`data`; those are the numbers this codebase controls. Total latency is a
property of Groq and should be tracked but not treated as a regression when it
moves.

---

## 8. Testing

```bash
pytest Tests/ -m "not live"    # 35 tests, offline, deterministic — gates CI
pytest Tests/ -m live          #  9 tests, real Ensembl/HGNC/Groq — may fail on upstream downtime
```

**Verified 2026-08-17: `35 passed, 9 deselected` in 34.5 s.**

| File | Tests | Covers |
|------|-------|--------|
| `test_gene_alias.py` | 9 (4 live) | local-index hits, mocked `httpx` parsing/cross-validation, live Ensembl/HGNC |
| `test_api.py` | 8 | 200/422/413/503/500 paths, request-ID echo, CORS preflight, `/health` with network blocked |
| `test_latency.py` | 8 | the full latency budget in §4 |
| `test_dataset_info.py` | 5 | exact/case-insensitive/substring/miss, all 7 entries valid |
| `test_agent.py` | 4 (live) | the four business-case walkthroughs end to end |
| `test_score_explainer.py` | 3 | known pair, unknown pair, Pydantic field types |
| `test_tool.py` | 1 (live) | tool smoke test |

The split matters for the write-up. The original suite was 100% live, and a
suite whose result depends on a third party's uptime on the day cannot back a
"all tests passed" claim — during the plan.md review, 4 of 7 gene-alias tests
failed purely because Ensembl was returning 500s. Unit tests now mock `httpx` so
they exercise parsing and cross-validation logic offline; live variants remain
for manual verification but never gate a result claim.

Latest live-suite state: 8/9 passed; the one failure is a Windows
`ProactorEventLoop` `RuntimeError: Event loop is closed` during async connection
teardown *after* the assertions had already run — an environment artifact, not a
tool or logic defect.

---

## 9. Known gaps and what is not yet done

| Item | Status | Impact |
|------|--------|--------|
| `score_explainer` reads mock data | `_MOCK_SCORES` holds one hardcoded `(ENSG00000141510, ACH-000001)` row | BC2 is demonstrable but not yet wired to real pipeline output. This is the largest outstanding functional gap; the replacement is a parquet read, and §7.4's concurrency note applies when it lands |
| Response cache is per-process | in-memory `TTLCache` | Fine at N=1; move to Redis before adding replicas |
| TLS termination | not configured | Infra, not app code — supplied by the hosting platform |
| nginx SSE config (§7.5) | documented, not deployed | Must be applied and verified per §7.5 after any deployment |
| Load testing | only the synthetic budget in `test_latency.py` | No concurrency/soak testing yet; the Groq quota is the ceiling anyway |
| API-key auth | implemented, disabled by default | Enable by populating `api_keys` before any public exposure |
| Python 3.14 dev environment | `langchain-core` Pydantic v1 warning | Cosmetic locally; the container pins 3.12 |

---

## 10. Design decisions, and why

| Decision | Choice | Reason |
|----------|--------|--------|
| Agent framework | LangChain single agent | LLM-agnostic, native `bind_tools()`, no hand-rolled routing |
| Tool selection | LLM reads the tool docstrings | Docstrings are the routing contract; no manual `if/else` to keep in sync |
| Gene source priority | local parquet first, REST as enrichment | Ensembl measured at 0% success during review; the local table has the same 19,213 genes at ~3 µs |
| Score explanation | pre-computed values + math reference in the prompt | The LLM explains, never computes — the arithmetic stays deterministic and auditable |
| Dataset descriptions | tool-returned, not prompt-embedded | Saves ~700 tokens on every LLM turn; the tool already returns the full entry |
| Upstream failure mode | `degraded[]` flag, never a 500 | A missing source must be visible to the scientist, not silent |
| LLM failure mode | 200 with `data` + `degraded`, not 5xx | The scientist still gets every number; a quota outage degrades rather than destroys |
| API surface | one business endpoint + `/health` | Four business cases share one route because the LLM does the routing |
| Streaming | content negotiation on one route, not a `/stream` route | Standard HTTP; keeps the surface minimal and the logic single-sourced |
| Response ordering | deterministic `data` before LLM `token`s | Answer usable in ~15 ms while prose streams for ~8 s |
| Health check | `/health` only, dependency-free | Checking Groq would let an upstream blip kill healthy containers |
| Versioning | `/v1` prefix from the start | Frontend contract stability; retrofitting versioning is painful |
| State | stateless; Redis only if N > 1 | Keeps horizontal scaling trivial |
| Config | `pydantic-settings`, secrets injected at runtime | Fails fast at startup; `.env` never in the image |
| Scaling limit | LLM quota, not app throughput | Groq free tier ~30 req/min is the binding constraint; replicas share it |
