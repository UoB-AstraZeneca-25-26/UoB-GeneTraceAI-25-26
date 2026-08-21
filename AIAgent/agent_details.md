# CellLineFinder AI Agent

## What it does

CellLineFinder's scoring pipeline ranks cancer cell lines for a given gene. The pipeline produces numbers — PIT ranks, correlation-penalised weights, core scores, confidence tiers — that are correct but not self-explanatory. Gene identifiers (symbols, Ensembl IDs, DepMap model IDs) are a well-known source of confusion in this domain.

The agent is a narration and lookup layer over deterministic data. It serves three purposes: resolving gene aliases, explaining how a cell line was scored, and describing the data sources behind the pipeline. In every case the agent retrieves structured data first, then uses an LLM to translate that data into plain language. The LLM explains — it never computes, never guesses, and never answers from its own training knowledge.

## Architecture

```
React UI
    │
    ▼
POST /v1/agent/query   ◄── one endpoint, all three business cases
    │
    ▼
┌─────────────────────────────────────────────────────────────┐
│               LangChain AgentExecutor                       │
│                                                             │
│   ┌─────────────────────┐    ┌────────────────────────┐     │
│   │   System prompt     │    │   LLM backend          │     │
│   │   Math reference    │    │   gpt-oss-120b (Groq)  │     │
│   │   (7 pipeline steps)│    │   Swappable via .env   │     │
│   └─────────────────────┘    └────────────────────────┘     │
│                                                             │
│                     .bind_tools()                           │
└──────────┬──────────────────┬──────────────────┬────────────┘
           │                  │                  │
           ▼                  ▼                  ▼
  ┌────────────────┐  ┌───────────────┐  ┌──────────────┐
  │gene_alias_lookup│ │score_explainer│  │ dataset_info │
  └───────┬────────┘  └──────┬────────┘  └──────┬───────┘
          │                  │                   │
          ▼                  ▼                   ▼
  ┌────────────────┐  ┌───────────────┐  ┌──────────────┐
  │ gene_lookup    │  │ Pipeline      │  │ datasets.json│
  │ .parquet       │  │ output        │  │ (7 sources)  │
  │ (19,213 genes) │  │ (parquet)     │  │              │
  │ + Ensembl/HGNC │  │               │  │              │
  │   REST (backup)│  │               │  │              │
  └────────────────┘  └───────────────┘  └──────────────┘
```

The agent uses one endpoint for everything. The scientist sends a natural-language query, the LLM picks the right tool from its docstrings, the tool returns structured JSON from a local data source, and the LLM formats that JSON into a readable explanation. The scientist sees the structured data immediately (~15 ms), then the prose explanation streams in over ~8 seconds.

## How to use it

The API exposes a single endpoint. All three business cases are invoked the same way — only the query string changes.

```
POST /v1/agent/query
Content-Type: application/json
Body: { "query": "..." }
```

### Business Case 1 — Gene alias lookup

A scientist types a gene name or Ensembl ID and needs to know its official symbol, synonyms, and whether the identifier is valid.

**When it triggers:** the search bar in the React UI, on submit.

**Example request:**
```javascript
const result = await fetch("http://localhost:8000/v1/agent/query", {
  method: "POST",
  headers: { "Content-Type": "application/json" },
  body: JSON.stringify({ query: "What are the aliases for TP53" })
});
const json = await result.json();
```

**Example response (json.data):**
```json
{
  "found": true,
  "query": "TP53",
  "symbol": "TP53",
  "ensembl_id": "ENSG00000141510",
  "full_name": "tumor protein p53",
  "previous_symbols": [],
  "synonyms": ["p53", "LFS1"],
  "cross_validated": true,
  "sources": ["ensembl", "hgnc"]
}
```

**Example response (json.answer):**
> TP53 (Ensembl ENSG00000141510) is also known as p53 and LFS1. These are the primary synonyms listed in Ensembl and HGNC for the tumor-protein p53 gene. No former symbols are recorded.

**How it resolves the gene:** the tool queries Ensembl and HGNC REST APIs concurrently first; a hit on either (or both) is returned immediately, cross-validated when both agree. Only if both network sources are unreachable does it fall back to a local parquet table (19,213 genes, ~3 µs lookup) reconciled from the same two sources upstream, so a fallback hit is still treated as cross-validated.

### Business Case 2 — Score explanation

A scientist sees a ranked cell line and wants to understand why it was ranked at that position — what the numbers mean and how they were computed.

**When it triggers:** the user clicks a "Why this rank?" button on a specific cell line row in the results table. The frontend already has the `ensg_id` and `model_id` from the ranking data.

**Example request:**
```javascript
const ensgId = "ENSG00000141510";  // from the current gene context
const modelId = "ACH-000001";      // from the row they clicked

const result = await fetch("http://localhost:8000/v1/agent/query", {
  method: "POST",
  headers: { "Content-Type": "application/json" },
  body: JSON.stringify({
    query: `Explain the score for ${ensgId} in ${modelId}`
  })
});
const json = await result.json();
```

**Two modes.** The tool reads the scoring intermediates from a separate scoring API, configured by `SCORING_API_URL`. No numbers are ever generated inside the agent.

- **API connected** (`SCORING_API_URL` set): the numeric fields carry the real pipeline values for that (gene, cell line) pair, and the answer walks through the 7 steps using them.
- **API not connected** (default): every numeric field is `null` and the answer is a methodology explanation of the 7 steps — no values are invented for the requested cell line.

**Example response (json.data) — methodology mode (default):**
```json
{
  "found": true,
  "ensg_id": "ENSG00000141510",
  "model_id": "ACH-000001",
  "raw_expression": null,
  "raw_proteomics": null,
  "pit_expression": null,
  "pit_proteomics": null,
  "rho_ep": null,
  "weight_expression": null,
  "weight_proteomics": null,
  "core_score": null,
  "confidence_tier": null,
  "is_tsg": true,
  "driver_gated": false,
  "message": "Methodology explanation — scoring API not connected. Explaining the 7-step pipeline methodology for this gene."
}
```

`is_tsg` still comes from the gene's role in `reference/gene_lookup.parquet` (TP53 is annotated `both`), so the explanation covers the TSG inversion that applies to this gene. `found` is `false` only when the ENSG id is outside the 19,213-gene panel.

With the API connected the same fields are populated from its response, and `message` becomes `"Scoring intermediates retrieved."`. Any field the API omits stays `null` rather than being filled in.

**Example response (json.answer) — methodology mode:**
> Ranking runs in seven steps. Raw expression (log2 TPM+1) and proteomics (log2 ratio) sit on incompatible scales, so each is converted to a PIT percentile — hypothetically, a pit_expression of 0.72 would mean the line expresses this gene more highly than 72% of the ~950-line panel. Expression and proteomics are correlated, so the layers are weighted to avoid double-counting redundant evidence, and the core score is their weighted mean. TP53 is a tumour suppressor, so the score is inverted (1 − core score): low abundance is the informative signal. Lines carrying a driver alteration are then placed ahead of non-driver lines, and a confidence tier records how much evidence was available. No scoring API is connected, so these are illustrative figures, not results for ACH-000001.

**How the explanation works:** the system prompt contains the full mathematical reference for the 7-step scoring pipeline (PIT normalisation, correlation-penalised weights, core score, TSG inversion, driver gating, confidence tiers). When the tool returns numbers, the LLM walks through each step using them and never computes anything independently — every number traces back to the tool output. When the tool returns nulls, it explains the method instead, labelling any illustrative figure as hypothetical.

### Business Case 3 — Dataset information

A scientist wants to know what a data source measures, where it comes from, and what caveats apply.

**When it triggers:** the user clicks an info icon next to a dataset name (e.g. "DepMap", "HPA") in the UI.

**Example request:**
```javascript
const datasetName = "depmap";  // from the UI element

const result = await fetch("http://localhost:8000/v1/agent/query", {
  method: "POST",
  headers: { "Content-Type": "application/json" },
  body: JSON.stringify({ query: `What is ${datasetName}` })
});
const json = await result.json();
```

**Example response (json.data):**
```json
{
  "name": "depmap",
  "full_name": "Cancer Dependency Map (DepMap)",
  "url": "https://depmap.org/portal/",
  "what_it_measures": "Genome-wide CRISPR knockout gene-effect scores (Chronos)...",
  "data_types": ["CRISPR gene effect (Chronos)", "RNAi", "expression (TPM)", "mutation", "copy number"],
  "cell_line_count": "~1100",
  "caveats": "Cell line panel is skewed toward well-characterised lines...",
  "reference": "Tsherniak et al., Cell 2017; Dempster et al., Nat Commun 2019"
}
```

**Example response (json.answer):**
> DepMap (the Cancer Dependency Map) is a public resource that profiles ~1,100 cancer cell lines. It provides genome-wide CRISPR knockout (Chronos) and RNAi gene-effect scores, plus matched multi-omics data such as RNA-seq expression (TPM), mutation calls, and copy-number profiles. Because the gene-effect scores are model-derived, they reflect predicted dependencies rather than direct biochemical measurements.

**How it retrieves the description:** the tool looks up the dataset name in a static JSON file containing 7 curated entries (depmap, hpa, geo, cellosaurus, ensembl, cosmic, hgnc). These descriptions were extracted once from AstraZeneca's project documentation and the primary papers. No web search or RAG is involved.

### Business Case 4 — Hallucination guard

When any of the three tools receives an identifier that does not exist, the agent says so explicitly rather than inventing data.

**Example request:**
```javascript
const result = await fetch("http://localhost:8000/v1/agent/query", {
  method: "POST",
  headers: { "Content-Type": "application/json" },
  body: JSON.stringify({
    query: "Explain the score for ENSG_FAKE in ACH-FAKE"
  })
});
```

**Example response (json.data):**
```json
{
  "found": false,
  "ensg_id": "ENSG_FAKE",
  "model_id": "ACH-FAKE",
  "raw_expression": null,
  "pit_expression": null,
  "core_score": null,
  "confidence_tier": null,
  "message": "No scoring data found for ENSG_FAKE in ACH-FAKE."
}
```

**Example response (json.answer):**
> The scoring pipeline couldn't find any data for the requested gene-cell-line pair. Neither expression nor proteomics measurements are available, so none of the intermediate values could be computed. If you need a score, please check that the gene identifier is correct and that the cell line is part of the CellLineFinder dataset.

This is not a separate feature — it is the required behaviour of the other three tools when a lookup misses. For a scientific tool, "says it doesn't know" is a deliverable, not an error.

## How hallucination prevention works

Three architectural layers ensure the agent never fabricates scientific claims.

**Layer 1 — Tool-first enforcement.** The system prompt instructs the LLM: "Always use a tool before answering. Never guess." Every factual claim must pass through a tool call against local data. The LLM cannot answer from its training knowledge.

**Layer 2 — Structured null responses.** When a tool cannot find data, it returns a fully structured Pydantic model with `found: false` and every numeric field set to `null`. The LLM sees explicit nulls — it cannot fabricate a `core_score` of 0.7 when the tool returned `null`.

**Layer 3 — Closed-world data sources.** Each tool reads from a specific, bounded data source: a parquet table (19,213 genes), pre-computed pipeline output, or a 7-entry JSON file. The LLM has no access to external knowledge, web search, or its own parametric memory for factual claims. Its only role is translating structured JSON into prose.

## Frontend integration

All three business cases use the same function. The query string is the only thing that changes.

```javascript
async function queryAgent(queryString) {
  const res = await fetch("http://localhost:8000/v1/agent/query", {
    method: "POST",
    headers: { "Content-Type": "application/json" },
    body: JSON.stringify({ query: queryString })
  });
  return await res.json();
}

// BC1 — search bar submit
queryAgent(`What are the aliases for ${geneInput}`)

// BC2 — "Why this rank?" button click
queryAgent(`Explain the score for ${ensgId} in ${modelId}`)

// BC3 — dataset info icon click
queryAgent(`What is ${datasetName}`)
```

The response always contains `answer` (prose explanation), `data` (structured tool output), `degraded` (any missing sources), and `cached` (whether this was a cache hit). The frontend should render `data` immediately and append the `answer` text as it arrives.

BC2 is the exception: for `score_explainer` the agent answers with the structured 7-step contract instead of prose, so the response carries `methodology` and `answer` is null. The frontend renders those steps as the interactive `MethodologyChart` pipeline. If the model returns prose anyway, `methodology` is null and `answer` holds the text — the chart is an enhancement, not a hard dependency.

For SSE streaming (progressive rendering), add the `Accept: text/event-stream` header. Events arrive in this order: `accepted` → `status` → `data` (render here, ~15 ms) → `status` → `token` × N (prose streams ~8 s) → `done`.

## Response structure

Every response from `/v1/agent/query` has the same shape:

```json
{
  "answer": "TP53 (ENSG00000141510) is also known as...",
  "data": { "found": true, "symbol": "TP53", ... },
  "tool_calls": [{ "name": "gene_alias_lookup", "args": {"query": "TP53"} }],
  "degraded": [],
  "methodology": null,
  "cached": false,
  "latency_ms": 8432,
  "request_id": "f484ab09-0b4c-44e0-9abb-c5cacfffcb58"
}
```

| Field | Type | What it contains |
|-------|------|-----------------|
| `answer` | string or null | LLM-generated prose explanation. Null if the LLM was unavailable (data is still present), or if the answer was the structured `methodology` contract |
| `data` | object | Structured output from whichever tool was called. Contains the actual scientific content |
| `methodology` | object or null | `score_explainer` only: `{"steps": [{key, value, formula, status}, ×7]}`. `status` is `active`, `skipped`, or `inverted`. Null for every other tool, and for a `score_explainer` answer that came back as prose |
| `tool_calls` | array | Which tool was called and with what arguments. For audit and debugging |
| `degraded` | array | Names of sources that were unavailable (e.g. `["ensembl"]`, `["llm_unavailable"]`). Empty when everything worked |
| `cached` | boolean | True if this was an exact repeat of a recent query, served from cache in <10 ms |
| `latency_ms` | integer | Total request time in milliseconds |
| `request_id` | string | UUID for tracing this request through logs |

When `degraded` contains `"llm_unavailable"`, the `answer` field is null but `data` still contains every number the scientist asked for. The LLM being down degrades the experience (no prose explanation) but does not destroy it (the structured data is still delivered).

## Failure behaviour

| Situation | HTTP status | What the scientist sees |
|-----------|-------------|------------------------|
| Normal operation | 200 | `data` + `answer` |
| BC2 (score_explainer) normal operation | 200 | `data` + `methodology` (7 chart steps), `answer` null |
| BC2 but the model answered in prose | 200 | `data` + `answer`, `methodology` null — the UI falls back to the paragraph |
| LLM rate-limited but tool succeeded | 200 | `data` present, `answer` null, `degraded: ["llm_unavailable"]` |
| Ensembl/HGNC unreachable | 200 | Answer from local table, `degraded: ["ensembl"]` |
| Invalid input format | 422 | Error message explaining what was wrong |
| LLM rate-limited and no tool output | 503 | Error with `Retry-After` header |

The design principle: an LLM outage after the deterministic data was obtained is not an error. The scientist still receives every number. Only a failure that prevents even the local lookup produces a 5xx.

## Latency design

An LLM takes ~8 seconds to generate a complete response. The agent cannot change that. What it can do is make the scientifically meaningful content available immediately and stream the explanation afterwards.

| What happens | When |
|-------------|------|
| Structured data available | ~15 ms (local lookup) |
| Prose explanation starts streaming | ~500 ms |
| Full explanation complete | ~8 s |
| Cached repeat query | <10 ms |

The frontend should render the structured data on arrival and fill in the explanation underneath as tokens stream in.

## Test results

44 tests total: 35 run offline (deterministic, no network dependency), 9 require live Ensembl/HGNC/Groq connections.

Offline suite verified 2026-08-17: **35 passed in 34.5 seconds.**

| Test area | Tests | What it covers |
|-----------|-------|----------------|
| Gene alias resolution | 9 | Local index hits, cross-validation logic, live API calls |
| API endpoints | 8 | Success/error paths, request ID echo, CORS, health check |
| Latency budget | 8 | Time-to-first-byte <50 ms, data-before-tokens ordering, cache <10 ms |
| Dataset info | 5 | All 7 datasets return valid JSON, case-insensitive matching |
| Agent (end-to-end) | 4 | All four business cases with real LLM calls |
| Score explainer | 3 | Known pair, unknown pair, Pydantic field validation |

## Technology stack

| Component | Choice | Why |
|-----------|--------|-----|
| Agent framework | LangChain (single AgentExecutor) | LLM-agnostic, native tool binding, production-tested |
| LLM | OpenAI gpt-oss-120b via Groq | 120B MoE model, strong tool-calling, free tier sufficient for demo |
| API framework | FastAPI | Async, Pydantic-native, SSE streaming support |
| Gene data | Ensembl/HGNC REST, primary + local parquet (19,213 genes) as fallback | Live APIs give current data; local table covers queries when both APIs are down |
| Dataset descriptions | Static JSON (7 entries) | Extracted once from project documentation, no runtime retrieval needed |
| Scoring pipeline reference | math_reference.md embedded in system prompt | The 7-step pipeline specification the LLM uses to explain scores |

The LLM backend is swappable via a single environment variable. The production code supports Groq; stubs exist for OpenAI and HuggingFace. Switching to a self-hosted model (Ollama on AWS EC2) requires changing one line in `.env`.

## Known gaps

| Gap | Impact | Status |
|-----|--------|--------|
| Response cache is per-process | Fine for single-replica demo; needs Redis for multi-replica deployment | Acceptable for current scope |
| Groq free tier rate limit | 8K tokens per minute; rapid sequential queries can hit the limit | Mitigated by max_tokens=512 and response caching |