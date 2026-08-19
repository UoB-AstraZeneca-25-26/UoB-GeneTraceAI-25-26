# CellLineFinder AI Agent — API Documentation

**Base URL:** `API URL (get it once hosted)`

This document describes the AI agent's API contract. The frontend team
needs only this document to integrate — no knowledge of the agent's
internals is required.

---

## Quick reference

| Method | Endpoint | Purpose |
|--------|----------|---------|
| POST | `/v1/agent/query` | All agent queries (alias, score explanation, dataset info) |
| GET | `/health` | Liveness check |

One endpoint serves all three business cases. The query string determines
which tool the agent invokes internally.

---

## POST /v1/agent/query

### Request

```
POST /v1/agent/query
Content-Type: application/json
```

```json
{
  "query": "string (required, max 2000 chars)"
}
```

The `query` field is a natural-language question. The agent determines
which tool to call based on the content.

### Response modes

The endpoint supports two response modes based on the `Accept` header:

| Accept header | Mode | Use case |
|--------------|------|----------|
| `application/json` (or omitted) | Buffered — single JSON response after completion | Postman, curl, simple fetch |
| `text/event-stream` | SSE streaming — progressive events | React UI with real-time rendering |

---

## Response shape (buffered mode)

Every response has the same structure regardless of which business case
was invoked:

```json
{
  "answer": "string | null",
  "data": { },
  "methodology": null,
  "tool_calls": [
    {
      "name": "gene_alias_lookup",
      "args": { "query": "BRAF" }
    }
  ],
  "degraded": [],
  "cached": false,
  "latency_ms": 8432,
  "request_id": "f484ab09-0b4c-44e0-9abb-c5cacfffcb58"
}
```

### Field reference

| Field | Type | Description |
|-------|------|-------------|
| `answer` | `string \| null` | LLM-generated prose explanation. Null if the LLM was unavailable or if `methodology` is populated instead |
| `data` | `object \| null` | Structured output from the tool that was called. Contains the scientifically meaningful content |
| `methodology` | `object \| null` | Structured 7-step pipeline explanation (BC2 only). When present, `answer` is null — the frontend renders the steps as a chart instead of prose |
| `tool_calls` | `array` | Which tool was called and with what arguments. For debugging and audit |
| `degraded` | `array` | Names of sources that were unavailable during this request. Empty when everything worked |
| `cached` | `boolean` | True if this was a cache hit (exact repeat of a recent query, served in <10 ms) |
| `latency_ms` | `integer` | Total request time in milliseconds |
| `request_id` | `string` | UUID for tracing — also returned in the `X-Request-ID` response header |

---

## Business Case 1 — Gene Alias Lookup

### When to call

When a user searches for a gene and the frontend needs the official
symbol, Ensembl ID, full name, and known synonyms.

### Query pattern

```json
{ "query": "What are the aliases for BRAF" }
```

Or with an Ensembl ID:

```json
{ "query": "What are the aliases for ENSG00000157764" }
```

### Response — `data` field

```json
{
  "found": true,
  "query": "BRAF",
  "symbol": "BRAF",
  "ensembl_id": "ENSG00000157764",
  "full_name": "B-Raf proto-oncogene, serine/threonine kinase",
  "previous_symbols": [],
  "synonyms": ["BRAF1", "RAFB1", "B-RAF1"],
  "cross_validated": true,
  "sources": ["local_index"],
  "degraded": []
}
```

### Field reference (data)

| Field | Type | Description |
|-------|------|-------------|
| `found` | `boolean` | Whether the gene was found in any source |
| `query` | `string` | The original query string |
| `symbol` | `string \| null` | HGNC-approved gene symbol |
| `ensembl_id` | `string \| null` | Ensembl gene ID (dot-version stripped) |
| `full_name` | `string \| null` | Full descriptive gene name |
| `previous_symbols` | `string[]` | Retired/former symbols |
| `synonyms` | `string[]` | Known aliases and alternative names |
| `cross_validated` | `boolean` | True if the result was confirmed against multiple sources |
| `sources` | `string[]` | Which sources contributed: `"local_index"`, `"ensembl"`, `"hgnc"` |
| `degraded` | `string[]` | Sources that were unreachable: `"ensembl"`, `"hgnc"` |

### When `found` is false

```json
{
  "found": false,
  "query": "FAKEGENE",
  "symbol": null,
  "ensembl_id": null,
  "full_name": null,
  "previous_symbols": [],
  "synonyms": [],
  "cross_validated": false,
  "sources": [],
  "degraded": []
}
```

The `answer` field will contain a message explaining that the gene
was not found and suggesting the user verify the identifier.

### Frontend integration

```javascript
// Call after the user submits a gene search
const result = await queryAgent(`What are the aliases for ${geneInput}`);

if (result.data?.found) {
  // Display: full_name under gene header, synonyms as tags
  setAliasData(result.data);
}
```

---

## Business Case 2 — Score Explanation (Methodology)

### When to call

When a user clicks "Explain scoring methodology" on a specific
gene × cell line pair.

### Query pattern

```json
{ "query": "Explain the score for ENSG00000157764 in ACH-000620" }
```

Use the Ensembl ID (not gene symbol) and the ACH model ID. Both are
available from the scoring API response.

### Response — `methodology` field

When the agent returns structured methodology steps, the `methodology`
field is populated and `answer` is null:

```json
{
  "answer": null,
  "methodology": {
    "steps": [
      {
        "key": "Raw measurements",
        "value": "Two independent assays: expression as log2(TPM+1) and proteomics as log2 protein-to-reference ratio, measured across the cell line panel.",
        "formula": "log₂(TPM + 1), log₂(ratio)",
        "status": "active"
      },
      {
        "key": "PIT normalisation",
        "value": "Each raw value converted to a percentile rank within this gene across all ~950 cell lines, making expression and proteomics directly comparable.",
        "formula": "F̂(x) = rank(x) / (n + 1)",
        "status": "active"
      },
      {
        "key": "Correlation penalty",
        "value": "Spearman correlation between expression and proteomics layers quantifies redundancy. Higher correlation means less independent evidence per layer.",
        "formula": "n_eff = 2 / (1 + |ρ_EP|)",
        "status": "active"
      },
      {
        "key": "Core score",
        "value": "Weighted sum of PIT ranks using correlation-penalised weights. For two exchangeable layers, optimal weights are equal: 0.5 each.",
        "formula": "core = w_E × E + w_P × P",
        "status": "active"
      },
      {
        "key": "TSG inversion",
        "value": "BRAF is an oncogene — inversion applies only to tumour suppressors where low abundance signals loss of function. This step is skipped.",
        "formula": "score_TSG = 1 − core_score",
        "status": "skipped"
      },
      {
        "key": "Driver gating",
        "value": "Cell lines carrying a known driver alteration matching the gene's role are promoted above all non-driver lines regardless of numeric score.",
        "formula": "sort = driver × 10⁶ + score",
        "status": "active"
      },
      {
        "key": "Confidence tier",
        "value": "Categorical assessment based on data completeness and evidence configuration. Not a calibrated probability — an auditable evidence ordering.",
        "formula": "HIGH | MEDIUM | CONTEXT | LOW",
        "status": "active"
      }
    ]
  },
  "data": { ... }
}
```

### Step field reference

| Field | Type | Description |
|-------|------|-------------|
| `key` | `string` | Short label (3–5 words) — displayed as the chart node title |
| `value` | `string` | Explanation (20–25 words) — displayed on click/expand |
| `formula` | `string \| null` | Mathematical formula — displayed below the key in the chart node |
| `status` | `string` | Visual state: `"active"` \| `"skipped"` \| `"inverted"` |

### Status values and their meaning

| Status | Meaning | When it occurs | Visual treatment |
|--------|---------|----------------|------------------|
| `active` | This step applied normally | Most steps for most genes | Default styling (teal on hover) |
| `skipped` | This step did not apply | TSG inversion for oncogenes | Dashed border, muted colours, "skipped" label |
| `inverted` | This step applied in reverse | TSG inversion for tumour suppressors (e.g. TP53) | Amber theme, swap icon |

### Fallback — prose instead of structured steps

If the LLM returns prose instead of structured JSON (e.g. during
rate limiting or model issues), `methodology` will be null and
`answer` will contain the prose explanation:

```json
{
  "answer": "The pipeline starts with raw measurements...",
  "methodology": null,
  "data": { ... }
}
```

The frontend should handle both cases:

```javascript
if (result.methodology?.steps) {
  // Render MethodologyChart component
  setMethodology(result.methodology);
} else if (result.answer) {
  // Fallback: render prose text
  setExplanation(result.answer);
}
```

### Frontend integration

```javascript
const [methodology, setMethodology] = useState(null);
const [explanation, setExplanation] = useState(null);
const [loading, setLoading] = useState(false);

async function handleExplain(ensgId, modelId) {
  setLoading(true);
  try {
    const result = await queryAgent(
      `Explain the score for ${ensgId} in ${modelId}`
    );
    if (result.methodology?.steps) {
      setMethodology(result.methodology);
    } else {
      setExplanation(result.answer || "Explanation unavailable.");
    }
  } catch {
    setExplanation("Agent unavailable.");
  }
  setLoading(false);
}
```

---

## Business Case 3 — Dataset Information

### When to call

When a user clicks an info icon next to a data source name (e.g.
"DepMap", "HPA", "GEO") to learn what that source measures and
where the data comes from.

### Query pattern

```json
{ "query": "What is depmap" }
```

Valid dataset names: `depmap`, `hpa`, `geo`, `cellosaurus`, `ensembl`,
`cosmic`, `hgnc`. Case-insensitive. Substring matching is supported
(e.g. "dep" matches "depmap").

### Response — `data` field

```json
{
  "found": true,
  "name": "depmap",
  "full_name": "Cancer Dependency Map (DepMap)",
  "url": "https://depmap.org/portal/",
  "what_it_measures": "Genome-wide CRISPR knockout gene-effect scores (Chronos) and RNAi dependency scores across a large panel of cancer cell lines, plus matched multi-omics (expression, mutation, copy number).",
  "data_types": [
    "CRISPR gene effect (Chronos)",
    "RNAi",
    "expression (TPM)",
    "mutation",
    "copy number"
  ],
  "cell_line_count": "~1100",
  "version_used": "DepMap public release (see pipeline config for exact quarter)",
  "caveats": "Cell line panel is skewed toward well-characterised lines; dependency scores are model-derived, not direct measurements.",
  "reference": "Tsherniak et al., Cell 2017; Dempster et al., Nat Commun 2019"
}
```

### Field reference (data)

| Field | Type | Description |
|-------|------|-------------|
| `name` | `string` | Short identifier |
| `full_name` | `string` | Full display name |
| `url` | `string` | Link to the source portal |
| `what_it_measures` | `string` | Plain-language description of what the data contains |
| `data_types` | `string[]` | List of measurement types available |
| `cell_line_count` | `string \| null` | Approximate number of cell lines covered |
| `version_used` | `string` | Which version/release is used in the pipeline |
| `caveats` | `string` | Known limitations or biases |
| `reference` | `string` | Primary citation(s) |

### When the dataset is not found

```json
{
  "error": "not found",
  "available": ["cellosaurus", "cosmic", "depmap", "ensembl", "geo", "hgnc", "hpa"]
}
```

### Frontend integration

```javascript
async function handleDatasetInfo(datasetName) {
  const result = await queryAgent(`What is ${datasetName}`);

  if (result.data?.name) {
    setDatasetModal({
      title: result.data.full_name,
      url: result.data.url,
      description: result.answer,       // prose explanation from LLM
      caveats: result.data.caveats,
      reference: result.data.reference,
      dataTypes: result.data.data_types,
    });
  }
}
```

---

## Business Case 4 — Hallucination Guard

Not a separate endpoint — this is the automatic behaviour when any
of the three tools receives an identifier that does not exist.

### Example: unknown gene

```json
{ "query": "What are the aliases for FAKEGENE123" }
```

Response:

```json
{
  "answer": "The gene 'FAKEGENE123' was not found in the CellLineFinder panel...",
  "data": {
    "found": false,
    "query": "FAKEGENE123",
    "symbol": null,
    "ensembl_id": null,
    "synonyms": [],
    "message": "Gene not found."
  }
}
```

### Example: unknown gene × cell line pair

```json
{ "query": "Explain the score for ENSG_FAKE in ACH-FAKE" }
```

Response:

```json
{
  "answer": "The scoring pipeline couldn't find any data for this pair...",
  "data": {
    "found": false,
    "ensg_id": "ENSG_FAKE",
    "model_id": "ACH-FAKE",
    "raw_expression": null,
    "core_score": null,
    "message": "No scoring data found for ENSG_FAKE in ACH-FAKE."
  },
  "methodology": null
}
```

### Frontend handling

Always check the `found` field before rendering:

```javascript
if (result.data?.found === false) {
  showErrorBanner(result.data.message || result.answer);
}
```

---

## SSE Streaming Mode

For real-time rendering in the React UI, set the `Accept` header to
`text/event-stream`:

```javascript
const response = await fetch("/v1/agent/query", {
  method: "POST",
  headers: {
    "Content-Type": "application/json",
    "Accept": "text/event-stream"
  },
  body: JSON.stringify({ query: "What are the aliases for BRAF" })
});
```

### Event sequence

```
event: accepted    → { request_id }
event: status      → { stage: "resolving", label: "Looking up gene aliases..." }
event: data        → { full structured tool output }           ← RENDER HERE (~15 ms)
event: status      → { stage: "explaining", label: "Writing explanation..." }
event: token       → { text: "BRAF" }                         ← append incrementally
event: token       → { text: " (ENSG00000157764)" }
event: token       → { text: " is also known as..." }
...
event: done        → { full QueryResponse }
```

### Frontend SSE handling

```javascript
const eventSource = new EventSource(url);  // or use fetch + ReadableStream

// Render the answer as soon as 'data' arrives (~15ms)
eventSource.addEventListener("data", (e) => {
  const payload = JSON.parse(e.data);
  setStructuredData(payload);    // render immediately
});

// Append prose tokens incrementally (~0.5–8s)
eventSource.addEventListener("token", (e) => {
  const payload = JSON.parse(e.data);
  setAnswer(prev => prev + payload.text);
});

// Final complete response
eventSource.addEventListener("done", (e) => {
  const payload = JSON.parse(e.data);
  setLoading(false);
  eventSource.close();
});
```

### Timing budget

| Event | Expected time | What to render |
|-------|--------------|----------------|
| `accepted` | 0 ms | Show spinner / skeleton |
| `data` | ~15 ms | Structured answer (render the real content) |
| `token` (first) | ~500 ms | Start appending prose |
| `token` (last) | ~8 s | Prose complete |
| `done` | ~8 s | Hide spinner, cache response |

**The key rule: render on `data`, not on `done`.** The structured data
arrives in 15 ms. The prose takes 8 seconds. Do not make the user wait
for prose to see the answer.

---

## Error Responses

Every error has the same envelope:

```json
{
  "type": "rate_limited",
  "title": "Upstream rate limit exceeded",
  "status": 503,
  "detail": "The LLM provider returned a rate limit error. Retry after the indicated period.",
  "request_id": "uuid"
}
```

### Error codes

| HTTP Status | `type` | When it occurs | Frontend action |
|-------------|--------|----------------|-----------------|
| 422 | `validation_error` | Missing or invalid `query` field | Show input validation error |
| 413 | `payload_too_large` | Query > 2000 chars or body > 64 KB | Show "query too long" message |
| 429 | `client_rate_limited` | Too many requests from this IP (>30/min) | Show "slow down" message, check `Retry-After` header |
| 503 | `rate_limited` | Groq quota exhausted, no data obtained | Show "service busy", retry after `Retry-After` header |
| 504 | `agent_timeout` | Agent exceeded 45s execution budget | Show "request timed out" |
| 500 | `internal_error` | Unhandled error (no details exposed) | Show generic error |

### Degraded responses (HTTP 200 but incomplete)

When the LLM is unavailable but the tool succeeded, the response is
still HTTP 200 with data present and the degradation flagged:

```json
{
  "answer": null,
  "data": { "found": true, "symbol": "BRAF", ... },
  "degraded": ["llm_unavailable"],
  "latency_ms": 45
}
```

The frontend should:
- Render the structured `data` normally
- Show a subtle banner: "Explanation unavailable — showing data only"
- Never show an error page for a 200 response

---

## GET /health

```
GET /health
```

Returns:

```json
{
  "status": "ok"
}
```

This endpoint makes zero outbound calls. It is safe to poll at any
frequency. Use it for container liveness probes, uptime monitors,
or frontend "is the agent running?" checks.

---

## CORS Configuration

The agent allows requests from `http://localhost:5173` (Vite dev server)
by default. If the frontend runs on a different origin, update the
`ALLOWED_ORIGINS` environment variable in `AIAgent/.env`:

```
ALLOWED_ORIGINS=http://localhost:5173,http://localhost:3000
```

---

## Rate Limits

### Client-side (per IP)

30 requests per minute. Exceeding this returns HTTP 429 with a
`Retry-After` header.

### LLM-side (Groq free tier)

8,000 tokens per minute for `openai/gpt-oss-120b`. The agent uses
~2,000 tokens per request (system prompt + tool output + response).
This allows approximately 3–4 requests per minute before hitting the
Groq rate limit.

**Practical implication:** rapid sequential agent calls (e.g. a user
clicking "explain" on 5 cell lines within 30 seconds) will trigger a
Groq 429 error. The frontend should:
- Debounce rapid clicks
- Show a "please wait" message if a Groq rate limit is returned
- Use the `cached` field — repeated identical queries return instantly

---

## Caching

Exact repeated queries are cached for 1 hour (1024 entries max).
A cached response has `cached: true` and returns in <10 ms.

Cache keys are normalised: `"What are the aliases for BRAF"` and
`"what are the aliases for braf"` are the same cache key.

---

## Quick-start examples

### curl

```bash
# Gene alias
curl -X POST /v1/agent/query \
  -H "Content-Type: application/json" \
  -d '{"query": "What are the aliases for BRAF"}'

# Score explanation
curl -X POST /v1/agent/query \
  -H "Content-Type: application/json" \
  -d '{"query": "Explain the score for ENSG00000157764 in ACH-000620"}'

# Dataset info
curl -X POST /v1/agent/query \
  -H "Content-Type: application/json" \
  -d '{"query": "What is depmap"}'

# Health check
curl /health
```

### JavaScript (fetch)

```javascript
const AGENT_API = "";

async function queryAgent(queryString) {
  const res = await fetch(`${AGENT_API}/v1/agent/query`, {
    method: "POST",
    headers: { "Content-Type": "application/json" },
    body: JSON.stringify({ query: queryString }),
  });
  if (!res.ok) {
    const err = await res.json();
    throw new Error(err.detail || `Agent error: ${res.status}`);
  }
  return res.json();
}

// BC1: Gene alias
const alias = await queryAgent("What are the aliases for BRAF");

// BC2: Score methodology
const methodology = await queryAgent("Explain the score for ENSG00000157764 in ACH-000620");

// BC3: Dataset info
const dataset = await queryAgent("What is depmap");
```