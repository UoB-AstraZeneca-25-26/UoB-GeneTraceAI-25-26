# Cell Line Finder — Backend Integration & Design

How the frontend is architected, how it talks to the backend today, and exactly what the
backend must provide (and in what shape) to run on live data. This reflects the current
codebase.

---

## 1. Purpose

Cell Line Finder ranks cancer cell-line models for a target gene (or a joint / selectivity
query) using multi-omics evidence, and lets a user drill into any line's evidence profile.
The frontend is a React + TypeScript + Vite single-page app. It is **presentation only** —
all scoring, ranking, and evidence live in the backend. This document is the contract
between the two.

---

## 2. Architecture at a glance

```
QueryBuilder ──(QueryParams)──▶ App.handleSearchSubmit
                                      │
                                      ▼
                         lib/api.ts  fetchRanking()            ← the single entry point
                          ├─ routes by query shape:
                          │    ≥2 targets            → GET /prod/genes        (multi)
                          │    1 target + exclusion  → GET /prod/exclude/many (selectivity)
                          │    1 target              → GET /prod/gene         (single)
                          ├─ transport: getJson / getJsonLenient
                          └─ adapters: toSingleRanked / toMultiRanked /
                                       toSelectivityRanked  (snake_case → domain types)
                                      │
                       { results: RankedCellLine[], meta: ResultMeta }
                                      │
                                      ▼
   ResultsTable ── Table / Evidence grid / Lineages ──▶ inspect(modelId)
        │                    │                                   │
        │        EvidenceMatrix fetches                          ▼
        │        GET /prod/gene/detail per (gene × line)   ProfileView fetches
        │                                                  GET /prod/gene/detail
        │                                                  (+ uses ranked row for context)
        ▼
   Assistant / "did you mean"  ──▶ lib/assistant.ts ──▶ POST /v1/agent/query
```

**Layered design (why it's easy to wire):**

| Layer | Files | Responsibility |
|---|---|---|
| Transport | `src/lib/api.ts` (`getJson`, `getJsonLenient`) | fetch, error shaping, JSON repair |
| Endpoints | `src/lib/api.ts` (`fetch*`) | URL building, query params, `USE_MOCK` switch |
| Adapters | `src/lib/api.ts` (`to*Ranked`, `toCellLineDetail`) | map API snake_case → domain camelCase, compute derived fields |
| Domain types | `src/types/index.ts` | `RankedCellLine`, `CellLineDetail`, `ResultMeta`, `QueryParams` |
| Mock | `src/lib/mockApi.ts`, `src/lib/mockData.ts` | mirror the exact API shapes behind `USE_MOCK` |
| Agent | `src/lib/assistant.ts`, `src/config.ts` | chat + gene-alias lookup |

Components never fetch on their own except through the `api.ts` fetchers. **All snake_case
JSON keys below are read verbatim by the adapters — they must match exactly.**

---

## 3. How it works today (mock vs. live)

- `src/lib/api.ts` exports **`USE_MOCK`** (currently `false`). When `true`, every fetcher
  returns local data from `mockApi.ts` in the identical shapes — the app runs with no
  backend. Set `false` for live.
- **Base URLs are environment-switched** (`import.meta.env.DEV`):
  - **Dev:** requests go to same-origin paths (`/api`, `/agent-api`) that **Vite proxies**
    server-side to the real APIs (`vite.config.ts`). This exists because neither API sends
    CORS headers yet, so the browser never makes a cross-origin call in dev.
  - **Prod build:** requests go directly to the deployed URLs.
- Scoring base (prod): `https://9368clqa34.execute-api.eu-north-1.amazonaws.com`, all under
  the **`/prod`** stage. Agent base (prod):
  `https://mhyvpwmjma3gqy44dk4wskkxje0gqfud.lambda-url.eu-west-2.on.aws`.
- `TOP_N = 30` is sent as `top_n`; the UI slices further to the "Number of recommendations"
  slider value.

**Config touchpoints** (the only places URLs/flags live):

| What | Where |
|---|---|
| `USE_MOCK` on/off | `src/lib/api.ts` |
| Scoring base URL + dev-proxy switch | `src/lib/api.ts` (top) and `vite.config.ts` (`/api`) |
| Agent base URL + dev-proxy switch | `src/config.ts` and `vite.config.ts` (`/agent-api`) |
| `top_n` | `src/lib/api.ts` (`TOP_N`) |

---

## 4. Endpoints the frontend calls

| Purpose | Method & path | Query / body | Parse |
|---|---|---|---|
| Single-gene ranking | `GET /prod/gene` | `gene=<SYM>&top_n=30` | lenient |
| Multi-gene joint ranking | `GET /prod/genes` | `genes=<A,B,…>&top_n=30` | lenient |
| Selectivity ranking | `GET /prod/exclude/many` | `gene_a=<HIGH>&exclude=<LOW1>&exclude=<LOW2>…&top_n=30` | **strict** |
| Cell-line detail | `GET /prod/gene/detail` | `gene=<SYM>&model_id=<ACH>` | **strict** |
| Agent (chat + alias) | `POST /v1/agent/query` | `{ "query": "…" }` | strict |

- **lenient** = bare `NaN`/`Infinity` are repaired to `null` before parsing (per-gene scores
  can legitimately be "no evidence").
- **strict** = parsed as-is; `/exclude/many` and `/gene/detail` **must not emit `NaN`/
  `Infinity`** — send `null`.
- Selectivity sends **one `exclude=` param per excluded gene** (repeated), not a CSV.

---

## 5. Ranking responses

### 5.1 `/prod/gene` and `/prod/genes` — `JointApiResponse`

Both endpoints return the **same** shape. Lines are expected **pre-sorted best-first** (the
UI numbers rank by array order).

```jsonc
{
  "genes": ["EGFR"],
  "total_passing": 1039,
  "floor": 0,
  "lines": [
    {
      "model_id": "ACH-000552",       // ACH id — PRIMARY identity everywhere
      "name": "A549",                 // display name; null or "None" → falls back to model_id
      "joint_score": 0.9847,          // 0–1 GLOBAL score (single: the gene's score; multi: joint)
      "scores": { "EGFR": 0.9847 },   // per-gene scores; null = no evidence for that gene
      "limiting_gene": "EGFR",        // (multi) gene capping the joint score → "▼" marker
      "lineage_score": 0.981,         // OPTIONAL 0–1 within-lineage score (§7)
      "metadata": { "lineage": "lung" }
    }
  ]
}
```

| Field | Req? | Consumed by |
|---|---|---|
| `genes` | yes | result header; multi per-gene columns |
| `total_passing` | yes | "Top N of M"; population context |
| `lines[].model_id` | **yes** | primary label (table, profile, grid, lineages), inspect key |
| `lines[].name` | rec | secondary alias under the ACH id |
| `lines[].joint_score` | **yes** | global score column; ordering; profile top card |
| `lines[].scores` | multi | "Per-gene scores" column; a `null` shows the PARTIAL badge and "▼" |
| `lines[].limiting_gene` | opt | limiting-gene marker |
| `lines[].lineage_score` | opt | in-lineage score (§7) |
| `lines[].metadata.lineage` | rec | Lineages view, lineage column + filter, in-lineage grouping |

### 5.2 `/prod/exclude/many` — `ExcludeManyApiResponse`

```jsonc
{
  "gene_a": "BRAF",
  "excluded_genes": ["EGFR"],
  "total_ranked": 912,
  "lines": [
    {
      "model_id": "ACH-000219", "name": "A-375",
      "score_a": 0.94,                    // target score (component-scores column, profile chip ↑)
      "exclusion_scores": { "EGFR": 0.08 }, // per-excluded-gene scores (chips ↓)
      "selectivity": 0.86,                // 0–1 GLOBAL score for this mode
      "lineage_score": 0.9,               // OPTIONAL (§7)
      "metadata": { "lineage": "skin" }
    }
  ]
}
```

Required: `gene_a`, `excluded_genes`, `total_ranked`; per line `model_id`, `score_a`,
`exclusion_scores`, `selectivity`. `name`, `metadata.lineage`, `lineage_score` as in 5.1.

---

## 6. Detail response — `/prod/gene/detail` (`CellLineDetailApiResponse`)

Fetched **per gene, per line** — by the profile (for the inspected gene) and by the Evidence
grid (for every gene in the query × every shown line, see §8).

```jsonc
{
  "gene": "EGFR",
  "ensg": "ENSG00000146648",
  "cell_line": { "model_id": "ACH-000552", "name": "A549" },
  "rank": 3, "total": 1523,
  "score": 0.9666,
  "tier": "HIGH",                  // HIGH | MEDIUM | LOW | CONTEXT
  "n_layers": 3,
  "driver_alteration": true,
  "p_mutation": 0.72,              // mutation PRESENT if non-null (value not shown)
  "p_fusion": null,                // fusion PRESENT if non-null
  "has_cna_alteration": true,
  "expression_level": 0.87,        // 0–1 → Low/Moderate/High band; null = not measured
  "proteomics_level": 0.64,        // 0–1 → band; null = not measured
  "expression_sources": ["DepMap","HPA"],  // OPTIONAL present/absent chips (of DepMap/HPA/GEO)
  "proteomics_sources": ["ProCan"],         // OPTIONAL present/absent chips (of ProCan/CCLE)
  "lineage_score": 0.98, "lineage_rank": 3, "lineage_total": 45,  // OPTIONAL (§7)
  "metadata": { "lineage": "lung", "primary_disease": "lung cancer", "sex": "male", "age": "58" },
  "rna_alternatives": [
    { "model_id": "ACH-000788", "name": "NCI-H1975", "similarity": 0.91,
      "lineage": "lung", "primary_disease": "lung cancer" }  // lineage/disease OPTIONAL
  ]
}
```

| Field | Req? | Consumed by |
|---|---|---|
| `gene`, `cell_line.{model_id,name}` | **yes** | profile header (ACH primary, name secondary) |
| `ensg` | opt | available in data; not shown in header currently |
| `rank`, `total`, `score` | rec | fallback cards when opened outside a ranking |
| `tier` | **yes** | tier pill — must be `HIGH`/`MEDIUM`/`LOW`/`CONTEXT` |
| `n_layers` | yes | "Evidence layers" card |
| `driver_alteration` | yes | "driver" tag on mutation |
| `p_mutation`, `p_fusion` | yes | Mutation / Fusion **Yes/No** (presence = non-null) |
| `has_cna_alteration` | yes | Copy-number Yes/No |
| `expression_level`, `proteomics_level` | rec | Low/Mod/High band + Evidence grid |
| `expression_sources`, `proteomics_sources` | opt | present/absent source chips |
| `lineage_score` / `lineage_rank` / `lineage_total` | opt | profile in-lineage card (§7) |
| `metadata` | rec | metadata block (recognised keys below) |
| `rna_alternatives[]` | opt | "RNA-similar alternatives" list |

**Recognised `metadata` keys** (any subset; unknown ignored; `"None"`/empty skipped):
`lineage`, `lineage_subtype`, `primary_disease`, `subtype`, `cellosaurus_ncit_disease`,
`primary_or_metastasis`, `sample_collection_site`, `default_growth_pattern`, `sex`, `age`.

**Source label vocabulary** (must match for chips to light up): expression → `DepMap`,
`HPA`, `GEO`; proteomics → `ProCan`, `CCLE`.

> **Not needed:** `expression_measurements` / `proteomics_measurements` (raw per-source bars)
> were removed from the UI. The types tolerate them, but the backend need not send them.

---

## 7. Design: the two scores (global + within-lineage)

Every line shows **two** scores. This is deliberate: the top lines cluster near the ceiling
globally, so "how good within its own tissue" is often the more actionable number.

- **Global** = `joint_score` / `selectivity` — always from the backend; drives the default
  ranking and the primary score everywhere.
- **In-lineage** = the app prefers the backend's `lineage_score` (0–1). **If the ranking
  response omits it, the table derives a fallback** by normalising each line against the
  other *shown* lines of the same lineage — honest but limited to the shown set. The Lineages
  filter and the "re-rank by in-lineage" table sort use this value.
- On the **profile**, in-lineage comes from `detail.lineage_score` + `lineage_rank` /
  `lineage_total` ("#3 of 45 in lung"); absent → shows `—`.

**To make in-lineage fully correct:** compute each line's score relative to all lines of the
same lineage in the full population and return it as `lineage_score` on ranking lines and on
`/gene/detail` (with `lineage_rank` / `lineage_total`). No frontend change needed — it's read
if present.

---

## 8. Design: the Evidence grid (multi-gene detail fan-out)

The Evidence grid shows cell lines (rows) × evidence layers (columns). For **selectivity and
multi-gene** queries it repeats the layer block **per gene**, so a user sees a line is high in
the target and low in the excluded gene on one screen.

- It fetches `/prod/gene/detail` for **every (gene × shown line)** pair. Example: a
  selectivity query showing 30 lines fires `2 × 30 = 60` detail calls (target + one excluded);
  a 3-gene joint query showing 20 lines fires `60`.
- Calls run in parallel with a progress counter and are abort-controlled on gene/line change.
- **Backend implication:** `/gene/detail` should be cheap and idempotent, and ideally cached,
  because it's called in bursts. If this fan-out is too heavy, the cleanest optimisation is a
  batch endpoint (e.g. `POST /prod/gene/detail/batch` taking `{gene, model_ids[]}` or
  `{pairs:[{gene,model_id}]}`) — but that's an enhancement; the app works with the per-pair
  endpoint today.

---

## 9. Design: reference tables + gene dropdown (client-side, to replace)

Three static datasets ship as a **preview** and should be swapped for full tables:

| Constant | File | Shape | Powers |
|---|---|---|---|
| `CELL_LINE_DIRECTORY` | `src/lib/mockApi.ts` | `{ ach, name, lineage, disease }[]` | Reference → Cell lines |
| `GENE_REF` | `src/lib/mockData.ts` | `{ symbol, name, ensg }[]` | Reference → Genes |
| `GENE_UNIVERSE` | `src/lib/mockData.ts` | `string[]` (from `GENE_REF`) | searchable gene dropdown |

Two ways to connect real data:

1. **Bundle** the full tables as JSON in `src/lib/` and point these constants at them
   (simplest; larger bundle).
2. **Fetch** them on load from new endpoints (e.g. `GET /prod/cell-lines`,
   `GET /prod/genes/list`) into app state (best for large panels). The dropdown accepts any
   typed symbol regardless, so it degrades gracefully.

The DepMap directory and the gene panel already exist backend-side; exposing a slim
projection (`{id,name,lineage}` / `{symbol,name,ensg}`) is all that's needed.

---

## 10. Agent endpoint (only if using assistant / alias features)

`POST {AGENT_API_URL}/v1/agent/query`, body `{ "query": "<text>" }`.

- **Assistant chat** (`sendMessage`) → expects `{ "answer": "<string>" }`.
- **"Did you mean?" gene resolution** (`lookupGeneSuggestion`) sends
  `query: "What are the aliases for <input>"` → expects
  `{ "data": { "found": bool, "symbol": "<SYM>", "full_name": "<name>" } }`.

Both fail soft (fallback message / silent "not found"); ranking never depends on the agent.
The agent's CORS is app-controlled (`ALLOWED_ORIGINS`), which is why dev proxies it.

---

## 11. Networking & JSON hygiene

- **Dev works only via the Vite proxies** (`/api`, `/agent-api`) until CORS is enabled.
- **Prod build** needs one of: the scoring API sends `Access-Control-Allow-Origin` for the
  app origin (or same-origin hosting); and the agent Lambda has the prod origin in
  `ALLOWED_ORIGINS`.
- The UI error "Network error contacting the ranking API: Failed to fetch" is almost always
  CORS or an unreachable base URL.
- `/gene` and `/genes` tolerate bare `NaN`/`Infinity` (repaired to `null`); prefer emitting
  `null`. `/exclude/many` and `/gene/detail` are **strict** — any `NaN`/`Infinity` throws.
- `name` may be `null` or the literal `"None"`; both fall back to the ACH id.

---

## 12. How integration should happen (ordered)

1. **Point at the API, keep the proxy.** Confirm `USE_MOCK = false`. Leave the dev proxies in
   place; run `npm run dev`. This is the fastest path to first live data (no CORS work).
2. **Verify the core shapes** against §5/§6 — especially `model_id`, the score field per mode,
   `scores`, `tier ∈ {HIGH,MEDIUM,LOW,CONTEXT}`, `metadata.lineage`. Fix any key mismatches
   backend-side (the adapters read exact names).
3. **Guarantee strict-JSON cleanliness** on `/exclude/many` and `/gene/detail` (no NaN/Inf).
4. **Add the optional fields** for full feature parity (§13).
5. **Replace the reference tables** (§9) — bundle or endpoints.
6. **Enable CORS** on both APIs, then remove the corresponding Vite proxy entries and confirm
   a production build calls the deployed URLs directly.
7. **Deploy** and run the smoke test (§14).

---

## 13. Checklist

**Must-have to run on live data**
- [ ] `USE_MOCK = false` for the live build.
- [ ] Four scoring endpoints reachable; CORS or dev-proxy in place.
- [ ] Shapes match §5/§6 (`model_id`, `joint_score`/`selectivity`, `scores`, `tier`, `metadata.lineage`).
- [ ] `/exclude/many` and `/gene/detail` emit no `NaN`/`Infinity`.

**For full feature parity**
- [ ] `lineage_score` on ranking lines; `lineage_score`/`lineage_rank`/`lineage_total` on detail (§7).
- [ ] `expression_sources` / `proteomics_sources` on detail (source chips).
- [ ] `expression_level` / `proteomics_level` populated (bands + grid).
- [ ] `rna_alternatives[].lineage` / `.primary_disease`.
- [ ] Replace `CELL_LINE_DIRECTORY`, `GENE_REF`, `GENE_UNIVERSE` with full tables (§9).
- [ ] Agent endpoint live if assistant / "did you mean" are wanted (§10).

**Performance / optional**
- [ ] Consider a batch detail endpoint if the grid fan-out (§8) is heavy under load.

---

## 14. Smoke test after wiring

1. `USE_MOCK = false`, `npm run dev`. Single gene (e.g. EGFR) → table populates; each row
   shows global + in-lineage scores; ACH id primary, name beneath.
2. Two targets → MULTI; per-gene columns; PARTIAL badge where a gene has no score.
3. One target + one exclusion → SELECTIVITY; component-score columns.
4. Evidence grid tab → layers load per gene (both genes for selectivity); progress counter clears.
5. Lineage filter on the table → re-ranks the chosen lineage by in-lineage score.
6. Inspect a line → header shows the query's gene(s) with scores; top cards match the clicked
   row's rank/score; alterations Yes/No; expression/proteomics bands + source chips.
7. Reference tab → search resolves a real ACH id / gene symbol / ENSG.
8. `npm run build` is clean; a production build calls the deployed URLs (not the proxy).
