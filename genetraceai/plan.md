# Wire GeneTraceAI Frontend to Real APIs

## Current Architecture Analysis

The entire frontend lives in a **single 1042-line file**: [App.jsx](file:///d:/MSc%20Data%20Science/Dissertation/AstraZeneca/UoB-GeneTraceAI-25-26/genetraceai/src/App.jsx). There are no separate components, pages, services, or API directories. The CSS is embedded as a template literal at [line 800](file:///d:/MSc%20Data%20Science/Dissertation/AstraZeneca/UoB-GeneTraceAI-25-26/genetraceai/src/App.jsx#L800-L1041).

---

## A. Search Flow (Current)

| Question | Answer |
|----------|--------|
| **Search bar component** | [`SearchBox`](file:///d:/MSc%20Data%20Science/Dissertation/AstraZeneca/UoB-GeneTraceAI-25-26/genetraceai/src/App.jsx#L537-L578) — renders on the landing page (big variant) and in the topbar (small variant) |
| **What happens on submit** | Calls [`submit(term)`](file:///d:/MSc%20Data%20Science/Dissertation/AstraZeneca/UoB-GeneTraceAI-25-26/genetraceai/src/App.jsx#L619-L622) which uppercases the term and sets `gene` state. A `useEffect` on `gene` (line 609) calls `fetchGene(gene)` |
| **Mock/static data source** | [`fetchGene(symbol)`](file:///d:/MSc%20Data%20Science/Dissertation/AstraZeneca/UoB-GeneTraceAI-25-26/genetraceai/src/App.jsx#L122-L126) fetches `GET /data/{SYMBOL}.json` — a **static JSON file** served from `public/data/`. There is also [`fetchIndex()`](file:///d:/MSc%20Data%20Science/Dissertation/AstraZeneca/UoB-GeneTraceAI-25-26/genetraceai/src/App.jsx#L116-L120) that loads `/data/index.json` for the autocomplete dropdown. |
| **Problem** | **There is no `public/data/` directory.** No JSON files exist. Both fetches will return 404. The frontend currently shows nothing. |

### Data shape the UI expects (from `adaptPayload`)

The [`adaptPayload`](file:///d:/MSc%20Data%20Science/Dissertation/AstraZeneca/UoB-GeneTraceAI-25-26/genetraceai/src/App.jsx#L87-L112) function transforms raw JSON into this shape:

```javascript
{
  symbol: "BRAF",
  ensg: "ENSG00000157764",
  name: "B-Raf Proto-Oncogene",
  role: "oncogene",         // from gene["COSMIC role"]
  chr: "",                  // from gene["locus type"]
  geneMeta: {},             // entire gene metadata block
  klass: "Tier 1",          // j.class
  nCandidates: 1523,        // total scored
  verdict: { state, headline, detail[], qualify_ranking },
  lines: [{
    name: "jhh-1",
    model_id: "ACH-000620",
    lineage: "Liver",       // capitalised
    subtype: "Hepatocellular",
    score: 0.9996,          // core_score, float [0,1]
    tier: "high",           // lowercase: high|moderate|low|prior|unknown
    reason: "...",
    nLayers: 2,
    nModalities: 2,
    metadata: {},
    tracks: {               // per-evidence-layer breakdown
      expression: { present: true, score: 0.95, detail: {...}, finding: "..." },
      proteomics: { present: true, score: 0.88, ... },
      cna: { present: true, score: 0.7, total_cn: 4, direction: "gain", ... },
      mutation: { present: true, score: 1, driver_alteration: true, ... },
      fusion: { present: false, score: 0, ... },
      signature: { present: false, score: 0, ... },
    },
  }],
}
```

> [!IMPORTANT]
> The API returns a **flat summary** (`rank, model_id, name, score, tier, n_layers, driver_alteration`). The frontend expects a **richer object** with `lineage`, `subtype`, `tracks` (6 evidence layers), `reason`, `metadata`, etc. This is the **key mismatch** that must be bridged.

---

## B. Results Display (Current)

| Component | Location | What it renders |
|-----------|----------|----------------|
| [`CellLineRow`](file:///d:/MSc%20Data%20Science/Dissertation/AstraZeneca/UoB-GeneTraceAI-25-26/genetraceai/src/App.jsx#L177-L200) | Ranked list view | Row with: name, model_id, lineage dot, subtype, 6-segment evidence spine, score bar |
| [`NetworkGraph`](file:///d:/MSc%20Data%20Science/Dissertation/AstraZeneca/UoB-GeneTraceAI-25-26/genetraceai/src/App.jsx#L242-L296) | Network view | Force-directed network with gene at centre, cell lines as nodes |
| [`OmicsMap`](file:///d:/MSc%20Data%20Science/Dissertation/AstraZeneca/UoB-GeneTraceAI-25-26/genetraceai/src/App.jsx#L301-L404) | Omics map view | Sankey-like diagram: Gene → 6 evidence tracks → cell lines |
| [`DetailPanel`](file:///d:/MSc%20Data%20Science/Dissertation/AstraZeneca/UoB-GeneTraceAI-25-26/genetraceai/src/App.jsx#L409-L495) | Right sidebar | Full evidence readout for selected cell line |
| [`VerdictBanner`](file:///d:/MSc%20Data%20Science/Dissertation/AstraZeneca/UoB-GeneTraceAI-25-26/genetraceai/src/App.jsx#L507-L532) | Above results | Abstention/warning banner for genes that can't be ranked |

**Detail view on click**: Yes — clicking a row sets `selected` state, which renders `DetailPanel` with the full evidence breakdown.

**Exclusion/selectivity view**: **Not implemented in the UI.** There is no two-gene search form or selectivity table in the current frontend.

---

## C. Existing API Calls

| Call | Location | URL |
|------|----------|-----|
| [`fetchIndex()`](file:///d:/MSc%20Data%20Science/Dissertation/AstraZeneca/UoB-GeneTraceAI-25-26/genetraceai/src/App.jsx#L116-L120) | On mount | `GET /data/index.json` (static file, **does not exist**) |
| [`fetchGene(symbol)`](file:///d:/MSc%20Data%20Science/Dissertation/AstraZeneca/UoB-GeneTraceAI-25-26/genetraceai/src/App.jsx#L122-L126) | On gene selection | `GET /data/{SYMBOL}.json` (static file, **does not exist**) |

No axios. No API client module. No third-party HTTP library. Uses native `fetch`.

---

## D. State Management

All state is `useState` inside the [`App`](file:///d:/MSc%20Data%20Science/Dissertation/AstraZeneca/UoB-GeneTraceAI-25-26/genetraceai/src/App.jsx#L590-L795) component. No Context, Redux, or Zustand.

| State variable | Line | Purpose |
|---------------|------|---------|
| `query` | [591](file:///d:/MSc%20Data%20Science/Dissertation/AstraZeneca/UoB-GeneTraceAI-25-26/genetraceai/src/App.jsx#L591) | Search input text |
| `gene` | [592](file:///d:/MSc%20Data%20Science/Dissertation/AstraZeneca/UoB-GeneTraceAI-25-26/genetraceai/src/App.jsx#L592) | Currently searched gene symbol (triggers fetch) |
| `geneData` | [600](file:///d:/MSc%20Data%20Science/Dissertation/AstraZeneca/UoB-GeneTraceAI-25-26/genetraceai/src/App.jsx#L600) | Transformed payload from `adaptPayload()` |
| `selected` | [593](file:///d:/MSc%20Data%20Science/Dissertation/AstraZeneca/UoB-GeneTraceAI-25-26/genetraceai/src/App.jsx#L593) | Currently clicked cell line (for detail panel) |
| `index` | [599](file:///d:/MSc%20Data%20Science/Dissertation/AstraZeneca/UoB-GeneTraceAI-25-26/genetraceai/src/App.jsx#L599) | Gene list for autocomplete |
| `loading` | [601](file:///d:/MSc%20Data%20Science/Dissertation/AstraZeneca/UoB-GeneTraceAI-25-26/genetraceai/src/App.jsx#L601) | Loading spinner state |
| `loadError` | [602](file:///d:/MSc%20Data%20Science/Dissertation/AstraZeneca/UoB-GeneTraceAI-25-26/genetraceai/src/App.jsx#L602) | Error message state |
| `lineageFilter` | [594](file:///d:/MSc%20Data%20Science/Dissertation/AstraZeneca/UoB-GeneTraceAI-25-26/genetraceai/src/App.jsx#L594) | Lineage dropdown filter |
| `minConf` | [595](file:///d:/MSc%20Data%20Science/Dissertation/AstraZeneca/UoB-GeneTraceAI-25-26/genetraceai/src/App.jsx#L595) | Minimum score slider |
| `view` | [596](file:///d:/MSc%20Data%20Science/Dissertation/AstraZeneca/UoB-GeneTraceAI-25-26/genetraceai/src/App.jsx#L596) | Active view tab (list/network/omics) |
| `showScoring` | [597](file:///d:/MSc%20Data%20Science/Dissertation/AstraZeneca/UoB-GeneTraceAI-25-26/genetraceai/src/App.jsx#L597) | Scoring methodology panel toggle |

---

## Proposed Changes

### 1. Create API client

#### [NEW] [client.js](file:///d:/MSc%20Data%20Science/Dissertation/AstraZeneca/UoB-GeneTraceAI-25-26/genetraceai/src/api/client.js)

Create `src/api/client.js` with 4 functions using native `fetch`. See exact contents in the execution step.

---

### 2. Replace `fetchGene()` and `fetchIndex()` — the core wiring change

#### [MODIFY] [App.jsx](file:///d:/MSc%20Data%20Science/Dissertation/AstraZeneca/UoB-GeneTraceAI-25-26/genetraceai/src/App.jsx)

> [!IMPORTANT]
> **The API response shape differs from what `adaptPayload` expects.** The API returns a flat list of `{ rank, model_id, name, score, tier, n_layers, driver_alteration }`. The UI expects `{ lineage, subtype, tracks, reason, metadata, ... }`. Two options:
>
> **Option A (recommended):** Write an `adaptApiResponse()` function that maps the API response to the shape `geneData` expects, filling in missing fields with reasonable defaults (lineage: "Unknown", tracks: empty, etc). The list view, network view, and omics map will render with degraded but functional data — the score and tier will be correct, but evidence spines will be empty and lineage dots will all be grey.
>
> **Option B:** Extend the scoring API to return the full evidence-layer breakdown. This is a backend change and out of scope.

**Specific changes to [App.jsx](file:///d:/MSc%20Data%20Science/Dissertation/AstraZeneca/UoB-GeneTraceAI-25-26/genetraceai/src/App.jsx):**

1. **Lines 1–2** — Add import for the API client:
   ```diff
   +import { searchGene, getGeneCellLine, getSelectivity, queryAgent } from "./api/client";
   ```

2. **Lines 65–126** — Replace `DATA_BASE`, `adaptPayload`, `fetchIndex`, `fetchGene` with a new `adaptApiResponse()` that maps the scoring API response to the shape `geneData` expects:
   ```javascript
   function adaptApiResponse(apiData) {
     return {
       symbol: apiData.gene,
       ensg: apiData.ensg,
       name: "",                    // not in API response
       role: "unknown",             // not in API response
       chr: "",
       geneMeta: {},
       klass: "",
       nCandidates: apiData.total,
       verdict: null,               // no verdict from API
       lines: (apiData.lines || []).map((l) => ({
         name: l.name || l.model_id,
         model_id: (l.model_id || "").toUpperCase(),
         lineage: "Unknown",        // not in API response
         subtype: "—",
         score: l.score,
         tier: (l.tier || "unknown").toLowerCase(),
         reason: "",
         nLayers: l.n_layers,
         nModalities: l.n_layers,
         metadata: {},
         tracks: buildTracksFromApi(l),
       })),
       rnaAlternatives: apiData.rna_alternatives_for_rank1 || [],
     };
   }
   ```

3. **Lines 604–617** — Replace the two `useEffect` blocks:
   - Remove `fetchIndex()` call on mount (no index.json to load; the search bar will work as a free-text input without autocomplete)
   - Replace `fetchGene(gene)` with `searchGene(gene)` and transform via `adaptApiResponse()`

4. **Lines 541–545** — In `SearchBox`, handle empty `index` gracefully (the API doesn't provide an index of available genes, so autocomplete will be empty until we add the agent alias lookup later)

5. **Lines 662–665** — Landing page "Try" buttons: replace `index.map(...)` with a hardcoded list of example genes: `["BRAF", "EGFR", "TP53", "KRAS", "PIK3CA"]`

---

### 3. Error handling

The existing `loadError` and `loading` states are already wired into the UI ([lines 675–691](file:///d:/MSc%20Data%20Science/Dissertation/AstraZeneca/UoB-GeneTraceAI-25-26/genetraceai/src/App.jsx#L675-L691)). The API client will throw on non-2xx responses, and the existing `catch` block at line 615 will set `loadError`. No additional error UI is needed.

However, the current "no result" state at [line 681](file:///d:/MSc%20Data%20Science/Dissertation/AstraZeneca/UoB-GeneTraceAI-25-26/genetraceai/src/App.jsx#L681) checks `!geneData` — this will also trigger on API errors. We should also display the actual error message:

```diff
-<h2>No models on record for "{gene}".</h2>
+<h2>{loadError ? "Unable to reach scoring API" : `No models on record for "${gene}".`}</h2>
+{loadError && <p className="method-note">{loadError}</p>}
```

---

## Open Questions

> [!IMPORTANT]
> **Tier mapping**: The API returns tiers as `"HIGH"`, `"MEDIUM"`, `"CONTEXT"`, `"LOW"` (uppercase). The frontend's [`confTier`](file:///d:/MSc%20Data%20Science/Dissertation/AstraZeneca/UoB-GeneTraceAI-25-26/genetraceai/src/App.jsx#L135-L141) function maps `"high"`, `"moderate"`, `"low"`, `"prior"`, `"unknown"` (lowercase, different names). Should I:
> - Map `HIGH→high`, `MEDIUM→moderate`, `CONTEXT→low`, `LOW→unknown`? (best guess)
> - Or leave the tier as-is and add new entries to `confTier` for the uppercase API tiers?

> [!IMPORTANT]
> **Network and Omics views**: These visualisations depend on `lineage` (for colour grouping) and `tracks` (for evidence flow widths). The API doesn't return these. With defaults, the network will show all grey nodes and the omics map will show no evidence flows. Is this acceptable as a first pass, or should we hide those view tabs until the data model is enriched?

> [!IMPORTANT]
> **Exclusion/selectivity search**: The current frontend has no UI for two-gene queries. The API supports `GET /exclude?gene_a=...&gene_b=...`. Should I add a selectivity search UI in this wiring pass, or defer to a later iteration?

> [!IMPORTANT]
> **Agent API integration**: The agent API (alias lookup, score explanation, dataset info) requires `localhost:8000` to be running. Should I wire these in this pass (with graceful fallback when the agent is offline), or defer?

---

## Verification Plan

### Manual Verification
1. Run `npm run dev` and open the Vite dev server
2. Search for **BRAF** — verify the results table populates with real data from the scoring API
3. Click a cell line row — verify the detail panel opens (with degraded evidence readout since tracks aren't available)
4. Search for a nonexistent gene (e.g. "ZZZZZ") — verify the error/no-results state displays cleanly
5. Verify the Network and Omics Map views render (even if with default colours/empty evidence flows)
6. Check the browser console for CORS errors from the AWS Lambda endpoint

### Automated Tests
- None currently exist. No test framework is set up.
