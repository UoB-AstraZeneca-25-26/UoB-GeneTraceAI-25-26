# Scoring Methodology Visualisation — Implementation Plan

Two changes: the AI agent returns structured steps (not prose), and the
frontend renders them as an interactive animated chart.

---

## Part 1 — AI Agent: Structured Step Output

### Goal

Replace the free-form prose explanation with a structured JSON response
that the frontend can render directly as chart nodes. Each step is a
key-value pair: the key is the chart label, the value is a 20–25 word
description shown on click.

### Change 1: Update the system prompt (Schema/Prompts.py)

In the rules section of the system prompt, add this instruction:

```
- For score_explainer results: respond ONLY with valid JSON. 
  No prose before or after. The JSON must be an object with a single
  key "steps" containing an array of exactly 7 objects.
  
  Each object has:
    - "key": short label (3-5 words, used as the chart node title)
    - "value": explanation (20-25 words, shown on click)
    - "formula": the mathematical formula for this step (optional, short)
    - "status": one of "active" | "skipped" | "inverted"
      Use "active" for steps that applied normally.
      Use "skipped" if the step did not apply (e.g. TSG inversion for an oncogene).
      Use "inverted" for TSG inversion when the gene IS a tumour suppressor.

  Example output format:
  {
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
        "value": "For tumour suppressors, low abundance signals loss of function. Score inverted so higher always means more relevant, regardless of gene class.",
        "formula": "score_TSG = 1 − core_score",
        "status": "skipped"
      },
      {
        "key": "Driver gating",
        "value": "Cell lines carrying a known driver alteration (mutation, fusion, or CNA matching the gene's role) are promoted above all non-driver lines.",
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
  }

  CRITICAL RULES for this output:
  - Return ONLY the JSON object. No markdown, no backticks, no preamble.
  - Always include all 7 steps in order.
  - Set "status" to "skipped" for TSG inversion if the gene is an oncogene.
  - Set "status" to "inverted" for TSG inversion if the gene IS a TSG.
  - Keep each "value" to 20-25 words. Not shorter, not longer.
  - Keep each "key" to 3-5 words.
  - If the gene name is known (from the ensg_id), mention it by name
    in relevant steps (e.g. "BRAF is an oncogene, so no inversion").
```

### Change 2: Update the agent response parsing (api/routes.py)

In `_agent_events()` or the response builder, when the tool is
`score_explainer`, attempt to parse the LLM's answer as JSON:

```python
import json

# After getting the LLM answer text:
if tool_name == "score_explainer":
    try:
        parsed = json.loads(answer_text)
        if "steps" in parsed:
            response.answer = None          # no prose needed
            response.methodology = parsed   # new field for the frontend
    except json.JSONDecodeError:
        response.answer = answer_text       # fallback to prose if parsing fails
```

### Change 3: Update response schema (api/schemas.py)

Add an optional `methodology` field to `QueryResponse`:

```python
class MethodologyStep(BaseModel):
    key: str
    value: str
    formula: str | None = None
    status: str = "active"  # active | skipped | inverted

class QueryResponse(BaseModel):
    answer: str | None = None
    data: dict | None = None
    methodology: dict | None = None  # {"steps": [MethodologyStep, ...]}
    tool_calls: list = []
    degraded: list = []
    cached: bool = False
    latency_ms: int = 0
    request_id: str = ""
```

### Change 4: Reduce token usage

The structured JSON response uses ~300 tokens (vs ~800 for prose).
Combined with max_tokens=512, this fits comfortably within one
request against the 8K TPM limit.

No changes needed to max_tokens — 512 is sufficient for the 7-step
JSON output.

---

## Part 2 — Frontend: Interactive Step Chart

### Goal

Render the 7 methodology steps as an animated vertical pipeline chart
with hover effects, click-to-expand descriptions, colour transitions,
and smooth animations. No external charting library needed — pure React
+ CSS animations + framer-motion.

### Install required package

```bash
npm install framer-motion
```

This is the only new dependency. It handles enter/exit animations,
hover states, and layout transitions natively in React.

### Create the component: src/components/MethodologyChart.jsx

```
MethodologyChart receives: { steps, geneName }

steps = [
  { key: "Raw measurements", value: "...", formula: "...", status: "active" },
  { key: "PIT normalisation", value: "...", formula: "...", status: "active" },
  ...7 items
]
```

### Visual design specification

```
LAYOUT: Vertical pipeline, each step is a node connected by animated lines

    ┌─────────────────────────────┐
    │  ① Raw measurements        │ ← step node
    │     log₂(TPM+1), log₂(r)   │ ← formula (always visible, muted)
    └─────────────────────────────┘
              │  (animated line, draws downward)
              ▼
    ┌─────────────────────────────┐
    │  ② PIT normalisation       │
    │     F̂(x) = rank(x)/(n+1)   │
    └─────────────────────────────┘
              │
              ▼
          ... 5 more ...
```

### Step node states

| State | Background | Border | Text | Icon |
|-------|-----------|--------|------|------|
| Default (idle) | white | #e5e7eb (gray-200) | #374151 (gray-700) | step number, gray |
| Hover | #f0fdfa (teal-50) | #0d9488 (teal-600) | #0d9488 | step number, teal |
| Active (clicked/expanded) | #0d9488 (teal-600) | #0d9488 | white | check icon, white |
| Skipped | #f9fafb (gray-50) | #d1d5db (gray-300) dashed | #9ca3af (gray-400) | skip icon, gray |
| Inverted | #fef3c7 (amber-50) | #f59e0b (amber-500) | #92400e (amber-800) | swap icon, amber |

### Animations

1. **Entrance**: Steps stagger in from left, 100ms apart, using
   framer-motion's `variants` with `staggerChildren`:
   ```jsx
   const container = {
     hidden: {},
     show: { transition: { staggerChildren: 0.1 } }
   };
   const item = {
     hidden: { opacity: 0, x: -30 },
     show: { opacity: 1, x: 0, transition: { type: "spring", stiffness: 300, damping: 24 } }
   };
   ```

2. **Connecting lines**: Each line draws from top to bottom using
   SVG `stroke-dashoffset` animation. Start as fully dashed (invisible),
   animate to fully drawn, timed to complete just before the next
   step node appears.

3. **Hover**: Scale up slightly (1.02), border colour transitions
   to teal, subtle box-shadow appears. Use framer-motion `whileHover`:
   ```jsx
   whileHover={{ scale: 1.02, boxShadow: "0 4px 12px rgba(13,148,136,0.15)" }}
   ```

4. **Click to expand**: The description panel slides down with
   `AnimatePresence` + `motion.div` using `layout` prop for smooth
   height transition:
   ```jsx
   <AnimatePresence>
     {expanded && (
       <motion.div
         initial={{ height: 0, opacity: 0 }}
         animate={{ height: "auto", opacity: 1 }}
         exit={{ height: 0, opacity: 0 }}
         transition={{ duration: 0.25 }}
       >
         <p>{step.value}</p>
       </motion.div>
     )}
   </AnimatePresence>
   ```

5. **Colour pulse on expand**: When a step is clicked, the border
   briefly pulses brighter (teal-400 → teal-600) using a CSS keyframe:
   ```css
   @keyframes pulse-border {
     0% { border-color: #2dd4bf; }
     50% { border-color: #0d9488; }
     100% { border-color: #0d9488; }
   }
   ```

6. **Skipped steps**: Render with dashed border, muted colours,
   and a small "skipped" label. On hover, show a tooltip explaining
   why (e.g. "BRAF is an oncogene — inversion applies only to TSGs").

### Expanded description panel (on click)

```
    ┌─────────────────────────────────────────────┐
    │  ④ Core score                    ✓ active   │
    │     core = w_E × E + w_P × P                │
    ├─────────────────────────────────────────────┤
    │  Weighted sum of PIT ranks using             │
    │  correlation-penalised weights. For two      │
    │  exchangeable layers, optimal weights are     │
    │  equal: 0.5 each.                            │
    └─────────────────────────────────────────────┘
```

The description panel shows:
- The "value" text (20-25 words)
- A subtle divider line between formula and description
- Click again to collapse (toggle behaviour)

### Where to place in the existing UI

The chart replaces the current plain-text "Scoring methodology" section
in the right sidebar panel. Find where the current explanation text
renders (the block that starts with "Step 1 – Raw measurements...") and
replace it with:

```jsx
{methodology ? (
  <MethodologyChart steps={methodology.steps} geneName={geneData?.gene} />
) : explaining ? (
  <LoadingSpinner label="Loading methodology..." />
) : (
  <button onClick={handleExplain}>Explain scoring methodology</button>
)}
```

### Loading state

While waiting for the agent response (~8s), show a skeleton version
of the chart: 7 gray placeholder nodes with a shimmer animation,
no text. This gives immediate visual feedback that something is loading
in the shape it will take.

```jsx
const SkeletonStep = () => (
  <motion.div
    style={{
      height: 56,
      borderRadius: 8,
      background: "linear-gradient(90deg, #f3f4f6 25%, #e5e7eb 50%, #f3f4f6 75%)",
      backgroundSize: "200% 100%",
    }}
    animate={{ backgroundPosition: ["200% 0", "-200% 0"] }}
    transition={{ duration: 1.5, repeat: Infinity, ease: "linear" }}
  />
);
```

### Error/fallback handling

If the agent returns prose instead of JSON (parsing fails), fall back
to rendering the raw text in a styled paragraph — same as the current
behaviour. The chart is an enhancement, not a hard dependency.

```jsx
{methodology?.steps ? (
  <MethodologyChart steps={methodology.steps} geneName={geneData?.gene} />
) : explanation ? (
  <div className="prose-fallback">{explanation}</div>
) : null}
```

---

## Integration Flow Summary

```
User clicks "Explain scoring methodology"
    │
    ▼
Frontend calls: POST /v1/agent/query
    { "query": "Explain the score for ENSG... in ACH-..." }
    │
    ▼
Agent calls score_explainer tool → returns structured data
    │
    ▼
LLM reads math_reference.md + tool output
    │
    ▼
LLM returns structured JSON:
    { "steps": [{ "key": "...", "value": "...", "formula": "...", "status": "..." }, ×7] }
    │
    ▼
api/routes.py parses JSON → sets response.methodology
    │
    ▼
Frontend receives { methodology: { steps: [...] } }
    │
    ▼
MethodologyChart renders 7 animated step nodes
    │
    ├── Entrance: stagger animation (700ms total)
    ├── Connecting lines: draw-down SVG animation
    ├── Hover: scale + colour + shadow transition
    ├── Click: expand description panel with slide animation
    ├── Skipped steps: dashed border + muted + tooltip
    └── Inverted steps: amber theme + swap icon
```

---

## File Changes Summary

| File | Change |
|------|--------|
| `Schema/Prompts.py` | Add structured JSON output rules to system prompt |
| `api/schemas.py` | Add `methodology` field to QueryResponse |
| `api/routes.py` | Parse score_explainer answer as JSON, populate methodology field |
| `genetraceai/src/components/MethodologyChart.jsx` | NEW — the animated step chart component |
| `genetraceai/src/components/MethodologyChart.css` | NEW — animations and step node styles |
| `genetraceai/package.json` | Add framer-motion dependency |
| Sidebar component (find by tracing) | Replace prose explanation with MethodologyChart |

---

## Verification

1. Search "BRAF" → click any cell line → click "Explain scoring methodology"
2. Skeleton loader should appear immediately (7 shimmer nodes)
3. After ~8s, real steps appear with stagger animation
4. Step 5 (TSG inversion) should show as "skipped" (BRAF is oncogene)
5. Hover any step → teal highlight + subtle shadow
6. Click any step → description expands with slide animation
7. Click again → collapses
8. Search "TP53" → Step 5 should show as "inverted" (amber theme)
9. If agent is down → fallback to prose text or "Agent unavailable"