# CellLineFinder AI Agent

**Model:** `openai/gpt-oss-120b` via Groq Cloud (free tier)
**Framework:** LangChain AgentExecutor with `bind_tools()`
**Status:** All 4 tests passed for the 4 business cases designed as of now

---

## Architecture

```
┌──────────────────────────────────────────────────────────────┐
│                        React UI                              │
└──────────────────────┬───────────────────────────────────────┘
                       │
                       ▼
┌──────────────────────────────────────────────────────────────┐
│                  FastAPI /agent/query                        │
└──────────────────────┬───────────────────────────────────────┘
                       │
                       ▼
┌ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─┐
│              LangChain AgentExecutor                         │
│                                                              │
│   ┌─────────────────────┐   ┌──────────────────────────┐     │
│   │    System prompt    │   │      LLM backend         │     │
│   │  Math ref + dataset │   │  ChatGroq (gpt-oss-120b) │     │
│   │  JSON (~1500 tokens)│   │  ↔ ChatOllama (swap)     │     │
│   └─────────────────────┘   └──────────────────────────┘     │
│                                                              │
│                      .bind_tools()                           │
└ ─ ─ ─ ─ ─ ─ ┬─ ─ ─ ─ ─ ─ ─┬─ ─ ─ ─ ─ ─ ─┬─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ┘
               │              │              │
               ▼              ▼              ▼
┌──────────────────┐ ┌────────────────┐ ┌──────────────────┐
│ gene_alias_lookup│ │ score_explainer│ │   dataset_info   │
│ Pydantic I/O     │ │ Pydantic I/O   │ │ Pydantic I/O     │
└────────┬─────────┘ └───────┬────────┘ └────────┬─────────┘
         │                   │                   │
         ▼                   ▼                   ▼
┌──────────────────┐ ┌────────────────┐ ┌──────────────────┐
│ Ensembl REST API │ │ Pipeline output│ │ Static JSON store│
│ HGNC REST API    │ │ Parquet files  │ │ Pre-extracted    │
│ gene_lookup table│ │ Intermediates  │ │ from AZ PDFs     │
└──────────────────┘ └────────────────┘ └──────────────────┘
    Free, no auth       Pre-computed       One-time extraction
```

**Key design decisions:**
- Gene alias resolution uses deterministic API calls, not web search
- Score explanation uses pre-computed values — LLM explains, never computes
- Dataset descriptions stored as static JSON (~1500 tokens, loaded in system prompt)
- LLM backend is swappable via `.env` (Groq → Ollama → OpenAI)

---

## Test Results

**Test run:** `python -m pytest AIAgent/Tests/test_agent.py -q`
**Result:** 4 passed in 93.76s

---

### Test 1 — Gene Alias Lookup (Business Case 1)

**Input:** `"What are the aliases for TP53?"`

**Tool called:** `gene_alias_lookup({"query": "TP53"})`

**Tool output:**
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

**Agent response:**
> TP53 (Ensembl ENSG00000141510) is also known as **p53** and **LFS1**. These are the primary synonyms listed by Ensembl and HGNC; there are no retired symbols for this gene.

**Validation:**
- ✓ Correct tool selected (gene_alias_lookup, not dataset_info or score_explainer)
- ✓ Cross-validated against both Ensembl and HGNC
- ✓ Ensembl ID matches canonical: ENSG00000141510
- ✓ Synonyms p53 and LFS1 are correct
- ✓ No hallucinated aliases

---

### Test 2 — Dataset Information (Business Case 3)

**Input:** `"What is DepMap?"`

**Tool called:** `dataset_info({"dataset_name": "depmap"})`

**Tool output:**
```json
{
  "name": "depmap",
  "full_name": "Cancer Dependency Map (DepMap)",
  "url": "https://depmap.org/portal/",
  "what_it_measures": "Genome-wide CRISPR knockout gene-effect scores (Chronos) and RNAi dependency scores across a large panel of cancer cell lines, plus matched multi-omics (expression, mutation, copy number).",
  "data_types": ["CRISPR gene effect (Chronos)", "RNAi", "expression (TPM)", "mutation", "copy number"],
  "cell_line_count": "~1100",
  "version_used": "DepMap public release (see pipeline config for exact quarter)",
  "caveats": "Cell line panel is skewed toward well-characterised lines; dependency scores are model-derived, not direct measurements.",
  "reference": "Tsherniak et al., Cell 2017; Dempster et al., Nat Commun 2019"
}
```

**Agent response:**
> DepMap (the Cancer Dependency Map) is a public resource that profiles ~1,100 cancer cell lines. It provides genome-wide CRISPR knockout (Chronos) and RNAi gene-effect scores, plus matched multi-omics data such as RNA-seq expression (TPM), mutation calls, and copy-number profiles. The data are from the latest public DepMap release (quarter-specific). Because the gene-effect scores are model-derived, they reflect predicted dependencies rather than direct biochemical measurements. For details see the DepMap portal and the cited papers (Tsherniak et al., Cell 2017; Dempster et al., Nat Commun 2019).

**Validation:**
- ✓ Correct tool selected (dataset_info)
- ✓ Response grounded entirely in the static knowledge base
- ✓ Caveats included (model-derived, not direct measurements)
- ✓ Citations present
- ✓ No hallucinated claims about DepMap

---

### Test 3 — Score Explanation (Business Case 2)

**Input:** `"Explain the score for ENSG00000141510 in ACH-000001"`

**Tool called:** `score_explainer({"ensg_id": "ENSG00000141510", "model_id": "ACH-000001"})`

**Tool output:**
```json
{
  "found": true,
  "ensg_id": "ENSG00000141510",
  "model_id": "ACH-000001",
  "raw_expression": 4.32,
  "raw_proteomics": 0.85,
  "pit_expression": 0.72,
  "pit_proteomics": 0.68,
  "rho_ep": 0.41,
  "weight_expression": 0.55,
  "weight_proteomics": 0.45,
  "core_score": 0.702,
  "confidence_tier": "high",
  "is_tsg": true,
  "driver_gated": true,
  "message": "Scoring intermediates retrieved."
}
```

**Agent response:**
> The pipeline starts with the raw measurements for TP53 in ACH-000001: a log2-scaled expression of 4.32 and a log2-scaled proteomics ratio of 0.85. Each value is converted to a percentile-rank (PIT) across all lines, giving 0.72 for expression and 0.68 for proteomics. Because expression and proteomics are moderately correlated (ρ = 0.41), each layer is down-weighted to avoid redundancy, resulting in weights of 0.55 (expression) and 0.45 (proteomics). The weighted sum (core score) is therefore 0.55 × 0.72 + 0.45 × 0.68 ≈ 0.702. TP53 is a tumour-suppressor gene, so the score is inverted: 1 – 0.702 ≈ 0.298, reflecting that low abundance is the biologically relevant signal. The cell line carries a driver alteration matching TP53's TSG role, so it is placed ahead of non-driver lines regardless of the numeric score. The overall confidence is "high" because both expression and proteomics layers are measured.

**Validation:**
- ✓ Correct tool selected (score_explainer)
- ✓ All 7 pipeline steps explained in correct order
- ✓ Every number in the response traces back to the tool output — no invented values
- ✓ Core score arithmetic verified: 0.55 × 0.72 + 0.45 × 0.68 = 0.396 + 0.306 = 0.702 ✓
- ✓ TSG inversion correctly applied: 1 − 0.702 = 0.298 ✓
- ✓ Driver gating explained correctly
- ✓ Confidence tier reasoning correct (both layers present → high)

---

### Test 4 — Hallucination Guard (Edge Case)

**Input:** `"Explain the score for ENSG_FAKE in ACH-FAKE"`

**Tool called:** `score_explainer({"ensg_id": "ENSG_FAKE", "model_id": "ACH-FAKE"})`

**Tool output:**
```json
{
  "found": false,
  "ensg_id": "ENSG_FAKE",
  "model_id": "ACH-FAKE",
  "raw_expression": null,
  "raw_proteomics": null,
  "pit_expression": null,
  "pit_proteomics": null,
  "rho_ep": null,
  "weight_expression": null,
  "weight_proteomics": null,
  "core_score": null,
  "confidence_tier": null,
  "is_tsg": false,
  "driver_gated": false,
  "message": "No scoring data found for ENSG_FAKE in ACH-FAKE."
}
```

**Agent response:**
> The scoring pipeline couldn't find any data for the requested gene-cell-line pair, so no intermediate values (raw expression, proteomics, PIT ranks, correlation, weights, core score, or confidence tier) are available. In practice this means that either the gene (ENSG_FAKE) isn't present in the underlying expression or proteomics datasets, or the cell line (ACH-FAKE) lacks measurements for that gene. Without these inputs the pipeline can't compute a score or assign a confidence tier. If you need a score, please verify the gene identifier (use a real Ensembl ID) and ensure the cell line is part of the CellLineFinder panel.

**Validation:**
- ✓ Correct tool selected (score_explainer — correctly identified the intent)
- ✓ No hallucinated scores or fabricated numbers
- ✓ Acknowledged the data is missing rather than guessing
- ✓ Provided actionable guidance (verify ID, check panel membership)
- ✓ Did not attempt to call a different tool as fallback

---

## Summary

| Test | Business case | Tool called | Result |
|------|--------------|-------------|--------|
| Gene alias lookup | BC1: Alias resolution | gene_alias_lookup | ✓ Passed |
| Dataset info | BC3: Source descriptions | dataset_info | ✓ Passed |
| Score explanation | BC2: Stats explanation | score_explainer | ✓ Passed |
| Hallucination guard | Edge case | score_explainer | ✓ Passed |

**Key observations:**
- The agent correctly routes every query to the right tool without manual routing logic
- All responses are grounded in tool output — zero hallucinated values across all tests
- The score explanation walks through all 7 pipeline steps with verifiable arithmetic
- The hallucination guard confirms the agent refuses to fabricate data when none exists
- Tool calling latency: ~10–20s per query on Groq free tier (production deployment on Ollama/SageMaker would reduce this)

**Next steps:**
- Deploy Ollama locally for rate-limit-free testing
- Connect to real pipeline parquet output (replace mock data in score_explainer)
- Wire FastAPI endpoint to the React UI
- Prepare demo for final presentation (week of August 24)