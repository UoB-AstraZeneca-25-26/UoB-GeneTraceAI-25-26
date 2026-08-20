# CellLineFinder Scoring Pipeline — Mathematical Reference

This document is embedded in the agent's system prompt. The
score_explainer tool has two modes. When it returns numeric values, walk
through each step using ONLY those numbers — never invent values. When it
returns null numerics, no scoring API is connected: explain the method
itself, and label any illustrative number as hypothetical rather than
presenting it as a result for the requested cell line.

---

## Pipeline overview (7 steps)

Input: one (gene, cell_line) pair.
Output: a ranked position with a confidence tier.

```
raw measurement → PIT normalisation → correlation-penalised weights
→ core_score → TSG inversion → driver gating → confidence tier
```

---

## Step 1 — Raw measurements

| Field | What it is | Units |
|-------|-----------|-------|
| `raw_expression` | Gene expression abundance | log₂(TPM + 1) |
| `raw_proteomics` | Protein abundance relative to reference | log₂ ratio |

These are on incompatible scales. A log₂(TPM+1) of 4.0 means ~15 TPM.
A log₂ ratio of 0.85 means protein is 2^0.85 ≈ 1.8× the reference.
Direct comparison is meaningless — Step 2 fixes this.

---

## Step 2 — Probability Integral Transform (PIT)

| Field | Formula | Range |
|-------|---------|-------|
| `pit_expression` | rank of this cell line among ALL cell lines for this gene, divided by (n+1) | [0, 1] |
| `pit_proteomics` | same, for proteomics | [0, 1] |

```
F̂ₙ(x) = rank(x) / (n + 1)      where n ≈ 950 cell lines
```

**What it means:** pit_expression = 0.72 means this cell line's expression
of this gene is higher than 72% of all ~950 cell lines in the panel.

**Key properties:**
- Monotone-invariant: any strictly increasing transform of the input
  gives the same percentile. This is why log₂-TPM and log₂-ratio are
  combinable — the PIT is scale-free.
- Marginally uniform on [0,1] by the probability integral transform.
- Finite-sample error: ±0.044 uniform 95% band at n ≈ 950 (DKW inequality).
  Any difference in core_score smaller than ~0.044 is indistinguishable
  from sampling error.

---

## Step 3 — Correlation-penalised weights

| Field | Formula |
|-------|---------|
| `rho_ep` | Spearman rank correlation between expression and proteomics PIT values across all cell lines for this gene |
| `weight_expression` | w_E = (1/(1 + \|ρ_EP\|)) / Σw_raw, normalised so weights sum to 1 |
| `weight_proteomics` | w_P = (1/(1 + \|ρ_EP\|)) / Σw_raw |

**Why penalise correlation:** Two highly correlated layers carry redundant
information. An idempotent weighted mean (A(x,x) = x) avoids inflating
the score when correlated layers agree — unlike Noisy-OR which inflates
by up to 0.25 at the midpoint.

**Critical structural fact:** For two exchangeable layers, the GLS
minimum-variance weights are (½, ½) regardless of ρ. Correlation does
not change the optimal weights — it changes the precision of the result:
Var = σ²(1+ρ)/2. The effective evidence count is:

```
n_eff = 2 / (1 + |ρ_EP|)
```

At ρ_EP = 0.46 (measured): n_eff = 1.37 independent layers of evidence.

---

## Step 4 — Core score

| Field | Formula |
|-------|---------|
| `core_score` | w_E × pit_expression + w_P × pit_proteomics |

For the two-layer case with equal weights: core_score = ½E + ½P.

**Properties:**
- Bounded in [0, 1] as a convex combination of values in [0, 1].
- Idempotent: if E = P = x, then core_score = x (agreement is neither
  rewarded nor penalised).
- Sensitivity: ∂core/∂E = w_E = 0.5, so a ±0.044 DKW error in E
  propagates to at most ±0.022 in core_score.

**Stratum-aware re-ranking:** The core_score is then re-ranked within
(gene, n_layers) strata because one-layer and two-layer scores have
different distributions (sd = 0.289 vs 0.247). A raw core_score of 0.7
is rarer in the two-layer stratum than in the one-layer stratum.
Re-ranking within stratum maps both back to U(0,1).

---

## Step 5 — Tumour suppressor gene (TSG) inversion

| Field | Effect when true |
|-------|-----------------|
| `is_tsg` | core_score_TSG = 1 − core_score; stratum_rank_TSG = 1 − stratum_rank |

**Why:** For a TSG, LOW abundance is the biologically relevant signal
(loss of the suppressor). The involution φ(x) = 1 − x is the unique
order-reversing affine involution on [0,1]. After inversion, all genes
can be sorted descending unconditionally.

Example: core_score = 0.702 → inverted = 0.298. Lower inverted values
now indicate stronger loss of the TSG.

---

## Step 6 — Driver gating

| Field | Meaning |
|-------|---------|
| `driver_gated` | True if the cell line carries a known driver alteration matching this gene's role |

**Driver flag formula:**
```
mut_driver        = any_driver ∨ oncogene_hit ∨ tsg_hit
fusion_driver     = (max_confidence == "high")
has_cna_alt       = (amplification ∧ role ∈ {oncogene, both})
                  ∨ (deletion ∧ role ∈ {tsg, both})
has_driver_alt    = mut_driver ∨ fusion_driver ∨ has_cna_alt
```

**Effect on ranking:** Driver-altered lines are placed AHEAD of
non-driver lines regardless of core_score, via a lexicographic sort:
```
sort_key = has_driver_alt × 10⁶ + core_score
```

**Directionality (LOF vs GOF):**
- LOF-altered driver lines (truncating mutation or deletion in TSG):
  ranked by ascending core_score (lower = better).
- GOF-altered driver lines (amplification or activating mutation in
  oncogene): ranked by descending core_score (higher = better).

---

## Step 7 — Confidence tier

| Field | Rule |
|-------|------|
| `confidence_tier` | Categorical: high, moderate, low, unknown |

**Assignment rules (applied in precedence order, last wins):**
```
"high"     — measured regime ∧ abundance_tracking ∧ n_layers == 2
"moderate" — default fallback
"low"      — activation_driven ∧ score_only (no driver evidence)
"unknown"  — regime not measured
CNA upgrade — mechanistic CNA hit present → promote tier
Chronos upgrade — chronos_check == "validated" → promote to moderate
```

**Not a probability.** The tiers are a categorical decision rule over
evidence configurations. No calibration curve P(correct | tier) exists.
The ordering is auditable and explicit rather than numerically precise.

---

## Explanation guidelines

When the score_explainer tool returns numeric values (API connected):
1. State the raw measurements and what they mean physically.
2. Explain PIT values as "higher than X% of all cell lines for this gene."
3. Show the weight calculation from rho_ep.
4. Verify the arithmetic: core_score = w_E × pit_E + w_P × pit_P.
5. If is_tsg: explain the inversion and compute 1 − core_score.
6. If driver_gated: explain that this line is promoted above non-driver lines.
7. State the confidence tier and what evidence produced it.
8. Note the ±0.044 DKW band.

When the score_explainer tool returns null numeric values (methodology mode):
1. Explain each of the 7 pipeline steps using the formulas above.
2. Use concrete hypothetical examples to illustrate each step
   (e.g. "if a cell line had pit_expression = 0.72, that would mean...").
3. Explain WHY each step exists — what problem it solves.
4. Do NOT invent specific values and present them as real results.
5. Make clear this is a methodology explanation, not a result for
   a specific cell line.
