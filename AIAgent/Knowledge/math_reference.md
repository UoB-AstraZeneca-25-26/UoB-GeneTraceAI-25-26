# Scoring Pipeline — Math Reference

Condensed from `docs/MATH_REFERENCE.md` (Stage 2, 4, 5). Use these exact
definitions when explaining a `score_explainer` result — never invent numbers,
only interpret the ones returned by the tool.

## 1. Raw inputs
- Expression: `log2(TPM + 1)`.
- Proteomics: `log2` ratio.

## 2. PIT / percentile-rank normalisation
Each raw value is converted to a percentile rank **within that gene, across
all cell lines**, so expression and proteomics (different native scales)
become comparable on `[0, 1]`:
```
PIT(x) = rank(x) / n      # percentile rank, per gene, na_option="keep"
```
`pit_expression` and `pit_proteomics` in the tool output are these ranks.


## 3. Correlation-penalised weights
`rho_ep` is the Pearson correlation between expression and proteomics PIT
ranks. A highly-correlated layer is redundant, so it is down-weighted:
```
w_raw_E = 1 / (1 + |rho_ep|)
w_raw_P = 1 / (1 + |rho_ep|)
weight_expression, weight_proteomics = normalise(w_raw_E, w_raw_P)   # sum to 1
```

## 4. Core score
Weighted sum of the PIT ranks:
```
core_score = weight_expression * pit_expression + weight_proteomics * pit_proteomics
```

## 5. TSG inversion
For a **tumour-suppressor gene** (`is_tsg = true`), low abundance is the
biologically relevant signal (loss of the suppressor), so the score is
inverted before ranking:
```
core_score_TSG = 1.0 - core_score
```

## 6. Driver gating
`driver_gated = true` means the cell line carries a known driver alteration
(mutation, fusion, or copy-number event matching the gene's role — oncogene
amplification or TSG deletion). Driver-gated rows are sorted ahead of
non-driver rows regardless of `core_score`, for gene classes where a driver
event is the dominant signal.

## 7. Confidence tiers
Categorical, deterministic — not a numeric 0-1 score:
- **high** — measured regime, abundance-tracking class, both expression and
  proteomics layers present (`n_layers = 2`).
- **medium** — moderate/default case, or a CNA mechanistic hit upgrades the tier.
- **low** — activation-driven class ranked on score alone, no driver evidence.
- **insufficient** — regime not measured / too little data to score confidently.
