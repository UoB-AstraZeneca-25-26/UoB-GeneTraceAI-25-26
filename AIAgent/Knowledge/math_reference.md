# CellLineFinder Scoring Pipeline — Mathematical Reference

This document is embedded in the agent's system prompt. The
score_explainer tool has two modes. When it returns numeric values, walk
through each step using ONLY those numbers — never invent values. When it
returns null numerics, no scoring API is connected: explain the method
itself, without any worked example.

The formulas below are taken from `final_pipeline/Scoring/core_score.py`,
`final_pipeline/utils/common.py` (`robust_z_matrix`, `stouffer_combine`)
and `final_pipeline/Scoring/confidence_tiers.py`. They are the algorithm
as implemented, not a simplified account of it.

---

## Pipeline overview (7 steps)

Input: one (gene, cell_line) pair.
Output: a ranked position with a confidence tier.

```
robust z-score → within-lineage inverse-normal transform
→ protein residualisation against RNA → evidence weights
→ Stouffer combination → core_score → gene-role direction
→ confidence tier
```

The two evidence arms are NOT symmetric. The RNA arm carries the
within-lineage rank position; the protein arm carries only the part of
the protein signal that RNA does not already explain.

---

## Step 1 — Is this level unusual for this gene?

Each raw measurement is expressed as how far it sits from the typical
cell line, using the median and MAD rather than the mean — so one
extreme line cannot distort the scale.

```
robust_z(x) = (x − median) / (MAD × 1.4826)
```

The 1.4826 scaling factor makes MAD comparable to a standard deviation
under normality. The denominator is floored at 0.384 so a gene with
almost no variation cannot produce an unbounded z-score, and a column
needs at least 5 peers within its (source, lineage) group to be scored
at all.

This is applied per source within lineage, then the three RNA sources
are combined by a source-count-weighted Stouffer sum into a single
`rna_z`; the two proteomics platforms are combined the same way into
`prot_z`. A pair measured by fewer than all three RNA sources is
shrunk toward zero by √(n_observed / n_total) — a deliberate discount
for low corroboration.

---

## Step 2 — Where does this line rank within its lineage?

Tissue of origin drives most expression differences, so the comparison
is made within the same lineage (e.g. among upper aerodigestive lines
only), not against the whole panel. The rank is converted to a normal
quantile, and the proportion is clipped just inside 0 and 1 so the
transform stays finite at both boundaries.

```
rna_pct     = R / n_lineage                    (rank within lineage)
rna_pct_z   = Φ⁻¹(clip(rna_pct, 0.001, 0.999))
```

A line ranked top of its lineage scores high here even if other tissues
express the gene more. The clip is what keeps the transform defined: at
rank n_lineage the proportion is exactly 1.0, and Φ⁻¹(1) = +∞, so the
value is clipped to 0.999 (giving ≈ +3.09) before the transform. The
same applies at the bottom of the lineage.

Because n_lineage is far smaller than the full panel (often tens of
lines, not ~950), the sampling error on this rank is correspondingly
larger than a whole-panel percentile would carry.

---

## Step 3 — Does the protein add anything the RNA didn't say?

Protein levels partly follow RNA levels. Counting both would count the
same evidence twice. Only the part of the protein signal that RNA does
not explain is carried forward — the residual after regressing protein
on RNA.

```
prot_resid = (prot_z − β_g · rna_z − median_resid) / SD_resid
```

This produces **e\***, the decorrelated protein evidence.

- `β_g` is the per-gene OLS slope, `β_g = ρ_g × SD(prot_z) / SD(rna_z)`.
  It is the regression coefficient, not the correlation: the two SDs
  differ substantially (measured 1.525 vs 0.951), so using ρ directly
  would overcorrect by roughly 1.9×.
- `ρ_g` is empirical-Bayes shrunk toward a prior of 0.353 with a
  pseudo-count of 30 shared lines, so genes with few shared lines lean
  on the prior and genes with many lean on their own data. Genes with
  fewer than 30 shared lines use the prior with a global SD ratio.
- The residual is then re-standardised per gene (subtract the median,
  divide by the SD of the residuals) so it has unit variance, which is
  what the weighting in Step 4 assumes. Note this step deliberately
  uses SD rather than MAD: unit variance is required here, and the
  robustness of Step 1's MAD scale is knowingly traded away.

Note the residualisation is against the raw `rna_z` from Step 1, not
against the within-lineage transform from Step 2.

---

## Step 4 — How much should each side count?

Three RNA sources that largely agree carry less independent information
than three unrelated ones. The weights reflect how much genuinely
independent evidence each side contributes.

```
W_RNA  = √1.431      W_PROT = √1.45
```

These are Kish effective sample sizes, not hand-tuned values:

- RNA: 3 sources at mean pairwise ρ = 0.548 → n_eff = 3 / (1 + 2×0.548) = 1.431
- Protein: 2 platforms at ProCAN–CCLE ρ = 0.373 → n_eff = 2 / (1 + 0.373) = 1.45

The 0.373 platform-agreement correlation is a different quantity from
the 0.353 RNA–protein prior used in Step 3; they are not
interchangeable.

---

## Step 5 — Combining the two into one score

The within-lineage RNA value from Step 2 and the RESIDUALISED protein
signal (**e\*** from Step 3, not raw `prot_z`) are combined using the
weights from Step 4, then mapped onto a 0–1 scale via the normal CDF.

```
core_z     = (W_RNA · rna_pct_z + W_PROT · e*) / √(W_RNA² + W_PROT²)
core_score = Φ(core_z)
```

**CRITICAL:** the protein term is `e*`, the output of Step 3, not raw
`prot_z`. If Step 5 used raw `prot_z`, the decorrelation in Step 3
would be wasted — the same RNA-correlated signal would be counted
twice. Likewise the RNA term is `rna_pct_z`, the output of Step 2, not
the raw `rna_z` of Step 1.

Not every pair has both arms. A gene with no proteomics coverage, or a
line with no protein data, falls back to `core_score = Φ(rna_pct_z)`
and is marked `n_layers = 1`. Pairs with both arms carry `n_layers = 2`.

`core_score` is then ranked within (gene, lineage) strata to give
`stratum_rank`.

---

## Step 6 — Does high or low count as good here?

For some targets high abundance supports selection; for others the
opposite. The direction depends on the gene's role — oncogene, tumour
suppressor, or unclassified.

```
If is_tsg:            score is inverted (1 − core_score)
If is_tsg is unknown: gene-role-blind ranking, no inversion applied
```

For a tumour suppressor, LOW abundance is the biologically relevant
signal (loss of the suppressor), so the ordering has to be reversed
before all genes can be sorted in one direction. For an oncogene or an
abundance-tracking gene the score is used as-is.

---

## Step 7 — How much evidence stands behind the score?

The tier is separate from the score. It reflects the breadth of
independent evidence, not the magnitude of the score.

```
tier = f(score_percentile_within_gene, has_driver_alteration, direction)
```

With a within-gene cutoff at the top 20%:

| Tier | Condition |
|------|-----------|
| HIGH | top 20% (direction-aware) AND a driver alteration is present |
| MEDIUM | top 20% AND no driver alteration |
| CONTEXT | not top 20% AND a driver alteration is present |
| LOW | not top 20% AND no driver alteration |

For a loss-of-function alteration in a tumour suppressor, "top" means
the BOTTOM 20% of the score distribution — low abundance is the
supporting evidence there. MEDIUM and LOW are assigned role-blind.

MEDIUM means the abundance evidence is strong but no genomic alteration
corroborates it. It is not a weaker score; it is a score with less
independent backing. The tiers are a categorical decision rule over
evidence configurations, not a calibrated probability — no
P(correct | tier) curve exists.

---

## Mapping the tool's fields onto these steps

The `score_explainer` payload uses the older field names. They bind to
the steps above as follows:

| Field | Step | Meaning |
|-------|------|---------|
| `raw_expression`, `raw_proteomics` | before Step 1 | raw measurements |
| `pit_expression` | Step 2 | the within-lineage transformed RNA value |
| `pit_proteomics` | Step 3 | the residualised protein value, e* |
| `rho_ep` | Step 3 | the shrunk per-gene RNA–protein correlation ρ_g |
| `weight_expression`, `weight_proteomics` | Step 4 | W_RNA, W_PROT |
| `core_score` | Step 5 | Φ(core_z) |
| `is_tsg` | Step 6 | gene role, drives the direction |
| `confidence_tier`, `driver_gated` | Step 7 | tier and driver evidence |

---

## Explanation guidelines

When the score_explainer tool returns numeric values (API connected):
1. State the raw measurements and what they mean physically.
2. Explain the Step 2 value as a rank position WITHIN THIS LINEAGE, not
   against the whole panel.
3. State the weights and where they come from (Kish n_eff).
4. Verify the arithmetic:
   core_z = (W_RNA × rna_pct_z + W_PROT × e*) / √(W_RNA² + W_PROT²).
5. If is_tsg: explain the direction and compute 1 − core_score.
6. If driver_gated: explain what the driver alteration contributes.
7. State the confidence tier and what evidence produced it.
8. Note that the within-lineage n is small, so the rank in Step 2
   carries meaningfully more sampling error than a whole-panel
   percentile would.

When the score_explainer tool returns null numeric values (methodology mode):
1. Explain each of the 7 pipeline steps using the formulas above.
2. Give NO worked example and NO illustrative numbers — not even ones
   labelled hypothetical. The UI renders these steps as the method
   applied to this specific pair, so a worked example reads as a result.
3. Explain WHY each step exists — what problem it solves.
4. Make clear this is a methodology explanation, not a result for a
   specific cell line.
