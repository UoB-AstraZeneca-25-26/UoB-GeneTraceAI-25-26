# CellLineFinder Scoring Pipeline — Mathematical Reference

This document is embedded in the agent's system prompt. The
score_explainer tool has two modes. When it returns numeric values, walk
through each step using ONLY those numbers — never invent values. When it
returns null numerics, no scoring API is connected: explain the method
itself, and label any illustrative number as hypothetical rather than
presenting it as a result for the requested cell line.

Verified 2026-08-21 against the live `final_pipeline/` source (not the
earlier `src/pipeline/` prototype some older docs describe) — the code that
actually produces the numbers behind this website. Canonical spec:
`docs/MATH_REFERENCE.md` §11 ("Final Pipeline: Mathematics and Audit").

---

## Pipeline overview

Input: one (gene, cell_line) pair.
Output: a `core_score` in [0,1], plus a separately-computed confidence tier.

```
per-source robust z-score  →  within-modality combination (RNA arm, Protein arm)
    →  cross-modality Stouffer combination  →  core_score
    →  (offline only) driver-alteration flag + confidence tier
```

The **ranking exposed by this website is core_score only** — sorted
descending, gene-role-blind. Tumour-suppressor direction and driver
alterations affect an *offline* confidence-tier label, not the score value
or the sort order a user sees. That distinction matters and is explained in
Steps 5–6 below — don't imply the live ranking re-orders on driver status
or inverts for TSGs, because it doesn't.

---

## Step 1 — Raw measurements → per-source robust z-scores

| Field | What it is | Units |
|-------|-----------|-------|
| `raw_expression` | RNA abundance, robust z-scored per source | log₂(TPM+1) → z |
| `raw_proteomics` | Protein abundance, robust z-scored per source | log₂ ratio → z |

Each source is standardised with the **median/MAD** z-score (not mean/SD),
robust to outliers:

```
robust_z(x) = (x − median(x)) / (MAD(x) × 1.4826)
```

RNA has 3 sources (DepMap, HPA, GEO); protein has 2 platforms (ProCAN,
CCLE). A pair missing some sources gets its z-score shrunk toward zero
(`z_t = z_raw / sqrt(N_SOURCES / n_observed)`) rather than dropped — most
protein pairs (~83%) have only 1 of 2 platforms, so this shrink is the
common case, not an edge case.

*Provenance:* `final_pipeline/utils/common.py:robust_z_matrix`,
`stouffer_combine`.

---

## Step 2 — Arm-level combination (RNA and Protein are NOT treated the same way)

| Field | Formula | Range |
|-------|---------|-------|
| `pit_expression` | RNA's rank **within lineage**, quantile-transformed | [0, 1] → z |
| `pit_proteomics` | Protein residual after regressing out RNA, then z-scored | z-score |

This is the step that most differs from a naive "percentile both layers"
description:
- **RNA** is percentile-ranked *within the cell line's lineage* (not
  globally across all ~950 lines), then converted to a normal quantile.
- **Protein is not percentile-ranked at all.** It's residualised against
  RNA first — `prot_resid = prot_z − β_g × rna_z`, where `β_g` is a
  per-gene regression slope shrunk toward a global prior via empirical
  Bayes (`RHO_PRIOR = 0.353`, shrinkage weight `K_RHO/(n+K_RHO)`) — then
  the residual is standardised `(x − median)/SD`. The intent: only the
  part of protein abundance RNA *doesn't already explain* contributes
  independent evidence.

**What it means (RNA arm):** pit_expression high means this line's
expression of this gene is high *relative to other lines of the same
lineage*, not the whole panel.

*Provenance:* `final_pipeline/Scoring/core_score.py` (`_percentile_within_lineage`,
the `prot_resid` block).

---

## Step 3 — Cross-modality weights (fixed constants, not a live correlation)

| Field | Value |
|-------|-------|
| `rho_ep` | Not computed live — no runtime RNA-vs-protein correlation feeds the weights |
| `weight_expression` | `W_RNA = √1.431 ≈ 1.196` |
| `weight_proteomics` | `W_PROT = √1.45 ≈ 1.204` |

These are **fixed constants**, each derived from agreement *within* one
modality, not between them:
- `W_RNA`: Kish effective-sample-size across the 3 RNA sources
  (`n_eff = 3/(1+2×ρ̄)`, measured `ρ̄ = 0.548` between DepMap/HPA/GEO) →
  `n_eff = 1.431`.
- `W_PROT`: same idea across the 2 protein platforms
  (measured cross-platform ρ = 0.373) → `n_eff = 1.45`.

The two arms end up almost exactly balanced (RNA 49.7% / Protein 50.3% of
variance) — RNA is marginally the *weaker* arm, the opposite of an earlier,
uncorrected version of this pipeline.

*Provenance:* `final_pipeline/Scoring/core_score.py:39-40`.

---

## Step 4 — Core score (Stouffer combination, not a weighted average)

| Field | Formula |
|-------|---------|
| `core_score` | `Φ(core_z)` where `core_z = (W_RNA·rna_pct_z + W_PROT·prot_resid_z) / √(W_RNA² + W_PROT²)` |

This is a **Stouffer Z-score combination**, not `w_E×E + w_P×P`: the two
arms are summed with fixed weights, normalised so `core_z ~ N(0,1)` under
independence, then mapped back to [0,1] via the standard normal CDF `Φ`
rather than staying a simple convex combination of two percentiles.

**Stratum-aware re-rank:** `core_score` is then re-ranked **within
(gene, lineage)**, not within (gene, evidence-layer-count) — lineage is
the stratification variable in this pipeline, so a line's rank always
accounts for which tissue type it's being compared against.

*Provenance:* `final_pipeline/Scoring/core_score.py` (Stouffer combination,
`stratum_rank` groupby).

---

## Step 5 — Gene-role direction (TSG vs oncogene): tier-only, NEVER a score inversion

| Field | Effect |
|-------|--------|
| `is_tsg` | Used only for the offline confidence-tier percentile threshold; does not change `core_score` or the ranking this website returns |

**THERE IS NO `1 − core_score` STEP IN THIS PIPELINE. Do not describe one,
even conditionally, even as "if TSG then invert."** This is the single
most common mistake to make when explaining this pipeline — the phrase
"TSG inversion" sounds like it should exist, and an *older, retired*
version of this pipeline did have exactly that (`core_score_TSG = 1 −
core_score`), but the code that serves this website was rewritten and no
longer does it. If you output a structured step for this stage, its
formula field must NOT contain `1 − core_score`, and its `status` MUST be
`"skipped"` for every gene, TSG or not — never `"active"` or `"inverted"`.

**What actually happens instead:** for a TSG, low abundance is the
biologically relevant signal, so the *offline* pipeline (not the live
API) flips which percentile band counts as "top" when assigning a
confidence tier to LOF-classified genes — a rank-*threshold* flip, not a
value inversion, computed once when `predictions_with_confidence.parquet`
is built.

**Say this plainly whenever gene role, TSG, or "inversion" comes up:**
the `core_score` and the ranking order this website returns are
gene-role-blind. Every gene, oncogene or tumour suppressor, is sorted by
the same unmodified `core_score`, descending. Only the separately-reported
`confidence_tier` field is TSG-aware — the score itself never is.

*Provenance:* `final_pipeline/Scoring/confidence_tiers.py`
(`is_lof_for_direction`, `directional_top`); confirmed absent from
`final_pipeline/api/ranking.py` (no `invert`, `1 -`, `tsg`, or `lof`
token anywhere in that file).

---

## Step 6 — Driver-alteration flag — metadata, not a ranking mechanism (here)

| Field | Meaning |
|-------|---------|
| `driver_gated` | Reported per line as `driver_alteration`; does not reorder this website's rankings |

**Driver flag formula (as implemented):**
```
mut_driver     = p_mutation ≥ 0.5
fusion_driver  = p_fusion ≥ 0.5
has_cna_alt    = (amplification ∧ role ∈ {oncogene, both})
              ∨ (deletion ∧ role ∈ {tsg, both})
has_driver_alt = mut_driver ∨ fusion_driver ∨ has_cna_alt
```

A `sort_key = has_driver_alt × 10⁶ + core_score` lexicographic promotion
does exist in the codebase (`final_pipeline/Ranking/cli.py`, an offline
analysis tool), but the deployed ranking API sorts by `core_score` alone
(`final_pipeline/api/ranking.py`, comment: "straight core_score ranking").
So on this website, `driver_alteration` is informational — it tells you
*why* a line might be biologically interesting, it does not move it up the
list.

*Provenance:* `final_pipeline/Scoring/driver_routing.py`,
`final_pipeline/api/ranking.py`.

---

## Step 7 — Confidence tier

| Field | Rule |
|-------|------|
| `confidence_tier` | Categorical: **HIGH, MEDIUM, CONTEXT, LOW** (no "UNKNOWN" tier exists) |

A 2×2 matrix of (top-20%-within-gene percentile) × (has a driver
alteration), with the percentile threshold flipped for LOF-classified
genes so "top" means low-expression for those:

```
HIGH    — top 20% (direction-aware) AND has driver alteration
MEDIUM  — top 20% AND no driver alteration
CONTEXT — not top 20% AND has driver alteration
LOW     — not top 20% AND no driver alteration
```

**Not a probability.** No calibration curve `P(correct | tier)` exists.
The ordering is auditable and explicit, not numerically precise.

*Provenance:* `final_pipeline/Scoring/confidence_tiers.py`
(`TOP_SCORE_PERCENTILE = 0.80`).

---

## Explanation guidelines

When the score_explainer tool returns numeric values (API connected):
1. State the raw per-source z-scores and what they mean physically.
2. Explain the RNA arm as a within-lineage percentile, and the protein arm
   as a residual after regressing out RNA — these are different
   transforms, not the same operation applied twice.
3. Give the fixed weights (W_RNA, W_PROT) — note they come from
   within-modality source agreement, not a live RNA-vs-protein
   correlation.
4. Show `core_score = Φ(core_z)`.
5. If asked whether TSG status or driver alterations affected the score
   or ranking shown: say clearly that they didn't — those only affect the
   separately-reported confidence tier.
6. State the confidence tier and which of the four evidence
   configurations produced it.

When the score_explainer tool returns null numeric values (methodology mode):
1. Explain each of the 7 pipeline steps above using the formulas given.
2. Use concrete hypothetical examples to illustrate each step.
3. Explain WHY each step exists — what problem it solves.
4. Do NOT invent specific values and present them as real results.
5. Make clear this is a methodology explanation, not a result for a
   specific cell line.
6. For the Step 5 entry specifically: never write a formula containing
   `1 − core_score`, and always set its `status` to `"skipped"` regardless
   of the gene's role. This applies to every gene in every response, not
   just ones you're unsure about — re-read Step 5 above before writing
   this entry if you're about to write "invert" or "TSG inversion".
