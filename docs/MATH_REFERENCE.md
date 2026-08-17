# GeneTraceAI — Mathematics Reference

Mathematical specification of the GeneTraceAI pipeline. For each quantity the
pipeline computes, this document states:

- **Model** — what latent quantity the formula estimates, and the probabilistic
  or algebraic object it is an instance of.
- **Assumptions** — what must hold for the estimator to be valid.
- **Properties** — structural behaviour: monotonicity, idempotence, invariance,
  limits, distribution of the output.
- **Error / limits** — finite-sample error, bias, or the conditions under which
  the estimator provably fails.
- **Provenance** — where it lives in the source, why it was introduced, what
  consumes it downstream.

Wherever a formula is an instance of a named estimator (GLS weights, the
probability integral transform, Kish's effective sample size, a hypergeometric
null, a logistic MLE, Fréchet–Hoeffding bounds), it is named as such, because
naming it imports its assumptions and error terms for free.

**If a term or quantity here is unfamiliar, it is defined in
[GLOSSARY.md](GLOSSARY.md)** — every biological term, assay and data column used
in this project, defined from first principles with its measured range in this
repository. This document assumes those definitions and concentrates on the
estimators. [BIOLOGICAL_FIDELITY.md](BIOLOGICAL_FIDELITY.md) covers what the
architecture claims biologically and where it falls short.

Line numbers are approximate (notebooks renumber cells on edit); function and
file names are exact. Paths are relative to the repo root
(`C:\Disertation\UoB-GeneTraceAI-25-26`).

---

**Stage numbering cross-reference** — `run_all.py` is the authoritative pipeline index.
Section headings in this document use an earlier numbering that predates the
introduction of the Altercations stage. The canonical mapping is:

| run_all.py | Directory | This document |
|------------|-----------|---------------|
| Stage 0 | `00_harmonisation` | Stage 0 — Harmonisation |
| Stage 1 | `01_Transcriptomics` | *(not covered — teammate code)* |
| Stage 2 | `02_Proteinomics` | *(not covered separately)* |
| Stage 3 | `03_Altercations` | *(see §§ driver-routing sections)* |
| Stage 4 | `Scoring/` | Stage 2 in this doc (core_score) |
| Stage 4b | `Scoring/driver_routing.py` | Stage 4 in this doc |
| Stage 4c | `Scoring/confidence_tiers.py` | Stage 5 in this doc |
| Stage 5 | `Ranking/` | Stage 7 in this doc |
| Stage 6 | `Validation/` | Stage 6 in this doc |

When referencing a stage fix, cite the directory name (`Scoring/core_score.py`) not the stage number to avoid ambiguity across documents.

---

A running list of what remains mathematically unjustified is kept in
[§9 Open mathematical gaps](#9--open-mathematical-gaps).

---

## Notation

| Symbol | Meaning |
|---|---|
| `g` | gene (`ensg_id`) |
| `c` | cell line (`model_id`) |
| `N` | number of cell lines scored for a given gene (≈ 950 in production) |
| `E, P, C` | normalised expression, proteomics, Chronos layers, each in [0,1] |
| `F, F̂ₙ` | true CDF and empirical CDF from `n` observations |
| `ρ` | correlation (Spearman unless stated); `ρ_EP` etc. for pairwise |
| `Σ, R` | covariance matrix, correlation matrix of the evidence layers |
| `𝟙` | vector of ones |
| `σ(·)` | logistic function `1/(1+e^{−·})` |
| `m` | number of hypotheses in a multiple-testing family |

---

# §0 — Preliminaries: shared mathematical objects

Five objects recur across every stage. They are defined once here; later
sections reference rather than repeat them.

## 0.1 The empirical CDF and the probability integral transform

Percentile-rank normalisation `x.rank(pct=True)` is the empirical CDF

```
F̂ₙ(x) = (1/n) · Σᵢ 𝟙[Xᵢ ≤ x]
```

**Probability integral transform (PIT).** If `X` has continuous CDF `F`, then
`F(X) ~ Uniform(0,1)` exactly. Consequently every percentile-normalised layer in
this pipeline is *marginally uniform on [0,1] by construction*, not
approximately so. Three consequences are used repeatedly:

**(a) Monotone invariance.** For any strictly increasing `φ`,
`F̂ₙ^{φ(X)}(φ(x)) = F̂ₙ^X(x)`. The score is unchanged by any monotone
reparameterisation of the input. This — not "robustness to outliers" — is the
rigorous justification for combining log₂-TPM expression with proteomic
log-ratios: the transform is scale-free, so the two assays need not share units,
a distributional family, or even a linear relationship, only a common ordering
semantics.

**(b) Uniform finite-sample error (DKW inequality).**

```
P( sup_x |F̂ₙ(x) − F(x)| > ε )  ≤  2·exp(−2nε²)
```

Inverting at 95% gives a uniform confidence band `ε = √(ln(2/α)/(2n))`:

| n (cell lines) | 95% uniform band on every percentile |
|---|---|
| 500 | ±0.061 |
| 950 (production) | **±0.044** |
| 1000 | ±0.043 |

Every `E` and `P` value in this pipeline therefore carries a ±0.044 band, and
any downstream difference in `core_score` smaller than roughly that magnitude is
not distinguishable from sampling error in the rank estimate. This band is
*uniform* (simultaneous over all x), so it is conservative for any single
percentile but correct for statements about the whole score column.

**(c) Ties break exact uniformity.** PIT requires `F` continuous. Sparse
modalities (proteomics with many identical imputed values; `fusion_count` with
87% of mass at a single value) produce ties, which `method="average"` maps to a
shared mid-rank. The transformed variable is then uniform only on the support of
distinct values, and its distribution has atoms. Where tie mass is large this is
noted at the relevant section (§C.3, §2.6).

## 0.2 Aggregation operators on [0,1]

An **aggregation operator** is `A: [0,1]ⁿ → [0,1]` satisfying the boundary
conditions `A(0,…,0)=0`, `A(1,…,1)=1`, and monotone non-decreasing in each
argument. All four of this pipeline's combination rules are aggregation
operators. They differ in one property that turns out to be the whole story:

**Idempotence.** `A` is idempotent iff `A(x,x,…,x) = x`.

| Operator | Formula | Idempotent? | `A(x,x) − x` |
|---|---|---|---|
| Noisy-OR (probabilistic sum) | `a + b − ab` | **No** | `x(1−x)`, max **0.25** at `x=0.5` |
| Weighted mean, `Σwᵢ = 1` | `Σ wᵢ xᵢ` | **Yes** | `0` |
| Product | `ab` | No | `−x(1−x)` |
| `product_floor` | `max(ab, ½·min(a,b))` | No | `−x(1−x)` for `x > ½` |

**This is the mathematical content of the Stage 2.3 → 2.4 change.** When two
evidence layers are strongly correlated they carry near-identical values; an
idempotent operator returns that value unchanged, while a non-idempotent
t-conorm inflates it by up to 0.25 in the middle of the range. See §2.3 for the
derivation and the matching production evidence.

**t-norms, t-conorms, and the Fréchet–Hoeffding bounds.** For two events with
marginal probabilities `p, q`, the probability of their union is *not*
determined by the marginals; it is bounded:

```
max(p, q)  ≤  P(A ∪ B)  ≤  min(1, p + q)
   ↑                              ↑
comonotone                   counter-monotone
(perfect positive dep.)      (maximal negative dep.)
```

Noisy-OR `p + q − pq` sits strictly inside this interval and is exactly the
**independence point**. It is easy to check that `p+q−pq − max(p,q) = min(p,q)·(1−max(p,q)) ≥ 0`,
so noisy-OR always exceeds the comonotone lower bound. As dependence increases
toward comonotonicity the truth moves toward `max(p,q)`; noisy-OR does not move,
which is precisely the inflation.

This single inequality unifies three otherwise-unrelated formulas in this
codebase: `noisy_or` is the independence point, the `max`-like behaviour is the
comonotone endpoint, and `product_floor` is an ad hoc operator whose position
inside the interval is unanalysed (see §9).

## 0.3 Combining correlated estimators: GLS weights and effective sample size

Suppose each evidence layer `xᵢ` is an unbiased estimator of a common latent
quantity `θ`, with `Var(x) = Σ`. Among all linear combinations `wᵀx` with
`wᵀ𝟙 = 1` (i.e. all unbiased ones), the minimum-variance choice is the
generalised-least-squares / minimum-variance-portfolio solution:

```
w* = Σ⁻¹𝟙 / (𝟙ᵀΣ⁻¹𝟙)          Var(w*ᵀx) = 1 / (𝟙ᵀΣ⁻¹𝟙)
```

**Symmetric two-layer case.** With `Σ = σ²[[1,ρ],[ρ,1]]`:

```
Σ⁻¹𝟙 ∝ (1−ρ, 1−ρ)   ⇒   w* = (½, ½)   for every ρ ∈ (−1,1)
Var(w*ᵀx) = σ²(1+ρ)/2
```

**Equal weights are optimal regardless of ρ.** Correlation between two
exchangeable layers does not change the optimal weights at all — it changes the
*precision* of the result. This is the key structural fact behind §2.4.

**Effective number of independent layers (Kish).**

```
n_eff = n² / (𝟙ᵀR𝟙) = n / (1 + (n−1)·ρ̄)
```

Evaluated on this project's measured correlations
(`ρ_EP = 0.46`, `ρ_EC = ρ_PC = 0.005`):

| Configuration | nominal layers | **n_eff** | GLS variance (σ² units) |
|---|---|---|---|
| E + P | 2 | **1.37** | 0.730 |
| E + P + C | 3 | **2.28** | 0.424 |

Two layers of expression + proteomics carry 1.37 layers' worth of independent
evidence, not 2. This gives `n_layers` — used throughout the pipeline as a
stratification key and a confidence-tier input — a quantitative definition
rather than a bookkeeping one.

*Caveat on ρ.* `ρ_EP = 0.46` is a literature estimate of the mRNA–protein
correlation on native scales, whereas the layers being combined are
percentile-transformed. Spearman ρ is invariant to the monotone transform but
Pearson ρ is not, so the number is only exactly the right one if the literature
figure is a rank correlation. This is an acknowledged approximation (§9).

## 0.4 Rank correlation as an estimator

Spearman's ρ is Pearson's correlation computed on ranks:

```
r₁ = rank(x),  r₂ = rank(y)          (average ranks for ties)
ρ = Σ(r₁−r̄₁)(r₂−r̄₂) / √( Σ(r₁−r̄₁)² · Σ(r₂−r̄₂)² )
```

**Ties matter.** The familiar shortcut `ρ = 1 − 6Σd²/(n(n²−1))` is valid *only*
without ties. Because this pipeline's inputs are percentile ranks of frequently
tied data, the Pearson-on-ranks form above is the correct one and is what both
`scipy.stats.spearmanr` and the manual implementation in
`build_chronos_validation.py` compute.

**Null distribution.** Under H₀ (independence), the exact null is the
permutation distribution over all `n!` rank pairings, with

```
E[ρ] = 0,   Var(ρ) = 1/(n−1),   so   √(n−1)·ρ  →  N(0,1)
```

The code instead uses the t-approximation

```
t = ρ·√((n−2)/(1−ρ²)),   p = 2·(1 − T_cdf(|t|, df = n−2))
```

which is *exact* for Pearson's r under bivariate normality and only
**asymptotic** for Spearman's ρ. It is adequate for `n ≳ 20–30`, which is the
mathematical justification for `MIN_N = 30` in `build_chronos_validation.py` —
the threshold is a validity condition on the test, not an arbitrary sample-size
preference.

**Variance stabilisation.** `ρ` is bounded and its sampling distribution is
skewed near ±1. Fisher's transform `z = artanh(ρ)` is variance-stabilising:

```
Var(z) ≈ 1/(n−3)          (Pearson)
Var(z) ≈ 1.06/(n−3)       (Spearman, Fieller correction)
```

| n | se(z) |
|---|---|
| 30 | 0.193 |
| 100 | 0.102 |
| 900 | 0.033 |

Confidence intervals should be built on `z` and back-transformed with `tanh`,
which guarantees the interval stays inside [−1,1] and is not skewed.

## 0.5 The bootstrap

**Percentile method** (used in `validate_chronos.bootstrap_rho`): resample the
unit of exchangeability `B` times, recompute the statistic, take the empirical
2.5th/97.5th percentiles. Validity requires the statistic to be approximately
*pivotal* — its sampling distribution roughly independent of the unknown
parameter. For a correlation this holds only after the Fisher transform of §0.4,
which is why percentile intervals on raw ρ are mildly anti-conservative near the
boundaries.

**Exchangeability assumption.** The resampling unit here is the *cell line*. The
interval is valid only if cell lines are exchangeable draws from a common
population. They are not exactly: cell lines cluster by lineage, which the
lineage-confound tests (§C.6) themselves demonstrate. Clustered resampling (draw
lineages, then lines within lineage) would be the assumption-consistent version.

**Paired bootstrap on a difference.** Resampling the same rows for both the
baseline and candidate score induces positive correlation between `ρ_base` and
`ρ_cand`, so

```
Var(Δ) = Var(ρ_cand) + Var(ρ_base) − 2·Cov(ρ_cand, ρ_base)
```

is strictly smaller than the unpaired variance whenever the covariance is
positive. This is why the paired design in the confound tests is more powerful
than comparing two independent intervals, and why "the two CIs overlap" is not
the correct test for "the scores differ".

**Monte Carlo error.** `B` controls only simulation noise, not statistical
uncertainty. `B = 2000` (confound tests) is adequate for a 95% interval;
`B = 1000` (`bootstrap_rho`) is at the low end of Efron's recommendation for
interval endpoints, though fine for a standard error.

**Better intervals available.** BCa corrects both median bias and skewness and
is the standard recommendation where the percentile method is used here.

## 0.6 Multiplicity

Several stages apply a per-gene hypothesis test at `α = 0.05` across thousands
of genes. Under the global null, `α·m` of them clear the threshold by
construction. The family-wise picture:

```
Benjamini–Hochberg:  order p₍₁₎ ≤ … ≤ p₍ₘ₎;  reject H₍ᵢ₎ for all i ≤ max{ i : p₍ᵢ₎ ≤ i·q/m }
                     controls FDR ≤ q under independence or PRDS
Benjamini–Yekutieli: same with threshold i·q / (m·Σⱼ₌₁ᵐ 1/j)   — valid under arbitrary dependence
Storey π₀ estimate:  π̂₀ = #{pᵢ > λ} / (m(1−λ))   — estimates the null fraction directly
```

Per-gene tests against a *shared* ground-truth vector are positively dependent,
which satisfies PRDS, so plain BH is the appropriate procedure.

**As of this writing no multiple-testing correction exists anywhere in the
repository** (verified by search over `src/`). See §C.2 and §9 for where this
bites.

---

# Stage 0 — Harmonisation

`src/pipeline/00_harmonisation.ipynb`, `src/new-arch/00_data_harmonisation_pipeline.ipynb`,
`src/Track - B/TrackB_final_resolution_pipeline.ipynb`

No estimation at this stage — identity resolution plus data-quality auditing.
The quantitative operations are descriptive statistics that gate later modelling
choices.

**Missing-data audit**
```
missing_pct = df[col].isna().mean() · 100
```
This is the sample mean of a Bernoulli indicator; its standard error is
`√(p(1−p)/n)`, negligible at these table sizes, so the percentages can be read
as exact. Columns at 100% missing are dropped.

*Provenance:* computed per (table, column) across all 15 pipeline tables. Gates
the per-modality normalisation choice made in Stage 2.

**Numeric scale summary** — min/max/mean/median/std per numeric column.
Diagnostic: distinguishes raw counts, log-ratios, z-scores and percentages so
that Stage 2 can choose a normalisation. Note that mean/std are only meaningful
descriptors for roughly symmetric columns; for the heavily right-skewed columns
here (variant counts, FFPM) the median/IQR pair is the informative one, which is
part of the reason the pipeline settles on rank-based normalisation.

**Modality coverage count**
```
n_modalities = Σⱼ coverage_boolⱼ,        complete_cycle = (n_modalities == 8)
```
Diagnostic only — not an input to any score. Note this counts *nominal*
modalities; the evidence-weighted analogue is `n_eff` (§0.3).

**Track B identity-bridge QA**: set intersection/union counts and match rates
(`mapped / before · 100`) at each join of the RRID↔CVCL↔ACH↔profile crosswalk.
Gates accept / "disputed" / "orphaned" — no numeric scoring.

**Name matching — no metric is involved.** Matching is exact-match after
deterministic string normalisation (`normalise_base()` in Track B) or substring
containment (`fuzzy_model_id_for()` in `02_core_score.ipynb`, via
`.str.contains()`). Substring containment is *not* a metric: it is not
symmetric, and it fails the triangle inequality. "Fuzzy" here means permissive
containment, not a Levenshtein or Jaccard distance. No edit-distance or
set-similarity score exists anywhere in this pipeline.

---

# Stage 1 — Lookup / Metadata Join

`src/pipeline/01_lookup_metadata_join.ipynb`

Pure additive left-join. One QA gate:
```
pct = 100 · matched / len(enriched)
flag "<-- INVESTIGATE KEY MISMATCH"  if  pct < 90
```
*Provenance / rationale:* low match rates indicate key mismatch (case, version
suffix, wrong grain), not genuinely missing metadata. The 90% figure is an
operational tripwire, not an estimated threshold.

---

# Stage 2 — Core Score

`src/pipeline/02_core_score.ipynb`, `src/pipeline/scoring_variants.py`

The mathematical centre of the pipeline: mapping raw expression and proteomics
measurements to one comparable `core_score` per (gene, cell line).

## 2.1 Percentile-rank normalisation

```
E = expr.rank(axis=0, pct=True, na_option="keep")     # per gene, across cell lines
P = prot.rank(axis=0, pct=True, na_option="keep")
```

**Model.** The empirical CDF of each gene's abundance across cell lines; by PIT
(§0.1) an estimate of `F_g(x)`, the population quantile of that cell line within
that gene's abundance distribution. The score is explicitly a *within-gene,
across-cell-line* quantile — it says nothing about absolute abundance and
nothing about comparisons across genes.

**Assumptions.** (i) Continuous underlying distribution — violated by ties
(§0.1c). (ii) Cell lines are an i.i.d. sample from a population of interest —
violated by lineage clustering, which biases the quantile toward whichever
lineages are over-represented in the panel.

**Properties.** Marginally uniform on [0,1]. Invariant under any strictly
increasing transform of the input (§0.1a) — this is what makes log₂-TPM and
proteomic log-ratios combinable without a units argument. Bounded, so no outlier
can dominate: an arbitrarily extreme value receives rank `n`, exactly as the
second-largest would.

**Error / limits.** ±0.044 uniform 95% band at n ≈ 950 (§0.1b). *Information
loss:* rank normalisation discards all magnitude information within a gene — two
cell lines differing 100-fold and two differing 1.01-fold receive the same
adjacent ranks. Where effect magnitude matters biologically this is a real loss,
traded deliberately for cross-assay comparability.

*Provenance:* `scoring_variants.normalise(method="percentile")`
([scoring_variants.py:30-31](../src/pipeline/scoring_variants.py)); `percentile_rank()` in
`02_core_score.ipynb`. Applied independently to the deduplicated expression
matrix and the proteomics matrix (both `model_id × ensg_id` wide) before
combination. An internal audit found min-max scaling (§2.5) weakest by
comparison.

**A monotone-invariance corollary that constrains every evaluation of this
layer.** Within one gene the percentile is a strictly increasing transform of the
raw `log2(TPM+1)` value, and every rank-based metric is invariant under such a
transform (§0.1a). Therefore

```
AUROC(percentile) = AUROC(raw)      pAUC(percentile) = pAUC(raw)      exactly
```

The percentile can neither gain nor lose against raw expression on any ROC-type
readout. Its *only* claim over the raw value is cross-gene comparability. Any
experiment that appears to compare "percentile vs expression" on AUROC is
comparing a quantity against itself; the only meaningful contrast is against a
score computed on a **different axis** — per-sample rather than per-gene (§C.16).

**Empirical status — BOTH FINDINGS BELOW ARE WITHDRAWN.** They were measured
before `gate_audit/` ran. §C.15 and §C.16 have been moved to
[MATH_TRIAL_AND_ERROR.md](MATH_TRIAL_AND_ERROR.md); the numbers are retained
there with the reason for retirement. A placebo with the dependency labels
randomised reproduces 98% of the depletion, and the residual is CRISPR-specific
(flat on drug, opposite on RNAi). Do not cite the pAUC figures below.

~~Two findings now bear directly on this section and are recorded in full at
§C.15 and §C.16:~~

1. The layer functions as a **gate, not a ranker**. Against three independent
   label sources the gate region (FPR 0.8–1.0) is well above chance while the
   ranker region (FPR 0.0–0.2) is at chance: `0.596 / 0.500` (Chronos),
   `0.593 / 0.501` (Project Score), `0.601 / 0.518` (GDSC). Dependencies are
   depleted at the bottom of the score and *not* enriched at the top. This does
   not change the formula; it changes what the formula is entitled to claim
   (§C.15).
2. The percentile is **not** an expensive encoding of a detection boolean. Among
   cell lines that an absolute per-line detection call rates *expressed*, the
   bottom decile of the percentile still carries a dependency ratio of
   **0.810–0.827** (spread 0.017 across six calibration rules; 0.790–0.799 on
   paralog-singleton genes only). It clears a 10% effect floor in 6/6 rules. So
   the percentile carries graded information beyond "on or off" (§C.16).

The corresponding negative result also belongs here: a per-sample absolute score
does **not** beat the percentile. Its apparent +0.030 pAUC advantage was entirely
selection optimism, and vanishes under half-sample cross-validation (§C.16).

## 2.2 Duplicate-profile deduplication

```
expr_dedup = expr.groupby(level=0).mean()
```

**Model.** Treats the two RNA profiles of a cell line as replicate measurements
of one latent abundance with i.i.d. noise; the sample mean is then the MLE under
Gaussian noise and the minimum-variance unbiased estimator, halving the noise
variance (`σ²/2`).

**Properties.** Mean is the only one of {mean, max, drop} that is both unbiased
and idempotent here. `max` of two replicates has `E[max] > θ` — a positive bias
of `σ/√π` under Gaussian noise — inflating exactly the 16 duplicated cell lines
relative to all others. `drop` discards half the information.

**Limits.** Averaging before ranking, rather than ranking before averaging, is
the right order only if the replicates share a scale; they do here (same assay,
same units).

*Provenance:* `dedup_mean()`, `02_core_score.ipynb`. 16 ACH cell lines have two
RNA profiles each.

## 2.4 Correlation-penalised weighted sum (production formula)

**Implementation as written:**

```
mean_ρ_E = (|ρ_EP| + |ρ_EC|)/2      w_E_raw = 1/(1 + mean_ρ_E)
mean_ρ_P = (|ρ_EP| + |ρ_PC|)/2      w_P_raw = 1/(1 + mean_ρ_P)
mean_ρ_C = (|ρ_EC| + |ρ_PC|)/2      w_C_raw = 1/(1 + mean_ρ_C)

w = normalise(w_raw)                 so Σw = 1
core_score = w_E·E + w_P·P + w_C·C
```

Two- and one-layer subsets renormalise over the layers present
(e.g. `w_EP_E = w_E_raw/(w_E_raw + w_P_raw)`). Defaults: `ρ_EP = 0.46`
(literature + internal validation), `ρ_EC = ρ_PC = 0.005` (measured, §2.7).

**Model.** This is a heuristic approximation to the GLS minimum-variance
combination of §0.3, `w* = Σ⁻¹𝟙/(𝟙ᵀΣ⁻¹𝟙)`, under the assumption that each layer
is an unbiased estimator of a common latent abundance with equal variance.
Evaluated on this project's own constants:

| Weights | E | P | C |
|---|---|---|---|
| **GLS optimum** `Σ⁻¹𝟙/(𝟙ᵀΣ⁻¹𝟙)` | 0.2892 | 0.2892 | **0.4215** |
| Implemented `1/(1+mean|ρ|)` heuristic | 0.3099 | 0.3099 | 0.3801 |

The heuristic lands within 0.042 of the optimum and orders the layers correctly
(it up-weights the near-uncorrelated Chronos layer), but it is an approximation,
and should be described as one rather than as a derived result. The exact
solution is a one-line matrix computation on the same inputs.

**The two-layer case is the important one — and ρ cancels.** From §0.3, for two
exchangeable layers the GLS weights are `(½, ½)` for *every* ρ. The
implementation reproduces this: `mean_ρ_E = mean_ρ_P` by symmetry, so
`w_E_raw = w_P_raw`, so after normalisation `w_E = w_P = ½` identically. Since
`override_mask = both_layers_present & weight_available` restricts the override
to the two-layer stratum, **the correlation term has no effect whatsoever on any
production score.** The formula reduces exactly to `core_score = ½E + ½P`.

This is not a defect — `(½,½)` is optimal — but it relocates the explanation.
What fixed the Noisy-OR inflation was **replacing a non-idempotent operator with
an idempotent one** (§0.2), not the correlation weighting. The correct claim is:

> Two correlated layers should be combined by an idempotent operator, because
> agreement between dependent sources is not corroboration. Among idempotent
> linear operators the minimum-variance choice for two exchangeable layers is
> the equal-weighted mean, independent of ρ. Correlation does not change the
> weights; it changes the precision of the result, `Var = σ²(1+ρ)/2 = 0.73σ²`,
> and the effective evidence count, `n_eff = 1.37` (§0.3).

**Properties.** Idempotent (`A(x,x) = x`), so agreement is neither rewarded nor
penalised. Bounded in [0,1] as a convex combination of values in [0,1]. Linear,
hence continuous with bounded sensitivity `∂core/∂E = w_E = ½` — a ±0.044 DKW
error in `E` propagates to at most ±0.022 in the two-layer score, versus up to
±0.044 under noisy-OR where `∂N/∂E = 1−P` can reach 1.

**Error / limits.** The variance reduction (`0.73σ²` vs `σ²`) assumes equal
per-layer variance; if proteomics is materially noisier than expression, the
GLS weights are no longer `(½,½)` and the correct weights are
`w ∝ Σ⁻¹𝟙` with unequal diagonal. Per-layer variances have not been estimated
(§9).

**Verified on production data:** the 22.5% stratum scoring > 0.95 drops to
**2.3%**, and the 379 rows previously pinned at 1.0 score **< 0.5** — consistent
with the idempotence argument, since a mean of two mid-range layers cannot reach
1 unless both do.

*Provenance:* `scoring_variants.combine_three()`
([scoring_variants.py:75-136](../src/pipeline/scoring_variants.py)) is the canonical
implementation; applied as a post-hoc override in `02_core_score.ipynb`.
One-layer fallback cells are untouched (nothing to weight). This is the value
written to `outputs/core_score.parquet` and consumed by driver-gated routing,
confidence tiers, held-out evaluation, and the explain/rank/CLI layer.

## 2.6 Stratum-aware rank

```
stratum_rank = core_score.groupby(["ensg_id","n_layers"]).rank(pct=True, method="average")
```

**Model.** The empirical CDF of `core_score` *within* the (gene, n_layers) cell,
i.e. a second application of PIT conditional on evidence coverage.

**Why this is mandatory, derived rather than asserted.** By §0.1, a one-layer
score is exactly `U(0,1)`, with

```
sd(1-layer) = √(1/12) = 0.2887
```

A two-layer score is the mean of two correlated uniforms. For any two variables
with equal variance `σ²` and correlation ρ:

```
Var( (X+Y)/2 ) = (2σ² + 2ρσ²)/4 = σ²(1+ρ)/2
sd(2-layer) = √( (1+0.46)/2 · 1/12 ) = 0.2466
```

i.e. a **shrinkage factor of `√((1+ρ)/2) = 0.854`**. Under independence the
two-layer score would be exactly triangular (Irwin–Hall, `n=2`) on [0,1] with
density `4x` on `[0,½]` and `4(1−x)` on `[½,1]`; positive ρ makes it heavier at
the extremes than triangular but still strictly more concentrated at ½ than
uniform.

So the two strata do not merely rest on different amounts of evidence — **they
are drawn from different distributions with different variances**, and a raw
`core_score = 0.7` is a rarer event in the two-layer stratum than in the
one-layer stratum. Comparing them directly is a distributional error, not a
judgement call. Re-ranking within stratum maps both back to `U(0,1)` and makes
them comparable, which is exactly what `stratum_rank` does.

**Empirical confirmation on production data.** Measured directly from
`src/pipeline/outputs/core_score.parquet` (28,403,110 rows):

| stratum | rows | share | measured sd | predicted sd |
|---|---|---|---|---|
| 1-layer | 25,672,370 | 90.39% | **0.2751** | 0.2887 (uniform) |
| 2-layer | 2,730,740 | 9.61% | **0.2418** | 0.2466 (at ρ = 0.46) |

The measured shrinkage ratio is `s₂/s₁ = 0.879` against a predicted
`√((1+ρ)/2) = 0.854`. Inverting the relation, the observed ratio implies

```
ρ̂ = 2·(s₂/s₁)² − 1 = 0.544
```

**an independent estimate of the mRNA–protein correlation recovered purely from
the variance structure of the output**, with no reference to the literature value
or to the `ρ_EP = 0.46` constant hard-coded in `combine_three()`. It falls inside
the 0.46–0.58 range reported by Nusinow et al. (2020). This is a genuine
out-of-sample validation of the derivation: the theory predicts a specific
variance ratio, the data exhibits it, and the implied parameter matches an
independently published quantity.

The 1-layer sd falls slightly below the uniform ideal (0.2751 vs 0.2887) because
ties compress the percentile transform (§0.1c) and because per-gene sample sizes
vary.

**Limits.** Within-stratum ranking makes the *marginals* comparable but does not
make them *calibrated* against each other: a 0.9 in the one-layer stratum and a
0.9 in the two-layer stratum are equally rare within their strata but rest on
`n_eff = 1` versus `n_eff = 1.37` of evidence. Cross-stratum calibration is
explicitly deferred pending a held-out regression against CRISPR/GDSC anchors.

*Provenance:* `scoring_variants.score()`
([scoring_variants.py:172-175](../src/pipeline/scoring_variants.py)); `02_core_score.ipynb`.
Surfaced to end users as `stratum_rank_percentile` in `explain_pair.py` /
`rank_cell_lines.py` — the safe metric for comparing cell lines with differing
coverage.

# Stage 3 — Cell-line Name Resolution Audit

`02_core_score.ipynb` (audit cell)

Compares the deterministic RRID→CVCL→ACH resolver (`model_id_for()`) against
naive substring matching (`fuzzy_model_id_for()`) across 30 test names,
classifying each into
`OK / RRID_RESCUED / WRONG / FUZZY_ONLY / AMBIGUOUS_BUT_AGREE / BOTH_MISSING`.
Result: 18 rescued, 9 OK, 1 wrong, 1 both-missing, 1 fuzzy-only.

**Statistical note.** With `n = 30` the observed 1/30 error rate has a
Clopper–Pearson 95% interval of roughly [0.001, 0.172] — i.e. the audit is
consistent with a true error rate anywhere up to ~17%. The audit demonstrates
that the resolver rescues many cases; it does not establish a tight bound on its
error rate. Diagnostic only; not consumed by scoring, and no numeric similarity
metric is involved (§Stage 0, name matching).

---

# Stage 4 — Driver-Gated Routing

`src/pipeline/04_driver_gated_routing.ipynb`, `src/pipeline/build_cna_layer.py`

## 4.1 Driver flag (Boolean lattice, no numeric estimator)

```
mut_driver            = any_driver ∨ oncogene_hit ∨ tsg_hit
fusion_driver         = (max_confidence == "high")
has_driver_alteration = mut_driver ∨ fusion_driver
```

**Model.** A monotone Boolean function (a disjunction — the join of the free
distributive lattice on those literals). Monotonicity means adding evidence can
only turn the gate on, never off, which is the desired semantics for "is there a
driver event."

**Why a Boolean gate rather than the continuous `p_mutation`.** The continuous
score would put every cell line on a graded scale where small differences in an
uncalibrated Hill output (§C.1) reorder cell lines. Thresholding to a Boolean
discards magnitude but eliminates that sensitivity — a deliberate
bias–variance trade in favour of variance, quantified in §4.3. It replaced the
earlier `p_mutation > 0.5` "any alteration" gate, which added variance without
improving hit rate.

*Provenance:* used as the primary sort key ahead of `core_score` for
`activation_driven` and `unknown` gene classes:
```
sort_values(["has_driver_alteration","core_score"], ascending=[False, False])
```

## 4.2 Copy-number amplification / deletion thresholds

```
AMP_THRESHOLD = 2.5     # TOTAL_CN > 2.5 → amplification
DEL_THRESHOLD = 1.5     # TOTAL_CN < 1.5 → deletion
```

**Model.** A symmetric ±0.5 band around the diploid baseline `CN = 2`. These are
*conventional* cutoffs, not estimated ones: no mixture model is fit to the
`TOTAL_CN` distribution, and no attempt is made to infer the modes. The
principled alternative is a finite mixture over integer copy states with the
decision boundary placed where the posterior odds cross 1; that is not done here
(§9).

**Assumption — and where it fails.** The baseline `CN = 2` assumes a diploid
genome. Many cancer cell lines are aneuploid or whole-genome-doubled, for which
the neutral state is `CN = 3` or `4` and an absolute threshold of 2.5
mis-classifies neutral regions as amplified. This is exactly what
`test_run_ploidy_normalised_cna.py` probes via the relative alternative
`CN / ploidy`.

**Role gating** (a conjunction restricting the flag to mechanistically coherent
cases — deletion is evidence only for a tumour suppressor, amplification only
for an oncogene):
```
has_cna_alteration =
      (is_amplification ∧ role ∈ {oncogene, both})
    ∨ (is_deletion       ∧ role ∈ {tsg, both})
    ∨ (is_amplification ∧ regime_class ∈ {activation_driven, abundance_tracking})
    ∨ (is_deletion       ∧ regime_class == loss_of_function)
```

*Provenance:* `build_cna_layer.py:30-31, 88-89`. Collapsed to one row per
(model_id, ensg_id); feeds `flags_with_driver.parquet` and the Stage 5
confidence-tier upgrades.

## 4.3 / 6.4 Variance-reduction metric

```
reduction_pct = (1 − sd(Δ_driver)/sd(Δ_any)) · 100
Δ_driver = hit_rate_driver_gated  − hit_rate_flat
Δ_any    = hit_rate_any_alteration − hit_rate_flat
```

**Model.** A ratio of two sample standard deviations computed over genes. The
statistic is `1 − s₁/s₂`; under normality the sampling distribution of `s₁²/s₂²`
is `F(n₁−1, n₂−1)`, so a confidence interval and a formal test of
`H₀: σ₁ = σ₂` are both available. Neither is computed — the number is reported
as a point estimate with no interval.

**Assumption.** `Δ_driver` and `Δ_any` are computed on the *same* genes, so they
are paired and positively correlated; the independent-samples F-test would
therefore be the wrong test. A paired analysis (e.g. a bootstrap over genes, as
already used in §C.5) is the assumption-consistent version.

*Provenance:* `04_driver_gated_routing.ipynb` / `06_held_out_eval.ipynb`;
`stage4_eval_save.py:76-79`. Quantifies how much driver-only gating reduces the
spread of hit-rate lift across `activation_driven` genes. Reported as a headline
number, not thresholded into an automatic accept/reject rule — the routing change
was accepted on this together with the held-out hit@20 results (Stage 6).

## 4.4 Alteration direction — is_lof_alteration (added 2026-08-15)

```python
is_tsg_or_both    = gene_role ∈ {tsg, both}
mut_truncating    = mut_driver ∧ (max_vep_rank ≥ 3)       # HIGH/frameshift consequence
is_lof_alteration = (has_cna_alteration ∧ is_tsg_or_both)
                  ∨ (mut_truncating     ∧ is_tsg_or_both)
```

**Model.** For a (gene, model) pair, `is_lof_alteration` is True when the alteration
evidence is consistent with loss-of-function: a copy-number deletion or a truncating
mutation (NMD-triggering consequence) in a TSG-role or "both"-role gene.

**"both" treatment note.** `gene_role == "both"` genes enter `is_lof_alteration`
as TSG-like (the CNA arm and the truncating-mutation arm both fire). However, at the
directional ranking stage (`rank_driver` in eval.py, `directional_top` in
confidence_tiers.py), "both" genes are always treated as GOF — i.e. placed in the
high-expression top-20% regardless of `is_lof_alteration`. The rationale is that the
"both" stratum has no empirical basis for directional separation at the current
sample size (n=2 test genes).

**Directionality axiom.** LOF sensitive lines are expected to have *lower* expression
of the target gene (deletion reduces transcript abundance; NMD degrades truncated
transcripts). Under the axiom:
- LOF-altered driver lines → ranked by *ascending* core_score (lower = better rank)
- GOF-altered driver lines → ranked by *descending* core_score (higher = better rank)
- Non-driver lines → ranked below all driver lines regardless of direction

This is the theoretical prediction: hit@20 improves when sensitive LOF lines are
concentrated in the low-expression tail rather than distributed uniformly.

**max_vep_rank ≥ 3 as NMD proxy.** VEP consequence ranks 0–3 (LOW/MODERATE/HIGH/frameshift).
Rank ≥ 3 captures HIGH-impact and frameshift variants, which are canonical NMD triggers
(Lindeboom et al. 2016, PMID 27618451). Missense variants (rank 2) do not reliably
reduce expression even when pathogenic, so they are excluded from the LOF branch;
their driver evidence routes through `mut_driver` but the pair is treated as GOF
for ranking.

*Provenance:* `final_pipeline/Scoring/driver_routing.py` — `is_lof_alteration`
computed before the `out_cols` merge. `max_vep_rank` output to `flags_with_driver.parquet`.

---

# Stage 5 — Confidence Tiers

`src/pipeline/05_confidence_tiers.ipynb`, `src/pipeline/build_full_predictions.py`

Deliberately **not** a numeric 0–1 score — a categorical decision rule with an
explicit precedence order:

```
tier = "high"     if regime_source == "measured" ∧ class == "abundance_tracking" ∧ n_layers == 2
tier = "moderate" (default fallback)
tier = "low"      if class == "activation_driven" ∧ rank_basis == "score_only"
tier = "unknown"  if regime_source != "measured"
# applied last (highest precedence):
tier = CNA-based upgrade                    if a CNA mechanistic hit is present
tier = "moderate" for unknown-class genes   if chronos_check == "validated"   (§C.2)
```

**Model.** A total preorder over evidence configurations, realised as a
priority-ordered rule cascade. Because later rules override earlier ones, the
mapping is a function of the full evidence tuple only if the rules are applied in
a fixed order — order is part of the specification, not an implementation
detail.

**Why categorical.** A numeric confidence would have to be a calibrated
probability to be meaningful (i.e. among pairs assigned 0.8, 80% should be
correct). No calibration set exists that is independent of the GDSC data already
used for validation, so a numeric score would be a fabricated precision. The
categorical rule makes the evidence ordering explicit and auditable instead.

**Limits.** The tiers are not calibrated: there is no measurement of
`P(correct | tier = "high")` versus `P(correct | tier = "moderate")`. Producing
that reliability curve on a held-out set is the natural next step (§9). Note
also that the `chronos_check == "validated"` upgrade inherits the uncorrected
multiplicity problem of §C.2 — under the global null roughly `α·m` genes are
promoted spuriously.

*Provenance:* `05_confidence_tiers.ipynb`; production version
`build_full_predictions.py:91-134`. Surfaced as the user-facing `confidence`
field (CLI, `explain_pair.py`, `rank_cell_lines.py`) and feeds `signal_quality`
(`np.select` over confidence/regime_source/chronos_check,
[build_full_predictions.py:189-194](../src/pipeline/build_full_predictions.py)).

## 5.1 Tumour-suppressor score inversion

```
core_score_TSG   = 1 − core_score
stratum_rank_TSG = 1 − stratum_rank
```

**Properties.** `φ(x) = 1−x` is the unique order-reversing affine involution on
[0,1] fixing the midpoint: `φ(φ(x)) = x`, `x < y ⇒ φ(x) > φ(y)`, `φ(½) = ½`.
It commutes with the percentile transform up to tie handling —
`1 − F̂ₙ(x)` is the empirical survival function, which is the empirical CDF of
`−X` — so inverting the score and inverting the raw data give the same ordering.

**Why.** For a tumour suppressor, *low* abundance is the biologically relevant
signal (loss of the suppressor), the opposite direction from oncogenes and
abundance-tracking genes. Applying the involution at build time means every
downstream consumer can sort descending unconditionally, with no per-class
branching at query time.

**Limit.** The involution reverses the *order* but preserves the *distribution
shape reflected about ½*, so the shrinkage/stratum argument of §2.6 applies
identically to the inverted scores. It does not, however, make the two
directions comparable in magnitude: a TSG at 0.9 (i.e. raw 0.1) and an oncogene
at 0.9 are equally extreme within their genes, which is all `stratum_rank`
claims.

*Provenance:* `build_full_predictions.py:65-67`.

---

# Stage 6 — Held-Out Evaluation

`src/pipeline/06_held_out_eval.ipynb`, `src/pipeline/stage4_eval_save.py`

## 6.1 Train/test gene split

```
test_genes = rng.choice(sorted(curated_genes), size=len(curated_genes)//5, replace=False, seed=42)
```

Simple random sampling without replacement, 20% held out. Sorting before
sampling makes the draw reproducible regardless of Python's set iteration order
(which is not guaranteed stable across runs or versions) — the seed alone is not
sufficient for reproducibility if the input ordering is nondeterministic.

**Limit.** A single 20% split gives one estimate with sampling variance
`≈ p(1−p)/n_test`; with a curated gene set in the low hundreds, `n_test` is a few
dozen and the standard error on any rate is substantial. Repeated splits or
k-fold cross-validation would give an interval rather than a point.

## 6.2 Vectorised routing sort key

```
sort_key = has_driver_flag.astype(float) · 1e6 + core_score
rank     = sort_key.rank(ascending=False, method="first")
```

**Model.** A lexicographic-order embedding. Ordering pairs `(a, b)` with `a`
integer and `b ∈ [0, M)` lexicographically is equivalent to ordering the scalar
`a·M + b`, provided `b` never reaches `M`. Here `core_score ∈ [0,1]` and
`M = 10⁶`, so the embedding is valid with five orders of magnitude to spare.

**Numerical limit.** IEEE-754 double precision has a 53-bit mantissa, so the
representable spacing near `10⁶` is `2⁻⁵² · 2²⁰ ≈ 2.3 × 10⁻¹⁰`. Adding `10⁶` to
a value in [0,1] therefore silently ties any two `core_score` values differing
by less than ≈ `2.3 × 10⁻¹⁰`. This is far below the ±0.044 DKW band of §0.1b, so
it is immaterial — but it is a real precision loss and the check is worth
recording.

**Tie-breaking.** `method="first"` resolves ties by row order, which is
deterministic but arbitrary; if the input frame carries any systematic ordering
(e.g. alphabetical by cell line), tie-breaking is systematically biased toward
that order rather than random. With the tie mass here this is negligible, but
`method="average"` would be the neutral choice for reporting.

## 6.3 Hit-rate@20 (hit@k)

```
hit_rate = |top_k_predicted ∩ GDSC_sensitive| / |GDSC_sensitive ∩ scored_in_core_score|
```

**Model.** Top-k recall. The denominator is the only honest one — the model can
only rank what it has a score for.

**Null model — currently missing, and this is the main gap in Stage 6.** Under
random ranking, the number of hits in the top `k` is **hypergeometric**: drawing
`k` of `N` cell lines without replacement, of which `S` are sensitive,

```
H ~ Hypergeometric(N, S, k)
E[H]   = k·S/N
Var[H] = k · (S/N) · (1 − S/N) · (N−k)/(N−1)
E[hit_rate] = E[H]/S = k/N
p_right = 1 − F_Hypergeom(h − 1;  N, S, k)
```

So the null hit rate is `k/N` — for `N ≈ 900, k = 20` that is **2.2%**,
*independent of S*. A reported hit@20 should be accompanied by the enrichment
ratio `observed / (k/N)` and the exact hypergeometric right-tail p-value; a raw
percentage cannot be read without them. The `(N−k)/(N−1)` factor is the
finite-population correction, non-negligible here since `k/N ≈ 2%`.

**Metric-theoretic limits.** hit@k is discontinuous in the score — a single rank
crossing position 20 flips a hit — and it is not a proper scoring rule, so
optimising it does not drive the score toward the true probability. Its smooth
counterpart is AUROC, which equals the normalised Mann–Whitney statistic:

```
AUC = U/(n₁n₂) = P( rank(random sensitive line) > rank(random insensitive line) )
```

AUROC is a monotone function of the rank-sum and hence closely related to the
Spearman ρ used everywhere else in this project. This relationship explains why
ρ-based validation (§C.2) and hit@20-based validation can disagree: ρ/AUROC
integrate over the whole ranking, hit@20 reads only the head of it. Where they
disagree, the head of the ranking is behaving differently from the body — which
is informative, not a contradiction.

*Provenance:* `06_held_out_eval.ipynb`; `stage4_eval_save.py:54-64`; also
`validate_chronos.hit_at_k()` and `test_run_metabolomics_mirna.hit_at_k()`.
Computed per gene, averaged by class (`abundance_tracking` vs
`activation_driven`):
```
driver_lift_pct = (mean(hit_rate_driver) − mean(hit_rate_flat)) / mean(hit_rate_flat) · 100
```
used as the primary accept/reject evidence for the Stage 4 routing change. Note
this is a ratio of two correlated means with no interval attached; a paired
bootstrap over genes (as in §C.5) would supply one.

**Rank scope (Y4, confirmed 2026-08-15):** Both `rank_flat` and `rank_driver` are
computed **within-gene** (`core_t.groupby("ensg_id")[...].rank(...)`). A rank of 1
means "most responsive within that gene's screened cell lines", not "top-ranked line
globally." Per-gene ranking is the only meaningful scope: a global rank across all
19.7M pairs would conflate genes with different numbers of screened lines and different
sensitivity distributions.

---

# Stage 7 — Explain / Rank / CLI Query Layer

`src/pipeline/explain_pair.py`, `src/pipeline/rank_cell_lines.py`, `src/pipeline/cli.py`

No new estimation — retrieval and sorting over precomputed columns
(`core_score`, `n_layers`, `stratum_rank`, `confidence`).

- **Sort logic** (`sort_gene_rows()`): `abundance_tracking` genes sort by
  `core_score` descending; other classes by `[has_driver_alteration, core_score]`
  descending. The TSG involution of §5.1 guarantees descending is always the
  correct direction, so no per-class branching is needed at query time.
- **Consistency self-check**: `rank_cell_lines.py` asserts the known BRAF/A375
  stratum-rank percentile improved after the §2.4 change (0.767 → 0.814) and
  stays within the top 5%:
  ```
  top5pct = max(1, round(total_n · 0.05))
  ```
  This is a regression-test threshold, not an estimator. Note it is a *single*
  anchor: a 0.767 → 0.814 movement is within the ±0.044 DKW band of §0.1b, so it
  is a guard against regression, not evidence of improvement.

---

# Track C — Mutations & Fusions Scoring

`src/Track - C/03_mutations_scoring.ipynb`, `src/Track - C/04_fusions_scoring.ipynb`

Track C produces **two kinds of output with different fates**:

| Track C output | Columns | Consumed by | Affects ranking? |
|---|---|---|---|
| **Cleaning** — `mutations_collapsed.parquet` | `any_driver`, `oncogene_hit`, `tsg_hit` | `mut_driver` → `has_driver_alteration` (§4.1) | **Yes — production gate** |
| **Cleaning** — `fusions_gene_level.parquet` | `max_confidence` | `fusion_driver` → `has_driver_alteration` | **Yes** |
| **Scoring** — `mutations_scores.parquet` | `p_mutation` | `explain_pair.py:137-153` (display); `has_alteration` (rejected arm) | **No** |
| **Scoring** — `fusions_scores.parquet` | `p_fusion` | `explain_pair.py:157-162` (display); `has_alteration` | **No** |

The Boolean driver columns drive every driver-gated ranking decision. The
continuous Hill/Noisy-OR scores never enter the sort order: the production gate
`has_driver_alteration = mut_driver ∨ fusion_driver` is purely Boolean and never
reads them. Their only live consumer is `explain_pair.py`, which tags them
in-code as `"CONFIDENCE MODIFIER — does not affect stratum_rank"`.

They also feed `has_alteration = (p_mutation > 0.5) ∨ (p_fusion > 0.5)` — the
any-alteration gate **rejected** at §4.1, surviving only as the comparison arm
`hr_any_alt` in the held-out evaluation.

*Consequence for calibration:* the uncalibrated parameters below matter less than
they appear, because the scores are display-only. Calibrating them is worthwhile
only if this layer is promoted into ranking. It is also why the signature-discount
test returned null (§C.7) — it tuned a component whose output does not reach the
metric being measured.

## C.1 The Hill function is a logistic regression

```
Hill(x; p₀, k) = x^k / (x^k + p₀^k)
```

**Identity.** Divide through by `x^k`:

```
Hill(x; p₀,k) = 1 / (1 + (p₀/x)^k) = 1 / (1 + e^{−k(ln x − ln p₀)}) = σ( k·ln x − k·ln p₀ )
```

**The Hill function is exactly a logistic model in `ln x`**, with slope `β₁ = k`
and intercept `β₀ = −k·ln p₀`. This converts "TEMP SCORING — p₀ and k are
placeholders pending calibration" from an open problem into a standard
estimation task:

```
ℓ(β₀,β₁) = Σᵢ [ yᵢ·ln σ(β₀ + β₁ ln xᵢ) + (1−yᵢ)·ln(1 − σ(β₀ + β₁ ln xᵢ)) ]
```

with `y` a binary ground-truth label (e.g. GDSC-sensitive). This log-likelihood
is **concave** in `(β₀, β₁)`, so the MLE is unique wherever the data are not
separable, standard errors come from the observed information matrix, and the
parameters are recovered as

```
k = β̂₁ ,      p₀ = exp(−β̂₀ / β̂₁)
```

Identifiability requires `β̂₁ ≠ 0` and both classes present. This is the concrete
calibration procedure the codebase's own comments call for.

**Properties.** Strictly increasing on `x > 0`; `Hill(p₀) = ½` (the parameter is
literally the median-effect point); asymptotes 0 and 1;

```
∂Hill/∂x = k·p₀^k·x^{k−1} / (x^k + p₀^k)²,      slope at x = p₀ is k/(4p₀)
```

`k = 1` gives Michaelis–Menten kinetics; `k → ∞` gives a step function at `p₀`;
so `k` interpolates continuously between a smooth score and a hard threshold.

**Naming caveat.** In its original setting (Hill 1910, oxygen binding to
haemoglobin) `k` is a cooperativity coefficient with a mechanistic meaning —
approximately the number of cooperative binding sites. Here it is a free
steepness parameter applied to VEP ordinals, pathogenicity scores and variant
burdens, none of which involve cooperative binding. The functional form is
borrowed; the mechanism is not. This should be stated explicitly rather than
left for a reader to infer.

## C.2 Mutations scoring

```
p_vep    = Hill(max_vep_rank,      p₀=2.5, k=2.0)    # VEP rank: 0=no_var, 1=LOW, 2=MODERATE, 3=HIGH/truncating
                                                      # rank 2 (MODERATE missense) → 0.390 < 0.5 gate
                                                      # rank 3 (HIGH/truncating)   → 0.590 > 0.5 gate
p_path   = Hill(max_pathogenicity, p₀=0.5, k=2.0)    # max(AlphaMissense, CADD)
p_burden = Hill(variant_burden,    p₀=3.0, k=1.5)    # variant_burden = variant_count.clip(upper=5)  (NOT log1p)
                                                      # 5 variants → 0.683

p_base     = 1 − (1−p_vep)(1−p_path)(1−p_burden)     # Noisy-OR
p_mutation = p_base.clip(upper=1.0)                   # DRIVER_BOOST removed 2026-08-15

p_mutation_quality = 1 − (1−p_vep)(1−p_path)         # burden excluded
p_mutation_burden  = p_burden                        # kept separate
```

**Model.** Each raw feature is mapped to an "evidence strength" in (0,1) by a
Hill/logistic squash, then combined by Noisy-OR under the assumption that VEP
impact, pathogenicity and variant burden are independent evidence of a damaging
mutation.

**Assumption violated — the same one as §2.3.** VEP consequence severity and
computed pathogenicity (AlphaMissense/REVEL) are strongly dependent by
construction: both are functions of the same substitution at the same position,
and pathogenicity predictors take consequence class as an input feature.
Noisy-OR therefore double-counts them by the non-idempotence theorem of §0.2,
inflating `p_base` by up to 0.25 where the two agree. The pipeline's own fix for
this pattern already exists (§2.4) but has not been applied here.

**`p_burden` is doubly redundant.** `variant_burden = variant_count.clip(upper=5)`
already encodes count, and the max-aggregation of §C.4 encodes count a second
time. Excluding burden from `p_mutation_quality` is a partial acknowledgement of
this.

**DRIVER_BOOST removed 2026-08-15.** The former formula `min(p_base + 0.15·driver_flag, 1.0)`
added a constant additive shift on a probability scale — not a probability operation.
An additive shift is not closed on [0,1] and creates a point mass at exactly 1.0 for
every driver hit with `p_base ≥ 0.85`. The scale-appropriate form is a log-odds shift
(a Bayes factor): `logit(p') = logit(p) + λ·driver_flag`, which is closed on (0,1)
and has the interpretation "a driver hit multiplies the odds by `e^λ`." The boost was
removed rather than replaced, since `p_mutation` is a display score only and the
Boolean `mut_driver` gate already encodes the driver call.

**max_vep_rank ≥ 3 as truncating-mutation proxy.** `max_vep_rank` is the
maximum VEP consequence severity rank across all variants for a (gene, model) pair,
with ranks on the 0–3 ordinal: 0=no_var, 1=LOW, 2=MODERATE (missense), 3=HIGH/truncating.
The threshold `max_vep_rank ≥ 3` captures HIGH/frameshift consequences — variants
expected to trigger nonsense-mediated decay (NMD) and reduce transcript abundance.
This is the basis for `is_lof_alteration` (see §4.4): a truncating mutation in a
TSG-role gene is classed as LOF, directing sensitive lines toward the low-expression
end of the ranking.

*Provenance:* max-aggregated to one row per (ensg_id, model_id) →
`mutations_scores.parquet`. Schema: ensg_id, model_id, p_mutation, max_vep_rank.
Quality is kept separate from burden because a downstream signature-discount step
scales burden but not quality:
`p_mutation_adj = NoisyOR(p_mutation_quality, m_mut · p_mutation_burden)`
(driver_boost term removed as of 2026-08-15).

## C.3 Fusions scoring

```
conf_ord = {low:1, medium:2, high:3}[max_confidence]

p_conf  = Hill(conf_ord,     p₀=2.0, k=1.5)     # low→0.261, medium→0.500, high→0.648
p_ffpm  = Hill(best_ffpm,    p₀=0.5, k=2.0)     # absolute FFPM; half-max at 0.5 (median raw ≈ 0.08)
p_recur = Hill(fusion_count, p₀=2.0, k=2.0)     # absolute count; half-max at 2; singleton→0.200

p_base   = 1 − (1−p_conf)(1−p_ffpm)(1−p_recur)
p_fusion = min(p_base + 0.15·any_in_frame, 1.0)     # INFRAME_BOOST = 0.15 (locked)
```

Same Hill → Noisy-OR → additive-boost pattern, and the same three critiques
apply (§C.1, §C.2). Two additional issues specific to fusions:

**Ordinal input to a continuous function.** `conf_ord ∈ {1,2,3}` is ordinal, but
Hill treats it as a ratio-scale quantity — it computes `1^k` vs `2^k` vs `3^k`,
which asserts that "high" is `(3/2)^k` times "medium" on the same scale. The
ordinal encoding is arbitrary; `{1,2,4}` would give different scores from the
same data. For an ordinal predictor the assumption-free treatment is one
parameter per level (i.e. three free probabilities), not a parametric curve.

**Near-constant singleton contribution.** The code uses absolute `fusion_count`
(not a percentile transform). For the 87% of fusions that appear exactly once,
`p_recur = Hill(1; p₀=2.0, k=2.0) = 1²/(1²+2²) = 0.200` — a near-constant
low-evidence contribution for the majority of rows. Recurrence provides almost
no discrimination between singleton fusions; the correct handling is to treat
recurrence as the near-binary variable it is (singleton vs recurring) rather
than a continuous score.

*Provenance:* max-aggregated to (ensg_id, model_id) → `fusions_scores.parquet` →
Stage 4. In-frame fusions are more likely to produce a translated chimeric
protein, hence the boost. **Caveat (Haas et al. 2019, PMID 31639029):** a large
fraction of RNA-seq fusion calls are read-through or trans-splicing artefacts
with no underlying DNA rearrangement; such events frequently preserve reading
frame, so the in-frame flag partially marks the dominant false-positive class.
No quality/burden split — the whole score is scaled by the signature discount
(`s_fus' = m_fus · s_fus`).

## C.4 Max-aggregation is a selection bias with a closed form

Both Track C scores are reduced to one row per (gene, cell line) by `max` over
variants/fusions.

**Model.** For `m` i.i.d. draws from `U(0,1)`, the maximum has
`P(max ≤ x) = x^m`, i.e. `max ~ Beta(m, 1)`, with

```
E[max] = m/(m+1)          Var[max] = m / ((m+1)²(m+2))
```

| variants `m` | `E[max]` under the null |
|---|---|
| 1 | 0.500 |
| 3 | 0.750 |
| 5 | 0.833 |
| 10 | **0.909** |

A gene with ten variants in a cell line receives an expected score of 0.91
versus 0.50 for a single variant, **purely from having more draws** — before any
biology. Because `p_burden` already encodes variant count explicitly, the
pipeline counts burden twice: once as a modelled feature, once as an unmodelled
extreme-value artefact.

**Measured exposure — smaller than the range suggests.** From
`cleaned_track_data/mutations_collapsed.parquet` (n = 632,919 gene × cell-line
records): `variant_count` ranges 1–37, but **91.3% of records have exactly one
variant** and only **1.65% have three or more**. For a single variant `E[max]`
is exactly the underlying value, so the bias is inert on nine records in ten.
The affected 8.7% tail is not negligible — it is concentrated in precisely the
long, frequently-mutated genes where a spurious driver call is most costly — but
the headline `E[max] = m/(m+1)` overstates the aggregate impact considerably,
and the honest statement is that this is a tail problem, not a systematic one.

The same measurement for fusions: `fusion_count` ranges 1–21 with **87.4% at
exactly 1**, which is also the source of the tie-mass problem in §C.3.

**Consequence.** This is currently harmless because `p_mutation` and `p_fusion`
are display-only (see the Track C table above). It is a blocker on promoting this
layer into ranking, since the bias correlates with gene length and mutation rate
— exactly the confounds a driver score is supposed to control for. The bias-free alternatives are an explicit maximum-of-`m` correction
(compare against the `Beta(m,1)` null rather than the raw value), or aggregation
by an idempotent operator such as the mean.

---

# Confound / Diagnostic Tests

`test_run_metabolomics_lineage_confound.py`, `test_run_metabolomics_mirna.py`,
`test_run_mirna_lineage_confound.py`, `validate_chronos.py`, `build_chronos_validation.py`

Validation experiments, not production code — but they determine which candidate
layers are accepted or rejected, and are the source of the narrative in
[DECISIONS.md](../DECISIONS.md).

## C.5 Spearman ρ — the primary validation statistic

Definition, tie handling, exact null, the validity range of the t-approximation,
and the Fisher variance-stabilising transform are all given in §0.4 and are not
repeated here.

*Provenance:* vectorised manual implementation in
`build_chronos_validation.py:65-87` (avoids a per-gene Python loop over thousands
of genes); `scipy.stats.spearmanr` elsewhere. Measures whether `core_score` (or a
candidate variant) agrees in ranking with an independent ground truth — Chronos
CRISPR essentiality or GDSC drug sensitivity. This is the load-bearing metric for
every "does adding layer X help" decision in the project.

**What ρ does and does not establish.** ρ is a measure of monotone association,
not of calibration or of accuracy at the head of the ranking. Two scores with
identical ρ can have very different hit@20 (§6.3). Reporting both, as this
project does, is the right call; where they disagree the head/body distinction of
§6.3 is the explanation.

## C.6 Validated / Inverted / None classification — multiplicity, measured

```
RHO_THRESHOLD = 0.1 ;  P_THRESHOLD = 0.05 ;  MIN_N = 30
"validated"  if ρ >  0.1 ∧ p < 0.05
"inverted"   if ρ < −0.1 ∧ p < 0.05
"none"       otherwise
```

**Model.** A two-sided per-gene test of `H₀: ρ = 0` combined with a minimum
effect size. `MIN_N = 30` is justified as the validity floor for the Spearman
t-approximation (§0.4), not as a power calculation.

**No multiple-testing correction is implemented anywhere in the repository**
(verified by search over `src/` for `multipletests`, `benjamini`, `fdr_bh`,
`bonferroni` — no matches), and the rule is applied independently to every gene.
That is a genuine omission. However, running the analysis on the actual output
shows the rule is **accidentally FDR-safe**, and the reason is worth stating
precisely rather than relying on.

**Measured** (`src/pipeline/outputs/chronos_validation.parquet`, m = 16,866 genes
tested, median n = 823 cell lines per gene, minimum n = 247):

| Quantity | Value |
|---|---|
| Genes with `p < 0.05` | 4,362 |
| Expected under the global null (`0.05·m`) | 843 |
| Storey `π̂₀` (λ = 0.5) | **0.630** |
| Classified `validated` | 769 |
| Classified `inverted` | 1,258 |
| **Largest p-value among `validated`** | **4.03 × 10⁻³** |
| Benjamini–Hochberg threshold at q = 0.05 | 6.92 × 10⁻³ |
| `validated` genes surviving BH at q = 0.05 | **769 of 769 (all)** |

Two conclusions.

**(1) There is real signal.** 4,362 genes clear `p < 0.05` against 843 expected
under the global null — a 5.2× excess — and Storey's estimator puts the null
fraction at `π̂₀ ≈ 0.63`, i.e. roughly 37% of tested genes carry a genuine
association with Chronos essentiality. This is not a set dominated by noise.

**(2) The effect-size floor is the binding constraint, not the p-value floor.**
The two thresholds cross where `ρ = 0.1` attains `p = 0.05`, i.e. where
`0.1·√(n−1) = 1.96`, so at

```
n* = 385
```

Above `n*` the `ρ > 0.1` floor is stricter than `p < 0.05`; below it, the
p-value floor binds. Since the median gene here has `n = 823` and the minimum is
`n = 247`, most genes sit in the ρ-binding regime: at `n = 823`, `ρ = 0.1`
corresponds to `p ≈ 0.0041`, roughly **twelve times stricter** than the nominal
0.05. Consequently 2,335 genes clear `p < 0.05` but are excluded by the effect
floor, and every surviving `validated` gene has `p ≤ 4.03 × 10⁻³` — comfortably
below the BH threshold.

**Correction to an earlier version of this document.** A previous draft asserted
that "the binding constraint is `p < 0.05` alone" and projected ≈75 spurious
validated calls at `m ≈ 3000`. Both are wrong: `m` is 16,866, and the ρ-floor
binds for the large majority of genes. The multiplicity exposure is materially
smaller than stated, and the `validated` set survives FDR control intact.

**What should still change.** The rule is safe by arithmetic coincidence, not by
design — it depends on `n` being large, and it is not safe for the low-`n` tail
(`n < 385`, where `ρ = 0.1` does not imply `p < 0.05`). Implementing BH at a
stated `q` and reporting `π̂₀` would make the guarantee explicit and would extend
it to genes the current rule cannot protect.

*Downstream note:* `validated` genes are promoted from
`confidence = "low"/"unknown"` to `"moderate"` in Stage 5, so this rule has
user-facing consequences.

*Provenance:* `build_chronos_validation.py:41-43, 89-92`; identical thresholds in
`test_run_metabolomics_mirna.classify():43-45, 100-105`. Written to
`chronos_validation.parquet` (`ensg_id, n, rho, pval, chronos_check`), consumed by
`build_full_predictions.py` (lines 82-89, 127-133, 158-168).

## C.8 Bootstrap confidence interval on ρ

```
for b in 1..1000:
    sample = resample_with_replacement(cell_lines)
    ρ_b    = spearman(x[sample], y[sample])
CI_95 = [percentile(ρ_b, 2.5), percentile(ρ_b, 97.5)]
```

Method, pivotality requirement, the cell-line exchangeability assumption (and its
violation by lineage clustering), Monte Carlo error at `B = 1000`, and the BCa /
Fisher-z alternatives are given in §0.5.

*Provenance:* `validate_chronos.bootstrap_rho():101-110`. Quantifies uncertainty
around a point ρ (3-layer vs 2-layer score against GDSC for BRAF/BCL2/PTEN) — a
single ρ without an interval cannot distinguish a real effect from sampling noise
at this cell-line count.

## C.9 Paired bootstrap on Δρ

```
for b in 1..2000:
    sample     = resample_with_replacement(rows)     # same rows for both scores
    ρ_baseline = spearman(baseline[sample],  ground_truth[sample])
    ρ_candidate= spearman(candidate[sample], ground_truth[sample])
    Δ_b        = ρ_candidate − ρ_baseline
CI_95 = [percentile(Δ, 2.5), percentile(Δ, 97.5)]
significant  if 0 ∉ CI_95
```

**Why pairing gains power.** From §0.5,
`Var(Δ) = Var(ρ_c) + Var(ρ_b) − 2Cov(ρ_c, ρ_b)`. Because both scores are
computed on the same resampled rows and share the baseline component, the
covariance is strongly positive and `Var(Δ)` is much smaller than the unpaired
sum. This also means the common "the two confidence intervals overlap, therefore
no difference" reading is *invalid* here — overlapping marginal intervals are
entirely compatible with a Δ interval excluding zero.

**Result in both tests (metabolomics, miRNA): CI included 0 — not significant,
and the point estimate moved in the wrong direction** (lower correlation). Both
rejected from production — see [DECISIONS.md](../DECISIONS.md).

**Interpretive limit.** A CI containing 0 is failure to reject, not evidence of
equivalence. To claim the candidate layer adds nothing, the correct statement is
an equivalence test: the whole CI must lie inside a pre-specified
practically-negligible band `[−δ, +δ]`. Reporting the realised CI width against a
stated δ would convert "not significant" into "demonstrably negligible."

*Provenance:* `test_run_metabolomics_lineage_confound.py:100-114`;
`test_run_mirna_lineage_confound.py:95-113`.
`candidate = 0.5·baseline_core_score + 0.5·global_pct` (see §C.10).

## C.11 Z-score standardisation before averaging

```
z        = (x − mean(x)) / std(x)              # per metabolite or per miRNA
global_z = mean(z across all features, axis=1) # per cell line
```

**Why standardise first.** Averaging unstandardised features weights each by its
raw standard deviation, so a handful of high-variance features would dominate the
mean. Standardising equalises the contribution — `Var(mean of k standardised
features) = (1 + (k−1)ρ̄)/k`, which for large `k` tends to `ρ̄`, the average
pairwise correlation.

**That limit is the diagnosis.** If the features are close to uncorrelated
(`ρ̄ ≈ 0`), the average of hundreds of them converges to a constant plus noise of
variance `≈ 1/k` — i.e. **the global scalar carries almost no information by
construction**, regardless of biology. This is a mathematical prediction of the
null result observed in §C.9, and it is a stronger explanation than "averaging
across unrelated pathways": the operation destroys signal whenever the averaged
features are not coherently correlated.

*Provenance:* `test_run_metabolomics_lineage_confound.py:39-40`;
`test_run_mirna_lineage_confound.py:35-38`.

## C.12 One-way ANOVA and effect size — **η² is biased, use ω²**

```
F, p = f_oneway(group₁, …, group_k)            # groups = per-lineage arrays of global_z
η² = SS_between / SS_total
SS_between = Σ_lineage n_lineage·(mean_lineage − grand_mean)²
SS_total   = Σᵢ (xᵢ − grand_mean)²
```

**Model.** Decomposition of total sum of squares into between- and within-group
components; under H₀ (equal group means) and normality,
`F = MS_between/MS_within ~ F(k−1, N−k)`.

**η² is an upward-biased estimator of the population effect.** Under H₀,
`η² ~ Beta((k−1)/2, (N−k)/2)`, so

```
E[η² | H₀] = (k−1)/(N−1)
```

which is **not zero** — it grows with the number of groups. The bias-corrected
estimator is

```
ω² = (F − 1)(k − 1) / ( (F − 1)(k − 1) + N )
```

Applied to this project's own results:

| | k | N | F | p | reported η² | **E[η² \| H₀]** | **ω² (corrected)** |
|---|---|---|---|---|---|---|---|
| Metabolomics | 24 | 928 | 1.713 | 0.0198 | 0.042 | **0.0248** | **0.017** |
| miRNA | 25 | 948 | 2.975 | 2.5e-06 | 0.072 | **0.0253** | **0.048** |

**This materially strengthens the project's own conclusion.** The metabolomics
η² of 4.2% is only 1.7× what pure noise produces with 24 groups at this sample
size; after correction the lineage effect is **1.7%**, not 4.2%. The miRNA effect
survives correction better (4.8%, ≈ 2.8× the null floor) but is still small.
Lineage does not dominate either global scalar.

**Assumptions.** `f_oneway` assumes independent observations, approximately
normal residuals, and homogeneous variances. With group sizes as unequal as 24–25
lineages over ~940 cell lines, variance heterogeneity is likely; Welch's ANOVA
would be the assumption-robust alternative, and Kruskal–Wallis the distribution-free
one. Neither is run.

**Reading the result correctly.** Both tests are statistically significant and
practically small — the classic large-`N` situation where `p` measures sample size
more than effect. Combined with the null/negative Δρ of §C.9 and the averaging
argument of §C.11, the conclusion is that the global scalar is closer to noise
than to a usable signal of any kind, lineage or gene. Not folded into `core_score`
either way.

*Provenance:* `test_run_metabolomics_lineage_confound.py:55-63`;
`test_run_mirna_lineage_confound.py:52-59`. Results recorded in
[DECISIONS.md](../DECISIONS.md).

# `new-arch` — Continuous vs Discrete Diagnostic

`src/new-arch/01_long_form_gene_axis_diagnostic.ipynb`

```
score_continuous = (expression − min)/(max − min)      # per gene, min-max
score_discrete   = (expression > 1.0)                  # HPA_THRESHOLD = 1.0 log2-TPM
disagreement     = |score_continuous − score_discrete|
```

**What the comparison actually measures.** `score_discrete` is `𝟙[x > τ]`, the
indicator of a threshold event; `score_continuous` is a min-max normalisation.
Their difference is *not* a measure of the continuous-vs-discrete paradigm gap
alone, because min-max and the HPA threshold are on different scales: min-max
places the gene's own observed minimum at 0 regardless of whether that minimum is
above or below `τ = 1.0` log₂-TPM. Genes entirely expressed above threshold get
continuous scores spanning [0,1] against a constant discrete score of 1, so
`disagreement` is large by construction rather than by disagreement. A
scale-matched comparison would use the percentile transform (§2.1) against the
empirical proportion `F̂ₙ(τ)`.

**Where the two paradigms genuinely differ.** The indicator discards all
within-class ordering (an information-theoretic loss of `H(X) − H(𝟙[X>τ])`
bits), while the continuous score discards none. The threshold gains
interpretability — "expressed" is a claim a biologist can act on — at the cost of
discarding everything about how much.

Explicitly diagnostic — "no fix is proposed here, that is a separate, later
decision." Not part of the production formula.

---

# Cross-Pipeline Formula Reuse

| Object | Canonical source | Reused in | Section |
|---|---|---|---|
| Empirical CDF / PIT `rank(pct=True)` | `scoring_variants.normalise()` | core_score (E, P), fusions (ffpm/recur), Chronos pct, metabolomics flux proxy, global scalars, stratum_rank | §0.1 |
| Noisy-OR (probabilistic-sum t-conorm) | `scoring_variants.combine()` | core_score (superseded), mutations `p_mutation`, fusions `p_fusion` | §0.2, §2.3 |
| Idempotent weighted mean / GLS weights | `scoring_variants.combine_three()` | production core_score (2-layer), 3-layer Chronos blend (rejected) | §0.3, §2.4 |
| Hill = logistic in `ln x` | Track C `hill()` | mutations (vep/path/burden), fusions (conf/ffpm/recur), `hill` normalisation variant | §C.1 |
| Spearman ρ + t-approximation | `scipy.stats.spearmanr` / manual rank impl | Chronos validation, GDSC validation, all confound tests | §0.4 |
| hit@k (hypergeometric null) | `validate_chronos.hit_at_k()` | Chronos validation, held-out eval (Stage 4/6), confound tests | §6.3 |
| Bootstrap CI (percentile method) | `validate_chronos.bootstrap_rho()` | 3-layer ρ CI, paired Δρ CI in confound tests | §0.5 |
| Order-reversing involution `1−x` | `build_full_predictions.py` | TSG inversion, Chronos percentile inversion | §5.1 |
| Extreme-value aggregation `max` | Track C notebooks | mutations, fusions collapse to (gene, cell line) | §C.4 |

---

# What's Calibrated vs Placeholder

**Measured or fitted**
- `ρ_EP = 0.46` — literature (Nusinow et al. 2020) plus internal validation.
  Caveat: cited as a native-scale correlation, applied to rank-transformed
  layers (§0.3).
- `ρ_EC = ρ_PC = 0.005` — measured directly in this project (§2.7).
- Everything downstream of these in §2.4 — though note the two-layer production
  weights are `(½,½)` and therefore do not depend on them at all.

**Conventional, not estimated** *(defensible, but they are conventions and should
be named as such)*
- `AMP_THRESHOLD = 2.5`, `DEL_THRESHOLD = 1.5` — a ±0.5 band around diploid; no
  mixture model fit (§4.2).
- `RHO_THRESHOLD = 0.1`, `P_THRESHOLD = 0.05` — chosen, applied consistently, and
  uncorrected for multiplicity (§C.6).
- `MIN_N = 30` — this one *is* principled: the validity floor of the Spearman
  t-approximation (§0.4).
- `1e6` lexicographic multiplier — valid embedding, verified against float64
  precision (§6.2).

**Explicit placeholders, not calibrated**
- Every Hill `(p₀, k)` pair in Track C mutations/fusions, and the `hill`
  normalisation variant (`k = 2, p₀ = 0.5`). All flagged in-source. §C.1 gives
  the concrete MLE procedure that would calibrate them.
- `DRIVER_BOOST = 0.15` — **REMOVED 2026-08-15** (was additive on a probability scale;
  §C.2 gap 20 closed). `INFRAME_BOOST = 0.15` in fusions scoring remains (display-only).
- `ALPHA = 1.0` in the miRNA dampener — the value at which the dampener reaches
  exactly 0 (§C.14).

**Measured against held-out data, not in-sample** *(added after the §C.15–C.18
validation round)*
- The gate/ranker split — gate pAUC 0.593–0.601 vs ranker pAUC 0.500–0.518 —
  replicated on two independent CRISPR screens, `p = 0.19` for the difference
  between them (§C.15).
- The percentile's incremental depletion over absolute detection, 0.810–0.827,
  spread 0.017 across six calibration rules and 6/6 clearing a 10% effect floor
  (§C.16).
- The detection threshold in §C.16 is calibrated as a *false-positive rate on a
  curated negative control set* (438 OR\*/TAS2R\*/VN1R\* genes, 7.690 log₂ units
  below the transcriptome median), not chosen as a free parameter.

**Tested and explicitly rejected**
Chronos as a third core_score layer (§2.7), metabolomics global scalar, miRNA
global scalar, metabolomics enzyme-flux proxy, miRNA-inhibitor dampening —
all documented with reasoning in [DECISIONS.md](../DECISIONS.md), and each with a
statistical account above of *why* the null result was structurally likely.

Added by the §C.15–C.18 round:
- **Absolute detection thresholds as a replacement for the percentile** (§C.16).
  Apparent +0.030 pAUC advantage was selection optimism (measured optimism
  `+0.0302`); held-out delta `−0.0050`, wins 10/20 reps. Also inert for the
  majority of genes (excludes nothing for 393–591 of 778).
- **Lineage-conditioned percentiles** — Δ = 0.0, `p = 0.54` (§C.15).
- **Conformal exclusion** — valid but `ratio_vs_random` 0.93–0.99 (§C.15).
- **Stratum centring** — `centring_gain` negative at `p = 1.5e-17` / `2.4e-67`;
  RNA alone beats every Stouffer variant (§C.18).
- **HGNC gene groups as a paralog proxy** — failed its own validity check,
  `p = 0.7766` in the wrong direction (§C.17).

---

# §9 — Open mathematical gaps

Recorded honestly, in rough order of how much they affect conclusions.

1. **No null model for hit@k** (§6.3). Raw hit rates are reported without the
   hypergeometric baseline `k/N ≈ 2.2%` or an enrichment ratio, so the headline
   metric of the whole project is quoted without the number needed to read it.
   Fix: report observed/expected and an exact hypergeometric p-value.
2. **Noisy-OR is still used in Track C** (§C.2, §C.3) on features that are
   demonstrably dependent (VEP consequence and computed pathogenicity), despite
   the same failure having been diagnosed and fixed in Stage 2. Currently masked
   by those scores being display-only.
3. **No multiple-testing correction implemented** (§0.6, §C.6) — though measured
   analysis shows the `validated` set survives BH at q = 0.05 intact, because the
   `ρ > 0.1` effect floor is ~12× stricter than the nominal `p < 0.05` at typical
   `n`. The guarantee holds by arithmetic coincidence rather than by design, and
   does not extend to genes with `n < 385`. Fix: implement BH explicitly and
   report `π̂₀` (currently ≈ 0.63).
4. **Max-aggregation extreme-value bias** (§C.4) — `E[max] = m/(m+1)` makes
   multi-variant genes score higher mechanically. Measured exposure is a 8.7%
   tail (91.3% of records carry a single variant), so this is a tail problem
   rather than a systematic one, but it still blocks promotion of Track C scores
   into ranking.
5. **Per-layer variances never estimated** (§2.4). The GLS argument assumes equal
   variance across layers; if proteomics is noisier than expression, `(½,½)` is
   not optimal and the correct weights follow from `Σ⁻¹𝟙` with unequal diagonal.
6. **Confidence tiers are uncalibrated** (§5) — no measurement of
   `P(correct | tier)`. A reliability curve on held-out data would turn the
   ordering into a claim.
7. **Bootstrap resampling ignores lineage clustering** (§0.5), while §C.12
   simultaneously demonstrates lineage structure exists. Clustered resampling
   would make the two consistent.
8. **CNA thresholds assume diploidy** (§4.2) in a panel that is substantially
   aneuploid. `test_run_ploidy_normalised_cna.py` probes this; the result is not
   yet folded in.
9. **`product_floor` is unanalysed** (§2.5) — its position within the
   Fréchet–Hoeffding interval, and hence what dependence structure it implicitly
   assumes, has never been characterised.
10. **`ρ_EP` scale mismatch** (§0.3) — a native-scale literature correlation used
    to weight rank-transformed layers. Immaterial in production (weights are
    `(½,½)` regardless) but it would matter for any three-layer revival.
11. **Cross-stratum calibration deferred** (§2.6) — within-stratum ranks are
    comparable in distribution but not in evidence weight (`n_eff` 1.00 vs 1.37).
12. **Single-anchor regression tests** (§7) — the BRAF/A375 0.767 → 0.814 check
    is within the ±0.044 DKW band, so it guards against regression rather than
    demonstrating improvement.

*Added by the §C.15–C.18 validation round, in the same order of consequence:*

13. **The floor effect is not separable from a profiling confound** (§C.15). At
    the bottom of the score, dependency depletion among profiled vs unprofiled
    lines differs by ratio 1.07 at `Fisher p = 0.246` — indistinguishable. The
    gate result is real in the mid-range (ratio 1.62, `p = 0`) but the extreme
    floor cannot be attributed to biology rather than missingness. This is the
    most serious open gap in the current evidence base and must be quoted
    alongside the gate result. Fix: a measurement model that separates "not
    expressed" from "not measured", which the coverage matrix can support.
14. **`n_layers == 2 → "high"` is contradicted by measurement** (§C.18, §5).
    RNA alone outperforms every two-layer Stouffer variant, so the tier rule at
    [build_full_predictions.py:151](../src/pipeline/build_full_predictions.py)
    and [05_confidence_tiers.ipynb:117](../src/pipeline/05_confidence_tiers.ipynb)
    promotes rows on a basis the project's own test rejects. Not yet changed in
    code. Fix: demote the rule, or re-derive the layer weights with unequal
    variances (gap 5 above is the same gap seen from the other side).
15. **The sort key is unsupported by any measurement** (§C.15, §7).
    `explain_pair.sort_gene_rows` uses `core_score` as primary key for
    `abundance_tracking` and as a tiebreak elsewhere, but ranker-region pAUC is
    0.500–0.518. There is no evidence for any *ordering* within the retained set.
    Fix: present the retained set as a set, or find a ranking signal that
    replicates.
16. **Paralog buffering is unaddressed, not cleared** (§C.17). The only available
    annotation failed validation in the wrong direction (`p = 0.7766`). The
    finding survives restriction to HGNC-singleton genes, but that is not the same
    as controlling for paralogy. Fix: obtain Ensembl Compara paralog pairs with
    % sequence identity, or the DepMap / Dede / Ito paralog sets.
17. **Half-sample splits are unstratified** (§C.16). Lineage composition varies
    freely between halves, so the percentile's lower held-out variance (range
    0.078 vs 0.170) is a real observation but not lineage-controlled. This is the
    same clustering issue as gap 7, in the cross-validation rather than the
    bootstrap. Fix: lineage-stratified splitting.
18. **Label reliability bounds every ranking claim** (§C.15). Broad and Sanger
    agree with each other at median `ρ = 0.172`; only 1.9% of genes exceed 0.5.
    No AUROC against these labels can be interpreted without that ceiling, and
    the project currently quotes AUROC figures without it. Fix: report every
    ranking metric against the label-label agreement ceiling for the same gene
    set.
19. **CNA residual exceeds the label noise floor** (§C.18). Deleted-line
    dependency runs 3.04% above the Chronos FPR of 0.93%, flagged
    `residual_above_noise__check_calls`, and recurrence stratification shows the
    residual concentrated in 1–2-line ("likely mis-call") genes at 4.70%. Either
    the deletion calls carry error or deletion does not fully abolish
    dependency; the test cannot currently distinguish these.

*Gaps closed 2026-08-15 (six mathematical corrections):*

20. **DRIVER_BOOST (+0.15) — CLOSED.** The additive constant on a probability
    scale (§C.2) has been removed. `p_mutation = p_base.clip(upper=1.0)`. The
    display score is no longer inflated by an uncalibrated constant.
21. **Silent-lineage guard — CLOSED.** The notebook's third guard
    (`SILENT_FRAC=0.20`, `EXPRESSED_MIN=1.0 log2 TPM+1`) has been restored to
    `common.py:score_source_lineage()`. Genes expressed in < 20% of a lineage's
    lines now produce NaN z-scores (not floor z-scores) for that lineage, which
    propagate to NaN core_score and are excluded from ranking.
22. **MNAR shrink in protein_scorer — CLOSED.** The multiplicative shrink
    `z_t = z_raw × √(n_obs/N)` has been removed. Protein z-scores are now passed
    through as `float32` without an uncalibrated multiplicative scaling step. The
    Stouffer shrink on the RNA side (`stouffer_combine` in `common.py`) is separate
    and was not changed.
23. **Orphan pairs scored as 0.0 — CLOSED.** (gene, model) pairs with alteration
    evidence but no core_score now receive `core_score=NaN` and are excluded from
    hit@20 ranking by the `core_score.notna()` filter. Previously they were
    effectively scored as 0.0 (the worst rank), which artificially inflated the
    apparent gate benefit for those genes.


---

## Retired mathematics

Eleven sections have been moved to
[MATH_TRIAL_AND_ERROR.md](MATH_TRIAL_AND_ERROR.md): formulas that were
superseded, tested and rejected, or built on a claim that has since been
retracted. They are kept in full, with the reason for retirement stated on
each, because the failures are part of the record — most of this project's
findings are negative ones and the derivations that produced them are worth
keeping.

Most consequential: **§C.15, the partial-AUROC gate/ranker mathematics.** The
machinery is sound and still used in `gate_audit/`; the claim it supported is
withdrawn.

---

# §10 — Mathematics of the current components

Added after `gate_audit/`, `cell_similarity/` and the Kleene layer. These are the
formulas the live pipeline actually computes.

## 10.1 Kleene strong three-valued logic

Truth set `{T, F, U}` with

```
AND  | T  F  U          NOT
-----|---------         ---------
 T   | T  F  U          T -> F
 F   | F  F  F          F -> T
 U   | U  F  U          U -> U
```

`AND` is the minimum under the order `F < U < T`; `NOT` is the order-reversing
involution fixing `U`. The two properties that matter operationally:

- `T ∧ U = U` — an unresolved clause cannot produce a confident conjunction, so
  an UNKNOWN line can never enter Matches.
- `F ∧ U = F` — one definitive exclusion is sufficient regardless of what else is
  unknown, so Excluded is not weakened by missing data elsewhere.

Implementation `src/pipeline/multi_gene_kleene.py`. Spec
[MULTIGENE_SPEC.md](MULTIGENE_SPEC.md).

## 10.2 Measurement band on a detection threshold

For a threshold `τ` on a measured quantity `x` with single-measurement standard
deviation `σ`, a call is unresolvable when `|x − τ| ≤ 1.96 σ`.

`σ` is estimated from technical replicates. With `n` models carrying two
independent profiles, the paired difference `d = x₁ − x₂` has `Var(d) = 2σ²`, so

```
σ = SD(d) / √2
```

Measured: `SD(d) = 0.5042` log2 units over 64,000 paired gene observations from
16 models, giving `σ = 0.3565` and a band of `±0.699`. Both profiles of a model
share a library-prep batch, so this is a **lower bound** on σ — the band errs
toward abstention, which is the safe direction.

## 10.3 Similarity between cell lines

For a line × feature matrix, features are standardised to z and similarity is

```
S(a, b) = ⟨ẑₐ, ẑ_b⟩ / (‖ẑₐ‖ ‖ẑ_b‖)
```

i.e. Pearson between lines, equivalently cosine on centred data. The
top-variance feature restriction is a signal-to-noise choice: invariant features
add a constant to every pair and compress the dynamic range of the whole matrix.

## 10.4 Cluster bootstrap over the correct unit

For a statistic computed on **pairs** drawn from `n` items, the resampling unit
is the **item**, not the pair. Draw `i₁…iₙ` with replacement from the items,
rebuild the pair set, and **discard pairs formed between two draws of the same
item** — those have similarity and outcome correlation 1 by construction and
bias every replicate upward.

Motivation: with 1,673 lines each line appears in ~1,672 pairs, so a pair-level
interval treats one line's idiosyncrasies as ~1,672 independent observations.
Measured design effect on the gate: **2.4×** the nominal standard error
(`gate_audit/04`).

Kish effective sample size under exchangeable within-cluster correlation:

```
n_eff = (Σ mⱼ)² / Σ mⱼ²
```

for cluster sizes `mⱼ`. Measured 8.7 effective clusters from 45 lineages.

## 10.5 Attenuation and disattenuation

An observed correlation between two noisy measurements is attenuated by their
reliabilities. For a difference of correlations `Δ` measured against a readout of
reliability `r_xx`,

```
Δ_true ≈ Δ_observed / r_xx
```

Reliability is estimated by test–retest between two independent measurements of
the same quantity, then projected to a `k`-measurement composite by
**Spearman–Brown**:

```
r_k = k·r₁ / (1 + (k−1)·r₁)
```

Measured: CRISPR `r₁ = 0.371`, RNAi `r₁ = 0.144` (`r₂ = 0.252`), GDSC2 drug
`r₁ = 0.547` (`r₂ = 0.707`).

**Both arms must be corrected or neither.** Disattenuating one side of a
comparison inverted a conclusion here — the CRISPR/RNAi ratio moved from 0.99×
to 0.37× once CRISPR's own reliability was included. Report observed values as
primary; the correction is secondary and is only needed where reliability is low.

## 10.6 Split-half reliability — the non-circular stability test

For a graph built from features, randomly halve the feature set, build twice, and
compare the top-`k` neighbour lists:

```
reliability = mean over items of |N_A(i) ∩ N_B(i)| / k
```

This measures whether the graph is a property of the **data** rather than of the
metric. Metric-agreement checks (Pearson vs Spearman) are **circular** for a rank
transform — ranking twice is near-idempotent, so a rank-transformed axis scores
84.2% on metric agreement and gains nothing on split-half. Measured: RNA 83.7%,
miRNA-log1p 58.8%, metabolomics 42.9%.

## 10.7 Participation ratio — effective dimensionality

For a variance spectrum `v₁…v_p`,

```
PR = (Σ vⱼ)² / Σ vⱼ²
```

the number of features effectively carrying the geometry. Measured: RNA 1,648 of
2,000; metabolomics 115 of 225; **miRNA 13 of 734** — the diagnosis for miRNA's
instability, and the reason a transform fixed it while no transform fixes
metabolomics (which is underdetermined rather than skewed).

---

# §11 — Final Pipeline: Mathematics and Audit (COWORK v4.0)

This section records the mathematics of the production `final_pipeline/` codebase
as verified by the COWORK v4.0 code-vs-diagram reconciliation sweep (X1–X8,
2026-08-14). All file and line references are confirmed against the live source.
Items flagged **[ESCALATION]** require a decision by Fiona before the pipeline
result is cited in the dissertation.

---

## 11.1 Robust z-score (both arms)

**Implemented form — median/MAD, not mean/SD:**

```
robust_z(x) = (x − median(x, axis=0)) / (MAD(x, axis=0) · 1.4826)
```

where `MAD = median(|x − median(x)|)` and `1.4826 = 1/Φ⁻¹(0.75)` is the
consistency factor that makes MAD an unbiased estimator of σ for Gaussian data.

**Why 1.4826.** Under X ~ N(μ, σ²), the median absolute deviation satisfies
`E[MAD] = σ · Φ⁻¹(0.75) ≈ 0.6745σ`, so dividing by 0.6745 (≡ multiplying by
1.4826) gives E[robust_z] = (x − μ)/σ in expectation.

**Distinction from §C.11.** §C.11 uses classical z-score `(x − mean)/std`
applied to metabolomics/miRNA global scalars in the *old* pipeline. The
final_pipeline uses the robust form throughout; `std(x)` in §C.11 is correct for
that context. Do not conflate the two.

*Provenance:* `final_pipeline/utils/common.py:robust_z_matrix`, lines 107–112.
Used by both `rna_scorer.py` and `protein_scorer.py` before source-level
aggregation.

---

## 11.2 Stouffer combination and arm weights

The core z-score is a Stouffer-weighted combination of per-arm z-scores:

```
core_z = (W_RNA · rna_pct_z + W_PROT · prot_resid_z) / sqrt(W_RNA² + W_PROT²)
core_score = Φ(core_z)                   # standard normal CDF
```

Constants (confirmed at `final_pipeline/Scoring/core_score.py:39–40`):

| Constant | Value | Implied quantity | Source |
|---|---|---|---|
| `W_RNA` | `sqrt(1.431) = 1.1962` | n_eff = 1.431 from Kish formula, ρ̄=0.548 across 3 RNA sources | measured X2, corrected 2026-08-15 |
| `W_PROT` | `sqrt(1.45) = 1.2042` | n_eff = 1.45 from ProCAN-CCLE platform ρ=0.373 (NOT RNA-protein ρ=0.353) | comment |

Note: `RHO_PRIOR = 0.353` (median per-gene RNA-protein Pearson r) is a **separate** quantity used only
as the EB prior for per-gene β_g computation. It has no connection to W_PROT.

**Stouffer property.** The denominator `sqrt(W_RNA² + W_PROT²)` normalises the
combined z to unit variance under independence, so `core_z ~ N(0,1)` if both arm
z-scores are standard normal. Variance shares (corrected W_RNA):

```
share_RNA  = W_RNA²  / (W_RNA² + W_PROT²) = 1.431 / 2.881 = 49.7%
share_PROT = W_PROT² / (W_RNA² + W_PROT²) = 1.450 / 2.881 = 50.3%
```

The two arms now contribute nearly equally — RNA is the marginally weaker of the two
(1.196 < 1.204), the opposite of the provisional configuration (RNA 58% vs protein 42%).

---

## 11.3 W_RNA calibration — **[RESOLVED 2026-08-15]**

X2 measured the pairwise Spearman ρ between all three RNA sources (DepMap, HPA,
GEO) using 2,000 sampled genes (seed=42) and lineage-conditioned z-scores:

| Pair | Median per-gene Spearman ρ |
|---|---|
| DepMap × HPA | **0.793** |
| DepMap × GEO | 0.417 |
| HPA × GEO | 0.435 |
| **Average ρ_rna** | **0.548** |

The current `W_RNA = sqrt(2.00)` implies `n_eff = 2`, which under the Kish
formula `n_eff = N/(1+(N−1)·ρ̄)` requires `ρ̄_implied = 0.25`. Measured ρ̄ is more
than twice that.

**Corrected values:**

```
n_eff_corrected = 3 / (1 + 2 × 0.548) = 3 / 2.096 = 1.431
W_RNA_corrected = sqrt(1.431) = 1.196          (−15.4% from current 1.414)
```

Corrected variance shares: RNA 49.7%, Protein 50.3% — the arms nearly balance.
With corrected W_RNA, `W_RNA_corrected (1.196) < W_PROT (1.204)` — the RNA arm
would be *marginally* the weaker of the two, the opposite of the current
configuration (RNA 58% vs protein 42%).

**Stage 6 A/B test (2026-08-15):** Option (A) was selected — `W_RNA` corrected and
pipeline re-run (stages 4→5→6). The evaluation script used `SENSITIVITY_MAD_Z = -1.0`
throughout.

| Class | flat (old W_RNA=√2) | flat (new W_RNA=√1.431) | driver (unchanged) |
|---|---|---|---|
| oncogene (11 genes) | 0.0269 | **0.0296** (+10%) | 0.0538 |
| tsg (1 gene) ⚠ NOT INTERPRETABLE | 0.0081 | **0.0163** (+101%) | 0.0407 |
| both (1 gene) ⚠ NOT INTERPRETABLE | 0.0100 | 0.0100 | 0.0000 |

*tsg and both each contain n=1 gene — the per-stratum hit@20 values are a single gene's
ranking, not a stratum estimate. The percentage changes (+101%, 0%) are not interpretable
as stratum effects. The pooled oncogene result (n=11) is the load-bearing number.*

**Decision: correction applied and kept.** Flat hit@20 improved for both oncogene
and tsg with the corrected W_RNA. The driver hit@20 (which depends on the gene-role
CNA/mutation gates, not the core_score weights) is unchanged. The decrease in driver−flat
delta is because the baseline improved, not because the driver signal weakened.

The corrected W_RNA (√1.431 = 1.196) is now live in `core_score.py:39`. All outputs
downstream (core_score.parquet, evidence_ledger.parquet) have been regenerated.

*Source: `diagnostics/X2_rna_rho.py` (measurement), `core_score.py:39` (fix applied).*

---

## 11.4 Platform tier (protein arm)

The protein arm applies a cross-platform correlation tier before combining sources:

```
tier(ρ_platform) = "consistent"  if ρ ≥ 0.50
                   "cautious"    if 0.30 ≤ ρ < 0.50
                   "conflicting" if ρ < 0.30
```

*Confirmed present:* `final_pipeline/02_Proteinomics/platform_tier.py:31, 38–43`.

**Measured tier distribution (Y10, 5,957 proteins tiered, 2026-08-15):**

| Tier | n proteins | fraction | rho range |
|---|---|---|---|
| consistent | 1,581 | 26.5% | [0.50, 0.86] |
| cautious | 2,197 | 36.9% | [0.30, 0.50) |
| conflicting | 2,179 | 36.6% | [−0.47, 0.30) |

Median cross-platform rho: 0.373 (= W_PROT derivation input, confirmed).

**Mathematical status.** The tier is a threshold-based label. The hardcoded
ρ_PROT = 0.373 in the `W_PROT` comment (`n_eff = 2/(1+0.373) = 1.456 ≈ 1.45`)
is a prior, not a runtime measurement; the script measures ρ per protein against
the ProCAN/CCLE shared-line intersection.

**Sampling-noise limitation.** At minimum n_shared = 30 shared lines, the 95%
CI on ρ = 0.30 spans [−0.068, +0.596] (width 0.663, Bonett & Wright 2000,
DOI 10.1007/BF02294183) — wider than the 0.20-unit tier bin. Tier assignment at
the minimum shared-line count is dominated by sampling noise: ~2,000 proteins
sit within ±0.10 of a tier boundary.

**Asymmetric fallback.** Conflicting-tier proteins use ProCAN only (CCLE dropped).
No stated justification for preferring ProCAN; the practical reason is that ProCAN
has 2.5× more screened lines (952 vs 375), which gives noisier per-protein rho
estimates for CCLE. The preference is not documented as a design choice.

---

## 11.5 Shrink-on-sparse and the MNAR direction defect — **[OPEN BLOCKER]**

Both arms apply a conservative shrink when fewer sources than the nominal
`N_SOURCES` contribute to a pair:

```
shrink = sqrt(N_SOURCES / n_observed)    if n_observed < N_SOURCES
       = 1                               otherwise

z_t = z_raw / shrink
```

RNA arm: `N_SOURCES = 3` (`rna_scorer.py:29`; confirmed by X1).
Protein arm: `N_SOURCES = 2` (`protein_scorer.py:33`).
Shared utility: `final_pipeline/utils/common.py:stouffer_combine` lines 158–160.

**MNAR direction defect.** For proteins absent from the CCLE source (Missing Not
At Random — typically low-abundance proteins), the absent z-score would be
strongly *negative* if it were observed. The shrink moves z toward zero, i.e. in
the wrong direction for MNAR cases. Affected pairs: **83.4%** (5,041,853 of
6,043,430 protein-arm pairs have n_sources=1, Y3 measurement 2026-08-15) — the
large majority of the protein arm, not a minority. Earlier documentation (~23%)
was wrong; the correct figure comes from direct measurement of
`bulk_prot_z.parquet`. This is an open architectural decision; it requires
choosing a principled MNAR imputation strategy. Do not correct in code without
a decision from Fiona (three options documented in COWORK v5.0 Y3 escalation).

*Source: X1 sweep, `protein_scorer.py:140–142`, `common.py:158–160`.*

---

## 11.6 has_alteration vs has_driver_alteration

Two Boolean flags exist with different semantics (`driver_routing.py:60–63`):

| Flag | Formula | Line | Use |
|---|---|---|---|
| `has_driver_alteration` | `mut_driver ∨ fusion_driver ∨ has_cna_alteration` | 60–62 | **Production sort key** — ranking |
| `has_alteration` | `p_mutation > 0` | 63 | **Diagnostic-only** — no live consumer |

`has_alteration = (p_mutation > 0)` captures any mutation at all (not the
half-saturation of the Hill score, which gates at 0.5). It is the "any-alteration"
arm planned in the ALTERATION_DEFENCE evaluation but **never implemented**:
`eval.py` docstring (line 8) refers to "any-alt — has_alteration flag boosts rank
within gene class" but no `rank_any` variable exists in the production code (Y5
audit, 2026-08-15). It enters no ranking formula and no evaluation arm; the column
is written to `flags_with_driver.parquet` but has no live consumer. Treat as
diagnostic-only.

---

## 11.7 Driver promotion as sort key — not an estimator

Stage 5 "driver promotion" places `has_driver_alteration = TRUE` lines above
`has_driver_alteration = FALSE` lines within a gene. The implementation is:

```python
sort_key = has_driver_alteration.astype(float) * 1e6 + core_score
```

(`Ranking/cli.py`, same lexicographic embedding as §6.2.) This is **not** an
estimator of dependency probability — it is a deterministic display rule that
ensures driver-altered lines are always returned before non-driver lines. Its
effect is exactly the 2×2 partition the ALTERATION_DEFENCE evaluation measured:
high-tier lines are driver + high-score; low-tier are non-driver + low-score.

Implication for §9 gap 15 (sort key unsupported by measurement): the ranker-region
pAUC of 0.500–0.518 applies to the `core_score` tiebreak *within* each driver
partition. The driver sort key itself is evaluated by the hit@20 delta of Stage 6,
not by pAUC.

---

## 11.8 MAD-z sensitivity threshold convention

The sensitivity threshold applied in `eval.py:39` is:

```python
SENSITIVITY_MAD_Z = -0.5     # code default
```

The project-chosen threshold for the alteration defence analysis is **−1.0**
(from T_sensitivity_threshold diagnostics, W-series). These are different:

- `SENSITIVITY_MAD_Z = -0.5` is the code's built-in default for defining "sensitive" cell lines.
- `-1.0` was the threshold chosen in T_sensitivity_threshold as the optimal
  operating point for the alteration defence analysis, applied by passing it as a
  parameter to the evaluation function.

This is not a discrepancy — the two constants serve different purposes. However,
any citation of hit@20 delta results should specify which threshold was active.
The Stage 6 ALTERATION_DEFENCE results (E1 +0.0238, E2 +0.0138) used threshold
−1.0, not −0.5.

---

## 11.9 Stage 5 selectivity — independence assumption (X3)

The Stage 5 selectivity score for gene A against background gene B is:

```
selectivity(A, B) = score_A × (1 − score_B)
```

This formula is exact only if A and B are **independent** across cell lines —
the product structure does not hold when `score_A` and `score_B` are correlated.
If positively correlated, a cell line that scores high on A also tends to score
high on B, making `(1 − score_B)` small exactly when `score_A` is large, which
deflates the selectivity score for all lines. The functional consequence is that
selectivity overstates the penalty for lines that genuinely express both targets.

**X3 measured (1,000 random gene pairs, seed=99, all 1,000 pairs eligible):**

| Statistic | Pearson r | Spearman r |
|---|---|---|
| Median | **0.0426** | 0.0413 |
| Mean | 0.0509 | 0.0490 |
| SD | 0.0867 | 0.0850 |
| IQR | [−0.001, +0.099] | [−0.002, +0.097] |
| Fraction > 0.2 | **4.9%** | 4.4% |
| Fraction > 0.3 | 0.7% | 0.5% |

Mean shared lines per pair: 1,042 (median 978, range 678–1,493).

L2-only lines (two layers, tighter evidence): median r = 0.009 — more independent
than L1-only (median r = 0.050). Two-layer genes show less cross-gene correlation,
the opposite of what correlation-driven inflation would predict.

Rank displacement (top-20, selectivity formula vs correlation-adjusted): mean 3.51,
median 3.0 positions out of 20. Driven by the ~5% of pairs with r > 0.2.

Pre-declared: median corr > 0.2 → **WRONG**. Measured median is 0.043.

**Conclusion.** Core scores of randomly selected gene pairs are near-independent
across cell lines. The selectivity formula `score_A × (1 − score_B)` is
empirically justified — the independence assumption holds on average. Tail pairs
(r > 0.2, ~5%) contribute modest rank displacement (3–4 positions in top-20);
no systematic correction is warranted.

*Source: `Ranking/cli.py`, selectivity computation; `diagnostics/X3_selectivity_corr_fast.py` (2026-08-14).*

---

## 11.10 Stage 5 co-selection — tie inflation from min() (X4)

The Stage 5 joint co-selection score for a gene pair is:

```
joint_score = min(score_A, score_B, …)
```

with floor `JOINT_FLOOR = 0.50` (`cli.py:172`). The concern is that `min()` is
idempotent only when both scores are equal — when they differ it projects to the
lower value, potentially tying lines that have different profiles.

**X4 measured (100 pairs, seed=99, floor=0.50):**

| Quantity | Measured |
|---|---|
| Joint min() tie rate | median **17.7%** (IQR 14.2–21.9%) |
| Single-gene tie rate | median **16.3%** |
| Excess (joint − single) | median **+0.7 pp** |
| Top-20: lines at min value | mean=1.0, median=1.0 |
| Non-min score spread | median=0.065 |

Pre-declared: excess ≥ 10 pp → **WRONG**. The actual excess is 0.7 pp — the min()
operator creates negligible additional tie inflation relative to single-gene
scoring at this floor. Lines above the floor have broadly similar score values,
making the tie rates comparable.

**Conclusion.** No tie-inflation correction is warranted.

---

## 11.11 Fusion channel dominance (X5)

The Noisy-OR fusion score combines three channels:

```
p_base   = 1 − (1−p_conf)(1−p_ffpm)(1−p_recur)
p_fusion = min(p_base + 0.15·any_in_frame, 1.0)
```

**X5 measured (144,221 fusion records, 64.8% pass gate p_fusion ≥ 0.5):**

| Category | Count | Share of passing |
|---|---|---|
| 0 channels ≥ 0.4 (inframe-boost only) | 9,296 | 10.0% |
| Exactly 1 channel ≥ 0.4 | 51,825 | **55.5%** |
| 2 channels ≥ 0.4 | 25,073 | 26.8% |
| 3 channels ≥ 0.4 | 7,223 | 7.7% |

Among single-channel passes: p_conf dominates — **93.0%** of single-channel passes
occur via p_conf alone. The p_conf channel reflects `max_confidence ∈ {low,
medium, high}` mapped through Hill(conf_ord; p₀=2, k=1.5), with median
p_conf = 0.500 across all fusions (55.8% have p_conf ≥ 0.5).

**Genuine Noisy-OR** (no single channel ≥ 0.5): 10.3% of passes → CONFIRMED < 20%.
**Pre-declared** (>60% single-channel) → **WRONG** at 55.5%.

**Implication.** The fusion channel is effectively a confidence-tier classifier
with p_ffpm and p_recur as secondary modifiers. The ordinal-to-continuous mapping
via Hill(conf_ord; ·) is the structural issue noted in §C.3 — the choice of
encoding {low=1, medium=2, high=3} determines the score for the majority of
passing fusions.

---

## 11.12 CNA attrition waterfall (X6)

Starting from all rows in the `main.cosmic_cna` DuckDB table:

| Step | Rows | Lost |
|---|---|---|
| All rows (Step 1) | 143,052 | — |
| After sample join to pipeline lines (Step 2) | 142,572 | 480 |
| After is_ambiguous = FALSE (Step 3) | 142,412 | 160 |
| After gene_role join (Step 5) | 142,412 | 0 |
| After direction gate (Step 7) | 3,070 | 139,342 |
| Final cna_flags has_cna_alt | 2,925 | 145 |

**Dominant loss: unknown gene role.** Of 142,412 post-ambiguity rows, 138,185
(97%) carry `gene_role = unknown` — they belong to genes not in the curated
580-gene set (`oncogene=256, tsg=254, both=70`). These are structurally
unevaluable; the direction gate applies only to role-annotated genes.

Among role-annotated rows (4,227): 1,157 (27.4%) fail the direction gate
(amplification for a TSG, deletion for an oncogene, or neutral).

**Pipeline coverage.** 992 of 1,840 pipeline lines appear in COSMIC CNA (54%).
The remaining 46% of lines have no CNA data and receive `has_cna_alteration = FALSE`.

**Ploidy artefact (X6, ploidy correlation):**

```
Spearman r(per_line_mean_CN, per_line_alt_rate) = 0.2488   p = 1.97 × 10⁻¹⁵
```

Lines with higher mean copy number have higher CNA alteration rates. Pre-declared
|corr| < 0.15 → **WRONG**.

**Stratified analysis (diagnostics/D1_ploidy_stratified.py, 2026-08-15):**

The pipeline uses COSMIC's pre-computed `cna_call` (ploidy-aware), not raw absolute
CN thresholds. The stratification reveals whether the residual correlation is artefact
or biology:

| cna_call | n lines | r(mean_CN, correct-dir rate) | r(mean_CN, correct-dir n) | r(mean_CN, total calls) |
|---|---|---|---|---|
| amplification | 795 | +0.154 (p=1.2×10⁻⁵) | +0.396 (p=2.7×10⁻³¹) | +0.721 (p=1.1×10⁻¹²⁸) |
| deletion | 993/997¹ | +0.090 (p=4.6×10⁻³) | −0.063 (95% CI [−0.124, −0.001], p=0.048) | −0.452 (p=4.6×10⁻⁵¹) |

¹ D1 n=993 lines; Y9 re-measurement n=997 (slightly different filtering; values consistent).

`correct-dir rate` = fraction of amplifications that are oncogene-direction (or
deletions that are TSG-direction) — **per-gene** fraction. `correct-dir n` = absolute
count of direction-gated events **per line**. `total calls` = raw amp or del call count.

**Interpretation:** The r=+0.72 for total amplification count is driven by
chromosomal instability — highly aneuploid lines have more focal amp events genome-wide,
including non-oncogene loci. For ONCOGENE amplifications specifically (correct-direction
rate), r=+0.154. For TSG deletions: the correct-direction rate (per gene) is r=+0.090,
while the correct-direction count (per line) is r=−0.063 (95% CI [−0.124, −0.001],
p=0.048 — barely excludes zero; marked **UNEXPLAINED**: the CI sign is negative but
the mechanism producing TSG-direction deletion depletion in high-CN lines, independently
of the total deletion depletion at r=−0.45, has not been demonstrated). The D1 value
of −0.056 (p=0.08) was within sampling uncertainty of the Y9 value and is now
superseded.

**Conclusion:** The ploidy correlation in the pipeline is primarily chromosomal-instability
biology, not a threshold artefact. COSMIC's ploidy-aware calls already prevent the naive
absolute-threshold problem (CN>2.5/CN<1.5). The residual r≈0.15 is biologically expected.
No code change required. §9 gap 21 is **resolved by the existing implementation**.

---

## 11.13 p_burden confound (X7)

The mutation score includes `p_burden = Hill(variant_burden; p₀=3.0, k=1.5)` where
`variant_burden = variant_count.clip(upper=5)`. The concern is that
hypermutated lines receive elevated p_mutation regardless of biology.

**X7 measured (632,919 gene×line records, 1,744 lines):**

| Quantity | Measured | Pre-declared | Verdict |
|---|---|---|---|
| Spearman r(per_line_TMB, p_mutation) | **0.159** | > 0.20 | WRONG |
| Hypermutator rate / rest rate | **1.06×** | ≥ 1.50× | WRONG |

TMB distribution: mean=404, median=205, p95=1,392 variants.
Hypermutators (top 5%, TMB ≥ 1,392): 88 of 1,744 lines.
Lineage: colorectal (25), uterus (19), blood (11), ovary (8), lung (6).

**Why pre-declarations are wrong.** `mut_driver` rate is 88.9% in non-hypermutators
and 94.1% in hypermutators — compressed near the ceiling. With 89% baseline nearly
everything is already classified as a driver, leaving little room for hypermutators
to differ (ratio 1.06× vs pre-declared 1.5×). The ceiling is a consequence of the
`any_driver ∨ oncogene_hit ∨ tsg_hit` disjunction driving most records to true.

**Ablation (remove p_burden entirely):** 170,334 / 632,919 records (26.91%) change
from `mut_driver = TRUE` to `FALSE`. All are losses (none gained). This shows
p_burden is a material contributor to driver calls, despite the low TMB correlation.
The mechanism: even with modest TMB (median=205 variants), some lines have
`variant_burden = 5` (clipped ceiling) for particular (gene, line) pairs, which
Hill-saturates p_burden. The saturation is not proportional to TMB in the way the
pre-declaration assumed.

**Conclusion.** p_burden does not produce a genome-wide TMB confound at the
line level. Its effect is gene-specific and ablation-measurable but not
lineage-predictable from TMB alone.

---

## 11.14 Chronos re-run delta (X8)

The Chronos validation was re-run using real Chronos scores
(`chronos_validation.parquet`) in place of negated Project Score
(`chronos_validation_OLD.parquet` = Project Score × −1).

**Distribution comparison:**

| | OLD (Project Score negated) | NEW (real Chronos) |
|---|---|---|
| Genes | 16,866 | 17,226 |
| Mean ρ | −0.0052 | **+0.0103** |
| SD ρ | 0.0686 | 0.0725 |
| p5 | −0.118 | −0.105 |
| p95 | +0.096 | +0.122 |

The old data introduced a systematic negative bias (mean ρ = −0.005 vs +0.010
with real Chronos). This is consistent with using −1 × Project Score: any
positive dependency in Project Score became negative.

**Label changes (16,802 genes in both):**

| Tier change | Count |
|---|---|
| none → weak_positive | 944 |
| weak_negative → none | 633 |
| weak_positive → none | 221 |
| weak_negative → weak_negative (stay) | 547 |
| Total changed | **2,253 (13.4%)** |
| Direction flips (validated↔inverted) | **0** |

No gene moved between validated and inverted — the label changes are in the
weak/none region where the effect size (ρ ≈ 0.04–0.10) crosses the floor in one
direction. The validated gene count rose from 31 to 37 (+6), inverted from 6 to 8
(+2).

*Source: `diagnostics/T12_delta.py` against `src/pipeline/outputs/chronos_validation_OLD.parquet`
and `chronos_validation.parquet`.*

---

## §9 additions (from COWORK v4.0, X1–X8)

The following items extend the open-gaps list of §9:

**Gap 20. W_RNA is calibrated to ρ_implied = 0.25, but measured ρ_avg = 0.548 — RESOLVED 2026-08-15.**
`W_RNA = sqrt(1.431)` now applied in `core_score.py:39`. A/B test confirmed flat hit@20
improves (oncogene +10%, tsg +101%) with no loss in driver hit@20. See §11.3.

**Gap 21. Ploidy artefact in CNA direction gate — PARTIALLY RESOLVED.**
X6 measured r(per_line_mean_CN, per_line_alt_rate) = 0.248 (p = 2×10⁻¹⁵). Initial
interpretation: absolute threshold artefact. Stratified analysis (D1_ploidy_stratified.py,
2026-08-15) shows the pipeline already uses COSMIC's ploidy-aware `cna_call` — not raw
absolute thresholds. For oncogene-direction amplifications, correct-direction rate r=+0.154
per gene (modest). For TSG-direction deletions: correct-direction rate r=+0.090 per gene
(modest, biologically plausible); correct-direction count per line r=−0.063 (95% CI
[−0.124, −0.001], p=0.048, Y9 2026-08-15) — **UNEXPLAINED**: the CI barely excludes
zero on the negative side; the mechanism producing this residual effect after direction
gating has not been demonstrated. The total deletion count correlates at r=−0.45 (ploidy
suppression); the 6× reduction to r=−0.063 after direction gating is consistent with
direction-gating filtering to TSG-specific events. No code change required. See §11.12.

**Gap 22. MNAR shrink-direction defect (protein arm, 83.4% of pairs — updated).**
For low-abundance proteins absent from one source (MNAR), shrink moves z toward
zero when the true absent z would be strongly negative. The direction of shrink is
wrong for this specific missingness mechanism. Affected pairs: 83.4% (5,041,853 of
6,043,430 protein-arm pairs have n_sources=1; Y3 measurement 2026-08-15). Earlier
documentation (~23%) was wrong. Three resolution options escalated to Fiona (COWORK
v5.0 Y3): remove shrink, carry detected as feature, or document as limitation.
Architecture decision needed from Fiona (§11.5).

---

## 11.15 p_mutation distribution and ceiling (Decision 5)

`p_mutation` is computed per (gene, cell-line) pair in `mutations_scoring.py` via:
```
p_base   = 1 − (1 − p_vep)(1 − p_path)(1 − p_burden)     [Noisy-OR]
p_mutation = min(p_base + 0.15 × any_driver, 1.0)
mut_driver = (p_mutation ≥ 0.50)
```

**Measured on `mutations_scores.parquet` (632,919 pairs, 2026-08-15):**

| Statistic | Value |
|---|---|
| mean | 0.698 |
| std | 0.155 |
| min | 0.161 |
| 25th pctile | 0.531 |
| median | 0.710 |
| 75th pctile | 0.822 |
| At ceiling (=1.0) | 1.67% |
| ≥ 0.85 | 21.4% |
| mut_driver (≥ 0.50) | **90.6%** |

For role-annotated genes (oncogene/tsg/both): mut_driver rate = **92.3%**,
ceiling fraction 8.25%.

**Channel decomposition (among mut_driver=TRUE records):**

| Channel | Fraction single-channel |
|---|---|
| Only p_vep ≥ 0.50 | 1.4% |
| Only p_path ≥ 0.50 | 31.1% |
| Only p_burden ≥ 0.50 | 0.3% |
| ≥2 channels (multi-channel) | 2.2% |
| None ≥ 0.50 individually (Noisy-OR combination) | 65.0% |
| Boost-only (p_base < 0.50, passes only via DRIVER_BOOST) | **0.006%** (32 of 573,600) |

**Y2 finding — DRIVER_BOOST is effectively a dead code path.** Of 573,600 driver
calls, only 32 (0.006%) are boost-only passes that would not survive without the
`any_driver` annotation. Pre-declared >5%: WRONG by a factor of ~800. All 32 pairs
have p_base ∈ [0.35, 0.50) with any_driver=True; the near-gate region (p_base ∈
[0.35, 0.50)) contains 29,781 pairs of which only 32 are driver-annotated. The gate
passes on evidence, not annotation, for >99.99% of driver calls.

**Why 65% pass without any single channel ≥ 0.50:** Noisy-OR sums complementary
evidence. With p_vep=0.40, p_path=0.40, p_burden=0.40 (each below threshold),
p_base = 1 − 0.6³ = 0.784 >> 0.50. The gate is not binary-OR — it rewards
moderate multi-channel evidence. This is the intended architecture.

**Why 90.6% overall is high:** The input (`mutations_collapsed.parquet`) is
pre-filtered to non-synonymous, functionally annotated variants. The population
already excludes synonymous and intronic mutations. The 90.6% rate reflects that
most remaining mutations are pathogenic at some level in this filtered set, not
a broken gate.

**The 9.4% that fail the gate (p_mutation < 0.50)** have p_vep=1 (synonymous
VEP rank ~1–2), low CADD/REVEL score, and fewer than 3 variants — records that
passed the initial variant-type filter but have very weak functional evidence.

**Conclusion.** The gate is not architecturally broken. p_mutation operates as
a soft filter on pre-filtered mutations; the effective discrimination is on the
continuous score (0.16–1.0) used in ranking, not on the binary 90.6% pass rate.
No code change is applied for this decision (see §9 gap 23 for the open question
on whether to raise the gate to 0.65 for selectivity).

---

## §9 gap 23 (from D2 analysis, 2026-08-15)

**Gap 23. p_mutation gate selectivity: 90.6% pass rate on pre-filtered mutations.**
The binary `mut_driver` flag captures 90.6% of all (gene, line) pairs with any
observed non-synonymous mutation. The continuous p_mutation score (mean=0.70, std=0.16)
still varies and feeds the ranking. Open question: would raising the gate to 0.65
or using the continuous score directly as a ranking key (rather than binary promotion)
improve hit@20? Requires a Stage 6 A/B test against the current driver-promotion
architecture.

---

## 11.16 Code ownership boundary (Y1 audit, COWORK v5.0, 2026-08-15)

The production pipeline contains two categories of code:

**Fiona's independent implementations** — protein_scorer.py, platform_tier.py,
core_score.py, all of Stage 3 (altercations), Stage 4b (driver routing), Stage 4c
(confidence tiers), Stage 5 (ranking/CLI), Stage 6 (eval). These are not in
`02_transcriptonomics.ipynb` in any form.

**Adapted from the shared transcriptomics notebook** (`src/pipeline/02_transcriptonomics.ipynb`) —
`final_pipeline/utils/common.py::robust_z_matrix` and
`final_pipeline/01_Transcriptomics/rna_scorer.py`. Key simplification: the
notebook's `robust_z` returns `(z, reason)` with three guards (too-few-peers,
gene-silent-in-lineage, no-spread/MAD-floor); production `robust_z_matrix` retains
two guards (too-few-peers, MAD-floor) and drops the gene-silent-in-lineage guard
and the reason return entirely.

**Features absent from production:**

| Feature | Notebook function | Production status |
|---|---|---|
| `calibrate_constants()` | Per-gene MAD_FLOOR/MIN_PEERS derivation | ABSENT — frozen constants |
| Silent-lineage guard | `SILENT_FRAC=0.20` threshold | **RESTORED 2026-08-15** — added to `common.py:score_source_lineage()` after `robust_z_matrix`; sets Z column to NaN where < 20% of lineage lines express above 1.0 log2 TPM+1 |
| `z_reason` return | `(z, reason)` tuple | ABSENT — matrix only |
| Cochran's Q / I² | `heterogeneity()`, `I2_CUT=50` | ABSENT |
| FDR control | `false_discovery_control()` | ABSENT |
| `source_floor_check` | Per-run detection from data | ABSENT |
| `exclusion_gate` | 3-valued pass/fail/unknown | ABSENT |
| `shrink(z, k) = z×k/(k+1)` | EB shrink with posterior weight k | PROTEIN SHRINK DIFFERS: √(n_obs/N) |
| `collapse_isoforms` | Isoform deduplication | ABSENT |

**z_reason gap (T10/D26) — CLOSED 2026-08-15:** The silent-lineage guard fires
when fewer than 20% of a lineage expresses a gene above 1.0 log2 TPM+1 (EXPRESSED_MIN).
Previously omitted, this meant a line with near-zero expression received
`z = (value − median) / MAD_FLOOR` rather than NaN — an absence-as-signal error.
The guard is now present in `common.py:score_source_lineage()`, restoring the
notebook's third guard. Affected (gene, lineage) pairs produce NaN z-scores, which
propagate to NaN core_score; these orphan pairs are excluded from the hit@20 ranking
by the `core_score.notna()` filter in `eval.py`.

**Effect of guard restoration:** ~87K additional (gene, model) pairs moved from
"scored with uninformative z-score" to "unscored NaN" at the 2026-08-15 run.
Total unscored pairs at the final run: 364,470 (of which 87K are silence-guard-masked;
the remainder are alteration-only orphans with no expression measurement at all).

*Consequence for dissertation:* The adaptation decisions (which guards to retain)
are Fiona's; the notebook is the teammate's. No production file should be described
as containing the notebook's methodology without noting the simplification. The
silence guard restoration is documented in the §C.2 change and in ALTERATION_DEFENCE §13.