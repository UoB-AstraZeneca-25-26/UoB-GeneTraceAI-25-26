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

Line numbers are approximate (notebooks renumber cells on edit); function and
file names are exact. Paths are relative to the repo root
(`C:\Disertation\UoB-GeneTraceAI-25-26`).

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

## 2.3 Noisy-OR combination (superseded — retained because its failure motivates §2.4)

```
core_score = 1 − (1−E)(1−P)     both layers present
core_score = E                   E only
core_score = P                   P only
core_score = NaN                 neither
```

**Model — derivation.** Let `A` be the latent event "the gene product is
functionally abundant in this cell line." Model each layer as an independent
noisy channel that fires with probability equal to its normalised score. `A`
occurs iff at least one channel fires, so

```
P(A) = 1 − P(no channel fires) = 1 − ∏ᵢ (1 − pᵢ)
```

The independence assumption enters at exactly one step: factorising
`P(no channel fires)` into a product.

**Assumptions.** (i) `E` and `P` are calibrated probabilities of `A`, not merely
scores in [0,1]. (ii) The channels are conditionally independent. (iii) `A` is
binary, i.e. abundance is a presence/absence event rather than a continuum.

**Failure — stated as a theorem.** Assumption (ii) is false: mRNA and protein
abundance correlate at ρ ≈ 0.46–0.58 (Nusinow et al. 2020). Two independent
consequences follow.

*Non-idempotence.* From §0.2, when the layers agree (`E = P = x`):

```
N(x,x) − x = 2x − x² − x = x(1−x)  >  0   for all x ∈ (0,1),   max 0.25 at x = ½
```

So the operator inflates *most* precisely where the score is least informative,
and inflates *nothing* at the endpoints. Agreement between correlated layers is
rewarded as though it were corroboration from independent sources.

*Fréchet position.* From §0.2, the true `P(A∪B)` under positive dependence lies
below the independence point and toward `max(p,q)`. Noisy-OR is stuck at the
independence point, so its bias is upward and grows with dependence.

**Observed on production data, matching the prediction.** 22.5% of two-layer
rows scored > 0.95; **379 rows pinned at exactly 1.0**. The pinning at exactly
1.0 is the boundary case of the same theorem: any cell with `E = 1` or `P = 1`
(the top-ranked cell line for that gene in either assay) is mapped to 1
regardless of the other layer, because `1 − (1−1)(1−P) = 1`. Noisy-OR is
absorbing at 1.

*Provenance:* `scoring_variants.combine(method="noisy_or")`
([scoring_variants.py:59-60](../src/pipeline/scoring_variants.py)); original
`02_core_score.ipynb` cell. Superseded by §2.4; retained in the validation grid.

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

## 2.5 Alternative normalisation / combination variants (validation harness only)

`scoring_variants.py` implements a 4×3 grid (12 variants) used *only* to compare
empirically against GDSC ground truth. Only percentile + idempotent-mean (§2.1 +
§2.4) is in production.

| Normalisation | Formula | Mathematical character |
|---|---|---|
| `minmax` | `(x − min)/(max − min)`, span 0 → NaN | Affine, **not** outlier-robust: a single extreme value compresses all others toward 0. Sensitivity to the extremes is unbounded as the span shrinks. Flagged weakest by audit — this is why. |
| `percentile` | `rank(pct=True)` | PIT; see §2.1 |
| `zscore` | `z = (x−μ)/σ`, then `σ(z) = 1/(1+e^{−z})` | Standardisation is affine (so location/scale-invariant but not monotone-transform-invariant); the logistic squash is a smooth monotone bijection ℝ→(0,1) with maximum slope ¼ at `z=0`. Preserves magnitude information that percentile discards, at the cost of assuming a roughly symmetric input. |
| `hill` | percentile, then `pct^k/(pct^k + p₀^k)`, `k=2, p₀=0.5` | A logistic in `ln(pct)` — see §C.1. Marked uncalibrated in-source. |

| Combination | Formula | Mathematical character |
|---|---|---|
| `noisy_or` | `1 − (1−E)(1−P)` | Probabilistic-sum t-conorm; independence point of the Fréchet interval; non-idempotent (§2.3) |
| `weighted_sum` | `½E + ½P` | Idempotent; equals the production formula (§2.4) |
| `product_floor` | `max(E·P, ½·min(E,P))` | Monotone (pointwise max of two monotone functions) and satisfies the identity `T(a,1) = max(a, ½a) = a`, but **not** idempotent (`T(x,x) = max(x², ½x) < x` for `x<1`) and **not** associative, so it is not a t-norm. It has a kink where `E·P = ½·min(E,P)`, i.e. at `max(E,P) = ½`: below that the floor binds, above it the product does. Its position in the Fréchet interval — and hence the dependence structure it implicitly assumes — is unanalysed. |

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

## 2.7 Chronos three-layer blend (tested, rejected)

```
chronos_pct   = 1 − rank(chronos_raw, pct=True)      # invert: more negative raw = more essential
core_score_3  = w_EP · core_score_2layer + w_C · chronos_pct
```

Built only on a 50/50 train split of cell lines (to avoid circularity), gated per
gene by `ρ(chronos_pct, ground_truth) > 0.05` (§C.4).

**Note on the inversion.** `1 − F̂ₙ(x)` is the survival function, and equals
`F̂ₙ` of the negated variable up to tie handling — so inverting the percentile is
identical to ranking `−chronos_raw`, and is an order-reversing involution
(§5.1).

**Note on the split.** The 50/50 split over *cell lines* controls circularity
only for cell-line-level leakage. Gene-level information (which genes were
inspected when choosing the blend) is not controlled by it.

**What the measurement showed.** `ρ_EC = ρ_PC = 0.005` — Chronos is essentially
uncorrelated with the abundance layers. By §0.3 this is exactly the condition
under which a third layer is *most* informative (it raises `n_eff` from 1.37 to
2.28 and cuts variance from 0.73σ² to 0.42σ²), which is why the blend was worth
testing and why the GLS optimum up-weights it to 0.42.

**Outcome: rejected.** BCL2's hit@20 collapsed to 0 when Chronos was blended as
a third layer. Note the tension this creates: variance reduction and predictive
accuracy point in opposite directions here, which means the equal-variance
unbiasedness assumption of §0.3 fails — Chronos is *not* an unbiased estimator
of the same latent quantity as E and P. It measures dependency, not abundance.
Not used in production.

*Provenance:* `src/pipeline/build_chronos_layer.py:113-137`; re-derived in
`validate_chronos.py:54-77`.

---

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
p_vep    = Hill(max_vep_rank,      p₀=3.0, k=2.0)    # 0-3 ordinal VEP impact
p_path   = Hill(max_pathogenicity, p₀=0.5, k=2.0)    # max(AlphaMissense, REVEL)
p_burden = Hill(variant_burden,    p₀=3.0, k=1.5)    # variant_burden = log1p(variant_count)

p_base     = 1 − (1−p_vep)(1−p_path)(1−p_burden)     # Noisy-OR
p_mutation = min(p_base + 0.15·driver_flag, 1.0)     # DRIVER_BOOST = 0.15 (locked)

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

**`p_burden` is doubly redundant.** `variant_burden = log1p(variant_count)`
already encodes count, and the max-aggregation of §C.4 encodes count a second
time. Excluding burden from `p_mutation_quality` is a partial acknowledgement of
this.

**The `+0.15` driver boost is not a probability operation.** Additive shift on a
probability scale is not closed on [0,1] — hence the `min(·, 1.0)` clamp, which
introduces a point mass at exactly 1.0 for every driver hit with
`p_base ≥ 0.85`. The scale-appropriate version is an additive shift on the
log-odds (a Bayes factor):
`logit(p') = logit(p) + λ·driver_flag`, which is closed on (0,1), needs no clamp,
and has the interpretation "a driver hit multiplies the odds by `e^λ`."

*Provenance:* max-aggregated to one row per (ensg_id, model_id) →
`mutations_scores.parquet`. Quality is kept separate from burden because a
downstream signature-discount step scales burden but not quality:
`p_mutation_adj = NoisyOR(p_mutation_quality, m_mut · p_mutation_burden) + driver_boost`.

## C.3 Fusions scoring

```
conf_ord   = {low:1, medium:2, high:3}[max_confidence]
ffpm_pctl  = rank(best_ffpm,    pct=True) within ensg_id
recur_pctl = rank(fusion_count, pct=True) within ensg_id

p_conf  = Hill(conf_ord,   p₀=2.0, k=2.0)     # midpoint at "medium"
p_ffpm  = Hill(ffpm_pctl,  p₀=0.5, k=1.5)
p_recur = Hill(recur_pctl, p₀=0.5, k=1.5)

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

**Severe tie mass invalidates the percentile transform.** In-source note: 87% of
`fusion_count` mass sits at percentile 1. By §0.1c the PIT does not apply —
`recur_pctl` is not uniform, it is a two-atom distribution. `Hill(1.0; p₀=0.5,
k=1.5)` = 0.739 for 87% of rows, so `p_recur` contributes almost no
discriminating information while contributing a near-constant 0.739 to the
Noisy-OR product — which by itself floors `p_base` at 0.739 for the great
majority of fusion rows, before any other evidence is considered. The gentle
`k=1.5` mitigates but does not fix this; the correct handling is to treat
recurrence as the near-binary variable it actually is.

*Provenance:* max-aggregated to (ensg_id, model_id) → `fusions_scores.parquet` →
Stage 4. In-frame fusions are more likely to produce a translated chimeric
protein, hence the boost. No quality/burden split — the whole score is scaled by
the signature discount (`s_fus' = m_fus · s_fus`).

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

**Consequence.** This is currently harmless because `p_mutation` and `p_fusion`
are display-only (see the Track C table above). It is a hard blocker on ever
promoting this layer into ranking, since the bias correlates with gene length
and mutation rate — exactly the confounds a driver score is supposed to control
for. The bias-free alternatives are an explicit maximum-of-`m` correction
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

## C.6 Validated / Inverted / None classification — **uncorrected multiplicity**

```
RHO_THRESHOLD = 0.1 ;  P_THRESHOLD = 0.05 ;  MIN_N = 30
"validated"  if ρ >  0.1 ∧ p < 0.05
"inverted"   if ρ < −0.1 ∧ p < 0.05
"none"       otherwise
```

**Model.** A two-sided per-gene test of `H₀: ρ = 0` combined with a minimum
effect size. `MIN_N = 30` is justified as the validity floor for the Spearman
t-approximation (§0.4), not as a power calculation.

**The multiplicity problem, quantified.** This rule is applied independently to
every unknown-class gene. Under the global null, the expected number of genes
called "validated" is

```
E[false validations] = m · α/2 ≈ 0.025·m      (one-sided, positive direction)
```

At `m ≈ 3000` unknown-class genes that is **≈ 75 spurious "validated" calls**,
each of which is then promoted from `confidence = "low"/"unknown"` to
`"moderate"` in Stage 5. The `ρ > 0.1` effect-size floor removes some of these
but is not a substitute for error control: at `n = 900`, `ρ = 0.1` has
`p ≈ 0.003`, so the effect floor and the p-floor are close to redundant and the
binding constraint is `p < 0.05` alone.

**No correction is implemented** (verified by search over `src/` for
`multipletests`, `benjamini`, `fdr_bh`, `bonferroni` — no matches). The
appropriate procedure is Benjamini–Hochberg at a stated FDR `q` (§0.6): the
per-gene tests share a ground-truth vector and are therefore positively
dependent, satisfying PRDS, so plain BH is valid. Reporting `π̂₀` alongside would
directly estimate what fraction of the current "validated" set is null.

*Provenance:* `build_chronos_validation.py:41-43, 89-92`; identical thresholds in
`test_run_metabolomics_mirna.classify():43-45, 100-105`. Written to
`chronos_validation.parquet` (`ensg_id, n, rho, pval, chronos_check`), consumed by
`build_full_predictions.py` (lines 82-89, 127-133, 158-168).

## C.7 Chronos correlation gate

```
use_chronos = ( ρ(chronos_pct, ground_truth) > 0.05 )      # per gene
```

A per-gene safety gate: only add Chronos where it aligns with drug sensitivity.
Note this is an *unthresholded-in-p* gate — it uses effect size only, so at
typical `n` a substantial fraction of genes pass by chance (`P(ρ > 0.05) ≈ 0.07`
at `n = 900` under H₀ using `√(n−1)ρ ~ N(0,1)`). It is a filter, not a test, and
should be described as such.

*Provenance:* `validate_chronos.chronos_correlation_gate():121-128`.

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

## C.10 Constant-term global-scalar blend (rejected candidate layers)

```
global_pct = rank(mean_z_across_all_metabolites_or_mirnas, pct=True)     # per cell line
candidate  = 0.5·baseline_core_score + 0.5·global_pct
```

**Model.** Adding a *gene-agnostic* term to a gene-specific score. Note the
structural consequence: because `global_pct` is constant across genes within a
cell line, it cannot change the within-cell-line ordering of genes at all, and it
shifts every gene's score for that cell line by the same amount. Its only effect
is on the *across-cell-line* ranking that the pipeline actually uses — which is
precisely why it is a lineage-confound risk rather than a biological signal: a
scalar that predicts drug response identically for every gene is, by
construction, not a gene-specific mechanism.

**The confound being probed.** A scalar can look predictive purely by proxying
tissue-of-origin, which is itself strongly linked to certain drug responses
(e.g. BCL2 inhibitors in blood/lymphoid lineages). §C.12 tests exactly this.

*Outcome:* rejected for both modalities.

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

## C.13 Metabolomics enzyme-flux proxy (tested, null)

```
proxy = product_metabolite_percentile − substrate_metabolite_percentile
```

**Model.** A difference of two percentiles as a proxy for net flux through an
enzyme — a steady-state argument: at higher enzyme activity, product accumulates
relative to substrate.

**Why the difference of percentiles is a weak construction.** The difference of
two uniform variables is *triangular* on [−1,1], symmetric about 0, with variance
`(1 + 1 − 2ρ)/12 = (1−ρ)/6`. Two consequences: (i) the proxy is not on the same
[0,1] scale as everything else in the pipeline, and (ii) when substrate and
product are positively correlated (which they are, being co-measured metabolites
in the same pathway), the variance collapses toward 0 and the proxy becomes
mostly noise. This predicts a null result on statistical grounds independent of
the biology.

*Outcome:* null — no gene cleared `ρ > 0.1, p < 0.05`. Not used in production.

*Provenance:* `test_run_metabolomics_mirna.py:244-250` (`METABOLITE_MAP`);
literature-curated substrate/product pairs per gene (e.g. NAMPT:
niacinamide → NAD⁺), tested against Chronos essentiality via the §C.6 rule.

## C.14 miRNA-as-inhibitor dampening (tested, mixed/weak)

```
burden      = mean(expression of curated inhibitory miRNAs)     # per cell line
burden_pct  = rank(burden, pct=True)
dampener    = 1.0 − ALPHA · fillna(burden_pct, 0.0)             # ALPHA = 1.0
E_effective = E_raw · dampener
```

**Model.** Multiplicative repression — the mass-action form for a competitive
inhibitor, where available message scales by the fraction not bound.

**Properties and a boundary problem.** With `ALPHA = 1.0` the dampener spans
`[0, 1]` and *reaches exactly 0* for the top-ranked cell line (`burden_pct = 1`),
annihilating its expression signal entirely regardless of `E_raw`. The
multiplicative form also means `E_effective` is no longer marginally uniform —
it is a product of two dependent uniforms, concentrated toward 0, so it is no
longer on the same scale as the untreated `E` it is compared against. Re-ranking
after dampening (or `ALPHA < 1`) would avoid both issues.

*Outcome:* mixed/weak — MET improved slightly, BCL2's already-odd baseline
direction shrank toward zero. Not adopted.

*Provenance:* `test_run_metabolomics_mirna.py:174-187`. Models curated,
experimentally-known miRNA regulators (e.g. BCL2 ← miR-34a/15a/16) for the two
genes (BCL2, MET) with both a curated miRNA map and GDSC ground truth.

---

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
- `DRIVER_BOOST = 0.15`, `INFRAME_BOOST = 0.15` — locked constants, additive on a
  probability scale (§C.2 argues for log-odds instead).
- `ALPHA = 1.0` in the miRNA dampener — the value at which the dampener reaches
  exactly 0 (§C.14).

**Tested and explicitly rejected**
Chronos as a third core_score layer (§2.7), metabolomics global scalar, miRNA
global scalar, metabolomics enzyme-flux proxy, miRNA-inhibitor dampening —
all documented with reasoning in [DECISIONS.md](../DECISIONS.md), and each with a
statistical account above of *why* the null result was structurally likely.

---

# §9 — Open mathematical gaps

Recorded honestly, in rough order of how much they affect conclusions.

1. **No multiple-testing correction anywhere** (§0.6, §C.6). At `m ≈ 3000`
   unknown-class genes, ~75 spurious "validated" calls are expected under the
   global null, each promoting a gene's confidence tier in Stage 5. Fix: BH at a
   stated FDR; report `π̂₀`.
2. **No null model for hit@k** (§6.3). Raw hit rates are reported without the
   hypergeometric baseline `k/N ≈ 2.2%` or an enrichment ratio. Fix: report
   observed/expected and an exact hypergeometric p-value.
3. **Noisy-OR is still used in Track C** (§C.2, §C.3) on features that are
   demonstrably dependent (VEP consequence and computed pathogenicity), despite
   the same failure having been diagnosed and fixed in Stage 2. Currently masked
   by those scores being display-only.
4. **Max-aggregation extreme-value bias** (§C.4) — `E[max] = m/(m+1)` makes
   multi-variant genes score higher mechanically. Blocks promotion of Track C
   scores into ranking.
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
