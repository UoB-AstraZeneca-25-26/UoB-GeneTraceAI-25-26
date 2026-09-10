# GeneTraceAI Scoring Pipeline — Full Methods Description

This document describes, in dissertation-methods-chapter register, the statistical
methodology implemented in the four core scoring stages of the GeneTraceAI pipeline:
RNA expression scoring, protein abundance scoring, their combination into a single
`core_score`, and the downstream driver-gated routing and confidence-tier assignment.
Every method described below was confirmed by direct reading of the live source files
listed at the head of each section; no step is included that could not be located and
quoted in code. Design choices that were tested empirically within this project are
described together with the evidence that supports them. Parameters that are
project-specific calibrations rather than externally derived quantities are labelled
explicitly, consistent with this project's existing PROVEN / CONTESTED / UNVALIDATED
discipline (see `docs/RISK_REGISTER.md`, `docs/MATH_REFERENCE.md`).

**Verification note.** Two claims underlying the diagrams that accompany this document
were independently re-checked against the live code before this document was written:
(i) that the protein platform-combination step is coverage-weighted (weighted by how
many cell lines a gene was measured in on each platform) as a mechanism *distinct from*
the confidence-weighting that follows it, and (ii) that the RNA scorer's minimum-group-size
gate (`MIN_PEERS`) excludes a lineage's contribution outright, independent of detection
fraction, rather than falling back to a pooled estimate. Both were confirmed accurate
against the code exactly as described in §1 and §2 below.

**Update (full mathematical audit + follow-up task).** A companion document,
`docs/WHY_REFERENCE.md`, now supersedes this document's parameter-by-parameter
justification for every clip bound, shrinkage constant, and threshold across all three
scoring files — this document remains the reference for *how the pipeline works*, that
one for *why each number is what it is*, with explicit PROVEN/CONTESTED/UNVALIDATED
tagging per constant. Two corrections from that task, recorded here so this document is
not left stale:
1. **`W_RNA`/`W_PROT` correlation estimator.** §3.2's `ρ̄ = 0.548` (RNA) and
   `ρ = 0.373` (protein) are **both measured via Spearman correlation**, not a
   Pearson/Spearman split across the two arms — confirmed directly against
   `diagnostics/X2_rna_rho.py` (RNA, one-time, not part of `run_all.py`) and
   `final_pipeline/02_Proteinomics/platform_tier.py:97` (protein, live, run every
   pipeline execution as part of Stage 2). Neither value is read back into `core_score.py`
   at runtime; both `W_RNA`/`W_PROT` remain frozen historical literals. The Pearson
   correlation that genuinely exists live in `core_score.py` (§3.2's `rho_raw`, feeding
   `RHO_PRIOR`/`beta_g`) is a *different* quantity — the RNA-vs-protein correlation used
   for residualization, unrelated to the Kish weights. See `WHY_REFERENCE.md` §3.5.
2. **`n_layers=1` fallback scale mismatch (finding "C7") — PROMOTED 2026-08-29.** See
   the updated §3 text above and `WHY_REFERENCE.md` §3.6 for the full verification and
   promotion record.

**Update (promotion task, 2026-08-29).** Both C7 and the `beta_g` SD-shrinkage fix
(finding "C3") are now live in `core_score.py`, following this project's established
backup/replace/run/regenerate-downstream/spot-check procedure. §3 above reflects both.
C3's full downstream check — left undone in the prior task — came back flat-to-slightly-
improved, not worse; see `WHY_REFERENCE.md`'s Part C summary table for the numbers and
for why this does not count as a fifth instance of this project's separate "upstream
improvement fails to compose downstream" pattern.

---

## 1. RNA Expression Scoring

**Source:** `final_pipeline/01_Transcriptomics/transcriptomics.py`

### 1.1 What it does

The RNA scorer produces one z-score per gene, per cell line, from three RNA sources:
DepMap RNA-seq (`log2(TPM+1)`), HPA RNA-seq (TPM), and GEO microarray data. GEO is not
treated as a single source: samples are split into one *method* per combination of
microarray platform and processing pipeline (`build_method_axes()`, L277–350), because
two samples processed on the same platform by, for example, MAS5 versus RMA are not on
a comparable numeric scale (`classify_processing()`, `PROC_RULES`, L74–86, a regular-expression
classifier over each sample's `data_processing` metadata field).

For each method, replicate samples belonging to the same cell line are collapsed to a
single value by the median (`collapse_replicates()`, L504–522), preventing a heavily
replicated line from dominating its lineage's summary statistics. Each method's
measurement scale (logarithmic or linear) is then confirmed empirically from its own
value distribution rather than assumed from metadata (`detect_scale_from_values()`,
L529–539); the code explicitly logs a warning when the declared SOFT-file scale
conflicts with what the data show (`profile_method()`, L562–567), and trusts the data.
A per-method detection floor is then set: RNA-seq sources use a fixed floor at the
established convention of ≥1 TPM (`TPM_FLOOR = 1.0`, attributed in-line to Uhlén et al.
(2015) — Human Protein Atlas's own detection convention), while microarray methods
receive an empirically derived floor from the lower quantile (`floor_q`) of that
method's own values (L570–577), because microarray intensities have no shared absolute
zero comparable across platforms.

Within each lineage, and for each gene, the fraction of that lineage's cell lines with a
value at or above the detection floor (`det`) determines which of three scoring regimes
applies (`method_z()`, L604–699):

- **Well detected** (`det ≥ upper_cut`): a robust z-score, `(x − centre) / scale`, where
  `centre` is the median and `scale` derives from the median absolute deviation (MAD),
  computed within the lineage.
- **Partially detected** (`lower_cut ≤ det < upper_cut`): a rank-based inverse-normal
  transform of the within-lineage ranks (`rankdata`, then `norm.ppf((R − 0.5)/n)`,
  L663–667) — the Blom/van der Waerden form of the normal-scores transform. Because
  `rankdata` with `nan_policy="omit"` gives tied values a shared mid-rank, cell lines
  censored at the detection floor receive a common rank rather than an arbitrarily
  broken tie, so a censored value reads as a lower bound rather than a spurious precise
  estimate.
- **Below detection** (`det < lower_cut`): the gene is left unscored for that lineage
  (`NaN`, status `below_detection`) rather than assigned `z = 0`, which would falsely
  assert "typical expression" for a gene that was not observed.

Within the well-detected regime, the lineage's own MAD-based scale is *shrunk* toward
the method-wide MAD by a weight `w = nb / (nb + n0)`, where `nb` is the number of
scored cell lines in that lineage and `n0` is a calibrated shrinkage constant
(L653–661); smaller lineages are pulled further toward the pooled, method-wide scale,
and the lineage's own median is used as the centre only when the lineage has at least
`MIN_PEERS` scored cell lines, falling back to the method-wide median otherwise
(`use_lineage_centre`, L659–660).

A distinct, unconditional gate sits beneath this branching: if a lineage's group size
for a given method, `nb`, falls below `MIN_PEERS` (imported from `utils/common.py`,
value 5), **every** gene in that lineage for that method is set to `NaN` with status
`no_data` (`STATUS[3]`), regardless of what its own detection fraction would otherwise
have implied (L673–697, quoted and verified in the Phase A note above). The code
comment is explicit that this is deliberate: rather than compute a lineage-level
statistic — even one heavily shrunk toward the global scale — from too few lines, the
gate discards the estimate entirely, because a genuinely elevated numerator measured
against a near-silent global background can still produce an extreme (`|z| ≥ 10`)
value from what is effectively a single data point. This gate is scoped per method, per
lineage: it removes that one method's contribution for that lineage, not necessarily
the gene's final RNA score, since another method with sufficient coverage of the same
lineage can still contribute through the cross-source combination described next.

Once every method has produced its per-gene, per-lineage z-scores, GEO's individual
platform × processing methods are combined into a single GEO-source z-score by a
sample-size-weighted average (`combine_sources()` stage 1, L746–762; `stouffer()`,
L702–714), with weight `√n` for each method's number of contributing samples. The three
top-level sources — DepMap, HPA, and combined GEO — are then combined at **equal**
weight (`combine_sources()` stage 2, L768–771: `stouffer(Zs)` with no weight argument,
so each finite source contributes weight 1). The module docstring is explicit about why
this two-stage design was chosen over a single flat combination across every method:
"With GEO split on platform × processing there may be 10+ GEO methods; flat Stouffer
would hand GEO ~80% of the weight." The output is the final, lineage-conditioned RNA
z-score per gene, per cell line (`z_lineage` in `chunk_to_frame()`, L875–922).

### 1.2 Why this method was chosen

Each design choice addresses a specific failure mode of treating heterogeneous RNA
measurements as directly comparable:

- **Splitting GEO by platform × processing** is necessary because microarray
  normalisation pipelines (MAS5, RMA, gcRMA, etc.) place values on genuinely different
  numeric scales; pooling them before scoring would mix incommensurable units.
- **Empirical, data-driven scale detection** (rather than trusting SOFT-file metadata)
  guards against mislabelled or inconsistently annotated GEO submissions, a known
  reliability problem for public microarray repositories.
- **The three-regime detection-fraction branch** exists because a single functional
  form cannot correctly handle both well-measured and largely-absent genes. A robust
  z-score assumes enough finite, non-censored values to estimate a location and scale;
  a rank-based transform tolerates a substantial fraction of tied/censored values at the
  cost of only ordinal (rather than metric) information; and a gene mostly below the
  detection floor in a lineage should not receive *any* numeric estimate, since neither
  a location nor a scale is estimable from mostly-censored data. Encoding all three as a
  single always-robust-z estimator would either produce spuriously extreme values for
  low-coverage genes, or silently discard information for well-covered ones.
- **Median/MAD as the robust location–scale pair**, rather than mean/SD, is the standard
  choice when outliers or heavy-tailed measurement noise are expected, which is the
  case for expression data with occasional extreme values; this project uses
  `1.4826 × MAD` as the scale estimator, the standard consistency correction for the
  MAD's efficiency at the Gaussian (Rousseeuw & Croux, 1993).
- **Shrinkage of the lineage-level scale toward the method-wide scale** (`w = nb/(nb+n0)`)
  is a partial-pooling design directly in the spirit of empirical-Bayes/James–Stein
  estimation: it lets small lineages borrow statistical strength from the whole method's
  distribution rather than trusting an unstable, high-variance estimate from a handful
  of lines, while allowing large lineages to rely mostly on their own data. The specific
  value of `n0` used in production was not hand-set but chosen by `calibrate_params()`
  (L815–868): a grid search over `(floor_q, n0, upper_cut, lower_cut)` combinations,
  scored by held-out cross-source rank concordance (Spearman ρ between the scored
  source and a held-out reference source, in the ambiguous detection band only), with
  the final choice made by the conservative end of the "1-SE rule" — selecting, among
  parameter combinations within one standard error of the best-scoring combination, the
  one favouring more rank-based (conservative) rather than more parametric scoring
  (L858–864). The 1-SE rule for model-selection-under-noise is a named technique from
  Breiman, Friedman, Olshen & Stone (1984).
- **Equal weighting of the three top-level sources**, rather than a flat combination
  across all methods, is a deliberate choice to prevent GEO's greater method count from
  dominating the combined estimate purely as an artefact of its finer partitioning, not
  because it carries more independent biological information than DepMap or HPA. This
  reflects a "one source, one vote" principle distinct from within-GEO combination,
  which does weight by sample size because, within GEO, sample count genuinely reflects
  more independent replication of the same measurement modality.

### 1.3 What it replaced, and why the change was made

The current scorer, `transcriptomics.py`, replaced an earlier, simpler scorer,
`01_Transcriptomics/rna_scorer.py` (retained in the repository; no longer invoked by
`run_all.py`). The earlier scorer:

- treated GEO as a single, undifferentiated source (`common.fetch_geo()`), with no
  platform × processing splitting;
- applied a single scoring regime throughout — a plain robust z-score via
  `common.robust_z_matrix()` — with no detection-fraction branching, no rank-based
  treatment of partially detected genes, and no per-method scale or floor calibration;
- combined the three sources with a single flat, sample-size-weighted Stouffer
  combination (`common.stouffer_combine()`), including a "conservative penalty for
  sparse-source pairs" shrinkage term explicitly noted in that function's own comment
  as "not yet theoretically derived" — the same flat-weighting design that
  `transcriptomics.py`'s own docstring identifies as handing GEO's many methods a
  disproportionate share of the combined weight.

The replacement therefore directly addresses two identified problems: incommensurable
GEO measurement scales being pooled as if comparable, and a single point-estimate
z-score being computed even for genes with too little detected signal to support one.
`utils/common.py`'s primitives (`robust_z_matrix()`, `stouffer_combine()`) remain live
and are still used elsewhere in the pipeline (e.g., inside `core_score.py`'s silence
guard, and by the now-superseded protein scorer described in §2.3); they were not
removed, only replaced as the RNA-scoring entry point.

### 1.4 Statistical basis and citations

| Technique | Citation |
|---|---|
| Median Absolute Deviation as a robust scale estimator, `1.4826× MAD` consistency correction | Rousseeuw, P.J. & Croux, C. (1993). "Alternatives to the Median Absolute Deviation." *Journal of the American Statistical Association*, 88(424), 1273–1283. |
| Stouffer's method for combining z-scores (unweighted form) | Stouffer, S.A., Suchman, E.A., DeVinney, L.C., Star, S.A. & Williams, R.M. (1949). *The American Soldier, Vol. 1: Adjustment during Army Life.* Princeton University Press. |
| Weighted (`√n`) generalisation of Stouffer's method | Mosteller, F. & Bush, R.R. (1954). "Selected quantitative techniques." In G. Lindzey (Ed.), *Handbook of Social Psychology*, Vol. 1 (pp. 289–334). Addison-Wesley. |
| Rank-based inverse-normal ("normal scores") transformation | Blom, G. (1958). *Statistical Estimates and Transformed Beta-Variables.* Wiley; van der Waerden, B.L. (1952). "Order tests for the two-sample problem and their power." *Indagationes Mathematicae*, 14, 453–458. |
| Shrinkage / partial-pooling toward a pooled estimate (empirical-Bayes family) | Efron, B. & Morris, C. (1975). "Data Analysis Using Stein's Estimator and Its Generalizations." *Journal of the American Statistical Association*, 70(350), 311–319. — *Note*: this project's own `docs/TRANSCRIPTOMICS_SPEC.md` documents a prior instance (an older, now-superseded shrinkage step) where this citation was misapplied to a shrinkage factor that depended only on a fixed count, not on the data's variance structure, and was **not** a genuine empirical-Bayes estimator. The RNA scorer's `n0`-based shrinkage described above is closer to the genuine partial-pooling family in that its weight depends on the actual per-lineage group size, but the constant `n0` itself was chosen by grid search against a validation metric rather than derived analytically from a variance ratio — see the "not independently validated" note below for the precise boundary of this claim. |
| "1-SE rule" for choosing among near-equivalent model-selection candidates | Breiman, L., Friedman, J.H., Olshen, R.A. & Stone, C.J. (1984). *Classification and Regression Trees.* Wadsworth. |
| HPA's own ≥1 TPM detection convention | Uhlén, M. et al. (2015). "Tissue-based map of the human proteome." *Science*, 347(6220), 1260419. (As cited in-line in `transcriptomics.py`.) |

### 1.5 Parameters not independently validated

- **`MIN_PEERS = 5`** (`utils/common.py`): the minimum number of scored cell lines
  required to trust a lineage-level statistic. This is a project-chosen defensive
  floor, not a value derived from a formal power calculation; it is used consistently
  across this file and the protein scorer, but its value itself has not been
  independently validated against an external criterion.
- **`MAD_FLOOR = 0.384`** (`utils/common.py`): an absolute lower bound on the MAD-based
  scale estimate, added to prevent division by a near-zero or degenerate scale. The
  constant is documented in-code as "calibrated against the shared transcriptomics
  notebook" — a project-internal calibration, not a value with external statistical
  justification.
- **`n0`, the lineage-shrinkage constant**: although its *value* was selected by the
  grid-search-plus-1-SE procedure described in §1.2 (a genuine, data-driven internal
  validation), the *grid itself*, the *held-out concordance metric*, and the *choice to
  tie-break toward the conservative end* are project-specific methodological decisions;
  a reviewer should treat this as internally validated rather than externally derived.
- **`upper_cut` / `lower_cut`** (the well-detected / partially-detected / below-detection
  boundaries): calibrated by the same grid-search procedure as `n0`, and subject to the
  same caveat.
- **`floor_q`** (the empirical quantile used to set microarray detection floors):
  likewise calibrated by grid search, not derived from an external detection-limit
  theory for microarray platforms.

---

## 2. Protein Abundance Scoring

**Source:** `final_pipeline/02_Proteinomics/06_protein_score_M2.py` (and its precondition
step, `final_pipeline/02_Proteinomics/platform_tier.py`)

### 2.1 What it does

The protein scorer produces one z-score per gene, per cell line, from two independent
proteomics platforms: ProCAN and CCLE. Neither platform is treated as the reference
platform; both are carried through the pipeline symmetrically (module docstring, "ProCAN
is NOT preferred over CCLE").

**Isoform handling.** Where a gene maps to more than one measured protein identifier
(isoform or paralogous peptide group), the scorer first tests whether all pairs of
those forms correlate with one another across cell lines, using the *minimum* pairwise
Spearman correlation among all isoform pairs (`_min_pairwise_spearman()`, L225–236) —
a conservative, "weakest-link" concordance rule requiring every pair to agree, not just
the best pair. If every pair clears a data-derived concordance threshold, the gene's
per-line value is the **maximum** abundance across its isoforms; if any pair disagrees,
or cannot be assessed for lack of shared cell-line overlap, the value is the **mean**
across isoforms (`collapse_isoforms()`, L239–279) — a symmetric fallback chosen so that
an unrelated or poorly correlated isoform cannot dominate the gene-level value the way
an unconditional maximum would.

**Threshold derivation.** Both the isoform-concordance threshold and the
platform-agreement thresholds described next are derived from the data rather than
fixed constants (`derive_thresholds()`, L415–468). The isoform threshold is set as a
high percentile (default: the 95th) of a null distribution built by correlating
deliberately mismatched pairs of proteins mapped to different genes, within the same
source table (`_sample_null_protein_rho()`, L397–412) — i.e., the correlation level
exceeded by only a small, chosen fraction of unrelated protein pairs is used as the bar
for calling two same-gene isoforms genuinely concordant. The minimum shared-line
overlap required before any correlation is trusted (`seed_overlap`) is itself set from
a classical significance-power calculation: the smallest sample size at which even a
strong correlation (ρ = 0.9) would be statistically significant at the chosen α
(`_power_min_n()`, L386–394, the standard *t*-test for the significance of a Pearson/
Spearman correlation coefficient).

**Platform-agreement assessment.** For each gene, ProCAN and CCLE abundance (after
isoform collapse, not yet z-scored) are correlated across the cell lines both platforms
measured (`gene_agreement()`, L341–358). This per-gene correlation is compared against
two data-derived reference distributions: a **null** distribution built from
deliberately mismatched gene pairs across the two platforms, and a **matched**
distribution built from the real, same-gene correlations (`derive_thresholds()`,
L429–450). The null's high percentile sets a lower agreement bound (`agree_low`,
"disagreement"), and the matched distribution's percentile sets an upper bound
(`agree_high`, "clear agreement"), with the gap between them constrained to a minimum
width. Each gene's correlation is then mapped, by linear interpolation between these two
bounds, to a confidence score in [0, 1] (`confidence_from_rho()`, L361–368) — full
trust at or above `agree_high`, no trust at or below `agree_low`, and a continuous ramp
between them. A gene without enough shared cell lines to compute a correlation at all
defaults to full confidence (`confidence_from_rho()`, L364–365: `if not np.isfinite(rho):
return 1.0`), so the absence of an agreement estimate is not treated as evidence of
disagreement.

**Per-platform scoring.** A robust (median/MAD) z-score is computed independently for
each platform, within each lineage (`zscore_platform()`, `robust_z_matrix()`,
L295–334), gated by the same `MIN_PEERS` minimum-group-size floor used in the RNA
scorer (imported as a local constant, `MIN_PEERS = 5`, L292).

**Combination.** The two platforms' z-scores are then combined by a coverage-weighted
average (`compute_from_sources()` step 4, L521–535), quoted exactly:

```python
# 4. SYMMETRIC coverage-weighted combination (identical treatment of both platforms)
src_n = all_z.groupby(["gene_id", "source"])["model_id"].nunique().rename("src_n").reset_index()
all_z = all_z.merge(src_n, on=["gene_id", "source"], how="left")
all_z["w"] = np.sqrt(all_z["src_n"].astype(float))
```

Here `src_n` is the number of distinct cell lines that gene was measured in on that
platform, i.e., a measure of that platform's breadth of coverage for that gene, and
`w = √src_n` is the standard √n weighting used elsewhere in this pipeline's Stouffer-style
combinations. This produces `z_raw`, the coverage-weighted combined z-score, for every
scored (gene, cell line) pair.

**Confidence application.** Only after this combination is the platform-agreement
confidence applied, and only where it is meaningful to do so — where a cell line
actually has a measurement from *both* platforms for that gene (`n_sources ≥ 2`):

```python
result["z_score"] = np.where(result["n_sources"] >= 2,
                             result["z_raw"] * result["platform_confidence"],
                             result["z_raw"])
```

Where only one platform measured a given cell line, there is no second measurement to
disagree with, so the row passes through unmodified; the confidence score itself is
computed at the gene level, from *other* cell lines where both platforms were present,
and would not correctly characterise this row's evidence if applied to it. This scope
rule was independently re-verified against the code as part of the verification step
preceding this document (§0).

### 2.2 Why this method was chosen

- **Testing isoform concordance rather than silently aggregating** avoids conflating
  biologically or technically distinct measurements under one gene identifier when they
  do not, in fact, track one another; the "weakest-link" minimum-pairwise-correlation
  rule is deliberately conservative, requiring agreement among *all* isoform pairs, not
  merely the best one.
- **Deriving thresholds from permutation-style null distributions**, rather than
  hardcoding a correlation cutoff, ties every threshold in this scorer to the actual
  joint distribution of the data at hand, rather than to a domain convention imported
  from an unrelated dataset. This is a permutation/resampling approach to null
  construction in the general tradition established by Fisher (1935); the specific
  application here — comparing real, matched-feature correlations against correlations
  between deliberately mismatched features to set a significance-informed threshold —
  does not correspond to one single, uniquely citable published method by name, and
  should be described in a dissertation as *this project's application of the general
  permutation-null principle*, not attributed to one specific named test. If a single
  canonical citation is required at viva, this is a genuine gap that would benefit from
  a further citation search rather than an invented attribution.
- **Combining platforms by coverage, and applying agreement as a separate
  confidence multiplier**, rather than either dropping a disagreeing platform outright
  or gating combination on a fixed correlation cutoff, allows every available
  measurement to contribute to the combined estimate while still discounting the
  confidence assigned to genes where the two platforms substantively disagree. This
  design explicitly generalises Stouffer-style combination (weight by amount of
  evidence) and confidence-weighting (discount by agreement) as two independent axes,
  rather than conflating them into a single weight.
- **Scoping the confidence discount to rows with both platforms present** avoids
  penalising a cell line for a disagreement it was never party to: applying a
  gene-level confidence score (derived from lines where both platforms overlapped) to a
  line measured by only one platform would import evidence about *other* cell lines
  into a row that has none of its own to disagree with.

### 2.3 What it replaced, and why the change was made

The current scorer, `06_protein_score_M2.py`, replaced an earlier scorer,
`02_Proteinomics/protein_scorer.py` (retained in the repository; no longer invoked by
`run_all.py`). Its own docstring states three changes over that base version, each
confirmed against the retained file:

1. **Isoforms were not handled explicitly.** The earlier scorer pivoted each platform's
   long-format data to wide form with `aggfunc="first"` (`protein_scorer.py`,
   `_zscore_platform()`, L66): where a gene had multiple mapped UniProt identifiers,
   whichever value happened to appear first in the input was used, with no concordance
   test and no principled aggregation rule.
2. **ProCAN was implicitly preferred over CCLE.** The earlier scorer's combination rule,
   read directly from its own docstring and code, was: for genes whose platform
   agreement fell in the `consistent` or `cautious` tier (fixed cutoffs of ρ ≥ 0.5 and
   ρ ≥ 0.3 respectively, computed by `platform_tier.py`), combine via a Stouffer mean of
   both platforms; for genes falling in the `conflicting` tier (ρ < 0.3), **use ProCAN
   only**, discarding CCLE's measurement entirely. This is an asymmetric treatment: on
   disagreement, one platform's data is deleted rather than both platforms' evidence
   being retained at reduced confidence.
3. **Thresholds were hardcoded domain constants**, not derived from the data:
   `platform_tier.py`'s tier cutoffs (0.5 / 0.3) are fixed values with no permutation-null
   or power-based justification documented in that file.

The replacement therefore directly addresses an asymmetric-preference problem (data
loss when platforms disagree) and an unprincipled-threshold problem (fixed cutoffs with
no data-derived justification), replacing both with the symmetric, data-driven design
described in §2.1. `platform_tier.py` itself was **not removed** — it still runs as
part of Stage 2 of the pipeline and still writes `protein_platform_tier.parquet` using
its original fixed-cutoff logic — but `06_protein_score_M2.py` does not read or use its
output in any of its own calculations; `core_score.py` (§3) checks only that the file
exists as a precondition before running, never its contents (confirmed directly in
`core_score.py`, and independently confirmed in `run_protein_m2_stage.py`'s own
docstring, which documents this exact fact after checking it directly: *"06_protein_score_M2.py
does not read protein_platform_tier.parquet or reference 'tier' anywhere (checked
directly) — it is not an input dependency for M2's own scoring."*).

### 2.4 Statistical basis and citations

| Technique | Citation |
|---|---|
| Median/MAD robust z-score | Rousseeuw & Croux (1993), as in §1.4. |
| √n (coverage) weighting in combination | Mosteller & Bush (1954), as in §1.4. |
| Spearman rank correlation | Spearman, C. (1904). "The proof and measurement of association between two things." *American Journal of Psychology*, 15(1), 72–101. |
| Significance/power calculation for a correlation coefficient (the *t*-test conversion used in `_power_min_n()` / `t_crit()`) | Standard result for the sampling distribution of the Pearson/Spearman correlation coefficient under the null; foundational treatment in Fisher, R.A. (1915). "Frequency distribution of the values of the correlation coefficient in samples from an indefinitely large population." *Biometrika*, 10(4), 507–521. |
| Permutation-based null-distribution construction (general principle) | Fisher, R.A. (1935). *The Design of Experiments.* Oliver & Boyd. — see the caveat in §2.2: the specific null-construction procedure used here (mismatched-feature-pair correlations) is this project's own application of the general principle, not a single citable named method. |

### 2.5 Parameters not independently validated

- **`MIN_PEERS = 5`** (local constant, `06_protein_score_M2.py`, mirroring
  `utils/common.py`): the same project-chosen floor as in §1.5, applied here to the
  protein scorer. Its addition is documented as fixing a concrete, measured problem —
  the module's own comment records that, before this gate, "severely under-observed
  columns produce a real-but-spurious median/scale instead of NaN, driving the extreme
  outlier z-scores flagged earlier (min ≈ −14378, max ≈ 80725)" — but the specific
  threshold value (5) was not itself derived from a formal power calculation.
- **`DEFAULT_ALPHA = 0.05`, `DEFAULT_NULL_PCTILE = 95.0`, `DEFAULT_MATCHED_PCTILE = 60.0`,
  `DEFAULT_MIN_GAP = 0.10`**: these are described in the module's own docstring as
  "statistical CONVENTIONS (not domain-tuned values)" — i.e., standard significance-level
  and percentile conventions, exposed as overridable CLI parameters rather than buried
  constants — but they remain conventional choices, not values derived specifically for
  this dataset.
- **The isoform "weakest-link" rule (requiring all pairwise correlations, not just the
  best, to clear the threshold)**: a deliberate conservative design choice, not one
  validated against an alternative aggregation rule within this project.

---

## 3. Combination into `core_score`

**Source:** `final_pipeline/Scoring/core_score.py`

### 3.1 What it does

`core_score.py` combines the pre-scored RNA (§1) and protein (§2) z-scores into a
single score per (gene, cell line), then ranks cell lines within lineage strata.

**RNA standardization.** RNA z-scores are first standardized within strata defined by
gene and the number of contributing sources, `n_sources` (`_standardize_rna_by_stratum()`,
L83–121). For each `(gene, n_sources)` stratum, the median and standard deviation of
the raw RNA z-score are computed; strata with fewer than `RNA_STRATUM_MIN_N` (5)
contributing lines have their location and scale blended toward the gene-global
location and scale by a shrinkage weight `w = n/(n + RNA_STRATUM_N0)`, the same
`n/(n+n0)` partial-pooling form used in the RNA scorer itself (§1.2). The standardized
value is winsorized at ±`RNA_WINSOR_BOUND` (5.0) as a safety cap.

**Protein residualisation.** Rather than standardizing the protein z-score
independently of RNA, it is first *residualised* against RNA: for each gene, the
Pearson correlation between protein and RNA z-scores across shared cell lines
(`rho_raw`) is blended with a project-wide prior (`RHO_PRIOR = 0.353`, the median
per-gene correlation across the panel) by an empirical-Bayes-style weight,
`lam = K_RHO / (n + K_RHO)`, decreasing as the number of shared lines `n` grows
(L159–176). This shrunk correlation, `rho_g`, is converted from a correlation to a
regression (OLS) slope, `beta_g = rho_g × SD(protein)/SD(RNA)` — the code is explicit
that using the correlation directly, rather than the slope, would overcorrect, since
the two arms' measured standard deviations differ substantially (SD_RNA ≈ 1.525,
SD_protein ≈ 0.951). **As of 2026-08-29 (finding "C3", promoted), the two per-gene SDs
feeding this ratio are themselves shrunk toward their panel-wide values** with the
same `w = n/(n + SD_SHRINK_N0)` partial-pooling form used throughout this file
(`SD_SHRINK_N0 = 25`) — previously only `rho_g` was protected against small-sample
noise this way, while the two SDs multiplied into `beta_g` were raw per-gene
estimates; see `docs/WHY_REFERENCE.md` §3.2 for the full verification record
(near `N_MIN_RHO=30`, `beta_g` shifted by a median 19.9% under this fix, with the
full-pipeline downstream spike improving slightly, not worsening). The protein
residual, `prot_z − beta_g × rna_z`, is then
standardized per gene: `(x − median) / SD` (not MAD — the code notes this is a
deliberate trade-off of robustness for the unit variance the subsequent Stouffer-style
combination requires).

**Orthogonal combination.** The standardized RNA and residualised-and-standardized
protein arms are combined as:

```
core_z = (W_RNA · rna_std_z + W_PROT · prot_resid_z) / √(W_RNA² + W_PROT²)
```

with `W_RNA = √1.431` and `W_PROT = √1.45`, both derived from Kish's effective sample
size formula, `n_eff = N / (1 + (N−1)ρ̄)`, applied to the RNA arm's three correlated
sources (`N = 3`, mean pairwise correlation ρ̄ = 0.548, measured) and the protein arm's
two platforms (`N = 2`, ProCAN–CCLE ρ = 0.373). `core_z` is then mapped through the
standard normal CDF, `Φ(core_z)`, to produce `core_score ∈ [0, 1]`, a monotone
relabelling that preserves rank order while giving the score a probability-like scale.

Cell lines with expression data but no protein data for a given gene are not dropped:
they receive an `n_layers = 1` fallback score (L225–242), and genes present in RNA but
absent from either proteomics platform receive the same RNA-only treatment across all
lines (L244–256). **As of 2026-08-29 (finding "C7", promoted), this fallback score is
`Φ((W_RNA/W_NORM) · rna_std_z)`, not `Φ(rna_std_z)` alone** — prior to this fix, the
fallback path omitted the same `W_RNA/W_NORM ≈ 0.7048` downweight that `n_layers = 2`
rows apply before the CDF, so identical RNA evidence scored systematically more
extreme via the fallback path than via the combined path with a neutral protein
residual, while both row types were ranked together within the same `(gene, lineage)`
stratum. See `docs/WHY_REFERENCE.md` §3.6 for the full verification record (full-panel
spike 2.8444%→1.4790%; `n_layers=1`-only spike 2.7252%→0.4894%). Finally, within each `(gene, lineage)`
stratum, cell lines are ranked by `core_score` in descending order, with ties broken
deterministically by `model_id` — the code documents that, before this explicit
tie-break, unstable Python dictionary/set iteration order made 37.1% of rows' rank
differ between two runs on identical data (L258–268).

### 3.2 Why this method was chosen

- **Standardizing RNA per `(gene, n_sources)` stratum**, rather than per gene globally,
  addresses a specific, measured problem: the standard deviation of the raw RNA z-score
  is genuinely heterogeneous by how many sources contributed to it (measured on the
  live output: SD ≈ 0.88 / 1.54 / 1.86 for `n_sources = 1 / 2 / 3` respectively), and a
  naïve per-gene-global standardization does not remove this heterogeneity — the code
  documents that, after a naïve global standardization, per-`n_sources` SD remained
  heterogeneous (0.65 / 1.08 / 1.35), because the effect is a *stratum-level* scale
  difference, not a *gene-level* one. Standardizing within each `(gene, n_sources)`
  stratum directly targets this. The correction was validated on a worked example
  (`ENSG00000178828`), confirmed to correct a genuine rank misordering in that gene's
  scores rather than merely relabelling the same distortion under a different name —
  a concrete, gene-level empirical check, not only a summary-statistic argument.
- **Residualising protein against RNA before combining them** is chosen because RNA and
  protein abundance for the same gene are correlated, non-independent measurements of
  related biology; combining them as if independent would double-count shared signal.
  Regressing protein on RNA and combining the *residual* isolates the information
  protein contributes beyond what RNA already explains, which is the orthogonal
  combination the section header names.
- **Shrinking the per-gene correlation estimate toward a panel-wide prior** before
  computing the regression slope stabilises the slope estimate for genes with few
  shared RNA/protein observations, following the same logic as the RNA scorer's
  lineage-shrinkage (§1.2): more data (larger `n`) trusts the gene's own estimate more;
  less data leans more on the panel-wide value.
- **The Kish-derived weights (`W_RNA`, `W_PROT`)** replace a naïve assumption that each
  arm's constituent sources are fully independent, which would overstate each arm's
  effective evidence. This correction was independently re-derived and applied within
  this project: `docs/ALTERATION_DEFENCE.md`'s change log records `W_RNA` being
  corrected from a provisional `√2.00` to the measured `√1.431` (Kish n_eff at
  ρ̄ = 0.548), with the corrected pipeline re-run end-to-end and a measured effect on
  held-out evaluation of +0.0238 → +0.0240 hit@20 (a +0.7% change) — a concrete,
  measured empirical consequence of this specific correction, not merely a theoretical
  justification. The project's own preregistration documents (`EXPR_FUSION_PREREGISTRATION.md`)
  further flag that Kish's design-effect formula is a cluster-sampling result being
  applied here *by analogy* to correlated evidence sources, and specify a fallback test
  (comparing the Kish-scalar weights against an explicit covariance-based combination,
  e.g. via the `poolr` package's generalised Stouffer methods) to be run if the two
  disagree by more than 10%.
- **The `Φ(core_z)` transform to `[0, 1]`** is the standard probability integral
  transform of a (approximately) standard-normal combined score; it changes the scale
  on which the score is read without altering the rank order it induces.
- **The RNA-only fallback for protein-free lines/genes** was added because, without it,
  48.6% of RNA-scored pairs for the shared-gene set would receive no score at all
  purely because protein data was unavailable for that line — a coverage decision
  documented directly in the code's own comment at the point the fallback is
  constructed.

### 3.3 What it replaced, and why the change was made

`core_score.py`'s own docstring documents that the RNA-standardization step described
in §3.1 "replaces the prior within-lineage percentile + clipped `norm.ppf` treatment,
which double-conditioned RNA (already lineage-conditioned upstream in `rna_scorer.py`)
and left it structurally suppressed relative to protein's uncapped per-gene treatment."
In other words, the earlier version of this combination step took RNA's already
lineage-conditioned z-score (itself the output of a lineage-level scoring step,
§1.3's superseded `rna_scorer.py`) and applied a *second* lineage-conditioning
transform on top of it — ranking cell lines into a within-lineage percentile and then
inverse-normal-transforming that percentile back onto a z-like scale, with clipping.
This second transform compounds the effect of the first, compressing RNA's usable
dynamic range relative to protein, which received only a single, uncapped
`(x − median)/SD` treatment. The corrected version standardizes RNA exactly once, in a
stratum chosen to remove genuine `n_sources`-driven scale heterogeneity rather than to
re-apply lineage conditioning a second time, putting RNA and protein through
comparable single-step standardizations before they are combined. `W_RNA` and `W_PROT`
were deliberately left unchanged by this fix; the code comment explicitly flags whether
they still make sense given the corrected RNA treatment as "an explicitly separate,
out-of-scope question for a later task" — this is not resolved by the present
methodology and should be described as an open question, not a settled one, in any
write-up.

### 3.4 Statistical basis and citations

| Technique | Citation |
|---|---|
| Median/MAD-family robust standardisation, and its use here for a location–scale pair | Rousseeuw & Croux (1993), as in §1.4. |
| Kish's effective sample size (design effect) for correlated evidence | Kish, L. (1965). *Survey Sampling.* John Wiley & Sons. — applied here by analogy to correlated measurement sources rather than to its original cluster-sampling context; see the caveat above. |
| Empirical-Bayes-style shrinkage of a per-gene correlation toward a panel-wide prior | Efron & Morris (1975), as in §1.4, with the same caveat: the shrinkage weight here is data-adaptive (a function of the actual number of shared observations `n`), and the prior itself is empirically estimated (the panel-wide median correlation), which distinguishes it from the misapplied instance documented in `docs/TRANSCRIPTOMICS_SPEC.md`; however, the pseudo-count `K_RHO` that sets the *strength* of shrinkage is a fixed, hand-set value, not one derived from a formal variance-ratio calculation — see §3.5. |
| Generalised (correlation-corrected) Stouffer combination, as a documented alternative check on the Kish-scalar weighting | Strube, M.J. (1985). "Combining and comparing significance levels from nonindependent hypothesis tests." *Psychological Bulletin*, 97(2), 334–341; Hartung, J. (1999). "A Note on Combining Dependent Tests of Significance." *Biometrical Journal*, 41(7), 849–855. |
| Standard normal cumulative distribution function as a probability integral transform | Standard result; no specific attribution required. |

### 3.5 Parameters not independently validated

- **`RNA_STRATUM_N0 = 25`**: the code's own comment states this value "mirrors
  [`transcriptomics.py`'s `n0`] as a FIRST-PASS choice only — not independently
  validated for this context." This is the clearest, most explicit example in the
  codebase of a parameter the project itself flags as unvalidated, and should be
  described as such without softening.
- **`K_RHO = 30`** (the empirical-Bayes pseudo-count for the RNA–protein correlation
  blend): described in-code only as "EB pseudo-count: prior/data get equal weight at
  n = K_RHO," with no documented calibration procedure comparable to
  `transcriptomics.py`'s grid search for its own `n0`. This should be treated as a
  hand-set, unvalidated value.
- **`SD_SHRINK_N0 = 25`** (added 2026-08-29, finding "C3"): the shrinkage weight
  applied to `sd_rna`/`sd_prot` before they enter the `beta_g` ratio reuses the same
  `n0=25` convention as `RNA_STRATUM_N0`/`PROT_RESID_N0` by deliberate choice, not
  independent derivation for this specific application (shrinking an SD that feeds a
  ratio inside a slope calculation, not a direct location/scale standardization) — the
  fourth confirmed instance of this project's "one half of a paired estimate gets
  shrinkage, the other doesn't" pattern, now fixed. See `WHY_REFERENCE.md` §3.2.
- **`RNA_WINSOR_BOUND = 5.0`**: a safety cap described as replacing "the old
  percentile step's implicit ±3.09 bound," observed in the live data to bind on fewer
  than 0.2% of rows post-standardization — an empirical observation about its effect,
  but not a validated choice of the specific bound.
- **`W_RNA`, `W_PROT`**: unlike the shrinkage constants above, these *were* derived from
  a measured quantity (the panel's own mean pairwise source correlation) rather than
  hand-set, and the correction to them was empirically evaluated end-to-end (§3.2).
  They should therefore be described as internally validated, with the important
  caveat, stated directly in the code, that their continued appropriateness under the
  corrected RNA-standardization treatment is an open, unresolved question.

### 3.6 Promoted refinements since §3.1–3.5 (2026-08-27 / 2026-08-28)

Two further, independently-verified fixes have been promoted to live `core_score.py`
on top of the RNA-standardization fix described above. Both were investigated,
tested against full-scale live data, and verified structurally (bit-identity on
every unaffected row) before promotion — the same discipline as §3.3's original fix.

**Protein-residual shrinkage + clip (2026-08-27).** The RNA-standardization fix
added shrinkage (`w = n/(n+n0)` toward a gene-global scale) and a winsorization
bound to the RNA arm, but left the protein-residual standardization step
untouched — confirmed, by diffing against the pre-fix backup
(`core_score_v1_pre_rna_fix_*.py.bak`), to be an inherited gap rather than a
considered asymmetry. `_standardize_prot_resid_by_gene()` closes it: each gene's
own residual SD is shrunk toward the panel-wide residual SD with the same
`n/(n+PROT_RESID_N0)` weight (`PROT_RESID_N0 = 25`, reused from `RNA_STRATUM_N0`
by convention, not independently re-derived), applied continuously rather than
RNA's hard `n < 5` cutoff, and the result winsorized at ±`PROT_RESID_WINSOR_BOUND`
(5.0). Measured effect: `core_score` spike (fraction with `core_score < 0.01` or
`> 0.99`) fell from 3.098% → 3.039% on `n_layers=2` rows and 2.870% → 2.847%
full-panel. A companion investigation found that shrinking the *earlier*
per-(gene, lineage) protein z-scoring stage instead
(`02_Proteinomics/06_protein_score_M2.py`'s `robust_z_matrix()`) substantially
improved the raw protein z-scores (|z|≥5 fraction: −68% full-panel) but, tested
end to end through this combination step, made `core_score`'s own spike
*slightly worse* — it interacts with the per-gene `rho_g`/`beta_g` regression
this file fits fresh from whatever protein input it receives. That fix was
evaluated and explicitly **not** promoted.

**Regime-3 ("genuine disagreement") substitution (2026-08-28).** An unsupervised-
free regime-discovery pass — simple, interpretable cross-tabulation of
`|rna_std_z|`/`|prot_std_z|` against a "signaled" threshold (1.0), rather than a
black-box clustering method — identified five candidate regimes in the
`n_layers=2` population. One of them, where both arms are confidently signaled
(`|z| > 1.0`) but point in **opposite** directions (~1.1% of `n_layers=2` rows,
75,117 of 6,762,902), showed Variant A's residualization producing scores
*more* extreme (further from 0.5) than a plain, no-residualization average of
the two arms' independently-standardized z-scores — the only regime where this
ordering reverses — and the lowest per-gene `rho_g` (empirical RNA–protein
correlation) of any regime, meaning the regression A relies on is least
trustworthy exactly where the two arms disagree. For rows matching this regime
only, `core_z` is computed as the direct average
`(W_RNA·rna_std_z + W_PROT·prot_std_z)/norm` instead of A's residualized default.
This substitution is mechanistically distinct from an earlier, broader attempt
to apply the same no-residualization formula pipeline-wide (which measurably
hurt the aggregate spike, because it let RNA "run unchecked" wherever protein
was uninformative) — regime 3 requires protein to be confidently signaled by
definition, so that failure mode cannot occur here. Verified before promotion:
bit-identical to the pre-fix scorer on all non-matching rows (max diff = 0.0);
population reproduced exactly (75,117 rows); spike improved 3.0389% → 3.0313%
(`n_layers=2`), 2.8473% → 2.8444% (full panel) — small, since the regime is
~1.1% of rows, but real and directionally consistent with the mechanism.
Substituted rows are marked in a `combination_variant` audit column
(`"A"` default, `"C_regime3"` where substituted) rather than silently changed.

**Regime 1 ("protein-neutral") — tested, not resolved.** The same
regime-discovery pass identified a second, much larger regime (~21.1% of
`n_layers=2` rows, ~1.43M): protein reads near-neutral while RNA is strongly
signaled — the pattern responsible for known cases like *KRT86*/ACH-000878,
where A's residualization manufactures a spuriously suppressed score
(0.234) despite an unambiguous RNA signal (`rna_std_z ≈ 4.0`). Three
independent substitution strategies were tested against this regime and
**none improved the aggregate spike**:
1. *Precision-targeted C-substitution* (8,500 rows matching
   `|prot_std_z| < 0.15 ∧ |prot_resid_z| > 1.5`, the exact rule underlying
   the *KRT86* pattern): washed — no net aggregate improvement, because the
   direct-average formula is not uniformly "safer" for this population; it
   simply defers to RNA unchecked, which is not correct for every row it
   touches.
2. *Whole-bucket RNA-only substitution* (`core_score = Φ(rna_std_z)`, i.e.
   treating "protein present but neutral" identically to true `n_layers=1`):
   made the aggregate spike substantially **worse** (n_layers=2: 3.0389% →
   4.6204%, +1.58 pp). Diagnosed mechanistically: 8.02% of this regime already
   has `|rna_std_z| > 2.33` (Φ's saturation threshold) on its own, and A's
   weighted combination (`W_RNA/W_NORM ≈ 0.705`) was already correctly
   damping the great majority of these rows before the CDF — most of the
   regime was not actually broken, and removing that damping wholesale
   exposed those undamaged rows to the same saturation the fix was meant to
   prevent for the genuinely broken minority.
3. *Precision-targeted RNA-only substitution* (the same 8,500-row rule as
   attempt 1, but with attempt 2's RNA-only formula instead of the direct
   average): fixes *KRT86* itself exactly as intended (0.234 → 1.000), but
   still made the aggregate spike worse (n_layers=2: 3.0313% → 3.1149%,
   +0.084 pp, measured against the regime-3-promoted baseline). The same
   8,500-row rule that flags a genuine suppression case is, by construction,
   enriched for rows where `rna_std_z` is already large (since `prot_resid_z`
   is large precisely when RNA is large and protein is neutral), so even this
   narrow substitution pushes most of the targeted rows toward Φ's saturated
   tails.

With three independently-motivated substitution strategies tested at two
different targeting granularities (whole-regime and precision-targeted) and
two different formulas (direct average, RNA-only), none of which improves the
aggregate spike, regime 1 is documented here as a genuine open limitation of
the current scoring design — individual cases like *KRT86* can be shown to be
mis-scored by Variant A, but no simple variant substitution fixes the regime
as a whole without a larger aggregate cost elsewhere in the same population.
Resolving it, if pursued, likely requires new design work (e.g. a formula that
distinguishes *why* protein reads neutral — genuinely uninformative vs.
genuinely absent-of-effect — rather than a threshold-tuned substitution
between two existing variants) rather than further threshold tuning.

---

## 4. Driver-Gated Routing and Confidence Tiers

**Sources:** `final_pipeline/Scoring/driver_routing.py`, `final_pipeline/Scoring/confidence_tiers.py`

These two stages were not the focus of this investigation and are described here more
briefly, for completeness.

### 4.1 What they do

**Driver-gated routing** (`driver_routing.py`) merges `core_score` with three
independent alteration-evidence channels computed upstream — mutation impact
(`p_mutation`, from a Noisy-OR combination of VEP consequence rank, AlphaMissense/CADD
pathogenicity, and variant burden, each passed through a Hill-equation-shaped monotone
transform), fusion impact (`p_fusion`, a similarly-constructed Noisy-OR combination),
and copy-number alteration (role-gated amplification/deletion flags derived from
COSMIC's own ploidy-aware calls). A gene/cell-line pair is flagged
`has_driver_alteration` if any one channel clears a fixed threshold
(`MUT_DRIVER_THRESHOLD = FUS_DRIVER_THRESHOLD = 0.5`), and `is_lof_alteration` if the
alteration is a loss-of-function type (a truncating mutation or a deletion) in a
tumour-suppressor-role gene. A recovery step adds back alteration-driver pairs that
have no expression score at all (because the gene was silent in that lineage, or the
line was never profiled) with `core_score = NaN`, rather than silently dropping them —
the code documents that, without this recovery, roughly 32% of mutation drivers, 37%
of fusion drivers, and 56% of CNA alterations would be lost.

**Confidence tiers** (`confidence_tiers.py`) then assigns each pair to one of four
tiers from the cross of a within-gene score percentile (top 20% or not) and the driver
flag: `HIGH` (top-scoring and driver-supported), `MEDIUM` (top-scoring, no driver
support), `CONTEXT` (driver-supported, not top-scoring), `LOW` (neither). For
loss-of-function alterations in oncogene/TSG-role genes, the "top-scoring" criterion is
inverted — such alterations are expected to associate with *low*, not high, expression,
so the directional check ranks by the bottom 20% instead. Genes with `gene_role == "both"`
are treated as gain-of-function for this directional check regardless of alteration
type, because the code notes there is no empirical basis in this project for a
directional split within that stratum. Three partition invariants are asserted in code
at the end of this stage (HIGH+CONTEXT equals exactly the driver-flagged set;
MEDIUM+LOW equals exactly the non-driver set; all four tiers sum to the total row
count), so any future change that breaks the partition logic fails loudly rather than
silently mis-tiering rows.

### 4.2 Why chosen, and unvalidated parameters

The Noisy-OR combination rule for independent binary-ish evidence channels is a
standard construction from probabilistic graphical models (Pearl, 1988); the Hill
equation used to squash each individual channel onto [0, 1] is borrowed, by
functional-form analogy rather than mechanistic claim, from the dose–response
literature (Hill, 1910) — its use here is a monotone transform with a tunable
half-maximal point and steepness, not an assertion that mutation burden or fusion
recurrence follows ligand-binding kinetics. The specific Hill parameters
(`p0`, `k`, per channel) and the 0.5 driver thresholds, the 80th-percentile cutoff for
"top-scoring," and the fixed 0.15 in-frame boost for fusions are all project-specific,
hand-set calibrations, documented in-code with the qualitative behaviour they were
tuned to produce (e.g., "rank 3 (HIGH/truncating) = 0.590 > gate"), not derived from an
external statistical criterion. These should be described in a dissertation as
domain-informed engineering choices, not as independently validated statistical
thresholds.

### 4.3 Citations

| Technique | Citation |
|---|---|
| Noisy-OR combination of independent evidence channels | Pearl, J. (1988). *Probabilistic Reasoning in Intelligent Systems: Networks of Plausible Inference.* Morgan Kaufmann. |
| Hill equation (functional form borrowed by analogy) | Hill, A.V. (1910). "The possible effects of the aggregation of the molecules of hæmoglobin on its dissociation curves." *Journal of Physiology*, 40 (Suppl.), iv–vii. |
| AlphaMissense pathogenicity scores | Cheng, J. et al. (2023). "Accurate proteome-wide missense variant effect prediction with AlphaMissense." *Science*, 381(6664), eadg7492. |
| CADD pathogenicity scores | Rentzsch, P., Witten, D., Cooper, G.M., Shendure, J. & Kircher, M. (2019). "CADD: predicting the deleteriousness of variants throughout the human genome." *Nucleic Acids Research*, 47(D1), D886–D894. |
| Per-drug MAD-z sensitivity labelling, used downstream in held-out evaluation of these tiers | Iorio, F. et al. (2016). "A Landscape of Pharmacogenomic Interactions in Cancer." *Cell*, 166(3), 740–754. |

---

## 5. Summary: Validated vs. Unvalidated Parameters

| Parameter | Stage | Status |
|---|---|---|
| `n0`, `upper_cut`, `lower_cut`, `floor_q` | RNA scorer | Internally validated (grid search + held-out cross-source concordance, 1-SE rule) |
| `MIN_PEERS = 5` | RNA scorer, protein scorer | Unvalidated defensive floor |
| `MAD_FLOOR = 0.384` | shared utility | Unvalidated, project-internal calibration |
| Isoform "weakest-link" concordance rule | Protein scorer | Design choice, not validated against alternatives |
| `DEFAULT_ALPHA`, `DEFAULT_NULL_PCTILE`, `DEFAULT_MATCHED_PCTILE`, `DEFAULT_MIN_GAP` | Protein scorer | Standard statistical conventions, not dataset-specific derivations |
| `W_RNA`, `W_PROT` (Kish n_eff weights) | `core_score` | Internally validated (measured ρ̄, end-to-end evaluated); their fit to the *corrected* RNA treatment is an explicitly open question |
| `RNA_STRATUM_N0 = 25` | `core_score` | Explicitly flagged in-code as unvalidated for this context |
| `K_RHO = 30` | `core_score` | Unvalidated, no documented calibration |
| `SD_SHRINK_N0 = 25` (C3, promoted 2026-08-29) | `core_score` | Reused `n0=25` convention, unvalidated for this specific application; downstream effect measured flat-to-slightly-improved |
| `RNA_WINSOR_BOUND = 5.0` | `core_score` | Empirically observed to rarely bind; the bound itself not independently validated |
| Driver thresholds (0.5), top-score percentile (0.80), Hill parameters, in-frame boost (0.15) | Driver routing / confidence tiers | Hand-set, domain-informed engineering choices |
| `SENSITIVITY_MAD_Z = −1.0` | Held-out evaluation | Project-chosen threshold, checked post hoc against a 5–40%-sensitive-per-drug guard rather than derived a priori |

---

## 6. Full Reference List

- Blom, G. (1958). *Statistical Estimates and Transformed Beta-Variables.* Almqvist & Wiksell / John Wiley & Sons.
- Breiman, L., Friedman, J.H., Olshen, R.A. & Stone, C.J. (1984). *Classification and Regression Trees.* Wadsworth.
- Cheng, J. et al. (2023). "Accurate proteome-wide missense variant effect prediction with AlphaMissense." *Science*, 381(6664), eadg7492.
- Efron, B. & Morris, C. (1975). "Data Analysis Using Stein's Estimator and Its Generalizations." *Journal of the American Statistical Association*, 70(350), 311–319.
- Fisher, R.A. (1915). "Frequency distribution of the values of the correlation coefficient in samples from an indefinitely large population." *Biometrika*, 10(4), 507–521.
- Fisher, R.A. (1935). *The Design of Experiments.* Oliver & Boyd.
- Hartung, J. (1999). "A Note on Combining Dependent Tests of Significance." *Biometrical Journal*, 41(7), 849–855.
- Hill, A.V. (1910). "The possible effects of the aggregation of the molecules of hæmoglobin on its dissociation curves." *Journal of Physiology*, 40 (Suppl.), iv–vii.
- Iorio, F. et al. (2016). "A Landscape of Pharmacogenomic Interactions in Cancer." *Cell*, 166(3), 740–754. PMID 27397505.
- Kish, L. (1965). *Survey Sampling.* John Wiley & Sons.
- Mosteller, F. & Bush, R.R. (1954). "Selected quantitative techniques." In G. Lindzey (Ed.), *Handbook of Social Psychology*, Vol. 1 (pp. 289–334). Addison-Wesley.
- Pearl, J. (1988). *Probabilistic Reasoning in Intelligent Systems: Networks of Plausible Inference.* Morgan Kaufmann.
- Rentzsch, P., Witten, D., Cooper, G.M., Shendure, J. & Kircher, M. (2019). "CADD: predicting the deleteriousness of variants throughout the human genome." *Nucleic Acids Research*, 47(D1), D886–D894.
- Rousseeuw, P.J. & Croux, C. (1993). "Alternatives to the Median Absolute Deviation." *Journal of the American Statistical Association*, 88(424), 1273–1283.
- Spearman, C. (1904). "The proof and measurement of association between two things." *American Journal of Psychology*, 15(1), 72–101.
- Stouffer, S.A., Suchman, E.A., DeVinney, L.C., Star, S.A. & Williams, R.M. (1949). *The American Soldier, Vol. 1: Adjustment during Army Life.* Princeton University Press.
- Strube, M.J. (1985). "Combining and comparing significance levels from nonindependent hypothesis tests." *Psychological Bulletin*, 97(2), 334–341.
- Uhlén, M. et al. (2015). "Tissue-based map of the human proteome." *Science*, 347(6220), 1260419.
- van der Waerden, B.L. (1952). "Order tests for the two-sample problem and their power." *Indagationes Mathematicae*, 14, 453–458.

**Not independently confirmed / needs a further citation search if required for
submission:** the specific permutation-null-threshold construction used in
`06_protein_score_M2.py`'s `derive_thresholds()` (§2.2) does not correspond to one
single, uniquely named published method beyond the general permutation-testing
tradition (Fisher, 1935); a targeted literature search for "empirical null" threshold
methods (e.g., Efron's empirical-null work in large-scale inference) may yield a more
precise citation, but none is asserted here that was not directly confirmed.
