# MATH_TRIAL_AND_ERROR — retired mathematics

Formulas that were **superseded**, **tested and rejected**, or **built on a
claim since retracted**. Moved out of
[MATH_REFERENCE.md](MATH_REFERENCE.md) so that document describes only what
the live pipeline computes.

Nothing here is deleted. In a project whose main findings are negative, the
derivations that produced them are the evidence.

| section | status |
|---|---|
| 2.3 Noisy-OR combination (superseded | superseded |
| 2.5 Alternative normalisation / combination variants (validation harness only) | tested, rejected |
| 2.7 Chronos three-layer blend (tested, rejected) | tested, rejected |
| C.7 Chronos correlation gate | tested, rejected |
| C.10 Constant-term global-scalar blend (rejected candidate layers) | tested, rejected |
| C.13 Metabolomics enzyme-flux proxy (tested, null) | tested, rejected |
| C.14 miRNA-as-inhibitor dampening (tested, mixed/weak) | tested, rejected |
| C.15 Partial AUROC and the gate/ranker distinction | RETRACTED |
| C.16 Detection call vs percentile | tested, rejected |
| C.17 Paralog stratification | tested, rejected |
| C.18 Stratum centring and the two-layer score (protein is a net negative) | tested, rejected |

---

## 2.3 Noisy-OR combination (superseded — retained because its failure motivates §2.4)

> **Why retired.** Superseded by §2.4. Retained because its failure mode (non-idempotent t-conorm inflating correlated evidence) is the derivation that motivated the production formula.


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


---

## 2.5 Alternative normalisation / combination variants (validation harness only)

> **Why retired.** Validation-harness variants only. Never in production; kept for reproducibility of the harness.


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


---

## 2.7 Chronos three-layer blend (tested, rejected)

> **Why retired.** Tested and rejected. The label it was blended against was later shown to be Project Score, not Chronos (`project_chronos_label_provenance`).


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


---

## C.7 Chronos correlation gate

> **Why retired.** Built on the mis-identified Chronos label; the correlation gate it defines is not used by any live path.


```
use_chronos = ( ρ(chronos_pct, ground_truth) > 0.05 )      # per gene
```

A per-gene safety gate: only add Chronos where it aligns with drug sensitivity.
Note this is an *unthresholded-in-p* gate — it uses effect size only, so at
typical `n` a substantial fraction of genes pass by chance (`P(ρ > 0.05) ≈ 0.07`
at `n = 900` under H₀ using `√(n−1)ρ ~ N(0,1)`). It is a filter, not a test, and
should be described as such.

*Provenance:* `validate_chronos.chronos_correlation_gate():121-128`.


---

## C.10 Constant-term global-scalar blend (rejected candidate layers)

> **Why retired.** Rejected candidate layers. Never adopted.


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


---

## C.13 Metabolomics enzyme-flux proxy (tested, null)

> **Why retired.** Tested, null. Metabolomics carries no gene axis; it survives only as a cell-line similarity axis, and `cell_similarity/05` shows it does not beat a tissue lookup.


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


---

## C.14 miRNA-as-inhibitor dampening (tested, mixed/weak)

> **Why retired.** Tested, mixed/weak. Same outcome as C.13 — miRNA is a similarity axis, not a gene-level layer.


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


---

## C.15 Partial AUROC and the gate/ranker distinction

> **Why retired.** **RETRACTED.** This is the mathematics of the abundance gate. `gate_audit/03` showed 98% of the published decile-1 depletion is reproduced by a placebo with randomised labels — it was an artefact of scoring never-screened lines as confirmed non-dependencies. The partial-AUROC machinery itself is sound and is still used in `gate_audit/`; the *claim* it supported is withdrawn.


`test_run_abundance_gate.py`, `test_run_gate_mechanism.py`,
`test_run_screen_replication.py`

**The measurement error this section corrects.** Earlier evaluation reported
AUROC ≈ 0.52–0.56 for every scoring rule tried and read it as a weak result.
AUROC integrates over all thresholds and concentrates its mass where the *top* of
the ranking is, but `core_score` is an abundance prior — a necessary-not-
sufficient condition. The claim it licenses is not "high-abundance lines are
enriched for dependency" but "low-abundance lines are depleted of it". A filter
measured with a ranking metric scores near chance while working perfectly.

**McClish-standardised partial AUROC.** Over an FPR interval `[lo, hi]`,

```
A     = ∫_lo^hi TPR(f) df
A_min = (hi² − lo²)/2          (the chance line over the same interval)
A_max = hi − lo
pAUC  = ½ · (1 + (A − A_min)/(A_max − A_min))          so 0.5 = chance
```

The gate region is `FPR ∈ [0.8, 1.0]` — a rule excluding the bottom decile sits
near FPR 0.9. The ranker region `[0.0, 0.2]` is reported alongside because the
contrast between them *is* the finding.

**Result.** Gate region above chance, ranker region at chance, across three
independent label sources:

| label | gate pAUC | ranker pAUC | decile-1 ratio | decile-10 ratio |
|---|---|---|---|---|
| GDSC | 0.601 | 0.518 | 0.807 | 1.243 |
| Project Score | 0.593 | 0.501 | 0.831 | 1.001 |
| Chronos | 0.596 | 0.500 | 0.839 | 0.986 |

NPV is the metric that moves: `0.905 → 0.931` at 5% exclusion (Chronos),
`0.865 → 0.897` (GDSC).

**Null for the decile-1 count is exact, not permuted.** Drawing `n₁` of a gene's
`N` lines without replacement holds prevalence and line count fixed, so the null
count is hypergeometric and closed-form — the same thing a within-gene label
permutation approximates, but exact and scalable to thousands of genes.

**Replication (`test_run_screen_replication.py`).** Decile-1 depletion holds on
two independent CRISPR screens: Broad 0.9950, Sanger 0.9960, difference
`p = 0.19` (n = 626 genes; 176 shared lines, 16,731 shared genes). On selective
genes (n = 66) both hold and both are significant.

*The number that sets the ceiling for every ranking claim in this project:* the
two ground truths agree with **each other** at median Spearman `ρ = 0.172`, with
only 1.9% of genes above 0.5 and 5.3% negative. No score can exceed the
reliability of the label it is measured against, so a ranker AUROC near 0.5 is
partly a statement about Chronos and GDSC, not only about `core_score`.

**The adversarial check, and the caveat it produced
(`test_run_gate_mechanism.py`).** Two controls were run against the gate result:

- *Label noise floor.* Control gene families expected to be non-essential
  (TAS2R\*, 24 genes, 21,744 pairs) show a dependency rate of 2.56%. The measured
  floor rate is 1.39% — *at or below* the label's own noise floor, so the
  depletion is complete to within what the label can resolve, and quoting it as
  a percentage understates it. (KRTAP\* and VN1R\* were excluded from the control
  set as paralogy-contaminated: multi-targeting guides in near-identical families
  inflate apparent dependency.)
- *Profiling confound.* Stratifying on whether the line was profiled at all:
  floor stratum detected 2.45% vs not-detected 2.29%, ratio 1.07,
  **Fisher p = 0.246**; mid stratum 23.4% vs 14.5%, ratio 1.62, `p = 0`.

Verdict: `no_floor_specific_effect__profiling_confound_only`. **At the very
bottom of the score, depletion is not separable from "this gene was never
measured in this line."** This caveat must be quoted with the gate result, not
after it.

**What was ruled out (`test_run_improvement_scan.py`).** Lineage-conditioned
percentiles are neutral (Δ = 0.0, `p = 0.54`, better in 49.0% of genes — a coin
flip). Conformal exclusion is valid but adds nothing over random
(`ratio_vs_random` 0.93–0.99, calibrated safe fraction 0.0).

---


---

## C.16 Detection call vs percentile — is the layer an encoded boolean?

> **Why retired.** Supporting analysis for the retracted gate. The detection-vs-percentile contrast was the gate's defence against the tautology objection; with the gate withdrawn it no longer supports a live claim.


`test_run_detection_vs_percentile.py`

**Question.** If all of the measured signal is at the bottom (§C.15), the
continuous percentile may be encoding nothing more than "is this gene on or off",
in which case an absolute detection threshold would replace it — and would be
panel-independent, which matters because the percentile's denominator is whichever
cell lines are in the build.

**Why the comparison must be against a per-sample axis.** By the monotone-
invariance corollary in §2.1, percentile and raw expression are identical under
pAUC. Only a score normalised on the *other* axis reorders cell lines within a
gene. So the contrast is against a per-line detection distance.

**Calibration by curated unexpressed reference.** The unexpressed population is
seeded from genes not expressed in cultured lines — OR\*, TAS2R\*, VN1R\* (438
genes with an expression column), the same families §C.15 uses as controls.
Verified per run rather than assumed: reference median `log2TPM = −4.644` against
a transcriptome median of `+3.046`, a separation of **7.690 log₂ units**. Per-line
detection cuts are then quantiles of that reference, so the threshold parameter is
a false-positive rate on a real negative control set:

```
expressed(g, ℓ)  ⟺  log2TPM(g, ℓ) > Q_q( reference | ℓ )
q ∈ {0.90, 0.95, 0.99, 0.999}      plus a TPM ≥ 1 rule and a Hart variant
```

**A pseudocount trap, recorded because it silently produced a wrong answer
first.** `depmap_expr_clean.parquet` is `log2(TPM+1)`, and the `+1` collapses the
unexpressed population onto a spike at exactly 0 (17.4% of protein-coding cells).
Hart et al. (2013) zFPKM estimates the unexpressed mode by kernel density and its
spread from the right-hand half-normal; run on this data the unconstrained mode
search returned `μ = +4.13` log₂TPM — the **expressed** peak. Inverting the
pseudocount recovers the spread, but the low-TPM tail is then smeared thinly and
never out-densities the expressed peak. The rule that version actually tested was
`TPM ≥ 0.232`, arrived at by accident. Constraining the mode search left of the
per-line median fixes it (`μ = −6.65`). *A KDE mode is not a safe estimator of a
mixture component when a pseudocount has destroyed that component's width.*

**Result — the percentile is not an encoded boolean.** The deciding quantity is
the dependency ratio of the cell `{bottom decile} ∩ {rule says EXPRESSED}`. If
detection explained the gate, this cell would sit at 1.00.

| rule | median % excluded | dec1 ∩ EXPRESSED ratio | n genes |
|---|---|---|---|
| ref_q0.90 | 0.0% | 0.818 | 710 |
| ref_q0.95 | 0.0% | 0.810 | 694 |
| ref_q0.99 | 0.0% | 0.810 | 662 |
| ref_q0.999 | 2.7% | 0.827 | 635 |
| hart_constrained | 0.0% | 0.815 | 702 |
| TPM ≥ 1 | 0.0% | 0.822 | 632 |

Spread **0.017** across all six rules — insensitive to the threshold — and clears
the `|1 − ratio| ≥ 0.10` effect floor in 6/6. The detection increment (the mirror
cell) is *threshold-sensitive*, spanning 0.729–0.978 and clearing the floor in
only 3/6.

**Simpson reversal — read per gene, not pooled.** Pooled over 1,150,662 pairs the
same cell reads 0.932–1.043, i.e. "the percentile adds nothing". Pooling weights
by pair count, so the few tissue-restricted low-prevalence genes with large
not-expressed blocks dominate. The per-gene median gives each gene one vote and is
the correct aggregation here.

**Effect-size floors are mandatory at this n.** Pooled, an odds ratio of 0.870
returns `p = 3.1e-60`. Significance alone is uninformative on 10⁶ pairs; every
claim in this section clears both `p < 0.05` and an effect floor
(`MIN_RATIO_EFFECT = 0.10`, `MIN_PAUC_DELTA = 0.01`, `MAX_OR_EFFECT = 0.80`).

**The rules are structurally nested, which bounds what can be compared.** Decile 1
is exactly 10% of a gene's lines by construction, so:

```
rule excludes ≪ 10%  →  its exclusions sit inside decile 1   →  mirror cell empty
rule excludes ≫ 10%  →  decile 1 sits inside its exclusions  →  deciding cell empty
```

Measured containment confirms it: the fraction of a rule's exclusions falling
inside decile 1 runs 0.667–**1.000**, and for `TPM ≥ 1` it is exactly 1.000 —
perfect nesting, so the number of genes supporting a paired test is exactly **0**.
Only `ref_q0.999` (2.7% median exclusion, containment 0.234) supports one, at
n = 535. Small paired n in this design is *geometric, not statistical*. A further
practical fact: for 393–591 of 778 genes the detection rule excludes **nothing at
all**, so an absolute gate is inert on the majority of genes and could not replace
the percentile wholesale even if it matched it.

**Half-sample cross-validation kills the per-line score's apparent advantage.**
The per-line score's gate pAUC depends on its quantile centre, and the paired
delta against the percentile *changes sign* across the sweep:

| centre | per-line gate pAUC | delta (percentile − per-line) | p |
|---|---|---|---|
| Q0.90 | 0.5940 | **+0.0134** | 0.0068 |
| Q0.95 | 0.6178 | −0.0056 | 0.060 |
| Q0.99 | 0.6486 | −0.0301 | 1.2e-10 |
| Q0.999 | 0.6392 | −0.0280 | 1.2e-06 |

(percentile: 0.6076, fixed — it has no tunable parameter). Selecting the centre
and quoting its pAUC on the same labels is selection on the test set. Under 20
random half-splits of the 1,479 lines — centre chosen on half A, evaluated on
half B, **percentile recomputed within half B** because it is panel-dependent:

```
selection pAUC on half A (winner)   0.6578
held-out pAUC on B, per-line        0.6276
held-out pAUC on B, percentile      0.6200
held-out delta                     -0.0050   (range -0.0897 to +0.0662)

OPTIMISM (A - B)                   +0.0302
per-line wins on held-out data      10/20 reps
clears the effect floor             NO
```

The optimism `+0.0302` is numerically the entire claimed advantage `−0.0301`. On
held-out lines the two scores are indistinguishable and the per-line score wins in
exactly half the reps. Centre selection is stable (Q0.99 in 12/20, Q0.999 in 6/20)
— it just buys nothing out of sample.

*A secondary result, unanticipated:* the percentile is the **lower-variance**
estimator across panel composition. Held-out pAUC range 0.078 (0.5933–0.6713)
versus 0.170 (0.5541–0.7240) for the per-line score — notable given that
panel-dependence was the objection motivating the absolute alternative.

**Floor definitions have no headroom left.** The current exact-minimum floor gives
a dependency rate of 0.0086 against a Chronos label FPR of 0.0093 → residual
0.0000. `ref_q0.90` sits at 0.0093, exactly at the noise floor. Neither definition
leaves anything measurable.

*Verdict:* `percentile_increment_survives_effect_floor_but_pauc_is_tied`. The
quantity to quote is the **0.81 depletion ratio in the incremental cell**, not a
pAUC gap — the pAUC gap is precisely what failed to replicate.

---


---

## C.17 Paralog stratification — an annotation that failed its own validity check

> **Why retired.** An annotation that failed its own validity check. Never adopted.


`test_run_detection_vs_percentile.py`, Step 7b

Paralog buffering attenuates single-knockout phenotypes: a gene with a close
functional paralog can be deleted at little fitness cost. That depresses the
Chronos label for those genes independently of expression, so it is a candidate
confound for §C.16's incremental cell.

**No paralog annotation exists in this repository.** `gene_lookup.parquet` has no
such column; the only prior handling is §C.15's exclusion of KRTAP\* and VN1R\* as
paralogy-contaminated control families. The closest source on disk is HGNC
`gene_group` co-membership (`data/gene_with_protein_product.txt`), used here as an
explicit **proxy, not sequence identity**. Groups above 30 members are dropped
because the large ones are domain superfamilies rather than paralog families —
"Zinc fingers C2H2-type" has 763 members, "CD molecules" 385, and co-membership
there implies no paralogy at all. 1,862 of 1,968 groups survive; 10,615 of 19,173
genes get ≥1 paralog.

**Validation, and its failure.** Buffering predicts *lower* dependency prevalence
among paralogous genes. Measured:

```
prevalence WITH paralogs     0.1072   (n = 466)
prevalence WITHOUT paralogs  0.0629   (n = 312)
Mann-Whitney (paralog < singleton)   p = 0.7766
```

The effect runs the **wrong way**. HGNC gene groups skew toward well-characterised
families — ribosomal proteins, proteasome subunits, spliceosome components — which
are pan-essential, and that confound swamps buffering. What group membership
proxies is "belongs to a named family", which tracks essentiality, not redundancy.
The stratification was therefore **not interpreted as a paralog correction**, by
the guard built into the test rather than by later judgement.

**What can still be said.** The strata do not differ (deciding cell 0.829–0.842
paralogous vs 0.790–0.799 singleton, all `p > 0.18`), and restricting to
singletons — removing the suspect genes rather than adjusting for them — leaves
§C.16's finding intact and clearing the effect floor in 6/6 rules. So whatever
HGNC groups capture, it does not drive the result.

**The confound remains open, not cleared.** A real test needs sequence-level
paralog data: Ensembl Compara paralog pairs with % identity, DepMap's paralog
annotation, or the Dede/Ito paralog sets. None is present locally. See §9.

---


---

## C.18 Stratum centring and the two-layer score (protein is a net negative)

> **Why retired.** Tested; protein was a net negative. Not in any live path.


`test_run_stratum_centring.py`, `test_run_cna_gate.py`

**Stratum centring (§2.6 follow-up).** The 2-layer stratum genuinely sits higher
— median shift 0.210, positive in 86.3% of 1,562 genes — but correcting for it
makes performance *worse*, and RNA alone beats every Stouffer variant:

| label | RNA alone | Stouffer Z | Z centred |
|---|---|---|---|
| GDSC (AUROC) | **0.563** | 0.556 | 0.547 |
| Chronos (AUROC) | **0.519** | 0.517 | 0.509 |

Centring gain `−0.0145` (`p = 1.5e-17`) and `−0.0042` (`p = 2.4e-67`). Verdict:
`centring_does_not_recover__protein_is_a_net_negative`.

This bears directly on §5: the tier rule `n_layers == 2 → "high"` promotes rows on
a basis this test contradicts. It also reframes §0.3 — the Kish `n_eff = 1.37` for
E+P was already saying two layers carry less than two layers' evidence; this says
the second layer is worse than absent for ranking purposes.

**CNA deletion is a second, independent gate (`test_run_cna_gate.py`).** Over
17,278 genes × 754 lines = 13,027,612 pairs:

```
deleted pairs         8,467  (0.065% of the grid)
dep rate deleted      3.97%   vs base 10.44%   →  ratio 0.380
residual              3.04%   vs label FPR 0.93%   →  "check_calls"

incremental:  floor-not-deleted   ratio 0.119
              deleted-not-floor   ratio 0.466      Fisher OR 0.373, p = 1.7e-64
```

CNA carries information the expression floor does not. But the yield is 0.49 lines
excluded per gene, so it is usable only restricted to the **69 genes with ≥10
deleted lines** (CDKN2A 224, CDKN2B 202, MTAP 116). Verdict:
`cna_adds_information_beyond_expression`, usable form
`restrict to genes with >=10 deleted lines`.

---


---
