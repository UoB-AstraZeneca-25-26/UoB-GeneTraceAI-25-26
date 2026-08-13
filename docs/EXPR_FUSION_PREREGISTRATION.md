# PRE-REGISTRATION — Diagnostic Battery T1–T8
**Component under test:** multi-source RNA + orthogonal protein `core_score` architecture
(`05_bulk_rna_scorer.py` → `06_bulk_protein_scorer.py` → `02_core_score.ipynb` cells 7g/7h/7i)
**Status:** PRE-DECLARED. No test in this document has been run at time of writing.

| Field | Value |
|---|---|
| Document version | 1.0 |
| Author | F. — GeneTraceAI |
| Branch | `chai` |
| Commit at declaration | *(fill: `git rev-parse HEAD`)* |
| Declared UTC | *(fill: `date -u +%Y-%m-%dT%H:%M:%SZ`)* |
| Amendment policy | §9 |

> **Commit this file before running any test in it.** Its evidential value is
> the timestamp. A threshold written after seeing the number is not a threshold.

---

## 0. Why this battery exists

The architecture under test reinstates two components previously retracted on
measured grounds:
- the **protein arm**, excluded on reliability-ceiling and MNAR-missingness grounds;
- **Stouffer combination**, falsified by the measured DepMap–HPA redundancy (ρ ≈ 0.81).

Both were reinstated on the strength of a cleaner derivation, not a new
measurement. This battery supplies the missing measurement. The default outcome
if the battery is not completed is that both components remain retracted.

**Three arithmetic corrections are owed regardless of any result below** and are
not contingent on any test (§8).

---

## 1. Frozen state

Record before the first run. Any change to these invalidates results and
requires a re-run, not a re-interpretation.

```
git commit                    : ________
data/depmap_expr.parquet      : sha256 ________
data/hpa_rna.parquet          : sha256 ________
data/geo_matrix.parquet       : sha256 ________
data/procan_prot.parquet      : sha256 ________
data/ccle_prot.parquet        : sha256 ________
ground_truth_pairs.parquet    : sha256 ________   n_pairs = ________
lineage_map.parquet           : sha256 ________   n_lineages = ________
numpy / scipy / pandas        : ________
RNG seed (global)             : 20260813
```

**Excluded from this battery:** `chronos_long.parquet` contains Project Score
data rather than DepMap Chronos. No test here may use it. Any Chronos-based
external validation (Test A of the external sequence) is blocked until that file
is rebuilt and re-checksummed.

---

## 2. Primary endpoint and criterion

**Primary criterion:** agreement with the manually curated cell line × gene
ground-truth pairs. Current shipped value: **80.9%**.

**Criterion is fixed for the whole battery.** In particular:
- pAUC against CRISPR essentiality is **not** an endpoint here. The product
  selects cell lines for experimental design; it does not score vulnerabilities.
  Criterion reliability for CRISPR is 0.371, which caps observable effects
  regardless of method quality.
- No endpoint may be substituted after unblinding.

**Reported alongside every performance figure:** criterion reliability, the
UNKNOWN/NaN rate, and the number of scored pairs. A gain in agreement bought by
abstaining more is not a gain.

---

## 3. The null model (placebo spec)

One null is used throughout, so that all tests are commensurable.

### N1 — line-label permutation within (gene × source)

```
for each source s in {DepMap, HPA, GEO, ProCAN, CCLE}:
    for each gene g:
        observed_mask = isfinite(X[s][g, :])
        vals = X[s][g, observed_mask]
        X_null[s][g, observed_mask] = permute(vals, rng)
```

**Preserves exactly:** the missingness mask, per-(gene, source) coverage,
per-(gene, source) marginal distribution, `n_sources` for every (gene, line)
pair, lineage sizes, and the entire downstream pipeline.

**Destroys:** the association between a cell line's identity and its values;
hence all cross-source concordance and all biology.

**Rationale.** Under N1 there is no biological signal, but every structural
property that the machinery could exploit survives untouched. Any dependence of
the output on coverage, lineage size, or floor status that persists under N1 is
manufactured by the pipeline.

**Known limitation, stated in advance:** N1 also destroys cross-source
concordance, so it cannot separate "machinery exploits coverage" from
"machinery exploits concordance." A concordance-preserving null is not defined
here and is out of scope. This limitation is to be reported, not worked around.

### Replicates and seeds

- **B = 200** N1 replicates for T1 and T7; **B = 200** for T2's null.
- Seeds: `20260813 + replicate_index`, recorded per replicate.
- Bootstrap resamples where specified: **10,000**, resampled **over genes**, not
  over lines (lines are non-independent through the validated similarity graph).

### Placebo pass rule (global)

> Any effect claimed anywhere in this battery must be reproduced by **≤ 5%** of
> N1 replicates. An effect reproduced by more than 5% of placebo runs is
> withdrawn, not caveated.

---

## 4. Gating tests

T1 and T2 are gates. T3–T8 are conditional on them (§7).

---

### T1 — Completeness artefact (gate)

**Question.** Does the pipeline assign systematically larger `core_z` magnitude
to (gene, line) pairs simply because more sources happen to cover them?

**Predicted fault magnitude** (from the audit, to be confirmed or refuted):
plain Stouffer over three sources at ρ̄ ≈ 0.66 inflates sd to ≈ 1.52; the
single-source branch divides by √3 ≈ 1.732; ratio ≈ **2.64×**.

**Statistics, computed under N1** (biology absent by construction):

| ID | Statistic | Null expectation |
|---|---|---|
| T1a | `R_sd = sd(core_z | n_sources=3) / sd(core_z | n_sources=1)` | 1.00 |
| T1b | `R_abs = mean(|core_z| | n=3) / mean(|core_z| | n=1)` | 1.00 |
| T1c | Enrichment of 3-source pairs in the top 1,000 `core_score` vs base rate | 1.00 |
| T1d | Same as T1c for the bottom 1,000 | 1.00 |

Report each as median over 200 N1 replicates with a 95% percentile interval.

**PASS** — all of: `R_sd ∈ [0.90, 1.10]`, `R_abs ∈ [0.90, 1.10]`,
`T1c ≤ 1.10`, `T1d ≤ 1.10`.

**FAIL** — any statistic outside its band.

**KILL SWITCH K1.** On FAIL, the current within-layer combination is withdrawn.
It may not be patched by tuning `shrink`. Permitted replacements, in order:
1. Route coverage to an **evidence state**, not the score: `n_sources` becomes a
   reported field and a Kleene input; the score is computed on a fixed
   denominator (single-source or intersection panel).
2. Covariance-based combination (`poolr`, Strube/Hartung) with the measured
   covariance from T3, if and only if T3 passes.

Whichever replacement is chosen, **T1 is re-run against it before anything else
proceeds.** A replacement that has not cleared T1 is not eligible.

**Additional reporting requirement.** Regardless of PASS/FAIL, report the
Spearman correlation between `n_sources` and cell-line study frequency. Prior
measurement in this project: completeness–study-frequency ρ ≈ 0.60 (0.72 for
mutation coverage). If T1 fails and that correlation holds, the artefact is a
fame proxy and must be described as such in the write-up.

---

### T2 — Protein residual reproducibility (gate)

**Question.** Is `prot_resid` post-transcriptional biology, or platform noise?

**Design.** On the 291 lines shared by ProCAN and CCLE, compute the residual
independently within each platform:

```
for each gene g with both platforms and RNA:
    ρ_g^ProCAN  = shrunk Pearson(rna_z, procan_z)      # as in Stage 5 Step 1
    resid_ProCAN = procan_z − ρ_g^ProCAN × rna_z
    ρ_g^CCLE    = shrunk Pearson(rna_z, ccle_z)
    resid_CCLE   = ccle_z  − ρ_g^CCLE   × rna_z
    ρ_resid[g]   = Spearman(resid_ProCAN, resid_CCLE)  # across the 291 lines
```

Direction matters: this is **across lines, within gene**. A within-line,
across-gene correlation is inflated by the shared abundance profile and is not
an acceptable substitute.

**Primary statistic.** `ρ̃_resid` = median of `ρ_resid[g]` over genes with
n_overlap ≥ 30. Report with 95% bootstrap CI over genes (10,000 resamples).

**Reference values fixed in advance:**
- Raw cross-platform protein agreement: ρ = 0.373. This is the **ceiling**; the
  residual cannot reproduce better than the raw signal it is derived from.
- Empirical null: `ρ̃_resid` recomputed under N1 (200 replicates). Expected ≈ 0.

**PASS (protein arm may proceed)** — all of:
`ρ̃_resid ≥ 0.20`, lower bound of 95% CI `> 0.10`, and `ρ̃_resid` exceeds the
97.5th percentile of the N1 null.

**FAIL** — `ρ̃_resid < 0.10`, or the CI includes the null.

**INDETERMINATE** — `0.10 ≤ ρ̃_resid < 0.20` with CI excluding the null. Treated
as FAIL for shipping purposes; may be reported as a weak effect. The protein arm
does not enter the score on an indeterminate result.

**KILL SWITCH K2.** On FAIL or INDETERMINATE:
- The protein arm is removed from production scoring **permanently for this
  thesis**, and the removal is recorded as *measured*, not assumed.
- The orthogonalisation step (Stage 5, Steps 1–3) is withdrawn entirely — it has
  no purpose without a protein term.
- `W_PROT`, `RHO_PRIOR_CROSS`, and `N_MIN_RHO` are deleted from the parameter
  set rather than left dormant.
- The grounds recorded are: **reliability ceiling and MNAR missingness**, with
  T2 as the direct measurement. The phrase "protein is a net negative" is not
  reinstated (§8.4).

**Note fixed in advance so it cannot be argued after the fact.** ProCAN–CCLE
agreement (0.373) is at or below the literature range for mRNA–protein
correlation (0.37–0.46). Two measurements of the same molecule agreeing no
better than two different molecules is a statement about measurement error. The
disattenuated ceiling on the protein layer's correlation with any external truth
is √0.373 ≈ 0.61, and that is an unreachable bound, not a target.

---

## 5. Conditional tests

Run only if the gates permit (§7).

### T3 — Effective sample size of the RNA layer

**Question.** Is `W_RNA = √2.0 = 1.414` supported by the measured correlation
structure?

**Method.** Per-gene pairwise Spearman across shared lines for all three RNA
pairs (DepMap–HPA, DepMap–GEO, HPA–GEO), across lines within gene. Report median
and IQR per pair. Compute `ρ̄` as the mean of the three medians, then
`n_eff = 3 / (1 + 2ρ̄)`.

**Pre-declared reference table:**

| ρ̄ | n_eff | √n_eff |
|---|---|---|
| 0.25 | 2.00 | 1.414 (the claimed value) |
| 0.40 | 1.67 | 1.29 |
| 0.50 | 1.50 | 1.22 |
| 0.60 | 1.36 | 1.17 |
| 0.66 | 1.29 | 1.14 |

**PASS** — `n_eff ≥ 1.60`, i.e. the claimed weight is within ~10%.

**FAIL** — `n_eff < 1.60`. `W_RNA` is recomputed from the measurement. Note that
at ρ̄ ≥ 0.50, `√n_eff < W_PROT = 1.204`, which **inverts the relative weighting
of the two layers**. If T2 has already failed this is moot; if T2 passed, the
inversion must be applied, not discussed.

**Secondary test T3b — is the Kish scalar adequate?** Compute the combination
weights two ways: (i) Kish scalar `ρ̄`; (ii) explicit covariance via
`poolr` (Strube / Hartung). **If the resulting weights differ by > 10%, the Kish
scalar is abandoned** in favour of the explicit covariance. Kish's design effect
is a cluster-sampling result being used here by analogy and assumes exchangeable,
uniformly correlated sources; measured ρ spans 0.373–0.812, so the assumption is
already known to be false. T3b decides whether that matters numerically.

---

### T4 — Calibration of Φ(core_z)

**Question.** Is `core_score = Φ(core_z)` a calibrated probability, or a
monotone relabelling with a false confidence claim attached?

**Statistics:**

| ID | Statistic |
|---|---|
| T4a | Empirical sd, skew, excess kurtosis of `core_z` |
| T4b | Kolmogorov–Smirnov D between ECDF(`core_z`) and Φ |
| T4c | Observed / expected exceedance ratio at z = 2, 3, 4 |

**PASS** — `D ≤ 0.05` **and** all three exceedance ratios ∈ [0.5, 2.0].

**FAIL** — otherwise. `Φ` is replaced by the **empirical CDF** of `core_z`. This
preserves the identical ordering (both are monotone) while dropping the
normality claim.

**Fixed in advance.** Prior measurement on comparable data gave sd ≈ 1.606,
skew ≈ −4.54, kurtosis ≈ 78. At sd 1.6 the document's claim of
Φ(4) = 0.99997 is actually ≈ 0.9938 — the complement is wrong by a factor of
~200. If T4 fails, the marketing claim "magnitude preserved: z=4 → 0.99997" is
struck from all documents, and the ECDF is described as what it is: the same
ordering as `rank(pct=True)`, honestly labelled.

---

### T5 — Floor artefact in the lineage z

**Question.** Do the pipeline's highest-confidence calls concentrate on
(gene, lineage) strata where `MAD_FLOOR` is doing the work rather than the data?

**Strata definition, fixed in advance:**
- `FLOORED` — the (gene, lineage) stratum where raw MAD < `MAD_FLOOR`, i.e. the
  floor was binding.
- `SILENT` — lineage median within one floor-unit of the source detection floor
  (per-source floor from `source_floor_check`).
- `NORMAL` — neither.

**Statistic.** Enrichment of the top 1,000 `core_score` pairs in `FLOORED` and
in `SILENT` relative to their base rates.

**PASS** — both enrichments ≤ 1.25.

**FAIL** — otherwise. The silence guard is wired in as a **gate on the z**, not
as an annotation: where the lineage is silent, the output is `UNKNOWN`
("not expressed in this lineage"), never a large positive z.

**Mechanism being tested, stated so the result is interpretable either way.**
For a gene silent across a lineage, median ≈ 0 and MAD → floor 0.384, so a line
at log₂ value 2.0 receives z ≈ +5.2 → Φ ≈ 1.0. If T5 fails, the pipeline's most
confident calls are a floor artefact operating precisely where the biology is
"not expressed."

---

### T6 — Parameter sweeps (MIN_PEERS, MAD_FLOOR)

**Question.** Are the two free parameters load-bearing?

Both are swept in full. `MAD_FLOOR` is treated as a **free empirical parameter**,
because its stated closed form does not evaluate to its stated value (§8.1).

| Parameter | Sweep |
|---|---|
| `MIN_PEERS` | 5, 10, 15, 20, 30 |
| `MAD_FLOOR` | 0.10, 0.20, 0.384, 0.50, 0.75 |

Run the full 5 × 5 grid. For each cell report: ground-truth agreement, UNKNOWN
rate, scored-pair count, and Spearman rank correlation of the output against the
most conservative cell (`MIN_PEERS=30`, `MAD_FLOOR=0.75`).

**PASS** — ground-truth agreement varies by ≤ 3.0 percentage points across the
whole grid, **and** no adjacent-cell rank correlation falls below 0.90, **and**
no direction of effect reverses anywhere in the grid.

**FAIL** — otherwise. The adopted setting is then the **most conservative cell
that passes**, not the best-performing cell. Selecting the maximum of a sweep is
selecting noise.

**Fixed in advance.** `MIN_PEERS = 5` means a MAD estimated from five points,
which has breakdown 0.5 (two aberrant lines destroy the reference) and very high
sampling variance. Instability below ~15 would be the expected result, not a
surprise.

---

### T7 — End-to-end placebo

**Question.** Does the complete pipeline, run on data with no biology, reproduce
any effect the architecture claims?

**Method.** Full pipeline — Stage 1 through `core_score.parquet` — on N1
permuted inputs, 200 replicates. Recompute every headline statistic the
architecture asserts, including the primary endpoint of T8.

**PASS** — no claimed effect is reproduced by more than 5% of replicates.

**FAIL** — otherwise. The specific claim reproduced by placebo is withdrawn.

**Reporting requirement.** The placebo distribution of the primary endpoint is
reported as a figure alongside the observed value in all write-ups. Not as an
appendix.

---

### T8 — Head-to-head against M0

**Question.** Does the redesign beat the thing it replaces?

**M0 (comparator, frozen):** single-source DepMap within-gene percentile,
`expr.rank(axis=0, pct=True, na_option="keep")` — the current shipped E term.

**Primary endpoint:** ground-truth agreement. Baseline **80.9%**.

**Pre-declared superiority margin: δ = +3.0 percentage points.**
Ship threshold: ≥ **83.9%**.

**Simultaneous constraints — all must hold:**

| Constraint | Threshold |
|---|---|
| UNKNOWN / NaN rate increase | ≤ +5.0 pp absolute |
| Scored (gene, line) pair loss | ≤ 10% relative to M0 |
| T7 placebo reproduction of the gain | ≤ 5% of replicates |

**Inference.** Paired bootstrap **over genes**, 10,000 resamples, one-sided.
Superiority requires the lower bound of the 95% interval on the difference to
exceed 0.

**PASS** — margin met and all constraints hold.

**FAIL** — the redesign does not ship. M0 remains the shipped E term.

**KILL SWITCH K3.** On FAIL, the redesign is recorded as a **measured negative
result**: specified, audited before running, failed a pre-declared threshold.
It is written up as such. It is not re-tuned and re-tested against the same
ground truth — that would burn the only criterion available.

---

## 6. Multiplicity and re-use of the criterion

The curated ground truth is used as an endpoint in **T6 and T8 only**. It is not
consulted during T1–T5, T7, or during the choice of any replacement under K1.

Any test re-run after a K1 or K2 replacement counts as a **new test** and is
reported with its own row in the results table, including the fact that it is a
re-run and what changed. Re-runs are capped at **two per test**. A component
requiring a third attempt is withdrawn.

---

## 7. Decision tree

```
                    ┌─────────────────┐
                    │  T1  (coverage) │
                    └────────┬────────┘
              FAIL ──────────┴────────── PASS
                │                          │
        Apply K1: withdraw within-     ┌────┴────┐
        layer combination; choose      │   T2    │
        replacement; RE-RUN T1         │(protein)│
                │                      └────┬────┘
        still FAIL → STOP.        FAIL ─────┴───── PASS
        Ship M0 + Kleene.           │                │
        Both components stay     Apply K2:        Run T3, T3b
        retracted.               protein arm      (n_eff)
                                 removed              │
                                 permanently.         │
                                     │                │
                                 RNA-only path        │
                                     └────────┬───────┘
                                              │
                                    Run T4, T5, T6 (any order)
                                    Apply each FAIL's remedy
                                              │
                                         Run T7 (placebo)
                                              │
                                     FAIL → withdraw claim
                                              │
                                         Run T8 vs M0
                                              │
                                    ┌─────────┴─────────┐
                                  PASS               FAIL
                                    │                  │
                              Ship redesign      Apply K3.
                              with all           Ship M0.
                              remedies applied   Write up as
                                                 measured negative.
```

**Default outcome if the battery is abandoned part-way:** M0 + Kleene ships;
protein and Stouffer remain retracted. Incomplete evidence does not reinstate a
retracted component.

---

## 8. Corrections owed unconditionally

These are arithmetic and documentation errors. They are not contingent on any
test result and are to be fixed before the battery is run, so that the tests are
run against a correctly described system.

### 8.1 `MAD_FLOOR` closed form is false

Claimed: `0.384 = 1/(2√2·Φ⁻¹(0.75))`. Actual: `1/(2√2 × 0.674490) = 0.524179`.
No arrangement of those terms yields 0.384.

**Action.** Either state the true derivation, or label 0.384 as empirically
calibrated. If empirical, it is a free parameter and T6 sweeps it.

### 8.2 `W_RNA` is overstated

`n_eff = 2.0` for three sources requires ρ̄ = 0.25, when one pair is already
measured at 0.812. Cannot hold. T3 supplies the correct value; the claim is
struck from all documents immediately regardless.

### 8.3 The Φ tail claim is wrong

"z = 4 → 0.99997" holds only at sd = 1. Struck now; T4 determines the
replacement.

### 8.4 "Protein is a net negative"

Still present in `STATE_REPORT.md` and `RISK_REGISTER.md`, resting on a global
AUROC the runbook declares irrelevant. Replace with:

> Protein was excluded on measurement-reliability grounds (cross-platform
> reproducibility bounds the achievable correlation) and MNAR coverage bias, not
> on a discrimination metric.

If T2 passes, this becomes: "excluded pending T2, which it passed."

### 8.5 Reporting template operating characteristics

Correct values: **α = 0.116**, **85% recovery at 3 MAD**. Template currently
cites otherwise.

### 8.6 Citation claims requiring correction

- "MAD is 4.5× more efficient than SD for contaminated Gaussian data
  (Hampel 1974)" — MAD's asymptotic efficiency at the pure Gaussian is ~37%,
  i.e. *lower* than SD; its virtue is a 0.5 breakdown point. Either state the
  exact contamination fraction the 4.5× refers to, or cite Leys et al. 2013 for
  robustness instead of efficiency.
- "Stouffer is BLUE when weights equal √n (Hedges & Olkin 1985)" — true for
  study-level effect estimates with sample sizes n_i. Here each source
  contributes one value per (gene, line); panel size does not govern that
  value's precision. The claim does not transfer.
- "Durbin 2002" — that is a variance-stabilising transform; `log₂(x+1)` is not
  it, and the +1 pseudocount is not variance-stabilising at low counts.

### 8.7 Cross-lineage comparability contradiction

The architecture states cross-lineage comparison is invalid, then ranks across
lineages with `stratum_rank` grouped by `n_layers` only. Either add lineage to
the stratum key, or state in the CLI that the score is interpretable only within
a lineage. Pick one before the battery runs.

---

## 9. Amendment policy

1. Amendments before the first run: edit freely, bump the version, re-commit.
2. Amendments after any test has been run: **append only**, in §10, with UTC
   timestamp, the reason, and which results were already visible at the time of
   the amendment.
3. Thresholds in §4–§5 may not be changed after the corresponding test has been
   run. A threshold that turns out to be badly chosen is reported as such and
   the test is reported as failed.
4. Endpoints may not be substituted. If the criterion proves unusable, the
   battery halts and the reason is recorded; no replacement endpoint is
   introduced mid-battery.

---

## 10. Results log

Fill only after running. One row per test. Do not edit §4–§5.

| Test | Run UTC | Commit | Statistic | Threshold | Observed | Verdict | Kill switch fired |
|---|---|---|---|---|---|---|---|
| T1 | 2026-08-13 | see t1_completeness_result.json | R_sd, R_abs, T1c, T1d (N1, B=200) | [0.90,1.10] / ≤1.10 | 1.771 / 1.980 / 2.462 / 2.532 | FAIL | K1 fired |
| T2 | | | | | | | |
| T3 | | | | | | | |
| T3b | | | | | | | |
| T4 | | | | | | | |
| T5 | | | | | | | |
| T6 | | | | | | | |
| T7 | | | | | | | |
| T8 | | | | | | | |

**Amendments:**
*(none at declaration)*

---

## Appendix A — Statistic definitions

```python
# T1a — coverage scale artefact under N1
R_sd = core_z[n_sources == 3].std() / core_z[n_sources == 1].std()

# T1c — top-1000 enrichment
base  = (n_sources == 3).mean()
top   = core_score.nlargest(1000).index
T1c   = (n_sources.loc[top] == 3).mean() / base

# T2 — residual reproducibility, ACROSS LINES within gene
rho_resid[g] = spearmanr(resid_ProCAN[g, shared_lines],
                         resid_CCLE[g, shared_lines]).correlation
rho_tilde    = median(rho_resid[n_overlap >= 30])

# T3 — effective n
rho_bar = mean([median(rho_DepMap_HPA), median(rho_DepMap_GEO),
                median(rho_HPA_GEO)])
n_eff   = 3 / (1 + 2 * rho_bar)

# T4b — calibration
D = ks_2samp(core_z, norm.rvs(size=len(core_z), random_state=seed)).statistic

# T8 — paired bootstrap over GENES
diff = agreement(redesign) - agreement(M0)
boot = [diff on genes resampled with replacement for _ in range(10_000)]
ship = percentile(boot, 2.5) > 0 and diff >= 0.030
```

## Appendix B — What each gate outcome means for the write-up

| Outcome | Chapter framing |
|---|---|
| T1 FAIL | Absence-as-negative recurs in a fourth location, this time via shrinkage asymmetry rather than zero-filling. Detected by permutation control before shipping. Strengthens the recurring-fault-class finding. |
| T2 FAIL | The protein exclusion becomes a *measured* decision with a direct test, not an inference from published reliability figures. Stronger than the current justification. |
| T1 PASS, T2 PASS, T8 FAIL | Principled math did not beat the simple baseline. Consistent with the F1 lineage-z null and with Franks et al. 2018 on z-score instability. A result, reported as one. |
| All PASS | The redesign ships, with every constant measured rather than asserted, and with the placebo distribution published beside the headline figure. |
