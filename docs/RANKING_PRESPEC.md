# RANKING_PRESPEC — adoption rule for the Part A re-run

**Written 2026-08-09, before any re-run number exists.** Standing rule 6. The only inputs
consulted before writing this were the structure and scale of
`data/DepMap_Chronos/CCLEGeneCopyNumber20Q2.hdf5` (shape, dtype, value distribution,
ploidy-normalisation status) — none of which is a Part A outcome.

This document is the guard against reading a favourable result out of a noisy one. The
project has a documented history of parameter direction reversing under honest evaluation,
and of a conclusion surviving only because the sample was selected on the outcome
(`TRANSCRIPTOMICS_VALIDATION_REPAIR.md` B3, where an apparent +1.273 effect became −0.043
on a random sample).

---

## 1. What is being tested

**Question.** With a gene-complete copy-number check wired in — so that the confirmation
count has **two** near-universally assessed, materially discriminating inputs instead of one —
does the count become a viable ordering basis?

**Why this and not something else.** The previous Part A returned negative for a specific,
diagnosed reason: on gate-retained pairs, `cna_gate` was not assessed for 99.72% of pairs,
`protein` for 80.97%, `geo_corrob` for 67.92%, and `alteration` passed 97.15% of the time.
58.0% of retained pairs had exactly one non-gate check assessed. The count could not
discriminate because there was nothing to count. The CN file is the only identified lever
that changes that, and it is already on disk.

**What would still not be fixed if this succeeds.** A viable count is necessary, not
sufficient. Adoption for the product additionally requires the evaluation denominator fix
and a positive held-out comparison against a random ordering, neither of which is in scope
here.

---

## 2. The CN check — parameters fixed a priori

`CCLEGeneCopyNumber20Q2.hdf5`: 908 cell lines × 27,639 genes, **zero missing values**, scale
`log2(relative_CN + 1)` with neutral = 1.0. Verified already normalised to each line's own
ploidy (Spearman(per-line mean, GDSC ploidy) = +0.115; per-line mean sd = 0.020), so no
external ploidy correction is applied and the 30%-of-panel ploidy gap does not bite.

| Parameter | Value | Basis |
|---|---|---|
| **deletion call (primary)** | relative CN < **0.5** | ≈ loss of at least half the copies relative to the line's own ploidy — one-copy loss from diploid. The standard CCLE/DepMap relative-CN loss convention |
| deletion call (sensitivity) | 0.25, 0.75 | reported alongside; **not** used to choose |
| amplification call | relative CN > 2.0 | reported, not used in the gate |

**The deletion threshold is fixed before measuring and is NOT selected on any Part A
metric.** Selecting it on the outcome would repeat the B3 error this project has already
made once. The two sensitivity values are reported so the reader can see the result is not
an artefact of the cut; if the verdict flips between 0.25 and 0.75, that is itself a
reportable instability and the verdict becomes *inconclusive*, not *adopt*.

---

## 3. The primary statistic

**Primary variant: the raw confirmation count, restricted to pairs with ≥2 assessed checks.**

Rationale, fixed now: the raw count is the quantity the design actually proposes; the
≥2-assessed restriction is what makes "count of independent confirmations" meaningful at all,
and with a gene-complete CN layer that restriction should cost almost nothing (which is
itself a pre-registered prediction — see §6).

**Secondary variants, reported but not decisive:** the unrestricted raw count, and the
fraction-of-assessable `n_pass / n_assessed`.

**Evaluation set, fixed now:** valid genes only (silence guard, expressed in ≥20% of the
panel), measured **within the gate-retained set** — the level at which the count is proposed
to operate. Per standing rule 3, all-genes and all-lines figures are reported alongside, but
the valid/retained figure is the one the decision is made on.

**The gate check is excluded from the count it is being used to subdivide.** A check cannot
confirm the tier it defined.

---

## 4. The decision rule

Three criteria, all measured on the primary statistic, valid genes, within the retained set.
For reference, the previous run's values are shown — these are what "no change" looks like.

| # | Criterion | **ADOPT if** | **ABANDON if** | Previous run |
|---|---|---|---|---|
| **A1** | \|Spearman(coverage_count, count)\| — the count must not be a data-availability proxy | **≤ 0.20** | **≥ 0.40** | +0.583 |
| **A2a** | median largest tied group, as a fraction of retained lines | **≤ 0.50** | **≥ 0.65** | 0.705 |
| **A2b** | median number of distinct count values achieved | **≥ 4** | ≤ 2 | 4 |
| **A3** | median Spearman between two random genes' orderings — the count must not return the same list for every gene | **≤ 0.25** | **≥ 0.40** | +0.531 |

**ADOPT** requires **all four** to clear the adopt threshold.
**ABANDON** if **any one** hits its abandon threshold.
**INCONCLUSIVE** otherwise — report, do not adopt, and state what further coverage would be
needed. An inconclusive result is not a licence to adopt with caveats.

### Why these numbers

- **A1 ≤ 0.20.** ρ = 0.20 is 4% of variance shared with data availability — small enough that
  coverage is not a leading explanation of the ordering. ρ = 0.40 is 16%, at which point
  coverage is a leading explanation and the count is a proxy wearing a biology label.
- **A2a ≤ 0.50.** If the largest tier is at most half the retained set, the count genuinely
  partitions. At ≥ 0.65 two thirds of lines share one value and it is two buckets, not an
  ordering. This threshold is about the *user-facing* claim: a list where 65% of entries are
  tied should not be presented as a list.
- **A2b ≥ 4.** With the gate excluded there are 4 remaining checks and 5 possible values. Four
  distinct values achieved is the minimum for the count to express more than a coarse
  trichotomy.
- **A3 ≤ 0.25.** Below this, different genes return substantially different orderings — the
  product answers the question asked. At ≥ 0.40 the user gets the same shortlist regardless
  of gene, which is a product failure even if the statistic is defensible.

These are stated as absolute thresholds rather than "better than last time" deliberately. An
improvement from +0.583 to +0.45 would be a real improvement and still a coverage proxy.

---

## 5. Sensitivity and stopping conditions, fixed now

- All four criteria are reported at each of the three CN deletion thresholds
  (0.25 / **0.5** / 0.75). The decision is read off **0.5**. If the verdict differs across the
  three, the result is **INCONCLUSIVE** regardless of what 0.5 says.
- Gene sample: 300 protein-coding genes, seed 42 — **the same sample as the previous run**, so
  the comparison is paired and the change is attributable to the CN layer rather than to
  resampling.
- Additionally reported: the four criteria on the previous 5-check configuration recomputed
  on the same sample, so the before/after is like-for-like.
- No other variant will be introduced after seeing the numbers. If none of the pre-specified
  variants clears the bar, the answer is that the count does not work, not that a fifth
  variant should be tried.

---

## 6. Pre-registered predictions

Recorded so that being wrong is visible rather than absorbed.

1. **`cna_gate` assessed-fraction will rise from 0.28% to above 55% of retained pairs.** The
   CN matrix covers 908 lines against a ~1,543-line panel, so the ceiling is roughly 59%, not
   100%. **This is the main risk to the whole exercise**: a check assessed on 59% of lines
   cannot make the ≥2-assessed restriction free, and the pairs it misses will not be a random
   subset of the panel.
2. **A1 will improve but may not clear 0.20.** Adding a second check that is itself unevenly
   covered can move the confound in either direction.
3. **A2a is the criterion most likely to fail.** Two binary checks yield at most three count
   values on pairs where only those two are assessed, and the modal value will be large.
4. **A3 will improve most**, because the CN check is gene-specific in a way that data
   availability is not.

If the result contradicts these, the contradiction is reported as a finding.

---

## 7. What happens on each verdict

| Verdict | Action |
|---|---|
| **ADOPT** | Proceed to Parts B–E of the ranking brief: specify the confirmation set with the independence rule, the protein bands, the alteration tier, and the assembled ordering. Adoption is still conditional on the denominator fix and a positive held-out comparison |
| **INCONCLUSIVE** | Report the numbers, do not build. State precisely which criterion failed and what coverage would move it |
| **ABANDON** | The count is not a viable ordering basis even with the best available lever. Report, and redirect to the fallback: gate tiers with honest ties, and a coverage roadmap |

In all three cases the coverage finding is written up as a standalone result
(`COVERAGE_RESULT.md`), and the evaluation denominator fix proceeds regardless — it is
required for every future comparison independent of this outcome.
