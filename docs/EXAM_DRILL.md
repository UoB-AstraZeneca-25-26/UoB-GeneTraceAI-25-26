# GeneTraceAI — Viva / Review Drill Sheet

> ## ⚠ RETRACTION — the gate-region pAUC in this document is withdrawn
>
> **Issued 2026-08-10.** This is an EMPIRICAL retraction: the formula in
> §C.15 is mathematically correct and is still used in `gate_audit/`. What
> is withdrawn is the *claim* the formula was used to support. No arithmetic
> error is involved.
>
> This document cites **gate-region pAUC 0.596–0.608**, p = 2.8e-76, as verified
> and replicated. **That result does not stand.** It was measured before
> `gate_audit/` ran.
>
> | | |
> |---|---|
> | decile-1 ratio, as published | 0.8322 |
> | decile-1 ratio, screened lines only | **0.9875** |
> | placebo with dependency labels randomised | 0.8350 — **98% of the published depletion** |
>
> The effect was an artefact of counting never-screened cell lines as confirmed
> non-dependencies (`gate_audit/03_denominator.py`). It reproduces at DepMap 24Q4,
> so it was never a small-sample effect (`gate_audit/08`). What survives is
> CRISPR-specific: flat on GDSC2 drug response, and **opposite** on RNAi at
> matched prevalence (`gate_audit/09`).
>
> **Do not cite this number.** The claim it supports is withdrawn. Statements
> below about the gate being "the only place this architecture has replicated
> above-chance evidence" are superseded — the validated result in this project is
> the cell-line similarity graph (+0.186 [+0.149, +0.227] drug-response
> concordance over a same-tissue baseline, `cell_similarity/08`) and the Kleene
> filter (80.9% on 53 published cell-line facts,
> `docs/MULTIGENE_GROUND_TRUTH.md`).
>
> Retained rather than deleted so the record shows what was believed, when, and
> what changed it.

28 questions most likely to be asked, with the answer and where it lives in the repo.
Use this **after** reading `MATH_REFERENCE.md` end-to-end — it's a self-test, not a substitute.

Marked **[DANGER]** = a question where the honest answer is uncomfortable and you must
be able to give it without hesitating. These are the ones that lose marks when fumbled.

Line references: `MR` = `docs/MATH_REFERENCE.md`, `SV` = `src/pipeline/scoring_variants.py`.

---

## A. Core score — the mathematical centre

**A1. Why percentile rank and not z-score or min-max?**
Expression is log2-TPM, proteomics is a log-ratio — different, non-comparable native
scales. Percentile uses only *within-gene ordering*, so it's robust to outliers and to
the two assays having different distributions. Min-max was measured as the weakest.
→ `MR:79-95`, `SV:22-43`

**A2. What exactly is being ranked, and across what axis?**
`rank(axis=0)` — per **gene**, across **cell lines**. So `core_score = 0.9` means "this
cell line is in the top 10% of cell lines *for this gene*", not "this gene is highly
expressed in this cell line." Getting this axis backwards is the classic stumble.
→ `SV:30-31`

**A3. Why did you abandon Noisy-OR?**
Noisy-OR assumes independent evidence sources. mRNA and protein are **not** independent
(ρ≈0.46–0.58, Nusinow et al. 2020), so it double-counts agreeing evidence and inflates.
Measured: **22.5%** of 2-layer rows scored >0.95, and 379 rows pinned at exactly 1.0.
→ `MR:106-126`

**A4. What replaced it, and did it actually fix it?**
Correlation-penalised weighted sum: `w_i = 1/(1+|ρ̄_i|)`, normalised to sum to 1 — a
layer that's redundant with the others gets down-weighted. Verified on production data:
the 22.5% stratum drops to **2.32%** (I re-verified this against `core_score.parquet`),
and the 379 pinned rows now score <0.5.
→ `MR:127-174`, `SV:75-136`

**A5. Where does ρ_EP = 0.46 come from — did you fit it?**
No — it's taken from prior literature/validation, not fitted on this data. ρ_EC = ρ_PC =
0.005 *were* measured directly in this project. Be honest that 0.46 is an assumed
constant; the defence is that it's a published, widely-replicated value and the formula
is monotonic in it (small changes barely move the weights).
→ `MR:159-161`

**A6. [DANGER] What fraction of your data actually uses this formula?**
**9.61%.** Only 2,730,740 of 28,403,110 rows have two layers, because proteomics covers
**369 of 1,485** scored cell lines (~25%). For the other **90.39%**, `core_score` is just
the expression percentile — the correlation penalty never fires. Say this plainly: the
formula is correct and necessary where both layers exist, but the dominant regime in this
dataset is single-layer. This is a **data-availability limit, not a method flaw**, and it
is exactly what `n_layers` and `stratum_rank` exist to make visible rather than hide.
→ verified from `core_score.parquet`; `MR:194-210`

**A7. Why `stratum_rank` instead of just comparing `core_score`?**
A cell line with `n_layers=2, core_score=0.7` rests on more evidence than one with
`n_layers=1, core_score=0.7` — they aren't comparable. `stratum_rank` is a percentile
computed *within* the (gene, n_layers) group, so you only ever compare like with like.
Cross-stratum calibration is explicitly deferred, not silently ignored.
→ `MR:194-210`, `SV:172-175`

**A8. Why is a missing layer NaN rather than 0?**
0 would mean "measured, and it's the lowest" — a strong negative claim. NaN means "not
measured." The code reindexes to the union of models/genes specifically so structural
absence stays NaN and contributes *no term* rather than a zero term.
→ `SV:8-13, 46-72`

**A9. Why mean (not max, not drop) for the 16 duplicate-profile cell lines?**
Least lossy: max would inflate, dropping loses data. → `MR:96-105`

---

## B. Track C — mutations & fusions (the "why these columns" zone)

**B1. What is the Hill function doing, and why is it needed at all?**
`Hill(x,p0,k) = x^k/(x^k+p0^k)` — an S-shaped map from a raw feature onto (0,1) with
midpoint at `x=p0`. It's needed because Noisy-OR requires its inputs to be
*probabilities*; the raw features (VEP rank 0–3, pathogenicity 0–1, variant count) aren't.
Hill bounds each into (0,1) with a tunable threshold and steepness.
→ `MR:413-421`

**B2. Why *these three* columns for `p_mutation`?**
They are three **independent kinds of evidence** that a mutation is damaging:
- `max_vep_rank` (p0=3, k=2) — predicted *consequence severity* (is it nonsense/missense/silent)
- `max_pathogenicity` (p0=0.5, k=2) — predicted *functional impact*, max(AlphaMissense, REVEL)
- `variant_burden` (p0=3, k=1.5) — *how many* variants, `log1p(variant_count)`

Severity, impact, and recurrence are conceptually distinct, which is what licenses the
Noisy-OR combination. → `MR:422-444`

**B3. Why is burden kept separate from quality (`p_mutation_quality` vs `p_mutation_burden`)?**
Because the signature-discount step scales **burden but not quality**:
`p_mutation_adj = NoisyOR(quality, m_mut·burden) + driver_boost`. A hypermutated line's
high burden should be discounted; the severity of a specific variant should not be.
→ `MR:422-444`. **[DANGER]** — see D4: this split is *not* in the shipped file.

**B4. Why the flat +0.15 driver boost — why additive, why 0.15?**
Additive so a known oncogene/TSG hit lifts the score regardless of where the continuous
components landed. **0.15 is locked but not fitted** — say so. It's a stated placeholder
pending calibration against CRISPR/GDSC anchors. → `MR:429`, `MR:670-686`

**B5. [DANGER] Are the Hill parameters calibrated?**
**No.** Every `(p0,k)` pair in Track C is an explicit placeholder, flagged in-source as
"TEMP SCORING … pending calibration against CRISPR/GDSC anchors." Don't let this be
*discovered* — state it as a known limitation with a named remedy. → `MR:670-686`

**B6. Why different `k` for fusions?**
`p_ffpm` and `p_recur` use gentler `k=1.5` because those distributions are heavily
right-skewed (87% of recurrence mass sits at percentile 1); a steep `k` would saturate
almost everything to 1. `p_conf` uses p0=2.0 to place the midpoint at "medium"
confidence. → `MR:445-465`

**B7. Why the in-frame boost for fusions?**
In-frame fusions preserve the reading frame, so they can actually be translated into a
functional chimeric protein — out-of-frame usually can't. → `MR:456-459`

**B8. [DANGER] Does any of this Track C scoring actually affect your rankings?**
**No — and you must separate two things when answering.** Track C's *cleaning* outputs
are fully load-bearing: `any_driver | oncogene_hit | tsg_hit` from
`mutations_collapsed.parquet` and `max_confidence=="high"` from `fusions_gene_level.parquet`
are what build `has_driver_alteration`, the production gate. Every driver-gated ranking
decision runs on those Booleans.

Track C's *scoring* outputs (`p_mutation`, `p_fusion`) do **not** enter the sort order.
The gate is `has_driver_alteration = mut_driver | fusion_driver` — purely Boolean, never
reads them. Their only live consumer is `explain_pair.py:137-162`, which shows them to the
user tagged in-code as `"CONFIDENCE MODIFIER — does not affect stratum_rank"`. They also
feed `has_alteration` (p>0.5), the rejected any-alteration arm.

So: not unused, but display-only. This is *why* B5's uncalibrated Hill parameters matter
less than they look — and why the signature-discount test (D4) came back null.
→ `MR` Track C section, `explain_pair.py:137-162`

---

## C. Routing, confidence, evaluation

**C1. [DANGER] Driver-gating has a *lower* hit rate than any-alteration gating. Why did you pick it?**
Correct — `hr_any_alt = 0.0252` vs `hr_driver = 0.0212` vs flat `0.0192`. Any-alteration
wins on the point estimate. The reason for driver-gating is **stability**: the standard
deviation of the per-gene lift drops from 0.0147 to 0.0082, a **44.6% variance
reduction**. Any-alteration's higher mean comes with nearly double the spread — it helps
some genes and hurts others. Driver-gating is the more conservative, more generalisable
choice. Own the trade-off; don't pretend the bigger number isn't there.
→ `MR:283-297`, `outputs/stage4_eval_summary.json`

**C2. Why is confidence categorical and not a 0–1 score?**
Because a numeric confidence would imply a calibrated probability that was never fitted.
The tiers are deterministic rules with an explicit precedence order over four distinct
evidence types (GDSC validation, COSMIC role prior, CNA mechanistic confirmation,
Chronos cross-check). "Confidence is categorical and deterministic — not a fabricated
0–1 score." → `MR:300-325`

**C3. [DANGER] 92.5% of your 28.4M predictions are confidence="unknown" and only 0.03% are "high". Is the output useful?**
Yes, and the distribution is the honest consequence of the design. Only 141 genes have
curated GDSC ground truth; `regime_source` is `measured` for 209K rows, `prior` for
772K, `unknown` for 27.4M. Coverage is genome-wide *by construction* (the score is
computable for any gene with expression), but **validated** confidence is bounded by
ground-truth availability. The tiering exists precisely so a user can tell those apart
rather than trusting all 28.4M rows equally.
→ verified from `predictions_with_confidence.parquet`

**C4. [DANGER] Why invert the score for tumour suppressors — and does it work?**
**The rationale is sound but the measurement contradicts it. Volunteer this; do not get
caught by it.**

*The rationale:* for a TSG, *loss* is the actionable signal — low abundance, not high.
Inverting at build time (`1 - core_score`, `build_full_predictions.py:65-67`) lets every
downstream consumer sort descending with no per-class branching. → `MR:326-337`

*The measurement* (`test_run_tsg_two_hit.py`), on the 6 genes actually inverted in
production (`gene_role="tsg"` → `loss_of_function`):

```
mean hit@20, HIGH abundance      : 0.022050
mean hit@20, LOW  (production)   : 0.006933
delta +0.015117   ratio 3.18x   95% CI [+0.004117, +0.028008]  -> EXCLUDES 0
high beats low in 5/6 genes, 1 tie, 0 against
```

On all 14 TSG/both genes with ground truth the effect is the same: 3.01×,
CI [+0.0099, +0.0238], 13/14 genes. **Inverting makes ranking roughly 3× worse.**

*Why:* GDSC measures response to a drug **targeting** that gene, and those drugs need the
target present — ATM/ATR/CHEK2 inhibitors need the kinase; TP53's GDSC drugs are MDM2
inhibitors (nutlin) requiring functional wild-type p53. Inversion ranks lines that have
*lost* the gene, which cannot respond. The rule conflates two questions:

- "which lines have lost this suppressor's function?" — inversion correct
- "which lines respond to a drug targeting this gene?" — inversion backwards, **and this is
  what GDSC measures and what the pipeline is validated against**

*Caveats to state in the same breath:* n=6 (n=14 unrestricted); **not held-out** — the
curated 141-gene table has **zero** `loss_of_function` genes, so the seed-42 split contains
no TSGs; and any TSG with GDSC ground truth is nearly by definition a drug target, which is
plausibly why target-presence wins. Exploratory, not confirmatory.

*Related fact worth knowing:* only 6 of the 14 are inverted at all — `gene_role="both"`
maps to `activation_driven`, not `loss_of_function`.
→ `outputs/test_run_tsg_two_hit_results.json`

**C9. Did the Knudson two-hit rule help?**
No. Two-hit gating (`tsg_hit` AND (VAF≥0.75 OR multi-hit OR deletion)) scored 0.012194 vs
the production gate's 0.015631 — paired delta **−0.003437, CI [−0.010597, +0.004487]**,
includes 0. Mechanically the second hit came almost entirely from VAF (4,119 rows) and
multi-hit (2,919); the deletion route contributed only **56** rows, because just 103 TSG
(model, gene) pairs have both mutation and CNA data. Note `lohfraction` was deliberately
**not** used — it is a genome-wide per-cell-line scalar, and using it would repeat the
gene-agnostic-scalar error that got metabolomics and miRNA rejected.

**C5. Why the `1e6` multiplier in the sort key?**
`sort_key = flag·1e6 + core_score` collapses a two-column sort into one numeric key that
vectorises across all genes at once. `1e6` guarantees the boolean always dominates, since
`core_score ∈ [0,1]`. It's a performance trick, not modelling. → `MR:353-362`

**C6. Why sort genes before `rng.choice` in the train/test split?**
Python set iteration order isn't guaranteed stable across runs/versions, so sampling from
an unsorted set would silently give different test genes on different machines even with a
fixed seed. Sorting makes seed=42 genuinely reproducible. → `MR:344-352`

**C7. Why is hit@20's denominator "GDSC-sensitive ∩ scored", not all GDSC-sensitive?**
"The only honest denominator — we can only rank what we have a score for." Including
unscored lines would penalise the ranker for coverage gaps rather than ranking errors.
Be ready for the follow-up: this *does* make the number look better, which is why the
denominator is stated explicitly in the output JSON. → `MR:363-380`

**C8. [DANGER] Your hit rates are ~2%. Isn't that terrible?**
Contextualise, don't dodge: it's top-20 recall out of ~1,485 candidate cell lines, so
random baseline is ~1.3%. The reported flat baseline is 0.0192 and driver-routing gives
0.0212 (+10.5%). The honest statement is that the *lift* is small and measured on 28
held-out genes — this is a feasibility result, not a deployable predictor.
→ `outputs/stage4_eval_summary.json`

---

## D. Rejected experiments & known gaps (examiners love these)

**D1. Why did you reject Chronos as a third layer?**
Blending CRISPR essentiality into `core_score` collapsed BCL2's hit@20 to **0**. It was
built on a 50/50 cell-line train split to avoid circularity and still failed. Chronos
survives in the pipeline only as a **read-only per-gene validation check**
(ρ>0.1, p<0.05 → "validated"), which promotes unknown-class genes to moderate confidence.
Note the useful nuance: ρ_EC ≈ 0.005 means Chronos is almost *uncorrelated* with
abundance — it's measuring something real but different. → `MR:211-227`, `MR:494-511`

**D2. Why reject metabolomics and miRNA?**
Both were tested as a gene-agnostic global scalar blended 50/50 with core_score. Paired
bootstrap (2,000 resamples) on Δρ: **CI included 0** in both cases and the point estimate
moved the *wrong* way. The lineage-confound ANOVA came back small (metabolomics R²=4.2%,
miRNA R²=7.2%) — so the scalar isn't even a useful lineage proxy; it's closer to noise,
most likely from averaging across unrelated pathways. → `MR:534-608`, `DECISIONS.md`

**D3. Why a *paired* bootstrap rather than two independent ones?**
The same resampled rows feed both baseline and candidate, which cancels shared sampling
noise and gives a more powerful test of the *difference*. → `MR:534-556`

**D4. Is the signature discount applied, and does it matter?**
It is **not** applied — and it was **tested and measured null**, which is the stronger
answer. `signature_discount.parquet` (m_mut for 1,955 lines) is fully computed but never
wired in. `test_run_signature_discount.py` rebuilt `p_mutation` from source
(reproducing the shipped file to **0.000e+00** max difference, confirming the Hill
parameters), applied the discount, and re-ran the seed-42 held-out harness:

- hit@20 any-alteration: 0.025247 → 0.025356
- paired bootstrap Δ: **+0.000110, 95% CI [0.000000, +0.000329] — includes 0**
- activation_driven Δ-std: 0.014727 → 0.014735 (**variance did not shrink**)

The stated hypothesis — that hypermutation drove the any-alteration arm's variance — was
**falsified**. Mechanism: m_mut averages 0.950, so only 5,083 of 632,919 rows (0.8%) cross
the 0.5 threshold, and just 13 rows within the held-out genes. And by B8 it could never
have moved production ranking anyway.
→ `outputs/test_run_signature_discount_results.json`

**D4b. [DANGER] Your CNA thresholds assume a diploid genome. Are your cell lines diploid?**
**No — and this is a live confound in a production path.** `build_cna_layer.py:30-31,88-89`
uses absolute cutoffs (`TOTAL_CN > 2.5` = amplification, `< 1.5` = deletion) with no
normalisation by each line's baseline ploidy. But measured from `signatures_model_level.parquet`:

- median ploidy of scored lines = **3.06** (mean 3.10), not 2.0
- **53.0%** of lines have ploidy > 3
- WGD lines get **1.31×** the CNA-alteration rate of non-WGD lines
  (0.0277 vs 0.0212, Mann-Whitney **p = 1.4e-05**)
- ρ(mean TOTAL_CN, ploidy) = +0.109, p = 0.002

So the 2.5 "amplification" threshold sits *below* the median baseline copy number — for
much of the panel, "amplified" partly means "this line has a large genome." Since
`cna_flags` feeds confidence-tier upgrades (Stage 5), this propagates into user-facing
confidence labels.

**Follow-up you must be ready for — "did normalising by ploidy fix it?" Partly, and the
dissociation is the interesting part.** `test_run_ploidy_normalised_cna.py` applied
ploidy-relative cutoffs (ratio >1.25 amp, <0.75 del, chosen so diploid lines keep today's
calls exactly):

- ρ(alteration rate, ploidy): +0.090 (p=0.012) → **+0.045 (p=0.212, no longer significant)**
  — the continuous ploidy artifact **is** removed
- WGD vs non-WGD ratio: 1.43× → **1.40× (p=4.6e-07)** — the binary WGD gap **survives**

Interpretation: the arithmetic artifact was real and is now gone, but the WGD excess is
largely **genuine biology** — genome-doubled lines really do accumulate more focal
copy-number events. Correcting the threshold does not (and should not) erase that.

Blast radius: 120 of 117,540 pairs change (0.1%), but **93 of them currently carry "high"
confidence and would be downgraded** — roughly 1% of all high-confidence predictions,
retracting claims that rested on ploidy-inflated amplification calls. Ranking is
mathematically unchanged (CNA never enters the gate; asserted in the test, not assumed).
So the fix buys **precision, not performance**. → `outputs/test_run_ploidy_normalised_cna_results.json`

**D5. What's calibrated vs. placeholder in this pipeline?**
- *Measured/assumed:* ρ_EP=0.46 (literature), ρ_EC=ρ_PC=0.005 (measured), CN thresholds
  2.5/1.5 (standard diploid-based), ρ>0.1 & p<0.05 (chosen, applied consistently)
- *Placeholder:* every Track C Hill `(p0,k)`, the `hill` normalisation variant, driver
  boost 0.15
- *Tested & rejected:* Chronos 3rd layer, metabolomics scalar, miRNA scalar, enzyme-flux
  proxy, miRNA dampening
→ `MR:670-686`

**D6. Why is there no fuzzy string matching in cell-line resolution?**
There deliberately isn't a similarity *score* anywhere. Matching is exact-match after
deterministic normalisation, or substring containment. "Fuzzy" in this repo means
*permissive substring*, not Levenshtein. The RRID→CVCL→ACH resolver rescued 18 of 30
test names that naive substring matching got wrong. → `MR:44-54`, `MR:230-238`

**D7. Your AUROC is 0.52–0.56. Isn't the whole score just noise?**
No — that was **the wrong metric**, and proving so is the strongest result in the project.
AUROC concentrates its mass at the *top* of a ranking, but `core_score` is an abundance
prior: a necessary-not-sufficient condition. The claim it licenses is "low-abundance lines
are **depleted** of dependency", not "high-abundance lines are enriched". Split the ROC:

| label | gate region (FPR 0.8–1.0) | ranker region (FPR 0.0–0.2) |
|---|---|---|
| GDSC | 0.601 | 0.518 |
| Project Score | 0.593 | **0.501** |
| Chronos | 0.596 | **0.500** |

The gate region is above chance at `p ≈ 1e-60` to `1e-76`; the ranker region is *at* chance.
So the layer is a **filter, not a ranker**. It replicates across two independent CRISPR
screens (Broad 0.9950 vs Sanger 0.9960 decile-1 depletion, difference `p = 0.19`).
→ `MR:§C.15`

**D8. What's the honest caveat on that gate result?** *(volunteer this — it is the
strongest thing you can do with it)*
Part of the floor effect is a **profiling confound**, and the project's own adversarial
test found it. Stratifying on whether the line was profiled at all: at the floor, detected
2.45% vs not-detected 2.29% — ratio 1.07, **Fisher `p = 0.246`**, indistinguishable. In the
mid-range the same split gives ratio 1.62 at `p = 0`. So the gate is real in the mid-range,
but at the extreme floor "not expressed" is not separable from "never measured". Also: the
measured floor rate (1.39%) sits *at or below* the label's own noise floor (2.56%, from
TAS2R\* control genes), so depletion is complete to within what the label can resolve.
→ `MR:§C.15`, open gap 13

**D9. Then why keep a continuous percentile? Why not just a detection boolean?**
Because it was tested and the boolean loses. Among cell lines an absolute per-line detection
call rates **expressed**, the bottom decile of the percentile still carries a dependency
ratio of **0.810–0.827** — spread only 0.017 across six independent calibration rules
(four quantiles of a curated negative-control reference, a mode-constrained Hart zFPKM
variant, and TPM≥1), clearing a 10% effect floor in **6/6**. On paralog-singleton genes
only: 0.790–0.799, still 6/6. The percentile carries graded information beyond "on or off".

Two structural points worth having ready:
- Within a gene the percentile is a **monotone transform of raw log2TPM**, so
  `pAUC(percentile) ≡ pAUC(raw)` *exactly*. The only meaningful comparison is against a
  score on the other axis (per-sample, not per-gene).
- The detection rule is **inert for most genes** — it excludes nothing at all for 393–591
  of 778 genes — so it could not replace the percentile wholesale even if it matched it.
→ `MR:§2.1`, `MR:§C.16`

**D10. Didn't the absolute per-line score beat the percentile at one point?**
It appeared to, by +0.030 pAUC, and that number was **selection optimism** — a good example
of catching your own error. The per-line score's pAUC depends on which quantile it is
centred on, and the paired delta against the percentile *changes sign* across the sweep
(+0.0134 at Q0.90, −0.0301 at Q0.99). Choosing the best centre and quoting its pAUC on the
same labels is selection on the test set. Under 20 random half-splits — centre chosen on
half A, evaluated on half B, percentile recomputed within half B because it is
panel-dependent:

```
OPTIMISM (selection on A − held-out on B)  +0.0302
held-out delta (percentile − per-line)     −0.0050
per-line wins on held-out data              10/20 reps
```

The optimism `+0.0302` **is** the claimed advantage `−0.0301`. On held-out lines it is a
coin flip. Bonus finding: the percentile is the *lower-variance* estimator across panel
composition (held-out range 0.078 vs 0.170) — notable because panel-dependence was the
objection that motivated the absolute alternative in the first place. → `MR:§C.16`

**D11. Did you control for paralog buffering?**
**No, and I can show you why not** — which is the honest answer. Paralog buffering depresses
single-knockout dependency, so it is a real candidate confound. The repo has no paralog
annotation, so HGNC gene-group co-membership was tried as a proxy and **failed its own
validity check**: buffering predicts *lower* prevalence among paralogous genes, but measured
prevalence was 0.1072 with paralogs vs 0.0629 without, `p = 0.7766` in the wrong direction.
Reason: HGNC groups skew toward well-characterised families (ribosomal, proteasome,
spliceosome) that are pan-essential, and that swamps buffering. The test refused to
interpret the stratification, by a guard written into it beforehand. The finding *survives*
restriction to singleton genes, but that is not the same as controlling for paralogy. A real
test needs Ensembl Compara % identity or the DepMap/Dede/Ito paralog sets. → `MR:§C.17`,
open gap 16

**D12. Is the two-layer score better than RNA alone?**
No — and this contradicts a rule still live in the code. RNA alone beats every Stouffer
variant (GDSC AUROC 0.563 vs 0.556 vs 0.547 centred; Chronos 0.519 vs 0.517 vs 0.509), and
stratum centring makes it worse, not better (`centring_gain` −0.0145 at `p = 1.5e-17`).
The 2-layer stratum does sit higher (median shift 0.210, positive in 86.3% of 1,562 genes)
but correcting for that does not recover anything. Consequence:
`n_layers == 2 → "high"` at `build_full_predictions.py:151` promotes rows on a basis the
project's own test rejects. **Not yet fixed in code** — say so. → `MR:§C.18`, open gap 14

**D13. What bounds any ranking claim you could make?**
Label reliability. Broad and Sanger — two independent CRISPR screens of the same biology —
agree with **each other** at median Spearman `ρ = 0.172`, with only 1.9% of genes above 0.5
and 5.3% negative. No score can exceed the reliability of the label it is scored against,
so a ranker AUROC near 0.5 is partly a statement about the ground truth, not only about
`core_score`. Every AUROC in this project should be read against that ceiling.
→ `MR:§C.15`, open gap 18

---

## E. Fast-recall table

| Number | Value | Where |
|---|---|---|
| Total scored rows | **29,781,274** | `core_score.parquet` (re-measured) |
| Genes / cell lines | **19,177 / 1,746** | distinct in the score grid (re-measured) |
| Lines per gene | median **1,483** (min 5, max 1,746) | bimodal at 1,479 / 1,746 |
| 2-layer rows | **17.72%** (5,277,064 / 769 lines) | post-ProCan (re-measured) |
| 1-layer rows | 82.28% (24,504,210) | the score is single-omics 4 rows in 5 |
| Protein coverage | CCLE 375 + ProCan 952 lines | `coverage_matrix_enriched` |
| 2-layer >0.95 (post-fix) | 2.32% (was 22.5%) | `MR:165-167` |
| ρ_EP / ρ_EC | 0.46 / 0.005 | `MR:159-161` |
| Curated genes / held-out | 141 / 28 (seed 42) | `stage4_eval_summary.json` |
| hit@20 flat / driver / any | 0.0192 / 0.0212 / 0.0252 | same |
| Driver lift / variance reduction | +10.5% / 44.6% | same |
| Confidence: unknown / high | 92.5% / 0.03% | `predictions_with_confidence.parquet` |
| CN amp / del thresholds | >2.5 / <1.5 | `MR:262-267` |
| Chronos validation thresholds | ρ>0.1, p<0.05, n≥30 | `MR:494-501` |
| Lineage R²: metabolomics / miRNA | 4.2% / 7.2% | `MR:599-602` |

**Validation round §C.15–C.18** *(the numbers to have on instant recall)*

| Number | Value | Where |
|---|---|---|
| Gate vs ranker pAUC (Chronos) | **0.596 / 0.500** | `MR:§C.15` |
| Gate replication Broad / Sanger | 0.9950 / 0.9960, diff `p=0.19` | `MR:§C.15` |
| **Label–label agreement ceiling** | median **ρ = 0.172** (1.9% > 0.5) | `MR:§C.15` |
| Profiling confound at the floor | ratio 1.07, **`p = 0.246`** | `MR:§C.15` |
| Label noise floor vs measured floor | 2.56% vs 1.39% | `MR:§C.15` |
| Percentile increment over detection | **0.810–0.827**, spread 0.017, 6/6 | `MR:§C.16` |
| — singletons only | 0.790–0.799, 6/6 | `MR:§C.17` |
| Half-sample selection optimism | **+0.0302** (= the claimed advantage) | `MR:§C.16` |
| Held-out delta / reps won | −0.0050 / **10 of 20** | `MR:§C.16` |
| Held-out variance (pct / per-line) | range 0.078 / 0.170 | `MR:§C.16` |
| Unexpressed reference separation | 7.690 log₂ units (438 genes) | `MR:§C.16` |
| Paralog proxy validation | `p = 0.7766`, **wrong direction** | `MR:§C.17` |
| RNA alone vs Stouffer vs centred | 0.563 / 0.556 / 0.547 (GDSC) | `MR:§C.18` |
| CNA deleted dep ratio | 0.380 (base 10.44%) | `MR:§C.18` |
| CNA usable gene count | 69 genes with ≥10 deleted lines | `MR:§C.18` |
| Lineage conditioning / conformal | Δ=0.0 `p=0.54` / ratio 0.93–0.99 | `MR:§C.15` |
