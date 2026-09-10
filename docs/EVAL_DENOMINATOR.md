# EVAL_DENOMINATOR — the fix, and what it changes

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

**Status: fixed.** `src/pipeline/build_eval_tested_denominator.py` (additive; no existing
pipeline file modified). Outputs `eval_tested_denominator.parquet` / `.json`.

**This corrects a claim I made in an earlier report.** See §4.

---

## 1. The defect

Three sites left-join GDSC2 sensitivity onto every scored pair and then collapse the missing
values:

| Site | Code |
|---|---|
| `src/pipeline/stage4_eval_save.py:41` | `core_t["is_sensitive"].fillna(False)` |
| `src/pipeline/06_held_out_eval.ipynb` cell 4 | same |
| `src/pipeline/build_gene_regime.py:150` | `is_sensitive.astype("boolean").fillna(False)` |

Every (gene, line) pair GDSC2 never tested therefore entered the denominator as a **confirmed
negative**. This is the three-state failure (`BUILD_SPEC.md` C3, item 6) applied to the label
column itself.

### Scale — `VERIFIED`

| | Broken denominator | **Fixed (tested pairs only)** |
|---|---|---|
| median lines per gene (N) | 1,563 | **702** |
| median sensitive per gene (K) | 228 | 228 |
| **median prevalence** | **0.135** | **0.295** |

**About 861 untested lines per gene were being scored as true negatives** — 55% of each
gene's denominator. Prevalence was understated by a factor of 2.2.

## 2. The fix

Inner-join, not left-join. A pair with no GDSC2 row is `not_tested` and is excluded from
**both** numerator and denominator. It is not a negative.

The split, seed and ranking are otherwise unchanged, so the before/after is paired: same 141
curated genes, same 28 held out, same `core_score` ordering.

## 3. What it changes

Held-out genes (n = 28), paired bootstrap over genes, 10,000 resamples.
`*` = 95% CI on (observed − chance) excludes zero.

| Metric | Denominator | Observed | Chance | Excess | 95% CI |
|---|---|---|---|---|---|
| recall@1 | broken | 0.0013 | 0.0006 | +0.0007 | [−0.0005, +0.0026] |
| recall@1 | **fixed** | 0.0027 | 0.0017 | +0.0010 | [−0.0003, +0.0028] |
| recall@5 | broken | 0.0042 | 0.0031 | +0.0010 | [−0.0011, +0.0037] |
| recall@5 | **fixed** | 0.0101 | 0.0083 | +0.0019 | [−0.0009, +0.0050] |
| **recall@20** | broken | 0.0140 | 0.0126 | +0.0014 | [−0.0027, +0.0062] |
| **recall@20** | **fixed** | **0.0420** | **0.0331** | **+0.0090** | **[+0.0022, +0.0163]** `*` |
| anyhit@1 | fixed | 0.4643 | 0.3372 | +0.1271 | [−0.0475, +0.3128] |
| anyhit@5 | fixed | 0.8214 | 0.8437 | −0.0222 | [−0.1747, +0.1130] |
| anyhit@20 | fixed | 1.0000 | **0.9983** | +0.0017 | [+0.0010, +0.0025] `*` |

Tuning genes (n = 113), fixed denominator: recall@5 +0.0041 `*`, recall@20 +0.0120 `*`,
anyhit@1 +0.1134 `*`, anyhit@5 +0.0615 `*`.

## 4. The correction to what I told you earlier

In `STATE_REPORT.md` B1.4 I reported that on the 28 held-out genes, **every** metric at every
k had a 95% CI on excess-over-chance that included zero, and I described that as the strongest
available statement of §9 gap 15.

**That was measured under the broken denominator. Under the corrected one it is no longer
true.** Held-out `recall@20` is above chance: observed **4.20%** against a chance baseline of
**3.31%**, excess **+0.90 percentage points**, CI **[+0.22, +1.63]**.

So the corrected finding is: **the score-based ordering carries a small but real above-chance
signal at k = 20 on held-out genes.** The broken denominator was masking it — diluting both
the observed rate and the baseline with 861 untested lines per gene, which compressed the
difference toward zero.

### What this does and does not license

- **Does not** rescue any previously quoted figure. The `.docx`'s "top-20 93.6%" is the
  any-hit metric, whose corrected chance baseline is **99.83%** — even more clearly at chance
  than under the broken denominator (90.8%). Observed any-hit@20 is 100.00%. The metric is
  useless at a 29.5% prevalence and should not be quoted at all.
- **Does not** revive the confirmation count. `RANKING_FEASIBILITY.md` returned ABANDON on
  coverage-confound grounds, which this does not touch.
- **Does** mean the retained set has some internal order. The effect is small: 4.20% vs 3.31%,
  a relative lift of about 27% on a base rate of 3.3%, on 28 genes. Quotable only with its
  interval.
- **Does** create a tension worth resolving. Ranker-region pAUC measured **0.5002** against
  CRISPR dependency labels; this measures above-chance recall@20 against GDSC2 drug
  sensitivity. Different labels, different answers. Neither is wrong; they are different
  questions, and the write-up must say which label set any ordering claim rests on.

## 5. What is unaffected

**The gate result does not pass through this denominator.** Gate-region pAUC 0.596–0.608
(p = 2.8e-76) is measured against CRISPR dependency labels — Chronos / Project Score — with
prevalence and negatives defined by the screen, not by GDSC2 coverage. It is unchanged, as is
the label-FPR ceiling of 0.0093 derived from DepMap's own reference non-essentials.

## 6. Required follow-on

1. **Repoint the three defect sites.** They still ship the broken denominator; this build is
   additive and parallel. Repointing changes `gene_regime.parquet`'s
   `abundance_tracking` / `activation_driven` assignment, since `build_gene_regime.py:150`
   feeds the `hit_at_20_train` cut at `LIFT_DEFAULT = 2.0`. That reclassification must be
   re-derived, not carried over.
2. **Re-run every @k figure in the docs** against the fixed denominator, each with its
   hypergeometric baseline and a paired interval.
3. **Retire any-hit@k entirely.** At 29.5% prevalence its chance baseline is 99.83% at k = 20
   and 84.4% at k = 5. It cannot distinguish anything.
4. **Add the assertion.** `BUILD_SPEC.md` C3 A13 — no @k metric may be emitted without its
   chance baseline in the same record — plus a new one: no label column may be `fillna`'d into
   a negative.
