# TRANSCRIPTOMICS_DECISIONS

Questions only you can answer, with options and costs. Nothing marked *pending decision* in
`TRANSCRIPTOMICS_SPEC.md` should be implemented until the relevant question here is closed.

Ordered by consequence. Q3 is the one that decides whether the layer ships at all.

---

## Q1. Split GEO by processing method (or series)? — **recommendation: no**

**Why it's yours.** The flow diagram commits to a `source × method × lineage` stratification.
Abandoning it is a design reversal, and the diagram may reflect an intention I cannot see.

**What the measurements say.**

- **There is nothing to split on.** No processing-method column exists anywhere. All 2,329
  GSMs with platform metadata are `gpl570`, and all 19 series are on the same linear,
  MAS5-like scale (inferred from value distributions — medians 39–171, maxima 3.7k–12k;
  RMA would be log2-scaled ≈2–14). One platform, one processing family.
- **Splitting by series would amplify the existing fault, not correct it.** 276 of 590 GEO
  lines (46.8%) appear in ≥2 series, 158 in ≥3, one in **9**. Series-specific z-scores
  correlate at median ρ = **0.448**.
- **Quantified:** 3-source combination overstates the combined Z by **50.3%**. Splitting GEO
  into three series takes that to **76.5%** — worse by **26.2 percentage points** — while
  buying 0.28 of an effective independent source (1.61 of 5, against 1.33 of 3).
- **7 of 19 series fall below `MIN_PEERS_FOR_LINEAGE = 15`**; two (n=2, n=3) are at or below
  `MIN_PEERS_FOR_Z = 3` and cannot support a spread estimate at all.
- **The notebook already handles this correctly.** `fetch_gene_expression` collapses GEO to
  `median(...) GROUP BY model_id` with `count(*) AS n_samples` as a √n precision weight —
  the notebook's own stated replicate principle, correctly applied.

| Option | Cost |
|---|---|
| **A. Do not split (recommended)** | Keep the current collapse. Costs nothing; it is already implemented and already right. Requires correcting the flow diagram |
| B. Split by series | +26.2 pp false confidence, 7 unusable strata, and abandons the notebook's own replicate principle |
| C. Split by method | Not implementable — no method data exists |
| D. Download the 18 GSE `!Sample_data_processing` fields from NCBI first | Small metadata fetch. Would confirm the submitters' wording, but the numeric signature is stronger evidence and **the answer changes no decision** — one method either way |

**Recommendation: A.** If you want D for the record, say so and I will fetch it — I did not,
because it is an external download that changes nothing.

**One correction either way:** the notebook (cells 4 and 37) describes GEO as "much of it
RMA-normalised". The data says MAS5-like. The floor conclusion survives — GEO has a hard
minimum at 3.03 and **zero** values below the detection threshold — but the stated mechanism
should be corrected.

---

## Q2. Rename the shrinkage, or estimate it properly?

**Why it's yours.** One option is a five-minute rename; the other opens a real piece of
statistical work that also closes an existing gap.

**The problem.** Cell 0 describes step 5 as *"empirical-Bayes shrinkage"* and cites Efron &
Morris (1975). The defining property of that family is that the shrinkage factor is
**estimated from the data**. Here it is `z * k/(k+1)` — a function of the source count only.
It never looks at the data. With `PRIOR_STRENGTH = 1.0`, one source keeps ½ of its z and
three keep ¾, regardless of whether the sources agree or how much signal exists.

The notebook's own limitations section is honest (*"a judgement, not a derivation"*). The
headline description and the citation are not.

| Option | Cost |
|---|---|
| **A. Rename honestly (recommended, cheap)** | Call it `source_count_discount` with `DISCOUNT_STRENGTH`; drop the Efron & Morris citation. Behaviour unchanged, claim becomes true. ~5 minutes |
| **B. Estimate it properly** | Real James–Stein: `1 − (k−2)·σ²_within / Σ(z_i − z̄)²`. Requires a per-layer variance estimate the project does not have — **which is open gap 5 in `MATH_REFERENCE.md §9`. One estimate closes both.** Days, not minutes |
| C. Leave as-is | Not recommended. A reviewer who knows the result will check whether the shrinkage adapts, find it does not, and reasonably wonder what else was described more favourably than it behaves |

**Recommendation: A now, B if the layer ships as a scorer** (which Q3 may decide against).

**If it stays a chosen value, standing rule 4 applies:** sweep `PRIOR_STRENGTH ∈
{0, 0.5, 1.0, 2.0, 5.0}` with held-out selection. **`0` must be in the grid** as the null —
if no shrinkage wins, drop the step.

---

## Q3. Is the layer scoped to dependency prediction, or to unusual-expression detection? — **the decisive question**

**Why it's yours.** This decides what the layer is for, and therefore whether the scoring
machinery ships at all.

**What the decisive experiment found (F1).** The guarded scorer was run through the same
harness the main pipeline used for the plain version, with real DepMap Chronos labels
(scale check passed: essentials −0.9733, non-essentials +0.0012, label FPR 0.0093):

| Comparison | Median Δ pAUC | 95% CI (paired bootstrap over genes) | Better in |
|---|---|---|---|
| **guarded − lineage** | **+0.0022** | **[−0.0072, +0.0136]** | 51.4% |
| guarded − panel | −0.0004 | [−0.0241, +0.0098] | 48.3% |
| lineage − panel | −0.0120 | [−0.0256, +0.0053] | 44.2% |

**Negative.** The guards do not recover an effect that plain conditioning does not. All
three arms sit at pAUC ≈ 0.51 — near chance. The `lineage − panel` arm independently
reproduces §C.15's neutral result, so the harness is behaving.

**The important caveat, which is why this is a decision and not a verdict.** The guards
retain a median of **100.0%** of lines per gene and drop all lines for only 22 of 438 genes.
**A guard that does not bind cannot change a ranking metric.** The harness measures ranking
quality among scored lines; the guards produce *abstention*. Gate-region pAUC cannot price
abstention, and the genes where the guards bind hardest are the ones excluded from the
harness entirely.

So F1 says the layer is not a dependency predictor. It does **not** say the guards are
worthless — F2 measures where their value actually is.

| Option | Cost |
|---|---|
| **A. Rescope to unusual-expression detection (recommended)** | Honest given F1. The layer answers "is this line unusual for this gene within its lineage", which is a real and useful question — but it is not dependency prediction and must not be described as such. Loses the headline claim |
| B. Keep it scoped to dependency prediction | Not supportable. Requires quoting a pAUC of 0.51 with a CI spanning zero |
| **C. Ship only the silence guard, drop the scorer (see Q4)** | Smallest surface, largest measured value. Loses the stratification work, which is genuinely good |
| D. Re-run F1 on a gene set where the guards actually bind | The right follow-up if you want to keep B alive. Select genes with low `frac_expressed` and re-test. **Note this selects on the guard's own criterion**, so it needs a held-out design to avoid the B3 error |

**Recommendation: A, plus C.** The scoring method is defensible as an unusual-expression
detector with corrected operating characteristics; it is not defensible as a dependency
predictor. D is worth doing before closing the question permanently.

**Whichever you choose, cell 37's reporting template must not be used as written.** It
quotes a 1.2% false-positive rate and 100% recovery at 3 MAD. Both are superseded: the true
operating rate is **α = 0.116** (B1) and recovery at 3 MAD is **85%** (B2).

---

## Q4. Adopt the silence guard as a validity check on the main pipeline?

**Why it's yours.** It would withhold roughly 30% of the main pipeline's current output.
That is a large product change, not a technical detail.

**What the measurement says (F2).** Applied to the main pipeline's own pooled DepMap layer,
19,173 protein-coding genes × 1,479 lines:

| | Value |
|---|---|
| **genes failing the guard at `SILENT_FRAC = 0.20`** | **5,847 of 19,173 (30.5%)** |
| **gene × line combinations affected** | **8,647,713 of 28,356,867 (30.5%)** |

**And a spread floor does not catch them.** Of those 5,847 genes, the pipeline's own
`signal_spread` calls **2,420 "wide"**, and **2,091 exceed `FLAT_FOLD_CEILING = 3.0`**. A
dispersion-based check does not merely miss these — it actively certifies them as
well-dispersed. This confirms the notebook's decisive observation (EGFR in blood: MAD 0.178,
above any sensible floor, 38% of values at exactly zero) across the whole gene universe.

This bears directly on **§9 gap 13**, the main pipeline's highest-priority open item, and it
maps onto abstention states the query layer already defines (`UNINFORMATIVE` in
`docs/QUERY_LAYER_STATES.md`) — no new output vocabulary needed.

| Option | Cost |
|---|---|
| **A. Adopt at `SILENT_FRAC = 0.20` (recommended)** | Withholds 30.5% of output as `UNINFORMATIVE`. Large visible change; the withheld scores were not meaningful, so the loss is nominal rather than real |
| B. Adopt at 0.05 | Withholds 18.6% (3,565 genes). More conservative, catches less |
| C. Adopt as a warning label, not a withholding | Keeps output volume. Weaker — a caveat next to a number is usually read as a number |
| D. Do not adopt | Leaves 8.6M scores on degenerate distributions with no guard, and gap 13 untouched |

**Recommendation: A, after sweeping the cut.** Standing rule 4: sweep
`SILENT_FRAC ∈ {0.05, 0.10, 0.20, 0.30}` with held-out selection before committing.
Measured flag counts at each: 3,565 / 4,607 / 5,847 / 6,623 genes. **Do not inherit 0.20
unexamined from a different codebase** — it was chosen there for a different stratum
(lineage) than the one it would be applied to here (the pooled panel).

---

## Q5. Which combination route becomes primary?

**Why it's yours.** It changes every number the layer reports, and one option produces zero
significant calls.

**State.** Three faults compound in the current chain:

1. **Wrong quantity.** Stouffer combines standardised statistics; the notebook feeds it
   `robust_z`, an effect size with measured **sd = 1.606**, skew −4.54, kurtosis 78.1, and
   **P(|z| > 1.96) = 0.1055** against 0.05 predicted. The p-value read off the normal CDF is
   not a p-value.
2. **Independence assumed.** Correlation inflation of **36.5%** (two sources) to **50.3%**
   (three).
3. **The valid route is computed and discarded.** `empirical_p` is calculated per
   (line, source), used only to build `fisher_p`, which is itself never read.

**Effect of switching to the empirical-p chain (TSPAN6, 1,580 lines):** 232 lines (14.7%)
change verdict, **uniformly downward**, and significant calls go from **115 → 0**.

| Option | Cost |
|---|---|
| **A. Empirical-p → standardised → Strube/Hartung with measured R (recommended)** | Valid end to end. Produces **zero** significant calls for TSPAN6 — which is the measurement of how much of the current output was an artefact, not a reason to prefer the invalid route |
| B. Keep robust z, add only the R correction | Fixes fault 2, leaves fault 1. Still not a p-value |
| C. Fisher on empirical p | Assumption-light, but **also assumes independence** — needs Brown's method (1975), the Fisher analogue of Strube/Hartung. Does not escape Part A |
| D. Report effect sizes only, abstain on significance | Radical but coherent, and consistent with the resolution floor below |

**Recommendation: A, with D as the honest fallback** where the resolution floor makes
significance unreachable.

**The resolution floor must be stated in the methods either way.** `empirical_p` has a
minimum attainable value of `2/(n+1)`. At `MIN_PEERS_FOR_LINEAGE = 15` that is **0.125 —
2.5× the 0.05 threshold, before any BH correction.** The smallest lineage in which any line
can reach p < 0.05 at all is **n = 39**. **14 of 31 lineages are below that, containing 270
of 1,479 cell lines (18.3%)** — no line in them can ever be called significant, at any
effect size. This is a real limit on the design and belongs in the methods, not in a
workaround. Do not route those lineages back through the parametric path; that is the path
fault 1 invalidates.

---

## Q6. Should the notebook be committed to the repository?

**Why it's yours.** It is a working-practice question, but it has blocked this audit once
already.

**State.** `02_transcriptonomics.ipynb` is **not in the repository**. The previous audit
recorded it as `UNVERIFIED — could not locate` and the C5 decisive experiment could not be
run. It was supplied by hand for this task, from `C:\Users\vigne\Downloads\`.

Consequences that have already materialised:

- The flow diagram describes a `source × method × lineage` version with a
  `harvest_geo_methods` step. **The supplied notebook has neither.** Two versions are in
  circulation and only one is measurable.
- The notebook opens `celllineselector.db` by bare relative path (`DB_PATH =
  "celllineselector.db"`), which resolves only if run from `src/pipeline/outputs/`. Nothing
  records that.
- Nothing imports it and it writes no artefact, so no build would have caught its absence.

| Option | Cost |
|---|---|
| **A. Commit it, with the DB path fixed and the version question resolved (recommended)** | Small. Makes the work reviewable and stops the diagram and the code diverging |
| B. Keep it out until Q3 is decided | Defensible if the layer may be dropped — but F2's guard is worth adopting regardless of Q3, so something is being committed either way |
| C. Commit only the silence guard as a pipeline utility | Narrowest. Loses the stratification work, which is the layer's genuine contribution even though F1 was negative |

**Recommendation: A.** Also resolve which version is current — if a
`source × method × lineage` implementation exists somewhere, it needs measuring too, and Q1
would need re-answering against it (though on the A3/A4 numbers the answer would not change).
