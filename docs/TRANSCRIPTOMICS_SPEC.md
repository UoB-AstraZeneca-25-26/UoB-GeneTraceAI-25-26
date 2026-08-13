# TRANSCRIPTOMICS_SPEC — corrected design, parameter register, source roles, scope

Parts E and F of the transcriptomics brief. Tags: `PROVEN` (published) · `VERIFIED`
(computed here) · `UNVALIDATED` (asserted, no measurement) · `NOVEL FORMALISATION`.

---

# PART E — The free parameters

Every constant governing this design, with value and location. All are in
`02_transcriptonomics.ipynb` cell 1 unless stated.

## E1. The three groups

### Group 1 — published conventions that should be cited

| Constant | Value | Location | Citation | Cited in notebook? |
|---|---|---|---|---|
| `EXPRESSED_MIN` | 1.0 | cell 1 | HPA / standard RNA-seq detection convention — **derivation below** | as "HPA-style detection floor", derivation absent |
| `I2_CUT` | 50 | cell 1 | Higgins & Thompson (2002), *Stat Med* 21:1539 | yes, cell 0 |
| `HPA_FOLD` | 5.0 | cell 1 | HPA's 5-fold enriched/enhanced rule | as "HPA's 5-fold rule" |
| `DELTA_CUT` | 0.5 | cell 1 | Cliff's delta, "medium or larger" | as "medium or larger" |
| `TAU_BROAD` / `TAU_SPECIFIC` | 0.30 / 0.60 | cell 1 | Yanai et al. τ bands | as "Yanai conventions", no reference |
| BH-FDR at q<0.05 | 0.05 | cell 8 | Benjamini & Hochberg (1995), *JRSS-B* 57:289 | yes, cell 0 |
| MAD scaling 1.4826 | implicit | `stats.median_abs_deviation(scale="normal")` | Leys et al. (2013) | yes, cell 0 |

`I2_CUT = 50` carries a caveat that the citation does not cover: Higgins & Thompson's
convention was not intended for k = 2–3, and 52.1% of the flags it raises are at k = 2.
See `TRANSCRIPTOMICS_CORRELATION.md` Part D.

### Group 2 — values with a written derivation in the notebook

| Constant | Value | Location | Derivation given |
|---|---|---|---|
| `MIN_PEERS_FOR_Z` | 3 | cell 1 | yes — "below this a lineage has no distribution at all"; MAD is 0 by construction, e.g. `adrenal_cortex` n=1 |
| `SILENT_FRAC` | 0.20 | cell 1 | yes — cell 7 gives the EGFR-in-blood case (n=104, median 0.12, 38% at exactly zero, MAD 0.178) and argues a MAD floor cannot catch it |
| `SILENT_SOURCES` | `[depmap_expr, hpa_rna]` | cell 1 | yes — and, better, it is *derived at runtime* by `source_floor_check()` rather than trusted |
| `TAU_SOURCES` | `[depmap_expr, hpa_rna]` | cell 1 | yes — floored source compresses the ratio τ depends on; GEO τ = 0.666 vs 0.748/0.741 |
| `NO_LINEAGE` | `"(no lineage recorded)"` | cell 1 | yes — must not collide with DepMap's literal `"unknown"` label |

`SILENT_FRAC = 0.20` has a *rationale* but not a *derivation* — the EGFR case shows a guard
is needed, not that 0.20 is the right cut. It belongs in Group 3 for sweeping purposes; it
is listed here because the argument for its existence is written down, which is more than
most.

### Group 3 — chosen without either

| Constant | Value | Location | Consequence of the choice |
|---|---|---|---|
| `MIN_PEERS_FOR_LINEAGE` | 15 | cell 1 **and** cell 8 (redeclared) | flags lines as `small_lineage`; **duplicated in two places and could drift** |
| `PRIOR_STRENGTH` | 1.0 | cell 1 | sets shrinkage; see E3 |
| `MAD_FLOOR` | 0.1 | cell 1 | caps z where spread is small but non-zero |
| `FLOOR_FRAC_LOW` | 0.02 | cell 1 | the floor-detection trigger |
| `replicate_warn` | 1.0 | cell 8 signature | log2 units above which replicates are flagged as possible misidentification |
| verdict cut `abs(z_shrunk) > 1.5` | 1.5 | cell 8 `verdict()` | the strong/weak boundary — **undeclared as a constant at all**, hardcoded in a function body |
| verdict cut `z_shrunk > 1.0` | 1.0 | cell 8 `verdict()` | the weak/typical boundary — same |
| `n_perm` | 200 / 50 | cell 12 / cell 27 | validation resolution |
| `n_planted` | 20 | cell 28 | recovery test size |

**The two verdict cuts (1.5 and 1.0) are the most exposed.** They decide every user-facing
label, they are buried inside `verdict()` rather than declared with the others, and they
are on the scale of a statistic that C1 shows has sd 1.606 rather than 1.

## E2. `EXPRESSED_MIN` — the derivation, written out

The notebook calls it "HPA-style detection floor" and stops there. It looks arbitrary. It
is not, and the derivation is exact.

All three sources are placed on a `log2(x + 1)` scale before scoring — GEO and HPA by
`np.log2(d["value"].clip(lower=0) + 1)` in `fetch_gene_expression`, DepMap because it
ships as `log2(TPM+1)` already. On that scale:

```
value > EXPRESSED_MIN
⇔  log2(TPM + 1) > 1
⇔  TPM + 1 > 2
⇔  TPM > 1
```

**`EXPRESSED_MIN = 1.0` in log2(x+1) space is exactly the TPM ≥ 1 detection convention** —
the threshold the Human Protein Atlas uses to call a gene detected in a cell line, and the
most widely used "expressed" cut in RNA-seq generally (`PROVEN` as a convention).

The pleasing part is that the transform is what makes 1.0 the right number: pick any other
log base or offset and the correspondence breaks. Write this into the methods — as it
stands a reader sees a bare `1.0` and assumes it was tuned.

**Note the threshold is applied to HPA `nTPM`, not TPM.** nTPM is HPA's cross-sample
normalised value, so the correspondence is to HPA's own detection convention on its own
scale, which is the intended reading. State that explicitly too.

## E3. `PRIOR_STRENGTH` — it is not empirical Bayes

Cell 0 lists step 5 as *"empirical-Bayes shrinkage `z·k/(k+1)`"* and cites **Efron & Morris
(1975)**, *JASA* 70:311.

**This is not an empirical-Bayes estimator, and the citation does not support it.** The
defining property of the James–Stein / empirical-Bayes family is that the shrinkage factor
is **estimated from the data** — from the ratio of between-group to within-group variance,
so that the amount of shrinkage adapts to how much real signal is present. Here:

```python
def shrink(z, k, prior_strength=1.0):
    return z * k/(k+prior_strength)
```

The factor depends only on `k`, the number of sources. It is a **fixed, hand-set discount
for having few sources**. It never looks at the data. With `prior_strength = 1.0` a
single-source line keeps 1/2 of its z and a three-source line keeps 3/4, regardless of
whether those sources agree, how variable the gene is, or how much signal exists in the
panel. The notebook's own limitations section is honest about this — *"shrinkage strength
is a judgement, not a derivation"* — but the headline description and the citation are not.

### Recommendation — rename it, or estimate it properly

**Option A (recommended, cheap): rename honestly.** Call it a `source_count_discount` with
`DISCOUNT_STRENGTH = 1.0`, drop the Efron & Morris citation, and describe it as "a fixed
penalty on evidence seen in few sources". Nothing about the behaviour changes; the claim
becomes true.

**Option B (better, more work): estimate it.** The proper empirical-Bayes form shrinks
toward the grand mean by a factor estimated from the variance decomposition:

```
shrinkage = 1 - (k-2)·σ²_within / Σ(z_i - z̄)²        [James-Stein, k ≥ 3]
```

This requires a per-layer variance estimate, which the project does not have — it is
open gap 5 in `MATH_REFERENCE.md §9`. **One estimate would close both this and that gap.**

**Do not keep the current name.** A reviewer who knows the Efron & Morris result will check
whether the shrinkage adapts, find that it does not, and reasonably wonder what else was
described more favourably than it behaves.

### If it stays a chosen value — the required sweep

Standing rule 4. `PRIOR_STRENGTH ∈ {0, 0.5, 1.0, 2.0, 5.0}`, selected on **held-out genes**
against the F1 harness (real Chronos labels, gate-region pAUC, paired bootstrap over
genes). Report the retained fraction at each. `prior_strength = 0` (no shrinkage) must be
in the grid as the null hypothesis — if it wins, the step should be dropped.

The same sweep discipline applies to every Group 3 constant:

| Constant | Grid | Selection |
|---|---|---|
| `MIN_PEERS_FOR_LINEAGE` | {10, 15, 20, 30, 40} | held-out; report coverage lost at each. **Note C1c: below 39 no line can reach p<0.05 at all**, so this interacts with the resolution floor |
| `SILENT_FRAC` | {0.05, 0.10, 0.20, 0.30} | held-out; measured genes flagged at each: 3,565 / 4,607 / 5,847 / 6,623 of 19,173 (F2) |
| `MAD_FLOOR` | {0, 0.05, 0.1, 0.25, 0.5} | held-out; report how often it binds |
| `FLOOR_FRAC_LOW` | {0.01, 0.02, 0.05} | report which sources are called floored at each |
| verdict cuts | {1.0, 1.5, 2.0} × {0.5, 1.0, 1.5} | held-out; these must be re-derived anyway once C1 changes the statistic's scale |
| `I2_CUT` | {25, 50, 75} | held-out; report fire rate at each |

---

# PART F — Integration with the main pipeline

## F1. The decisive experiment — **the guards do not recover an effect**

### What was run

The notebook's **guarded** lineage scorer was run through the same harness the main
pipeline used for the plain version — `src/pipeline/test_run_improvement_scan.py` PART 1,
with real DepMap Chronos labels from `GeneFitnessEffect_Chronos_Achilles.hdf5`.

Scale check passed before scoring (the check whose absence let a negated Project Score
matrix pass as Chronos for the whole project): reference essentials median **−0.9733**
(expect ≈ −1), reference non-essentials **+0.0012** (expect ≈ 0), label FPR **0.0093**.

Three arms, identical labels, genes and lines:

| Arm | Score |
|---|---|
| `panel` | van der Waerden normal scores across the whole panel (pipeline baseline) |
| `lineage` | vdW normal scores within tissue, **no guards** (what §C.15 tested) |
| `guarded` | robust z within tissue **with** all three guards (`MIN_PEERS_FOR_Z`, silent-lineage, `MAD_FLOOR`) |

438 genes scored, median 731 lines/gene, median prevalence 0.124, 21 lineages with ≥15 lines.

### Result — VERIFIED

The guard makes some lines unscoreable, so pAUC over different line sets is not
like-for-like. Both are reported.

**(a) All lines — each arm on whatever it can score**

| Arm | Median gate pAUC | Genes |
|---|---|---|
| panel | 0.5222 | 438 |
| lineage | 0.5110 | 438 |
| guarded | 0.5104 | 412 |

**(b) Matched — every arm restricted to the lines the guard retains** (389 genes)

| Arm | Median gate pAUC |
|---|---|
| panel | 0.5278 |
| lineage | 0.5090 |
| guarded | 0.5062 |

**Paired bootstrap over genes, 10,000 resamples:**

| Comparison | Median Δ | 95% CI | p | Better in |
|---|---|---|---|---|
| **guarded − lineage** | **+0.0022** | **[−0.0072, +0.0136]** | 0.332 | 51.4% |
| guarded − panel | −0.0004 | [−0.0241, +0.0098] | 0.586 | 48.3% |
| lineage − panel | −0.0120 | [−0.0256, +0.0053] | 0.162 | 44.2% |

### Verdict

**Negative.** The guards do not recover an effect that plain lineage conditioning does not.
`guarded − lineage` is +0.0022 with a CI comfortably spanning zero and a better-in rate of
51.4% — a coin flip. The `lineage − panel` arm reproduces §C.15's neutral result
independently (−0.0120, CI spans zero, better in 44.2%), so the harness is behaving.

All three arms sit at pAUC ≈ 0.51, which is near chance.

### Why — and the important caveat

**The guards almost never fire on this gene sample.** The guard retains a median of
**100.0%** of lines per gene, and drops all lines for only 22 of 438 genes. A guard that
does not bind cannot change a ranking metric.

That is not a defect in the guards. It is a mismatch between what the guards do and what
this harness measures:

- The harness measures **ranking quality** among scored lines.
- The guards produce **abstention** — they refuse to score degenerate cases.

A guard's contribution is in a different currency, and gate-region pAUC cannot price it.
The genes on which the guards *do* bind are exactly the genes excluded from the harness,
because a gene with no scoreable lines produces no pAUC.

**So the honest conclusion is two-part, and the second part is not a consolation prize:**

1. **For dependency prediction, the guarded lineage scorer is not better than the plain
   one, and neither is better than the panel-wide baseline.** The layer must be rescoped —
   see F4 and `TRANSCRIPTOMICS_DECISIONS.md` Q3.
2. **The silence guard is nonetheless valuable, and F2 measures where.** It is just not
   valuable as a ranking improvement.

---

## F2. The salvageable part — the silence guard as a validity check on the main pipeline

The guard asks a question that does not require the group to be a lineage: **is this
distribution so degenerate that a relative position within it is meaningless?** The main
pipeline's pooled percentile has the identical failure mode and no guard against it.

### Specification

A validity check on the main pipeline's own scores, independent of lineage:

```
For each gene g, over the scoring stratum S (the pooled panel today):
    frac_expressed(g, S) = |{lines in S with log2(TPM+1) > EXPRESSED_MIN}| / |S|

    if frac_expressed(g, S) < SILENT_FRAC:
        score_validity(g, ·) = 'UNINFORMATIVE — gene silent in stratum'
        percentile is withheld, not reported with a caveat
```

This maps directly onto the abstention states the query layer already defines
(`docs/QUERY_LAYER_STATES.md`: `NO_EVIDENCE` / `UNINFORMATIVE`). It is an
`UNINFORMATIVE` producer, and it needs no new output vocabulary.

### How many gene × line combinations would be flagged — VERIFIED

Measured over the main pipeline's full DepMap layer: 19,173 protein-coding genes ×
1,479 lines = 28,356,867 gene×line combinations.

| Statistic | Value |
|---|---|
| median fraction of lines expressing a gene | 0.803 |
| IQR | 0.111 – 0.999 |
| **genes failing the guard at `SILENT_FRAC = 0.20`** | **5,847 of 19,173 (30.5%)** |
| **gene × line combinations affected** | **8,647,713 of 28,356,867 (30.5%)** |

Sensitivity to the cut:

| `SILENT_FRAC` | Genes flagged | % | Gene × line pairs |
|---|---|---|---|
| 0.05 | 3,565 | 18.6% | 5,272,635 |
| 0.10 | 4,607 | 24.0% | 6,813,753 |
| **0.20** | **5,847** | **30.5%** | **8,647,713** |
| 0.30 | 6,623 | 34.5% | 9,795,417 |

### A spread floor does not catch these — the notebook's central claim, confirmed at scale

Cross-checked against the pipeline's own `gene_dispersion.parquet`. Of the 5,847 genes
failing the silence guard, the pipeline's existing `signal_spread` label says:

| `signal_spread` | Genes |
|---|---|
| narrow | 2,549 |
| **wide** | **2,420** |
| moderate | 878 |

**2,420 of them are labelled "wide"** — a spread-based check does not merely miss them, it
actively certifies them as well-dispersed. And on dynamic range:

| | Value |
|---|---|
| median `dynamic_range` of the flagged genes | 1.947 |
| IQR | 0.496 – 3.863 |
| **exceeding `FLAT_FOLD_CEILING = 3.0`** | **2,091 genes** |

**2,091 genes have a dynamic range above the pipeline's own flatness ceiling while being
expressed in under 20% of the panel.** A floor alone does not catch this — exactly the
notebook's decisive observation (EGFR in blood: MAD 0.178, comfortably above any sensible
floor, 38% of values at exactly zero), now measured across the whole gene universe rather
than on one example.

**This is the transcriptomics layer's most valuable contribution to the project**, and it
is worth more than the scoring method it was built to support. It bears directly on §9
gap 13, the main pipeline's highest-priority open item.

### Required before adoption

`SILENT_FRAC` must be swept — {0.05, 0.10, 0.20, 0.30}, held-out selection — and the
abstention rate reported at each. Flagging 30.5% of the product's output is a large
change and the cut should not be inherited unexamined from a different codebase.

---

## F3. Source roles — settled, and for the right reason

### The specification

| Source | Role | Permitted to |
|---|---|---|
| **DepMap** + **HPA** | **one combined RNA opinion** | contribute to the score; judge whether a gene is off; compute τ |
| **GEO** | **held out** | corroborate (categorical), display |
| | | **barred** from the score, **barred** from judging whether a gene is off, **barred** from τ |

**DepMap and HPA must never be counted as two confirmations.** Measured ρ = 0.863 on
z-scores (0.829 on raw values). Counting them as two overstates the combined Z by
**36.5%**; they are worth **1.07** independent observations. They form one opinion, combined
under the Strube/Hartung correction with the measured R.

**GEO is held out because it is floored — not because it is redundant.**

This distinction is the whole point and it is easy to get backwards. On the measurements:

| Source pair | z-score ρ |
|---|---|
| DepMap–HPA | **0.863** |
| DepMap–GEO | 0.495 |
| HPA–GEO | 0.530 |

**GEO is the *least* redundant of the three.** It carries more information not already
present in the other two than either of them carries relative to the other. The tempting
design — keep the two sources that agree, discard the one that disagrees — is **exactly
backwards on information content**. It would discard the most independent source and keep
the near-duplicate pair.

GEO is held out for a different and specific reason: **it has a detection floor and
therefore cannot distinguish "not expressed" from "below detection".** Measured directly
by the notebook's own `source_floor_check()`:

| Source | n | min | p01 | frac below `EXPRESSED_MIN` | floored |
|---|---|---|---|---|---|
| depmap_expr | 1,495 | 0.000 | 0.000 | 0.1545 | False |
| **geo_expr** | **591** | **3.030** | **3.554** | **0.0000** | **True** |
| hpa_rna | 1,106 | 0.000 | 0.000 | 0.1799 | False |

GEO has **zero** values below the detection threshold and a hard minimum at 3.03. It
reports background as signal.

So the correct statement of the roles is:

> GEO is excluded from the score and from silence judgements because its detection floor
> makes it unable to represent absence, not because it is redundant. It is in fact the
> least redundant of the three sources, and its exclusion is a real cost paid for a
> specific reason. Where GEO's coverage is its own — 37 cell lines not covered by DepMap or
> HPA — that cost is a coverage loss, and it should be reported rather than absorbed.

A second cost, worth stating: GEO adds only 37 lines (2.3%) to the union of 1,580, so the
coverage loss is small. The information loss is not small, and it is the one to name.

**One correction to the notebook's own framing.** Cell 4 argues the floor "matters less
than it first appears" because z is computed within source and so the floor cancels, and
concludes "the floor breaks exactly one thing: judging whether a gene is off". That
argument is sound for the *z-score* but incomplete for the *combination*: a floored source
also compresses variance in silent tissues, which changes its z distribution and hence its
contribution to the Stouffer sum, not just its silence vote. Given GEO is being held out
of the score anyway, this is moot in the recommended design — but the notebook's stated
justification is stronger than the evidence supports and should be softened if GEO ever
re-enters the score.

---

## F4. What this layer can and cannot deliver — for the record

Stated plainly, as the brief requires.

1. **It produces a relative position within a group, not a quantity.** Every output —
   `z_combined`, `z_shrunk`, `pct_expressed`, `p_emp` — is a statement about where a line
   sits relative to its lineage peers. Change the peer set and the number changes. Nothing
   in it is a property of the cell line alone.
2. **It cannot answer "how much of this protein or transcript is present."** No absolute
   quantity is computed or recoverable. The `log2(x+1)` transform is applied to three
   sources on three different scales (`log2(TPM+1)`, `nTPM`, MAS5-like linear intensity),
   and the z-scoring then removes what remained of the level information. This is a
   deliberate trade for cross-assay comparability and it is the right trade — but it is a
   trade, and magnitude is what was given up.
3. **Its ordering has never been validated against any external outcome.** Both validation
   tests are internal — a permutation of its own data, and planted signal of its own
   construction. Neither uses Chronos, GDSC, Project Score, or any measured biology. Until
   F1 was run in this audit, the layer had been checked only against its own definition of
   unusual. F1 is now the first external check, and it is negative (pAUC ≈ 0.51 across all
   three arms).

### Does any downstream use assume otherwise? — VERIFIED

**No, because there is no downstream use.** `02_transcriptonomics.ipynb` is not in the
repository, is imported by nothing, and writes no artefact that any pipeline stage reads.
The layer is currently standalone.

Two forward-looking warnings, since the point of this task is integration:

- **The main pipeline's `core_score` is also a relative position** (a within-gene
  percentile across the panel) and is subject to the same three limits. If the two are ever
  combined, they are not independent evidence about different things — they are two
  relative positions computed from overlapping data. The A2 correlation applies.
- **Point 3 is the one that will be challenged.** "Validated" appears in the notebook's
  verdict vocabulary (`z_shrunk` thresholds producing `HIGH`/`LOW`) and in the reporting
  template in cell 37, which quotes a false-positive rate of 1.2% and 100% recovery at
  3 MAD. Both figures are superseded by this audit — the true operating rate is α = 0.116
  (B1) and recovery at 3 MAD is 85% (B2). **The reporting template in cell 37 should not be
  used as written.**
