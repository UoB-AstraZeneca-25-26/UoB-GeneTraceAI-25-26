> **ARCHIVED 2026-08-14 — NUMBERS NOT REPRODUCIBLE**
>
> Original figures (hr_flat=0.014025, hr_driver=0.017735, 28 genes) were computed by
> `src/pipeline/stage4_eval_save.py`, which no longer exists in the repository.
> The flags_with_driver file at the time had 885,844 rows; the current file has
> 19,704,388 rows (22x). The generating code and data state are both absent.
> These numbers cannot be regenerated from any source in the current repo.
>
> Superseded by: `docs/ALTERATION_DEFENCE.md` (pipeline state 2026-08-14,
> 19.7M flag rows, 26 test genes, three-method bootstrap).
>
> Retained for audit trail. Do not cite these numbers in any output.

---

# Alteration — what is measured, from what, and what it actually does

Defence notes for the alteration layer. Every number is measured on this
repository (2026-08-06) and cited to the file that produces it.

Companions: [TRACK_C_EDA.md](TRACK_C_EDA.md) (evidence quality),
[TRACK_C_PROVENANCE.md](TRACK_C_PROVENANCE.md) (how the data was collected),
[MATH_REFERENCE.md](MATH_REFERENCE.md) (formulae).

---

## 1. What exactly is being measured

**The unit is the (gene, cell line) pair.** Not a gene, not a cell line — the
pair. The claim the system makes is:

> in cell line *L*, gene *G* carries a somatic alteration of a class we believe
> changes its function

"Gene X is altered" and "cell line Y is altered" are both category errors in this
architecture. A gene is altered *in a line*; a line is altered *at a gene*.

Three separate quantities are produced per pair, and they have different fates:

| quantity | type | where it lives | affects ranking? |
|---|---|---|---|
| `has_driver_alteration` | boolean | `flags_with_driver.parquet` | **yes — this is the production gate** |
| `has_alteration` | boolean | `flags_with_driver.parquet` | no — rejected comparison arm |
| `p_mutation` / `p_fusion` | float [0,1] | `flags_with_driver.parquet` | **no — display only** |

---

## 2. Which datasets, which columns, which values

### Resources actually consumed

| resource | file | rows | columns used |
|---|---|---:|---|
| DepMap somatic mutations | `cleaned_track_data/mutations_collapsed.parquet` | 632,919 | `any_driver`, `oncogene_hit`, `tsg_hit` (**gate**); `max_vep_rank`, `max_pathogenicity`, `variant_count` (score) |
| DepMap mutation detail | `cleaned_track_data/mutations_variant_detail.parquet` | 704,713 | `af`, `vepimpact`, `revelscore`, `ampathogenicity`, `likelylof`, `hotspot` |
| DepMap fusions | `cleaned_track_data/fusions_gene_level.parquet` | 144,221 | `max_confidence` (**gate**); `best_ffpm`, `fusion_count`, `any_in_frame` (score) |
| COSMIC CNA | `outputs/cna_flags.parquet` | 117,540 | `max_cn`, `min_cn`, `cna_type` → `has_cna_alteration` |
| DepMap signatures | `cleaned_track_data/signatures_model_level.parquet` | 1,955 | `msiscore`, `ploidy` → `signature_discount` (**computed, never wired in**) |
| COSMIC CGC | warehouse `main.cosmic_cgc` | 768 | `gene_role`, `tier` — the source of "driver" |

Scale: 18,058 genes × 1,744 lines for mutations, 16,017 × 1,699 for fusions.

### Values in the shipped artefact

`flags_with_driver.parquet` — 885,844 rows × 13 columns:

| column | True | False |
|---|---:|---:|
| `has_alteration` | 539,099 | 346,745 |
| `has_driver_alteration` | **160,580** | 725,264 |
| `has_cna_alteration` | 3,217 | 882,627 |

`p_mutation` mean 0.423 (median 0.452); `p_fusion` mean 0.140 (median **0.000** —
the majority of pairs have no fusion evidence at all).

---

## 3. "How can you tell a mutation was present or not?"

**Present** means, precisely: a row exists in `mutations_collapsed` for that
`(ensg_id, model_id)`, which means at least one variant in that gene in that line
survived DepMap's Mutect2 + `FilterMutectCalls` pipeline **at VAF ≥ 0.15**, in a
sequencing profile that resolved to that model.

That is a chain of four conditions, and each can fail independently:

1. the line was sequenced (WES/WGS) and the profile resolved to a model
2. the variant was called by Mutect2 without a matched normal
3. it passed the filters, including the **VAF ≥ 0.15 floor** (measured: minimum
   `af` is exactly 0.150)
4. it survived Track C's universe filter

**Absent** currently means "no row" — and that is the weakest part of the design,
because no row conflates *at least* three different situations:

- the line was sequenced and nothing was found → genuine absence
- the line was never sequenced → unknown
- a variant exists below 15% VAF → undetectable

This is exactly what decision Q1 in [TRACK_C_EDA.md](TRACK_C_EDA.md) fixes, by
splitting `detected` into `True` / `False` / `NA`.

---

## 4. "How are you managing both fusion and mutation together?"

They are scored **separately on the same key**, then combined by boolean OR at
two different levels:

```
mutations_collapsed ──► mut_driver    = any_driver ∨ oncogene_hit ∨ tsg_hit
                                                                    │
                                                                    ├─► has_driver_alteration   (PRODUCTION GATE)
                                                                    │
fusions_gene_level  ──► fusion_driver = (max_confidence == "high") ─┘


mutations_scores ──► p_mutation ┐
                                ├─► has_alteration = (p_mut > 0.5) ∨ (p_fus > 0.5)   (REJECTED ARM)
fusions_scores   ──► p_fusion  ─┘
```

Both routes are **OR, not sum**. A pair with a driver mutation and a
high-confidence fusion is not ranked above a pair with only one — the gate is
binary and saturates at the first piece of evidence.

The continuous scores use Noisy-OR over Hill-transformed components:

```
p_mutation = 1 − (1−p_vep)(1−p_path)(1−p_burden)  +  0.15·driver_flag     (clipped to 1)
p_fusion   = 1 − (1−p_conf)(1−p_ffpm)(1−p_recur)  +  INFRAME_BOOST·inframe_flag

Hill(x; p₀,k) = x^k / (x^k + p₀^k)
   VEP (3.0, 2.0) · pathogenicity (0.5, 2.0) · burden (3.0, 1.5)
```

Noisy-OR assumes the components are **conditionally independent**. They are not —
`max_vep_rank`, `max_pathogenicity` and `variant_count` are all functions of the
same variant set, and the assumption has never been tested (Gap 2 in
`STATE_REPORT.md`).

**The clean answer if challenged:** the independence assumption is violated in
principle but inconsequential in practice, because `p_mutation` and `p_fusion`
are display annotations that do not enter the sort key. The ranking is unaffected
by this assumption. It would need to be resolved before this layer could be
promoted into ranking — and that is the condition under which it becomes a real
defect rather than a documented one.

---

## 5. "On what factor do you say the pair is altered?"

**One boolean, from four columns.**

```
has_driver_alteration = any_driver ∨ oncogene_hit ∨ tsg_hit ∨ (max_confidence == "high")
```

- `any_driver` — OR of five DepMap driver-annotation flags
- `oncogene_hit` / `tsg_hit` — direction-preserving, from DepMap's
  `oncogenehighimpact` / `tumorsuppressorhighimpact`, which are themselves keyed
  to COSMIC CGC's 768 genes
- `max_confidence == "high"` — STAR-Fusion's own confidence tier

**Neither VAF, nor pathogenicity score, nor FFPM enters this decision.** The
"how much" work in Track C's scoring notebooks produces `p_mutation`/`p_fusion`,
and those never reach the sort order. This is the single most important thing to
be able to say out loud in a viva: *the production system uses alteration as a
binary driver flag; the continuous scores are shown to the user but do not rank.*

### Where it enters the ranking

`stage4_eval_save.py:88-96` and `build_full_predictions.py`:

```
sort_driver = where(class == "abundance_tracking",
                    core_score,                                  # gate bypassed
                    has_driver_alteration · 1e6 + core_score)    # gate applied
```

The `1e6` makes it a strict lexicographic sort — every driver-flagged pair
outranks every unflagged one, and `core_score` only breaks ties within each
block. So alteration is not a weighted feature; it is a **hard partition** of the
ranking, and only for genes classed `activation_driven`.

---

## 6. What alteration actually contributes — measured

From `stage4_eval_summary.json` (28 held-out genes, seed 42, 113 train / 28 test):

| arm | hit rate | vs flat |
|---|---:|---|
| `hr_flat` (core_score only) | 0.014025 | — |
| `hr_driver` (**production**) | 0.017735 | +26.5% |
| `hr_any_alt` (rejected arm) | **0.020103** | **+43.3%** |

Three things a viva will push on, so they should be said first:

1. **The rejected arm scores higher than the production gate.** 0.0201 vs 0.0177.
   The driver gate was chosen for **variance reduction** (δ std 0.0126 → 0.0092,
   a 26.7% reduction), not for higher recall. That is a defensible choice — but
   it is a choice, and it should be stated as one. §6.1 shows the recall gap is
   not statistically distinguishable, which is what makes the choice defensible.
2. **The absolute numbers are tiny.** A 26.45% "lift" is a movement from 0.0140
   to 0.0177 — 0.0037 in absolute hit rate. Interval now attached, §6.1.
3. **The gate does nothing for `abundance_tracking` genes.** Measured in the
   per-class breakdown: `hr_flat = hr_any_alt = hr_driver = 0.046717`, delta
   exactly 0.0. The gate is bypassed by construction for that class, which is
   2 of the 28 test genes. Diagnosed in §6.2.

### 6.1 Confidence intervals — paired bootstrap over genes

`test_run_driver_gate_ci.py`, 2,000 resamples, seed 42. The bootstrap unit is the
**gene**, matching how the hit rate is aggregated; resampling pairs within a gene
would understate the gene-to-gene variance that actually matters.

| comparison | point | 95% CI | verdict |
|---|---:|---|---|
| driver − flat | +0.003710 | **[+0.000445, +0.007018]** | excludes zero |
| any-alt − flat | +0.006079 | [+0.001736, +0.010492] | excludes zero |
| **driver − any-alt** | −0.002369 | **[−0.005205, +0.000293]** | **straddles zero** |
| lift % | +26.45% | **[+2.68%, +61.21%]** | excludes zero, very wide |

Three readings, and all three should be volunteered rather than extracted:

- **The driver gate's improvement over flat ranking survives**, but only just —
  the lower bound is +0.000445, a hair above zero. The 26.45% headline has a
  95% interval of **[2.7%, 61.2%]**. Quote the interval, never the point estimate
  alone.
- **The driver-vs-any-alteration gap straddles zero.** The rejected arm looks
  better on the point estimate, but the two arms are **not statistically
  distinguishable** at n=28. This is what makes the variance-reduction argument
  the deciding one: when two arms tie on recall within noise, choosing the
  lower-variance one is the correct call, not a rationalisation.
- **Sign test:** 14 of 28 genes improved under the gate, 8 unchanged, 6 worsened.
  A majority, not a sweep.

### 6.2 The `abundance_tracking` bypass — diagnosed

The two bypassed genes, with the counterfactual: what would the hit rate have
been if the gate *had* been applied to them?

| gene | CGC role | lines driver-flagged | hr shipped (bypassed) | hr if gated | delta |
|---|---|---:|---:|---:|---:|
| **KDM6B** | unknown | 19 / 1,482 (1.28%) | 0.056000 | 0.024000 | **−0.032000** |
| **PIM1** | **oncogene, tier 1** | 5 / 1,582 (0.32%) | 0.037433 | 0.042781 | **+0.005348** |

The suspicion that the bypass might be hiding a gate that would have helped is
**half right, and worth saying so**:

- **KDM6B** is genuinely not driver-informative — no CGC role, 1.28% of lines
  flagged, and applying the gate would have *cut its hit rate by more than half*.
  The bypass is correct here.
- **PIM1 is a tier-1 CGC oncogene classed as `abundance_tracking`.** Applying the
  gate would have improved it, by +0.005348. So the class assignment is arguably
  wrong for this gene — but only 5 of its 1,582 lines carry a driver flag at all,
  so the gate has almost nothing to act on.

**Net across both: −0.0267.** Applying the gate to `abundance_tracking` genes
would have made the class worse overall, so the bypass is defensible as shipped.
The honest caveat is that this rests on n=2 genes, and PIM1 shows the regime
classifier and the CGC roster can disagree about the same gene.

### The fabricated-negative problem

`stage4_eval_save.py:34-35` does `.fillna(False)` on both alteration booleans
after a left join. Measured: `predictions_with_confidence` has **29,781,274**
rows, `flags_with_driver` has **885,844**. So **28,895,430 pairs (97.0%)** are
assigned `has_alteration = False` **by imputation, not by measurement**.

Those fabricated negatives are inside the evaluation harness that produces the
headline KPI. This is the same defect as the gate-denominator artefact already
recorded for the depletion headline, and it is the strongest argument for the
three-state `detected` field.

---

## 7. KPIs, and what each is actually sensitive to

| KPI | value | source | sensitive to |
|---|---:|---|---|
| hit rate @ top-5% (`hr_driver`) | 0.017735 | `stage4_eval_summary.json` | driver flags + core_score |
| driver lift | 26.45% **[2.68, 61.21]** | `test_run_driver_gate_ci_results.json` | ratio of two small numbers over n=28 |
| driver − flat (absolute) | +0.003710 **[+0.000445, +0.007018]** | same | excludes zero, barely |
| driver − any-alt | −0.002369 **[−0.005205, +0.000293]** | same | **straddles zero** — arms not distinguishable |
| variance reduction | 26.7% | `stage4_eval_summary.json` | the actual reason the gate was kept |
| driver / any ratio | 160,580 / 539,099 = **0.298** | `04_driver_gated_routing.ipynb` | CGC's 768-gene roster |
| pairs with any alteration | 539,099 | `flags_with_driver.parquet` | the p > 0.5 threshold |
| CNA contribution | 3,217 of 885,844 (**0.36%**) | same | COSMIC CNA coverage |

The last row needs its explanation attached, or it reads as a pipeline failure.
It is not one: **CNA coverage is limited by COSMIC's curation scope, not by a
defect in the join.** `cosmic_cna` is a curated event table listing called
copy-number regions — roughly 143 genes per model — rather than a genome-wide
copy-number matrix. Genome-wide CNA integration would require the output of a
separate copy-number caller (PureCN/PURPLE segmentation, or DepMap's own
`OmicsCNGene`), which is out of scope for this build. The 0.36% is the honest
reach of the curated source, and the join itself resolves correctly.

**No multiple-testing correction exists anywhere in the pipeline** (Gap 3).

---

## 8. "If I say Gene A excluding Gene B — is alteration a factor?"

Two readings, and the answer differs.

### (a) As a query — "A altered AND B not altered"

Expressible in principle: `has_alteration` is per (gene, line), so a boolean
exclusion is just a filter. **But it is not currently sound**, for one specific
reason.

The `NOT` clause is where the missing third state does maximum damage. "Gene B is
not altered" evaluates `has_alteration == False` — and 97.0% of pairs hold that
value by `.fillna(False)`, not by measurement. So a query for *"lines where A is
altered but B is not"* silently returns every line where **B was simply never
assayed**, presented identically to lines where B was assayed and found clean.

An exclusion query is strictly more dangerous than an inclusion query here,
because inclusion only ever returns measured positives, while exclusion returns
the entire unmeasured space. **This should be fixed before any exclusion syntax
is exposed to users.** With three-state `detected`, the correct semantics are:

```
A altered ∧ B assayed ∧ B not altered      -- sound
A altered ∧ ¬(B altered)                   -- currently returns the unmeasured space
```

### (b) As biology — mutual exclusivity

**No. Not modelled at all, and it is a real gap.**

Every (gene, line) pair is scored independently. There is no cross-gene term
anywhere in the architecture — no co-occurrence, no exclusivity, no pathway
grouping. The system cannot currently represent "KRAS and EGFR are rarely mutated
in the same lung line because they occupy the same pathway node", which is one of
the most robust findings in cancer genomics.

This is a genuine limitation, not an oversight to hide. If asked what it would
take:

- test each gene pair for co-occurrence/exclusivity — Fisher exact or a
  permutation test on the (gene, line) alteration matrix
- correct for the confounds already measured in
  [TRACK_C_EDA.md](TRACK_C_EDA.md) §2.4, since gene length (ρ = 0.464),
  hypermutators (top 1% of lines = 11.7% of pairs) and lineage (22.4× spread)
  all manufacture spurious co-occurrence
- correct for multiple testing across ~163M gene pairs — which the pipeline
  currently has no machinery for at all
- established methods to cite: DISCOVER, CoMEt, MEMo

Honest position: **mutual exclusivity is a known, unimplemented extension**, and
the confound correction is a prerequisite, not an afterthought.

---

## 9. One-line answers

- **What are you measuring?** Whether a (gene, cell line) pair carries a somatic
  alteration believed to change gene function.
- **What datasets?** DepMap somatic mutations (Mutect2, WES/WGS), DepMap fusions
  (STAR-Fusion on RNA), COSMIC CNA, COSMIC CGC for driver roles.
- **How do you know a mutation is present?** A row exists — meaning ≥1 variant
  passed Mutect2 filtering at VAF ≥ 0.15 in a resolved sequencing profile.
- **How do you know it is absent?** You don't, cleanly. That is the known defect.
- **Which columns decide it?** `any_driver`, `oncogene_hit`, `tsg_hit`,
  `max_confidence == "high"`. Four columns, OR'd.
- **Which values?** Booleans. Not VAF, not pathogenicity, not FFPM.
- **How are mutation and fusion combined?** Boolean OR, at both the gate and the
  score level. Never summed.
- **What does it contribute?** A hard lexicographic partition of the ranking for
  `activation_driven` genes: +0.0037 absolute hit rate, 95% CI
  [+0.000445, +0.007018], and a 26.7% variance reduction — which is the real
  reason it was kept.
- **Is exclusion supported?** Syntactically yes, soundly no, until `detected`
  has three states.
- **Is mutual exclusivity modelled?** No. Known gap, prerequisites identified.

---

## 10. Say these unprompted

Three statements that, volunteered before being asked, demonstrate command of the
system's own weaknesses:

1. **"97% of absent alteration flags are imputed by `fillna`, not measured."**
   `.fillna(False)` inside the evaluation harness assigns `has_alteration = False`
   to 28,895,430 of 29,781,274 pairs. The fix is a three-state `detected` field —
   assayed-and-absent, assayed-and-present, never-assayed. Until that exists,
   exclusion queries return the entire unmeasured space, not confirmed negatives.
2. **"The driver gate was kept for variance reduction, not recall."** The simpler
   any-alteration arm scores higher on the point estimate (0.0201 vs 0.0177) —
   but the difference straddles zero (95% CI [−0.0052, +0.0003]), so the arms are
   not statistically distinguishable at n=28. When recall ties within noise,
   picking the lower-variance arm is the correct call. That was a documented
   choice, not a post-hoc rationalisation.
3. **"The lift is real, but narrow."** Driver versus flat is +0.003710, 95% CI
   **[+0.000445, +0.007018]** — the lower bound is a hair above zero. Lead with
   that modesty, not with the 26.45% headline, whose own interval is
   [2.7%, 61.2%]. The effect excludes zero, so it is real; n=28 genes cannot
   support a precise claim about its size.
