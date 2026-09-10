# RANKING_FEASIBILITY — can a confirmation count order anything?

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

Part A of the ranking brief. **Verdict: negative, twice.**

> **Re-run, 2026-08-09.** §7 identified the unused copy-number file as the only available
> lever. It has since been wired in (`src/pipeline/build_cn_layer_complete.py`, a 141×
> increase in CN pairs) and Part A re-run against a decision rule pre-specified beforehand in
> `docs/RANKING_PRESPEC.md`. **Verdict: ABANDON**, unanimous across all three pre-specified
> thresholds — and every criterion moved in the wrong direction. See **§10** for the re-run
> and `docs/COVERAGE_RESULT.md` for why more data made it worse.
>
> Separately, the evaluation denominator has been fixed
> (`src/pipeline/build_eval_tested_denominator.py`). It does not affect the verdict below,
> but it does change a different earlier claim — see `docs/EVAL_DENOMINATOR.md` §4.

Standing rule 7 applies and the task stops here. Parts B–F were not attempted; §7 states what
would have to change first, and §10 reports what happened when it did.

Measured 2026-08-09 against the repository as committed on branch `chai`. Analysis scripts:
`docs/audit_scripts/rank_checks.py`, `partA.py`, `partA_repair.py`, `partA_within_gate.py`,
`partA_v2.py`.

**Code added (§10 re-run only, both additive — no existing pipeline file was modified):**
`src/pipeline/build_cn_layer_complete.py` and `src/pipeline/build_eval_tested_denominator.py`.
Sections 1–9 below were produced with no production code written at all.

Tags: `VERIFIED` (computed here) · `UNVERIFIED — could not locate`.

---

## 1. Verdict

**A confirmation count over the five candidate checks cannot order cell lines.** Two
independent failures, either of which is disqualifying:

| Failure | Raw count | Fraction-of-assessed variant |
|---|---|---|
| **Coverage confound** (A1) — Spearman with gene-independent data availability | **+0.583** | **−0.391** |
| **Ties** (A2) — largest tied group as a fraction of retained lines | **70.5%** | **83.8%** |
| **Same answer every gene** (A3) — Spearman between two random genes' orderings | **+0.531** | +0.181 |

*(All figures on valid genes, measured within the gate-retained set — the level at which the
count is actually proposed to operate.)*

The raw count ranks data availability. The fraction variant removes that confound but
inverts it — it penalises lines for having more checks assessed — and ties 84% of lines.

**The underlying cause is not an arithmetic problem and cannot be fixed by choosing a
different normalisation: there is almost nothing to count.** On gate-retained pairs, three
of the four remaining checks have no data for 68–99.7% of pairs, and the fourth passes
97.15% of the time. **58.0% of retained pairs have exactly one other check assessed**, so
for the majority of the panel a "count of independent confirmations" is a count over a
single observation.

Building the ordering on this would be worse than having no ordering, exactly as the brief
anticipated.

---

## 2. What was built and measured

### 2.1 The five checks, three-state

Each check is `passed` / `failed` / `not assessed` per (gene, line), per standing rule 2.
Definitions used (`rank_checks.py`):

| Check | Assessed when | Passes when | Data source |
|---|---|---|---|
| `rna_gate` | DepMap or HPA covers the line | consensus `log2(TPM+1) > 1` (i.e. TPM > 1, the HPA detection convention) | warehouse `depmap_expr`, `hpa_rna` |
| `cna_gate` | COSMIC CNA context exists for the pair | `cna_type` is not `loss`/`both` | `cna_flags.parquet` |
| `protein` | protein quantified for that line | above the protein's panel median | `procan_proteomics.parquet` via `procan_protein_map` + `gdsc_models` bridge |
| `alteration` | line has mutation or fusion coverage | no alteration present (**demote direction assumed**, see §6) | `flags_with_driver.parquet`, `coverage_matrix_enriched.parquet` |
| `geo_corrob` | GEO covers the pair | GEO within-gene percentile within 0.25 of the RNA-consensus percentile | warehouse `geo_expr` |

300 randomly sampled protein-coding genes (seed 42), 1,543 cell lines, **462,900 (gene, line)
pairs**. Per standing rule 3 every measurement is reported twice: **all genes** and **valid
genes only** (valid = passes the silence guard, expressed in ≥20% of the panel — 231 of 300,
77.0%).

### 2.2 State distribution — the result that explains everything else

Over all 462,900 pairs (`VERIFIED`):

| Check | passed | failed | **not assessed** |
|---|---|---|---|
| `rna_gate` | 0.6431 | 0.3569 | **0.0000** |
| `cna_gate` | 0.0021 | 0.0016 | **0.9963** |
| `protein` | 0.0650 | 0.0653 | **0.8697** |
| `alteration` | **0.9721** | 0.0156 | 0.0123 |
| `geo_corrob` | 0.1981 | 0.1114 | **0.6906** |

Only `rna_gate` is both universally assessed and materially discriminating. `cna_gate` is
**not assessed for 99.63% of pairs** — it is inert. `alteration` passes 97.21% of the time —
it is inert in the other direction. `protein` and `geo_corrob` are absent for the majority
of pairs.

This is the C3 three-state finding from the earlier audit reappearing as a design
constraint: once "not measured" is honestly separated from "measured negative", there is
very little measured anything.

---

## 3. A1 — the coverage confound (decisive)

Gene-independent evidence-type availability per line (0–5), from
`coverage_matrix_enriched.parquet`:

| Evidence types available | Lines |
|---|---|
| 1 | 9 |
| 2 | 599 |
| 3 | 153 |
| 4 | 389 |
| 5 | 393 |

### Correlation between coverage and confirmation count — `VERIFIED`

| Measurement | All genes | **Valid genes only** |
|---|---|---|
| pooled Spearman(coverage_count, n_pass), all lines | +0.331 | **+0.445** |
| per-gene Spearman, median | +0.470 | **+0.540** (IQR +0.406 to +0.642) |
| within-lineage pooled Spearman, median over 17 lineages | +0.307 | **+0.421** |
| **within the gate-retained set** (where Level 2 operates) | +0.581 | **+0.583** |

The confound **strengthens** on valid genes and strengthens again inside the retained set —
the opposite of what a real biological signal would do. Per standing rule 3, the
valid-subset number is the one that counts: **+0.583**.

### The sharpest form of the same measurement

`Spearman(n_assessed, n_pass) = +0.545`, Pearson **+0.603**. Mean confirmation count by
number of checks with data:

| checks assessed | pairs | mean n_pass |
|---|---|---|
| 1 | 3,363 | 0.63 |
| 2 | 286,859 | 1.59 |
| 3 | 142,588 | 2.24 |
| 4 | 29,929 | 3.12 |
| 5 | 161 | 4.34 |

**Each additional check that happens to have data adds ≈0.85 to the count**, near-linearly,
irrespective of biology. That is the definition of a coverage proxy.

### Coverage is itself confounded, so the proxy is actively misleading — `VERIFIED`

- **Mean RNA:** Spearman(coverage_count, per-line mean RNA) = **+0.139**, p = 4.37e-08,
  n = 1,543. Consistent with the previously measured proteomics coverage bias
  (Cliff's δ = +0.251).
- **Lineage:** Kruskal–Wallis H = 165.3, **p = 7.05e-25** across 21 lineages. Per-lineage
  mean coverage ranges from **2.18 (biliary tract) to 4.15 (breast)** — a near two-fold
  spread.
- **Curated set:** 1 of the 300 sampled genes is in the 141-gene curated set, so curated-set
  membership could not be assessed on this sample. `UNVERIFIED — could not locate` a
  meaningful overlap; it would need a stratified sample to test.

So the ranking would systematically favour well-studied lineages and higher-expressing
lines, and would present that as evidence about the queried gene.

---

## 4. A2 — the tie structure

The count is not an ordering; it is a small number of buckets. `VERIFIED`.

Pooled over all (gene, line) pairs:

| Confirmation count | pairs | % |
|---|---|---|
| 0 | 3,532 | 0.76% |
| **1** | **140,653** | **30.39%** |
| **2** | **236,328** | **51.05%** |
| 3 | 72,534 | 15.67% |
| 4 | 9,771 | 2.11% |
| 5 | 82 | 0.02% |

**81.4% of all pairs sit at count 1 or 2.**

Per gene:

| | All genes | Valid genes only |
|---|---|---|
| largest tied group, as a fraction of lines (median) | 0.668 | **0.622** (IQR 0.569–0.728) |
| distinct count values achieved (median, of 6 possible) | 4 | 4 |
| lines per gene (median) | 1,543 | 1,543 |

**Within the gate-retained set**, which is where Level 2 would operate:

| Level-2 count (4 checks, gate excluded) | pairs | % |
|---|---|---|
| 0 | 5,873 | 1.97% |
| **1** | **209,899** | **70.51%** |
| 2 | 72,074 | 24.21% |
| 3 | 9,769 | 3.28% |
| 4 | 82 | 0.03% |

**94.7% of retained pairs fall into two values.** Largest tied group: median **70.5%** of
retained lines (raw), **83.8%** (fraction variant).

**What this means for the user.** For the typical gene, roughly 900 of 1,280 retained lines
carry an identical count. This is not a list; it is two tiers with a long tail. Presenting it
as a ranked list of 1,280 would imply an ordering that does not exist.

---

## 5. A3 — is the count a property of the line rather than the gene?

`VERIFIED`.

| | All genes | Valid genes only |
|---|---|---|
| fraction of variance **between lines** | 0.191 | **0.281** |
| per-gene ordering vs the gene-independent line-mean ordering (median Spearman) | +0.561 | **+0.628** |
| **Spearman between two random genes' orderings** (all lines) | +0.341 | **+0.371** |
| **the same, within the gate-retained set** | +0.500 | **+0.531** |

Within the retained set, **two unrelated genes produce orderings correlated at +0.531**. The
user would receive substantially the same list whichever gene they queried. That is a product
failure independent of whether the count is statistically defensible — and it follows
directly from A1, because the shared component is data availability.

---

## 6. Does the obvious repair rescue it? — no

The brief's B3 offers three ways to handle unmeasured checks. All three were tested before
concluding. `VERIFIED`, valid genes, measured over all lines:

| Variant | Coverage Spearman | Pairs retained |
|---|---|---|
| (i) raw count | **+0.445** | 100% |
| (ii) fraction of assessable, `n_pass / n_assessed` | **−0.208** | 100% |
| (iii) minimum assessed ≥ 2 | +0.435 | 99.3% |
| (iii) minimum assessed ≥ 3 | +0.144 | **40.0%** |
| (iii) minimum assessed ≥ 4 | −0.028 | **8.2%** |

- **(ii) over-corrects.** It moves the confound from +0.445 to −0.208 — and inside the
  retained set to **−0.391**. A score that penalises lines for having *more* data assessed is
  as misleading as one that rewards it, just harder to notice. It also worsens the ties:
  largest tied group **83.8%**.
- **(iii) works only by discarding the panel.** The confound is acceptable at ≥4 assessed
  checks (−0.028), but that retains **8.2% of pairs**. A rule that can order 8% of the panel
  is not an ordering rule for the product.
- To its credit, (ii) does fix the A3 product failure over all lines — between-gene ordering
  agreement drops from +0.371 to **+0.065**. But inside the retained set it is back to +0.181,
  and the tie structure makes the point moot.

### The count is largely the RNA gate wearing four other labels

| | All genes | Valid genes only |
|---|---|---|
| Spearman(n_pass, `rna_gate` alone) | +0.781 | **+0.629** |
| Spearman(fraction, `rna_gate` alone) | +0.861 | **+0.746** |

Mean n_pass is 1.20 when the RNA gate fails and 2.29 when it passes (valid genes) — the gate
alone accounts for most of the movement. The other four checks are not adding independent
confirmations; they are mostly not there.

### What is left to count, on retained pairs

| Check | passed | failed | **not assessed** |
|---|---|---|---|
| `cna_gate` | 0.0024 | 0.0004 | **0.9972** |
| `protein` | 0.0952 | 0.0951 | **0.8097** |
| `alteration` | **0.9715** | 0.0162 | 0.0123 |
| `geo_corrob` | 0.2198 | 0.1010 | **0.6792** |

Number of non-gate checks assessed per retained pair: **0 → 0.71%, 1 → 57.99%,
2 → 31.80%, 3 → 9.45%, 4 → 0.05%.**

For 58% of retained pairs the count has one input. Standing rule 1 (correlated checks count
once) never even gets to bite, because there is rarely more than one check to correlate.

**Note the assumed direction.** `alteration` was scored with alterations *demoting*
(unaltered passes), per the still-open decision in `DECISIONS_REQUIRED.md` Q1. Reversing it
would move 1.56% of pairs and would not change any conclusion above.

---

## 7. What this means, and what survives

### The exclusion stage is unaffected and remains the defensible part

None of this touches Level 1. The gate result stands: gate-region pAUC **0.596–0.608**,
p = 2.8e-76, replicated on two screens. What Part A shows is that **there is no second level
to put underneath it** — not that the first level is wrong.

### The fallback in F3 also needs redefining — `VERIFIED`

The brief's fallback is "the retained set, grouped into tiers by confirmation count". That
fallback depends on the count and therefore falls with it. Measuring the count-free
alternative — gate tier alone, using the three exclusion checks:

| | All genes | Valid genes only |
|---|---|---|
| survived all **assessed** exclusions (median lines/gene) | 1,165 | **1,202** (IQR 972–1,420) |
| excluded by exactly 1 (median) | 378 | 340 |
| excluded by 2 or more (median) | **0** | **0** |
| tiers actually populated, of 3 (median) | 2 | 2 |

And the reason the third tier is empty: **86.67% of pairs have only ONE of the three
exclusion checks assessed**; 13.26% have two; 0.07% have all three. A pair cannot be
"excluded by 2+" if only one exclusion was ever run on it.

**So the fallback as currently constituted returns ~1,200 lines in 2 tiers.** That is not the
"shortlist of twenty with honest ties" the brief's appendix describes as the goal. It is the
panel with about 22% removed.

This is the actionable finding: **the binding constraint is not the ordering method, it is
coverage.** The architecture cannot produce a short list because it cannot assess most
exclusions on most pairs.

### What would have to change before this task can be re-attempted

1. **CNA coverage.** `cna_gate` is inert at 99.72% not-assessed. COSMIC CNA covers 997 of
   2,145 models but only 0.31% of (gene, line) pairs. A gene-complete CNA source — e.g. the
   unused `data/DepMap_Chronos/CCLEGeneCopyNumber20Q2.hdf5` (202 MB, currently an orphan
   input) — would make this check real. **This is the single highest-value change** and it is
   already on disk.
2. **A second discriminating check that is near-universally assessed.** Today only the RNA
   gate qualifies. Until a second one exists, "count of independent confirmations" has a
   maximum meaningful value of about 1.
3. **The evaluation denominator fix.** Independently required; nothing here can be validated
   until untested pairs stop being recorded as confirmed negatives.

Note that (1) and (2) are the same requirement seen twice.

---

## 8. Decisions this leaves you

Standing rule 7 says stop, so `RANKING_SPEC.md`, `RANKING_EVALUATION_PLAN.md` and
`RANKING_DECISIONS.md` were **not** produced. Four questions remain, and none of them is
answerable by more measurement of the count.

**Q1. Wire in a gene-complete CNA source?** `CCLEGeneCopyNumber20Q2.hdf5` is on disk and
unused. It would move `cna_gate` from 99.72% not-assessed to near-complete, and is the only
identified route to a second real check. Cost: a build step and a re-derivation of the CNA
thresholds against ploidy (which `BUILD_SPEC.md` C6 already requires for other reasons —
one correction serves both).

**Q2. Accept ~1,200 lines in 2 tiers as the product, or narrow the gate?** The current
exclusion set removes ~22% of the panel. Narrowing it (stricter thresholds, more exclusions)
would produce a shorter list but needs a swept, held-out justification per standing rule 4 —
and the gate's measured pAUC of 0.596 is the evidence budget being spent.

**Q3. Do the protein bands proceed independently?** Part C was not attempted because
standing rule 7 stopped the task. The protein tier is Level 3 and does not depend on the
count, so it *could* be specified separately. But note the coverage arithmetic from this
report: `protein` is not assessed for 80.97% of retained pairs, so a protein band would order
at most a fifth of each gene's retained set. Worth doing only if Q1 also proceeds.

**Q4. Alteration direction** — still open from `DECISIONS_REQUIRED.md` Q1, and now lower
stakes: at 1.56% of pairs failing, the tier is inert either way for the typical gene.

---

## 9. The honest framing

The brief's appendix says a shortlist of twenty with a stated basis and honest ties beats an
ordered list of fourteen hundred whose tail carries no information. That is right, and this
measurement says the project cannot currently produce either one — not because the ordering
method is wrong, but because **on most (gene, line) pairs there is only one thing measured**,
and one measurement cannot be corroborated by itself.

Part A was designed to find that out cheaply. It did, in about 460,000 pairs and no
production code, and the answer is more useful than a built ranking would have been: the next
piece of work is a coverage problem with a 202 MB file already sitting in `data/`.

---

## 10. The re-run, with the copy-number lever pulled

§7 named `CCLEGeneCopyNumber20Q2.hdf5` as "the single highest-value change". It was built in
and Part A re-run. This section records what happened.

### 10.1 What changed

`src/pipeline/build_cn_layer_complete.py` — additive, does not modify `build_cna_layer.py`.

| | COSMIC CNA (shipped) | Gene-complete CN (new) |
|---|---|---|
| models × genes | 969 × 15,545 | 908 × 18,276 |
| **(gene, line) pairs** | **117,540** | **16,594,608** (141×) |
| missing values in source | — | **zero** |
| assessed share of gate-retained pairs | **0.28%** | **52.29%** |
| ploidy correction | assumes diploidy (gap 8) | **not needed** — already relative to each line's own ploidy |

The ploidy point is worth keeping: the matrix is `log2(relative_CN + 1)` with neutral = 1.0,
verified by Spearman(per-line mean, GDSC ploidy) = +0.115 with per-line mean sd 0.020 over a
panel spanning ploidy 1.52–5.40. **Gap 8 is closed for this layer**, and the 30%-of-panel
missing-ploidy problem does not arise.

### 10.2 The pre-specified rule and the outcome

Thresholds fixed in `docs/RANKING_PRESPEC.md` before any re-run number existed. Primary
statistic: raw count, ≥2 assessed checks, gate excluded, valid genes, within the retained set.
Paired with the first run — identical 300 genes, seed 42, identical other checks.

| Criterion | Adopt if | Abandon if | **Before** | **After (rel<0.5)** | Moved |
|---|---|---|---|---|---|
| A1 coverage confound \|ρ\| | ≤ 0.20 | ≥ 0.40 | +0.182 | **+0.531** | worse |
| A2a largest tied group | ≤ 0.50 | ≥ 0.65 | 0.587 | **0.654** | worse |
| A2b distinct values | ≥ 4 | ≤ 2 | 4 | 4 | — |
| A3 between-gene ρ | ≤ 0.25 | ≥ 0.40 | +0.092 | **+0.475** | worse |
| **Verdict** | | | INCONCLUSIVE | **ABANDON** | |

Sensitivity across the three pre-specified CN thresholds — all **ABANDON**, so the prespec's
"verdicts differ ⇒ inconclusive" clause does not apply:

| CN deletion cut | A1 | A2a | A3 | Verdict |
|---|---|---|---|---|
| rel < 0.25 | +0.544 | 0.657 | +0.491 | ABANDON |
| **rel < 0.5 (primary)** | **+0.531** | **0.654** | **+0.475** | **ABANDON** |
| rel < 0.75 | +0.479 | 0.624 | +0.372 | ABANDON (on A1 alone) |

Usable pairs did rise — 33.9% → 55.4% of retained pairs clear the ≥2-assessed bar, and median
orderable lines per gene rose from 544 to 1,011. The count got *bigger*. It did not get
*better*.

### 10.3 Why more data made it worse

Two measured mechanisms:

1. **CN availability is a line-level indicator, not gene-level evidence.** 840 of 1,543 lines
   carry CN data, and every covered line has data for *all* sampled genes (min = max = 288).
   Its presence therefore encodes only "is this line in the DepMap release", and that
   correlates with other-evidence coverage at **ρ = +0.242, p = 4.9e-22**.
2. **The check almost never fails.** Among assessed retained pairs it passes **98.49%**.
   Median 2 deleted lines per gene out of 686 assessed; **104 of 300 genes have zero deleted
   retained lines**.

Adding a near-constant "pass" present for a specific half of the panel is arithmetically
adding a coverage indicator. That is why A1 rose from +0.182 to +0.531.

Per-check state on retained pairs after the change:

| Check | assessed | pass (of assessed) | fail | discrimination |
|---|---|---|---|---|
| `alteration` | 98.77% | 98.36% | 1.64% | 0.016 |
| `cna_gate` (new, rel<0.5) | 52.29% | 98.49% | 1.51% | 0.015 |
| `geo_corrob` | 32.08% | 68.51% | 31.49% | 0.315 |
| `protein` | 19.03% | 50.03% | 49.97% | 0.500 |

The generalisation — **coverage and discrimination are inversely related, ρ = −0.80** — is
written up as a standalone result in `docs/COVERAGE_RESULT.md`.

### 10.4 Pre-registered predictions, scored

Three of four were wrong, all in the same direction.

| # | Prediction | Outcome |
|---|---|---|
| 1 | assessed-fraction rises above 55% | **FAILED** — 52.29% |
| 2 | A1 improves but may not clear 0.20 | **FAILED** — worsened to +0.531 |
| 3 | A2a most likely to fail | partly right — it failed, but so did A1 and A3 |
| 4 | A3 improves most | **FAILED badly** — +0.092 → +0.475 |

### 10.5 Where this leaves it

The lever named in §7 has been pulled and the count is not viable. The conclusion is stronger
than before, not weaker: it is no longer "there is not enough data to tell", it is **"the best
available additional data makes it worse, for a reason that generalises"**. Q1 of §8 is now
answered — wiring in the CN source was worth doing, but not for the confirmation count.

The CN layer is still worth keeping. At `rel < 0.75` its discrimination is 0.117 and 210 of
300 genes gain ≥10 deleted retained lines, which makes it a usable **exclusion** signal —
the role CNA deletion already had (Fisher ratio 0.380 on 69 genes), now on 141× the pairs and
without the diploidy assumption. It is a gate input, not a confirmation.

**What would still change the answer:** raising `protein` (discrimination 0.500) from 19%
coverage, or `geo_corrob` (0.315) from 32%. Both already clear the discrimination bar and
fail only on coverage. No further ordering method should be attempted before one of them
does.
