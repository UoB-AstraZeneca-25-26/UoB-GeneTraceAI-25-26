# Proteomics — what is measured, from what, and what it actually does

Defence notes for the protein layer. Every number is measured on this repository
(2026-08-11) and cited to the file that produces it.

Companion to [ALTERATION_DEFENCE.md](ALTERATION_DEFENCE.md), same structure.
Formulae in [MATH_REFERENCE.md](MATH_REFERENCE.md); the retired experiments in
[MATH_TRIAL_AND_ERROR.md](MATH_TRIAL_AND_ERROR.md) §C.18.

---

## 1. What exactly is being measured

**The unit is the (gene, cell line) pair**, same as alteration. The claim is:

> in cell line *L*, the protein product of gene *G* was present at this relative
> abundance

Three things that are *not* being measured, and are easy to conflate:

- **Not absolute concentration.** Both platforms produce relative abundance —
  ratios against a reference channel (TMT) or against the run (DIA). A protein
  cannot be compared across genes, only across cell lines within a gene.
- **Not protein activity.** Abundance says nothing about phosphorylation state,
  localisation, or whether the protein is in a functional complex.
- **Not the same quantity as mRNA.** Measured on this repo: per-gene
  mRNA–protein Spearman ρ has **median 0.459** across 8,561 genes
  (`validation/prepared/expr_protein_correlation_per_gene.parquet`), range
  −0.211 to 0.894. The layers agree moderately and disagree substantially.

---

## 2. Which datasets, which columns, which values

### Resources actually consumed

| resource | file / table | shape | what it is |
|---|---|---|---|
| CCLE / Gygi proteomics | `cleaned_track_data/proteomics.parquet` | 375 × 10,129 proteins | **TMT** 10-plex isobaric labelling, MS2/MS3 |
| ProCan proteomics | `outputs/procan_proteomics.parquet` | 952 × 8,453 proteins | **DIA-MS** (data-independent acquisition) |
| UniProt → symbol (CCLE) | `main.protein_map` | 12,558 | accession resolution |
| UniProt → symbol (ProCan) | `main.procan_protein_map` | 8,453 | accession resolution |
| mRNA–protein ρ per gene | `validation/prepared/expr_protein_correlation_per_gene.parquet` | 8,561 | the layer-agreement estimate |

Columns are **lowercase UniProt accessions**, not gene symbols — the join to the
gene axis runs through `uniprot_ids` in the gene dimension.

### Coverage — the headline reach

From `HARMONISATION_V2.md` §4.1, over the 2,145-line hub:

| | lines | % of hub |
|---|---:|---:|
| CCLE/Gygi only | 84 | 3.9% |
| ProCan only | 661 | 30.8% |
| **both platforms** | **291** | 13.6% |
| **either** | **1,036** | **48.3%** |

### Reach into the score

`core_score.parquet` — 29,781,274 pairs:

| | value |
|---|---:|
| pairs with `n_layers = 2` (RNA + protein) | **5,277,064 (17.7%)** |
| pairs with `n_layers = 1` (RNA only) | 24,504,210 (82.3%) |
| genes with any 2-layer pair | 10,823 of 19,177 |
| lines with any 2-layer pair | **769 of 1,746** |

So the protein layer touches under a fifth of the scored space.

---

## 3. "How can you tell a protein was detected or not?"

**Detected** means a non-null value in the abundance matrix for that
(accession, line) — the peptide was identified and quantified in that MS run.

**Absent** is where this layer has the same defect as alteration, in a sharper
form. A missing protein value conflates at least four situations:

1. the protein is genuinely not expressed
2. it is expressed below the instrument's limit of detection
3. it was not selected for fragmentation (**TMT/DDA only** — stochastic
   precursor selection means the same protein can be quantified in one run and
   missed in the next, from identical material)
4. no tryptic peptide of that protein is detectable — some proteins are
   effectively invisible to shotgun MS regardless of abundance

Case 3 is specific to the CCLE/Gygi arm and is not a property of the sample at
all — it is a property of the run. Case 4 is a fixed property of the protein's
sequence. **Neither is biology, and neither is recoverable from the matrix.**

The abundance-stratified agreement below is direct evidence that this matters:
low-abundance proteins are exactly where the two platforms stop agreeing.

---

## 4. "How are you managing two proteomics platforms?"

**They are stored as separate layers and never concatenated.** That constraint
was imposed after measurement, not by preference — `test_run_proteomics_platform_overlap.py`
compared them on the 291 lines and 6,121 proteins they share, *before* any merge
was permitted, precisely because Stage 2 percentile-ranks each protein across
cell lines and a naive merge would let a protein look "high" merely because it
was measured on the more sensitive instrument.

### What the comparison found

| measure | value |
|---|---|
| per-protein ρ, median | **0.373** |
| proteins above ρ = 0.5 | 26.7% |
| proteins with negative ρ | 4.3% |
| platform tiers | 1,632 consistent (26.7%) · 2,243 cautious (36.6%) · **2,246 conflicting (36.7%)** |
| **per-line agreement, median cross-protein ρ** | **0.097** |
| **lines classified "trusted"** | **0 of 291** |

Two results deserve to be said out loud:

- **More proteins are classified conflicting (36.7%) than consistent (26.7%).**
- **Not one of the 291 dual-measured lines passed the trust threshold.** Median
  cross-protein correlation within a line is 0.097 — close to no relationship.

### Agreement depends on abundance

| abundance quartile | median ρ | % consistent |
|---|---:|---:|
| Q1 (low) | 0.235 | 5.3% |
| Q2 | 0.322 | 11.8% |
| Q3 | 0.430 | 33.0% |
| Q4 (high) | **0.511** | **52.9%** |

Agreement more than doubles from the lowest to the highest abundance quartile.
The two platforms substantially agree about abundant proteins and substantially
disagree about scarce ones — which is the expected signature of detection-limit
effects, and is why §3's missing-value ambiguity is not a pedantic point.

Shipped verdict: `merge_with_caution` — "merge only with per-dataset
normalisation, and report the batch effect as a caveat".

---

## 5. On what factor does protein enter the score?

```
core_score = w_E·E + w_P·P + w_C·C      with w = normalise(1/(1 + mean|ρ|))
```

where E, P, C are **percentile ranks within gene across cell lines**, not raw
abundances. Two- and one-layer subsets renormalise over the layers present.

**For the two-layer case the weighting collapses to a simple mean.** By symmetry
`mean ρ_E = mean ρ_P`, so `w_E = w_P = ½` identically, for every ρ:

```
core_score = ½·E + ½·P
```

This is worth being able to state precisely, because it is a common trap: the
correlation term is in the code and has **no effect whatsoever** on any
production two-layer score. That is not a defect — (½, ½) is the GLS optimum for
two exchangeable layers at any ρ — but the correlation weighting cannot be
claimed as the thing that improved the score.

**Two layers are not two layers' worth of evidence.** Kish's effective sample
size on the measured correlations (ρ_EP = 0.46):

```
n_eff = n / (1 + (n−1)·ρ̄)  =  2 / 1.46  =  1.37
```

So `n_layers = 2` carries **1.37** independent layers, not 2. `n_layers` is used
throughout the pipeline as a stratification key and a confidence-tier input, and
this gives it a quantitative meaning rather than a bookkeeping one.

---

## 6. What proteomics actually contributes — measured

This is the section to lead with, because it contains a genuine tension.

### 6.1 The retired verdict: "protein is a net negative"

`test_run_stratum_centring.py`, reported in §C.18. The 2-layer stratum genuinely
sits higher (median shift 0.210, positive in 86.3% of 1,562 genes), but
correcting for that made things worse, and on **global AUROC** RNA alone beat
every combination tried:

| arm | GDSC AUROC (n=110) | Chronos AUROC (n=1,030) |
|---|---:|---:|
| **RNA alone** | **0.5630** | **0.5193** |
| production (2-layer) | 0.5544 | 0.5151 |
| Stouffer Z | 0.5560 | 0.5174 |
| Z centred | 0.5473 | 0.5088 |

Verdict recorded: `centring_does_not_recover__protein_is_a_net_negative`.

### 6.2 The same test on pAUC says the opposite — and it was not carried forward

The same results file also reports **partial AUC**, which measures only the top
of the ranking. Both tables are **medians of the per-gene statistic**, not pooled
values, and `cur` is labelled in the script itself
(`test_run_stratum_centring.py:334`) as `"current 0.5E+0.5P (equal weight)"` —
i.e. the production two-layer formula:

| arm | GDSC pAUC | Chronos pAUC |
|---|---:|---:|
| RNA alone | 0.6415 | 0.6133 |
| **production (2-layer)** | **0.6649** | **0.6304** |
| Stouffer Z | 0.6289 | 0.6047 |
| Z centred | 0.6002 | 0.5899 |

**On pAUC the production two-layer score beats RNA alone on both arms** — +0.023
(GDSC) and +0.017 (Chronos). The Stouffer and centred variants still lose on
both metrics, so §C.18's rejection of *those* stands. But the blanket phrase
"protein is a net negative" is carried in
[STATE_REPORT.md](STATE_REPORT.md), [RISK_REGISTER.md](RISK_REGISTER.md) and
[REBUILD_RUNBOOK.md](REBUILD_RUNBOOK.md) as though it applied to the production
score, and on the metric that matches the product it does not.

**Why this matters for the defence.** The product is a top-N recommender — it
shows a user the highest-ranked cell lines for a gene. pAUC is the metric that
matches that use; global AUROC weights the entire ranking including the tail
nobody sees. So the defensible statement is:

> The protein layer costs ~0.004–0.009 AUROC globally and gains ~0.017–0.023
> pAUC at the top of the ranking. For a top-N product, that trade favours
> keeping it.

That is a substantially stronger position than the one currently written down,
and it is honest. It should be verified and then propagated to the three
documents that carry the flattened claim.

### 6.3 The tier rule (Gap 14)

`build_full_predictions.py:151` and `05_confidence_tiers.ipynb` cell 5:

```
(regime_source == "measured") & (class == "abundance_tracking") & (n_layers == 2)  →  "high"
```

Having a protein layer is what promotes a row to **high** confidence. Under §6.1
alone this is indefensible — the layer was measured as harmful. Under §6.2 it is
defensible on pAUC but still overstated, because `n_layers = 2` means 1.37
effective layers, not 2.

Honest position: **the rule is directionally justified by pAUC but has never
been calibrated.** Nobody has checked whether "high"-tier rows are actually more
often correct than "moderate"-tier ones. That is the missing test, and it is a
cheap one.

---

## 7. KPIs, and what each is actually sensitive to

| KPI | value | source | sensitive to |
|---|---:|---|---|
| protein coverage (lines) | 1,036 / 2,145 = **48.3%** | `HARMONISATION_V2.md` §4.1 | ProCan's 661-line contribution |
| dual-platform lines | 291 (13.6%) | same | the only cross-check available |
| 2-layer pairs in score | 5,277,064 / 29,781,274 = **17.7%** | `core_score.parquet` | UniProt→ENSG mapping |
| cross-platform ρ (median) | **0.373** | `test_run_proteomics_platform_overlap_results.json` | abundance quartile (0.235 → 0.511) |
| lines passing platform trust | **0 / 291** | same | the trust threshold itself |
| mRNA–protein ρ (median) | **0.459** | `expr_protein_correlation_per_gene.parquet` | literature comparator is 0.58 |
| effective layers (Kish) | **1.37** | `MATH_REFERENCE.md` §0.3 | ρ_EP = 0.46 |
| ΔAUROC vs RNA alone | −0.0086 GDSC / −0.0042 Chronos | `test_run_stratum_centring_results.json` | metric choice — see §6.2 |
| ΔpAUC vs RNA alone | **+0.0234 GDSC / +0.0171 Chronos** | same | metric choice — see §6.2 |

Note the mRNA–protein correlation: **0.459 measured here versus ρ ≈ 0.58 quoted
from Nusinow 2020** in the Stage 2 notebook. Our layers agree *less* than the
literature figure the weighting scheme was parameterised with. Since the two-layer
weights collapse to (½, ½) regardless of ρ, this does not change any score — but
it does mean the `n_eff = 1.37` figure is optimistic. At ρ = 0.459 the arithmetic
is nearly identical; the caveat is that ρ_EP was taken from literature on native
scales while the layers combined are percentile-transformed.

---

## 8. The two hard questions

### (a) "mRNA and protein disagree for this gene — which do you believe?"

**Currently: neither, and that is the honest answer.** The production score
averages them, `½E + ½P`, with no disagreement detection. A gene where RNA is at
the 90th percentile and protein at the 10th produces the same core_score (0.50)
as a gene where both sit at the 50th — two completely different evidential
situations collapsed to one number.

Measured, this is not rare: per-gene mRNA–protein ρ has median 0.459 and reaches
**−0.211** at the low end. Genes with negative ρ are ones where the layers
systematically contradict each other.

What should exist and does not: a **disagreement flag** on the pair, so that
`n_layers = 2` with concordant layers and `n_layers = 2` with contradictory
layers are distinguishable. Today they are not, and the tier rule promotes both
to "high" identically. This is the proteomics analogue of the three-state
`detected` field — the architecture cannot currently represent "measured twice,
got conflicting answers".

### (b) "Which platform is right?"

**Unanswerable with what is held, and the data says so unusually clearly.**
Zero of 291 dual-measured lines passed the agreement threshold; median
within-line cross-platform ρ is 0.097. There is no third platform to arbitrate.

What can be said: agreement is a function of abundance (0.235 → 0.511 across
quartiles), so for **abundant proteins** the platforms broadly concur and either
can be trusted; for **scarce proteins** they do not, and neither should be. The
`protein_platform_tier.parquet` sidecar already carries the per-protein tier —
it is computed and available but is not consumed by the score. Wiring it in as a
reliability weight is the obvious next step and would be a defensible use of an
existing measurement.

---

## 9. One-line answers

- **What are you measuring?** Relative protein abundance for a (gene, cell line)
  pair — not concentration, not activity.
- **What datasets?** CCLE/Gygi TMT (375 lines × 10,129 proteins) and ProCan
  DIA-MS (952 × 8,453). Stored separately, never concatenated.
- **How do you know a protein was detected?** A non-null value in the abundance
  matrix for that accession and line.
- **How do you know it was absent?** You don't — missing conflates not-expressed,
  below-detection, not-selected-for-fragmentation, and not-detectable-by-MS.
- **Which columns?** Lowercase UniProt accessions, joined to the gene axis via
  `uniprot_ids`.
- **How are the two platforms combined?** They are not. Measured first
  (median ρ 0.373, 0/291 lines trusted), and kept as separate layers on that basis.
- **How does protein enter the score?** As a percentile rank, `core_score = ½E + ½P`.
  The ρ-weighting in the code has no effect on any two-layer score.
- **How much evidence is two layers?** 1.37 effective layers (Kish), not 2.
- **What does it contribute?** −0.004 to −0.009 AUROC globally, **+0.017 to
  +0.023 pAUC** at the top of the ranking. For a top-N product the second number
  is the relevant one.
- **What is the known defect?** Layer disagreement is invisible — concordant and
  contradictory 2-layer pairs are scored and tiered identically.

---

## 10. Say these unprompted

1. **"The protein layer reaches 17.7% of scored pairs and 769 of 1,746 lines."**
   Coverage is the first thing to concede — this is not a layer that applies
   broadly, and any claim about it generalises only to the half of the panel that
   has proteomics at all.
2. **"Zero of 291 dual-measured lines passed cross-platform agreement."** Median
   within-line ρ is 0.097, and more proteins are classified conflicting (36.7%)
   than consistent (26.7%). Agreement is abundance-dependent — 0.235 in the
   lowest quartile, 0.511 in the highest. That is why the platforms are never
   merged.
3. **"'Protein is a net negative' is true on AUROC and false on pAUC."** The
   layer costs 0.004–0.009 AUROC globally but gains 0.017–0.023 pAUC at the top
   of the ranking, which is the region a top-N recommender actually serves. The
   flattened claim in `STATE_REPORT.md` should be corrected, and the confidence
   tier rule it undermines is more defensible than currently documented — though
   still uncalibrated.
4. **"Two layers are 1.37 layers."** Kish effective sample size on ρ_EP = 0.46.
   The pipeline treats `n_layers` as a count; it is not one.
