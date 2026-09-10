# `gate_audit/` — the abundance gate, re-examined

A self-contained folder. Nothing here writes to `src/pipeline/outputs/`, to the
DuckDB catalog, or to any production table. Every script is read-only with
respect to the pipeline and writes only into `gate_audit/outputs/`.

It exists because the gate result — the project's only solid finding — carried
a series of open challenges, and answering them inside `src/pipeline/` would have
meant nine more `test_run_*.py` files in a directory that already has twenty-two.

---

## Run it

From the repo root:

```bash
python gate_audit/run_all.py
```

Individually, in order:

```bash
python gate_audit/00_build_panel.py --n-genes 20000
```

```bash
python gate_audit/01_tie_mass.py
```

```bash
python gate_audit/02_gate_valid_genes.py
```

```bash
python gate_audit/03_denominator.py
```

```bash
python gate_audit/04_cluster_bootstrap.py --n-boot 1000 --n-genes 3900
```

```bash
python gate_audit/05_region_sweep.py
```

```bash
python gate_audit/06_screen_replication.py
```

```bash
python gate_audit/07_independent_labels.py
```

```bash
python gate_audit/08_release_refresh.py
```

```bash
python gate_audit/09_perturbation_types.py
```

`00` is the slow step (~15 min) and every later step reads its output, so
`run_all.py --skip-panel` reuses it. `run_all.py --quick` runs a smaller version
end to end.

---

## Architecture

| File | What it is |
|---|---|
| `common.py` | Every loader, every transform, every metric. The single definition of "the gate". |
| `00_build_panel.py` | Builds `outputs/panel.parquet` — one row per (gene, line) with score, label, lineage — and `outputs/genes.parquet` with per-gene validity and tie structure. |
| `01_tie_mass.py` | Tie mass and tie-handling in the expression percentile. |
| `02_gate_valid_genes.py` | The gate recomputed on valid genes only. |
| `03_denominator.py` | The profiling confound, in its direct form. **Read this one first if you read only one.** |
| `04_cluster_bootstrap.py` | Lineage and two-way cluster bootstrap for an honest interval. |
| `05_region_sweep.py` | Sweeps the three post-hoc free parameters. |
| `06_screen_replication.py` | Broad vs Sanger on the real Chronos release. |
| `07_independent_labels.py` | GDSC2 as an independent perturbation type; coverage per label per gene. |
| `08_release_refresh.py` | Rebuild on DepMap 24Q4 (1,103 usable lines vs 754). |
| `09_perturbation_types.py` | CRISPR vs RNAi vs drug at matched prevalence. The tiebreak. |
| `run_all.py` | Runs the above in order. |

### The DepMap release on disk

`data/DepMap_Chronos/` is the Zenodo **14067047** benchmark bundle, built on
**DepMap 20Q2** (`DepMapSampleInfo20Q2.csv`, 907 lines). Achilles Chronos in it
carries 767 lines, 754 after intersecting with expression. Current public DepMap
releases carry ~1,100–1,200 lines in `CRISPRGeneEffect.csv`, so a refresh would
plausibly take 754 → ~1,000. That is a power gain, not a fix: it does not touch
the audit-3 denominator finding and it does not change the GDSC2 null. If the
CRISPR label is refreshed, the expression matrix should be refreshed from the
same release, or model-ID drift between vintages becomes a new confound.

The bundle also contains MAGeCK, CERES, BAGEL2 and Naive processings of the
*same* screens. Those are not independent labels — same perturbation, same reads,
different estimator — and are not used here.

The point of `common.py` is that the original diagnostics each re-loaded and
re-derived everything independently, which is five chances to disagree about what
is being measured. Here the panel is built once and the six audits argue about
statistics, not inputs.

**The label is the real DepMap Chronos release**
(`data/DepMap_Chronos/GeneFitnessEffect_Chronos_Achilles.hdf5`), not
`validation/prepared/chronos_long.parquet`, which is Project Score scaled BF
negated. `common.load_chronos()` runs the reference-set scale check and **aborts**
rather than warns — that check's absence is what let the misnamed file pass as
Chronos project-wide.

---

## What the audits found

Panel: 3,866 genes × 754 lines = 2.9M scored (gene, line) pairs. 3,514 genes
valid, 352 invalid.

### 03 — the denominator. This is the finding.

`src/pipeline/test_run_abundance_gate.py` builds each gene's line set from the
expression and proteomics indices, then labels positives by set membership:

```python
idx = sorted(set(zE.index) | set(zP.index))
d["pos"] = d.index.isin(pos_by_gene[g])
```

`pos_by_gene[g]` contains only CRISPR-screened lines. Any line with expression
but no screen therefore gets `pos = False` — counted as a confirmed
non-dependency rather than excluded as unmeasured.

| | lines |
|---|---|
| expression | 1,479 |
| Chronos (Achilles) | 767 |
| overlap | 754 |
| **never screened** | **725 (49% of the denominator)** |

Those guaranteed-negative lines are **not** spread evenly across the score. They
concentrate at the bottom: decile 1 is 55.7% unscreened against 49.0% overall.
Decile 1 is enriched for lines that *cannot* be positive.

| construction | decile-1 ratio | depletion |
|---|---|---|
| as published | 0.8322 | 16.8% |
| screened lines only | 0.9875 | 1.3% |
| **placebo — labels randomised, no biology at all** | **0.8350** | **16.5%** |

The placebo keeps the construction and discards the dependency information
entirely. It reproduces **98% of the published depletion**. Gate-region pAUC
falls from 0.6125 to 0.5228 under the same correction.

The published 0.839 is reproduced exactly by the original construction (0.8322
here), so this is not a disagreement about implementation — it is the same
number, and it is an artefact of the denominator.

This is the profiling confound, one step closer to the outcome than the protein
version: not "profiled lines differ biologically" but "unprofiled lines were
scored as negatives".

> **A correction to how this was framed.** The protein-detection result in
> `test_run_gate_mechanism.py` reads the other way round from the summary that
> prompted this audit. Ratio 1.07 with Fisher p = 0.246 **at the floor** means
> profiled and unprofiled lines are *indistinguishable there* — the floor is
> clean of that confound. Ratio 1.62 with p = 0 **in the mid stratum** means the
> confound lives in the middle. The script's own verdict is
> `no_floor_specific_effect__profiling_confound_only`. So "report the gate as a
> mid-range result and exclude the floor" inverts it: the mid-stratum is the
> contaminated one. The denominator finding above supersedes both readings.

### 01 — tie mass: not the problem

Two tie decisions, in series, and they are different:

1. **The percentile.** `vdw()` uses `rank(method="min")`. Ties *collapse* to one
   shared z at the block minimum — so a floor block sits at the bottom of the
   score, not smeared through the middle.
2. **The decile.** `rankdata(method="ordinal")` — ties *spread*, in cell-line
   array order.

So the answer to "average rank or arbitrary tiebreak" is: both, at different
stages. The score is not contaminated; the stratification is.

| | median | p90 |
|---|---|---|
| tie mass (share of lines tied) | 0.182 | 0.798 |
| — valid genes | 0.162 | |
| — invalid genes | **0.941** | |

Tie mass is concentrated almost entirely in the invalid genes, which is what the
silence guard already says. And it does not move the answer:

| tie rule | decile-1 ratio |
|---|---|
| ordinal (production) | 0.9799 |
| average | 0.9800 |
| drop ties | 0.9819 |

Spread 0.0021. Paired shift `average − ordinal`: +0.0000, p = 0.45. **The
headline does not depend on array order.** Partial AUROC was never exposed to
this anyway — equal scores collapse into one ROC step.

### 02 — valid genes: also not the problem

| | ALL | VALID | INVALID |
|---|---|---|---|
| genes | 3,866 | 3,514 | 352 |
| decile-1 ratio | 0.9799 | 0.9835 | 0.9031 |
| gate pAUC | 0.5232 | 0.5254 | 0.5118 |

Restricting to the valid genes moves the decile-1 ratio by +0.0036 and pAUC by
+0.0021. **This doubt is closed.** The degenerate distributions the silence guard
flags are not what the gate was measuring. Report on valid genes anyway — it
costs nothing and removes the question.

Note the small absolute magnitudes are *not* about validity. They come from 03.

### 04 — cluster bootstrap: the effect survives, the p-value does not

Do not quote p = 2.8e-76. It treats 3,866 genes × 754 lines as independent units
when lines cluster by lineage (project gaps 7 and 17) and genes cluster by
co-expression.

1,000 replicates, 45 lineage clusters, 3,866 genes:

| | SE | 95% CI | bootstrap p | design effect |
|---|---|---|---|---|
| nominal | 0.0043 | — | 1.4e-10 | 1.0× |
| lineage cluster bootstrap | 0.0092 | [0.513, 0.548] | 0.000 | **2.14×** |
| two-way (lineage × gene) | 0.0103 | [0.510, 0.551] | **0.002** | **2.39×** |

The nominal standard error understates the real one by a factor of 2.4. Also
worth stating plainly: with 45 lineage clusters the **Kish effective sample size
is 8.7 clusters**, not 754 lines — line-level claims have an effective *n* an
order of magnitude below the nominal one.

Quote this:

> gate-region pAUC = 0.523, 95% CI [0.510, 0.551], p = 0.002 (two-way cluster
> bootstrap over 45 lineage clusters and 3,866 genes, 1,000 replicates)

The interval excludes chance, so the effect survives clustering — but note where
it lands. The nominal p was 1.4e-10; the honest one is 0.002. That is still a
real effect, and it is defensible, but it is not overwhelming and it should not
be described as such.

### 05 — the post-hoc region: a plateau, not a spike

`FPR ∈ [0.8, 1.0]` was chosen after full AUROC came back flat at 0.52. That is
post-hoc selection and this sweep does not undo it — it shows whether the choice
was load-bearing.

| lower bound | 0.70 | 0.75 | 0.80 | 0.85 | 0.90 |
|---|---|---|---|---|---|
| pAUC (valid) | 0.5189 | 0.5200 | 0.5254 | 0.5301 | 0.5351 |

Range across the full `[lo, hi]` grid: 0.0177. **0.8 is not doing any work.**
Bottom-bin depletion is likewise stable across 4–20 bins (0.9915 → 0.9777, sign
and story unchanged). NPV uplift over base decays smoothly from +0.017 at a 1%
cut to +0.003 at 25%.

Report the range, not the point — and state the pre-specification defence
explicitly rather than letting a panel find it: the region was fixed before the
screen replication in 06, and screen identity was not used to choose it.

### 06 — Broad vs Sanger on real Chronos

**Correcting the record:** the gate *has* already been run on real Chronos —
`test_run_abundance_gate_chronos_results.json` carries `label_fpr = 0.00928` and
the reference-set medians, which only the real-Chronos branch produces. What had
not been done is the two-screen comparison under the corrections above.

| | Broad | Sanger | paired p |
|---|---|---|---|
| decile-1 ratio | 0.9927 | 0.9960 | 0.338 |
| gate pAUC | 0.5077 | 0.5051 | 0.14 |
| decile-1 ∩ EXPRESSED | 0.9948 | 0.9960 | 0.376 |

Per-gene Spearman between screens: 0.32 (decile-1), 0.37 (pAUC).

**A selection effect that must be stated.** Sanger's panel is 317 lines against
Broad's 767. Requiring ≥10 dependent and ≥10 non-dependent lines in *both*
screens selects genes essential in a large share of lines: the paired set's
median base rate is 0.435 against the panel's 0.118. A gate cannot act on a gene
essential in 44% of lines, so the replication must be read by prevalence band:

| band | genes | Broad d1 | Sanger d1 | Broad pAUC | Sanger pAUC |
|---|---|---|---|---|---|
| <2% | 67 | 0.9019 | 0.9731 | 0.5077 | 0.5200 |
| 2–5% | 117 | 0.9921 | 0.9222 | 0.5072 | 0.5688 |
| **5–10%** | **115** | **0.8858** | **0.8903** | **0.5677** | **0.5944** |
| 10–25% | 233 | 0.9771 | 0.9486 | 0.5047 | 0.5347 |
| >25% | 992 | 0.9962 | 1.0082 | 0.5067 | 0.4809 |

**The 5–10% prevalence band is where the gate actually lives, and it replicates
across two independent screens** — 11% depletion and pAUC ~0.58 in both. That is
a real, narrow, replicated result.

On the tautology question — "isn't this just saying you can't knock out a gene
that isn't expressed?" — the intersection of bottom decile with
rule-says-EXPRESSED still depletes in Broad (0.9867, p = 1e-12) but not
significantly in Sanger (0.9960, p = 0.08). The percentile is not *merely* an
encoded detection boolean, but on the corrected denominator that defence is much
thinner than the 0.810–0.827 previously reported.

### 07 — biology or assay artefact? The gate is CRISPR-specific.

Broad and Sanger are the same *perturbation type*. Their agreement rules out
screen-specific noise, not a property of CRISPR knockout as an assay. GDSC2
(small-molecule dose-response, different institution, different panel) breaks
that — and it is already in the repo, verified here as 100% GDSC2 with no GDSC1
contamination.

Comparison made on the 140 genes present in GDSC2 *and* both CRISPR screens, each
label on its own denominator:

| label | genes scored | base rate | decile-1 ratio | p | gate pAUC | p |
|---|---|---|---|---|---|---|
| Achilles (CRISPR, Broad) | 63 | 0.117 | **0.5851** | 3.6e-07 | **0.7132** | 3.9e-07 |
| Score (CRISPR, Sanger) | 34 | 0.173 | **0.5369** | 2.9e-05 | **0.8303** | 2.0e-04 |
| GDSC2 (drug) | 140 | 0.298 | 0.9960 | 0.304 | 0.5012 | 0.476 |

**GDSC2's `sensitive` flag has a median base rate of 0.30, so every CRISPR band
below 10% is empty for it** — the band where the gate lives could not be compared
at all. A null there would have been an artefact of the binarisation. So GDSC2
was re-binarised per gene on continuous response (minimum AUC across the gene's
drugs), taking the most-sensitive *q* fraction as positive:

| q | base rate | decile-1 ratio | p | gate pAUC | p |
|---|---|---|---|---|---|
| 0.02 | 0.020 | 0.9901 | 0.179 | 0.5365 | 0.658 |
| 0.05 | 0.050 | 0.9037 | 0.761 | 0.5254 | 0.588 |
| 0.10 | 0.100 | 1.0678 | 0.229 | 0.4802 | 0.437 |
| 0.25 | 0.250 | 1.0278 | 0.198 | 0.4875 | 0.417 |

Flat at every matched prevalence. Against Achilles on the *same 140 genes* at the
same prevalence: 0.000 / 0.602 / 0.311. Paired within-gene, Spearman between
per-gene gate strength under CRISPR and under drug is **−0.09** (decile-1) and
**−0.14** (pAUC) — not attenuated correlation, no correlation.

**The gate does not reproduce on a different perturbation type.** Given GDSC2's
known causal distance — a drug must reach the cell, hit the annotated target
rather than an off-target, and killing must follow; clinical kinase inhibitors
are extensively polypharmacological (Klaeger et al. 2017) — this is not a
refutation. But it means the claim cannot be stated as a general
abundance–dependency law. It is a statement about CRISPR-derived dependency until
a third perturbation type says otherwise.

**A second finding, and it is a positive one.** On these 140 druggable target
genes the CRISPR gate is far stronger than panel-wide: decile-1 0.585 and pAUC
0.713, against 0.98 and 0.523 across the full panel. Drug targets are
selectively-essential oncogenes and kinases, and that is exactly where an
abundance prior should work. The gate's real domain may be selective
dependencies rather than all genes — worth testing directly.

**Coverage is reported per label per gene, never pooled** (`7A`,
`outputs/07_per_label_per_gene.parquet`):

| label | lines measured, median [p10, p90] | lines used after ∩ expression |
|---|---|---|
| Achilles | 767 [755, 767] | 754 [743, 756] |
| Score | 317 [317, 317] | 267 [249, 312] |
| GDSC2 | 932 [685, 960] | 702 [533, 920] |

### 08 — DepMap 24Q4 refresh: more lines, same artefact, different band

| | screened | expression | usable |
|---|---|---|---|
| 20Q2 bundle (Zenodo 14067047) | 767 | 1,479 | 754 |
| **24Q4 release** | **1,178** | **1,673** | **1,103** |

A 1.46× gain. Both CRISPR and expression taken from the same release, so no
vintage drift. **570 of 1,673 expression lines still have no screen** — the
audit-3 hazard shrinks but does not disappear and still has to be handled.

The denominator artefact reproduces at the larger *n*: decile-1 0.8631 →
**0.9788** on screened lines only (+0.1157, p = 1.05e-22 paired). It was never a
small-sample effect. Validity again does nothing (0.9788 all vs 0.9809 valid).

**The banded result moves, and this supersedes audit 6.** With 1,002 scored genes
on the corrected denominator:

| band | genes | decile-1 ratio | p | gate pAUC | p |
|---|---|---|---|---|---|
| **<2%** | 216 | **0.7098** | 0.00066 | 0.5970 | 0.0049 |
| **2–5%** | 158 | **0.8003** | 0.00139 | 0.5972 | 0.0033 |
| 5–10% | 89 | 0.9056 | 0.148 | 0.5506 | 0.184 |
| 10–25% | 120 | 0.9757 | 0.295 | 0.5146 | 0.554 |
| >25% | 313 | 1.0041 | 0.278 | 0.4931 | 0.071 |

The gate lives in the **rare-dependency bands (<5%)**, not the 5–10% band that
audit 6 identified on the smaller, Sanger-filtered gene set. That is the
selective-dependency regime, and it agrees with audit 7's finding that the gate
is strongest on selectively-essential drug-target genes. Quote the <5% bands.

### 09 — three perturbation types: the gate is CRISPR-only

|  | genetic? | complete? | acts on the gene product? |
|---|---|---|---|
| CRISPR knockout | yes | yes | yes |
| RNAi knockdown | yes | **no** (partial) | yes |
| Drug (GDSC2) | no | no | **no** (annotated target may be wrong) |

DEMETER2 RNAi is the discriminating case: if the gate reproduces there, the drug
null was causal distance; if not, the gate is specific to CRISPR.

**A scale problem had to be solved first.** `DEP_THRESHOLD = -0.5` is a Chronos
quantity and does not transfer — curated essentials sit at **−0.98** under CRISPR
24Q4 but **−0.35** under DEMETER2, because knockdown is partial. Applying one
absolute cut would have compared three different prevalences and called the
difference biology. So all three labels are binarised identically: per gene, the
most-dependent *q* fraction of that gene's **own measured lines**.

Decile-1 ratio at matched prevalence (< 1 = depletion = gate working):

| q | CRISPR | p | RNAi | p | drug |
|---|---|---|---|---|---|
| 0.02 | **0.9034** | 3.1e-10 | 0.9980 | 5.3e-05 | 0.9067 (ns) |
| 0.05 | **0.9034** | 7.3e-16 | **1.1976** | 1.3e-30 | 1.1537 (ns) |
| 0.10 | 0.9937 | 3.2e-09 | **1.0694** | 1.3e-44 | 1.1252 (ns) |
| 0.25 | 0.9721 | 1.5e-06 | **1.0334** | 1.5e-45 | 1.0018 (ns) |

**RNAi does not merely fail to reproduce the gate — it runs in the opposite
direction.** At q = 0.05 the bottom expression decile is *enriched* 1.20× for RNAi
dependency, with gate-region pAUC 0.477 (below chance, p = 4.8e-29). Paired
within gene, CRISPR and RNAi correlate at Spearman **+0.08 to +0.14** — the two
genetic perturbations essentially do not agree about which genes show the effect.

Coverage, per label per gene (never pooled):

| label | lines measured, median [p10, p90] | lines used |
|---|---|---|
| CRISPR 24Q4 | 1,178 [1,178, 1,178] | 1,103 |
| RNAi DEMETER2 | 538 [336, 699] | 509 [310, 657] |
| GDSC2 | 898 [548, 960] | 665 [426, 709] |

**Verdict: two of three perturbations say no, and one says the opposite.** The
gate is a property of CRISPR knockout, not a general abundance–dependency law.

The most likely mechanism is mechanical rather than biological. CRISPR cutting
toxicity scales with copy number, copy number correlates with expression, and
that is precisely the confound Chronos and CERES were built to correct — the
residual would look exactly like this. **That is the next test, and the data for
it is already on disk** (`data/DepMap_24Q4/OmicsCNGene.csv` was not downloaded,
but `data/DepMap_Chronos/CCLEGeneCopyNumber20Q2.hdf5` is present, and
`src/pipeline/build_cn_layer_complete.py` already builds a CN layer): recompute
the banded gate with copy number partialled out. If the CRISPR gate survives CN
adjustment it is biology; if it vanishes, it is cutting toxicity.

---

## What this means for the thesis

Three findings compound, and the order matters.

**1. The published gate is an artefact of its denominator (03, confirmed at 24Q4 in 08).**
16% decile-1 depletion, pAUC 0.61, p = 2.8e-76. A placebo with the dependency
labels randomised reproduces 98% of that depletion. On screened lines only it is
1-2%. This reproduces at 1,103 lines, so it was never small-sample.

**2. What survives is real but narrow (08).** On the corrected denominator at
24Q4, depletion of 20-29% in the rare-dependency bands (<5% prevalence),
p < 0.002, pAUC ~0.597. Nothing detectable above 10% prevalence. The robustness
work all holds — tie rules (01), gene validity (02), bin granularity and region
boundary (05) all move the answer by less than 0.02 — it is just robustness
around a smaller, narrower number.

**3. That surviving effect does not generalise beyond CRISPR (07, 09).** Drug
sensitivity is flat. RNAi runs *opposite*, significantly. Two genetic
perturbations of the same genes, at matched prevalence, correlate at Spearman
+0.08. The gate is a property of CRISPR knockout.

The defensible claim is therefore narrow and conditional:

> Within CRISPR knockout screens, low RNA abundance is associated with reduced
> dependency for genes that are selectively essential (<5% of lines), with a
> depletion of roughly 20-29% and gate-region pAUC ~0.60 on a correctly
> constructed denominator. The effect does not appear under RNAi knockdown or
> drug sensitivity, and is therefore not established as a general
> abundance-dependency relationship.

That is a real, defensible, publishable-as-a-negative finding. It is much smaller
than what was written down, and it is the version that survives a viva.

### What to do next, in order

1. **Copy-number adjustment (09's open question).** CRISPR cutting toxicity
   scales with copy number; copy number correlates with expression. Recompute the
   banded gate with CN partialled out. If the gate survives, it is biology; if it
   vanishes, it is cutting toxicity and the finding becomes a methodological one.
   This is now the single highest-value test.
2. **Re-run every downstream result that consumed the gate.** Anything derived
   from `test_run_abundance_gate.py` inherits the denominator, including the
   NPV/abstention thresholds feeding the confidence tiers.
3. **Fix the construction at source.** `test_run_abundance_gate.py` and anything
   sharing its `d.index.isin(pos_by_gene[g])` pattern must intersect with the
   screened line set rather than defaulting unmeasured lines to negative.
4. **Migrate the pipeline to 24Q4** if the gate is retained in any form — 754 to
   1,103 lines is a 1.46x power gain for free, now that the files are on disk.

### Data added by this audit

| path | what | source |
|---|---|---|
| `data/DepMap_24Q4/` | CRISPRGeneEffect, expression TPM, Model, controls | figshare 27993248 (DepMap 24Q4 Public) |
| `data/DEMETER2/` | D2_combined_gene_dep_scores, CL data | figshare 6025238 (DEMETER2 data v6) |

`depmap.org` itself is behind a Cloudflare Turnstile human-verification challenge
and cannot be fetched programmatically; figshare hosts the same releases without
one. Note that figshare's newest DepMap release is 24Q4 — 25Q+ releases appear to
be portal-only.

BioGRID ORCS was deliberately skipped: 106 extra lines under heterogeneous
protocols is not worth the harmonisation cost, and after 09 it would not change
the conclusion.
