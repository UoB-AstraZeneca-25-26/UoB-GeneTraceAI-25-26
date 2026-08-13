# External validation sources — what's available, what's used, what's missing

Reference doc, not a decisions doc. Citations verified against PubMed/journal
records 2026-08-11 (one correction below). Companions:
[ALTERATION_DEFENCE.md](ALTERATION_DEFENCE.md),
[gate_audit/README.md](../gate_audit/README.md) (the deepest existing use of
these datasets), [TRANSCRIPTOMICS_DECISIONS.md](TRANSCRIPTOMICS_DECISIONS.md).

---

## What "validation" means for this pipeline

The pipeline answers: *for gene G, which cell lines are most relevant to
study it?* So the validation question is: *do the cell lines ranked highly
for G actually behave differently with respect to G than the ones ranked
low?* That defines the requirement for external data — an independent,
experimentally measured readout of a gene–cell line relationship, collected
by a different lab, a different method, and not used to build the ranking.

---

## Status at a glance

| Dataset | Perturbation type | **On disk?** | **Wired into pipeline?** | Verified citation |
|---|---|---|---|---|
| DepMap Chronos (Broad) | CRISPR KO | ✅ `data/DepMap_Chronos/GeneFitnessEffect_Chronos_Achilles.hdf5` | ❌ confidence tier reads a mislabelled substitute instead (see below) | Dempster et al. 2021, *Genome Biology*, DOI 10.1186/s13059-021-02540-7, PMID 34930405 |
| Project Score (Sanger) | CRISPR KO | ✅ `data/GDSC/Project_score_combined_Sanger_v2_Broad_21Q2_*`, `GeneFitnessEffect_Chronos_Score.hdf5` | ✅ `gate_audit/06_screen_replication.py` | Behan et al. 2019, *Nature* 568:511–516, PMID 30971826 |
| GDSC2 | Drug sensitivity | ✅ `validation/prepared/gdsc_*.parquet` | ✅ `05_confidence_tiers.ipynb`, `stage4_eval_save.py`, `gate_audit/07,09` | Iorio et al. 2016, *Cell* 166:740–754, PMID 27397505 |
| DEMETER2 (Achilles RNAi) | RNAi | ✅ `data/DEMETER2/D2_Achilles_gene_dep_scores.csv` | ✅ `gate_audit/09_perturbation_types.py` | McFarland et al. 2018, *Nat Commun* 9:4610, DOI 10.1038/s41467-018-06916-5 |
| DRIVE (Novartis RNAi) | RNAi | ✅ `data/DEMETER2/D2_DRIVE_gene_dep_scores.csv` | ❓ **not confirmed used anywhere** — present but unchecked | McDonald et al. 2017, *Cell* 170:577–592, PMID 28753431 |
| Broad/Sanger cross-study harmonisation | — | referenced, not a raw dataset | used conceptually in `06`/`08` | Pacini et al. 2021, *Nat Commun* 12:1661, DOI 10.1038/s41467-021-21898-7 |
| CTRPv2 | Drug sensitivity | ❌ not downloaded | ❌ | Seashore-Ludlow et al. 2015, *Cancer Discovery* 5:1210–1223, PMID 26482930; Rees et al. 2016, *Nat Chem Biol* 12:109–116, PMID 26656090 |
| PRISM Repurposing | Drug sensitivity | ❌ not downloaded | ❌ | Corsello et al. 2020, *Nature Cancer* 1:235–248, DOI 10.1038/s43018-019-0018-6, PMID 32613204 |
| Perturb-seq (Replogle) | CRISPRi, single-cell | ❌ not downloaded | ❌ | Replogle et al. 2022, *Cell* 185:2559–2575, DOI 10.1016/j.cell.2022.05.013, PMID 35688146 |
| Gonçalves et al. 2024 dependency map | curated marker associations | ❌ not downloaded | ❌ | Gonçalves et al. 2024, *Cancer Cell* 42:301–316, DOI 10.1016/j.ccell.2023.12.016, PMID 38215750 |

**The one citation correction from the original list:** Chronos is Dempster
et al. 2021 in **Genome Biology**, not *Nature Methods* — the PMID and DOI
originally attached to it belonged to a different paper. Every other
citation below checked out against PubMed/journal records.

**Reading the table:** five of the ten sources are already downloaded. Two
of those five are confirmed wired into production analysis
(`gate_audit/`), one (GDSC2) is wired into three different places, one
(real Chronos) is on disk but **not** what the confidence tier actually
reads, and one (DRIVE) is on disk with no confirmed consumer at all. Only
four sources are genuine net-new acquisitions.

---

## Category 1 — Independent CRISPR essentiality screens

### 1a. DepMap Chronos / Achilles (Broad)

Genome-wide CRISPR-Cas9 knockout, ~1,000+ lines. Chronos is the scoring
model (not the screen itself) — it corrects for copy-number and
guide-efficiency confounds during inference.

**What it validates:** for gene G, do the top-ranked lines show lower
Chronos essentiality (more essential) than the bottom-ranked ones — a
direct functional dependency readout.

**Where this project actually stands, verified this session:**
`build_chronos_validation.py` — the script that sets `essentiality_check`
in every pair explanation, including the CD86 example already traced —
reads `validation/prepared/chronos_long.parquet`.
[`gate_audit/README.md`](../gate_audit/README.md) states directly: *"The
label is the real DepMap Chronos release
(`data/DepMap_Chronos/GeneFitnessEffect_Chronos_Achilles.hdf5`), not
`validation/prepared/chronos_long.parquet`, which is Project Score scaled
BF negated."* So every "tested against CRISPR essentiality" statement the
CLI currently prints for unknown-class genes was checked against Sanger
Project Score under a Chronos label, not against Chronos.

Open questions this repair would answer:
- Of the ~18,500 unknown-class genes, how many change `chronos_check`
  once pointed at the real file?
- Does the pre-declared ρ≥0.30 / weak 0.10 floor (chosen specifically to
  avoid the GAPDH false-positive documented in the script's own docstring)
  still hold on real Chronos, or does the housekeeping-gene inflation
  reappear at a different threshold?
- Chronos coverage is ~1,000–1,700 lines depending on release; what
  fraction of the panel does that leave untested (`chronos_check =
  "untested"`)?

### 1b. Project Score (Sanger)

Independent CRISPR screen, ~320–1,100 lines depending on release,
Bayesian-Factor scoring rather than Chronos's population-dynamics model.
Already on disk in two forms (`data/GDSC/Project_score_combined_*` and
`data/DepMap_Chronos/GeneFitnessEffect_Chronos_Score.hdf5`) and already the
comparator in `gate_audit/06_screen_replication.py`, which found:

| | Broad | Sanger | paired p |
|---|---:|---:|---:|
| decile-1 ratio | 0.9927 | 0.9960 | 0.338 |
| gate pAUC | 0.5077 | 0.5051 | 0.14 |

Agreement is high on core-essential genes and much weaker in the rare-
dependency band the gate actually depends on — per-gene Spearman between
screens is 0.32–0.37 in that regime (`gate_audit` §06). This is already
measured, not a gap.

---

## Category 2 — Independent drug sensitivity screens

### 2a. GDSC2 — already wired, already returned a null

Used three ways today: (i) `05_confidence_tiers.ipynb`'s curated-gene
intersection, (ii) `stage4_eval_save.py`'s held-out hit-rate eval (the
ALTERATION_DEFENCE numbers), (iii) `gate_audit/07,09` as an independent
perturbation type. Result on (iii), matched-prevalence against the CRISPR
gate: **flat, pAUC ≈ 0.50, p = 0.48–0.66**, and paired within-gene
correlation between CRISPR-gate strength and drug-gate strength of
**−0.09 to −0.14** — no relationship. This is a legitimate negative
control result, not a failure: CRISPR dependency and drug sensitivity are
different causal chains (a drug also has to reach the cell and hit the
annotated target), so their disagreement is informative rather than
alarming — provided it's reported as such rather than left silent.

### 2b. CTRPv2 — not on disk, would replicate 2a's negative result

Different lab (Broad Center for Science of Therapeutics vs Sanger/EMBL-EBI
for GDSC), different, largely non-oncology compound set. If GDSC2 is a
genuine negative control, CTRPv2 is the obvious replication: a second null
on a differently-sourced drug panel would substantially strengthen "drug
sensitivity is decoupled from this ranking" as a claim, rather than
leaving it resting on one dataset. A second *positive* result would be the
more interesting (and more work to explain) outcome.

Before downloading: check gene-target coverage overlap with the current
panel first — CTRPv2's compound set skews toward non-oncology chemical
probes, so the matched-target gene set may be small.

### 2c. PRISM Repurposing — not on disk, widest compound coverage

4,518 compounds (77% non-oncology) across 578 lines, pooled/barcoded assay
— noisier per-well than GDSC2/CTRPv2's individual dose-response curves.
Its value here is coverage, not precision: for genes with no GDSC or CTRP
match, PRISM may still have a compound against the same target or pathway.
Given the higher noise floor, any PRISM-derived signal should carry a
larger effect-size bar before being trusted than the CRISPR-based results
above — consistent with how this project has treated weak correlations
elsewhere (the 0.10→0.30 ρ floor revision in `build_chronos_validation.py`
is the direct precedent).

---

## Category 3 — RNAi (orthogonal perturbation type)

### 3a. DEMETER2 (Achilles-integrated) — already used, already reversed

`gate_audit/09_perturbation_types.py`: RNAi does not fail to reproduce the
CRISPR gate, it runs **opposite** — bottom-decile expression is *enriched*
1.20× for RNAi dependency at q=0.05, gate-region pAUC 0.477 (below chance,
p=4.8e-29). CRISPR/RNAi agreement on which genes show any effect at all:
Spearman +0.08 to +0.14. This is explainable (RNAi knockdown is partial,
off-target effects are substantial, and it's a genetically different
perturbation from complete knockout) and is already the strongest evidence
in the repo that the CRISPR gate is CRISPR-specific rather than a general
abundance–dependency law — already correctly reported as such in
`gate_audit/README.md`'s final verdict.

### 3b. DRIVE (Novartis RNAi) — on disk, consumer unconfirmed

`data/DEMETER2/D2_DRIVE_gene_dep_scores.csv` exists. DEMETER2's own
integration (McFarland 2018) already folds Achilles RNAi, DRIVE, and a
76-line breast screen into one harmonised score — so it's possible DRIVE's
signal is already present *inside* DEMETER2's numbers rather than absent.
Whether `gate_audit/09` used the DEMETER2-integrated score or the
Achilles-only RNAi file specifically determines whether DRIVE is a real
gap or already represented. **Not yet checked** — worth confirming before
treating this as new work.

---

## Category 4 — Perturb-seq (Replogle et al. 2022)

Genome-scale CRISPRi with single-cell transcriptomic readout — K562
genome-wide (>2.5M cells, 9,867 genes), K562 Essential, RPE1 Essential.
Only 41% of perturbations produce a measurable transcriptomic effect;
typical perturbations shift ~45 genes, essential-gene knockdowns shift
>500.

**What's different about this one:** every other source here validates
the *dependency* claim (does abundance track essentiality/sensitivity).
Perturb-seq is the only source that could validate the **expression
layer** itself — whether the transcriptomic state Perturb-seq measures
after perturbing a gene matches what this pipeline's E layer implies about
a cell line's molecular state for that gene. That's a genuinely different
question from everything else in this doc.

**The hard limit:** only two cell lines (K562, RPE1). This validates gene
biology in depth and cell-line breadth almost not at all — most of the
panel's ~1,700 lines have no Perturb-seq counterpart. Only useful for
genes where K562 or RPE1 themselves rank highly in this pipeline's output.

---

## Category 5 — Gonçalves et al. 2024 (curated dependency-marker map)

930 cell lines, multi-omic, and a curated table of 370 clinically-prioritised
dependency-marker associations across 27 cancer types. Not a raw
perturbation dataset — a second analytical pipeline's *conclusions*,
independently derived (different feature engineering, different marker
selection) from the same underlying DepMap CRISPR data this project also
draws on.

**What it validates:** agreement with their curated associations is a
genuine external replication only where the two pipelines' feature choices
are independent. Where both use expression as an input feature, agreement
partly reflects using the same data twice, not independent confirmation —
worth checking their published methodology before treating any overlap as
clean external validation rather than partial self-agreement.

---

## Priority, given what's actually on disk

1. **Fix `build_chronos_validation.py`'s data source.** Zero new downloads
   — the real Chronos file is already present. This is the single highest
   fully-free action in this whole list, and it directly changes what
   every `[UNINFORMATIVE]`/`unknown`-confidence CLI line currently claims.
2. **Confirm whether DRIVE is already represented inside the DEMETER2
   score `gate_audit/09` used**, or genuinely unconsumed. Determines
   whether Category 3b is a real gap or already closed.
3. **CTRPv2 before PRISM**, if pursuing new downloads — CTRPv2 is a
   cleaner replication of the already-established GDSC2 negative result;
   PRISM's coverage advantage matters most once a null is already
   established twice, not before.
4. **Perturb-seq and the Gonçalves map are lowest priority** — both
   answer narrower or partially-overlapping questions relative to their
   acquisition cost, and neither closes a gap the existing `gate_audit`
   work has left open.
