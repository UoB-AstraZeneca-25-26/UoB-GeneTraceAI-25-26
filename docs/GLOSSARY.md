# GeneTraceAI — Glossary of Terms and Quantities

Every term and every measured quantity used anywhere in this project, defined
from first principles. Written to be **lifted directly into a dissertation**:
each entry stands alone, defines the thing before using it, and states what its
values actually mean rather than assuming the reader already knows.

Each entry follows the same shape:

> **What it is** — the definition, in one or two sentences, assuming no prior knowledge.
> **How it is measured** — the assay, algorithm or database it comes from.
> **Range and typical value in this project** — measured from this repo's own tables, not asserted.
> **Why it matters here** — what the pipeline does with it.

Companion documents: [MATH_REFERENCE.md](MATH_REFERENCE.md) (the estimators and
their properties) and [BIOLOGICAL_FIDELITY.md](BIOLOGICAL_FIDELITY.md) (what the
architecture claims biologically and where it falls short).

**All distributional figures below were measured from this repository's tables**
on the dates the pipeline was last run. Source table is named in each entry.

---

## A. Biological entities

### Cell line

**What it is.** An immortalised population of cells originally taken from a
single patient's tumour and propagated indefinitely in culture. Because every
cell in the line descends from that one sample, a cell line behaves as a single
reproducible genotype that any laboratory can obtain and experiment on. Cell
lines are the standard preclinical model system: before a drug is tested in an
animal or a person, it is tested against panels of them.

**How it is identified here.** Three identifier systems, which is why the
pipeline needs a resolver at all:

| System | Format | Issued by |
|---|---|---|
| DepMap model ID | `ACH-000123` | Broad Institute DepMap |
| Cellosaurus accession | `CVCL_0031` | Cellosaurus (ExPASy) |
| RRID | `RRID:CVCL_0031` | Research Resource Identification Initiative |

Plus free-text names (`A375`, `MCF7`), which are ambiguous — the same name can
denote different lines, and the same line appears under many spellings.
`model_id` (the DepMap `ACH-` form) is this project's canonical key.

**Scale in this project.** 1,485 cell lines carry a `core_score`
(`src/pipeline/outputs/core_score.parquet`).

**Why it matters.** The pipeline's entire output is a ranking *of cell lines*:
given a gene, which lines are the best experimental models for targeting it.

### Gene, and the ENSG identifier

**What it is.** A stretch of DNA encoding a product — usually a protein. An
**ENSG identifier** (`ENSG00000157764` = BRAF) is Ensembl's stable accession for
a gene; unlike gene *symbols* (`BRAF`), which get renamed by committee and are
reused, ENSG IDs are permanent. This is why the pipeline joins on ENSG and treats
symbols as display labels only.

**Scale in this project.** 19,213 protein-coding genes in
`reference/gene_lookup.parquet`; 19,176 of them receive a `core_score`.

### Oncogene, tumour suppressor gene (TSG), and gene role

**What they are.** Two opposite mechanisms by which a gene can drive cancer.

- An **oncogene** causes cancer through *gain of function*: it becomes
  overactive. Mechanisms are activating mutation, amplification (extra DNA
  copies), or fusion to a partner that drives its expression. The therapeutic
  logic is to inhibit it — you need the protein to be *present* to drug it.
- A **tumour suppressor gene (TSG)** causes cancer through *loss of function*: it
  normally restrains growth, and cancer arises when it is switched off.
  Mechanisms are damaging mutation and deletion. There is nothing to inhibit —
  the therapeutic logic is indirect (synthetic lethality, or restoring the
  pathway).
- A gene marked **both** acts as either depending on context and tissue.

**How it is assigned here.** From COSMIC's Cancer Gene Census, joined into
`reference/gene_lookup.parquet` as `gene_role`.

**Distribution measured here** (`reference/gene_lookup.parquet`, n = 19,213):

| `gene_role` | count | share |
|---|---|---|
| `unknown` | 18,633 | 97.0% |
| `oncogene` | 256 | 1.3% |
| `tsg` | 254 | 1.3% |
| `both` | 70 | 0.4% |

**Why it matters.** Role determines *direction*. For an oncogene, high abundance
is the signal of interest; for a TSG, the pipeline inverts the score on the
grounds that low abundance indicates loss. §3.4 of
[BIOLOGICAL_FIDELITY.md](BIOLOGICAL_FIDELITY.md) shows that inversion is
measurably backwards for the question the project is validated against. Role also
gates which copy-number events count as evidence (amplification for oncogenes,
deletion for TSGs).

### Driver mutation vs passenger mutation

**What they are.** A tumour genome typically carries thousands of mutations, but
only a handful caused the cancer. A **driver** confers a growth advantage and was
positively selected during tumour evolution. A **passenger** is incidental — it
happened to occur in a cell that was expanding for other reasons, and carries no
functional consequence. The ratio is heavily skewed toward passengers, which is
the single most important statistical fact in cancer genomics: treating every
observed mutation as meaningful is the classic error.

**How it is assigned here.** `any_driver | oncogene_hit | tsg_hit` in
`cleaned_track_data/mutations_collapsed.parquet`, derived from COSMIC Cancer Gene
Census membership plus variant consequence.

**Why it matters.** This Boolean is the only DNA-level quantity that reaches the
production ranking, as a tie-break ahead of `core_score`.

### Lineage

**What it is.** The tissue of origin of a cell line — blood, lung, skin, colon
and so on. Cell lines from the same lineage share large amounts of biology, so
they are *not* independent samples.

**Why it matters.** Drug response is strongly lineage-dependent. This creates the
**lineage confound**: a variable that predicts drug response only because it is
secretly a proxy for tissue type looks predictive without being biologically
informative about the gene in question. Several candidate layers in this project
were rejected on exactly this basis (see `DECISIONS.md`), and the statistical
test for it is in [MATH_REFERENCE.md](MATH_REFERENCE.md) §C.12.

---

## B. Data sources

### DepMap / CCLE

The Broad Institute's **Cancer Dependency Map** and **Cancer Cell Line
Encyclopedia**: the reference multi-omic characterisation of ~1,500 cancer cell
lines — expression, mutation, copy number, fusion, and CRISPR knockout screens,
all on a common `ACH-` identifier. Supplies most of this pipeline's input.

### GDSC — Genomics of Drug Sensitivity in Cancer

Large-scale drug screens: hundreds of compounds tested against hundreds of cell
lines, each yielding a dose–response curve. The summary measure is **IC50**, the
drug concentration required to inhibit growth by 50% — *lower IC50 means more
sensitive*. Because each drug has a known intended protein target, GDSC can be
reindexed from (drug, cell line) to (gene, cell line), which is what makes it
usable as ground truth here.

**Why it matters.** GDSC is this project's **only ground truth**. Every accept/
reject decision reduces to "does this change improve agreement with GDSC?" This
is also a limitation worth stating in the report: the pipeline is validated
against *drug response to compounds targeting the gene*, not against the broader
question of whether a cell line is a good model for studying the gene.

### Chronos (CRISPR gene-effect score)

**What it is.** A measure of how much a cell line *depends* on a gene. In a CRISPR
knockout screen, every gene is disabled in turn and the resulting change in cell
proliferation is measured. Chronos is DepMap's algorithm for converting those raw
screen readouts into a per-(gene, cell line) effect score, correcting for copy
number and screen-quality artefacts.

**Scale.** Calibrated so that **0 = knocking the gene out has no effect**, and
**−1 = the median effect of a pan-essential gene** (one every cell needs). More
negative means stronger dependency. Positive values mean knockout slightly
*helped* the cells grow.

**Why it matters.** Dependency is arguably what the project is really trying to
predict — a good model cell line is one that *needs* the gene. Chronos was tested
as a third ranking layer and rejected; it survives as a per-gene validation check.
Its near-zero correlation with abundance (ρ ≈ 0.005) is the project's central
finding.

### COSMIC

The **Catalogue Of Somatic Mutations In Cancer** (Wellcome Sanger Institute).
Supplies two things here: the **Cancer Gene Census** (the curated oncogene/TSG
role assignments above) and the **copy-number calls** used to build the CNA layer.

### Human Protein Atlas (HPA)

A resource mapping protein expression across human tissues and cell lines. Its
convention is *binary*: a gene is "expressed" if its level exceeds 1.0 on the
log₂-TPM scale. This project scores continuously instead, and the `new-arch`
diagnostic compares the two conventions.

### Ensembl VEP — Variant Effect Predictor

A tool that takes a DNA variant and predicts its consequence for the gene:
`synonymous_variant` (no amino-acid change), `missense_variant` (one amino acid
substituted), `stop_gained` (protein truncated), and so on. The consequences form
a rough severity ordering, which this pipeline compresses to an ordinal
`max_vep_rank` in 0–3.

### AlphaMissense and REVEL

Two computational **pathogenicity predictors**. Given a specific amino-acid
substitution, each returns a score in [0,1] interpretable as the predicted
probability that the substitution is damaging to protein function. AlphaMissense
is DeepMind's structure-informed predictor; REVEL is an ensemble of earlier
methods. This pipeline takes the maximum of the two available scores as
`max_pathogenicity`.

### Cellosaurus

A cell-line knowledge resource that assigns each line a stable `CVCL_` accession
and records synonyms, misidentifications and derivation history. This project
uses it as the bridge between free-text names and DepMap `ACH-` identifiers.

---

## C. Molecular measurement quantities

### TPM and log₂(TPM + 1) — RNA expression

**What it is.** **TPM** (Transcripts Per Million) measures how much messenger RNA
of a gene is present in a sample. Raw sequencing gives read *counts*, which are
not comparable between genes (longer genes collect more reads) or between samples
(deeper sequencing collects more reads overall). TPM corrects both: it divides by
transcript length, then rescales so that the TPM values across all genes in a
sample sum to one million. A TPM of 50 therefore means "50 out of every million
transcripts in this cell are from this gene."

**Why the log.** Expression spans several orders of magnitude, so the raw scale is
dominated by a handful of very high genes. DepMap distributes **log₂(TPM + 1)**;
the `+1` pseudocount is there so that genes with zero expression map to 0 rather
than to −∞.

**Range measured here** (`data/parquet/2_DepMap_...TPMLogp1Profile.parquet`,
1,495 cell lines). Per-gene ranges vary widely by gene:

| Example gene | min | median | max |
|---|---|---|---|
| DPM1 (housekeeping) | 3.66 | 6.51 | 9.18 |
| TNMD (rarely expressed) | 0.00 | 0.00 | 5.25 |
| FGR (lineage-restricted) | 0.00 | 0.04 | 8.01 |

A value of **0 means not detected**; **1.0 is the Human Protein Atlas threshold
for "expressed"**; values above ~8 are strongly expressed.

**Why it matters.** This is ranking layer 1, and for 90.4% of scored pairs it is
the *entire* score.

### Proteomic log-ratio — protein abundance

**What it is.** Mass-spectrometry measurement of how much of a protein is present.
Critically, it is **relative, not absolute**: the value is a log₂ ratio of the
protein's signal in this cell line against a common reference channel. A value of
0 means "this cell line has the panel-average amount of this protein"; +1 means
twice the average; −1 means half.

**Consequence for interpretation.** You **cannot** compare two different proteins'
values to say which is more abundant — only compare the same protein across cell
lines. This is precisely why the pipeline percentile-ranks *within gene, across
cell lines* rather than across genes.

**Range measured here** (`cleaned_track_data/proteomics.parquet`): 375 cell lines
× 10,129 proteins. Per-protein distributions are centred near 0 with a typical
range of about ±2 to ±3.75.

**Coverage — the key limitation.** Only 375 of 1,485 cell lines have proteomics
at all, which is why the two-layer stratum covers just **9.6%** of scored pairs.

### TOTAL_CN — total copy number

**What it is.** The number of physical copies of a stretch of DNA — and therefore
of any gene on it — present in the cell's genome. A normal human somatic cell is
**diploid**: two copies of every autosomal gene, one inherited from each parent,
so `TOTAL_CN = 2` is the normal baseline. Cancer genomes routinely depart from
this: a region may be **amplified** to 5, 10 or more copies (typically the
mechanism by which an oncogene becomes overactive) or **deleted** down to 1 or 0
copies (typically the mechanism by which a tumour suppressor is lost).

**How it is measured.** Not by counting molecules directly. Copy number is
*inferred* from sequencing data by two signals: **read depth** (a region present
in four copies yields roughly twice the reads of a region present in two) and
**B-allele frequency** (the ratio in which the two parental alleles appear, which
shifts when one is gained or lost). Algorithms such as ASCAT combine these into an
integer copy-number call per segment. COSMIC's `TOTAL_CN` field, used here, is
such a call.

**Values and what they mean.**

| `TOTAL_CN` | Interpretation |
|---|---|
| 0 | Homozygous deletion — the gene is entirely absent |
| 1 | Heterozygous loss — one copy remains |
| 2 | Normal diploid |
| 3–4 | Gain (or normal, in a genome-doubled line — see below) |
| ≥ 5 | Amplification |

**The critical caveat — the baseline is not always 2.** `TOTAL_CN` is an
*absolute* count. Whether a given count is abnormal depends on the cell's overall
**ploidy** (next entry). In a genome-doubled line whose baseline is 4 copies, a
`TOTAL_CN` of 3 is a *loss*, not a gain. This project's thresholds
(`> 2.5` = amplification, `< 1.5` = deletion) assume a diploid baseline, and
[BIOLOGICAL_FIDELITY.md](BIOLOGICAL_FIDELITY.md) §3.3 measures what that assumption
costs.

**Why it matters.** Copy number feeds `has_cna_alteration`, which upgrades
confidence tiers. It never enters the ranking.

### VAF — variant allele frequency

**What it is.** When a mutation is called from sequencing, some reads covering
that position carry the mutant base and some carry the reference base. **VAF is
the fraction carrying the mutant.** It is a direct readout of what proportion of
the DNA in the sample carries the mutation.

**How to read it.**

| VAF | Interpretation (in a pure, diploid sample) |
|---|---|
| ≈ 0.5 | Heterozygous — one of the two copies is mutated |
| ≈ 1.0 | Either homozygous (both copies mutated) **or** the wild-type copy has been lost |
| < 0.5 | Subclonal — only a fraction of the cells carry it |

**Why the ≈1.0 case matters.** A VAF near 1 at a tumour-suppressor locus is
evidence that *both* functional copies are gone — the "second hit" of Knudson's
two-hit hypothesis. This makes `max_vaf` the correct **per-locus** proxy for loss
of heterozygosity, in contrast to `lohfraction`, which is genome-wide.

**Range measured here** (`cleaned_track_data/mutations_collapsed.parquet`,
n = 632,919, 100% populated): 0.15 to 1.00, median **0.452** — i.e. the typical
called mutation is heterozygous, as expected. The floor at 0.15 is the variant
caller's detection threshold, not biology.

### FFPM — fusion fragments per million

**What it is.** A gene fusion is detected in RNA sequencing by finding read
fragments that span the junction between two different genes. **FFPM counts those
junction-spanning fragments, normalised per million total fragments** — so it is
an expression-normalised measure of how strongly the fusion transcript is being
made. A high FFPM means the fusion is abundantly transcribed, which is evidence
both that it is real (rather than a sequencing artefact) and that it is
functionally relevant.

**Range measured here** (`cleaned_track_data/fusions_gene_level.parquet`,
n = 144,221): 0 to 26.28, **median 0.081** — an extremely right-skewed
distribution in which most detected fusions are weakly transcribed.

### In-frame vs out-of-frame fusion

**What it is.** DNA is read in triplets (codons). When two genes fuse, the
question is whether the second gene's codons still line up correctly after the
junction. **In-frame** means they do, and the cell translates a single chimeric
protein combining domains from both genes — this is how *BCR-ABL1* produces a
constitutively active kinase. **Out-of-frame** means the reading frame shifts, so
translation produces nonsense downstream of the junction and usually terminates
early, yielding no functional product.

**Measured here.** 15.3% of gene-level fusion records have `any_in_frame = True`.

**Why it matters.** Only in-frame fusions can produce a functional fusion protein,
which is the rationale for the `INFRAME_BOOST` term in fusion scoring.

---

## D. Genome-level signature quantities

These six describe the *state of the whole genome* of a cell line, not any
individual gene. All are per-cell-line scalars from
`cleaned_track_data/signatures_model_level.parquet` (1,955 lines; 1,622 with
copy-number-derived fields).

**A structural note that matters for the report.** Because these are
gene-agnostic — one number per cell line, identical for every gene — they can
never distinguish between genes within a cell line. Adding one as a scoring
*layer* would shift every gene's score for that line by the same amount, which is
structurally identical to the metabolomics and miRNA global scalars this project
rejected as lineage confounds. Their correct role is as **modifiers on how other
evidence is read**, not as evidence in their own right.

### Ploidy

**What it is.** The average number of copies of each chromosome across the whole
genome — the cell's overall genome size. Normal human somatic cells are diploid,
ploidy = 2. Cancer cells frequently carry more.

**Measured here.** Range 1.673 – 4.995, **median 2.974**, mean 3.000. **40.2% of
lines have ploidy > 3.**

**Why it matters.** Ploidy is the correct baseline against which `TOTAL_CN` should
be judged. That the median is ~3 rather than 2 is what makes the pipeline's
absolute amplification threshold of 2.5 problematic — it sits *below* the typical
baseline.

### WGD — whole-genome doubling

**What it is.** A single catastrophic event in which the entire chromosome
complement duplicates at once, taking a diploid genome to tetraploid. It is one of
the commonest large-scale events in cancer evolution and typically precedes
further chaotic gains and losses.

**Measured here.** `wgd` is Boolean: **1,060 True, 562 False**, 333 missing —
i.e. **65% of lines with data have undergone whole-genome doubling.**

**Why it matters.** In a WGD line the neutral copy number is 4, not 2, so every
absolute copy-number threshold is mis-calibrated for the majority of this panel.

### CIN — chromosomal instability

**What it is.** The *rate* at which a cell line gains and loses whole chromosomes
or large segments during division — a measure of ongoing genomic churn, as
distinct from the accumulated damage that ploidy and aneuploidy record.

**Measured here.** Range 0 – 0.859, median 0.533.

### Aneuploidy score

**What it is.** A count of how many chromosome arms deviate from the cell's modal
ploidy. The standard score counts **39 arms** (22 autosomes × 2 arms, less the
5 acrocentric short arms, which carry no unique genes), so the theoretical maximum
is 39 and the observed maximum here is exactly that.

**Measured here.** Range 0 – 39, median 19 — i.e. the typical line in this panel
has abnormal copy number on roughly *half* of all chromosome arms.

### LOH fraction

**What it is.** **Loss of heterozygosity** means that at a given locus, one of the
two parental alleles has been lost, leaving only one version present.
`lohfraction` is the proportion of the genome in that state.

**Measured here.** Range 0 – 0.93, median 0.186 — the typical line has lost one
parental allele across about 19% of its genome.

**Why it matters, and its limitation.** LOH is mechanistically central to tumour
suppressor inactivation. But `lohfraction` is a **genome-wide** figure: it cannot
say whether the *specific* locus of interest is affected. For a per-locus LOH
signal, `max_vaf` is the correct quantity.

### MSI score — microsatellite instability

**What it is.** Microsatellites are short repeated DNA motifs (e.g. `CACACACA`)
that DNA polymerase copies unreliably. Normally a **mismatch repair** system fixes
the resulting slippage errors. When mismatch repair fails, these repeats
accumulate insertions and deletions rapidly, and the genome becomes
**hypermutated** — carrying orders of magnitude more mutations than a stable
genome. MSI status is clinically important: MSI-high tumours respond well to
immunotherapy.

**Measured here.** Range 0.30 – 93.22, **median 2.19**, mean 6.06 — extremely
right-skewed. Most lines are microsatellite-stable; a small hypermutator tail
drags the mean to nearly three times the median.

**Why it matters.** Hypermutated lines accumulate mutations in essentially every
gene by chance. Any score that counts mutations will therefore rank MSI-high
lines highly *for every gene*, which is a confound rather than a signal. Measured
in [BIOLOGICAL_FIDELITY.md](BIOLOGICAL_FIDELITY.md) §3.3b.

---

## E. Pipeline-derived quantities

### `core_score`

**What it is.** The pipeline's central output: a number in [0,1] for each
(gene, cell line) pair expressing **how abundant that gene's product is in that
cell line, relative to all other cell lines**. It is a *within-gene,
across-cell-line* quantity — 0.9 means "this line is in the top 10% for this
gene", and says nothing about how this gene compares to any other gene.

**How it is computed.** Expression and proteomics are each converted to
percentile ranks within gene, then averaged where both are available (see
[MATH_REFERENCE.md](MATH_REFERENCE.md) §2.4).

**Measured here** (`src/pipeline/outputs/core_score.parquet`): 28,403,110 rows
covering 1,485 cell lines × 19,176 genes. Mean 0.499, median 0.477, full [0,1]
range.

### `n_layers`

**What it is.** How many evidence layers contributed to a given `core_score`: 1
(expression only, or proteomics only) or 2 (both).

**Measured here.** **90.39% of scored pairs are 1-layer; 9.61% are 2-layer** —
a direct consequence of proteomics covering only 375 of 1,485 cell lines.

**Why it matters, and the subtlety.** The two strata are not on the same scale.
A 1-layer score is a single percentile and is uniformly distributed; a 2-layer
score is an average of two percentiles and is therefore concentrated toward 0.5.
Measured standard deviations bear this out exactly:

| stratum | rows | measured sd | predicted sd |
|---|---|---|---|
| 1-layer | 25,672,370 | **0.2751** | 0.2887 (uniform) |
| 2-layer | 2,730,740 | **0.2418** | 0.2466 (at ρ = 0.46) |

The measured ratio 0.879 implies a between-layer correlation of **ρ ≈ 0.54**,
which sits inside the literature range of 0.46–0.58 for mRNA–protein correlation
— an independent confirmation that the two-layer stratum behaves as the theory
predicts. Full derivation in [MATH_REFERENCE.md](MATH_REFERENCE.md) §2.6.

### `stratum_rank`

**What it is.** `core_score` re-ranked into a percentile **within its
(gene, n_layers) group**, precisely to correct the scale mismatch above. This is
the number surfaced to users, and the only one that is safe to compare between
cell lines with different data coverage.

### `p_mutation`, `p_fusion`

**What they are.** Continuous scores in (0,1) summarising how much evidence there
is that a cell line carries a functionally significant mutation (respectively,
fusion) in a gene. Built by mapping each input feature through an S-shaped
(Hill) curve and combining.

**Critically: these do not affect the ranking.** The production gate is Boolean
(`has_driver_alteration`) and never reads them. They are displayed to users as
context, tagged in-code as `"CONFIDENCE MODIFIER — does not affect
stratum_rank"`. Their parameters are explicitly uncalibrated.

### `variant_count` and `variant_burden`

**What they are.** `variant_count` is the number of distinct mutations found in a
gene in a cell line; `variant_burden = log₁p(variant_count)` compresses that
count onto a gentler scale.

**Measured here** (n = 632,919 gene × cell-line records): range 1–37, but
**91.3% of records have exactly one variant** and only 1.65% have three or more.
The distribution is far more concentrated than the range suggests, which limits
how much any burden-based term can contribute.

### `max_vep_rank`

**What it is.** The severity of the worst variant consequence found in that gene,
compressed to an ordinal 0–3 from VEP's consequence classes (roughly: synonymous
→ low-impact → missense → truncating).

**Measured here.** 0–3, mean 2.17 — most records sit at the missense end.

### `max_pathogenicity`

**What it is.** `max(AlphaMissense, REVEL)` — the higher of two computational
predictions that a missense substitution damages the protein, in [0,1].

**Measured here.** 0–1, median 0.322. **Populated for 508,054 of 632,919 records
(80.3%)** — the missing 20% are variants for which no missense predictor applies
(e.g. truncating or non-coding variants).

### `max_confidence` (fusions)

**What it is.** The fusion caller's own confidence that the detected fusion is
real rather than an artefact, as an ordinal `low` / `medium` / `high`.

**Measured here.** low 44.2%, medium 23.0%, **high 32.7%**. Since
`fusion_driver = (max_confidence == "high")`, roughly a third of fusion records
pass the driver gate.

### `regime_class` and `regime_source`

**What they are.** `regime_class` records *how a gene's biology is expected to
relate to abundance*:

- `abundance_tracking` — more protein means more drug target; abundance should
  predict response directly.
- `activation_driven` — what matters is whether the gene is *activated* (by
  mutation or fusion), not how much of it there is.
- `loss_of_function` — the gene matters through being lost.

`regime_source` records whether that classification was **measured** (empirically
derived from GDSC) or assumed. Only measured genes are eligible for the highest
confidence tier.

**Why it matters.** `regime_class` selects the routing strategy: abundance-tracking
genes rank on `core_score` alone; activation-driven genes rank on the driver flag
first, `core_score` second.

### `confidence` tier

**What it is.** A categorical label — `high` / `moderate` / `low` / `unknown` —
summarising how much the pipeline's own evidence supports a given prediction. It
is deliberately **not** a 0–1 number: a numeric confidence would have to be
calibrated against held-out outcomes to mean anything, and no such calibration set
exists that is independent of the GDSC data already used for validation.

### `chronos_check`

**What it is.** Per gene, whether `core_score` ranking agrees with Chronos
dependency: `validated` (agrees), `inverted` (runs backwards), or `none` (no
detectable relationship).

**Measured here** (`src/pipeline/outputs/chronos_validation.parquet`, 16,866
genes tested, median 823 cell lines per gene):

| `chronos_check` | count | share |
|---|---|---|
| `none` | 14,839 | 88.0% |
| **`inverted`** | **1,258** | **7.5%** |
| `validated` | 769 | 4.6% |

**Why this table is one of the project's most important results.** Among genes
where abundance and dependency are related at all, **more run backwards than
forwards, by 1.6 to 1**. See [BIOLOGICAL_FIDELITY.md](BIOLOGICAL_FIDELITY.md) §3.1.

---

## F. Evaluation quantities

### `sensitive` (GDSC ground truth)

**What it is.** A Boolean per (gene, cell line): whether GDSC found this cell line
responsive to a drug targeting this gene. Derived by thresholding the dose–response
summary. This is the target the pipeline is scored against.

### hit@20 (hit rate at k = 20)

**What it is.** Take a gene; rank all cell lines by the pipeline's prediction;
look at the top 20. **hit@20 is the fraction of the cell lines *known* to be
sensitive that appear in that top 20.** It answers the question a user actually
asks: "if I can only run 20 experiments, how many of the truly good models do I
catch?"

**The denominator matters.** It is the number of GDSC-sensitive lines that also
have a `core_score` — the pipeline can only rank what it has scored, so any other
denominator would flatter or penalise it for coverage rather than accuracy.

**What counts as a good value.** Not obvious without a baseline. Under random
ranking the expected hit rate is `k/N` — about **2.2%** for a typical N ≈ 900
scored lines. Observed values in this project are ~1.9–2.5% against a ~1.3%
baseline on the evaluated subset, i.e. a real but modest improvement over chance.
A raw hit@20 figure reported without its baseline is uninterpretable.

### Spearman's ρ (rank correlation)

**What it is.** A measure of whether two rankings agree, from −1 to +1. It replaces
each value by its rank and measures correlation between the ranks. **+1** means
the two orderings are identical; **0** means no relationship; **−1** means exactly
reversed. Unlike Pearson correlation it does not assume a straight-line
relationship, only a consistent direction — which is appropriate here, since the
pipeline makes ordering claims, not magnitude claims.

**How to read the magnitudes in this project.** ρ ≈ 0.005 (abundance vs
dependency) means *no relationship whatsoever*. ρ ≈ 0.46–0.58 (mRNA vs protein)
is a moderate relationship — the two agree substantially but are far from
interchangeable. This project's `RHO_THRESHOLD = 0.1` is the floor below which an
association is treated as too weak to act on.

### Bootstrap confidence interval

**What it is.** A way of putting error bars on a statistic without assuming a
distribution for it. The dataset is resampled with replacement many times (1,000–
2,000 here), the statistic recomputed each time, and the middle 95% of those
recomputed values reported as the interval. **If the interval excludes 0, the
effect is unlikely to be sampling noise.**

**Why it appears everywhere here.** Most of this project's comparisons are
"does this candidate change improve correlation with ground truth?" — a
difference between two correlations, which has no simple closed-form standard
error. The bootstrap gives one directly.

### FDR — false discovery rate

**What it is.** When thousands of hypothesis tests are run at once, a 5%
significance threshold produces 5% false positives *by construction*: at 16,866
genes, roughly 843 would clear p < 0.05 even if nothing were real. **FDR control**
(the Benjamini–Hochberg procedure) adjusts the threshold so that a stated fraction
of the *declared* discoveries are expected to be false.

**Applied to this project's chronos validation:** 4,362 of 16,866 genes clear
p < 0.05, against 843 expected under a global null — so there is genuine signal.
Storey's estimate puts the null fraction at **π̂₀ ≈ 0.62**, i.e. about 38% of
tested genes carry a real association. All 769 `validated` genes have
p ≤ 4.0 × 10⁻³ and survive Benjamini–Hochberg at FDR 5% (threshold
6.9 × 10⁻³). See [MATH_REFERENCE.md](MATH_REFERENCE.md) §C.6 for why the
effect-size floor, not the p-value floor, is what makes this hold.

---

## G. Identifier and join keys

| Key | Format | Grain | Notes |
|---|---|---|---|
| `model_id` | `ACH-000123` | cell line | Canonical key; uppercased everywhere |
| `ensg_id` | `ENSG00000157764` | gene | Canonical key; bare, uppercase, no version suffix |
| `hgnc_symbol` | `BRAF` | gene | Display only — symbols are renamed and reused |
| `CVCL_xxxx` | `CVCL_0132` | cell line | Cellosaurus accession, used in the resolver |
| `profile_id` | DepMap-internal | assay run | One cell line may have several; 16 have two RNA profiles |
| `depmap_id` | `ACH-000123` | cell line | Same space as `model_id`; naming differs by table |

**The version-suffix trap.** Ensembl IDs appear in the wild both bare
(`ENSG00000157764`) and versioned (`ENSG00000157764.13`). Joining one form to the
other silently produces zero matches, which is why Stage 1 flags any join falling
below a 90% match rate: a low rate almost always indicates a key-format bug rather
than genuinely missing data.
