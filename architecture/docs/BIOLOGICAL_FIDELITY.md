# How Faithfully Does This Architecture Model the Biology?

A structural audit of what biological claim the GeneTraceAI pipeline actually
makes, where that claim is well founded, and where biology enters as a label
rather than as structure.

**Written to be usable as report text.** Every biological mechanism is explained
from first principles before it is relied on, so that a reader who does not
already know what loss of heterozygosity or whole-genome doubling is can follow
the argument. Terms are defined at first use; the full definitions, with measured
ranges from this repository's own tables, are in [GLOSSARY.md](GLOSSARY.md). The
statistical machinery behind each number is in
[MATH_REFERENCE.md](MATH_REFERENCE.md).

**All numbers here are measured from this repo's own outputs, not asserted.**
Where a claim was later tested and failed, the failure is recorded in place
rather than removed — including one case where an earlier version of this
document was wrong.

---

## 1. The claim the architecture makes

### 1.1 The problem the pipeline is trying to solve

A researcher wants to study or drug a particular gene. Before any experiment in
an animal or a person, the work is done in **cell lines** — immortalised
populations of cancer cells, each descended from one patient's tumour, that can
be grown indefinitely and shared between laboratories. There are roughly 1,500
well-characterised human cancer cell lines available. Choosing which ones to use
matters: a line that does not depend on the gene in question will produce a null
result no matter how good the experiment is.

The pipeline's job is to rank cell lines for a given gene, best model first.

### 1.2 The claim, in one sentence

> **A cell line is a good model for targeting gene X if gene X's product is
> abundant in it.**

Everything else the pipeline computes — driver flags, copy-number alterations,
confidence tiers — modifies *how much you should trust* that statement. Nothing
else changes the ranking.

### 1.3 Where each data type sits on the central dogma

The **central dogma** of molecular biology describes the flow of information in a
cell: DNA is transcribed into RNA, RNA is translated into protein, and protein
does the functional work. Cancer genomics data types map onto that chain, and it
is illuminating to ask where in the chain this architecture actually operates.

| Biological level | What it measures | Data in this repo | Role in the architecture |
|---|---|---|---|
| **DNA** — mutation, fusion, copy number | What is written in the genome | Track C, COSMIC CNA | **Confidence modifier only.** Never ranks, except as a Boolean tie-break. |
| **RNA** — transcript abundance | How much message is being made | expression (log₂ TPM) | **Ranking layer 1** — and for **90.4%** of scored pairs, the *entire* score |
| **Protein** — protein abundance | How much functional product exists | proteomics (375 of 1,485 lines) | **Ranking layer 2** — active on only **9.6%** of pairs |
| **Function** — gene dependency | Whether the cell actually *needs* the gene | Chronos CRISPR screens | **Validation only.** Tested as a ranking layer, rejected. |
| **Phenotype** — drug response | Whether the cell dies when the gene is drugged | GDSC | **Ground truth only.** |

**The structural observation: the pipeline ranks on the middle of the central
dogma.** It uses the beginning (DNA) as an annotation and the end (function,
phenotype) only to check itself. Abundance is the one thing that actually orders
cell lines.

This is not automatically wrong — abundance is a reasonable first-order proxy,
and it is by far the best-covered data type. But it is a specific and
consequential architectural commitment, and it should be stated as one rather
than left implicit.

---

## 2. Where the biology is modelled correctly

These are genuine strengths and should be defended in the report as deliberate
choices rather than accidents.

**1. Per-gene, across-cell-line normalisation asks the biologically meaningful
question.** Expression values for different genes are not comparable — a
housekeeping gene sits at log₂ TPM ≈ 6.5 in every cell, while a
lineage-restricted gene sits at 0 in most and 8 in a few. Asking "is this gene
highly expressed?" in absolute terms is therefore meaningless. Ranking each gene
separately across cell lines (`rank(axis=0)`) asks instead: *is this line high
**for this gene***. That is the question a researcher actually has.

**2. Two-layer RNA + protein encodes a real biological fact.** Messenger RNA is
an imperfect proxy for protein: transcripts are degraded at different rates,
translated at different efficiencies, and proteins have their own half-lives. The
measured correlation is moderate — around ρ = 0.46–0.58 in the literature, and
**ρ ≈ 0.54 as implied by this project's own two-layer variance
([MATH_REFERENCE.md](MATH_REFERENCE.md) §2.6)**. The pipeline combines the two
layers with an operator that refuses to treat agreement between them as
independent corroboration, which is the correct response to that correlation.

**3. ~~TSG score inversion gets directionality right.~~ WITHDRAWN — measured
false.** The architecture does distinguish tumour suppressors from oncogenes
rather than treating them identically, which is the right *instinct*. But the
direction it chose is backwards for the question it is validated against:
inverting makes ranking about three times worse (§3.4a). **Credit the
distinction, not the direction.**

**4. Gene-role-gated copy number is correct mechanism-matching.** Copy-number
change means something different depending on the gene. An oncogene drives cancer
by becoming overactive, so *amplification* — extra DNA copies producing more
protein — is evidence. A tumour suppressor drives cancer by being lost, so
*deletion* is evidence. An amplified tumour suppressor is not evidence of
anything. The pipeline encodes exactly this asymmetry.

**5. Driver-versus-passenger separation reflects the central concept in cancer
genomics.** A tumour genome carries thousands of mutations but only a handful
caused the cancer; the rest are **passengers**, incidental damage carried along
by a cell that was expanding for other reasons. Gating on
`any_driver | oncogene_hit | tsg_hit` rather than "any mutation" is the single
most important correction available, and the pipeline makes it.

**6. The in-frame fusion boost is mechanistically correct.** DNA is read in
three-base codons. When two genes fuse, the chimeric transcript only produces a
working protein if the second gene's codons still line up after the junction —
**in frame**. Out of frame, translation produces nonsense and terminates early.
Only in-frame fusions can yield a functional fusion protein such as BCR-ABL1, and
only those receive the boost. (Measured here: 15.3% of fusion records are
in-frame.)

That is a defensible amount of biology. **The architecture is not naive.**

---

## 3. Where biology is missing — ranked by how much it costs

### 3.1 Abundance is not dependency — and this project measured that itself

**The biology.** There are two distinct questions one can ask about a gene in a
cell line:

- *Is the gene's product present?* — an **abundance** question, answered by
  expression and proteomics.
- *Does the cell need the gene to survive and grow?* — a **dependency** question,
  answered by knocking the gene out and seeing what happens.

These come apart constantly. Housekeeping genes are abundantly expressed in every
cell and many are entirely dispensable. Conversely, a transcription factor
present at low copy number can be absolutely required — a cell can depend
critically on something it makes very little of. Abundance tells you a gene is
*available*; dependency tells you it is *load-bearing*.

**How dependency is measured.** In a genome-wide **CRISPR knockout screen**, each
gene in turn is disabled and the effect on cell proliferation is recorded.
DepMap's **Chronos** algorithm converts those readouts into a per-(gene, cell
line) score calibrated so that 0 means knockout had no effect and −1 means the
effect typical of a gene every cell needs. More negative means greater
dependency.

**The measurement that undercuts the central claim.**

```
ρ(expression, Chronos essentiality) ≈ 0.005
```

Abundance and functional dependency are **essentially uncorrelated in this
dataset**. The pipeline ranks on abundance while holding direct evidence that
abundance does not track dependency.

**A second, stronger measurement pointing the same way.** The per-gene validation
table (`chronos_validation.parquet`, 16,866 genes tested, median 823 cell lines
per gene) classifies each gene by whether abundance-based ranking agrees with
Chronos dependency:

| `chronos_check` | count | share | meaning |
|---|---|---|---|
| `none` | 14,839 | 88.0% | no detectable relationship |
| **`inverted`** | **1,258** | **7.5%** | abundance ranking runs **backwards** against dependency |
| `validated` | 769 | 4.6% | abundance ranking agrees with dependency |

**Among the genes where abundance and dependency are related at all, more run
backwards than forwards — by 1.6 to 1.** This is a stronger statement than
ρ ≈ 0.005 on its own, because ρ ≈ 0 is consistent with "no signal anywhere",
whereas this table shows there *is* signal (Storey's estimate puts the null
fraction at π̂₀ ≈ 0.62, so roughly 38% of genes carry a real association) and it
predominantly points the wrong way.

**The architecture's response was defensible.** Chronos was tried as a third
ranking layer, collapsed BCL2's hit@20 to zero, and was demoted to a per-gene
validation check. Given that outcome, not shipping it was correct.

**But the consequence stands.** This is an **expression-availability model, not a
dependency model**, and the hit rates reflect it: hit@20 ≈ 0.019–0.025 against a
random baseline of ~1.3% on the evaluated subset.

**This is the honest headline framing of the whole project. It is not a weakness
to hide; it is the finding.** A quantified demonstration that the abundance
heuristic underpinning much cell-line selection does not track functional
dependency — and more often than not points the wrong way where it points at all
— is a more valuable scientific output than the ranking itself. Say it before you
are asked.

### 3.2 No pathway or network structure — the "why these columns" gap

**The biology.** Genes do not act alone. They operate in **pathways**: chains of
proteins that pass a signal from a receptor at the cell surface down to the
nucleus, where it changes which genes are switched on. The canonical example for
this project is the **MAPK pathway**, in which BRAF is one link in a chain
running RAS → RAF → MEK → ERK → proliferation.

This matters because *the reason a cell line responds to a BRAF inhibitor is
usually not that it has a lot of BRAF protein.* It is that a mutation such as
BRAF-V600E has locked the pathway into the "on" position, so the entire chain
downstream is constitutively firing and the cell has become addicted to that
signal. The relevant state is the **pathway's**, not the gene's.

**What the pipeline does.** Every gene is scored independently. There is no
pathway membership, no epistasis (genes modifying each other's effects), no
synthetic lethality (pairs where losing both is fatal but losing either alone is
survivable), and no upstream/downstream relationship anywhere in the
architecture.

**Why this is the largest omission.** The pipeline can represent "BRAF is
present" but has no way to represent "the MAPK pathway is switched on." This is
precisely why the feature columns feel unmotivated when you try to justify them
in writing: **they describe the gene's state, never the pathway's state.** Any
single-gene abundance model has this ceiling, and the ~2% hit@20 is what the
ceiling looks like in practice.

**What would fix it.** Pathway-level aggregation — scoring the pathway and
attributing that score to its member genes. This is a new architecture rather
than a patch, which is why it sits last in the priority list despite being the
biggest gap.

### 3.3 Mutational process is ignored — measured consequences

**The biology.** Two genomes can carry the same mutation and mean entirely
different things by it, because of *how* each genome became damaged. A genome
with functioning repair machinery that carries a mutation in a cancer gene has
probably been selected for it. A **hypermutated** genome carrying tens of
thousands of mutations has probably acquired that one by accident. The
**mutational process** is therefore a modifier on how any individual variant
should be read.

The six signature quantities in this repo (`msiscore`, `lohfraction`, `cin`,
`ploidy`, `wgd`, `aneuploidy` — all defined in [GLOSSARY.md](GLOSSARY.md) §D)
describe exactly this. They cover **95.4% of scored lines** and are **currently
used nowhere in scoring.** Two consequences have been measured.

#### (a) Ploidy-naive copy-number thresholds — live in production

**The biology.** **Copy number** is how many physical copies of a gene the cell
carries. Normal human cells are **diploid** — two copies, one from each parent —
so `TOTAL_CN = 2` is the baseline, and departures from it are meaningful:
amplification to 5+ copies is how oncogenes become overactive, deletion to 0 or 1
is how tumour suppressors are lost.

**But the baseline is not always 2.** **Ploidy** is the cell's *average* copy
number across the whole genome. **Whole-genome doubling (WGD)** — a single
catastrophic event that duplicates the entire chromosome complement — takes a
diploid genome to tetraploid in one step, and is one of the commonest large-scale
events in cancer evolution. In a genome-doubled line the neutral copy number is 4,
and a `TOTAL_CN` of 3 is a *loss*, not a gain.

**What the pipeline does.** `build_cna_layer.py` uses absolute cutoffs —
`TOTAL_CN > 2.5` for amplification, `< 1.5` for deletion — assuming a diploid
baseline throughout.

**Measured reality of this panel:**

```
median ploidy of scored lines        : 3.06   (not 2.0)
lines with ploidy > 3                : 53.0%
whole-genome doubled (full table)    : 1,060 of 1,622 with data = 65%
CNA-alteration rate, WGD vs non-WGD  : 0.0277 vs 0.0212  (1.31x, p = 1.4e-05)
ρ(mean TOTAL_CN, ploidy)             : +0.109, p = 0.002
```

*(The full signatures table, over all 1,622 lines with copy-number-derived
fields, gives median ploidy 2.97 and 40.2% above 3; the higher figures above are
for the subset that actually receives a score. Both are reported because the
difference itself is informative — scored lines skew toward higher ploidy.)*

**The amplification threshold sits below the median baseline copy number.** For
much of the panel, "amplified" partly means "this line has a large genome."
Because `cna_flags` feeds confidence-tier upgrades, this propagates into
user-facing labels. Standard practice in the field — GISTIC, ABSOLUTE, PureCN —
normalises to ploidy (`TOTAL_CN / ploidy`) precisely to avoid this.

**TESTED — `test_run_ploidy_normalised_cna.py`. The result splits in two, and an
earlier version of this document claiming it was "the single most fixable defect"
does not survive it.** Applying ratio cutoffs (`> 1.25` amplification, `< 0.75`
deletion, chosen so that diploid lines keep today's calls exactly):

```
ρ(alteration rate, ploidy) : +0.090 (p=0.012)  ->  +0.045 (p=0.212, n.s.)   ARTIFACT REMOVED
WGD vs non-WGD ratio       : 1.43x (p=3.0e-08) ->  1.40x (p=4.6e-07)        GAP PERSISTS
```

**Interpretation.** The arithmetic artefact was real and is removed — after
normalisation, alteration rate no longer tracks ploidy. But the WGD excess
survives almost untouched, which means it is largely **genuine biology**:
genome-doubled lines really do accumulate more focal copy-number events, because
having spare copies relaxes the selective constraint against losing one. A
correct threshold should *not* erase that, and does not. This is a useful
demonstration that removing a confound and removing a signal are different
operations, and that a well-chosen correction distinguishes them.

**Blast radius.** 120 of 117,540 pairs change (0.1%) — but **93 of those
currently carry "high" confidence and would be downgraded**, about 1% of all
high-confidence predictions, retracting claims that rested on ploidy-inflated
amplification calls. Ranking is mathematically unchanged, because copy number
never enters the ranking gate (asserted in the test).

**So the fix is still worth making, but for precision, not performance.** It
corrects which predictions you are entitled to call "high confidence", and
improves no metric.

#### (b) Hypermutation inflates mutation burden — real but small

**The biology.** **Microsatellites** are short repeated DNA motifs that the
copying machinery handles unreliably. A **mismatch repair** system normally fixes
the resulting errors; when that system fails, the genome becomes
**microsatellite-unstable (MSI-high)** and accumulates mutations at an
enormously elevated rate. The measured distribution here is extremely skewed —
`msiscore` has median 2.19 but a maximum of 93.2, so a small tail of hypermutator
lines sits far from the rest.

**The confound.** A hypermutated line acquires mutations in essentially every
gene by chance. Any score that counts mutations will therefore rank MSI-high
lines highly **for every gene**, which is a property of the cell line rather than
evidence about the gene.

**Measured:**
```
ρ(mean p_mutation,  msiscore) = +0.121, p = 4.4e-07
ρ(n_mutated_genes,  msiscore) = +0.172, p = 4.9e-13
```

**Status.** The `m_mut` signature discount was built to correct this and measured
null on held-out hit@20 (confidence interval included 0). But note carefully
*why* it was null: **`p_mutation` does not reach the ranking at all** — the
production gate is Boolean. The confound is real; the correction was applied to a
component that does not affect the metric being measured. That is a testing-design
lesson worth recording in the report, not evidence that the confound is harmless.

### 3.4 No two-hit logic for tumour suppressors — textbook biology, data available

**The biology.** **Knudson's two-hit hypothesis** is foundational to cancer
genetics. Because we carry two copies of every gene, disabling one copy of a
tumour suppressor usually changes nothing — the remaining copy still produces
enough functional protein. Loss of function generally requires **both** copies to
be inactivated. This is why inherited cancer syndromes work the way they do: a
person born with one defective copy needs only one further somatic event, while
someone with two good copies needs two.

**The practical consequence.** A tumour suppressor with a damaging mutation in
one copy *and* loss of the other copy is a genuinely inactivated gene. Either
alone is much weaker evidence. Distinguishing these two situations requires
per-locus information about both alleles.

**What the pipeline has.** Three of the necessary ingredients:

- `tsg_hit` — is this gene a tumour suppressor, and is it mutated
- `max_vaf` — **variant allele frequency**, the fraction of sequencing reads
  carrying the mutant base. In a pure diploid sample, ≈0.5 means one copy is
  mutated; **≈1.0 means either both copies are mutated or the wild-type copy has
  been lost**. This makes it the correct *per-locus* proxy for the second hit.
  (100% populated on TSG rows; median 0.452 overall.)
- `lohfraction` — the fraction of the genome showing **loss of heterozygosity**,
  i.e. where one parental allele is gone (83% coverage)

**What it does with them: nothing.** For tumour suppressors it simply inverts
`core_score`.

**TESTED — `test_run_tsg_two_hit.py`. The two-hit rule is null. But testing it
uncovered that the existing TSG inversion is itself backwards, which matters far
more.**

#### (a) The TSG inversion is measurably wrong for the validated question

On the 6 genes actually inverted in production (`gene_role = "tsg"` →
`loss_of_function`):

```
mean hit@20, HIGH abundance    : 0.022050
mean hit@20, LOW  (production) : 0.006933
delta +0.015117    ratio 3.18x    95% CI [+0.004117, +0.028008]   EXCLUDES 0
high beats low in 5/6 genes, 1 tie, 0 against
```

Unrestricted to all 14 TSG/`both` genes with ground truth: same direction, 3.01×,
CI [+0.0099, +0.0238], 13/14 genes.

**This is the only significant effect any diagnostic in this project has
produced, and it contradicts a live production rule.**

**The mechanism, which is the interesting part.** GDSC measures response to a
drug **targeting** the gene — and drugs need their target to be present. ATM,
ATR and CHEK2 inhibitors need the kinase to exist in order to inhibit it. TP53's
GDSC compounds are MDM2 inhibitors (nutlins), which work by releasing *functional
wild-type p53* from suppression — they require the tumour suppressor to be
intact. Inverting the score ranks the cell lines that have **lost** the gene,
which are exactly the lines that cannot respond.

The rule conflates two different questions:

| Question | Is inversion right? |
|---|---|
| Which lines have lost this suppressor's function? | **Yes** |
| Which lines respond to a drug targeting this gene? | **No** — and this is what GDSC measures, and what the pipeline is validated and used for |

**Caveats, stated plainly.** n = 6 (14 unrestricted). **Not held out** — the
curated 141-gene table contains zero `loss_of_function` genes, so the seed-42
split has no tumour suppressors in its test set. And any tumour suppressor with
GDSC ground truth is nearly by definition a drug target, which is plausibly why
target-presence wins. This is **exploratory, not confirmatory**. Note also that
only 6 of 14 are inverted at all — `gene_role = "both"` maps to
`activation_driven` and is left un-inverted.

#### (b) The two-hit rule itself: null

Gating on `tsg_hit AND (max_vaf ≥ 0.75 OR multi_hit_high_impact OR min_cn < 1.5)`
scored 0.012194 against the production gate's 0.015631 — paired delta
**−0.003437, CI [−0.010597, +0.004487]**, which includes 0.

**Why it failed is more informative than that it failed.** The second hit came
almost entirely from VAF (4,119 rows) and multi-hit (2,919); the deletion route
contributed only **56** rows, because just **103** tumour-suppressor
(cell line, gene) pairs carry both mutation and copy-number data. The per-locus
data needed for a proper two-hit call is far thinner than the headline coverage
figures suggest. This is a data-availability result, not a refutation of
Knudson's hypothesis.

**`lohfraction` was deliberately excluded**, and the reason generalises: it is a
genome-wide per-cell-line scalar, not a per-locus measurement. Using it would
repeat exactly the gene-agnostic-scalar error that got metabolomics and miRNA
rejected as lineage confounds — a number that is identical for every gene in a
cell line cannot be evidence about a particular gene. `max_vaf` is the correct
per-locus LOH proxy.

### 3.5 Lineage is not conditioned on by default

**The biology.** Cell lines derive from different tissues — blood, lung, skin,
colon — and tissue of origin determines a great deal of a cell's biology. Drug
response is correspondingly lineage-dependent: BCL2 inhibitors work in blood and
lymphoid lineages far more than elsewhere, because those lineages depend on BCL2
for survival in a way that solid tumours generally do not.

**The consequence for ranking.** Ranking across the whole panel means a gene whose
relevance is lineage-specific will surface its lineage rather than the best models
within it — and a score that appears predictive may be predicting tissue type
rather than anything about the gene.

**Status: in progress, not absent.** `DECISIONS.md` records that within-stratum
scoring has cleared its feasibility gates — **2,421 usable gene × stratum cells,
140 genes with ≥ 3 strata** — so the data supports it. Worth saying so explicitly
in the report rather than listing this as a plain gap.

### 3.6 Fusion partner identity is discarded

**The biology.** A **gene fusion** joins two previously separate genes into one
transcript. Whether that is oncogenic depends almost entirely on *which two
genes*. BCR-ABL1 — the Philadelphia chromosome, the target of imatinib — is
oncogenic because the BCR portion drives constitutive activation of the ABL1
kinase. ABL1 fused to something irrelevant is not the same event and carries no
therapeutic implication.

**What the pipeline does.** `fusions_gene_level.parquet` carries
`top_partner_ensg_id`, but fusion scoring uses only confidence, FFPM, recurrence
and in-frame status — never the partner's identity. The partner is often the whole
story, and it is discarded.

### 3.7 Genetic redundancy — paralog buffering is unmodelled and untested

**The biology.** Gene families arise by duplication, and duplicates frequently
retain overlapping function. Knock out one member and its paralog covers the loss,
so the cell shows little or no fitness defect — the gene is *functionally*
important but *experimentally* silent in a single-knockout screen. This is not a
rare edge case: synthetic-lethal paralog pairs are a recognised class of drug
target precisely because single knockouts miss them.

**Why it matters here specifically.** Paralog buffering attacks the pipeline from
both ends at once:

- It **depresses the ground truth.** Chronos and Project Score are single-knockout
  screens, so a buffered gene reads as non-dependent regardless of its abundance.
  Every validation figure in `MATH_REFERENCE.md` §C.15–C.16 is measured against a
  label that systematically under-calls this class of gene.
- It **breaks the architecture's central inference in the other direction.** The
  pipeline's core claim is that abundance is a necessary condition for dependency
  (§3.1). For a buffered gene the reverse holds: the gene can be abundantly
  expressed and still non-essential, because the paralog makes it dispensable. Such
  genes are systematic false positives for the abundance prior — and unlike the
  "abundance is not sufficient" gap of §3.1, this one has a known mechanism and a
  known gene list.

Multi-targeting CRISPR guides in near-identical families push the opposite way,
inflating apparent dependency — which is why `test_run_gate_mechanism.py` had to
drop KRTAP\* and VN1R\* from its negative-control set.

**What the pipeline does.** Nothing. There is no paralog annotation anywhere:
`reference/gene_lookup.parquet` has no such column, and no scoring or gating step
references gene family membership.

**What was attempted, and why it does not count as addressed.** HGNC `gene_group`
co-membership was tried as a proxy and **failed its own validity check** — see
`MATH_REFERENCE.md` §C.17. Paralog buffering predicts lower dependency prevalence
among paralogous genes; measured prevalence was 0.1072 with paralogs against 0.0629
without, `p = 0.7766` in the wrong direction. HGNC groups skew toward
well-characterised families — ribosomal proteins, proteasome subunits, spliceosome
components — which are pan-essential, so group membership proxies "belongs to a
named family" and tracks essentiality rather than redundancy.

The §C.16 finding does survive restriction to HGNC-singleton genes (0.790–0.799,
6/6 rules), but that is a robustness check, not a control for paralogy.

**What would close it.** Sequence-level paralog data — Ensembl Compara paralog
pairs with % identity, DepMap's paralog annotation, or the published Dede/Ito
paralog sets — then either stratify the validation by paralog status or exclude
buffered genes from the abundance claim. The biological payoff is high and the
effort is a single annotation join; the obstacle is purely that the data is not
local.

---

## 4. Conclusion

**The architecture is a well-engineered abundance model with mechanism-aware
annotation — but it is not yet a mechanism model.**

Biology enters this pipeline in three legitimate ways:

- as **direction** — inverting the score for tumour suppressors (though §3.4a
  shows the chosen direction is wrong for the validated question);
- as **gating** — driver versus passenger, role-matched copy number;
- as **labels** — confidence tiers.

It does not enter as **structure**. There is no representation of pathways, of
allele state, of mutational process, or of the difference between a gene being
*present* and a gene being *required*.

That is a coherent and defensible design for a first system. It is auditable,
every step is traceable, and the rejected experiments are documented rather than
buried — which is itself a contribution, since most of the value in this project
lies in the negative results. But it sets a ceiling on performance, and the ~2%
hit@20 against a ~1.3% baseline is what that ceiling looks like in practice.

**The project's most valuable scientific output is arguably not the ranking but
the measurement that abundance does not track dependency**: ρ ≈ 0.005 overall,
and — among the genes where any relationship is detectable at all — **1,258
running backwards against 769 running forwards**. That is a direct, quantified
demonstration that the abundance heuristic underpinning much of this style of
cell-line selection is not merely weak but frequently inverted. Framing the
dissertation around that finding is considerably stronger than framing it around
a 10.5% driver lift.

### Ranked next steps, by biological payoff per unit of effort

| # | Change | Status | Why | Effort |
|---|---|---|---|---|
| 1 | **Reverse or remove the TSG inversion** | **TESTED** (§3.4a) | 3.18× effect, CI excludes 0 — the only significant diagnostic result in the project, contradicting a live production rule. Exploratory (n=6, not held out), so investigate before flipping | Low |
| 2 | Ploidy-normalise CNA (`TOTAL_CN / ploidy`) | **TESTED** (§3.3a) | Removes a real arithmetic artefact; the WGD gap that survives is genuine biology. Buys precision — 93 high-confidence retractions — not performance | Low |
| 3 | TSG two-hit rule (mutation + LOH + VAF) | **TESTED** (§3.4b) | Null (CI includes 0), and the per-locus data is too thin to support it — only 103 TSG pairs have both mutation and CNA | Medium |
| 4 | Make within-lineage percentile the default | **OPEN** | Feasibility gates already passed (`DECISIONS.md`); drug response is lineage-driven | Medium |
| 5 | Use fusion partner identity | **OPEN** | The partner determines whether a fusion is oncogenic at all | Medium |
| 6 | Pathway-level aggregation | **OPEN** | Addresses the structural ceiling of §3.2 — but is a new architecture, not a patch | High |

Every item marked OPEN is genuinely untested — no diagnostic has been run against
any of them.

### One thing not to do

**Do not add signature scalars (CIN, aneuploidy, ploidy) as `core_score`
layers.** They are gene-agnostic per-cell-line values: one number per cell line,
identical for every gene. Adding one as a layer shifts every gene's score for that
line by the same amount and cannot distinguish between genes at all — which makes
them structurally identical to the metabolomics and miRNA global scalars this
project already tested and rejected as lineage confounds.

Their correct role is exactly what this document argues for throughout:
**modifiers on how other evidence is interpreted** — normalising copy number to
ploidy, discounting burden for hypermutation, qualifying a tumour-suppressor hit
with allele state — not independent evidence of their own.
