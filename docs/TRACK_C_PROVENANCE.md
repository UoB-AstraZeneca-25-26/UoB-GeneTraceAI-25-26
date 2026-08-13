# Track C — data provenance and collection limitations

Where the four alteration datasets came from, why they were generated, how they
were measured, and what those choices constrain. Companion to
[TRACK_C_EDA.md](TRACK_C_EDA.md), which measures the consequences in our copy.

Every claim here is either cited to upstream documentation or measured on this
repository (marked **measured**). Where the two disagree, both are stated.

---

## 0. A provenance gap in the repository itself

**No DepMap release quarter is recorded anywhere in this repo.** The raw files
are named `6_DepMap_OmicsSomaticMutationsProfile.parquet`,
`5_OmicsFusionFilteredSupplementary.parquet`, `14_OmicsGlobalSignatures.parquet`
with no version marker, and a grep for a `NNqN` release tag across `docs/` and
`src/Track - C/` returns nothing.

This matters more than it looks. DepMap re-calls every quarter; gene models,
filters and even the cell line roster change between releases. Two analyses run
against "DepMap mutations" a year apart are not comparable, and right now we
cannot say which release we hold. Only `CCLE_miRNA_20181103` carries a date, and
that is in the filename by luck rather than design.

**Recommendation:** record the release quarter and download date for all four
sources before the re-run, and treat that as part of the output contract.

---

## 1. Mutations — `OmicsSomaticMutationsProfile`

### What and why

DepMap's somatic variant calls across the CCLE cell line panel. The purpose is
to let a dependency signal (CRISPR/RNAi knockout, drug response) be interpreted
against genotype — "this line dies when you knock out X, and it carries a
mutation in Y". Our use is the same shape: alteration evidence per (gene, line).

### How collected

- **Caller:** GATK **Mutect2**, run per sample, followed by `FilterMutectCalls`
  ([DepMap mutation pipeline](https://storage.googleapis.com/shared-portal-files/Tools/25Q2_Mutation_Pipeline_Documentation.pdf)).
- **Sequencing:** whole-genome and whole-exome, keyed by *sequencing profile*,
  which is why the raw table is profile-grain and needs the `profileid → modelid`
  bridge that Track C applies.
- **Annotation:** Funcotator plus OpenCravat, which is where `revelscore`,
  `ampathogenicity`, `hotspot`, the OncoKB/CGC driver flags and the 22 external
  database columns Track C dropped all originate.

### The central limitation: no matched normal

Mutect2 is designed to separate somatic from germline variants by comparing a
tumour against a **matched normal** from the same patient. Cell lines have no
matched normal — the patient is long gone and there is no germline reference for
that individual. DepMap substitutes a panel of normals and population allele
frequencies, plus a custom `fix_mutect2.py` post-processing step the pipeline
documentation explicitly describes as handling a germline-site discrepancy.

That substitution is population-level, not individual-level. Consequences:

- **Residual germline contamination.** Rare private germline variants that are
  not in the population resource cannot be filtered and will appear as somatic.
- **The bias is not uniform across ancestries.** Population allele-frequency
  resources are best-powered for European-ancestry variation, so lines from
  under-represented ancestries retain proportionally more unfiltered germline
  calls. Any per-line burden statistic inherits that bias.

This is not a defect Track C can fix downstream. It is a property of the
measurement and belongs in the report as a stated limitation.

### Other collection artefacts

- **A detection floor at VAF 0.15.** **Measured**: minimum `af` is exactly 0.150
  across all 704,713 variants. Sub-threshold variants were removed upstream, so
  "not present" in this table means "not detected at ≥15% VAF".
- **Mitochondrial calls are included and are a different quantity.**
  **Measured**: 2,243 variants on `chrm` across 8 genes and 1,210 lines. Nuclear
  VAF is the fraction of 2–4 chromosome copies; mitochondrial VAF is heteroplasmy
  across hundreds of mtDNA copies per cell. See EDA §6 — they share a column and
  mean different things.
- **Cell lines drift.** Passage number, clonal selection and culture conditions
  change the genome over time. A variant called from one lab's stock at passage
  40 may be absent from another's at passage 12. Nothing in the data records
  passage.

---

## 2. Fusions — `OmicsFusionFilteredSupplementary`

### How collected

- **Caller:** **STAR-Fusion v1.6.0** on RNA-seq, run against the GENCODE v29
  plug-n-play resource, with `--no_annotation_filter --min_FFPM 0` so that
  filtering happens downstream rather than inside the caller
  ([depmap_omics pipeline docs](https://github.com/broadinstitute/depmap_omics/blob/master/documentation/DepMap_processing_pipeline.md)).
- **Documented post-filters:** remove fusions involving mitochondrial, HLA or
  immunoglobulin genes; remove STAR-Fusion "red herring" annotated calls; remove
  `SpliceType=INCL_NON_REF_SPLICE AND LargeAnchorSupport=No AND FFPM<0.1`; and
  **remove FFPM < 0.05**.

### Limitations

- **Fusion detection is expression-dependent.** STAR-Fusion reads RNA, so a
  genuine genomic rearrangement in a gene that is not transcribed is invisible.
  Absence of a fusion call is not absence of a rearrangement — it is absence of a
  *transcribed* rearrangement. This directly limits what a `detected = False`
  can be claimed to mean.
- **FFPM is an abundance measure, not a functional one.** Fusion fragments per
  million scales with the expression of the fusion transcript, so it inherits
  every expression confound — library depth, transcript length, and the tissue's
  baseline expression of the partner genes.
- **The confidence tier is a detection-reliability signal, not a consequence
  signal.** **Measured** in EDA §5: median FFPM separates 5.6× across confidence
  tiers, while in-frame rate is flat (0.137 / 0.180 / 0.158). They are orthogonal.
- **Our copy does not match the documented filter.** **Measured**: 45,217 of
  144,221 gene-level rows (31.4%) have `best_ffpm` below the documented 0.05
  floor, and 27 are exactly 0. Since `best_ffpm` is the *maximum* over
  contributing events, a value below the floor means every contributing event was
  below it. Either the file is the unfiltered supplementary companion rather than
  the filtered release, or the release we hold predates that filter. **This needs
  resolving before the re-run** — it changes what a fusion call means.

---

## 3. Signatures — `OmicsGlobalSignatures`

### How collected

Genome-wide instability phenotypes derived from DepMap's WGS pipeline: MSI score,
LOH fraction, CIN, ploidy, WGD, aneuploidy. Per the
[Sanger DepMap documentation](https://depmap.sanger.ac.uk/documentation/cell-models/msi-ploidy-mutational-burden/),
MSI is estimated with **MSIsensor-pro** (score ≥ 7 called MSI) and ploidy comes
from copy-number segmentation (PureCN / PURPLE).

These are **not** SBS-style mutational signature decompositions, despite the
table name. Track C's doc already flags this; worth repeating because the name
misleads.

### Limitations

- **Derived, not measured.** Every column is the output of an inference
  algorithm over copy-number or microsatellite data. Ploidy in particular is a
  model fit with known degeneracy — a ploidy-2 and ploidy-4 solution can fit the
  same segmentation, and callers pick between them with heuristics. Treat ploidy
  as an estimate with unquantified error, which is exactly why EDA §Q3 keeps it
  out of the primary VAF column.
- **A structural 17% hole.** **Measured**: 333 of 1,955 models (17.0%) are null
  across all five WGS-derived columns as one block; `msiscore` alone is complete.
  These are lines without adequate WGS. The missingness is structural, not
  random — do not impute.
- **MSI is rare in this panel.** **Measured**: 128 of 1,955 lines (6.5%) exceed
  the MSIsensor-pro threshold of 7. Small, but those lines carry enormous
  mutation burdens and dominate gene-level counts (EDA §4b).
- **WGD is the majority state.** **Measured**: 1,060 of 1,622 scored lines (65%)
  are whole-genome doubled. This is why assuming diploidy anywhere in the
  pipeline is wrong — mean ploidy is 3.00.

---

## 4. miRNA — `CCLE_miRNA_20181103`

### How collected

**NanoString nCounter** miRNA expression assay, normalised with NanoString's
nSolver software. Published as part of the CCLE next-generation characterisation
([Ghandi et al., *Nature* 2019](https://www.nature.com/articles/s41586-019-1186-3)),
which reports miRNA across 1,072 lines run in **14 batches**, each including two
replicates of K-562 as a control.

The `nmir00001.1`-style identifiers in our copy are NanoString probe IDs, which
is the fingerprint of this platform — **measured**: 734 probes, values ranging
0.47 to 283,656 (median 22.0), the heavily-skewed profile characteristic of
nSolver-normalised counts.

### Limitations — the most constrained of the four

- **Probe-based, not sequencing.** The assay measures a **fixed, pre-selected
  panel** of 734 miRNAs. Nothing outside the codeset can be detected, and no
  novel miRNA discovery is possible. Absence from this table means "not on the
  panel" at least as often as "not expressed".
- **Batch structure exists upstream and is lost in our copy.** The publication
  describes 14 batches with K-562 controls precisely because batch effects were
  expected. Our cleaned table has **no batch column** — the wide source
  (`mirna_clean.parquet`, 734 × 956) carries only `name`, `description` and cell
  line columns. So the batch effect is present in the values and unmodellable
  from what we hold. Recovering the batch assignment from the original CCLE
  supplement should be on the re-run list.
- **Oldest of the four sources.** Dated 2018-11-03, against DepMap omics that
  are several years newer. The cell line roster and identifiers have moved on,
  which is why Track C's column-header resolution was the hardest of the four.
- **No gene axis, by construction.** miRNA→target mapping is a separate inference
  problem (miRTarBase, TargetScan) and is not in this file. Track C's Option B
  decision to keep miRNA as its own entity rather than exploding onto genes is
  the right call given the source.

---

## 5. Cross-cutting: what none of these can tell you

Three limitations apply to all four and should be stated once in the report
rather than four times:

1. **Cell lines are not tumours.** They are clonally selected, immortalised,
   cultured lines that have adapted to plastic and serum for decades. Alteration
   frequencies here do not estimate alteration frequencies in patients.
2. **The panel is not a random sample.** Cell lines exist for cancers that were
   easy to culture and commercially interesting. **Measured** in EDA §4c: median
   mutation burden varies 22× across tissues, and the lineage composition of the
   panel is itself a historical artefact rather than a sampling design.
3. **Absence is ambiguous in every modality.** Mutations have a VAF floor,
   fusions need transcription, signatures need WGS, miRNA needs codeset
   membership. In no case does "no row" cleanly mean "not altered" — which is the
   direct justification for the three-state `detected` field in the EDA decisions.

---

## Sources

- [DepMap 25Q2 Mutation Pipeline Documentation](https://storage.googleapis.com/shared-portal-files/Tools/25Q2_Mutation_Pipeline_Documentation.pdf)
- [broadinstitute/depmap_omics — processing pipeline](https://github.com/broadinstitute/depmap_omics/blob/master/documentation/DepMap_processing_pipeline.md)
- [GATK — Mutect2 and germline filtering](https://gatk.broadinstitute.org/hc/en-us/articles/360050722212-FAQ-for-Mutect2)
- [Sanger DepMap — MSI, ploidy and mutational burden](https://depmap.sanger.ac.uk/documentation/cell-models/msi-ploidy-mutational-burden/)
- [Sanger DepMap — fusions dataset](https://depmap.sanger.ac.uk/documentation/datasets/fusions/)
- [Ghandi et al. (2019), Next-generation characterization of the Cancer Cell Line Encyclopedia, *Nature*](https://www.nature.com/articles/s41586-019-1186-3)
- [NanoString nCounter miRNA Expression Assay](https://nanostring.com/wp-content/uploads/2022/09/MAN-C0009-08-nCounter-miRNA-Expression-Assay-User-Manual.pdf)
