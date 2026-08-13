# Track C — gene-axis EDA

Re-examination of the four Track C alteration tables against the Stage 0 gene
dimension, in two passes: what is on the gene axis, and what would break a model
built on it.

Reproduce with:

```bash
python src/pipeline/track_c_eda.py
```

Read-only. Writes `src/pipeline/outputs/track_c_eda_stats.json` and six figures
to `src/pipeline/outputs/figures/track_c/`. **Every number below comes from that
JSON** — none is transcribed by hand. Provenance and collection limitations are
in [TRACK_C_PROVENANCE.md](TRACK_C_PROVENANCE.md).

Inputs: `cleaned_track_data/` (the live pipeline copy, not the stale
`src/Track - C/outputs/` one) joined to `main.gene_enriched` in the warehouse.
Keys lowercased on read, following HARMONISATION_V2 §5.2.

---

## Part 1 — Decisions settled

The twelve design questions, resolved. Rationale is in the sections cited.

| # | Question | Decision | Basis |
|---|---|---|---|
| 1 | `False` vs `NA` sub-threshold | **Three states**: `True` / `False` (assayed, nothing ≥0.15) / `NA` (outside assayed universe) | §1.1 — coverage is 18,058 of 20,163 genes and 1,744 of 2,145 lines, so no-row conflates two things |
| 2 | Surface the 0.15 floor? | **In negative claims only** | §2.2 — a floor changes how a negative reads; irrelevant to a positive |
| 3 | Ship ploidy-adjusted VAF? | **Yes, as a secondary named column**, `vaf` stays primary | gene-level CN covers 0.3% of pairs; whole-line ploidy is a weak proxy |
| 4 | Impute missing ploidy? | **No — leave null** | structural 17% WGS block, not MAR |
| 5 | Drop `vep_rank_fallback`? | **Drop it; route on mechanism instead** | §2.1 — the fallback assigns near-constant severity to the most severe variants |
| 6 | REVEL vs AlphaMissense | **Keep both** + concordance flag; combine on percentiles, not raw values | r = 0.684, 19.6% discordance, asymmetric 2.3:1 |
| 7 | Low-confidence fusions | **Include with a reliability discount**, reusing the Layer 4 mechanism | §2.5 — confidence separates FFPM 5.6×, so the discount is calibrated |
| 8 | Out-of-frame fusions | **Enter scoring**, frame as a role-dependent magnitude modifier | §2.5 — frame is independent of confidence |
| 9 | Normalisation axis | **Gene-wise across lines**, with shrinkage below n = 20 | §2.3 — 33.3% of genes fall below that |
| 10 | Logit VAF user-facing? | **No** — raw VAF displayed, logit internal only | 342 variants sit at exactly 1.0 where logit is undefined |
| 11 | Collapse Layer 1 and 2? | **Never in storage**; combine only at query time | 21.5% of rows are "detected, magnitude unavailable" |
| 12 | Gate or additive? | **Gate** — magnitude is undefined, not zero, when nothing is detected | §1.2 — at 2.01% density, fabricated zeros would dominate every gene statistic |

Q11 and Q12 together make this a **hurdle model** (Cragg 1971): a binary presence
process plus a continuous magnitude process conditional on presence. That is the
standard construction for semicontinuous data and can be cited as such.

---

## Part 2 — Simple EDA: what is on the gene axis

![gene-axis coverage](../src/pipeline/outputs/figures/track_c/01_gene_axis_coverage.png)

### 1.1 Universe reconciliation — clean

| | genes |
|---|---:|
| Stage 0 gene dimension | 20,163 |
| HGNC protein-coding | 19,215 |
| mutations | 18,058 |
| fusions | 16,017 |
| **outside the dimension** | **0 / 0** |
| **not protein-coding** | **0 / 0** |
| never altered by either | 1,497 |

The universe filter Track C added holds perfectly — zero drift in either
direction. This is the one part of Track C that needs no rework.

Gene-axis overlap: 15,409 genes carry both mutation and fusion evidence, 2,649
mutation-only, 608 fusion-only.

**One anomaly:** a gene with an empty `hugo_symbol` ranks 6th most-mutated (703
lines) and 3rd most-fused (516 lines). It is inside the dimension and passes the
protein-coding filter but has no symbol to display, so it will surface in any
top-N ranking as a blank. Worth resolving before the ranking layer ships.

### 1.2 Density — the number that drives everything else

| | genes × lines | possible | observed | dense |
|---|---|---:|---:|---:|
| mutations | 18,058 × 1,744 | 31,493,152 | 632,919 | **2.01%** |
| fusions | 16,017 × 1,699 | 27,212,883 | 144,221 | **0.53%** |

Rows exist only where something was detected. This is the single most important
structural fact about Track C, and it is what forces the gate in Q12: writing a
zero into the 98% of empty cells would make every gene-wise mean and MAD a
statement about fabricated data.

### 1.3 Coverage per gene

| | median lines/gene | genes seen in < 20 lines |
|---|---:|---:|
| mutations | 26 | 35.1% |
| fusions | 5 | **90.5%** |

Fusions are far thinner than mutations. Any per-gene normalisation of fusion
evidence is estimating a distribution from a median of five observations —
shrinkage is not optional there, it is the only defensible option.

### 1.4 Alteration rate by curated role — the sanity check passes

| role | genes | median lines mutated | median lines fused | % ever mutated |
|---|---:|---:|---:|---:|
| both onco+TSG | 70 | 46.0 | 9.5 | 97.1 |
| TSG | 254 | 49.5 | 10.0 | 99.2 |
| oncogene | 257 | 36.0 | 8.0 | 98.4 |
| unknown | 19,582 | 23.0 | 3.0 | 89.3 |

Cancer Gene Census genes are altered roughly twice as often as background, and
TSGs more than oncogenes — which is what tumour-suppressor biology predicts,
since loss of function can be achieved by many more distinct variants than
gain of function can. The data behaves.

### 1.5 Top genes — and why the ranking is not yet trustworthy

Most mutated: **MUC4** (1,099), **TP53** (1,096), **TTN** (968), **MUC3A** (908),
**MUC16** (724), *(blank)* (703), **MT-ND5** (569), **MT-CYB** (497).

Only TP53 is a real driver. MUC4, TTN, MUC3A and MUC16 are the classic
long-and-repetitive artefact genes, and two mitochondrial genes follow. A ranking
built on raw counts would present four artefacts and two different-quantity genes
in its top eight. §2.4a and §2.6 quantify why.

---

## Part 3 — In-depth EDA: what would break a naive model

### 2.1 Missingness predicts outcome, strongly

Sweeping every nullable column against every outcome flag:

| null column | null rate | outcome | rate when null | rate when present | lift |
|---|---:|---|---:|---:|---:|
| `ampathogenicity` | 29.5% | `likelylof` | 0.552 | 0.003 | **216×** |
| `revelscore` | 22.9% | TSG high-impact | 0.042 | 0.0002 | **191×** |
| `revelscore` | 22.9% | oncogene high-impact | 0.018 | 0.0001 | **160×** |
| `revelscore` | 22.9% | `likelylof` | 0.700 | 0.006 | **125×** |
| `proteinchange` | 3.2% | `likelylof` | 0.985 | 0.138 | 7.1× |
| `revelscore` | 22.9% | `hotspot` | 0.005 | 0.005 | 1.1× |
| `ampathogenicity` | 29.5% | `hessdriver` | 0.002 | 0.004 | 0.5× |

Severity-missing rate by VEP impact class: **high 98.5%**, modifier 92.8%,
low 96.0%, moderate 6.4%.

AlphaMissense and REVEL are *missense* predictors; nonsense, frameshift and
splice variants have no missense score by construction. So the "severity hole" is
inverted — the unscored rows are the **most** consequential variants, enriched
190-fold for tumour-suppressor high-impact hits.

Note the two rows at the bottom that do *not* show lift: `hotspot` is flat (1.1×)
and `hessdriver` is actually *depleted* (0.5×). Both are missense-oriented
annotations themselves, so their behaviour is consistent with the same mechanism
rather than contradicting it.

**Consequence:** `mutation_scores.parquet`'s `severity_source = vep_rank_fallback`
assigns a near-constant severity to precisely this class. Route on mechanism
(`likelylof`, `vepimpact`, the role flags — all 100% populated) instead.

### 2.2 The detection floor

![detection floor](../src/pipeline/outputs/figures/track_c/02_detection_floor.png)

VAF minimum exactly **0.150**, maximum 1.000, median 0.444. 1.09% sit at the
floor; 0.07% at the ceiling, of which **342 are exactly 1.0** — where the logit
transform is undefined and will produce infinities unless squeezed.

The right panel is the control: VAF distributions are nearly identical for
severity-scored and severity-absent variants (median 0.447 vs 0.444), so the
severity gap in §2.1 is not a sequencing-depth or quality artefact.

### 2.3 Sparsity versus the normalisation axis

![MAD degeneracy](../src/pipeline/outputs/figures/track_c/03_mad_degeneracy.png)

**6,016 genes (33.3%)** have fewer than 20 variants to normalise against.

One correction to an earlier concern: **MAD collapses to zero for only 92 genes
(0.5%)**, not the widespread degeneracy I expected. Robust z is therefore viable
on the median gene, and switching to biweight midvariance is not necessary —
low n, not MAD failure, is the real constraint. Empirical-Bayes shrinkage toward
the pooled distribution handles it and degrades gracefully, where a hard n ≥ 20
cutoff would discard a third of the gene axis.

### 2.4 Three confounds

![confounds](../src/pipeline/outputs/figures/track_c/04_confounds.png)

**(a) Gene length.** Spearman correlation between log observed genomic span and
log lines mutated is **0.464** across 17,483 genes. Longer genes accumulate more
mutations for purely mechanical reasons. This is why MUC4, TTN and MUC16 top the
raw ranking. Any per-gene alteration rate needs a length offset — the standard
treatment is to model counts with `log(length)` as an offset term rather than
comparing raw counts.

**(b) Hypermutators.** The top **1%** of cell lines contribute **11.7%** of all
(gene, line) mutation pairs; the top 5% contribute **32.9%**. Median genes
mutated per line is 197; the maximum is **7,495** — 38× the median. These are the
MSI and POLE-deficient lines (6.5% of the panel exceed the MSIsensor-pro
threshold, per provenance §3). Without per-line normalisation, gene-level
frequencies are substantially a statement about which hypermutators are in the
panel.

**(c) Lineage.** Median burden ranges from 75 (eye) to 1,677 (endometrium) — a
**22.4× spread**. Endometrial lines are heavily MMR-deficient, so (b) and (c) are
partly the same effect. Since the panel's lineage composition is a historical
artefact rather than a sampling design, tissue must be a stratifying variable,
not something averaged over.

### 2.5 Fusion: confidence and frame are orthogonal

![fusion axes](../src/pipeline/outputs/figures/track_c/05_fusion_axes.png)

| confidence | n | median FFPM | in-frame rate |
|---|---:|---:|---:|
| low | 63,799 | 0.0495 | 0.137 |
| medium | 33,194 | 0.0940 | 0.180 |
| high | 47,228 | 0.2773 | 0.158 |

FFPM separates 5.6× across tiers; in-frame rate is flat. Confidence answers "is
this call real" and frame answers "is this call consequential" — two axes, two
layers. This is the measurement behind decisions 7 and 8.

**Filter discrepancy:** 45,217 rows (31.4%) fall below the FFPM ≥ 0.05 threshold
DepMap documents for the filtered fusion release, and 27 are exactly zero. Since
`best_ffpm` is a maximum over contributing events, every contributing event must
also have been below the floor. Resolve before the re-run — see provenance §2.

### 2.6 Mitochondrial variants are a different quantity

![mitochondrial VAF](../src/pipeline/outputs/figures/track_c/06_mitochondrial_vaf.png)

**2,243 variants** across **8 genes** in **1,210 lines** (0.32% of variants).

| | median VAF | VAF > 0.95 |
|---|---:|---:|
| nuclear | 0.444 | 4.7% |
| mitochondrial | **0.973** | **59.3%** |

Nuclear VAF is the fraction of 2–4 chromosome copies carrying an allele.
Mitochondrial VAF is heteroplasmy across hundreds of mtDNA copies per cell, and
most mtDNA variants in a clonal line have drifted to homoplasmy. They share a
column, a name and a [0,1] range, and they are not the same measurement.

Two of these genes rank 7th (MT-ND5) and 8th (MT-CYB) among the most-mutated in
the whole panel, so under a VAF-based magnitude they would rank near-maximal
essentially by construction. Mitochondrial variants need either a separate scale
or explicit exclusion — and ploidy scaling is meaningless for them regardless.

---

## What this changes

1. **`vep_rank_fallback` should go**, and severity should route on mechanism.
   This is the most consequential finding — the current design systematically
   down-weights truncating tumour-suppressor variants.
2. **Mitochondrial variants need separating** before any VAF-based magnitude.
3. **The fusion FFPM discrepancy needs resolving** — it changes what a fusion
   call means.
4. **Gene length, hypermutator and lineage offsets** belong in the model, not in
   a caveat paragraph.
5. **The stale `src/Track - C/outputs/` copy should be deleted** so the track
   folder stops disagreeing with the pipeline input, and `track-C.md`'s §1 table
   should be corrected to 18,058 / 16,017.

The universe filter, the key resolution and the collapse logic all hold up. The
rework is in the scoring layer, not the cleaning layer.
