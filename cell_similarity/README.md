# `cell_similarity/` — the cell-line similarity graph

Self-contained folder, read-only with respect to the pipeline, writes only into
`cell_similarity/outputs/`. Three weeks of work as four scripts, plus two
robustness scripts added after review.

It reuses `gate_audit/common.py` for the DepMap 24Q4 loaders and their scale
checks — imported by explicit path under a distinct module name, because both
folders have a `common.py`. Reuse rather than a copy, so the scale guards cannot
drift.

---

## Run it

```bash
python cell_similarity/run_all.py
```

Then query it:

```bash
python cell_similarity/04_query.py --gene TP53 --line ACH-000552
```

| File | What it is |
|---|---|
| `common.py` | Loaders, the similarity metric, neighbour machinery. One definition of "similar". |
| `01_build_axes.py` | **Week 1.** Builds the three axes separately. |
| `02_validate_dependency.py` | **Week 2.** Do similar lines behave similarly? |
| `03_validate_tissue.py` | **Week 2.** Do same-tissue lines come out close? |
| `04_query.py` | **Week 3.** Library + CLI for the query path. |
| `05_robustness.py` | Cluster bootstrap over lines, RNAi arm, common essentials stripped. |
| `06_mirna_stability.py` | miRNA metric instability: diagnose, fix, or demote. |
| `07_reliability.py` | RNAi reliability, disattenuation, power control. Power or confound? |
| `08_drug_arm.py` | Third validation arm: GDSC2 drug response. |

---

## Week 1 — the graph

RNA first, then metabolomics, then miRNA. Per axis: drop features missing in
>20% of lines, median-impute, keep the top 2,000 most variable, z-score, Pearson
between lines (= cosine on centred data).

| axis | source | lines | features |
|---|---|---|---|
| RNA | DepMap 24Q4 TPM log2+1 | **1,673** | 18,623 → 2,000 |
| metabolomics | CCLE metabolomics | **928** | 225 |
| miRNA | CCLE miRNA | **948** | 734 |

Overlap: **900** lines have all three; **726** have RNA only.

**They are not fused, and week 1 is what justifies that.** Agreement between
axes, measured on the 900 shared lines:

| pair | pair-level Spearman | top-25 neighbour overlap |
|---|---|---|
| RNA vs metabolomics | 0.32 | 9.7% |
| RNA vs miRNA | 0.23 | 14.6% |
| metabolomics vs miRNA | 0.11 | 6.9% |

Three genuinely different views. A fused score would also have been RNA alone for
the 726 lines that have nothing else, with no way to tell afterwards which lines
got which.

**Resolution** — a similarity that is 0.98 for everything ranks nothing:

| axis | median pair sim | sd | nearest-neighbour minus median |
|---|---|---|---|
| RNA | −0.044 | 0.201 | **0.725** |
| metabolomics | −0.006 | 0.226 | 0.584 |
| miRNA | −0.010 | 0.130 | 0.476 |

**Is it an artefact of the imputation?** Pearson vs Spearman (rank-based, far
less sensitive to imputed cells): pair-level ρ = 0.995 / 0.990 / 0.925 and top-25
overlap 90% / 83% / **61.5%** for RNA / metab / miRNA. RNA and metabolomics are
stable. **miRNA is the shakiest axis** and its neighbour lists move noticeably
with the metric — worth knowing before trusting a miRNA top-1.

---

## Week 2, test 1 — do similar lines behave similarly?

Dependency = DepMap 24Q4 CRISPR gene effect, restricted to **1,345 selective
genes** (sd > 0.20, not pan-essential). Including pan-essentials would make every
pair correlate for reasons unrelated to the graph and the test would pass no
matter what.

**Pooled** (the stated test) — median dependency-profile correlation:

| axis | top 0.1% | top 1% | all pairs | ρ(sim, dep) |
|---|---|---|---|---|
| RNA | 0.500 | 0.468 | 0.362 | 0.355 |
| metabolomics | 0.408 | 0.391 | 0.359 | 0.254 |
| miRNA | 0.379 | 0.395 | 0.359 | 0.105 |

All three pass. **That result is close to worthless on its own**, and here is
why: the top 1% of RNA pairs is **65.0% same-lineage** while the panel is 5.4%
same-lineage. Close pairs are mostly same-tissue pairs, and same-tissue lines
share dependencies for reasons a tissue label already gives you for free.

**So every comparison was run again within lineage** — close pairs against random
pairs *from the same lineage*:

| axis | same-lineage pairs | top 1% | all same-lineage | delta | p |
|---|---|---|---|---|---|
| RNA | 29,308 | 0.5173 | 0.4120 | **+0.1053** | 7.4e-72 |
| metabolomics | 13,910 | 0.4556 | 0.4108 | +0.0448 | 2.8e-09 |
| miRNA | 13,792 | 0.4604 | 0.4107 | +0.0497 | 1.2e-12 |

And against the free baseline directly:

| axis | random pair | same lineage (free) | graph top 1% | **graph − tissue** |
|---|---|---|---|---|
| RNA | 0.362 | 0.412 | 0.468 | **+0.056** |
| metabolomics | 0.359 | 0.411 | 0.391 | **−0.020** |
| miRNA | 0.359 | 0.411 | 0.395 | **−0.015** |

**This is the headline for week 2.** RNA beats knowing the tissue. Metabolomics
and miRNA do not — as a global ranking they are *worse* than simply picking a
same-tissue line at random, even though both carry real, highly significant
within-lineage signal.

That is not a reason to discard them. It is a reason to use them differently:
they resolve *within* a tissue, they do not rank *across* tissues. Week 3
presents them that way.

---

## Week 2, test 2 — the boring test with a known answer

Lineage from `cleaned_track_data/sample_info.parquet` (30-value controlled
vocabulary, not parsed free text).

| axis | top-25 purity | base rate | enrichment | AUROC (same vs diff lineage) |
|---|---|---|---|---|
| RNA | 56.4% | 5.8% | **9.8×** | **0.789** |
| metabolomics | 20.2% | 7.4% | 2.7× | 0.619 |
| miRNA | 23.8% | 7.2% | 3.3× | 0.602 |

The graph passes. It does **not** put lung next to blood at any meaningful rate.

**Per lineage (RNA)** — some tissues are recovered far better than others:

| well recovered | | poorly recovered | |
|---|---|---|---|
| blood, lymphocyte | ~81–86% | thyroid | 9% |
| skin, peripheral nervous system | ~81–83% | prostate, cervix | ~12–14% |
| kidney, colorectal | ~71–75% | bile duct, gastric | ~15% |
| lung | 58% | esophagus, urinary tract | ~19% |

**The failures are mostly real biology, which is evidence *for* the graph.**
27.2% of lines have a nearest neighbour from another lineage; the commonest
confusions are:

| actual | confused with | n |
|---|---|---|
| uterus | ovary | 11 |
| bile duct | pancreas | 9 |
| esophagus | upper aerodigestive | 7 |
| lymphocyte | blood | 7 |

Müllerian-duct tissues, pancreatobiliary tissues, contiguous squamous
epithelium, and two haematopoietic labels. These are developmental and
anatomical adjacencies, not errors. The label is more discrete than the biology.

`03` also names the genuinely bad cases: lines whose nearest neighbour is both a
lineage mismatch *and* low similarity (~0.25–0.35). Those are lines with **no
real neighbour**, and the query path flags them rather than confidently
returning a top-1.

> **Read this section against week 2 test 1.** A high score here is not
> straightforwardly good — an axis that reproduces the tissue label perfectly and
> adds nothing beyond it scores 1.0 and is worth nothing over a lookup. RNA
> earns its place by scoring well here *and* beating tissue in test 1.

---

## Week 3 — the query path

One route: **gene -> what is known -> candidate lines -> alternatives.**

```bash
python cell_similarity/04_query.py --gene EGFR --line ACH-000552 --k 10
```

Importable: `explore(gene, model_id, k)`, `describe_gene(q)`,
`candidate_lines(q, k, lineage=None)`, `similar_lines(model_id, k)`. `--json` for
machine output.

### What the tool claims, and what it does not

**Claims:** for a cell line, other lines that behave similarly *under drug
perturbation*. Validated at **+0.186 [+0.149, +0.227]** over a same-tissue
baseline, cluster-bootstrapped over cell lines, 649 lines and 150 compounds.
The claim is printed with the answer, not buried in a paper.

**Does not claim:** anything about whether knocking a gene out will kill a line.
That claim died in `gate_audit/` and is absent here on purpose.

So the gene view is **descriptive**. "Candidate lines" means *lines where this
gene is expressed* - an observation. Ranking them by predicted dependency would
reintroduce the dead claim under a new name, and the output says so explicitly:

```
2. CANDIDATE LINES  (EGFR)
  Basis: measured RNA expression, descending
  NOT a prediction: these lines express the gene; this is not a claim that
                    they depend on it

  1,220 of 1,673 measured lines express this gene. Top 5:
    line              log2(TPM+1)  lineage
    ach-001523              10.53  skin
    ach-000741              10.08  urinary_tract
```

### Axis roles are printed, not implied

```
3. ALTERNATIVE LINES  ach-001523   lineage: skin
  Validated on GDSC2 drug response (150 compounds): +0.186 [0.149, 0.227]
  over a same-tissue lookup, 649 lines, cluster-bootstrapped over cell lines.

  [rna]  RNA expression -- VALIDATED RANKER
        +0.186 drug-response concordance over a same-tissue baseline
    rank line                 sim  lineage
       1 ach-001695         0.723   unknown
       2 ach-000832         0.703   upper_aerodigestive

  [metab]  NOT COVERED -- this line has no metabolomics data (928 lines do have it)
  [mirna]  NOT COVERED -- this line has no miRNA data (948 lines do have it)
```

| axis | role in the UI | why |
|---|---|---|
| RNA | **validated ranker** | +0.186 vs drug response over tissue; split-half 83.7% |
| metabolomics | within-tissue context | +0.008 n.s.; split-half 42.9% |
| miRNA | within-tissue context | +0.019 n.s.; split-half 58.8% |

Four rules, each inherited from something measured:

1. **Absence is reported as absence.** `NOT COVERED - this line has no
   metabolomics data (928 lines do have it)`, never an empty list. This is the
   friendly form of the mistake `gate_audit` found, where unmeasured lines were
   counted as confirmed negatives.
2. **Axes stay separate** - top-25 overlap is 7-15%, so a merged list would be a
   fiction. There is a "found by more than one axis" section instead.
3. **Each axis is labelled with what it is worth**, so the two weak ones cannot
   be mistaken for the validated one.
4. **No dependency prediction anywhere.**

Abstention follows `docs/QUERY_LAYER_STATES.md`: `NO_EVIDENCE`, `UNINFORMATIVE`
(verified: `OR2T1` -> 0.0% expressed, `TAS2R43` -> 0.1%), `UNKNOWN_GENE`.

---

## Robustness (`05`, `06`) - intervals, RNAi, and the housekeeping floor

Everything above was point estimates on non-independent pairs. Three checks.

### Cluster bootstrap over LINES, not pairs

Each line appears in ~1,600 pairs, so pair-level statistics treat one line's
idiosyncrasies as ~1,600 observations. The resampling unit is the **cell line**;
pairs formed between two draws of the same line are dropped (they would be
similarity 1.0 by construction).

### All three checks, one table

`graph - tissue`, with 95% CI. `*` = interval excludes zero.

| axis | dependency | lines | random | tissue | graph | **graph - tissue** | 95% CI |
|---|---|---|---|---|---|---|---|
| RNA | CRISPR selective | 1,046 | 0.362 | 0.412 | 0.469 | **+0.0570** | [+0.044, +0.066] * |
| RNA | CRISPR, no common ess. | 1,046 | 0.280 | 0.337 | 0.417 | **+0.0802** | [+0.069, +0.092] * |
| RNA | RNAi selective | 659 | 0.527 | 0.585 | 0.591 | +0.0053 | [-0.014, +0.032] |
| RNA | **RNAi, no common ess.** | 659 | 0.245 | 0.308 | 0.328 | **+0.0201** | [+0.001, +0.042] * |
| metab | CRISPR selective | 649 | 0.359 | 0.411 | 0.391 | -0.0203 | [-0.033, -0.007] * |
| metab | CRISPR, no common ess. | 649 | 0.279 | 0.338 | 0.326 | -0.0117 | [-0.023, +0.002] |
| metab | RNAi, no common ess. | 574 | 0.248 | 0.317 | 0.253 | -0.0631 | [-0.083, -0.042] * |
| miRNA | CRISPR selective | 646 | 0.359 | 0.411 | 0.412 | +0.0007 | [-0.014, +0.012] |
| miRNA | RNAi, no common ess. | 573 | 0.249 | 0.317 | 0.277 | -0.0401 | [-0.064, -0.015] * |

**1. The +0.056 survives.** +0.0570, CI [+0.044, +0.066]. It was not an artefact
of pair pseudo-replication.

**2. Stripping common essentials makes it bigger, as suspected.** The random-pair
floor drops from 0.362 to 0.280 - that floor was shared housekeeping - and the
advantage grows to **+0.0802** [+0.069, +0.092]. The graph now moves 19.1% of the
available headroom from floor toward 1.0, against 16.7% before.

**3. The RNAi test - the deciding one - is positive but weak.** At 5,000
replicates: **+0.0201, 95% CI [+0.0011, +0.0423], bootstrap p = 0.020**. It
excludes zero, but only just; the CI straddled zero at 400 replicates. Two things
are true and both must be reported:

- It is **not purely** the expression-CRISPR artefact. The advantage is still
  positive on a perturbation that lacks that confound.
- But it is **0.25x the CRISPR magnitude**. Only a quarter of the CRISPR-measured
  advantage replicates, so the +0.080 figure is substantially inflated by
  something CRISPR-specific.

**The defensible estimate of what the graph adds is the RNAi number, +0.020, not
the CRISPR one.** The RNAi arm has 659 lines against 1,046 and a noisier readout,
so this is underpowered rather than negative - but it is not established, and this
component should not be described as having cleanly shipped on it.

### `06` - miRNA, and a surprise about metabolomics

Diagnosis of miRNA's 61.5% Pearson/Spearman instability: **participation ratio
13.0**. Of 734 miRNAs, ~13 effectively carry the similarity, with the top 10
holding 72.7% of the variance (RNA: 1,648 and 1.7%). Median skew 3.54, 65.9% of
features with |skew| > 2. The graph was resting on about a dozen numbers.

The Pearson-vs-Spearman check turned out to be **circular** for the obvious fix:
rank-within-line scores 84.2% there because ranking twice is near-idempotent. So
`06` adds **split-half feature reliability** - halve the features, build the graph
twice, compare neighbour lists. No transform can game it:

| axis / variant | split-half top-25 agreement |
|---|---|
| RNA | **83.7%** |
| **metabolomics** | **42.9%** |
| miRNA as-built | 51.9% |
| miRNA **log1p** | **58.8%** |
| miRNA rank-within-line | 51.4% (confirms the 84.2% was circular) |

Two consequences:

- **miRNA is rebuilt with log1p.** It improves the axis materially: `graph -
  tissue` moves -0.0155 to +0.0007 (CI [-0.014, +0.012], now indistinguishable
  from zero rather than negative), and tissue purity 23.8% to 29.2%, AUROC 0.602
  to 0.649.
- **Metabolomics is the least reproducible axis (42.9%), worse than miRNA.**
  Week 1's reassuring 82.8% Pearson/Spearman agreement was measuring metric
  invariance, not graph reliability - different things, and only the second
  matters. This was not asked for and is the more concerning of the two.

Neither secondary axis should carry a ranking.

### Why metabolomics and miRNA behave this way - the mechanism

Better within lineage, worse across it, is coherent rather than contradictory.
Global metabolic geometry is dominated by growth rate, culture medium and batch -
factors that vary across the whole panel and swamp tissue identity. Within a
tissue those are more nearly constant, so real biology shows through. That is why
both axes carry significant within-lineage signal (+0.045, +0.050, p < 1e-8) and
still lose to a plain tissue lookup as global rankers.

**They refine; they do not rank.** That is a mechanism rather than an
observation, and it is why the query path labels them `within-tissue context`.

---

## Reliability, and a third arm (`07`, `08`)

The 0.25x CRISPR/RNAi ratio has an exact statistical reading rather than a
rhetorical one. If each dependency profile carries reliability `r_xx`, the
observed correlation - and so `graph - tissue` - is attenuated by it, and
`delta_true = delta_observed / r_xx`. The question "power or confound" therefore
reduces to a measurable quantity.

### `07A` - RNAi reliability, measured

DEMETER2 is fitted to separate source screens. 226 lines were screened in both
Achilles shRNA and Novartis DRIVE, giving a direct test-retest of the *same*
perturbation type - which is the clean comparison, since CRISPR-vs-RNAi
confounds reliability with mechanism.

| | value |
|---|---|
| per-line profile test-retest, single screen | **r = 0.144** (IQR 0.109-0.180) |
| Spearman-Brown to the 2-screen combined estimate | **r_xx = 0.252** |

RNAi profiles are extremely noisy. Disattenuating the RNAi arm alone gives
+0.0201 / 0.252 = **+0.0797**, or 0.99x the CRISPR +0.0802 - which looks like a
complete explanation.

### `07D` - the symmetry check, which overturns that

Correcting one arm and not the other rigs the comparison: CRISPR's reliability
is not 1.0 either. Measured the same way, on 176 lines screened in both Broad
Achilles and Sanger Project Score:

| arm | r_xx | observed | **disattenuated** |
|---|---|---|---|
| CRISPR (24Q4, single screen) | 0.371 | +0.0802 | **+0.2162** |
| RNAi (DEMETER2, 2 screens combined) | 0.252 | +0.0201 | **+0.0797** |

**Ratio of disattenuated estimates: 0.37x.** The gap shrinks under correction but
does not close. Attenuation explains a large part of it; it does not explain all
of it.

### `07B` - the dose-response, which is the strongest evidence for attenuation

DEMETER2 publishes per-observation posterior SDs. Restricting to
precisely-estimated genes should raise the delta if attenuation is the
explanation, and leave it flat if a confound is:

| gene subset | genes | delta | 95% CI |
|---|---|---|---|
| all | 6,640 | +0.0201 | [-0.000, +0.047] |
| lowest 50% posterior SD | 3,320 | +0.0308 | [+0.011, +0.057] |
| lowest 25% posterior SD | 1,660 | +0.0405 | [+0.020, +0.067] |
| lowest 10% posterior SD | 664 | **+0.0556** | [+0.030, +0.080] |

Monotone, and nearly triples. Measurement error is demonstrably suppressing the
RNAi arm.

### `07C` - the power control

| arm | lines | delta | 95% CI |
|---|---|---|---|
| CRISPR, all lines | 1,046 | +0.0820 | [+0.068, +0.093] |
| CRISPR, restricted to the RNAi lines | 541 | **+0.0870** | [+0.068, +0.103] |
| RNAi | 659 | +0.0201 | [+0.001, +0.047] |

Sample size is **not** the explanation - CRISPR on the same lines is if anything
slightly stronger.

### `08` - the drug arm, and it is the strongest result here

GDSC2, 880 lines x 150 compounds, AUC, **centred per drug** so that "resistant to
everything" does not read as similar (the drug analogue of the common-essentials
floor). Profile = a line's response across the compound panel.

| axis | lines | random | tissue | graph | **graph - tissue** | 95% CI |
|---|---|---|---|---|---|---|
| **RNA** | 649 | -0.005 | 0.093 | 0.279 | **+0.1856** | [+0.150, +0.222] * |
| metabolomics | 566 | 0.000 | 0.097 | 0.105 | +0.0082 | [-0.032, +0.053] |
| miRNA | 567 | -0.002 | 0.096 | 0.115 | +0.0187 | [-0.027, +0.072] |

Split-half over drugs gives r_xx = 0.707, so the drug readout is far more
reliable than either genetic one.

**Quote it as: +0.186 observed, +0.262 after correcting for a readout reliability
of 0.707.** Always paired, and the observed number carries the claim. At r_xx =
0.707 the correction changes little and disattenuated correlations are an easy
thing to attack, so there is no reason to rest any weight on a correction this
argument does not need. The correction earns its place in the RNAi comparison
(`07`), where reliability is 0.252 and ignoring it would be the error - not here.

### The three arms on a common scale

| readout | failure mode | r_xx | **observed** | disattenuated | correction needed? |
|---|---|---|---|---|---|
| CRISPR | expression / copy-number cutting toxicity | 0.371 | +0.080 | +0.216 | yes - low reliability |
| **Drug (GDSC2)** | **polypharmacology** | **0.707** | **+0.186** | +0.262 | **no - observed stands** |
| RNAi | knockdown efficiency, seed effects | 0.252 | +0.020 | +0.080 | yes - and still short |

The disattenuated column is secondary throughout. It exists to put three readouts
of very different quality on one scale, which is the only way the CRISPR/RNAi
comparison can be read at all. It is not needed for the drug arm and no claim
here depends on it.

**This changes the reading.** CRISPR and drug response - two readouts with
completely different failure modes - agree closely once reliability is accounted
for (+0.216 vs +0.262). RNAi is the outlier, not CRISPR. If the CRISPR advantage
were the expression artefact from `gate_audit`, drug response should not
reproduce it, because drug response is not expression-confounded in that way. It
does reproduce it, and on the most reliable readout of the three.

The residual RNAi shortfall is most likely that disattenuation assumes classical
random measurement error, while RNAi's dominant error is *systematic* (seed-based
off-target effects), which correction of this form cannot repair.

### `08C` - the headline survives the lineage-shortcut checks

Per-drug centring handles "resistant to everything". The remaining worry was that
drug response is strongly lineage-structured, so the top 1% might win on a
haematological/solid split rather than on resolution.

| subset | lines | tissue | graph | **graph - tissue** | 95% CI |
|---|---|---|---|---|---|
| all lines (headline) | 649 | 0.093 | 0.279 | **+0.1856** | [+0.149, +0.227] |
| solid tumours only, haem dropped | 542 | 0.069 | 0.259 | **+0.1901** | [+0.139, +0.244] |
| within-lineage pairs only (15,713 pairs) | - | 0.093 | 0.345 | **+0.2519** | - |

The effect is *strongest* within lineage, so the split is not carrying it.

**Provenance, confirmed:** the random / tissue / graph triple all come from one
call on the same matrices, `D` built from the 880 x 150 drug panel and sliced to
the same line set. Nothing is carried from the CRISPR arm - visible in the
numbers themselves, since the drug tissue baseline is 0.0932 against the CRISPR
arm's 0.337.

**The CI was already a cluster bootstrap over lines**, not pair-level: the
estimator resamples cell lines with replacement and drops pairs formed between
two draws of the same line. Re-run at 2,000 replicates: [+0.1492, +0.2270].

### Metabolomics - a fair attempt, and a clean reason to hold it

| variant | split-half reliability |
|---|---|
| as-built | 43.4% |
| log1p (shifted) | 42.7% |
| rank-within-line | 42.5% |
| winsorised 1-99% | 43.6% |

No transform moves it. And the diagnosis differs from miRNA's: participation
ratio is **115.2 of 225** features with the top 10 holding only 18.0% of variance,
so it is *not* the low-dimensional fragility that log1p fixed for miRNA. 225
metabolites is simply not enough signal to estimate a stable 928-line neighbour
graph. That is an intrinsic limit rather than a preprocessing choice, which makes
it a defensible reason to hold the axis rather than an unexplained one.

---

## What this is worth, stated plainly

**RNA similarity beats the tissue baseline on three readouts with three
different failure modes.** Observed, which is what the claim rests on: drug
response **+0.186 [+0.149, +0.227]**, CRISPR +0.080 [+0.069, +0.092], RNAi
+0.020 [+0.001, +0.047].

Secondarily, and only so three readouts of very different quality can be compared
on one scale: after correcting each for its own measured reliability (0.707,
0.371, 0.252) they become +0.262, +0.216, +0.080. **No claim here needs that
correction** - the drug arm's reliability is high enough that the observed
number stands on its own.

The two most reliable readouts agree, and they fail in unrelated ways - CRISPR
through expression-linked cutting toxicity, drug response through
polypharmacology. No single artefact explains both. The earlier worry that the
CRISPR advantage was the `gate_audit` expression confound wearing a new hat does
not survive the drug arm, which is not expression-confounded in that way and
gives the largest effect of the three.

**Quote the drug arm, observed: +0.186 [+0.149, +0.227].** It is the most
reliable (r_xx = 0.707 against 0.371 and 0.252), it is the largest, and it is the
one that maps onto the actual product claim: choose a backup line from this graph
and a compound is more likely to behave the same way in it. Pair it with the
corrected value if the reliability question is raised; do not lead with it.

**The lineage confusions are the most persuasive evidence here, and they need no
p-value.** Uterus/ovary is Müllerian. Bile duct/pancreas is pancreatobiliary.
Oesophagus/upper aerodigestive is contiguous squamous. Lymphocyte/blood is one
thing labelled twice. A graph that fails only where the ontology is more discrete
than the biology is a graph that is working.

**Metabolomics and miRNA are not global rankers, and there is now a mechanism for
why.** Both are exposed as within-tissue context. On split-half reliability
metabolomics is the weakest axis at 42.9%, below miRNA's rebuilt 58.8% and far
below RNA's 83.7% - so neither should carry a ranking, and the concern miRNA's
instability raised turned out to apply more strongly to metabolomics.

**Nothing is fused, and the case for fusing has not been made.** The axes agree
weakly and cover different lines. If fusion is attempted later, the test it must
pass is the one in `02`: does the fused score beat RNA alone on `graph − tissue`?
Anything less and it's added complexity for nothing.

### Wired into the pipeline query layer

`src/pipeline/cell_similarity_bridge.py` makes the validated claim reachable from
the pipeline's own query surface, next to `rank_cell_lines()` / `explain_pair()`.
Additive and read-only: it does not modify `core_score`, does not re-sort any
ranking, and writes nothing.

```bash
python src/pipeline/cell_similarity_bridge.py --gene EGFR --top-n 10
```

```
GENE ensg00000146648   class=activation_driven   candidates=1674
  RANKING BASIS : core_score -- ABUNDANCE ONLY
    an abundance / expressibility prior (RNA + protein percentile)
    NOT a dependency estimate (gate_audit/03, gate_audit/09)
  ALTERNATIVES  : VALIDATED
    +0.186 [0.149, 0.227] over a same-tissue lookup on GDSC2 drug response

  #3    ach-001328     abundance=0.987  conf=high   skin squamous cell ca
        alternatives: ach-001690(0.72), ach-000606(0.71), ach-001694(0.70)
  #4    ach-000452     abundance=0.981  conf=high   esophageal squamous c
        alternatives: ach-000873(0.39), ...  << low confidence
```

**The provenance problem this had to solve.** `rank_cell_lines()` sorts by
`core_score`, and `gate_audit/` withdrew the *dependency* interpretation of that
column. Attaching a validated claim to an invalidated one without labelling would
let the pair read as equally supported. So every result carries a `claims` block
naming what each half rests on:

| part of the answer | status | rests on |
|---|---|---|
| ranking order (`core_score`) | **ABUNDANCE ONLY** | an expressibility prior; the dependency claim is withdrawn |
| `alternatives` per line | **VALIDATED** | +0.186 [+0.149, +0.227] vs drug response |

`core_score` is still an abundance ordering, which is what it always measured -
it is the dependency reading that is withdrawn, and the caveat says that rather
than implying the column is meaningless.

Only RNA is wired in. Metabolomics and miRNA are named in
`alternatives_coverage.context_axes_not_wired` with the reason (+0.008 and
+0.019, both n.s.), so their absence is explicit rather than silent. Uncovered
lines return `covered: false` with the axis size, never an empty list.

API: `attach_alternatives(result, k, top_n)`, `similar_lines_for(model_id, k)`,
`graph_available()`.

### Not done

- **Similarity is computed on all genes, not gene-conditioned.** "Which lines are
  similar *for this gene's pathway*" is a different and probably more useful
  question. This graph answers the global version.
- **No confidence intervals on the neighbour lists.** A top-1 at 0.727 and a
  top-1 at 0.402 are presented with the same structure; only a hard threshold
  flag distinguishes them.
- **The 27.2% nearest-neighbour lineage mismatch is not decomposed** into "real
  biology" versus "genuinely wrong". The confusion table suggests most is the
  former, but that's an eyeball judgement, not a measurement.
