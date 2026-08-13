# TRANSCRIPTOMICS_CORRELATION — measured redundancy and the disagreement flag

Parts A and D of the transcriptomics brief. Artefact under review: `02_transcriptonomics.ipynb`
(supplied 2026-08-09 from `C:\Users\vigne\Downloads\`; **still not committed to the repo**).
Measurements run against `src/pipeline/outputs/celllineselector.db`, the database the
notebook itself opens.

Scripts: `docs/audit_scripts/tx_*.py`, `a2_source_corr.py`, `a3_within_geo.py`, `a4_inflation.py`.

Tags: `VERIFIED` (computed here) · `PROVEN` (published) · `UNVERIFIED — could not locate`.

---

## Credit where it is due

Three things in this notebook are right, and one of them is right in a way the main
pipeline is not.

1. **The stratification solves a real problem the main pipeline had not identified.**
   Pooling every lineage mixes populations where a gene is on with populations where it
   is off. The main pipeline's `core_score` does exactly this and has no guard against it.
2. **The replicate principle is correctly stated and correctly applied at the profile
   level.** Several DepMap profiles of one cell line are collapsed to a median with
   `n_samples` carried as a √n precision weight, explicitly because "counting them as two
   sources would inflate the Stouffer combination, which assumes sources are INDEPENDENT"
   (cell 8 docstring). That is the correct principle, stated in the notebook's own words.
3. **The floored-source distinction is directly useful to the main pipeline.** Excluding a
   detection-floored source from the *silence* judgement while keeping its z-scores is a
   real and non-obvious insight, and it bears on §9 gap 13 — the main pipeline's most
   serious open item.

The fault in Part A is that principle (2) is applied one level down and not one level up.
The fix is the notebook's own stated principle, applied consistently.

---

## A1. Version confirmation — the flow diagram describes a notebook that does not exist

| Feature the flow diagram describes | Present? | Evidence |
|---|---|---|
| Per-source detection floor check | **PRESENT** | `source_floor_check()`, cell 5; invoked cell 18 |
| GEO processing-method column | **ABSENT** | `fetch_gene_expression` (cell 3) reads `geo_expr` and collapses `median(gid) … GROUP BY model_id`. No method column is read, derived or stored |
| Stratum `source × method × lineage` | **ABSENT** | the stratum is `source × lineage` — `for (src, lin), g in d.groupby(["source", "lineage"])`, cell 8 |
| `harvest_geo_methods` step | **ABSENT** | `UNVERIFIED — could not locate`. No function of that or any similar name exists in the notebook |
| Combination step taking a correlation matrix | **ABSENT** | `stouffer(zs, weights)` computes `Z = (w*z).sum()/np.sqrt((w**2).sum())`, cell 8 — the independence formula. No `R` argument, no covariance term |

**The supplied notebook is the `source × lineage` version.** Standing rule 1 applies: the
measurement wins over the diagram. Everything below measures the notebook that exists.

### The method split has nothing to split on — VERIFIED

There is no processing-method column anywhere, and the underlying GEO holding is
homogeneous:

- `geo_info.platform_id` is populated for 2,329 of 3,267 GSMs and is **`gpl570`
  (Affymetrix HG-U133 Plus 2.0) for every one of them**. `type` is `rna` for all 2,329.
- Processing method is not recorded in any field. Inferred instead from the value
  distributions, per series (60-gene sample):

| gse_id | n GSM | min | median | max | inferred scale |
|---|---|---|---|---|---|
| (no gse_id) | 938 | 3.49 | 90.06 | 11,983 | linear, MAS5-like |
| gse57083 | 627 | 5.15 | 70.21 | 8,545 | linear, MAS5-like |
| gse50811 | 241 | 5.22 | 56.17 | 5,023 | linear, MAS5-like |
| gse34211 | 230 | 6.61 | 64.85 | 8,341 | linear, MAS5-like |
| gse10843 | 206 | 13.28 | 109.73 | 7,938 | linear, MAS5-like |
| …all 19 series… | | 4.17–18.95 | 39.21–171.37 | 3,771–11,983 | linear, MAS5-like |

**All 19 series are on the same linear, MAS5-like scale.** RMA output would be log2-scaled
(≈2–14). One platform, one processing family, **no heterogeneity to stratify on**.

Two consequences:

- The "MAS5 stratum in the single digits" the brief anticipates does not exist as a small
  stratum. MAS5-like is essentially 100% of the GEO holding.
- **A correction to the notebook's own text.** Cells 4 and 37 state GEO is "much of it
  RMA-normalised microarray, which floors at background". The data says MAS5-like. The
  *floor conclusion survives* — MAS5 background-corrects to a positive floor, and one is
  measured (GEO min 3.03, **0.00%** of values below the expressed threshold, against
  15.45% for DepMap and 17.99% for HPA) — but the stated mechanism is wrong and should be
  corrected to "background-corrected single-channel array (MAS5-like)".

A download of the 18 GSE `!Sample_data_processing` fields from NCBI would confirm the
submitters' own wording. It was **not performed**: the numeric signature is stronger
evidence than a free-text field, and the answer changes no decision — whether it is MAS5
or gcRMA, there is one method.

---

## A2. The source correlation matrix

Measured within gene, across shared cell lines, on 400 randomly sampled genes present in
all three warehouse tables (seed 42), using the notebook's own transforms
(`log2(x+1)` on GEO and HPA nTPM; DepMap is already `log2(TPM+1)`).

### Raw values — confirms the main pipeline audit

| Pair | Shared lines | Genes | Median ρ | IQR | Frac ρ > 0.5 |
|---|---|---|---|---|---|
| DepMap–HPA | 1,041 | 400 | **0.829** | 0.743 – 0.912 | 0.935 |
| DepMap–GEO | 525 | 400 | 0.459 | 0.263 – 0.605 | 0.428 |
| HPA–GEO | 513 | 399 | 0.508 | 0.296 – 0.640 | 0.514 |

Against the main-pipeline audit's 0.817 / 0.467 / 0.523 — agreement to within 0.015.
Both `VERIFIED`, independently, on different line sets and different transforms.

### Z-scores — the quantity the Stouffer combination actually assumes is zero

Robust z within `(source, lineage)`, with the notebook's own guards applied
(`MIN_PEERS_FOR_Z = 3`, `SILENT_FRAC = 0.20`, `MAD_FLOOR = 0.1`).

| Pair | Shared lines | Genes | Median ρ | IQR | Frac ρ > 0.5 |
|---|---|---|---|---|---|
| DepMap–HPA | 1,041 | 351 | **0.863** | 0.777 – 0.918 | **0.992** |
| DepMap–GEO | 525 | 351 | 0.495 | 0.410 – 0.601 | 0.490 |
| HPA–GEO | 513 | 343 | 0.530 | 0.443 – 0.631 | 0.583 |

**R, the matrix the corrected combination needs:**

```
          DepMap     HPA     GEO
DepMap    1.0000  0.8631  0.4953
HPA       0.8631  1.0000  0.5303
GEO       0.4953  0.5303  1.0000
```

### The difference — and a correction to the brief's expectation

The brief predicts the z correlation will be *higher* than the raw correlation, because
stratified z removes platform-specific variance. **On the marginal comparison it is; on
the paired comparison it is not.** Both are reported because the difference between them
is a selection effect, not a contradiction.

| Pair | Genes | Raw median | Z median | Raw IQR | Z IQR | Median paired Δ | Frac z higher |
|---|---|---|---|---|---|---|---|
| DepMap–HPA | 351 | 0.842 | **0.863** | 0.772–0.920 | 0.777–0.918 | **−0.024** | 23.7% |
| DepMap–GEO | 351 | 0.495 | 0.495 | 0.358–0.625 | 0.410–0.601 | −0.011 | 45.0% |
| HPA–GEO | 343 | 0.544 | 0.530 | 0.403–0.660 | 0.443–0.631 | −0.028 | 36.4% |

Gene-for-gene, most genes see a *small decrease*. The marginal median rises because the
guards remove degenerate genes, which are the low-correlation ones. What the z transform
reliably does is **narrow the IQR** — DepMap–GEO's lower quartile moves from 0.358 to
0.410 — and push `frac ρ > 0.5` for DepMap–HPA from 0.994 to **0.992** (already saturated).

**The finding that matters: stratified z-scoring does not reduce the redundancy at all.**
Median cross-source correlation is essentially unchanged, and for DepMap–HPA it is
*higher* on the surviving gene set. The redundancy is in the biology, not in the platform
scaling. It cannot be normalised away, which is precisely why it has to be corrected for
in the combination step.

---

## A3. Within-GEO correlations — the sharp one

There is no method column, so the honest proxy for a processing stratum is the GEO
**series** (`gse_id`), since processing is a per-submission property. This measures what
the proposed split would actually do.

### One cell line, many GEO series — VERIFIED

| Series a line appears in | Lines |
|---|---|
| 1 | 314 |
| 2 | 118 |
| 3 | 83 |
| 4 | 44 |
| 5 | 15 |
| 6 | 9 |
| 7 | 6 |
| **9** | **1** |

**276 of 590 GEO cell lines (46.8%) appear in ≥2 series. 158 appear in ≥3. One appears in 9.**

### Series stratum sizes — VERIFIED

| gse_id | lines | | gse_id | lines |
|---|---|---|---|---|
| gse57083 | 362 | | gse50830 | 16 |
| (no gse_id) | 237 | | gse65216 | 14 |
| gse34211 | 144 | | gse7127 | 14 |
| gse10843 | 120 | | gse23806 | 13 |
| gse12790 | 51 | | **gse84557** | **8** |
| gse10890 | 44 | | **gse14315** | **5** |
| gse15329 | 35 | | **gse50451** | **3** |
| gse50811 | 26 | | **gse30240** | **2** |
| gse53798 | 23 | | | |
| gse50831 | 20 | | | |
| gse41445 | 18 | | | |

**7 of 19 series fall below `MIN_PEERS_FOR_LINEAGE = 15`. Four are below 10. Two
(gse30240 n=2, gse50451 n=3) are at or below `MIN_PEERS_FOR_Z = 3`** — they cannot
support a spread estimate at all. The brief's concern about single-digit strata is
correct; it simply attaches to series rather than to methods.

### Correlation between series-specific z-scores — VERIFIED

Pooled over (line, gene) pairs for lines shared by both series; 40 series pairs with
≥5 shared lines and ≥50 pairs.

| Statistic | Value |
|---|---|
| Series pairs measured | 40 |
| **Median Spearman** | **0.448** |
| IQR | 0.388 – 0.510 |
| Pair-count-weighted mean | 0.471 |
| Fraction of pairs with ρ > 0.5 | 32.5% |
| Highest pair | gse10843–gse12790, ρ = 0.835 (27 lines) |
| Lowest pair | (no gse_id)–gse15329, ρ = 0.149 (25 lines) |

### Are these independent observations, or one observation counted several times?

**They are repeat measurements of the same cell line.** Stated plainly:

The notebook already answers this for the level below. Its own cell-8 docstring says
several DepMap profiles of one line "are technical replicates, not independent
evidence … counting them as two sources would inflate the Stouffer combination". Two GEO
series containing the same cell line are the **same category of thing** — repeat
measurements of one biological entity, differing in submitter, passage and batch, not in
what they are measuring.

A ρ of 0.448 is *lower* than DepMap–HPA's 0.863, and it is tempting to read that as "more
independent". That reading is wrong. Low correlation between two measurements of the same
entity means **high measurement noise**, not independent evidence. Stouffer would credit
two noisy readings of one cell line as two confirmations.

**The current version already handles this correctly** — `fetch_gene_expression` collapses
GEO to `median(value) … GROUP BY model_id` with `count(*) AS n_samples` carried as a √n
precision weight. That is the notebook's own stated principle, correctly applied.
**Splitting by method or series would abandon it.** The split is an amplification of the
existing fault, not a correction. Quantified in A4.

---

## A4. Quantifying the inflation

```
current      Z = Σ w z / √(Σ w²)        [assumes independence]
corrected    Z = Σ w z / √(w' R w)      [Strube 1985 / Hartung 1999]   PROVEN
```

The ratio `Z_current / Z_corrected = √(w'Rw) / √(w'w)` **depends only on which sources are
present and their weights — not on the data**. So it is exact per source-presence pattern.

### Exact inflation by pattern (equal weights) — VERIFIED

| Sources | k | Inflation ratio | Effective independent sources | Overstated by |
|---|---|---|---|---|
| DepMap + HPA | 2 | **1.365** | **1.07** of 2 | **36.5%** |
| DepMap + GEO | 2 | 1.223 | 1.34 of 2 | 22.3% |
| HPA + GEO | 2 | 1.237 | 1.31 of 2 | 23.7% |
| **DepMap + HPA + GEO** | 3 | **1.503** | **1.33** of 3 | **50.3%** |

### The five-input case, if GEO were split by series — VERIFIED

Using the measured within-GEO series correlation ρ = 0.448 and three GEO series:

| Configuration | Inflation ratio | Effective sources | Overstated by |
|---|---|---|---|
| 3-source (GEO unsplit, current) | 1.503 | 1.33 of 3 | 50.3% |
| **5-input (GEO split into 3 series)** | **1.765** | **1.61 of 5** | **76.5%** |

**Splitting GEO into three series makes the overstatement worse by 26.2 percentage
points.** It buys 0.28 of an effective independent source while adding 26 points of
false confidence. This is the decisive number for the method-split decision.

### Empirical distribution across real cells — VERIFIED

200 genes, real `√n_samples` weights, 1,106 lines with ≥2 sources:

| Pattern | Lines | Median ratio | IQR |
|---|---|---|---|
| DepMap+HPA | 561 | 1.365 | 1.365–1.365 |
| DepMap+HPA+GEO | 478 | 1.469 | 1.379–1.503 |
| DepMap+GEO | 39 | 1.223 | 1.223–1.223 |
| HPA+GEO | 28 | 1.237 | 1.208–1.237 |

**Overall: median inflation 1.365, IQR 1.365–1.434, max 1.503. Median effective
independent sources: 1.10.**

### How much of the reported confidence is evidence?

Stated explicitly, as the brief asks:

> For the typical scored line — two sources, DepMap and HPA — the combined Z is
> **36.5% larger** than the evidence supports. The two sources are worth **1.07**
> independent observations, not 2. For a three-source line the Z is **50.3%** too large
> and the three sources are worth **1.33**, not 3.
>
> Across all lines with ≥2 sources, the median line's evidence base is **1.10 independent
> observations**. Almost everything the second and third source appear to add is the same
> measurement counted again.

This does **not** mean the scores are worthless — the ordering within a gene is barely
affected, because the inflation factor is near-constant within a source-presence pattern.
It means the **p-values, q-values and the "Conflicting Sources"/significance verdicts
built on them are not calibrated**, which Part B measures directly.

---

# PART D — The disagreement flag

`heterogeneity()` computes Cochran's Q and I², and `sources_agree = I2 < I2_CUT` with
`I2_CUT = 50`, driving the `"Conflicting Sources"` verdict.

## How often does it fire, and on how few sources? — VERIFIED

Eight genes, all 1,580 lines each:

| Gene | Lines with ≥2 sources | Conflicting | % of scoreable | on 2 sources | on 3 sources |
|---|---|---|---|---|---|
| TSPAN6 | 914 | 16 | 1.75% | 11 (68.8%) | 5 |
| EGFR | 914 | 14 | 1.53% | 3 (21.4%) | 11 |
| TP53 | 1,106 | 31 | 2.80% | 10 (32.3%) | 21 |
| MYC | 1,106 | 15 | 1.36% | 3 (20.0%) | 12 |
| BRAF | 1,106 | 40 | 3.62% | 24 (60.0%) | 16 |
| PTEN | 1,106 | 65 | 5.88% | 43 (66.2%) | 22 |
| KRAS | 1,106 | 35 | 3.16% | 18 (51.4%) | 17 |
| CDKN2A | 1,106 | 51 | 4.61% | 27 (52.9%) | 24 |

**Pooled: 267 of 8,464 scoreable lines (3.15%) receive the verdict. 139 of those (52.1%)
have only TWO sources.** More than half of all disagreement verdicts rest on a statistic
computed from a single degree of freedom.

## How biased is I² at these k? — VERIFIED

200,000 draws of k independent N(0,1) z's — i.e. **no true heterogeneity whatsoever**:

| k | Mean I² | Median I² | P(I² > 50) |
|---|---|---|---|
| **2** | **15.1%** | 0.0% | **15.75%** |
| **3** | **14.9%** | 0.0% | **13.51%** |
| 5 | 13.4% | 0.0% | 9.02% |
| 7 | **12.5%** | 0.0% | 6.18% |

The k=7 figure reproduces the ~12 percentage-point overstatement the brief cites
(`PROVEN`, Higgins & Thompson 2002 and the subsequent small-k literature), and confirms
it is **worse below that** — 15.1% at k=2.

Two further observations:

- The **median** I² is 0.0% at every k while the **mean** is 12–15%. The statistic is not
  merely biased, it is bimodal at small k: usually exactly zero, occasionally large. A
  fixed threshold on such a statistic is a coin-flip on the occasions it fires.
- The observed fire rate (3.15%) is **far below** the 15.75% a k=2 null predicts. That is
  not reassurance — it is another measurement of the redundancy. Real sources agree far
  more than independent normals because ρ = 0.863, so I² rarely clears 50. The flag is
  simultaneously **biased upward** by small k and **suppressed** by source correlation,
  and those two errors do not cancel in any controlled way.

## Recommendation

Of the three defensible options, take the **second, with a fallback to the third**:

1. **Use the underlying test's p-value instead of the point statistic.** `heterogeneity()`
   already computes `pq = 1 - chi2.cdf(Q, k-1)` and **discards it** — only Q and I² are
   returned into the row. Flag on Cochran's Q p-value, which at least has a defined null
   distribution, and report it. This is a two-line change to what is already computed.
2. **If the flag is kept as-is, relabel it.** `"Conflicting Sources"` reads as a test
   result. It is a heuristic warning. Rename to `"sources may disagree (heuristic)"` and
   state `I2_CUT = 50` and k in the output.
3. **Do not put an interval on I² at k=2.** An interval on a statistic with one degree of
   freedom will span nearly the whole range and will mislead more than the point estimate.
   The interval option is defensible at k≥5 and not here.

**A silent threshold is not defensible**, and that is what currently ships.

**Free parameter.** `I2_CUT = 50` is Higgins & Thompson's conventional cut (`PROVEN` as a
convention) but is applied at k=2–3 where the convention was never intended to operate. If
it stays, sweep it over {25, 50, 75} with held-out selection, and report the fire rate at
each.
