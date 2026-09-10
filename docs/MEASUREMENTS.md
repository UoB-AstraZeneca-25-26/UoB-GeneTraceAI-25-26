# MEASUREMENTS — the six B3 computations

All numbers below were computed on 2026-08-09 against the repository as committed on branch `chai`. Every result is tagged **VERIFIED** (computed here) unless stated otherwise. Scripts are reproduced in full so each measurement can be re-run independently.

**Environment:** `C:\Python314\python.exe` — Python 3.14.3, pandas 2.3.3, scipy 1.17.1, pyarrow, duckdb 1.5.4. Run from repo root.

**Shared conventions.**
- Gene universe: `reference/gene_lookup.parquet` filtered to `biotype == "protein_coding" & hgnc_status == "Approved"` → **19,213** genes. This is filter A of the three described in `STATE_REPORT.md` B2.5.
- Cell-line identity: all sources mapped to `model_id` through `src/pipeline/outputs/harmonised_enriched.parquet` (`profile_ids` for DepMap, `geo_accessions` for GEO, `hpa_cell_lines` for HPA). Multiple source records per `model_id` are averaged.

---

## B3.1 Source correlation matrix

**Question.** Pairwise Spearman ρ between DepMap, HPA and GEO expression, computed *within gene across shared cell lines*, then summarised across genes. Determines whether an independence-assuming combination is valid.

### Result — 3×3 matrix of median per-gene Spearman ρ

|  | DepMap | HPA | GEO |
|---|---|---|---|
| **DepMap** | 1.000 | **0.817** | 0.467 |
| **HPA** | **0.817** | 1.000 | 0.523 |
| **GEO** | 0.467 | 0.523 | 1.000 |

### Full distributions

| Pair | Shared lines | Genes | Median ρ | IQR | Mean ρ | Frac ρ>0 |
|---|---|---|---|---|---|---|
| DepMap–HPA | 1,040 | 19,171 | **0.817** | 0.700 – 0.903 | 0.765 | 0.995 |
| DepMap–GEO | 523 | 15,962 | 0.467 | 0.266 – 0.622 | 0.439 | 0.946 |
| HPA–GEO | 512 | 15,962 | 0.523 | 0.290 – 0.664 | 0.471 | 0.943 |

### Interpretation

**DepMap and HPA are near-redundant measurements of the same quantity**, not two independent observations. The 95%-of-genes-positive figure with a median of 0.817 leaves very little independent information in the second source. Treating them as independent in a Stouffer or noisy-OR combination will inflate apparent agreement — the same error class as §2.3, and larger in magnitude than the ρ_EP = 0.46 the pipeline already acknowledges for expression-vs-protein.

**GEO is materially less correlated with both** (0.467 / 0.523), consistent with its detection floor and platform heterogeneity. Note the direction: GEO correlates *better* with HPA than with DepMap, which is not what a pure "GEO is noisy" story predicts and is worth a sentence of explanation before GEO is relied on for anything.

Required correction, if DepMap and HPA are ever combined: see `BUILD_SPEC.md` C1.

### Script — `b3_1_2_sources.py`

```python
"""B3.1 source correlation matrix + B3.2 line overlap."""
import os, json, ast, re, numpy as np, pandas as pd
import pyarrow.parquet as pq
os.chdir(r"C:\Disertation\UoB-GeneTraceAI-25-26")

h = pd.read_parquet("src/pipeline/outputs/harmonised_enriched.parquet",
                    columns=["model_id", "profile_ids", "geo_accessions", "hpa_cell_lines"])

def explode_map(col):
    out = {}
    for mid, v in zip(h.model_id, h[col]):
        if v is None: continue
        try: lst = ast.literal_eval(v) if isinstance(v, str) else list(v)
        except Exception: continue
        for x in lst:
            if x: out[str(x).lower()] = mid
    return out

prof2m, gsm2m, hpa2m = (explode_map(c) for c in
                        ["profile_ids", "geo_accessions", "hpa_cell_lines"])

gl = pd.read_parquet("reference/gene_lookup.parquet",
                     columns=["ensg_id", "biotype", "hgnc_status"])
uni = set(gl[(gl.biotype == "protein_coding") &
            (gl.hgnc_status == "Approved")].ensg_id.str.lower())

# ---- DepMap: headers are 'SYMBOL (ENSG…)', rows are profiles
DEP_PATH = "data/parquet/2_DepMap_OmicsExpressionAllGenesTPMLogp1Profile.parquet"
dcols = pq.ParquetFile(DEP_PATH).schema_arrow.names
col2ensg = {c: m.group(1).lower() for c in dcols[1:]
            if (m := re.search(r"(ENSG\d+)", str(c)))}
keep = [dcols[0]] + [c for c, e in col2ensg.items() if e in uni]
dep = pd.read_parquet(DEP_PATH, columns=keep).reset_index()
dep = dep.rename(columns={dep.columns[0]: dcols[0]})
dep[dcols[0]] = dep[dcols[0]].astype(str).str.lower()
dep["model_id"] = dep[dcols[0]].map(prof2m)
dep = dep.dropna(subset=["model_id"]).drop(columns=[dcols[0]])
dep.columns = [col2ensg.get(c, c) for c in dep.columns]
DEP = dep.groupby("model_id").mean()

# ---- HPA: long form, use nTPM
hpa = pd.read_parquet("data/parquet/1_4_hpa_rna_celline.parquet")
gcol = [c for c in hpa.columns if c.lower() == "gene"][0]
ccol = [c for c in hpa.columns if "cell" in c.lower()][0]
lower2real = {c.lower(): c for c in hpa.columns}
vcol = lower2real.get("ntpm") or lower2real["tpm"]
hpa["model_id"] = hpa[ccol].astype(str).str.lower().map(hpa2m)
hpa[gcol] = hpa[gcol].astype(str).str.lower()
hpa = hpa[hpa.model_id.notna() & hpa[gcol].isin(uni)]
HPA = hpa.groupby(["model_id", gcol])[vcol].mean().unstack()
del hpa

# ---- GEO: genes x samples, transpose
GEO_PATH = "data/parquet/3_GEOexpression.parquet"
gcols = pq.ParquetFile(GEO_PATH).schema_arrow.names
gsm_keep = [c for c in gcols[1:] if c.lower() in gsm2m]
geo = pd.read_parquet(GEO_PATH, columns=[gcols[0]] + gsm_keep)
geo[gcols[0]] = geo[gcols[0]].astype(str).str.lower()
geo = geo[geo[gcols[0]].isin(uni)].set_index(gcols[0])
GEO = geo.T
GEO.index = [gsm2m[str(i).lower()] for i in GEO.index]
GEO = GEO.groupby(level=0).mean()
del geo

# ---- B3.2 Venn
D, H, G = set(DEP.index), set(HPA.index), set(GEO.index)

# ---- B3.1 per-gene Spearman, vectorised
def per_gene_spearman(A, B, min_lines=20):
    lines = sorted(set(A.index) & set(B.index))
    genes = sorted(set(A.columns) & set(B.columns))
    a, b = A.loc[lines, genes], B.loc[lines, genes]
    ra, rb = a.rank(axis=0, na_option="keep"), b.rank(axis=0, na_option="keep")
    mask = ra.notna() & rb.notna()
    ra, rb = ra.where(mask), rb.where(mask)
    ra, rb = ra - ra.mean(), rb - rb.mean()
    rho = ((ra * rb).sum() / np.sqrt((ra**2).sum() * (rb**2).sum())
           ).replace([np.inf, -np.inf], np.nan)
    return rho[mask.sum() >= min_lines], len(lines)

for A, B, nm in [(DEP, HPA, "DepMap-HPA"), (DEP, GEO, "DepMap-GEO"), (HPA, GEO, "HPA-GEO")]:
    rho, nl = per_gene_spearman(A, B)
    print(f"{nm}: lines={nl} genes={len(rho):,} median={rho.median():.3f} "
          f"IQR={rho.quantile(.25):.3f}-{rho.quantile(.75):.3f} "
          f"mean={rho.mean():.3f} frac>0={(rho>0).mean():.3f}")
```

**Caveat.** Scales differ across sources (DepMap `log2(TPM+1)`, HPA `nTPM`, GEO mixed). Spearman is rank-based so this does not bias ρ, but it does mean the matrix says nothing about *agreement in level* — only about agreement in ordering. See C1.

---

## B3.2 Line overlap

**Question.** How many cell lines are shared by all three sources, and by each pair?

| Region | Count |
|---|---|
| DepMap total | 1,479 |
| HPA total | 1,104 |
| GEO total | 588 |
| DepMap ∩ HPA | 1,040 |
| DepMap ∩ GEO | 523 |
| HPA ∩ GEO | 512 |
| **All three** | **484** |
| DepMap only | 400 |
| HPA only | 36 |
| GEO only | 37 |
| DepMap ∩ HPA, not GEO | 556 |
| DepMap ∩ GEO, not HPA | 39 |
| HPA ∩ GEO, not DepMap | 28 |
| **Union** | **1,580** |

VERIFIED. Script as B3.1 above.

**Interpretation.** GEO adds only **37 lines** not covered by DepMap or HPA — 2.3% of the union. Its value is corroboration on the 523 lines it shares with DepMap, not coverage extension. Conversely, DepMap alone contributes 400 unique lines. A design that requires all three sources would collapse the panel to 484 lines, a 69% loss against the union.

---

## B3.3 Proteomics coverage bias (the MNAR problem)

**Question.** For lines with proteomics vs without, compare the distribution of mean RNA expression.

**Method.** Mean `log2(TPM+1)` across 19,173 protein-coding genes per cell line (DepMap, n=1,479 lines). Split by proteomics coverage from `coverage_matrix_enriched.parquet`. Effect sizes: Cliff's δ (rank-based, no distributional assumption) and Cohen's d. Test: Mann–Whitney U.

| Platform | n with | n without | mean RNA with | mean RNA without | **Cliff's δ** | Cohen's d | MWU p |
|---|---|---|---|---|---|---|---|
| ProCan DIA | 687 | 792 | 2.6614 | 2.5819 | **+0.2158** | +0.3887 | 7.61e-13 |
| CCLE/Gygi TMT | 369 | 1,110 | 2.6743 | 2.6003 | **+0.2123** | +0.3592 | 9.55e-10 |
| **Either platform** | **769** | **710** | **2.6622** | **2.5718** | **+0.2511** | **+0.4445** | **6.68e-17** |

VERIFIED.

**Interpretation.** Lines with proteomics are drawn preferentially from higher-expressing lines, at a small-to-moderate effect size that is overwhelmingly significant. Cliff's δ = +0.251 means that for a randomly chosen covered/uncovered pair, the covered line has higher mean RNA about 63% of the time. **Both platforms show the same bias at nearly the same magnitude**, so this is a property of which lines get selected for mass spectrometry, not of a particular assay.

This is §9 gap 13 in numeric form. Any protein-derived quantity — including absolute copies per cell under C6 — inherits it. Precision on the covered subset must not be reported as representativeness over the panel.

### Script — `b3_3456.py` (B3.3 section)

```python
import numpy as np, pandas as pd
from scipy.stats import mannwhitneyu
# mean_rna: per-model_id mean of DepMap log2(TPM+1) over the protein-coding universe
# (built exactly as in b3_1_2_sources.py, then .mean(axis=1))
cov = pd.read_parquet("src/pipeline/outputs/coverage_matrix_enriched.parquet").set_index("model_id")
for plat, col in [("ProCan DIA", "procan_proteomics"), ("CCLE/Gygi TMT", "proteomics"),
                  ("either platform", None)]:
    has = (cov["procan_proteomics"] | cov["proteomics"]) if col is None else cov[col]
    j = pd.concat([mean_rna, has.rename("has")], axis=1).dropna()
    a, b = j[j.has].mean_rna, j[~j.has].mean_rna
    u, p = mannwhitneyu(a, b)
    cliff = 2 * u / (len(a) * len(b)) - 1
    d = (a.mean() - b.mean()) / np.sqrt(
        ((len(a)-1)*a.var() + (len(b)-1)*b.var()) / (len(a)+len(b)-2))
    print(f"{plat}: with={len(a)} without={len(b)} {a.mean():.4f} vs {b.mean():.4f} "
          f"Cliff={cliff:+.4f} d={d:+.4f} p={p:.3e}")
```

---

## B3.4 Ploidy distribution

**Question.** Does ploidy exist, and what is its distribution? Bounds the accuracy of any proteomic-ruler estimate.

**Source.** `src/pipeline/outputs/gdsc_models.parquet`, columns `ploidy_wes` / `ploidy_wgs`, originating from `data/GDSC/model_list_20260709.csv`. **PROVEN present** (A1.3).

| Estimate | n | min | Q1 | median | Q3 | max | frac > 2.5 | frac near-diploid (1.9–2.1) |
|---|---|---|---|---|---|---|---|---|
| `ploidy_wes` | 1,248 | 1.64 | 2.050 | **2.671** | 3.268 | 5.40 | **55.8%** | 26.6% |
| `ploidy_wgs` | 254 | 1.52 | 2.000 | 2.590 | 2.995 | 5.00 | 52.8% | 24.0% |
| either (WGS preferred) | **1,502** | — | — | **2.653** | — | — | — | — |

Coverage: 1,502 of 2,145 panel models = **70.0%**. VERIFIED.

**Interpretation.** The panel is heavily aneuploid: the median line has ploidy **2.65**, and **only about a quarter of lines are near-diploid**. Two consequences:

1. **For C6.** A proteomic ruler applied with an assumed ploidy of 2.0 would misestimate DNA per cell by a median factor of 2.65/2.0 = **1.33×**, and by up to 2.7× at the extreme. Ploidy correction is mandatory, exactly as the brief states — and it is available for 70% of the panel, which also means **30% of lines cannot be ploidy-corrected at all** and must be reported as not-estimable rather than silently assumed diploid.
2. **For gap 8.** `AMP_THRESHOLD = 2.5` and `DEL_THRESHOLD = 1.5` on absolute copy number are applied to a panel whose median ploidy is 2.65. On a median line, the amplification threshold sits *below* the genome-wide average copy number. The same correction closes both items.

---

## B3.5 Lineage size distribution

**Question.** Lines per lineage; how many lineages have ≥15 lines, and what fraction of the panel that covers.

**Source.** `gdsc_models.parquet:tissue` (30 distinct labels). A second, unreconciled vocabulary exists at `sample_info.lineage` — see `STATE_REPORT.md` A1.5.

| Tissue | n | | Tissue | n |
|---|---|---|---|---|
| haematopoietic and lymphoid | 382 | | bladder | 41 |
| lung | 317 | | biliary tract | 39 |
| large intestine | 223 | | cervix | 33 |
| skin | 135 | | endometrium | 32 |
| central nervous system | 119 | | liver | 30 |
| esophagus | 119 | | thyroid | 25 |
| ovary | 106 | | prostate | 18 |
| head and neck | 103 | | eye | 13 |
| pancreas | 92 | | uterus | 10 |
| bone | 82 | | unknown | 7 |
| breast | 82 | | placenta | 5 |
| soft tissue | 73 | | testis | 5 |
| kidney | 64 | | small intestine | 3 |
| stomach | 53 | | vulva | 3 |
| peripheral nervous system | 51 | | adrenal gland | 1 |

### Coverage at `MIN_PEERS = 15`

| Denominator | Labelled | Unlabelled | Lineages ≥15 | Models in those lineages |
|---|---|---|---|---|
| Coverage panel (2,145 models) | 1,781 (83.0%) | **364 (17.0%)** | 21 of 30 | 1,725 = **80.4%** of panel |
| Scored in `core_score` (1,746 models) | 1,634 (93.6%) | 112 (6.4%) | 21 of 30 | 1,587 = **90.9%** |

VERIFIED. (This is close to `test_run_improvement_scan.py:17-18`'s note of "21 of 30 tissues have ≥15 lines (1,726 of 1,781 models)" — that figure uses labelled models as denominator; the table above uses the panel and the scored set, which is the honest denominator for a coverage claim.)

**Interpretation.** Lineage statistics are stable for **90.9% of scored lines**. The cost of a ≥15 rule is modest. But **17% of panel models carry no tissue label at all**, and those must be routed to a pooled fallback, not dropped and not assigned a default lineage. Nine lineages (eye, uterus, unknown, placenta, testis, small intestine, vulva, adrenal gland, and prostate at n=18 only marginally above) are too small for stable within-lineage statistics.

---

## B3.6 Alteration prevalence

**Question.** For mutation, fusion and CNA separately: what fraction of (gene, line) pairs carry an alteration? Determines whether requirement C4 demotes a handful of pairs or a large fraction.

**Denominator.** All 29,781,274 scored (gene, line) pairs in `core_score.parquet` — 19,177 genes × 1,746 lines. Alteration flags left-joined from `flags_with_driver.parquet`; unmatched treated as *no flag present* (see the caveat below, which is the point).

| Alteration | Pairs | % of scored pairs |
|---|---|---|
| Mutation (`p_mutation > 0.5`) | 349,044 | **1.172%** |
| Fusion (`p_fusion > 0.5`) | 129,046 | 0.433% |
| **Any alteration** | **475,695** | **1.597%** |
| Driver alteration | 143,523 | 0.482% |
| CNA alteration | 2,717 | **0.009%** |

### Distribution across genes

| Statistic | Value |
|---|---|
| Per-gene fraction of lines altered — median | 0.0116 |
| — IQR | 0.0061 – 0.0202 |
| — max | 0.6229 |
| Genes with >10% of lines altered | **87** of 19,177 |
| Genes with **0%** of lines altered | **596** of 19,177 |

VERIFIED.

### Interpretation — the direct answer to C4's sizing question

**A global "alterations demote" rule would touch 1.6% of scored pairs.** It demotes a handful, not a large fraction. The distribution is extremely skewed: the median gene has ~1.2% of its lines altered, while 87 genes exceed 10% and one reaches 62%. So a global rule is nearly a no-op for most genes and a substantial re-ordering for a few — which argues against a global continuous demotion and in favour of a gate applied only where the alteration is actually prevalent enough to matter.

### The CNA caveat — this is a C3 finding, not a C4 one

**CNA context is available for only 91,879 of 29,781,274 pairs (0.31%).** The 0.009% "CNA alteration" rate is therefore not a biological prevalence — it is 2,717 alterations found within the 0.31% of pairs that were *looked at*. For the other 99.69% of pairs, `has_cna_alteration` is `False` because `build_full_predictions.py:122` and `build_evidence_ledger.py:68` apply `.fillna(False)`, collapsing **not-measured into measured-negative**.

Any statement of the form "this line has no CNA in this gene" is, for 99.69% of pairs, a statement that no CNA data was consulted. This is the single clearest instance of the three-state failure C3 is about, and it is listed first in `BUILD_SPEC.md` C3.

### Script — `b3_3456.py` (B3.6 section)

```python
import pandas as pd
OUT = "src/pipeline/outputs"
fl = pd.read_parquet(f"{OUT}/flags_with_driver.parquet",
                     columns=["model_id","ensg_id","p_mutation","p_fusion","has_alteration",
                              "has_driver_alteration","has_cna_alteration","has_cna_context"])
core = pd.read_parquet(f"{OUT}/core_score.parquet", columns=["model_id","ensg_id"])
core["model_id"] = core["model_id"].str.lower()
fl["model_id"] = fl["model_id"].str.lower()
key = core.merge(fl, on=["model_id","ensg_id"], how="left")
n = len(key)
for label, mask in [("mutation (p>0.5)", key.p_mutation > 0.5),
                    ("fusion (p>0.5)",   key.p_fusion > 0.5),
                    ("any alteration",   key.has_alteration.fillna(False)),
                    ("driver alteration", key.has_driver_alteration.fillna(False)),
                    ("CNA alteration",   key.has_cna_alteration.fillna(False))]:
    print(f"{label:22s} {int(mask.sum()):>10,}  {100*mask.mean():6.3f}%")
pg = key.assign(alt=key.has_alteration.fillna(False)).groupby("ensg_id")["alt"].mean()
print(f"per-gene: median={pg.median():.4f} IQR={pg.quantile(.25):.4f}-{pg.quantile(.75):.4f} "
      f"max={pg.max():.4f}; >10%: {int((pg>0.10).sum()):,}; ==0: {int((pg==0).sum()):,}")
print(f"CNA context pairs: {int(key.has_cna_context.fillna(False).sum()):,} "
      f"({100*key.has_cna_context.fillna(False).mean():.2f}%)")
```

---

## Supporting measurements (not in B3, used elsewhere in the report)

### Chronos sign diagnostic — `b11_chronos_sign.py` / `b11b_followup.py`

| Quantity | Value |
|---|---|
| median `essentiality`, DepMap ReferenceEssentials (n=1,043) | −5.198 |
| median `essentiality`, ReferenceNonEssentials (n=558) | +3.492 |
| Mann–Whitney p | 3.29e-229 |
| median `essentiality`, RPL*/RPS*/POLR2A/PSMA*/EIF* (n=147) | −10.682 |
| mean `chronos_pct` — essentials / non-essentials | 0.4995 / 0.4995 |
| control-set median ρ (1,149 genes) | −0.0655 (79.3% negative) |
| Spearman(per-line screen strength, per-line mean `core_score`) | −0.019, p=0.58, n=891 |
| shipped `chronos_validation` significant genes | 4,332 → 2,468 neg / 1,864 pos = 57.0%, binomial p = 4.46e-20 |

### Metric recomputation — `b13_b14_metrics.py`

Full table in `STATE_REPORT.md` B1.4. Chance baselines used:
- `E[recall@k] = k/N`
- `P[any-hit@k] = 1 − C(N−K, k)/C(N, k)` (hypergeometric), computed per gene via `gammaln` and averaged.
- Paired bootstrap: 10,000 resamples over genes of `mean(observed − chance)`, seed 42.
- Median N (scored lines per gene) = 1,563; median K (sensitive lines per gene) = 228.

### Proteomics characterisation — `a12_proteomics.py`, `a12b_histones.py`, `a12c_ruler.py`

| Quantity | ProCan DIA | CCLE/Gygi TMT |
|---|---|---|
| value median / % negative | 3.44 / 1.40% | −0.042 / 52.6% |
| per-line median: sd, range | 0.212, [2.667, 4.811] | 0.088, [−0.636, 0.109] |
| per-line median spread p1–p99 | 0.98 log2 units (≈2× linear) | — |
| proteins per line (median) | 5,218 | 9,072 |
| histone proteins detected | 17 | 23 |
| histone/total linear-intensity ratio | median 0.0606, IQR 0.0454–0.0816, **CV 0.411** | n/a (ratios) |
| Spearman(histone ratio, depth) | −0.066, p=0.042 | — |

Histone ratio by depth quintile (ProCan): 0.0687 / 0.0604 / 0.0561 / 0.0579 / 0.0625 at median depths 4,665 / 5,032 / 5,220 / 5,421 / 5,677. No monotone depth trend, but a **41% coefficient of variation** on the ruler's internal standard. See `BUILD_SPEC.md` C6.
