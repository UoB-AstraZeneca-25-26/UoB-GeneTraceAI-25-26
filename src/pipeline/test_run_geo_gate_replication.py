"""
test_run_geo_gate_replication.py
--------------------------------
Read-only diagnostic. Writes no production file.

QUESTION
--------
The gate is the project's central finding: a gene that a cell line does not
express behaves as non-essential there, so the bottom of the abundance ranking
can be abstained on. It was established in test_run_abundance_gate.py and it
rests entirely on ONE expression source -- DepMap RNA-seq. The gate, the floor
result and the scope condition are all downstream of that single matrix.

test_run_screen_replication.py already removed the equivalent single-source risk
on the LABEL side: it re-ran the gate against Sanger's independent CRISPR screen
and it held. This script does the same thing on the EXPRESSION side, which is
the side that was never checked.

    Recompute the gate using a GEO series as the expression source, against the
    same Chronos labels, restricted to the cell lines that series covers.

  reproduces  -> the gate is a property of GENE EXPRESSION, not an artefact of
                 DepMap's alignment/normalisation pipeline. The strongest
                 available statement, because GEO is a different assay
                 (Affymetrix GPL570 arrays) processed by different people.
  fails       -> a platform dependence exists and must be DISCLOSED as a scope
                 condition: the gate is calibrated on RNA-seq and does not
                 transfer to array expression.

Both outcomes are publishable. Neither is a nuisance result.

WHY AN ARRAY IS A HARD TEST, AND WHY THAT IS THE POINT
------------------------------------------------------
GPL570 measures hybridisation intensity, not transcript counts. It has a
compressed dynamic range, a non-zero background, and no true zero -- an unex-
pressed gene reads as background noise, not as 0 TPM. That is precisely the
region the gate operates in, so if the gate survives being asked to find a floor
that the assay does not physically have, it is not an artefact of RNA-seq's zeros.

CONFOUND THIS SCRIPT EXISTS TO KILL
-----------------------------------
A GEO series covers far fewer lines than DepMap. A weaker gate on GEO could
therefore mean "arrays disagree" OR simply "fewer lines, noisier estimate".
These are not the same finding and the naive version of this test cannot tell
them apart. So the gate is computed in THREE arms:

  A  GEO expression        on the series' lines          the replication
  B  DepMap expression     on THE SAME LINES, same genes the matched control
  C  DepMap expression     on the full DepMap cohort     the power reference

A vs B is the platform comparison at identical statistical power -- same lines,
same genes, same labels, same code path, so the ONLY thing that differs is which
instrument measured the RNA. B vs C is the cost of shrinking to 200-odd lines.
Quoting A against C, which is what a single-arm test would do, silently charges
the platform for the loss of power.

RANK INVARIANCE
---------------
Every score is van-der-Waerden transformed within gene across lines before use,
so any monotone difference in units between platforms -- log-TPM against MAS5
intensity, different normalisation constants -- cancels exactly. What survives
is disagreement about the ORDER of cell lines, which is the only thing the gate
depends on. This is why no cross-platform rescaling is attempted or needed.

THE FLOOR, ON AN ASSAY WITH NO FLOOR
------------------------------------
test_run_abundance_gate.py defines the floor as the tie-block at the gene's
minimum -- the lines with literally zero reads. Arrays produce no such block
(measured here: essentially every value distinct), so that definition returns
nothing and cannot be compared across platforms. This script therefore uses a
platform-neutral ABSENT CALL instead, applied identically to both arms: within
each cell line, a gene is called absent if its value falls at or below that
line's own ABSENT_Q-th percentile across all genes. That is the array field's
standard notion of "background", it needs no zeros, and because the threshold is
per line it is invariant to between-line scaling. DepMap's true-zero floor is
reported alongside for arm B/C as the link back to the original result.

READING THE RESULT
------------------
  A and B both show decile-1 depletion, similar magnitude
        -> replicates. The gate is about expression, not about DepMap.
  B depleted, A flat
        -> platform dependence. Disclose; the gate is an RNA-seq calibration.
  Neither A nor B depleted, but C is
        -> not a platform result at all, just a power result. Say so, and stop
           drawing a platform conclusion from 200 lines.

Outputs (suffixed by series):
    src/pipeline/outputs/test_run_geo_gate_replication_<gse>_results.json
    src/pipeline/outputs/test_run_geo_gate_replication_<gse>_per_gene.parquet

Run:
    python src/pipeline/test_run_geo_gate_replication.py
    python src/pipeline/test_run_geo_gate_replication.py --gse gse34211
    python src/pipeline/test_run_geo_gate_replication.py --gse all
"""
import argparse
import json
import warnings
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
import pyarrow.parquet as pq
from scipy.stats import hypergeom, norm, rankdata, spearmanr, wilcoxon

warnings.filterwarnings("ignore", category=RuntimeWarning)

OUTPUTS = Path("src/pipeline/outputs")
REF_DIR = Path("reference")
DATA_CLEAN = Path("data/parquet/data_clean")
CHR_DIR = Path("data/DepMap_Chronos")
BROAD = CHR_DIR / "GeneFitnessEffect_Chronos_Achilles.hdf5"

DEP_THRESHOLD = -0.5     # Chronos gene effect at or below this = dependent
MIN_LINES = 100          # deciles need enough lines to mean anything
MIN_POS = 8
MIN_NEG = 8
ABSENT_Q = 10            # per-line percentile defining the array background
SEED = 42

ap = argparse.ArgumentParser()
ap.add_argument("--gse", default="gse57083",
                help="GEO series to use as the expression source, or 'all' to "
                     "pool every series with per-series batch removal")
ap.add_argument("--n-genes", type=int, default=0,
                help="0 = every eligible gene; otherwise a random sample")
ap.add_argument("--min-lines", type=int, default=MIN_LINES)
args = ap.parse_args()
GSE = args.gse.lower()
MIN_LINES = args.min_lines

results = {"series": GSE, "dep_threshold": DEP_THRESHOLD,
           "absent_call_percentile": ABSENT_Q}


# ============================================================ helpers
def vdw(s: pd.Series) -> pd.Series:
    """Van der Waerden: rank -> normal quantile. Monotone, so unit-free."""
    n = int(s.notna().sum())
    return pd.Series(norm.ppf(s.rank(method="min") / (n + 1)), index=s.index)


def partial_auc(score, pos, lo, hi):
    """McClish-standardised partial AUROC over FPR in [lo, hi]. 0.5 = chance."""
    ok = ~np.isnan(score)
    score, pos = score[ok], pos[ok]
    if pos.sum() < 2 or (~pos).sum() < 2:
        return np.nan
    order = np.argsort(-score, kind="mergesort")
    s, p = score[order], pos[order]
    distinct = np.r_[np.diff(s) != 0, True]
    fpr = np.r_[0, np.cumsum(~p)[distinct] / (~p).sum()]
    tpr = np.r_[0, np.cumsum(p)[distinct] / p.sum()]
    grid = np.linspace(lo, hi, 501)
    area = np.trapezoid(np.interp(grid, fpr, tpr), grid)
    a_min, a_max = (hi ** 2 - lo ** 2) / 2.0, hi - lo
    return 0.5 * (1.0 + (area - a_min) / (a_max - a_min))


def auroc(score, pos):
    ok = ~np.isnan(score)
    score, pos = score[ok], pos[ok]
    n1, n0 = int(pos.sum()), int((~pos).sum())
    if n1 < 1 or n0 < 1:
        return np.nan
    r = rankdata(score)
    return (r[pos].sum() - n1 * (n1 + 1) / 2.0) / (n1 * n0)


def boot_ci(v, n_boot=2000, stat=np.median):
    v = np.asarray(pd.Series(v).dropna(), dtype=float)
    if len(v) < 5:
        return (np.nan, np.nan)
    rg = np.random.default_rng(SEED)
    b = [stat(rg.choice(v, len(v), replace=True)) for _ in range(n_boot)]
    return (float(np.percentile(b, 2.5)), float(np.percentile(b, 97.5)))


def gate_metrics(score, pos, absent, prefix, rec):
    """Every gate-shaped readout for one (score, label) pair, written into rec.

    score   van-der-Waerden score, high = more abundant
    pos     boolean dependency label
    absent  boolean absent-call mask, same length
    """
    ok = ~np.isnan(score)
    s, p, ab = score[ok], pos[ok], absent[ok]
    n, npos = len(s), int(p.sum())
    if n < MIN_LINES or npos < MIN_POS or (n - npos) < MIN_NEG:
        return False
    base = npos / n
    rec[f"n_{prefix}"] = n
    rec[f"base_{prefix}"] = float(base)

    # deciles by ordinal rank so tie blocks distribute deterministically
    r = rankdata(s, method="ordinal")
    q = np.minimum((r - 1) * 10 // n, 9)
    rates = np.array([p[q == i].mean() if (q == i).sum() else np.nan
                      for i in range(10)])
    rec[f"profile_{prefix}"] = (rates / base).tolist()
    rec[f"d1_{prefix}"] = float(rates[0] / base)
    rec[f"d10_{prefix}"] = float(rates[9] / base)
    # exact null for the decile-1 count: draw n1 of n lines without replacement
    # holding prevalence fixed -> hypergeometric, one-sided for depletion
    n1 = int((q == 0).sum())
    rec[f"d1_p_{prefix}"] = float(hypergeom.cdf(int(p[q == 0].sum()), n, npos, n1))

    rec[f"pauc_gate_{prefix}"] = partial_auc(s, p, 0.8, 1.0)
    rec[f"pauc_rank_{prefix}"] = partial_auc(s, p, 0.0, 0.2)
    rec[f"auc_{prefix}"] = auroc(s, p)

    for pctl in (5, 10, 20):
        k = max(1, int(round(n * pctl / 100)))
        below = np.argsort(s, kind="mergesort")[:k]
        rec[f"npv{pctl}_{prefix}"] = float(1.0 - p[below].mean())
    rec[f"npv_base_{prefix}"] = float(1.0 - base)

    # absent call floor -- the platform-neutral version of the floor test
    if ab.sum() >= 15:
        rec[f"n_absent_{prefix}"] = int(ab.sum())
        rec[f"absent_rate_{prefix}"] = float(p[ab].mean())
        rec[f"absent_ratio_{prefix}"] = float(p[ab].mean() / base)
    return True


# ============================================================ Step 1
print("=" * 78)
print("STEP 1 -- Chronos labels (Broad), with the scale check")
print("=" * 78)

gl = pd.read_parquet(REF_DIR / "gene_lookup.parquet",
                     columns=["ensg_id", "hgnc_symbol", "biotype", "hgnc_status"])
gl["e"] = gl.ensg_id.astype(str).str.split(".").str[0].str.lower()
s2e = dict(zip(gl.hgnc_symbol.astype(str).str.upper(), gl.e))
sym_of = {v: k for k, v in s2e.items()}
pc_genes = set(gl[(gl.biotype == "protein_coding") &
                  (gl.hgnc_status == "Approved")].e)

with h5py.File(BROAD, "r") as fh:
    D = fh["data"][:]
    ch_lines = [x.decode().lower() for x in fh["dim_0"][:]]
    ch_genes = [x.decode() for x in fh["dim_1"][:]]
ki, ke = [], []
for i, c in enumerate(ch_genes):
    e = s2e.get(c.split(" (")[0].strip().upper())
    if e:
        ki.append(i)
        ke.append(e)
CH = pd.DataFrame(D[:, ki], index=ch_lines, columns=ke)
CH = CH.T.groupby(level=0).mean().T
del D

ess = set(pd.read_csv(CHR_DIR / "ReferenceEssentials.csv").iloc[:, 0]
          .astype(str).str.split(" (", regex=False).str[0].str.upper())
non = set(pd.read_csv(CHR_DIR / "ReferenceNonEssentials.csv").iloc[:, 0]
          .astype(str).str.split(" (", regex=False).str[0].str.upper())
e_ess = {s2e[g] for g in ess if g in s2e} & set(CH.columns)
e_non = {s2e[g] for g in non if g in s2e} & set(CH.columns)
m_ess = float(np.nanmedian(CH[list(e_ess)].values))
m_non = float(np.nanmedian(CH[list(e_non)].values))
LABEL_FPR = float(np.nanmean(CH[list(e_non)].values <= DEP_THRESHOLD))
print(f"  Chronos Achilles: {CH.shape[0]} lines x {CH.shape[1]:,} genes")
print(f"  scale check: median {np.nanmedian(CH.values):+.4f} (expect ~0);  "
      f"essentials {m_ess:+.4f} (~-1);  non-essentials {m_non:+.4f} (~0)")
if not (m_ess < -0.5 < m_non + 0.5):
    raise SystemExit("Chronos scale check FAILED. Aborting.")
print(f"  LABEL FPR at {DEP_THRESHOLD}: {LABEL_FPR:.4f}  "
      f"({len(e_non)} curated non-essential genes)")
results["label_fpr"] = LABEL_FPR

# ============================================================ Step 2
print()
print("=" * 78)
print(f"STEP 2 -- GEO expression source ({GSE})")
print("=" * 78)

info = pd.read_parquet(DATA_CLEAN / "geo_info_clean.parquet")
info["geo_accession"] = info.geo_accession.astype(str).str.lower()
info["cvcl"] = info.cellosaurus_id.astype(str).str.lower()
info["gse_id"] = info.gse_id.astype("string").str.lower()
si = pd.read_parquet(DATA_CLEAN / "sample_info_clean.parquet",
                     columns=["depmap_id", "rrid"])
cvcl2model = dict(zip(si.rrid.astype(str).str.lower(),
                      si.depmap_id.astype(str).str.lower()))
info["model_id"] = info.cvcl.map(cvcl2model)

geo_cols = set(pq.ParquetFile(DATA_CLEAN / "geo_expr_clean.parquet")
               .schema_arrow.names) - {"gene"}
sel = info[info.model_id.notna() & info.geo_accession.isin(geo_cols)].copy()
if GSE != "all":
    sel = sel[sel.gse_id == GSE]
    if not len(sel):
        raise SystemExit(f"series {GSE} has no usable samples")
sel["series"] = sel.gse_id.fillna("(no_gse)")
print(f"  usable GSM samples : {len(sel):,}  "
      f"across {sel.series.nunique()} series")
print(f"  distinct cell lines: {sel.model_id.nunique():,}")

# ------------------------------------------------------------ gene set
geo_gene_col = pq.read_table(DATA_CLEAN / "geo_expr_clean.parquet",
                             columns=["gene"]).to_pandas()["gene"].astype(str)
expr_names = pq.ParquetFile(DATA_CLEAN / "depmap_expr_clean.parquet").schema_arrow.names
dm_col = {c.split(".")[0].lower(): c for c in expr_names if c != "index"}
genes = sorted(set(geo_gene_col) & set(dm_col) & set(CH.columns) & pc_genes)
if args.n_genes and args.n_genes < len(genes):
    genes = sorted(np.random.default_rng(SEED).choice(genes, args.n_genes,
                                                      replace=False))
print(f"  genes in GEO n DepMap n Chronos n protein-coding : {len(genes):,}")

# ------------------------------------------------------------ GEO matrix
gsm_cols = sorted(sel.geo_accession.unique())
GX = pq.read_table(DATA_CLEAN / "geo_expr_clean.parquet",
                   columns=["gene"] + gsm_cols).to_pandas()
GX = GX.set_index("gene").reindex(genes)
# GPL570 values are linear intensities. log2 makes replicate averaging a
# geometric mean, which is the right summary for intensity data; every score is
# rank-transformed later so this affects nothing downstream except that average.
GX = np.log2(GX.clip(lower=1e-3))
print(f"  GEO matrix: {GX.shape[0]:,} genes x {GX.shape[1]:,} samples")

gsm2model = dict(zip(sel.geo_accession, sel.model_id))
gsm2series = dict(zip(sel.geo_accession, sel.series))

if GSE == "all" and sel.series.nunique() > 1:
    # Pooling series would let a between-series batch shift masquerade as
    # between-line biology. Standardise WITHIN series, per gene, before pooling:
    # each series contributes only its own ordering of its own lines.
    print("  pooling mode: standardising within series per gene before pooling")
    parts = []
    for s_name, sub in sel.groupby("series"):
        cols = sorted(set(sub.geo_accession))
        if len(cols) < 5:
            continue
        blk = GX[cols]
        by_line = blk.T.groupby(pd.Series(cols).map(gsm2model).values).mean().T
        if by_line.shape[1] < 5:
            continue
        # per gene, rank lines within this series, mapped to [0,1]
        rk = by_line.rank(axis=1, pct=True, method="average")
        parts.append(rk)
    G_line = pd.concat(parts, axis=1)
    G_line = G_line.T.groupby(level=0).mean().T          # line seen in >1 series
else:
    G_line = GX.T.groupby(pd.Series(GX.columns).map(gsm2model).values).mean().T

G_line = G_line.loc[:, ~G_line.columns.duplicated()]
print(f"  GEO collapsed to cell lines: {G_line.shape[1]:,}")

# ============================================================ Step 3
print()
print("=" * 78)
print("STEP 3 -- DepMap expression, and the three line sets")
print("=" * 78)

prof = pd.read_parquet(REF_DIR / "depmap_profiles.parquet")
rp = prof[prof.datatype == "rna"][["profileid", "modelid"]].copy()
rp["model_id"] = rp.modelid.astype(str).str.lower()
rp["profileid"] = rp.profileid.astype(str)
E = pq.read_table(DATA_CLEAN / "depmap_expr_clean.parquet",
                  columns=["index"] + [dm_col[g] for g in genes]).to_pandas()
if "index" not in E.columns:
    E = E.reset_index()
E["profileid"] = E["index"].astype(str)
E = E.drop(columns=["index"]).merge(rp[["profileid", "model_id"]],
                                    on="profileid", how="inner")
E = E.drop(columns=["profileid"]).groupby("model_id").mean()
E.columns = [c.split(".")[0].lower() for c in E.columns]
E = E[genes]
print(f"  DepMap expression: {E.shape[0]:,} lines x {E.shape[1]:,} genes")

L_match = sorted(set(G_line.columns) & set(E.index) & set(CH.index))
L_full = sorted(set(E.index) & set(CH.index))
print(f"  arm A/B matched line set (GEO n DepMap n Chronos) : {len(L_match):,}")
print(f"  arm C   full DepMap cohort (DepMap n Chronos)     : {len(L_full):,}")
results["lines"] = {"matched": len(L_match), "full_depmap": len(L_full),
                    "geo_lines_total": int(G_line.shape[1])}
if len(L_match) < MIN_LINES:
    raise SystemExit(f"only {len(L_match)} matched lines; need {MIN_LINES}")

GM = G_line[L_match]                     # genes x lines
EM = E.loc[L_match].T                    # genes x lines
EF = E.loc[L_full].T
CM = CH.loc[L_match, genes].T
CF = CH.loc[L_full, genes].T

# ------------------------------------------------------------ absent calls
# Per LINE, the ABSENT_Q-th percentile across all genes is that line's
# background. A gene is absent in a line if it sits at or below it. Per-line, so
# it is immune to between-line scaling; identical rule on both platforms.
def absent_mask(M):
    thr = np.nanpercentile(M.values, ABSENT_Q, axis=0)
    return pd.DataFrame(M.values <= thr[None, :], index=M.index, columns=M.columns)


AB_G, AB_M, AB_F = absent_mask(GM), absent_mask(EM), absent_mask(EF)
if GSE == "all":
    print("  NOTE: in pooling mode GEO values are within-series percentile ranks,")
    print("  not intensities, so the arm-A absent call is a rank-space background")
    print("  and is NOT the same construct as arm B's. Read the absent-call rows")
    print("  from a single-series run; the decile and pAUC rows are unaffected,")
    print("  since those are rank-based on both arms either way.")
zero_floor = (EM.values == np.nanmin(EM.values, axis=1, keepdims=True))
print(f"  absent calls at the {ABSENT_Q}th per-line percentile; "
      f"GEO {AB_G.values.mean():.3f} of entries, DepMap {AB_M.values.mean():.3f}")
tie_geo = float(np.mean(GM.values == np.nanmin(GM.values, axis=1, keepdims=True)))
print(f"  true tie-at-minimum blocks: DepMap {zero_floor.mean():.4f} of entries, "
      f"GEO {tie_geo:.4f}")
print(f"    -> the array has no floor to speak of, which is why the absent call")
print(f"       above is the only floor definition comparable across the two.")
results["floor_definition"] = {
    "depmap_tie_at_min_frac": float(zero_floor.mean()),
    "geo_tie_at_min_frac": tie_geo,
    "absent_frac_geo": float(AB_G.values.mean()),
    "absent_frac_depmap": float(AB_M.values.mean())}

# ============================================================ Step 4
print()
print("=" * 78)
print("STEP 4 -- per-gene gate, three arms")
print("=" * 78)
print("  A = GEO on matched lines | B = DepMap on the SAME lines | "
      "C = DepMap, full cohort")

rows = []
gv, ev = GM.values, EM.values
fv, cmv, cfv = EF.values, CM.values, CF.values
abg, abm, abf = AB_G.values, AB_M.values, AB_F.values
for i, g in enumerate(genes):
    rec = {"ensg_id": g, "symbol": sym_of.get(g, "?")}

    # cross-platform agreement on the ordering of the matched lines: the
    # scope condition for everything below it
    a, b = gv[i], ev[i]
    ok = ~np.isnan(a) & ~np.isnan(b)
    if ok.sum() >= 30 and np.std(a[ok]) > 0 and np.std(b[ok]) > 0:
        rec["rho_platform"] = float(spearmanr(a[ok], b[ok]).statistic)

    lab_m = cmv[i]
    pos_m = lab_m <= DEP_THRESHOLD
    okA = gate_metrics(vdw(pd.Series(a)).values, pos_m, abg[i], "A", rec)
    okB = gate_metrics(vdw(pd.Series(b)).values, pos_m, abm[i], "B", rec)

    lab_f = cfv[i]
    okC = gate_metrics(vdw(pd.Series(fv[i])).values, lab_f <= DEP_THRESHOLD,
                       abf[i], "C", rec)
    if okA or okB or okC:
        rows.append(rec)

G = pd.DataFrame(rows)
out_pg = OUTPUTS / f"test_run_geo_gate_replication_{GSE}_per_gene.parquet"
G.to_parquet(out_pg)
both = G.dropna(subset=["d1_A", "d1_B"])
print(f"  genes scored: A {G.d1_A.notna().sum():,}  B {G.d1_B.notna().sum():,}  "
      f"C {G.d1_C.notna().sum():,}  |  A and B both: {len(both):,}")

# ============================================================ Step 5
print()
print("=" * 78)
print("STEP 5 -- DO THE TWO PLATFORMS EVEN AGREE ABOUT EXPRESSION?")
print("=" * 78)
print("  Per-gene Spearman(GEO, DepMap) over the matched lines. This bounds")
print("  everything else: a gate cannot transfer through a measurement the two")
print("  platforms disagree about.")
RP = G.rho_platform.dropna()
qs = RP.quantile([.05, .25, .5, .75, .95])
print(f"    n genes {len(RP):,}   median {RP.median():+.4f}   mean {RP.mean():+.4f}")
print(f"    p05 {qs[.05]:+.3f}  p25 {qs[.25]:+.3f}  p50 {qs[.5]:+.3f}  "
      f"p75 {qs[.75]:+.3f}  p95 {qs[.95]:+.3f}")
print(f"    rho > 0.5 : {100*(RP > 0.5).mean():5.1f}%     "
      f"rho > 0.2 : {100*(RP > 0.2).mean():5.1f}%     "
      f"rho < 0 : {100*(RP < 0).mean():5.1f}%")
results["platform_agreement"] = {
    "n_genes": int(len(RP)), "median_rho": float(RP.median()),
    "mean_rho": float(RP.mean()),
    "quantiles": {str(k): float(v) for k, v in qs.items()},
    "frac_above_0.5": float((RP > 0.5).mean()),
    "frac_above_0.2": float((RP > 0.2).mean()),
    "frac_negative": float((RP < 0).mean())}

# ============================================================ Step 6
print()
print("=" * 78)
print("STEP 6 -- THE DECILE PROFILE, side by side")
print("=" * 78)
print("  positive rate / base rate, median across genes. decile 1 = lowest")
print("  abundance. 1.00 = exactly the base rate; below 1 = DEPLETED.")
print()


def med_profile(col):
    v = [p for p in G[col].dropna() if isinstance(p, (list, np.ndarray))]
    return np.nanmedian(np.vstack(v), axis=0) if v else np.full(10, np.nan)


pA, pB, pC = (med_profile(f"profile_{t}") for t in "ABC")
print(f"  {'decile':>7} {'A  GEO':>10} {'B  DepMap':>11} {'C  DepMap':>11}   "
      f"{'(A)':<24}")
print(f"  {'':>7} {'matched':>10} {'matched':>11} {'full':>11}")
for i in range(10):
    bar = "#" * max(0, int(round(pA[i] * 20)))
    print(f"  {i+1:>7} {pA[i]:>10.3f} {pB[i]:>11.3f} {pC[i]:>11.3f}   {bar}")
for t, p in zip("ABC", (pA, pB, pC)):
    results.setdefault("decile_profiles", {})[t] = [float(x) for x in p]

print()
print(f"  {'':<40} {'A GEO':>10} {'B DepMap':>11} {'C DepMap':>11}")
print(f"  {'':<40} {'matched':>10} {'matched':>11} {'full':>11}")
summary = {}
for name, col, ref in [("decile-1 ratio (median)", "d1", 1.0),
                       ("decile-10 ratio", "d10", 1.0),
                       ("gate-region pAUC (FPR .8-1)", "pauc_gate", 0.5),
                       ("ranker-region pAUC (FPR 0-.2)", "pauc_rank", 0.5),
                       ("full AUROC", "auc", 0.5),
                       ("absent-call ratio (floor)", "absent_ratio", 1.0)]:
    vals, ps = [], []
    for t in "ABC":
        c = f"{col}_{t}"
        v = G[c].dropna() if c in G else pd.Series(dtype=float)
        vals.append(v.median() if len(v) else np.nan)
        ps.append(wilcoxon(v - ref).pvalue if len(v) > 10 else np.nan)
    print(f"  {name:<40} {vals[0]:>10.4f} {vals[1]:>11.4f} {vals[2]:>11.4f}")
    print(f"  {'  vs ' + str(ref) + ', Wilcoxon p':<40} "
          f"{ps[0]:>10.2e} {ps[1]:>11.2e} {ps[2]:>11.2e}")
    summary[col] = {"A_geo": float(vals[0]), "B_depmap_matched": float(vals[1]),
                    "C_depmap_full": float(vals[2]),
                    "p_A": float(ps[0]), "p_B": float(ps[1]), "p_C": float(ps[2])}
results["summary"] = summary

# ============================================================ Step 7
print()
print("=" * 78)
print("STEP 7 -- A vs B PAIRED, the platform comparison at equal power")
print("=" * 78)
print("  Same genes, same lines, same labels. Any difference here is the")
print("  INSTRUMENT and nothing else.")
print()
print(f"  {'metric':<32} {'A GEO':>9} {'B DepMap':>10} {'A-B':>9} "
      f"{'paired p':>10} {'A 95% CI':>18}")
paired = {}
for name, col in [("decile-1 ratio", "d1"),
                  ("gate-region pAUC", "pauc_gate"),
                  ("absent-call floor ratio", "absent_ratio"),
                  ("full AUROC", "auc")]:
    s = G.dropna(subset=[f"{col}_A", f"{col}_B"])
    if len(s) < 20:
        continue
    a, b = s[f"{col}_A"], s[f"{col}_B"]
    p = wilcoxon(a - b).pvalue
    lo, hi = boot_ci(a)
    print(f"  {name:<32} {a.median():>9.4f} {b.median():>10.4f} "
          f"{(a - b).median():>9.4f} {p:>10.2e} [{lo:>7.4f},{hi:>7.4f}]")
    paired[col] = {"n": int(len(s)), "A": float(a.median()), "B": float(b.median()),
                   "diff": float((a - b).median()), "p": float(p),
                   "A_ci": [lo, hi]}
results["paired_A_vs_B"] = paired

# per-gene concordance: do the SAME genes gate on both platforms?
s = G.dropna(subset=["d1_A", "d1_B"])
if len(s) > 20:
    cc = spearmanr(s.d1_A, s.d1_B)
    agree = float(np.mean((s.d1_A < 1) == (s.d1_B < 1)))
    print()
    print(f"  per-gene concordance of the decile-1 ratio: "
          f"Spearman {cc.statistic:+.4f} (p {cc.pvalue:.2e})")
    print(f"  same direction (both depleted or both not) in "
          f"{100*agree:.1f}% of genes")
    print(f"  depleted in A: {100*(s.d1_A < 1).mean():.1f}%   "
          f"in B: {100*(s.d1_B < 1).mean():.1f}%")
    results["concordance"] = {"spearman_d1": float(cc.statistic),
                              "p": float(cc.pvalue), "sign_agreement": agree,
                              "frac_depleted_A": float((s.d1_A < 1).mean()),
                              "frac_depleted_B": float((s.d1_B < 1).mean())}

# ============================================================ Step 8
print()
print("=" * 78)
print("STEP 8 -- GATE-ELIGIBLE GENES ONLY")
print("=" * 78)
print("  Pan-essential genes are expressed everywhere and dependent everywhere;")
print("  no expression gate can act on them, and pooling them in dilutes the")
print("  only regime under test. Restrict to prevalence < 10% in both arms --")
print("  the same restriction test_run_screen_replication.py judges on.")
elig = G[(G.base_A < 0.10) & (G.base_B < 0.10)].copy()
print(f"  gate-eligible genes: {len(elig):,}")
sel_out = {"n": int(len(elig))}
if len(elig) > 20:
    print()
    print(f"  {'metric':<32} {'A GEO':>9} {'B DepMap':>10} {'paired p':>10} "
          f"{'A vs null p':>12}")
    for name, col, ref in [("decile-1 ratio", "d1", 1.0),
                           ("gate-region pAUC", "pauc_gate", 0.5),
                           ("absent-call floor ratio", "absent_ratio", 1.0)]:
        s2 = elig.dropna(subset=[f"{col}_A", f"{col}_B"])
        if len(s2) < 10:
            continue
        a, b = s2[f"{col}_A"], s2[f"{col}_B"]
        pp = wilcoxon(a - b).pvalue
        pn = wilcoxon(a - ref).pvalue
        print(f"  {name:<32} {a.median():>9.4f} {b.median():>10.4f} "
              f"{pp:>10.2e} {pn:>12.2e}")
        sel_out[col] = {"A": float(a.median()), "B": float(b.median()),
                        "n": int(len(s2)), "paired_p": float(pp),
                        "A_vs_null_p": float(pn),
                        "B_vs_null_p": float(wilcoxon(b - ref).pvalue)}
results["gate_eligible"] = sel_out

# ------------------------------------------------------------ by agreement
print()
print("  BY CROSS-PLATFORM AGREEMENT. If the gate is real but the array measures")
print("  some genes badly, the gate should appear in A exactly where the two")
print("  platforms agree, and fade where they do not. That pattern is the")
print("  signature of a MEASUREMENT limit rather than a false finding.")
sub = elig.dropna(subset=["rho_platform", "d1_A", "d1_B"])
band_rows = []
if len(sub) > 60:
    sub["band"] = pd.qcut(sub.rho_platform, 3,
                          labels=["low agree", "mid agree", "high agree"])
    print()
    print(f"  {'band':<12} {'genes':>6} {'rho':>7} {'d1 A':>8} {'d1 B':>8} "
          f"{'pAUC A':>8} {'pAUC B':>8}")
    for bnd in sub.band.cat.categories:
        s3 = sub[sub.band == bnd]
        pa = s3.pauc_gate_A.dropna()
        pb = s3.pauc_gate_B.dropna()
        print(f"  {str(bnd):<12} {len(s3):>6} {s3.rho_platform.median():>7.3f} "
              f"{s3.d1_A.median():>8.4f} {s3.d1_B.median():>8.4f} "
              f"{pa.median():>8.4f} {pb.median():>8.4f}")
        band_rows.append({"band": str(bnd), "n": int(len(s3)),
                          "median_rho": float(s3.rho_platform.median()),
                          "d1_A": float(s3.d1_A.median()),
                          "d1_B": float(s3.d1_B.median()),
                          "pauc_A": float(pa.median()),
                          "pauc_B": float(pb.median())})
results["by_platform_agreement"] = band_rows

# ------------------------------------------------------------ by prevalence
print()
print("  BY PREVALENCE. The gate is known to decay as prevalence rises.")
elig2 = G.dropna(subset=["d1_A", "d1_B", "base_B"]).copy()
elig2["pband"] = pd.cut(elig2.base_B, [0, .02, .05, .10, .25, 1.01],
                        labels=["<2%", "2-5%", "5-10%", "10-25%", ">25%"])
print()
print(f"  {'band':>8} {'genes':>6} {'d1 A':>8} {'d1 B':>8} {'A 95% CI':>18}")
pb_rows = []
for bnd in elig2.pband.cat.categories:
    s4 = elig2[elig2.pband == bnd]
    if len(s4) < 10:
        continue
    lo, hi = boot_ci(s4.d1_A)
    print(f"  {str(bnd):>8} {len(s4):>6} {s4.d1_A.median():>8.4f} "
          f"{s4.d1_B.median():>8.4f} [{lo:>7.4f},{hi:>7.4f}]")
    pb_rows.append({"band": str(bnd), "n": int(len(s4)),
                    "d1_A": float(s4.d1_A.median()),
                    "d1_B": float(s4.d1_B.median()), "d1_A_ci": [lo, hi]})
results["by_prevalence"] = pb_rows

# ============================================================ Step 9
print()
print("=" * 78)
print("STEP 9 -- NEGATIVE PREDICTIVE VALUE, the number the abstention rule uses")
print("=" * 78)
print(f"  {'cut':<22} {'NPV A':>9} {'uplift A':>10} {'NPV B':>9} {'uplift B':>10}")
npv_out = {}
for pctl in (5, 10, 20):
    r = {}
    line = f"  bottom {pctl:>2}% of score  "
    for t in "AB":
        s5 = G.dropna(subset=[f"npv{pctl}_{t}", f"npv_base_{t}"])
        npv = s5[f"npv{pctl}_{t}"].median()
        up = (s5[f"npv{pctl}_{t}"] - s5[f"npv_base_{t}"]).median()
        line += f"{npv:>9.4f} {up:>+10.4f} "
        r[t] = {"npv": float(npv), "uplift": float(up)}
    print(line)
    npv_out[f"bottom_{pctl}pct"] = r
results["npv"] = npv_out

# ============================================================ Step 10
print()
print("=" * 78)
print("STEP 10 -- VERDICT")
print("=" * 78)
ge = results["gate_eligible"]
d1A = ge.get("d1", {}).get("A", np.nan)
d1B = ge.get("d1", {}).get("B", np.nan)
paA = ge.get("pauc_gate", {}).get("A", np.nan)
paB = ge.get("pauc_gate", {}).get("B", np.nan)
pA_null = ge.get("d1", {}).get("A_vs_null_p", np.nan)
pB_null = ge.get("d1", {}).get("B_vs_null_p", np.nan)

A_gates = (d1A < 0.95 and pA_null < 0.05)
B_gates = (d1B < 0.95 and pB_null < 0.05)
A_pauc = paA > 0.52
B_pauc = paB > 0.52

if A_gates and B_gates:
    verdict = ("gate_REPLICATES_on_independent_expression_platform"
               if A_pauc and B_pauc else
               "gate_replicates_partially__depletion_yes_pauc_weak")
    note = ("Decile-1 depletion is present on GEO array expression as well as "
            "DepMap RNA-seq, on the same lines and labels. The gate is a "
            "property of gene expression, not of DepMap's pipeline.")
elif B_gates and not A_gates:
    verdict = "PLATFORM_DEPENDENCE__gate_does_not_transfer_to_array"
    note = ("On identical lines and labels the gate is present in DepMap "
            "RNA-seq and absent in GEO arrays. This is a scope condition and "
            "must be disclosed: the gate is calibrated to RNA-seq.")
elif not B_gates:
    verdict = "UNDERPOWERED__no_platform_conclusion_available"
    note = ("The gate does not appear even in the DepMap arm restricted to "
            "these lines, so this line set cannot resolve it either way. "
            "Nothing about GEO can be concluded from this run.")
else:
    verdict = "inconclusive"
    note = "Mixed evidence; read the per-band tables before writing anything."
results["verdict"] = verdict

print(f"  cross-platform expression agreement : rho = {RP.median():.3f} (median gene)")
print(f"  gate-eligible genes                 : {ge.get('n', 0):,}")
print(f"  decile-1 ratio   A GEO {d1A:.4f} (p {pA_null:.2e})   "
      f"B DepMap {d1B:.4f} (p {pB_null:.2e})")
print(f"  gate pAUC        A GEO {paA:.4f}                B DepMap {paB:.4f}")
print(f"  power reference  C DepMap full cohort, decile-1 "
      f"{summary['d1']['C_depmap_full']:.4f} over {len(L_full):,} lines")
print()
print(f"  VERDICT: {verdict}")
print(f"  {note}")

out = OUTPUTS / f"test_run_geo_gate_replication_{GSE}_results.json"
out.write_text(json.dumps(results, indent=2, default=float), encoding="utf-8")
print(f"\nwrote {out}")
print(f"wrote {out_pg}")
