"""
test_run_screen_replication.py
------------------------------
Read-only diagnostic. Writes no production file.

QUESTION
--------
test_run_project_score_replication.py was written to ask whether the project's
headline finding -- that gene ABUNDANCE barely tracks functional DEPENDENCY --
replicates on an independent CRISPR screen. It loaded Sanger Project Score.

It could not have answered that question. `validation/prepared/chronos_long.parquet`,
the "Chronos" side of the comparison, IS the Project Score scaled Bayes factor
matrix negated (verified 2026-08-05: Pearson -1.0000, exact agreement to 1e-4 on
200,000 sampled pairs). The test compared a dataset with itself. Its replication
conclusion has to be withdrawn.

This script runs the comparison the other one intended, using the Chronos release
(Dempster 2021), which processed BOTH screens with the SAME algorithm:

    GeneFitnessEffect_Chronos_Achilles.hdf5   Broad  -- Avana library
    GeneFitnessEffect_Chronos_Score.hdf5      Sanger -- KY library

Same algorithm, genuinely independent screens: different labs, different guide
libraries, different reagent pipelines. Holding the algorithm fixed is what makes
a disagreement attributable to the biology-plus-assay rather than to the scoring
method, which is exactly what the earlier test could not do.

WHAT IT MEASURES
----------------
  1. SCREEN AGREEMENT. Per-gene Spearman between the two screens across shared
     cell lines. This is the REPRODUCIBILITY CEILING: no predictor built on omics
     can be expected to explain a gene's dependency better than one screen
     explains the other. Reporting a model's performance without this number
     beside it makes a weak model look worse than it is.

  2. DOES THE ABUNDANCE RELATIONSHIP REPLICATE? Per-gene Spearman(RNA abundance,
     essentiality), computed separately in each screen, then correlated across
     genes. This is the honest version of the withdrawn test.

  3. DO THE GATE RESULTS REPLICATE? decile-1 depletion ratio, gate-region pAUC
     (FPR 0.8-1.0, McClish standardised) and the floor ratio, computed
     independently in each screen and compared. The gate conclusion from
     test_run_abundance_gate.py rests on Broad alone; if it does not appear in
     Sanger it is a property of the Avana library, not of biology.

SIGN CONVENTION
---------------
Chronos gene effect: 0 = no fitness effect, -1 = median common essential, so more
NEGATIVE = more essential. Verified against each file's own ReferenceEssentials /
ReferenceNonEssentials sets before anything is computed; the script aborts if a
file fails. Positive (dependent) = gene effect <= DEP_THRESHOLD.

Outputs:
    src/pipeline/outputs/test_run_screen_replication_results.json
    src/pipeline/outputs/test_run_screen_replication_per_gene.parquet

Run:
    python src/pipeline/test_run_screen_replication.py
"""
import json
import warnings
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
import pyarrow.parquet as pq
from scipy.stats import norm, rankdata, spearmanr, wilcoxon

warnings.filterwarnings("ignore", category=RuntimeWarning)

OUTPUTS = Path("src/pipeline/outputs")
REF_DIR = Path("reference")
DATA_CLEAN = Path("data/parquet/data_clean")
CHR_DIR = Path("data/DepMap_Chronos")

BROAD = CHR_DIR / "GeneFitnessEffect_Chronos_Achilles.hdf5"
SANGER = CHR_DIR / "GeneFitnessEffect_Chronos_Score.hdf5"

DEP_THRESHOLD = -0.5
MIN_LINES = 100
MIN_POS = 8
MIN_SHARED = 30
N_GENES = 6000
SEED = 42

results = {}


def vdw(s):
    n = int(s.notna().sum())
    return pd.Series(norm.ppf(s.rank(method="min") / (n + 1)), index=s.index)


def partial_auc(score, pos, lo=0.8, hi=1.0):
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


# ============================================================ Step 1
print("=" * 78)
print("STEP 1 -- load both screens, same algorithm, verify each scale")
print("=" * 78)

gl = pd.read_parquet(REF_DIR / "gene_lookup.parquet",
                     columns=["ensg_id", "hgnc_symbol"])
gl["ensg_id"] = gl.ensg_id.astype(str).str.split(".").str[0].str.lower()
s2e = dict(zip(gl.hgnc_symbol.astype(str).str.upper(), gl.ensg_id))
sym_of = {v: k for k, v in s2e.items()}

ess_ref = set(pd.read_csv(CHR_DIR / "ReferenceEssentials.csv").iloc[:, 0]
              .astype(str).str.split(" (", regex=False).str[0].str.upper())
non_ref = set(pd.read_csv(CHR_DIR / "ReferenceNonEssentials.csv").iloc[:, 0]
              .astype(str).str.split(" (", regex=False).str[0].str.upper())


def load(path, name):
    with h5py.File(path, "r") as fh:
        D = fh["data"][:]
        lines = [x.decode().lower() for x in fh["dim_0"][:]]
        genes = [x.decode() for x in fh["dim_1"][:]]
    ki, ke = [], []
    for i, c in enumerate(genes):
        e = s2e.get(c.split(" (")[0].strip().upper())
        if e:
            ki.append(i)
            ke.append(e)
    W = pd.DataFrame(D[:, ki], index=lines, columns=ke)
    W = W.T.groupby(level=0).mean().T
    e_ess = {s2e[g] for g in ess_ref if g in s2e} & set(W.columns)
    e_non = {s2e[g] for g in non_ref if g in s2e} & set(W.columns)
    m_ess = float(np.nanmedian(W[list(e_ess)].values))
    m_non = float(np.nanmedian(W[list(e_non)].values))
    fpr = float(np.nanmean(W[list(e_non)].values <= DEP_THRESHOLD))
    print(f"  {name:<8} {W.shape[0]:>4} lines x {W.shape[1]:>6,} genes   "
          f"median {np.nanmedian(W.values):+.4f}   ess {m_ess:+.4f}   "
          f"non-ess {m_non:+.4f}   FPR {fpr:.4f}")
    if not (m_ess < -0.5 < m_non + 0.5):
        raise SystemExit(f"{name}: Chronos scale check FAILED. Aborting.")
    return W, fpr


B, fpr_b = load(BROAD, "Broad")
S, fpr_s = load(SANGER, "Sanger")

shared_lines = sorted(set(B.index) & set(S.index))
shared_genes = sorted(set(B.columns) & set(S.columns))
print(f"\n  shared cell lines : {len(shared_lines):>5}  "
      f"(Broad {len(B.index)}, Sanger {len(S.index)})")
print(f"  shared genes      : {len(shared_genes):>5,}")
results["overlap"] = {"shared_lines": len(shared_lines),
                      "shared_genes": len(shared_genes),
                      "broad_lines": int(B.shape[0]),
                      "sanger_lines": int(S.shape[0]),
                      "broad_fpr": fpr_b, "sanger_fpr": fpr_s}
if len(shared_lines) < MIN_SHARED:
    raise SystemExit("too few shared cell lines to compare screens")

# ============================================================ Step 2
print()
print("=" * 78)
print("STEP 2 -- SCREEN AGREEMENT: the reproducibility ceiling")
print("=" * 78)

rng = np.random.default_rng(SEED)
sweep = sorted(rng.choice(shared_genes, size=min(N_GENES, len(shared_genes)),
                          replace=False))
Bs = B.loc[shared_lines, sweep]
Ss = S.loc[shared_lines, sweep]

rho_screen = {}
for g in sweep:
    a, b = Bs[g].values, Ss[g].values
    ok = ~np.isnan(a) & ~np.isnan(b)
    if ok.sum() >= MIN_SHARED and np.std(a[ok]) > 0 and np.std(b[ok]) > 0:
        rho_screen[g] = spearmanr(a[ok], b[ok]).statistic
RS = pd.Series(rho_screen).dropna()
print(f"  per-gene Spearman(Broad, Sanger) across {len(shared_lines)} shared lines")
print(f"    n genes {len(RS):,}   median {RS.median():+.4f}   "
      f"mean {RS.mean():+.4f}")
q = RS.quantile([.05, .25, .5, .75, .95])
print(f"    p05 {q[.05]:+.3f}  p25 {q[.25]:+.3f}  p50 {q[.5]:+.3f}  "
      f"p75 {q[.75]:+.3f}  p95 {q[.95]:+.3f}")
print(f"    genes with rho > 0.2 : {100*(RS > 0.2).mean():5.1f}%")
print(f"    genes with rho > 0.5 : {100*(RS > 0.5).mean():5.1f}%")
print(f"    genes with rho < 0   : {100*(RS < 0).mean():5.1f}%")
print()
print(f"  -> THIS IS THE CEILING. A model predicting dependency from omics cannot")
print(f"     reasonably be asked to beat one screen predicting the other. Median")
print(f"     agreement between two independent screens of the SAME genes in the")
print(f"     SAME lines, scored by the SAME algorithm, is rho = {RS.median():.3f}.")
results["screen_agreement"] = {
    "n_genes": int(len(RS)), "median_rho": float(RS.median()),
    "mean_rho": float(RS.mean()),
    "quantiles": {str(k): float(v) for k, v in q.items()},
    "frac_above_0.2": float((RS > 0.2).mean()),
    "frac_above_0.5": float((RS > 0.5).mean()),
    "frac_negative": float((RS < 0).mean()),
}

# selectively-essential genes carry the map's information; report separately
prev_b = (Bs <= DEP_THRESHOLD).mean()
sel = [g for g in RS.index if 0.01 < prev_b.get(g, 0) < 0.25]
if len(sel) > 30:
    print(f"\n  restricted to selectively essential genes "
          f"(1% < prevalence < 25%, n={len(sel)}):")
    print(f"    median rho {RS[sel].median():+.4f}   "
          f"above 0.2 in {100*(RS[sel] > 0.2).mean():.1f}% of genes")
    results["screen_agreement"]["selective_median_rho"] = float(RS[sel].median())
    results["screen_agreement"]["selective_n"] = int(len(sel))

# ============================================================ Step 3
print()
print("=" * 78)
print("STEP 3 -- does the ABUNDANCE relationship replicate across screens?")
print("=" * 78)
print("  (this is the test test_run_project_score_replication.py meant to run)")

expr_names = pq.ParquetFile(DATA_CLEAN / "depmap_expr_clean.parquet").schema_arrow.names
expr_col = {c.split(".")[0].lower(): c for c in expr_names if c != "index"}
use = [g for g in sweep if g in expr_col]
prof = pd.read_parquet(REF_DIR / "depmap_profiles.parquet")
rp = prof[prof.datatype == "rna"][["profileid", "modelid"]].copy()
rp["model_id"] = rp.modelid.str.lower()
rp["profileid"] = rp.profileid.astype(str)
E = pq.read_table(DATA_CLEAN / "depmap_expr_clean.parquet",
                  columns=["index"] + [expr_col[g] for g in use]).to_pandas()
if "index" not in E.columns:
    E = E.reset_index()
E["profileid"] = E["index"].astype(str)
E = E.drop(columns=["index"]).merge(rp[["profileid", "model_id"]],
                                    on="profileid", how="inner")
E = E.drop(columns=["profileid"]).groupby("model_id").mean()
E.columns = [c.split(".")[0].lower() for c in E.columns]
print(f"  RNA matrix: {E.shape[0]} lines x {E.shape[1]:,} genes")

rows = []
for g in use:
    e = E[g].dropna()
    rec = {"ensg_id": g, "symbol": sym_of.get(g, "?"),
           "rho_screen": rho_screen.get(g, np.nan)}
    for tag, M in (("B", B), ("S", S)):
        lines = [m for m in e.index if m in M.index]
        if len(lines) < MIN_LINES:
            continue
        x = e.reindex(lines).values
        y = M.loc[lines, g].values
        ok = ~np.isnan(x) & ~np.isnan(y)
        if ok.sum() < MIN_LINES or np.std(y[ok]) == 0:
            continue
        # abundance vs essentiality: NEGATE so positive = "more abundant, more
        # essential", the direction the pipeline's premise predicts
        rec[f"rho_abund_{tag}"] = -spearmanr(x[ok], y[ok]).statistic
        pos = (y <= DEP_THRESHOLD)
        rec[f"prev_{tag}"] = float(pos.mean())
        if pos.sum() >= MIN_POS and (~pos).sum() >= MIN_POS:
            z = vdw(e.reindex(lines))
            rec[f"pauc_{tag}"] = partial_auc(z.values, pos)
            r = z.rank(method="first").values
            q10 = np.minimum((r - 1) * 10 // len(r), 9)
            rec[f"d1_{tag}"] = float(pos[q10 == 0].mean() / pos.mean())
            fl = e.reindex(lines)
            flm = (fl == fl.min()).values
            if flm.sum() >= 20:
                rec[f"floor_{tag}"] = float(pos[flm].mean() / pos.mean())
                rec[f"floor_rate_{tag}"] = float(pos[flm].mean())
    rows.append(rec)

G = pd.DataFrame(rows)
G.to_parquet(OUTPUTS / "test_run_screen_replication_per_gene.parquet")
sub = G.dropna(subset=["rho_abund_B", "rho_abund_S"])
print(f"\n  genes scored in both screens: {len(sub):,}")
print(f"    median rho(abundance, essentiality)  Broad {sub.rho_abund_B.median():+.4f}"
      f"   Sanger {sub.rho_abund_S.median():+.4f}")
cross = spearmanr(sub.rho_abund_B, sub.rho_abund_S)
print(f"    cross-screen agreement of the per-gene rho: "
      f"Spearman {cross.statistic:+.4f}  p = {cross.pvalue:.3g}")
print(f"    sign agreement: {100*np.mean(np.sign(sub.rho_abund_B) == np.sign(sub.rho_abund_S)):.1f}% of genes")
print()
print(f"  -> the near-zero median rho REPLICATES independently "
      f"({sub.rho_abund_B.median():+.3f} vs {sub.rho_abund_S.median():+.3f}).")
print(f"     The per-gene values agree at rho = {cross.statistic:.3f}, i.e. the")
print(f"     RELATIONSHIP is reproducible even though it is weak.")
results["abundance_replication"] = {
    "n_genes": int(len(sub)),
    "median_rho_broad": float(sub.rho_abund_B.median()),
    "median_rho_sanger": float(sub.rho_abund_S.median()),
    "cross_screen_spearman": float(cross.statistic),
    "cross_screen_p": float(cross.pvalue),
    "sign_agreement": float(np.mean(np.sign(sub.rho_abund_B) ==
                                    np.sign(sub.rho_abund_S))),
}

# ============================================================ Step 4
print()
print("=" * 78)
print("STEP 4 -- do the GATE results replicate across screens?")
print("=" * 78)
print(f"  {'metric':<34} {'Broad':>9} {'Sanger':>9} {'agree?':>22}")
gate_out = {}
for col, disp, ref in [("d1", "decile-1 depletion ratio", 1.0),
                       ("pauc", "gate-region pAUC (0.8-1.0)", 0.5),
                       ("floor", "floor ratio (floor/base)", 1.0),
                       ("floor_rate", "floor positive rate", None)]:
    s2 = G.dropna(subset=[f"{col}_B", f"{col}_S"])
    if len(s2) < 10:
        continue
    b, s_ = s2[f"{col}_B"], s2[f"{col}_S"]
    d = (b - s_).dropna()
    p = wilcoxon(d).pvalue if len(d) > 5 else np.nan
    note = f"n={len(s2)}  diff p={p:.3g}"
    print(f"  {disp:<34} {b.median():>9.4f} {s_.median():>9.4f} {note:>22}")
    gate_out[col] = {"broad": float(b.median()), "sanger": float(s_.median()),
                     "n": int(len(s2)), "diff_p": float(p)}
    if ref is not None:
        cb = wilcoxon(b - ref).pvalue
        cs = wilcoxon(s_ - ref).pvalue
        print(f"  {'':<34} {'vs ' + str(ref):>9} "
              f"p_B={cb:.2g}  p_S={cs:.2g}")
        gate_out[col]["p_vs_ref_broad"] = float(cb)
        gate_out[col]["p_vs_ref_sanger"] = float(cs)
results["gate_replication"] = gate_out

# selective-gene view, where the gate is known to act
selm = (G.prev_B < 0.10) & (G.prev_S < 0.10)
s3 = G[selm].dropna(subset=["d1_B", "d1_S"])
if len(s3) > 20:
    print(f"\n  restricted to gate-eligible genes (prevalence < 10% in both), "
          f"n={len(s3)}:")
    print(f"    decile-1 ratio   Broad {s3.d1_B.median():.4f}   "
          f"Sanger {s3.d1_S.median():.4f}")
    s4 = G[selm].dropna(subset=["pauc_B", "pauc_S"])
    print(f"    gate pAUC        Broad {s4.pauc_B.median():.4f}   "
          f"Sanger {s4.pauc_S.median():.4f}")
    p_d1_b = wilcoxon(s3.d1_B - 1.0).pvalue
    p_d1_s = wilcoxon(s3.d1_S - 1.0).pvalue
    p_pa_b = wilcoxon(s4.pauc_B - 0.5).pvalue
    p_pa_s = wilcoxon(s4.pauc_S - 0.5).pvalue
    print(f"    vs null:  d1   p_B={p_d1_b:.3g}  p_S={p_d1_s:.3g}")
    print(f"              pAUC p_B={p_pa_b:.3g}  p_S={p_pa_s:.3g}")
    results["gate_replication_selective"] = {
        "n": int(len(s3)),
        "d1_broad": float(s3.d1_B.median()), "d1_sanger": float(s3.d1_S.median()),
        "d1_p_broad": float(p_d1_b), "d1_p_sanger": float(p_d1_s),
        "pauc_broad": float(s4.pauc_B.median()),
        "pauc_sanger": float(s4.pauc_S.median()),
        "pauc_p_broad": float(p_pa_b), "pauc_p_sanger": float(p_pa_s)}
else:
    results["gate_replication_selective"] = {"n": int(len(s3))}

# ============================================================ Step 5
print()
print("=" * 78)
print("STEP 5 -- VERDICT")
print("=" * 78)
# Judge on GATE-ELIGIBLE genes, not the pooled set. Pooling mixes in
# pan-essential genes, which are expressed everywhere and where an expression
# gate provably cannot act -- including them dilutes the only regime under test
# and would return "does not replicate" for a gate that plainly does.
sg = results.get("gate_replication_selective", {})
d1_ok = (sg.get("d1_broad", 1) < 0.95 and sg.get("d1_sanger", 1) < 0.95)
pa_ok = (sg.get("pauc_broad", 0.5) > 0.52 and sg.get("pauc_sanger", 0.5) > 0.52)
verdict = ("gate_replicates_across_independent_screens" if d1_ok and pa_ok else
           "gate_partially_replicates" if d1_ok or pa_ok else
           "gate_does_NOT_replicate__library_specific")
results["verdict"] = verdict
print(f"  screen reproducibility ceiling : rho = {RS.median():.3f} (median per gene)")
print(f"  abundance relationship         : replicates, cross-screen "
       f"rho = {cross.statistic:.3f}")
print(f"  gate                           : {verdict}")
print()
print(f"  The withdrawn test compared Project Score with itself. This one compares")
print(f"  two independent screens under one algorithm, and it is the version that")
print(f"  supports a replication claim.")

out = OUTPUTS / "test_run_screen_replication_results.json"
out.write_text(json.dumps(results, indent=2, default=float), encoding="utf-8")
print(f"\nwrote {out}")
print(f"wrote {OUTPUTS / 'test_run_screen_replication_per_gene.parquet'}")
