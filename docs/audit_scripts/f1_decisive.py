"""F1 -- THE DECISIVE EXPERIMENT.

Runs the notebook's GUARDED lineage scorer through the same evaluation harness
the main pipeline used for the plain (unguarded) version:
src/pipeline/test_run_improvement_scan.py PART 1, real DepMap Chronos labels.

Three arms, identical labels, identical genes, identical lines:
  panel    van der Waerden normal scores across the whole panel   (pipeline baseline)
  lineage  vdW normal scores within tissue, NO guards             (what C.15 tested)
  guarded  robust z within tissue WITH the notebook's three guards
           (MIN_PEERS_FOR_Z, SILENT_FRAC silent-lineage check, MAD_FLOOR)

The guard makes some lines unscoreable. pAUC over different line sets is not a
like-for-like comparison, so both are reported:
  (a) all-lines  -- each arm on whatever it can score
  (b) MATCHED    -- every arm restricted to the lines the guarded arm retains
Paired bootstrap over genes on the matched comparison.
"""
import os, json, numpy as np, pandas as pd, h5py
import pyarrow.parquet as pq
from pathlib import Path
from scipy.stats import norm, wilcoxon
from scipy import stats
os.chdir(r"C:\Disertation\UoB-GeneTraceAI-25-26")
W = os.path.dirname(os.path.abspath(__file__))
import tx_lib as T

REF_DIR = Path("reference"); DATA_CLEAN = Path("data/parquet/data_clean")
DEP_THRESHOLD, MIN_LINES, MIN_POS, MIN_LINEAGE = -0.5, 100, 10, 15
N_GENES, SEED = 2000, 42


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


# ---------------- labels: real Chronos, with the release's own scale check ------
cp = Path("data/DepMap_Chronos/GeneFitnessEffect_Chronos_Achilles.hdf5")
with h5py.File(cp, "r") as fh:
    D = fh["data"][:]
    lines_h = [x.decode() for x in fh["dim_0"][:]]
    genes_h = [x.decode() for x in fh["dim_1"][:]]
gl0 = pd.read_parquet(REF_DIR / "gene_lookup.parquet", columns=["ensg_id", "hgnc_symbol"])
gl0["ensg_id"] = gl0.ensg_id.astype(str).str.split(".").str[0].str.lower()
s2e = dict(zip(gl0.hgnc_symbol.astype(str).str.upper(), gl0.ensg_id))
keep_i, keep_e = [], []
for i, c in enumerate(genes_h):
    e = s2e.get(c.split(" (")[0].strip().upper())
    if e:
        keep_i.append(i); keep_e.append(e)
Wm = pd.DataFrame(D[:, keep_i], index=[m.lower() for m in lines_h], columns=keep_e)
Wm = Wm.T.groupby(level=0).mean().T
ess = set(pd.read_csv(cp.parent / "ReferenceEssentials.csv").iloc[:, 0]
          .astype(str).str.split(" (", regex=False).str[0].str.upper())
non = set(pd.read_csv(cp.parent / "ReferenceNonEssentials.csv").iloc[:, 0]
          .astype(str).str.split(" (", regex=False).str[0].str.upper())
e_ess = {s2e[g] for g in ess if g in s2e} & set(Wm.columns)
e_non = {s2e[g] for g in non if g in s2e} & set(Wm.columns)
m_ess = float(np.nanmedian(Wm[list(e_ess)].values))
m_non = float(np.nanmedian(Wm[list(e_non)].values))
label_fpr = float(np.nanmean(Wm[list(e_non)].values <= DEP_THRESHOLD))
print(f"Chronos scale check: essentials median {m_ess:+.4f} (expect ~-1), "
      f"non-essentials {m_non:+.4f} (expect ~0), label FPR {label_fpr:.4f}")
assert m_ess < -0.5 < m_non + 0.5, "Chronos scale check FAILED"

ch = Wm.stack().rename("essentiality").reset_index()
ch.columns = ["model_id", "ensg_id", "essentiality"]
ch["dep"] = ch.essentiality <= DEP_THRESHOLD

si = pd.read_csv(cp.parent / "DepMapSampleInfo20Q2.csv", low_memory=False)
si["model_id"] = si.DepMap_ID.astype(str).str.lower()
tissue_of = dict(zip(si.model_id, si.lineage))
vc = pd.Series(tissue_of).value_counts()
big = set(vc[vc >= MIN_LINEAGE].index)
print(f"lineages: {len(vc)} total, {len(big)} with >= {MIN_LINEAGE} lines")

# ---------------- expression ----------------------------------------------------
expr_names = pq.ParquetFile(DATA_CLEAN / "depmap_expr_clean.parquet").schema_arrow.names
expr_col = {c.split(".")[0].lower(): c for c in expr_names if c != "index"}
rng = np.random.default_rng(SEED)
cands = sorted(set(ch.ensg_id.unique()) & set(expr_col))
sweep = sorted(rng.choice(cands, size=min(N_GENES, len(cands)), replace=False))
prof = pd.read_parquet(REF_DIR / "depmap_profiles.parquet")
rp = prof[prof.datatype == "rna"][["profileid", "modelid"]].copy()
rp["model_id"] = rp.modelid.str.lower(); rp["profileid"] = rp.profileid.astype(str)
E = pq.read_table(DATA_CLEAN / "depmap_expr_clean.parquet",
                  columns=["index"] + [expr_col[g] for g in sweep]).to_pandas()
if "index" not in E.columns:
    E = E.reset_index()
E["profileid"] = E["index"].astype(str)
E = E.drop(columns=["index"]).merge(rp[["profileid", "model_id"]], on="profileid", how="inner")
E = E.drop(columns=["profileid"]).groupby("model_id").mean()
E.columns = [c.split(".")[0].lower() for c in E.columns]
print(f"RNA matrix: {E.shape[0]} lines x {E.shape[1]} genes")

dep_by_gene = {g: s.set_index("model_id")["dep"]
               for g, s in ch[ch.ensg_id.isin(sweep)][["ensg_id", "model_id", "dep"]]
               .groupby("ensg_id")}


def guarded_z(e, tis):
    """The notebook's robust_z + three guards, applied within lineage."""
    out = pd.Series(np.nan, index=e.index)
    for lg, idx in e.groupby(tis).groups.items():
        v = e.loc[idx].to_numpy(float)
        if len(v) < T.MIN_PEERS_FOR_Z:                       # guard 1
            continue
        if (v > T.EXPRESSED_MIN).mean() < T.SILENT_FRAC:     # guard 2 (silent lineage)
            continue
        mad = stats.median_abs_deviation(v, scale="normal")
        if mad <= 0:                                          # guard 3
            continue
        out.loc[idx] = (v - np.median(v)) / max(mad, T.MAD_FLOOR)
    return out


rows = []
for g in sweep:
    if g not in dep_by_gene or g not in E.columns:
        continue
    dep = dep_by_gene[g]; e = E[g].dropna()
    lines = [m for m in e.index if m in dep.index and tissue_of.get(m) in big]
    if len(lines) < MIN_LINES:
        continue
    e = e.reindex(lines); y = dep.reindex(lines).astype(bool).values
    if y.sum() < MIN_POS or (~y).sum() < MIN_POS:
        continue
    tis = pd.Series([tissue_of[m] for m in lines], index=lines)
    z_panel = vdw(e)
    z_lin = e.groupby(tis).transform(
        lambda s: pd.Series(norm.ppf(s.rank(method="min") / (len(s) + 1)), index=s.index))
    z_grd = guarded_z(e, tis)

    keep = z_grd.notna().to_numpy()
    rec = {"ensg_id": g, "n_lines": len(lines), "n_pos": int(y.sum()),
           "base": float(y.mean()), "n_kept": int(keep.sum()),
           "frac_kept": float(keep.mean())}
    # (a) each arm on whatever it can score
    rec["pauc_panel"] = partial_auc(z_panel.values, y)
    rec["pauc_lineage"] = partial_auc(z_lin.values, y)
    rec["pauc_guarded"] = partial_auc(z_grd.values, y)
    # (b) MATCHED: all arms on the guard-retained lines only
    if keep.sum() >= MIN_LINES and y[keep].sum() >= MIN_POS and (~y[keep]).sum() >= MIN_POS:
        rec["m_panel"] = partial_auc(z_panel.values[keep], y[keep])
        rec["m_lineage"] = partial_auc(z_lin.values[keep], y[keep])
        rec["m_guarded"] = partial_auc(z_grd.values[keep], y[keep])
    rows.append(rec)

R = pd.DataFrame(rows)
print(f"\ngenes scored: {len(R)}   median lines/gene {R.n_lines.median():.0f}   "
      f"median prevalence {R.base.median():.3f}")
print(f"guard retains median {100*R.frac_kept.median():.1f}% of lines per gene; "
      f"drops all lines for {int((R.frac_kept==0).sum())} genes")


def boot(d, n=10000, seed=42):
    d = np.asarray(pd.Series(d).dropna())
    rg = np.random.default_rng(seed)
    m = d[rg.integers(0, len(d), size=(n, len(d)))].mean(axis=1)
    return np.percentile(m, [2.5, 97.5]), d


print("\n" + "=" * 74)
print("(a) ALL-LINES  -- each arm on whatever it can score (NOT like-for-like)")
print("=" * 74)
for nm in ("panel", "lineage", "guarded"):
    s = R[f"pauc_{nm}"].dropna()
    print(f"  {nm:<10} median pAUC {s.median():.4f}   n_genes {len(s)}")

print("\n" + "=" * 74)
print("(b) MATCHED -- every arm on the lines the GUARD retains (like-for-like)")
print("=" * 74)
M = R.dropna(subset=["m_panel", "m_lineage", "m_guarded"])
print(f"  genes with a valid matched comparison: {len(M)}")
for nm in ("panel", "lineage", "guarded"):
    print(f"  {nm:<10} median pAUC {M[f'm_{nm}'].median():.4f}")
print()
for a, b in [("guarded", "panel"), ("guarded", "lineage"), ("lineage", "panel")]:
    d = (M[f"m_{a}"] - M[f"m_{b}"]).dropna()
    ci, dd = boot(d)
    try:
        p = wilcoxon(dd).pvalue
    except Exception:
        p = np.nan
    star = "" if (ci[0] <= 0 <= ci[1]) else "  *"
    print(f"  {a:<8} - {b:<8} : {d.median():+.4f}  95% CI [{ci[0]:+.4f}, {ci[1]:+.4f}]"
          f"  p={p:.4g}  better in {100*(d>0).mean():.1f}% of genes{star}")
print("\n  * = 95% paired-bootstrap CI over genes excludes zero")

out = {"label_fpr": label_fpr, "n_genes": int(len(R)), "n_genes_matched": int(len(M)),
       "median_frac_lines_kept": float(R.frac_kept.median()),
       "all_lines": {nm: float(R[f"pauc_{nm}"].median()) for nm in
                     ("panel", "lineage", "guarded")},
       "matched": {nm: float(M[f"m_{nm}"].median()) for nm in
                   ("panel", "lineage", "guarded")}}
for a, b in [("guarded", "panel"), ("guarded", "lineage"), ("lineage", "panel")]:
    d = (M[f"m_{a}"] - M[f"m_{b}"]).dropna(); ci, dd = boot(d)
    out[f"{a}_minus_{b}"] = {"median": float(d.median()),
                             "ci_lo": float(ci[0]), "ci_hi": float(ci[1]),
                             "frac_better": float((d > 0).mean())}
json.dump(out, open(os.path.join(W, "f1_results.json"), "w"), indent=1)
R.to_parquet(os.path.join(W, "f1_per_gene.parquet"), index=False)
print("\nwritten f1_results.json")
