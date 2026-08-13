"""A4: quantify the inflation.

Current rule      Z = sum(w z) / sqrt(sum(w^2))          [independence]
Correlation-aware Z = sum(w z) / sqrt(w' R w)            [Strube 1985 / Hartung 1999]

The RATIO Z_current / Z_corrected = sqrt(w'Rw)/sqrt(w'w) depends only on which
sources are present and their weights -- NOT on the data. So it is exact per
source-presence pattern, and its distribution across genes is just the mixture
of those patterns.
"""
import os, json, itertools, numpy as np, pandas as pd
from scipy import stats
os.chdir(r"C:\Disertation\UoB-GeneTraceAI-25-26")
W = os.path.dirname(os.path.abspath(__file__))
import tx_lib as T

A2 = json.load(open(os.path.join(W, "a2_results.json")))
R3 = A2["R_zscore"]
SRC = ["depmap_expr", "hpa_rna", "geo_expr"]
NAME = {"depmap_expr": "DepMap", "hpa_rna": "HPA", "geo_expr": "GEO"}
R = np.array([[R3[NAME[a]][NAME[b]] for b in SRC] for a in SRC])
print("measured z-score correlation matrix R:")
print(pd.DataFrame(R, index=[NAME[s] for s in SRC], columns=[NAME[s] for s in SRC]).round(4).to_string())


def infl(idx, w=None, Rm=R):
    """Z_current / Z_corrected for the given subset of sources."""
    k = len(idx)
    w = np.ones(k) if w is None else np.asarray(w, float)
    Rk = Rm[np.ix_(idx, idx)]
    return float(np.sqrt(w @ Rk @ w) / np.sqrt(w @ w))


print("\n=== A4 exact inflation by source-presence pattern (equal weights) ===")
rows = []
for k in (2, 3):
    for combo in itertools.combinations(range(3), k):
        r = infl(list(combo))
        rows.append({"sources": " + ".join(NAME[SRC[i]] for i in combo), "k": k,
                     "inflation_ratio": round(r, 4),
                     "effective_n_sources": round(k / r ** 2, 3),
                     "pct_overstated": round(100 * (r - 1), 1)})
print(pd.DataFrame(rows).to_string(index=False))

# ---- the five-input case, IF GEO were split by series ----
print("\n=== A4 five-input case: GEO split into 3 series (rho_within_GEO = 0.448) ===")
A3 = json.load(open(os.path.join(W, "a3_results.json")))
rho_geo = float(np.median([r["spearman"] for r in A3["series_pair_corr"]]))
print(f"  median within-GEO series-pair correlation = {rho_geo:.4f}")
n_geo = 3
m = 2 + n_geo
R5 = np.eye(m)
R5[0, 1] = R5[1, 0] = R[0, 1]                       # DepMap-HPA
for j in range(2, m):
    R5[0, j] = R5[j, 0] = R[0, 2]                   # DepMap-GEOseries
    R5[1, j] = R5[j, 1] = R[1, 2]                   # HPA-GEOseries
for i, j in itertools.combinations(range(2, m), 2):
    R5[i, j] = R5[j, i] = rho_geo                   # GEOseries-GEOseries
r5 = infl(list(range(m)), Rm=R5)
print(f"  5-input inflation ratio        = {r5:.4f}   ({100*(r5-1):.1f}% overstated)")
print(f"  effective independent sources  = {m / r5**2:.2f} of {m}")
r3 = infl([0, 1, 2])
print(f"  vs 3-source (GEO unsplit)      = {r3:.4f}   ({100*(r3-1):.1f}% overstated)")
print(f"  -> splitting GEO into {n_geo} series makes the overstatement "
      f"{'WORSE' if r5 > r3 else 'BETTER'} by {100*(r5-r3):.1f} pp")

# ---- empirical distribution across real (gene, line) cells, real weights ----
print("\n=== A4 empirical distribution across real cells (sqrt(n_samples) weights) ===")
con = T.connect()
N_GENES, SEED = 200, 42
dep_cols = set(T.cols_of(con, "depmap_expr")); geo_cols = set(T.cols_of(con, "geo_expr"))
hpa_g = set(x[0] for x in con.execute("SELECT DISTINCT lower(gene) FROM hpa_rna").fetchall())
gall = sorted(g for g in dep_cols if g.startswith("ensg") and g in geo_cols and g in hpa_g)
genes = sorted(np.random.default_rng(SEED).choice(gall, size=N_GENES, replace=False))
lin = T.lineage_map(con)

sel = ",".join(f'"{g}"' for g in genes)
DEP = con.execute(f"SELECT model_id,{sel} FROM depmap_expr").df().set_index("model_id")
depw = DEP.groupby(level=0).size().rename("n")          # replicate profiles per line
DEP = DEP.groupby(level=0).median()
GEO = con.execute(f"SELECT model_id,{','.join(f'median(x.\"{g}\") AS \"{g}\"' for g in genes)}, "
                  f"count(*) AS n_gsm FROM geo_expr x WHERE x.model_id IS NOT NULL "
                  f"GROUP BY model_id").df().set_index("model_id")
geow = GEO.pop("n_gsm")
GEO = np.log2(GEO.clip(lower=0) + 1)
hp = con.execute(f"SELECT model_id, lower(gene) AS gene, avg(ntpm) AS v, count(*) AS n "
                 f"FROM hpa_rna WHERE lower(gene) IN "
                 f"({','.join(chr(39)+g+chr(39) for g in genes)}) GROUP BY 1,2").df()
HPA = np.log2(hp.pivot(index="model_id", columns="gene", values="v").clip(lower=0) + 1)
hpaw = hp.groupby("model_id")["n"].median()


def zmat(M):
    lg = pd.Series(lin.reindex(M.index).fillna(T.NO_LINEAGE).values, index=range(len(M)))
    Z = pd.DataFrame(np.nan, index=M.index, columns=M.columns)
    for ln in lg.unique():
        pos = np.flatnonzero((lg == ln).to_numpy()); sub = M.iloc[pos]
        n = sub.notna().sum(); med = sub.median()
        mad = (sub - med).abs().median() * 1.4826
        fe = (sub > T.EXPRESSED_MIN).sum() / n.replace(0, np.nan)
        ok = (n >= T.MIN_PEERS_FOR_Z) & (fe >= T.SILENT_FRAC) & (mad > 0)
        Z.iloc[pos, :] = ((sub - med) / mad.clip(lower=T.MAD_FLOOR)).to_numpy() * \
                         np.where(ok.to_numpy(), 1.0, np.nan)
    return Z


ZD, ZH, ZG = zmat(DEP), zmat(HPA), zmat(GEO)
for _M in (ZD, ZH, ZG):
    _M.drop(index=[i for i in _M.index if not isinstance(i, str)], inplace=True, errors="ignore")
lines = sorted(set(ZD.index) | set(ZH.index) | set(ZG.index))
present = pd.DataFrame({
    "depmap_expr": ZD.reindex(lines).notna().sum(axis=1) > 0,
    "hpa_rna":     ZH.reindex(lines).notna().sum(axis=1) > 0,
    "geo_expr":    ZG.reindex(lines).notna().sum(axis=1) > 0}).fillna(False)
wts = pd.DataFrame({
    "depmap_expr": np.sqrt(depw.reindex(lines).fillna(1.0)),
    "hpa_rna":     np.sqrt(hpaw.reindex(lines).fillna(1.0)),
    "geo_expr":    np.sqrt(geow.reindex(lines).fillna(1.0))}).fillna(1.0)

recs = []
for mid in lines:
    idx = [i for i, s in enumerate(SRC) if present.loc[mid, s]]
    if len(idx) < 2:
        continue
    w = wts.loc[mid, [SRC[i] for i in idx]].to_numpy(float)
    recs.append({"model_id": mid, "k": len(idx),
                 "pattern": "+".join(NAME[SRC[i]] for i in idx),
                 "ratio": infl(idx, w)})
E = pd.DataFrame(recs)
print(f"lines with >=2 sources: {len(E):,}")
print(E.groupby("pattern")["ratio"].agg(
    n="size", median="median", q25=lambda s: s.quantile(.25),
    q75=lambda s: s.quantile(.75)).round(4).to_string())
print(f"\noverall: median inflation {E.ratio.median():.4f}  "
      f"IQR {E.ratio.quantile(.25):.4f}-{E.ratio.quantile(.75):.4f}  max {E.ratio.max():.4f}")
print(f"median effective sources: {(E.k / E.ratio**2).median():.2f}")

json.dump({"R": R.tolist(), "patterns": rows, "rho_within_geo": rho_geo,
           "five_input_ratio": r5, "three_source_ratio": r3,
           "empirical_median_ratio": float(E.ratio.median()),
           "empirical_by_pattern": E.groupby("pattern")["ratio"].median().round(4).to_dict()},
          open(os.path.join(W, "a4_results.json"), "w"), indent=1)
E.to_parquet(os.path.join(W, "a4_per_line.parquet"), index=False)
print("\nwritten a4_results.json")
