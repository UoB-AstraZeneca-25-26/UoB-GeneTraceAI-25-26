"""A2: correlation between RAW values and between the notebook's Z-SCORES.

The notebook's stratified z removes platform-specific variance, which is the
main thing making the sources differ -- so the z correlation is expected to be
HIGHER than the raw correlation. That is the quantity the Stouffer combination
actually assumes is zero.

Vectorised over many genes; applies the notebook's own guards (MIN_PEERS_FOR_Z,
SILENT_FRAC, MAD_FLOOR) so the z's measured here are the z's it combines.
"""
import os, json, numpy as np, pandas as pd, duckdb
from scipy import stats
os.chdir(r"C:\Disertation\UoB-GeneTraceAI-25-26")
W = os.path.dirname(os.path.abspath(__file__))
import tx_lib as T

N_GENES = 400
SEED = 42
con = T.connect()

dep_cols = set(T.cols_of(con, "depmap_expr"))
geo_cols = set(T.cols_of(con, "geo_expr"))
genes_all = sorted(g for g in dep_cols if g.startswith("ensg") and g in geo_cols)
hpa_genes = set(x[0] for x in con.execute(
    "SELECT DISTINCT lower(gene) FROM hpa_rna").fetchall())
genes_all = [g for g in genes_all if g in hpa_genes]
print(f"genes present in all three warehouse tables: {len(genes_all):,}")
rng = np.random.default_rng(SEED)
genes = sorted(rng.choice(genes_all, size=min(N_GENES, len(genes_all)), replace=False))
print(f"sampled {len(genes)} genes (seed {SEED})")

lin = T.lineage_map(con)
print(f"lineage labels: {lin.nunique()} distinct, {len(lin):,} models")

# ---------------- build the three (model_id x gene) matrices, notebook transforms
sel = ",".join(f'"{g}"' for g in genes)

DEP = con.execute(f"SELECT model_id,{sel} FROM depmap_expr").df().set_index("model_id")
DEP = DEP.groupby(level=0).median()                       # replicate collapse (cell 8)

GEO = con.execute(f"SELECT model_id,{','.join(f'median(\"{g}\") AS \"{g}\"' for g in genes)} "
                  f"FROM geo_expr GROUP BY model_id").df().set_index("model_id")
GEO = np.log2(GEO.clip(lower=0) + 1)                      # log_geo=True

hp = con.execute(
    f"SELECT model_id, lower(gene) AS gene, avg(ntpm) AS v FROM hpa_rna "
    f"WHERE lower(gene) IN ({','.join(chr(39)+g+chr(39) for g in genes)}) "
    f"GROUP BY 1,2").df()
HPA = hp.pivot(index="model_id", columns="gene", values="v")
HPA = np.log2(HPA.clip(lower=0) + 1)                      # log_hpa=True

for nm, M in [("DepMap", DEP), ("GEO", GEO), ("HPA", HPA)]:
    print(f"  {nm}: {M.shape[0]:,} lines x {M.shape[1]:,} genes")


def robust_z_matrix(M, lineage):
    """Notebook robust_z + guards, applied per (lineage, gene), vectorised."""
    lg = pd.Series(lineage.reindex(M.index).fillna(T.NO_LINEAGE).values, index=range(len(M)))
    Z = pd.DataFrame(np.nan, index=M.index, columns=M.columns)
    for lname in lg.unique():
        pos = np.flatnonzero((lg == lname).to_numpy())
        sub = M.iloc[pos]
        n = sub.notna().sum()
        med = sub.median()
        mad = (sub - med).abs().median() * 1.4826
        frac_expressed = (sub > T.EXPRESSED_MIN).sum() / n.replace(0, np.nan)
        ok = (n >= T.MIN_PEERS_FOR_Z) & (frac_expressed >= T.SILENT_FRAC) & (mad > 0)
        denom = mad.clip(lower=T.MAD_FLOOR)
        z = (sub - med) / denom
        # NB: DataFrame.where aligns a Series on the INDEX, not the columns, so
        # broadcast the per-gene guard mask explicitly across rows.
        Z.iloc[pos, :] = z.to_numpy() * np.where(ok.to_numpy(), 1.0, np.nan)
    return Z


ZD, ZG, ZH = (robust_z_matrix(M, lin) for M in (DEP, GEO, HPA))
print("z matrices built")


def per_gene_corr(A, B, min_lines=20, method="spearman"):
    lines = sorted(set(A.index) & set(B.index))
    genes_ = sorted(set(A.columns) & set(B.columns))
    a, b = A.loc[lines, genes_], B.loc[lines, genes_]
    if method == "spearman":
        a = a.rank(axis=0, na_option="keep"); b = b.rank(axis=0, na_option="keep")
    m = a.notna() & b.notna()
    a, b = a.where(m), b.where(m)
    a, b = a - a.mean(), b - b.mean()
    rho = ((a * b).sum() / np.sqrt((a ** 2).sum() * (b ** 2).sum())
           ).replace([np.inf, -np.inf], np.nan)
    return rho[m.sum() >= min_lines].dropna(), len(lines)


def summarise(rho, label, nl):
    return {"pair": label, "n_lines_shared": nl, "n_genes": int(len(rho)),
            "median": round(float(rho.median()), 4),
            "q25": round(float(rho.quantile(.25)), 4),
            "q75": round(float(rho.quantile(.75)), 4),
            "frac_gt_0.5": round(float((rho > 0.5).mean()), 4)}


PAIRS = [("DepMap", "HPA", DEP, HPA, ZD, ZH),
         ("DepMap", "GEO", DEP, GEO, ZD, ZG),
         ("HPA", "GEO", HPA, GEO, ZH, ZG)]

raw_rows, z_rows, delta_rows = [], [], []
for a, b, RA, RB, ZA, ZB in PAIRS:
    rr, nl = per_gene_corr(RA, RB)
    zz, nlz = per_gene_corr(ZA, ZB)
    raw_rows.append(summarise(rr, f"{a}-{b}", nl))
    z_rows.append(summarise(zz, f"{a}-{b}", nlz))
    common = rr.index.intersection(zz.index)
    d = (zz[common] - rr[common])
    delta_rows.append({"pair": f"{a}-{b}", "n_genes": int(len(d)),
                       "median_delta": round(float(d.median()), 4),
                       "q25": round(float(d.quantile(.25)), 4),
                       "q75": round(float(d.quantile(.75)), 4),
                       "frac_z_higher": round(float((d > 0).mean()), 4)})

print("\n=== A2 RAW-VALUE correlation (within gene, across shared lines) ===")
print(pd.DataFrame(raw_rows).to_string(index=False))
print("\n=== A2 Z-SCORE correlation -- the quantity Stouffer assumes is zero ===")
print(pd.DataFrame(z_rows).to_string(index=False))
print("\n=== A2 DIFFERENCE (z minus raw) ===")
print(pd.DataFrame(delta_rows).to_string(index=False))

M = pd.DataFrame(1.0, index=["DepMap", "HPA", "GEO"], columns=["DepMap", "HPA", "GEO"])
for r in z_rows:
    a, b = r["pair"].split("-")
    M.loc[a, b] = M.loc[b, a] = r["median"]
print("\n=== 3x3 Z-SCORE correlation matrix R (median per-gene Spearman) ===")
print(M.to_string())

json.dump({"n_genes_sampled": len(genes), "seed": SEED,
           "raw": raw_rows, "z": z_rows, "delta": delta_rows,
           "R_zscore": M.to_dict()}, open(os.path.join(W, "a2_results.json"), "w"), indent=1)
for a, b, RA, RB, ZA, ZB in PAIRS:
    zz, _ = per_gene_corr(ZA, ZB)
    zz.rename("rho").to_frame().to_parquet(os.path.join(W, f"a2_z_{a}_{b}.parquet"))
print("\nwritten a2_results.json")

# ---- like-for-like: raw restricted to the SAME genes the z survived on -------
print("\n=== A2 LIKE-FOR-LIKE (raw and z on the identical surviving gene set) ===")
lfl = []
for a, b, RA, RB, ZA, ZB in PAIRS:
    rr, _ = per_gene_corr(RA, RB)
    zz, _ = per_gene_corr(ZA, ZB)
    c = rr.index.intersection(zz.index)
    lfl.append({"pair": f"{a}-{b}", "n_genes": len(c),
                "raw_median": round(float(rr[c].median()), 4),
                "z_median": round(float(zz[c].median()), 4),
                "raw_IQR": f"{rr[c].quantile(.25):.3f}-{rr[c].quantile(.75):.3f}",
                "z_IQR": f"{zz[c].quantile(.25):.3f}-{zz[c].quantile(.75):.3f}",
                "raw_frac_gt_.5": round(float((rr[c] > .5).mean()), 4),
                "z_frac_gt_.5": round(float((zz[c] > .5).mean()), 4),
                "median_paired_delta": round(float((zz[c] - rr[c]).median()), 4),
                "frac_z_higher": round(float((zz[c] > rr[c]).mean()), 4)})
lfl = pd.DataFrame(lfl)
print(lfl.to_string(index=False))
json.dump(lfl.to_dict("records"), open(os.path.join(W, "a2_like_for_like.json"), "w"), indent=1)
