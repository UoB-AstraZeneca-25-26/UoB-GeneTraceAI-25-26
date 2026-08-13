"""A3: within-GEO correlations.

There is no processing-method column and one platform (gpl570), so the honest
proxy for a "processing stratum" is the GEO SERIES (gse_id) -- processing is a
per-submission property. This measures whether splitting GEO by series would
create independent observations or re-count one observation several times.
"""
import os, json, itertools, numpy as np, pandas as pd
from scipy import stats
os.chdir(r"C:\Disertation\UoB-GeneTraceAI-25-26")
W = os.path.dirname(os.path.abspath(__file__))
import tx_lib as T

N_GENES, SEED = 400, 42
con = T.connect()

dep_cols = set(T.cols_of(con, "depmap_expr"))
geo_cols = set(T.cols_of(con, "geo_expr"))
hpa_genes = set(x[0] for x in con.execute("SELECT DISTINCT lower(gene) FROM hpa_rna").fetchall())
genes_all = sorted(g for g in dep_cols if g.startswith("ensg") and g in geo_cols and g in hpa_genes)
rng = np.random.default_rng(SEED)
genes = sorted(rng.choice(genes_all, size=N_GENES, replace=False))
lin = T.lineage_map(con)

# ---- long form: one row per (gse_id, model_id, gene)
aggs = ",".join(f'median(ge."{g}") AS "{g}"' for g in genes)
q = (f"SELECT coalesce(gi.gse_id,'(no gse_id)') AS gse, ge.model_id AS model_id, "
     f"count(*) AS n_gsm, {aggs} "
     f"FROM geo_expr ge LEFT JOIN geo_info gi ON lower(ge.sample)=lower(gi.geo_accession) "
     f"WHERE ge.model_id IS NOT NULL GROUP BY 1,2")
d = con.execute(q).df()
print(f"rows (series x line): {len(d):,}")

# ---- how many lines appear under >1 series
per_line = d.groupby("model_id")["gse"].nunique()
print("\n=== lines by number of GEO series they appear in ===")
vc = per_line.value_counts().sort_index()
for k, v in vc.items():
    print(f"  in {k} series: {v:5d} lines")
print(f"  lines in >=2 series: {int((per_line >= 2).sum()):,} of {len(per_line):,} "
      f"({100*(per_line>=2).mean():.1f}%)")
print(f"  lines in >=3 series: {int((per_line >= 3).sum()):,}")

print("\n=== series stratum sizes (distinct cell lines per series) ===")
sz = d.groupby("gse")["model_id"].nunique().sort_values(ascending=False)
print(sz.to_string())
print(f"\n  series with <15 lines (below MIN_PEERS_FOR_LINEAGE): "
      f"{int((sz < 15).sum())} of {len(sz)}")
print(f"  series with <10 lines: {int((sz < 10).sum())};  <5 lines: {int((sz < 5).sum())}")

# ---- robust z within (series, lineage), notebook guards
d["lineage"] = d["model_id"].map(lin).fillna(T.NO_LINEAGE)
val = np.log2(d[genes].clip(lower=0) + 1)
d = pd.concat([d[["gse", "model_id", "lineage", "n_gsm"]], val], axis=1)

Z = pd.DataFrame(np.nan, index=d.index, columns=genes)
for (gse, lg), idx in d.groupby(["gse", "lineage"]).groups.items():
    pos = d.index.get_indexer(idx)
    sub = d.loc[idx, genes]
    n = sub.notna().sum()
    med = sub.median()
    mad = (sub - med).abs().median() * 1.4826
    fe = (sub > T.EXPRESSED_MIN).sum() / n.replace(0, np.nan)
    ok = (n >= T.MIN_PEERS_FOR_Z) & (fe >= T.SILENT_FRAC) & (mad > 0)
    z = (sub - med) / mad.clip(lower=T.MAD_FLOOR)
    Z.iloc[pos, :] = z.to_numpy() * np.where(ok.to_numpy(), 1.0, np.nan)
zl = pd.concat([d[["gse", "model_id"]], Z], axis=1).melt(
    id_vars=["gse", "model_id"], var_name="gene", value_name="z").dropna(subset=["z"])
print(f"\nnon-null series-level z values: {len(zl):,}")

# ---- pairwise correlation between series-specific z, pooled over (line, gene)
print("\n=== A3 correlation between SERIES-specific z, for lines shared by both ===")
rows = []
for a, b in itertools.combinations(sorted(zl.gse.unique()), 2):
    A = zl[zl.gse == a].set_index(["model_id", "gene"])["z"]
    B = zl[zl.gse == b].set_index(["model_id", "gene"])["z"]
    common = A.index.intersection(B.index)
    nl = len({m for m, _ in common})
    if len(common) < 50 or nl < 5:
        continue
    r, p = stats.spearmanr(A.loc[common], B.loc[common])
    rows.append({"series_a": a, "series_b": b, "n_lines": nl,
                 "n_pairs": len(common), "spearman": round(float(r), 4),
                 "p": float(p)})
rr = pd.DataFrame(rows).sort_values("spearman", ascending=False)
print(rr.to_string(index=False))
if len(rr):
    w = rr.n_pairs
    print(f"\n  series pairs measured        : {len(rr)}")
    print(f"  median  Spearman             : {rr.spearman.median():.4f}")
    print(f"  IQR                          : {rr.spearman.quantile(.25):.4f}-{rr.spearman.quantile(.75):.4f}")
    print(f"  pair-count-weighted mean     : {(rr.spearman*w).sum()/w.sum():.4f}")
    print(f"  fraction of pairs rho > 0.5  : {(rr.spearman>0.5).mean():.3f}")

# ---- benchmark: the same quantity BETWEEN genuinely different platforms
print("\n=== benchmark: DepMap-vs-HPA z correlation, same pooled statistic ===")
print("  (from a2_results.json: median per-gene Spearman 0.863)")

json.dump({"lines_by_n_series": {int(k): int(v) for k, v in vc.items()},
           "n_lines_ge2_series": int((per_line >= 2).sum()),
           "n_lines_ge3_series": int((per_line >= 3).sum()),
           "series_sizes": {k: int(v) for k, v in sz.items()},
           "series_pair_corr": rows},
          open(os.path.join(W, "a3_results.json"), "w"), indent=1)
print("\nwritten a3_results.json")
