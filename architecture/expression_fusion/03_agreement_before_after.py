"""
expression_fusion/03_agreement_before_after.py
-------------------------------------------------
Re-runs the D1 cross-source agreement check from
01_denominator_and_agreement.py, swapping the old pooled/mis-scaled
`geo_expr` (warehouse) for the new per-series scale-corrected
`geo_expr_v2` (02_rebuild_geo_scale_corrected.py), on the SAME 250-gene
sample (seed=42) -- so the two numbers are directly comparable and the only
thing that changed is the GEO scale bug.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy import stats

import common as C

N_GENES = 250

C.banner(f"D1 agreement -- old (mis-scaled) vs new (per-series corrected) GEO, {N_GENES} genes")

con = C.connect()
genes = C.sample_genes(con, N_GENES)
gene_ids = genes.gene_id.tolist()

dm = C.fetch_depmap(con, gene_ids, verbose=False)
hp = C.fetch_hpa(con, gene_ids, verbose=False)
geo_old = C.fetch_geo(con, gene_ids, verbose=False)          # old: pooled detect_scale (the bug)
geo_new = C.fetch_geo_v2(gene_ids, verbose=False)             # new: per-series corrected
con.close()

print(f"old geo_expr (warehouse, pooled scale-detect): "
      f"{geo_old.model_id.nunique()} model_ids, {len(geo_old):,} rows, "
      f"value range [{geo_old.value.min():.2f}, {geo_old.value.max():.2f}]")
print(f"new geo_expr_v2 (per-series corrected):        "
      f"{geo_new.model_id.nunique()} model_ids, {len(geo_new):,} rows, "
      f"value range [{geo_new.value.min():.2f}, {geo_new.value.max():.2f}]\n")


def agreement(geo_df, label):
    long = pd.concat([dm.assign(source="depmap_expr"), hp.assign(source="hpa_rna"),
                      geo_df.assign(source="geo")], ignore_index=True)
    rows = []
    for gid, g in long.groupby("gene_id"):
        w = g.pivot_table(index="model_id", columns="source", values="value", aggfunc="median")
        for pair in [("depmap_expr", "geo"), ("hpa_rna", "geo")]:
            if pair[0] not in w or pair[1] not in w:
                continue
            inter = w[[pair[0], pair[1]]].dropna()
            if len(inter) >= 15:
                rho, _ = stats.spearmanr(inter[pair[0]], inter[pair[1]])
                rows.append({"gene_id": gid, "pair": f"{pair[0]}_{pair[1]}",
                            "n": len(inter), "rho": rho})
    df = pd.DataFrame(rows)
    print(f"--- {label} ---")
    for pair, g in df.groupby("pair"):
        s = g["rho"].dropna()
        n_med = g["n"].median()
        print(f"  {pair:22s} n_genes={len(s):4d}  median_panel_n={n_med:5.0f}  "
              f"median rho={s.median():.3f}  IQR=[{s.quantile(.25):.3f} - {s.quantile(.75):.3f}]")
    return df


old = agreement(geo_old, "OLD -- geo_expr (pooled scale-detect, likely double-logged)")
print()
new = agreement(geo_new, "NEW -- geo_expr_v2 (per-series scale-corrected)")

old.to_parquet(C.OUT / "agreement_geo_old.parquet", index=False)
new.to_parquet(C.OUT / "agreement_geo_new.parquet", index=False)
print(f"\nwrote {C.OUT / 'agreement_geo_old.parquet'} and agreement_geo_new.parquet")
