"""
expression_fusion/04_isolated_scale_bug_test.py
---------------------------------------------------
03_agreement_before_after.py compared two GEO sources that differ in TWO ways
at once -- the scale-detection granularity (the thing under test) AND the
underlying sample population (warehouse geo_expr vs the 15 freshly-fetched
series) -- so its result (agreement went DOWN after the "fix") cannot be
attributed to the scale bug alone.

This script isolates the one variable. Both arms are built from the exact
same 15-series raw data (data/geo/), the exact same probe->gene mapping, and
the exact same gsm->model_id resolution -- the ONLY difference is:

  NAIVE      pool every series' raw values for a gene together FIRST, run
             detect_scale ONCE on the pool, apply that one decision to
             everything -- reproducing the bug in both the shared notebook's
             fetch_gene_expression and this project's own first-draft
             fetch_geo.
  CORRECTED  transform each series with its own manifest-declared scale
             BEFORE pooling (what 02_rebuild_geo_scale_corrected.py does).

Same genes, same model_ids, same everything else. Whatever difference shows
up now is the scale-bug's effect alone.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

import common as C

N_GENES = 250
DATA_GEO = C.ROOT / "data" / "geo"

C.banner(f"Isolated test -- naive pooled-then-detect vs per-series-then-pool, {N_GENES} genes")

con = C.connect()
genes = C.sample_genes(con, N_GENES)
gene_ids = set(genes.gene_id.tolist())

gene = con.execute("SELECT gene_id, gene_names FROM main.gene").df()
pairs = gene.explode("gene_names").dropna(subset=["gene_names"])
pairs["name_lower"] = pairs["gene_names"].astype(str).str.lower()
counts = pairs.groupby("name_lower")["gene_id"].nunique()
unambiguous = set(counts[counts == 1].index)
symbol_to_ensg = (pairs[pairs.name_lower.isin(unambiguous)]
                  .drop_duplicates("name_lower")
                  .set_index("name_lower")["gene_id"].to_dict())

gsm_map = con.execute(
    "SELECT geo_accession, model_id FROM main.geo_info WHERE model_id IS NOT NULL"
).df().drop_duplicates("geo_accession").set_index("geo_accession")["model_id"]

dm = C.fetch_depmap(con, list(gene_ids), verbose=False)
hp = C.fetch_hpa(con, list(gene_ids), verbose=False)
con.close()

# =============================================================== load raw, per series, unscaled
raw_frames = []
for d in sorted(p for p in DATA_GEO.iterdir() if p.is_dir()):
    acc = d.name
    mat_p, map_p, man_p = d / f"{acc}_matrix.parquet", d / f"{acc}_probe_map.parquet", d / f"{acc}_manifest.json"
    if not (mat_p.exists() and map_p.exists() and man_p.exists()):
        continue
    manifest = json.loads(man_p.read_text())
    mat = pd.read_parquet(mat_p)
    pmap = pd.read_parquet(map_p)
    long = mat.reset_index().rename(columns={"ID_REF": "probe_id"}).melt(
        id_vars="probe_id", var_name="gsm", value_name="raw_value").dropna(subset=["raw_value"])
    long["gsm"] = long["gsm"].str.lower()
    long = long.merge(pmap, on="probe_id", how="inner")
    long["name_lower"] = long["gene_symbol"].astype(str).str.lower()
    long["gene_id"] = long["name_lower"].map(symbol_to_ensg)
    long = long.dropna(subset=["gene_id"])
    long = long[long.gene_id.isin(gene_ids)]
    if not len(long):
        continue
    long = long.groupby(["gene_id", "gsm"], as_index=False)["raw_value"].median()
    long["series"] = acc
    long["series_detected_scale"] = manifest["scale_detected"]["detected_scale"]
    raw_frames.append(long)

raw = pd.concat(raw_frames, ignore_index=True)
raw["model_id"] = raw["gsm"].map(gsm_map)
raw = raw.dropna(subset=["model_id"])
print(f"raw (gene, gsm) rows across {raw.series.nunique()} series, resolved to model_id: {len(raw):,}\n")


def detect_scale(v):
    v = np.asarray(v, float)
    v = v[np.isfinite(v)]
    if not len(v):
        return "unknown"
    if v.min() >= -10 and np.nanpercentile(v, 99) < 25 and v.max() <= 50:
        return "log2"
    return "linear"


# =============================================================== NAIVE: pool per gene, detect once
naive_rows = []
for gid, g in raw.groupby("gene_id"):
    scale = detect_scale(g["raw_value"])
    v = np.log2(np.clip(g["raw_value"], 0, None) + 1) if scale == "linear" else g["raw_value"]
    naive_rows.append(pd.DataFrame({"gene_id": gid, "model_id": g["model_id"], "value": v}))
naive = pd.concat(naive_rows, ignore_index=True)
naive = naive.groupby(["gene_id", "model_id"], as_index=False)["value"].median()

# =============================================================== CORRECTED: transform per series first
raw["value"] = raw["raw_value"]
is_linear = raw["series_detected_scale"] == "linear"
raw.loc[is_linear, "value"] = np.log2(np.clip(raw.loc[is_linear, "raw_value"], 0, None) + 1)
corrected = raw.groupby(["gene_id", "model_id"], as_index=False)["value"].median()

print(f"naive:     {naive.model_id.nunique()} model_ids, {len(naive):,} (gene, model_id) rows, "
      f"value range [{naive.value.min():.2f}, {naive.value.max():.2f}]")
print(f"corrected: {corrected.model_id.nunique()} model_ids, {len(corrected):,} (gene, model_id) rows, "
      f"value range [{corrected.value.min():.2f}, {corrected.value.max():.2f}]\n")


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
            if len(inter) >= 15 and inter[pair[0]].nunique() > 1 and inter[pair[1]].nunique() > 1:
                rho, _ = stats.spearmanr(inter[pair[0]], inter[pair[1]])
                rows.append({"gene_id": gid, "pair": f"{pair[0]}_{pair[1]}", "n": len(inter), "rho": rho})
    df = pd.DataFrame(rows)
    print(f"--- {label} ---")
    for pair, g in df.groupby("pair"):
        s = g["rho"].dropna()
        print(f"  {pair:22s} n_genes={len(s):4d}  median_panel_n={g['n'].median():5.0f}  "
              f"median rho={s.median():.3f}  IQR=[{s.quantile(.25):.3f} - {s.quantile(.75):.3f}]")
    return df


agreement(naive, "NAIVE -- pool across series, detect scale once (the bug)")
print()
agreement(corrected, "CORRECTED -- detect+transform per series, pool after (the fix)")

# =============================================================== direct value comparison
same = naive.merge(corrected, on=["gene_id", "model_id"], suffixes=("_naive", "_corrected"))
same["diff"] = (same.value_naive - same.value_corrected).abs()
print(f"\n{len(same):,} identical (gene, model_id) pairs scored by both.")
print(f"fraction where naive and corrected values differ (>0.01): "
      f"{(same['diff'] > 0.01).mean():.1%}")
print(f"median |difference| where they differ: "
      f"{same.loc[same['diff'] > 0.01, 'diff'].median():.3f}")
