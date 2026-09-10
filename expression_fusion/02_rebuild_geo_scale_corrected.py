"""
expression_fusion/02_rebuild_geo_scale_corrected.py
------------------------------------------------------
Rebuilds a GEO expression layer from the 15 series fetched by
`src/scripts/geo_fetch.py` into `data/geo/`, fixing the bug found in the
previous diagnostic: BOTH the shared notebook's `fetch_gene_expression` (pools
across series per gene) and this project's own `expression_fusion/common.py`
`fetch_geo` (pools across the whole gene sample) call `detect_scale` on a pile
that mixes RMA (log2, 12 of 15 series, 94.9% of samples) and MAS5 (linear, 3
of 15 series) together. A few large MAS5 values are enough to make the whole
pile read as "linear," so the already-log2 RMA majority gets log2(x+1)
applied a second time.

THE FIX
-------
Apply the scale transform PER SERIES, using that series' own manifest
(`data/geo/<acc>/<acc>_manifest.json`, `scale_detected.detected_scale`) --
which was itself computed on that series alone, not pooled -- BEFORE any
series are combined with each other or with other genes. Only after every
series is on the same scale is anything pooled.

WHAT THIS DOES NOT FIX
-----------------------
The GSM -> model_id resolution gap (only 774 of 3,267 rows in `geo_info` map
to a model_id at all, and GSE57083 itself is only 19.6% resolved) is a
SEPARATE problem -- name/accession matching, not scale -- and is untouched
here. This script can only correctly scale the samples that were already
resolvable; it does not add coverage on that axis.

Output: expression_fusion/outputs/geo_expr_v2.parquet
        long: gene_id, model_id, value (log2 scale, correctly single-logged),
        n_samples (GSM count collapsed into this model_id/gene, across all
        series that measured it).
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd

import common as C

DATA_GEO = C.ROOT / "data" / "geo"

C.banner("REBUILD -- per-series scale-corrected GEO layer")

con = C.connect()

# =============================================================== symbol -> ensg
gene = con.execute("SELECT gene_id, gene_names FROM main.gene").df()
pairs = gene.explode("gene_names").dropna(subset=["gene_names"])
pairs["name_lower"] = pairs["gene_names"].astype(str).str.lower()
counts = pairs.groupby("name_lower")["gene_id"].nunique()
unambiguous = set(counts[counts == 1].index)
n_ambiguous_names = int((counts > 1).sum())
symbol_to_ensg = (pairs[pairs.name_lower.isin(unambiguous)]
                  .drop_duplicates("name_lower")
                  .set_index("name_lower")["gene_id"].to_dict())
print(f"gene-name -> ensg lookup: {len(symbol_to_ensg):,} unambiguous names "
      f"({n_ambiguous_names:,} ambiguous names excluded)")

# =============================================================== gsm -> model_id
# main.geo_expr's OWN sample/model_id/is_ambiguous columns are the better
# resolution (2,310 of 3,269 samples, 70.7%) -- not geo_info.geo_accession
# (only 774). The notebook's own fetch_gene_expression relies on this richer
# mapping ('SELECT model_id, ... FROM geo_expr GROUP BY model_id'), which is
# why its printed run reached ~590 GEO lines while a geo_info-based join
# reaches far fewer. Use the same resolution here for a fair rebuild.
gsm_map = con.execute(
    "SELECT sample, model_id FROM main.geo_expr "
    "WHERE model_id IS NOT NULL AND is_ambiguous = FALSE"
).df().drop_duplicates("sample").set_index("sample")["model_id"]
print(f"gsm -> model_id map: {len(gsm_map):,} resolved samples "
      f"(via geo_expr's own resolution, is_ambiguous=FALSE; scale is the only "
      f"thing this script changes)\n")

# =============================================================== per-series load
series_dirs = sorted(p for p in DATA_GEO.iterdir() if p.is_dir())
frames = []
series_report = []

for d in series_dirs:
    acc = d.name
    mat_p, map_p, man_p = d / f"{acc}_matrix.parquet", d / f"{acc}_probe_map.parquet", d / f"{acc}_manifest.json"
    if not (mat_p.exists() and map_p.exists() and man_p.exists()):
        continue
    manifest = json.loads(man_p.read_text())
    detected = manifest["scale_detected"]["detected_scale"]

    mat = pd.read_parquet(mat_p)          # probe_id (index) x gsm
    pmap = pd.read_parquet(map_p)         # probe_id, gene_symbol (already exploded, one gene/row)

    long = mat.reset_index().rename(columns={"ID_REF": "probe_id"}).melt(
        id_vars="probe_id", var_name="gsm", value_name="value").dropna(subset=["value"])
    long["gsm"] = long["gsm"].str.lower()
    long = long.merge(pmap, on="probe_id", how="inner")
    long["name_lower"] = long["gene_symbol"].astype(str).str.lower()
    long["gene_id"] = long["name_lower"].map(symbol_to_ensg)
    long = long.dropna(subset=["gene_id"])

    # --- THE FIX: transform per series, using that series' own detected scale
    if detected == "linear":
        long["value"] = np.log2(np.clip(long["value"], 0, None) + 1)
        applied = "log2(x+1) applied"
    elif detected in ("log_like", "log_like_signed"):
        applied = "left as-is"
    else:
        applied = f"left as-is (detected={detected!r}, ambiguous/insufficient -- not blindly converted)"

    # probe -> gene: collapse multiple probes per (gene, gsm) by median
    long = long.groupby(["gene_id", "gsm"], as_index=False)["value"].median()

    # gsm -> model_id
    long["model_id"] = long["gsm"].map(gsm_map)
    n_before = long["gsm"].nunique()
    long = long.dropna(subset=["model_id"])
    n_resolved = long["model_id"].nunique()

    frames.append(long[["gene_id", "model_id", "value"]])
    series_report.append({
        "series": acc, "detected_scale": detected, "transform": applied,
        "n_samples_total": n_before, "n_samples_resolved_to_model_id": n_resolved,
        "n_gene_value_rows": len(long),
    })
    print(f"  {acc:10s} detected={detected:14s} -> {applied:60s} "
          f"{n_resolved}/{n_before} samples resolved")

report_df = pd.DataFrame(series_report)
report_df.to_parquet(C.OUT / "geo_v2_series_report.parquet", index=False)

# =============================================================== combine across series
all_long = pd.concat(frames, ignore_index=True)
# a (gene, model_id) pair can be measured by more than one series -- collapse
# by median, weight (n_samples) is the GSM count contributing, same replicate
# principle the shared notebook already uses for GEO.
combined = all_long.groupby(["gene_id", "model_id"]).agg(
    value=("value", "median"), n_samples=("value", "size")).reset_index()

out_path = C.OUT / "geo_expr_v2.parquet"
combined.to_parquet(out_path, index=False)

print(f"\ncombined: {len(combined):,} (gene, model_id) rows, "
      f"{combined.gene_id.nunique():,} genes, {combined.model_id.nunique():,} model_ids")
print(f"wrote {out_path}")
print(f"wrote {C.OUT / 'geo_v2_series_report.parquet'}")
con.close()
