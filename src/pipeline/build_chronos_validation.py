"""
build_chronos_validation.py
----------------------------
Independent validation check for "unknown"-class genes (no curated GDSC
regime, no COSMIC oncogene/TSG role) using genome-wide CRISPR essentiality
(Chronos) as ground truth for whether abundance (core_score) actually
tracks functional dependency for that gene.

This does NOT feed Chronos into the ranking score itself (that was tried
and rejected earlier -- BCL2's hit@20 collapsed to 0 when Chronos was
blended into core_score as a 3rd layer). This is a read-only diagnostic:
for each gene, does the cell-line ordering implied by core_score agree
with the cell-line ordering implied by real essentiality data?

Per-gene classification (Spearman rho between core_score and Chronos
essentiality percentile, across cell lines both cover):
    validated : rho >  0.1, p < 0.05  -- abundance ranking tracks real
                dependency for this gene; safe to trust
    inverted  : rho < -0.1, p < 0.05  -- abundance ranking is likely
                backwards for this gene (e.g. low abundance = high
                dependency, a threshold/buffering effect)
    none      : no significant relationship -- stays honestly unknown,
                but now a CHECKED unknown rather than an assumed one

Role in pipeline: read by build_full_predictions.py to adjust confidence
tier and signal_quality for unknown-class genes only. Curated (measured)
and COSMIC-prior genes are untouched by this file.

Output: src/pipeline/outputs/chronos_validation.parquet
    (ensg_id, n, rho, pval, chronos_check)
"""
import time
import numpy as np
import pandas as pd
from pathlib import Path
from scipy.stats import t as tdist

OUTPUTS = Path("src/pipeline/outputs")
VAL     = Path("validation/prepared")

RHO_THRESHOLD = 0.1
P_THRESHOLD   = 0.05
MIN_N         = 30

t0 = time.time()

pred = pd.read_parquet(OUTPUTS / "predictions_with_confidence.parquet",
                        columns=["model_id", "ensg_id", "core_score", "class"])
pred["model_id"] = pred["model_id"].str.upper()
unknown = pred[pred["class"] == "unknown"][["model_id", "ensg_id", "core_score"]]
print(f"[{time.time()-t0:.1f}s] unknown-class rows: {len(unknown):,}  genes: {unknown.ensg_id.nunique():,}")

ch = pd.read_parquet(VAL / "chronos_long.parquet")
ib = pd.read_parquet(VAL / "id_bridge.parquet")
ch = ch.merge(ib, on="sanger_model_id", how="inner")
ch["model_id"] = ch["model_id"].str.upper()
ch["ensg_id"]  = ch["ensg_id"].str.upper()
# More negative raw essentiality = more essential -> invert so high chronos_pct = essential
ch["chronos_pct"] = 1.0 - ch.groupby("ensg_id")["essentiality"].rank(pct=True, method="average")
print(f"[{time.time()-t0:.1f}s] chronos_pct built for {ch.ensg_id.nunique():,} genes")

merged = unknown.merge(ch[["model_id", "ensg_id", "chronos_pct"]], on=["model_id", "ensg_id"], how="inner")
print(f"[{time.time()-t0:.1f}s] merged: {len(merged):,} rows, {merged.ensg_id.nunique():,} genes with any overlap")

# Vectorised Spearman (Pearson-on-ranks via groupby aggregation -- no per-gene python loop)
merged["r1"] = merged.groupby("ensg_id")["core_score"].rank(method="average")
merged["r2"] = merged.groupby("ensg_id")["chronos_pct"].rank(method="average")

g = merged.groupby("ensg_id")
merged["dx"]   = merged["r1"] - g["r1"].transform("mean")
merged["dy"]   = merged["r2"] - g["r2"].transform("mean")
merged["dxdy"] = merged["dx"] * merged["dy"]
merged["dx2"]  = merged["dx"] ** 2
merged["dy2"]  = merged["dy"] ** 2

agg = merged.groupby("ensg_id").agg(
    n=("r1", "size"),
    sum_dxdy=("dxdy", "sum"),
    sum_dx2=("dx2", "sum"),
    sum_dy2=("dy2", "sum"),
)
agg["rho"] = agg["sum_dxdy"] / np.sqrt(agg["sum_dx2"] * agg["sum_dy2"])
agg = agg[agg["n"] >= MIN_N].copy()

with np.errstate(invalid="ignore", divide="ignore"):
    tstat = agg["rho"] * np.sqrt((agg["n"] - 2) / (1 - agg["rho"] ** 2))
agg["pval"] = 2 * tdist.sf(np.abs(tstat), df=agg["n"] - 2)

agg["chronos_check"] = np.where(
    (agg.rho > RHO_THRESHOLD) & (agg.pval < P_THRESHOLD), "validated",
    np.where((agg.rho < -RHO_THRESHOLD) & (agg.pval < P_THRESHOLD), "inverted", "none")
)
agg = agg.reset_index()[["ensg_id", "n", "rho", "pval", "chronos_check"]]

print(f"[{time.time()-t0:.1f}s] classified {len(agg):,} genes")
print(agg["chronos_check"].value_counts().to_string())

agg.to_parquet(OUTPUTS / "chronos_validation.parquet", index=False)
print(f"\n[{time.time()-t0:.1f}s] Written: {OUTPUTS / 'chronos_validation.parquet'}  ({len(agg):,} rows)")
