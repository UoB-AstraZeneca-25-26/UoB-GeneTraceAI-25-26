"""
W4_rna_pctz_lineage.py
Compute mean(rna_pct_z) per lineage × n_layers and check whether the
lineages over-represented in L2 also have higher rna_pct_z.

Compositional contribution = sum_lineage (share_L2 - share_L1) * mean_rna_pct_z_lineage
Compare against T4's measured mean component (1.20 pp).

Run: python diagnostics/W4_rna_pctz_lineage.py
"""
import sys
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import stats
import duckdb

REPO = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(REPO / "architecture"))
from config import RNA_Z, CORE_SCORE, DB

# T4 mean component (from session summary / T4_decomp.py run)
T4_MEAN_COMPONENT_PP = 1.20  # pp, measured

print("=" * 70)
print("W4 — mean(rna_pct_z) per lineage × n_layers")
print("=" * 70)

# Load RNA z-scores and lineage
rna = pd.read_parquet(RNA_Z).rename(columns={"z_t": "rna_z", "gene_id": "ensg_id"})
print(f"RNA rows: {len(rna):,}  genes: {rna.ensg_id.nunique():,}  models: {rna.model_id.nunique():,}")

con = duckdb.connect(str(DB), read_only=True)
lin = con.execute(
    "SELECT model_id, lineage FROM main.sample_info WHERE lineage IS NOT NULL"
).df()
con.close()
lin = lin.drop_duplicates("model_id")
lineage_map = lin.set_index("model_id")["lineage"]
print(f"Lineage map: {len(lineage_map)} models, {lineage_map.nunique()} lineages")

# Compute rna_pct_z (identical logic to core_score.py)
rna_wide = rna.pivot(index="model_id", columns="ensg_id", values="rna_z")
rna_wide = rna_wide.loc[rna_wide.index.isin(lin["model_id"])]

# Percentile within lineage, per gene
rna_wide_with_lin = rna_wide.copy()
rna_wide_with_lin["lineage"] = lineage_map.reindex(rna_wide_with_lin.index)

pct = rna_wide_with_lin.groupby("lineage", group_keys=False).apply(
    lambda g: g.rank(pct=True, ascending=True, na_option="keep"),
    include_groups=False,
)
# pct has same shape as rna_wide (numeric columns only)
rna_pct_z_wide = pd.DataFrame(
    stats.norm.ppf(np.clip(pct.values, 0.001, 0.999)),
    index=pct.index, columns=pct.columns,
)
print(f"rna_pct_z computed: shape {rna_pct_z_wide.shape}")

# n_layers per model from CORE_SCORE
cs = pd.read_parquet(CORE_SCORE, columns=["model_id", "ensg_id", "n_layers"])
model_layers = cs.groupby("model_id")["n_layers"].max().reset_index()
model_layers = model_layers.merge(lin, on="model_id", how="left")
model_layers = model_layers[model_layers["lineage"].notna()]

print(f"L1 models: {(model_layers['n_layers']==1).sum()}  "
      f"L2 models: {(model_layers['n_layers']==2).sum()}")

# Mean rna_pct_z per (lineage × n_layers)
# Melt rna_pct_z_wide to long format
rna_pct_z_long = rna_pct_z_wide.stack().reset_index()
rna_pct_z_long.columns = ["model_id", "ensg_id", "rna_pct_z"]
rna_pct_z_long = rna_pct_z_long.merge(model_layers[["model_id","lineage","n_layers"]],
                                        on="model_id", how="left")
rna_pct_z_long = rna_pct_z_long[rna_pct_z_long["lineage"].notna()]

# Lineage × n_layers mean
lineage_mean = rna_pct_z_long.groupby(["lineage", "n_layers"])["rna_pct_z"].agg(
    mean="mean", sd="std", n="count"
).reset_index()
lineage_mean_wide = lineage_mean.pivot_table(
    index="lineage", columns="n_layers", values="mean"
)
lineage_mean_wide.columns = [f"L{int(c)}_mean_rna_pct_z" for c in lineage_mean_wide.columns]

# Lineage × n_layers n_models (for composition weights)
lin_n = rna_pct_z_long.drop_duplicates(["model_id","lineage","n_layers"]).groupby(
    ["lineage","n_layers"]
)["model_id"].count().reset_index()
lin_n_wide = lin_n.pivot_table(index="lineage", columns="n_layers", values="model_id").fillna(0)
lin_n_wide.columns = [f"L{int(c)}_n" for c in lin_n_wide.columns]
# Share
l1_n = lin_n_wide.get("L1_n", 0)
l2_n = lin_n_wide.get("L2_n", 0)
total_l1 = l1_n.sum()
total_l2 = l2_n.sum()
lin_n_wide["L1_share"] = l1_n / max(total_l1, 1)
lin_n_wide["L2_share"] = l2_n / max(total_l2, 1)
lin_n_wide["diff_pp"]  = (lin_n_wide["L2_share"] - lin_n_wide["L1_share"]) * 100

combined = lineage_mean_wide.join(lin_n_wide, how="outer")
combined = combined.sort_values("diff_pp", ascending=False)

print(f"\nmean(rna_pct_z) per lineage × n_layers (sorted by L2 over-representation):")
print(combined.round(4).to_string())

# Compositional contribution
if "L1_mean_rna_pct_z" in combined.columns and "L2_mean_rna_pct_z" in combined.columns:
    # Use mean rna_pct_z per lineage (averaging L1 and L2 as global lineage trait)
    combined["mean_rna_pct_z_lineage"] = combined[
        ["L1_mean_rna_pct_z","L2_mean_rna_pct_z"]
    ].mean(axis=1)
elif "L2_mean_rna_pct_z" in combined.columns:
    combined["mean_rna_pct_z_lineage"] = combined["L2_mean_rna_pct_z"]
else:
    combined["mean_rna_pct_z_lineage"] = np.nan

# Compositional contribution = Σ (share_L2 - share_L1) * mean_rna_pct_z
# This is in rna_pct_z units; convert to core_score pp using W_RNA / W_NORM weight
combined_valid = combined.dropna(subset=["mean_rna_pct_z_lineage"])
comp_contrib_z = (combined_valid["diff_pp"] / 100 * combined_valid["mean_rna_pct_z_lineage"]).sum()

# core_score is Phi((W_RNA * rna_pct_z + W_PROT * prot_resid_z) / W_NORM)
# Approximate pp contribution: if we shift rna_pct_z by delta at the mean (norm.cdf'(0)=0.399)
# delta_core_score ≈ 0.399 * W_RNA/W_NORM * delta_rna_pct_z
W_RNA  = np.sqrt(2.00)
W_PROT = np.sqrt(1.45)
W_NORM = np.sqrt(W_RNA**2 + W_PROT**2)
rna_share = W_RNA / W_NORM
phi_prime_at_0 = stats.norm.pdf(0)  # ≈ 0.399

comp_contrib_pp = comp_contrib_z * rna_share * phi_prime_at_0 * 100  # in pp

print(f"\nCompositional contribution:")
print(f"  Σ (share_L2 - share_L1) × mean_rna_pct_z = {comp_contrib_z:.4f} rna_pct_z units")
print(f"  RNA weight in core_z: W_RNA/W_NORM = {rna_share:.4f}")
print(f"  ΔΦ per unit z at mean: {phi_prime_at_0:.4f}")
print(f"  ≈ compositional contribution to L2-L1 gap: {comp_contrib_pp:.2f} pp")
print(f"  T4 mean component (measured): {T4_MEAN_COMPONENT_PP:.2f} pp")
if T4_MEAN_COMPONENT_PP > 0:
    frac = comp_contrib_pp / T4_MEAN_COMPONENT_PP
    print(f"  Fraction explained by composition: {frac*100:.1f}%  "
          f"(pre-declared ≥50%)")
    print(f"  PRE-DECLARED: {'PASS' if frac >= 0.5 else 'FAIL — composition accounts for <50%'}")

# V7 link: are the 5 most over-represented L2 lineages above the global mean?
global_mean_z = rna_pct_z_long["rna_pct_z"].mean()
print(f"\nGlobal mean(rna_pct_z): {global_mean_z:.4f}")
top5_over = combined.nlargest(5, "diff_pp")
print(f"\nTop 5 lineages over-represented in L2 vs global mean:")
for lin_name, row in top5_over.iterrows():
    lin_mean = row.get("mean_rna_pct_z_lineage", float("nan"))
    above = "ABOVE" if lin_mean > global_mean_z else "BELOW"
    print(f"  {lin_name:20s}  diff={row.get('diff_pp',0):+.2f}pp  "
          f"mean_rna_pct_z={lin_mean:.4f}  {above} global mean")

print(f"\nWell-studied lineages check (lung, breast, colorectal):")
for lname in ["lung", "breast", "large_intestine"]:
    if lname in combined.index:
        r = combined.loc[lname]
        print(f"  {lname}: diff_pp={r.get('diff_pp',float('nan')):+.2f}  "
              f"mean_rna_pct_z={r.get('mean_rna_pct_z_lineage',float('nan')):.4f}")
