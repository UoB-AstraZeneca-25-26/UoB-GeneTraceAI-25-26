"""
T4_decomp.py — Three checks after T4 post-fix widening.

Check 1: SD(prot_resid_z) and SD(rna_pct_z) on L2 subset.
          Predicted: ~1.43 and ~1.00.
          If MAD underestimates spread (non-normal), "z-scores" come out with SD > 1.

Check 2: Mean vs variance decomposition of the -4.84 pp gap.
          Counterfactual A: give L2 L1's SD, keep L2's mean → P(top20%|L2)?
          Counterfactual B: give L2 L1's mean, keep L2's SD → P(top20%|L2)?
          Splits structural (selection) vs fixable (standardisation artefact).

Check 3: How many pairs changed mut_driver status with P0_VEP 3.0→1.5?
          Is p_vep contributing anything to the driver call?

Run from repo root: python diagnostics/T4_decomp.py
"""
import sys
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import stats

REPO = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(REPO / "architecture"))
from config import RNA_Z, PROT_Z, CORE_SCORE, PREDICTIONS

W_RNA  = np.sqrt(2.00)
W_PROT = np.sqrt(1.45)
W_NORM = np.sqrt(W_RNA**2 + W_PROT**2)
RHO_PRIOR = 0.353
K_RHO = 30
N_MIN_RHO = 30

print("=" * 70)
print("T4_decomp — SD decomposition, mean/variance split, p_vep audit")
print("=" * 70)

# ── Load predictions (has n_layers, score_rank_pct, core_score) ───────────
core = pd.read_parquet(PREDICTIONS, columns=["model_id","ensg_id","core_score","n_layers","score_rank_pct"])
core["core_z"] = stats.norm.ppf(np.clip(core["core_score"], 1e-9, 1-1e-9))
core["top20"]  = core["score_rank_pct"] >= 0.80

l1 = core[core["n_layers"] == 1]
l2 = core[core["n_layers"] == 2]

mean_z1 = l1["core_z"].mean();  sd_z1 = l1["core_z"].std(ddof=1)
mean_z2 = l2["core_z"].mean();  sd_z2 = l2["core_z"].std(ddof=1)
p_top_l1 = l1["top20"].mean()
p_top_l2 = l2["top20"].mean()

print(f"\n[CHECK 2] core_z statistics by n_layers:")
print(f"  L1: mean={mean_z1:.4f}  SD={sd_z1:.4f}  P(top20%)={p_top_l1:.6f}  n={len(l1):,}")
print(f"  L2: mean={mean_z2:.4f}  SD={sd_z2:.4f}  P(top20%)={p_top_l2:.6f}  n={len(l2):,}")
print(f"  Gap: {(p_top_l2-p_top_l1)*100:+.4f} pp  mean diff={mean_z2-mean_z1:+.4f}  SD ratio={sd_z2/sd_z1:.4f}")

# Counterfactuals — applied within-gene (vectorized via groupby rank).
# Rescale L2 core_z to have L1's SD/mean, then compute new within-gene rank.
# CF_A: L2 gets L1's SD, keeps L2's mean.  CF_B: L2 gets L1's mean, keeps L2's SD.
l2_idx   = core["n_layers"] == 2
z2_unit  = (core.loc[l2_idx, "core_z"] - mean_z2) / sd_z2  # mean=0, SD=1

cfA_col = z2_unit * sd_z1 + mean_z2   # L2 mean, L1 SD
cfB_col = z2_unit * sd_z2 + mean_z1   # L1 mean, L2 SD

cfA_df = core[["ensg_id","core_z","n_layers"]].copy()
cfB_df = core[["ensg_id","core_z","n_layers"]].copy()
cfA_df.loc[l2_idx, "core_z"] = cfA_col.values
cfB_df.loc[l2_idx, "core_z"] = cfB_col.values

def within_gene_pct80(df):
    df = df.copy()
    df["rk"] = df.groupby("ensg_id")["core_z"].rank(pct=True, ascending=True)
    return df

cfA_df = within_gene_pct80(cfA_df)
cfB_df = within_gene_pct80(cfB_df)

p_cfA = (cfA_df.loc[l2_idx, "rk"] >= 0.80).mean()
p_cfB = (cfB_df.loc[l2_idx, "rk"] >= 0.80).mean()

print(f"\n[CHECK 2] Counterfactual decomposition of {(p_top_l2-p_top_l1)*100:+.2f} pp gap:")
print(f"  P(top20% | L2_observed)         = {p_top_l2:.6f}")
print(f"  P(top20% | L2, L1_SD, L2_mean) = {p_cfA:.6f}  delta from obs: {(p_cfA-p_top_l2)*100:+.2f} pp  ← SD effect")
print(f"  P(top20% | L2, L1_mean, L2_SD) = {p_cfB:.6f}  delta from obs: {(p_cfB-p_top_l2)*100:+.2f} pp  ← mean effect")
print(f"  P(top20% | L1_observed)         = {p_top_l1:.6f}")
frac_sd   = abs(p_cfA - p_top_l2) / abs(p_top_l2 - p_top_l1) if abs(p_top_l2-p_top_l1) > 0 else np.nan
frac_mean = abs(p_cfB - p_top_l2) / abs(p_top_l2 - p_top_l1) if abs(p_top_l2-p_top_l1) > 0 else np.nan
print(f"  SD accounts for  {100*frac_sd:.1f}% of the gap")
print(f"  Mean accounts for {100*frac_mean:.1f}% of the gap")

# ── Check 1: recompute prot_resid_z and rna_pct_z directly ────────────────
print(f"\n[CHECK 1] Recomputing prot_resid_z and rna_pct_z for L2 subset ...")

rna  = pd.read_parquet(RNA_Z).rename(columns={"gene_id":"ensg_id","z_t":"rna_z"})
prot = pd.read_parquet(PROT_Z).rename(columns={"gene_id":"ensg_id","z_t":"prot_z"})

genes_both = set(rna.ensg_id) & set(prot.ensg_id)
merged = prot[prot.ensg_id.isin(genes_both)].merge(
    rna[rna.ensg_id.isin(genes_both)][["ensg_id","model_id","rna_z"]],
    on=["ensg_id","model_id"], how="inner"
)

# Per-gene beta (same logic as fixed core_score.py)
rows = []
for g, grp in merged.groupby("ensg_id"):
    if len(grp) < N_MIN_RHO:
        continue
    rho_raw, _ = stats.pearsonr(grp["prot_z"].fillna(0), grp["rna_z"].fillna(0))
    n = len(grp); lam = K_RHO / (n + K_RHO)
    rho_g = lam * RHO_PRIOR + (1 - lam) * rho_raw
    sd_rna  = grp["rna_z"].std(ddof=1)
    sd_prot = grp["prot_z"].std(ddof=1)
    beta_g  = rho_g * (sd_prot / sd_rna) if sd_rna > 0 else rho_g
    rows.append({"ensg_id": g, "beta_g": beta_g})
rho_df = pd.DataFrame(rows)

sd_rna_g  = merged["rna_z"].std(ddof=1)
sd_prot_g = merged["prot_z"].std(ddof=1)
fallback  = RHO_PRIOR * (sd_prot_g / sd_rna_g)

merged = merged.merge(rho_df, on="ensg_id", how="left")
merged["beta_g"]     = merged["beta_g"].fillna(fallback)
merged["prot_resid"] = merged["prot_z"] - merged["beta_g"] * merged["rna_z"]

# Standardise per gene: (x − median) / SD — matches core_score.py post-fix.
prot_resid_wide = merged.pivot(index="model_id", columns="ensg_id", values="prot_resid")
_prot_med = prot_resid_wide.median(axis=0)
_prot_sd  = prot_resid_wide.std(ddof=1, axis=0).clip(lower=1e-6)
prot_resid_z_wide = prot_resid_wide.subtract(_prot_med, axis=1).div(_prot_sd, axis=1)

# rna_pct_z: need lineage-conditioned percentile → probit
import duckdb, sys as _sys
from config import DB
con = duckdb.connect(str(DB), read_only=True)
lin = con.execute("SELECT model_id, lineage FROM main.sample_info WHERE lineage IS NOT NULL").df()
lin = lin.drop_duplicates("model_id")
con.close()
lin_map = lin.set_index("model_id")["lineage"]

rna_wide = rna.pivot(index="model_id", columns="ensg_id", values="rna_z")
rna_wide["lineage"] = lin_map.reindex(rna_wide.index)
rna_pct = rna_wide.groupby("lineage", group_keys=False).apply(
    lambda g: g.rank(pct=True, ascending=True, na_option="keep"), include_groups=False
)
rna_pct_z_wide = pd.DataFrame(
    stats.norm.ppf(np.clip(rna_pct.values, 0.001, 0.999)),
    index=rna_pct.index, columns=rna_pct.columns
)

# L2 models are those with protein data
l2_models = list(set(prot_resid_z_wide.index) & set(rna_pct_z_wide.index))
shared_genes = list(set(prot_resid_z_wide.columns) & set(rna_pct_z_wide.columns))

prz_l2 = prot_resid_z_wide.loc[l2_models, shared_genes].values.ravel()
rpz_l2 = rna_pct_z_wide.loc[l2_models, shared_genes].values.ravel()

prz_l2 = prz_l2[np.isfinite(prz_l2)]
rpz_l2 = rpz_l2[np.isfinite(rpz_l2)]

sd_prz = np.std(prz_l2, ddof=1)
sd_rpz = np.std(rpz_l2, ddof=1)

print(f"  SD(prot_resid_z | L2) = {sd_prz:.4f}  (predicted ≈1.43; 1.0 = perfectly standardised)")
print(f"  SD(rna_pct_z    | L2) = {sd_rpz:.4f}  (expected ≈1.00; probit of ranks ≈ unit variance)")

# Implied Var(core_z | L2) from components
r_cross, _ = stats.pearsonr(prz_l2[:len(rpz_l2)], rpz_l2[:len(prz_l2)])
var_core_pred = (W_RNA**2 * sd_rpz**2 + W_PROT**2 * sd_prz**2 + 2*W_RNA*W_PROT*r_cross*sd_rpz*sd_prz) / W_NORM**2
print(f"  r(rna_pct_z, prot_resid_z | L2) = {r_cross:.4f}")
print(f"  Predicted Var(core_z|L2) from components = {var_core_pred:.4f}  (measured SD²={sd_z2**2:.4f})")
prot_var_share = (W_PROT**2 * sd_prz**2) / (W_NORM**2 * var_core_pred)
rna_var_share  = (W_RNA**2  * sd_rpz**2) / (W_NORM**2 * var_core_pred)
print(f"  RNA contribution to Var(core_z):  {100*rna_var_share:.1f}%")
print(f"  Prot contribution to Var(core_z): {100*prot_var_share:.1f}%  (expected ≈42% if both SD=1)")

# ── Check 3: p_vep contribution to mut_driver ─────────────────────────────
print(f"\n[CHECK 3] p_vep contribution to driver call ...")

# mutations_collapsed.parquet has the raw inputs (max_vep_rank, max_pathogenicity, variant_count)
mut = pd.read_parquet(REPO / "cleaned_track_data" / "mutations_collapsed.parquet")
print(f"  mutations_collapsed rows: {len(mut):,}")

def hill(x, p0, k):
    return x**k / (x**k + p0**k)

# Reproduce components (matching mutations_scoring.py)
vb  = mut["variant_count"].clip(upper=5).astype(float)
pth = mut["max_pathogenicity"].fillna(0).astype(float)
vrk = mut["max_vep_rank"].fillna(0).astype(float)

p_path   = hill(pth,  0.5, 2.0)
p_burden = hill(vb,   3.0, 1.5)
p_vep_old = hill(vrk, 3.0, 2.0)   # original P0_VEP
p_vep_cur = hill(vrk, 2.5, 2.0)   # current P0_VEP (rank2=0.390<gate, rank3=0.590>gate)
p_vep_no  = pd.Series(np.zeros(len(vrk)))

p_base_old    = 1 - (1-p_vep_old) * (1-p_path) * (1-p_burden)
p_base_cur    = 1 - (1-p_vep_cur) * (1-p_path) * (1-p_burden)
p_base_no_vep = 1 - (1-p_path)    * (1-p_burden)

# Pass rate at variant level — the guard statistic
pass_rate_cur = (p_vep_cur > 0.5).mean()

driver_old    = (p_base_old > 0.5).sum()
driver_cur    = (p_base_cur > 0.5).sum()
driver_no_vep = (p_base_no_vep > 0.5).sum()
changed       = ((p_base_old > 0.5) != (p_base_cur > 0.5)).sum()
decisive_vep  = driver_cur - driver_no_vep

saturated = (p_base_cur >= 0.999).sum()
print(f"  p_vep pass rate at variant level (P0=2.5):   {100*pass_rate_cur:.2f}%  (guard: must be <30%)")
print(f"  Pairs with p_base > 0.5 (P0_VEP=3.0 orig): {driver_old:,}")
print(f"  Pairs with p_base > 0.5 (P0_VEP=2.5 curr): {driver_cur:,}")
print(f"  Pairs with p_base > 0.5 (no p_vep):         {driver_no_vep:,}")
print(f"  Pairs that crossed 0.5 threshold old→curr:  {changed:,}  ({100*changed/len(mut):.3f}%)")
print(f"  Pairs where p_vep is decisive (curr):        {decisive_vep:,}  ({100*decisive_vep/len(mut):.3f}%)")
print(f"  Pairs saturated at p_base ≥ 0.999:           {saturated:,}  ({100*saturated/len(mut):.1f}%)")

print("\n" + "=" * 70)
print("SUMMARY:")
print(f"  SD(prot_resid_z | L2) = {sd_prz:.4f}  vs predicted 1.43")
print(f"  SD(rna_pct_z    | L2) = {sd_rpz:.4f}  vs expected  1.00")
print(f"  Prot arm variance share: {100*prot_var_share:.1f}%  (expected 42%)")
print(f"  Gap decomposition: {100*frac_sd:.1f}% from SD, {100*frac_mean:.1f}% from mean")
print("=" * 70)
