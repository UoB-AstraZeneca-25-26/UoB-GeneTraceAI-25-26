"""
T_beta_verify.py — quick diagnostic: does β = ρ × SD(prot)/SD(rna) move
r(rna_pct_z, prot_resid_z) from −0.100 toward 0?

Runs entirely from the existing parquet files; does NOT rerun core_score.py.
Simulates the corrected residual on a sample of genes to confirm direction.

Pre-declared outcome: corrected r should be within ±0.03 of 0.
"""
import sys
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import stats

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "final_pipeline"))
from config import RNA_Z, PROT_Z, CORE_SCORE

RHO_PRIOR_OLD = 0.373
RHO_PRIOR_NEW = 0.353
K_RHO = 30
N_MIN_RHO = 30

print("=" * 70)
print("T_beta_verify — β correction: r(rna_pct_z, prot_resid_z) before/after")
print("=" * 70)

rna  = pd.read_parquet(RNA_Z).rename(columns={"gene_id": "ensg_id", "z_t": "rna_z"})
prot = pd.read_parquet(PROT_Z).rename(columns={"gene_id": "ensg_id", "z_t": "prot_z"})
core = pd.read_parquet(CORE_SCORE, columns=["model_id","ensg_id","core_score"])

genes_both = set(rna.ensg_id) & set(prot.ensg_id)
print(f"[COUNT] Genes in both RNA + prot: {len(genes_both):,}")

merged = prot[prot.ensg_id.isin(genes_both)].merge(
    rna[rna.ensg_id.isin(genes_both)][["ensg_id","model_id","rna_z"]],
    on=["ensg_id","model_id"], how="inner"
)
print(f"[COUNT] Merged rows: {len(merged):,}")

# Build per-gene rho and beta (OLD and NEW) ─────────────────────────────────
old_rows, new_rows = [], []
for g, grp in merged.groupby("ensg_id"):
    if len(grp) < N_MIN_RHO:
        continue
    rho_raw, _ = stats.pearsonr(grp["prot_z"].fillna(0), grp["rna_z"].fillna(0))
    n = len(grp)
    lam = K_RHO / (n + K_RHO)

    rho_old = lam * RHO_PRIOR_OLD + (1 - lam) * rho_raw
    rho_new = lam * RHO_PRIOR_NEW + (1 - lam) * rho_raw

    sd_rna  = grp["rna_z"].std(ddof=1)
    sd_prot = grp["prot_z"].std(ddof=1)
    beta_new = rho_new * (sd_prot / sd_rna) if sd_rna > 0 else rho_new

    old_rows.append({"ensg_id": g, "beta": rho_old})   # OLD: beta = rho (wrong)
    new_rows.append({"ensg_id": g, "beta": beta_new})   # NEW: beta = rho * SD ratio

old_df = pd.DataFrame(old_rows)
new_df = pd.DataFrame(new_rows)

print(f"[COUNT] Genes with per-gene β: {len(old_df):,}")
print(f"[MEASURED] OLD β: median={old_df.beta.median():.4f}  mean={old_df.beta.mean():.4f}")
print(f"[MEASURED] NEW β: median={new_df.beta.median():.4f}  mean={new_df.beta.mean():.4f}")
print(f"[MEASURED] Overcorrection factor OLD/NEW: "
      f"{old_df.beta.median()/new_df.beta.median():.3f}×  (expect ≈1.88)")

# Compute residuals with OLD and NEW beta ─────────────────────────────────────
merged = merged.merge(old_df.rename(columns={"beta":"beta_old"}), on="ensg_id", how="left")
merged = merged.merge(new_df.rename(columns={"beta":"beta_new"}), on="ensg_id", how="left")
# fallback for genes < N_MIN_RHO
sd_rna_global  = merged["rna_z"].std(ddof=1)
sd_prot_global = merged["prot_z"].std(ddof=1)
fallback_old = RHO_PRIOR_OLD
fallback_new = RHO_PRIOR_NEW * (sd_prot_global / sd_rna_global)
merged["beta_old"] = merged["beta_old"].fillna(fallback_old)
merged["beta_new"] = merged["beta_new"].fillna(fallback_new)

merged["resid_old"] = merged["prot_z"] - merged["beta_old"] * merged["rna_z"]
merged["resid_new"] = merged["prot_z"] - merged["beta_new"] * merged["rna_z"]

# ── Join rna_pct_z from core_score (it's not directly output, so use probit of lineage rank)
# Simpler: just measure r(rna_z, resid) directly — same signed direction as r(rna_pct_z, resid)
r_old, p_old = stats.pearsonr(merged["rna_z"].fillna(0), merged["resid_old"].fillna(0))
r_new, p_new = stats.pearsonr(merged["rna_z"].fillna(0), merged["resid_new"].fillna(0))

print()
print(f"[MEASURED] r(rna_z, resid_OLD) = {r_old:+.4f}  p={p_old:.2e}  "
      f"  (expect negative: overcorrected)")
print(f"[MEASURED] r(rna_z, resid_NEW) = {r_new:+.4f}  p={p_new:.2e}  "
      f"  (pre-declared: within ±0.03 of 0)")

# ── Analytical check: what should r(rna_z, resid) be? ─────────────────────
# cov(Y - β·X, X) = cov(Y,X) - β·var(X)
# For OLD: β = rho_old ≈ 0.373 → cov = r·SD_prot·SD_rna - 0.373·SD_rna^2
#   = 0.318×0.951×1.525 - 0.373×1.525^2 = 0.461 - 0.868 = -0.407  (negative)
# For NEW: β = rho_new × SD_prot/SD_rna → cov = r·SD_prot·SD_rna - β·SD_rna^2
#   = r·SD_prot·SD_rna - rho_new·SD_prot/SD_rna·SD_rna^2 = 0  (exactly zero)
sd_prot_m = merged["prot_z"].std()
sd_rna_m  = merged["rna_z"].std()
cov_pr    = r_old * sd_prot_m * sd_rna_m   # approx from pooled data
print()
print("[ANALYTICAL] Expected r(rna_z, resid) under OLD β:")
cov_old_expected = cov_pr - old_df.beta.median() * sd_rna_m**2
r_old_pred = cov_old_expected / (sd_rna_m * merged["resid_old"].std())
print(f"  cov(Y,X)={cov_pr:.4f}, β_old·var(X)={old_df.beta.median()*sd_rna_m**2:.4f} "
      f"→ predicted r≈{r_old_pred:.4f}")

print()
if abs(r_new) <= 0.03:
    print("VERDICT: PRE-DECLARED OUTCOME MET — r(rna_z, resid_NEW) within ±0.03 of 0")
else:
    print(f"VERDICT: PRE-DECLARED OUTCOME NOT MET — r = {r_new:+.4f} "
          f"(tolerance ±0.03; SD-ratio explanation may be incomplete)")
print("=" * 70)
