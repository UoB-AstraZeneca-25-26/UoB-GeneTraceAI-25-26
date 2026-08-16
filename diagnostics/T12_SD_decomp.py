"""
T12 supplement — SD(prot_z) and SD(rna_z) directly — item 12.
Explains why SD(core_score | L2) = 1.13 instead of the predicted 0.95.
Run from repo root: python diagnostics/T12_SD_decomp.py
"""
import sys
from pathlib import Path
import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "final_pipeline"))
from config import RNA_Z, PROT_Z, CORE_SCORE

print("=" * 72)
print("T12 SUPPLEMENT — SD DECOMPOSITION FOR n_layers VARIANCE")
print("=" * 72)

core = pd.read_parquet(CORE_SCORE)
rna  = pd.read_parquet(RNA_Z)
prot = pd.read_parquet(PROT_Z)

print(f"[COUNT] core: {len(core):,}  rna: {len(rna):,}  prot: {len(prot):,}")
print(f"[DOCUMENTED] core cols: {list(core.columns)}")
print(f"[DOCUMENTED] rna  cols: {list(rna.columns)}")
print(f"[DOCUMENTED] prot cols: {list(prot.columns)}")

# SD of core_score by n_layers
for nl in [1, 2]:
    sub = core[core["n_layers"] == nl]
    print(f"\n[MEASURED] n_layers={nl}: SD(core_score) = {sub.core_score.std():.6f}  "
          f"mean={sub.core_score.mean():.6f}  n={len(sub):,}")

# Rename to common key
rna_r  = rna.rename(columns={"gene_id":"ensg_id","z_t":"rna_z"})[["ensg_id","model_id","rna_z"]]
prot_r = prot.rename(columns={"gene_id":"ensg_id","z_t":"prot_z"})[["ensg_id","model_id","prot_z"]]

# SD of RNA z-score (marginal)
print(f"\n[MEASURED] RNA z (bulk_rna_z.z_t): SD = {rna_r.rna_z.std():.6f}  "
      f"mean={rna_r.rna_z.mean():.6f}")
print(f"  p5={rna_r.rna_z.quantile(0.05):.4f}  p95={rna_r.rna_z.quantile(0.95):.4f}")

# SD of protein z-score (marginal)
print(f"\n[MEASURED] Prot z (bulk_prot_z.z_t): SD = {prot_r.prot_z.std():.6f}  "
      f"mean={prot_r.prot_z.mean():.6f}")
print(f"  p5={prot_r.prot_z.quantile(0.05):.4f}  p95={prot_r.prot_z.quantile(0.95):.4f}")

# Join core (L2) with RNA and protein for paired SD analysis
print("\n[STATUS] Joining L2 core_score pairs with rna_z and prot_z...")
l2 = core[core["n_layers"]==2][["ensg_id","model_id","core_score"]].copy()
l1 = core[core["n_layers"]==1][["ensg_id","model_id","core_score"]].copy()

j2 = l2.merge(rna_r,  on=["ensg_id","model_id"], how="inner") \
        .merge(prot_r, on=["ensg_id","model_id"], how="inner") \
        .dropna()
print(f"[COUNT] L2 joined triplets (core, rna, prot): {len(j2):,}")

if len(j2) > 100:
    sd_r  = j2.rna_z.std()
    sd_p  = j2.prot_z.std()
    sd_cs = j2.core_score.std()
    corr_rp = j2[["rna_z","prot_z"]].corr().iloc[0,1]

    print(f"\n[MEASURED] L2 SD components:")
    print(f"  SD(rna_z)           = {sd_r:.6f}")
    print(f"  SD(prot_z)          = {sd_p:.6f}")
    print(f"  r(rna_z, prot_z)    = {corr_rp:.6f}")
    print(f"  SD(core_score | L2) = {sd_cs:.6f}")

    # Expected SD from Stouffer combination
    # core_score = (W_RNA * rna_z + W_PROT * prot_z) / sqrt(W_RNA^2 + W_PROT^2)
    # SD(core_score) = sqrt((W_RNA*sd_r)^2 + (W_PROT*sd_p)^2 + 2*W_RNA*W_PROT*r*sd_r*sd_p)
    #                  / sqrt(W_RNA^2 + W_PROT^2)
    W_RNA  = np.sqrt(2.00)
    W_PROT = np.sqrt(1.45)
    W_NORM = np.sqrt(W_RNA**2 + W_PROT**2)
    sd_expected = np.sqrt(
        (W_RNA*sd_r)**2 + (W_PROT*sd_p)**2 + 2*W_RNA*W_PROT*corr_rp*sd_r*sd_p
    ) / W_NORM
    print(f"\n[MEASURED] Stouffer expected SD(core_score|L2) = {sd_expected:.6f}")
    print(f"[MEASURED] Observed  SD(core_score|L2)          = {sd_cs:.6f}")
    print(f"[MEASURED] Ratio observed/expected              = {sd_cs/sd_expected:.4f}")

    # Compare L1
    j1 = l1.merge(rna_r, on=["ensg_id","model_id"], how="inner").dropna()
    sd_r1  = j1.rna_z.std()
    sd_cs1 = j1.core_score.std()
    print(f"\n[MEASURED] L1 SD components:")
    print(f"  SD(rna_z)           = {sd_r1:.6f}")
    print(f"  SD(core_score | L1) = {sd_cs1:.6f}")
    print(f"[MEASURED] SD ratio observed L2/L1 = {sd_cs/sd_cs1:.6f}")
    print(f"[MEASURED] SD ratio expected L2/L1 ~ 0.95 (pre-declared)")

    print(f"\n[DOCUMENTED] WHY SD(L2) > SD(L1):")
    if corr_rp > 0:
        print(f"  r(rna_z, prot_z) = {corr_rp:.4f} (POSITIVE)")
        print(f"  SD(prot_z) = {sd_p:.4f}  vs  SD(rna_z) = {sd_r:.4f}")
        if sd_p > sd_r:
            print(f"  Protein z-scores have LARGER spread than RNA z-scores.")
            print(f"  With positive correlation, Stouffer AMPLIFIES variance")
            print(f"  instead of cancelling it. Pre-declared model assumed")
            print(f"  SD(prot) ≈ SD(rna) with ~zero correlation → shrinkage.")
        else:
            print(f"  Despite positive correlation, SD(prot) < SD(rna).")
            print(f"  The Stouffer numerator adds positively correlated terms,")
            print(f"  preventing the variance reduction expected at r=0.")
    else:
        print(f"  r(rna_z, prot_z) = {corr_rp:.4f} (NEGATIVE or zero)")
        print(f"  Negative correlation would cause shrinkage. Discrepancy")
        print(f"  is likely explained by differences in scale/calibration.")

print("\n" + "=" * 72)
print("SUMMARY:")
print("  The pre-declared SD(L2)/SD(L1) ~ 0.95 assumed:")
print("    (1) SD(prot_z) ~ SD(rna_z)")
print("    (2) r(prot_z, rna_z) ~ 0 across cell lines")
print("  Measured values above show which assumption is violated.")
print("=" * 72)
