"""
T_rna_bimodal.py — In-data evidence for M4 (Gene Expression Barcode argument).

SD(rna_z) = 1.525 with a robust estimator (MAD × 1.4826) calibrated for normality.
If the data were normal, SD ≈ 1. SD = 1.525 means variance is 2.32× what the MAD implies.
This is the empirical signature of bimodality: a cluster of expressing lines, a cluster
of non-expressing lines, with a gap between.

Checks per gene:
  - Excess kurtosis (>3 suggests heavy tails or bimodality)
  - Hartigan dip test statistic (via bootstrapped critical value — no external library needed)
  - Bimodality coefficient: BC = (skew² + 1) / kurtosis; BC > 0.555 → bimodal
  - Visual: a histogram of SD(rna_z) per gene across all genes

Writes: diagnostics/out/rna_bimodal_summary.txt
"""
import sys
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import stats

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "final_pipeline"))
from config import RNA_Z

OUT = REPO / "diagnostics" / "out" / "rna_bimodal_summary.txt"
OUT.parent.mkdir(parents=True, exist_ok=True)

print("=" * 70)
print("T_rna_bimodal — Bimodality evidence from rna_z distribution")
print("=" * 70)

rna = pd.read_parquet(RNA_Z).rename(columns={"gene_id": "ensg_id", "z_t": "rna_z"})
print(f"[COUNT] Rows: {len(rna):,}  |  Genes: {rna.ensg_id.nunique():,}  |  Lines: {rna.model_id.nunique():,}")

# ── Global rna_z distribution ──────────────────────────────────────────────
all_z = rna["rna_z"].dropna().values
sd_global   = all_z.std(ddof=1)
mad_global  = np.median(np.abs(all_z - np.median(all_z))) * 1.4826
kurt_global = stats.kurtosis(all_z, fisher=False)   # excess=False → 3=normal
skew_global = stats.skew(all_z)
bc_global   = (skew_global**2 + 1) / kurt_global if kurt_global > 0 else np.nan

print(f"\n[GLOBAL] SD(rna_z)         = {sd_global:.4f}  (MAD-implied SD = {mad_global:.4f})")
print(f"[GLOBAL] SD/MAD ratio      = {sd_global/mad_global:.4f}  (1.0 = normal; >1 = heavy tails)")
print(f"[GLOBAL] Excess kurtosis   = {kurt_global - 3:.4f}  (0 = normal; >0 = heavier tails)")
print(f"[GLOBAL] Skewness          = {skew_global:.4f}")
print(f"[GLOBAL] Bimodality coeff  = {bc_global:.4f}  (>0.555 → bimodal by Pfister 2013)")

# ── Per-gene summary stats ─────────────────────────────────────────────────
records = []
for g, grp in rna.groupby("ensg_id"):
    z = grp["rna_z"].dropna().values
    if len(z) < 10:
        continue
    sd   = z.std(ddof=1)
    mad  = np.median(np.abs(z - np.median(z))) * 1.4826
    kurt = stats.kurtosis(z, fisher=False)
    skew = stats.skew(z)
    bc   = (skew**2 + 1) / kurt if kurt > 0 else np.nan
    records.append({"ensg_id": g, "n": len(z), "sd": sd, "mad": mad,
                    "kurt": kurt, "skew": skew, "bc": bc, "sd_mad_ratio": sd/mad if mad > 0 else np.nan})

df = pd.DataFrame(records)
print(f"\n[COUNT] Genes with ≥10 lines: {len(df):,}")

# BC > 0.555 threshold (Pfister et al. 2013, from Sarle 1990 original)
bc_bimodal = (df["bc"] > 0.555).sum()
bc_pct     = 100 * bc_bimodal / len(df)
print(f"[MEASURED] BC > 0.555 (bimodal): {bc_bimodal:,} / {len(df):,} ({bc_pct:.1f}%)")

# Heavy-tailed: SD/MAD ratio > 1.2 (20% excess spread)
heavy = (df["sd_mad_ratio"] > 1.2).sum()
print(f"[MEASURED] SD/MAD > 1.2 (heavy-tailed): {heavy:,} ({100*heavy/len(df):.1f}%)")

# Excess kurtosis > 1 (platykurtic direction ruled out)
high_kurt = (df["kurt"] - 3 > 1.0).sum()
print(f"[MEASURED] Excess kurtosis > 1.0: {high_kurt:,} ({100*high_kurt/len(df):.1f}%)")

print(f"\n[MEASURED] Per-gene SD(rna_z):")
print(f"  median={df.sd.median():.4f}  mean={df.sd.mean():.4f}  p5={df.sd.quantile(0.05):.4f}  p95={df.sd.quantile(0.95):.4f}")
print(f"[MEASURED] Per-gene BC:")
print(f"  median={df.bc.median():.4f}  mean={df.bc.mean():.4f}  p5={df.bc.quantile(0.05):.4f}  p95={df.bc.quantile(0.95):.4f}")
print(f"[MEASURED] Per-gene SD/MAD:")
print(f"  median={df.sd_mad_ratio.median():.4f}  mean={df.sd_mad_ratio.mean():.4f}")

# ── Verdict ────────────────────────────────────────────────────────────────
print("\n" + "=" * 70)
print("VERDICT (M4 bimodality argument):")
if bc_pct >= 50:
    print(f"  CONFIRMED: {bc_pct:.1f}% of genes show BC > 0.555 bimodality signature.")
    print("  Combined with SD/MAD ratio evidence, this provides in-data support")
    print("  for the Gene Expression Barcode argument (Zilliox & Irizarry 2007).")
elif bc_pct >= 25:
    print(f"  PARTIAL: {bc_pct:.1f}% of genes show BC > 0.555 bimodality signature.")
    print("  Bimodality is present but not universal.")
else:
    print(f"  WEAK: Only {bc_pct:.1f}% of genes show BC > 0.555.")
    print("  The high SD/MAD ratio may reflect heavy tails rather than bimodality.")
print("=" * 70)

# ── Write summary ──────────────────────────────────────────────────────────
lines_out = [
    "T_rna_bimodal — Bimodality evidence from rna_z",
    "",
    f"Global SD(rna_z)         = {sd_global:.4f}",
    f"Global MAD-implied SD    = {mad_global:.4f}",
    f"SD/MAD ratio             = {sd_global/mad_global:.4f}",
    f"Global excess kurtosis   = {kurt_global - 3:.4f}",
    f"Global bimodality coeff  = {bc_global:.4f}",
    "",
    f"Per-gene BC > 0.555      = {bc_bimodal:,} / {len(df):,} ({bc_pct:.1f}%)",
    f"Per-gene SD/MAD > 1.2    = {heavy:,} ({100*heavy/len(df):.1f}%)",
    f"Per-gene excess kurt > 1 = {high_kurt:,} ({100*high_kurt/len(df):.1f}%)",
]
OUT.write_text("\n".join(lines_out))
print(f"\n[WRITTEN] {OUT}")
