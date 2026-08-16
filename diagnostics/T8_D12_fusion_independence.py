"""
T8 — D12: Noisy-OR independence assumption for fusions.
Measures correlation between p_ffpm and p_recur post-fix (absolute scales).
Run from repo root: python diagnostics/T8_D12_fusion_independence.py
"""
import sys
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import stats

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "final_pipeline" / "03_Altercations"))
from fusions_scoring import hill

print("=" * 72)
print("T8 — D12: NOISY-OR INDEPENDENCE (FUSION COMPONENTS)")
print("=" * 72)

raw = pd.read_parquet(REPO / "cleaned_track_data/fusions_gene_level.parquet")
print(f"[COUNT] fusions_gene_level rows: {len(raw):,}")

raw["conf_ord"] = raw["max_confidence"].map({"low":1,"medium":2,"high":3}).astype(float)
raw["p_conf"]  = hill(raw["conf_ord"],        p0=2.0, k=1.5)
raw["p_ffpm"]  = hill(raw["best_ffpm"],        p0=0.5, k=2.0)
raw["p_recur"] = hill(raw["fusion_count"],     p0=2.0, k=2.0)

# Drop rows missing any component
df = raw.dropna(subset=["p_conf","p_ffpm","p_recur"]).copy()
print(f"[COUNT] rows with all three components: {len(df):,}")

# Component distributions
for c in ["p_conf","p_ffpm","p_recur"]:
    print(f"[MEASURED] {c}: mean={df[c].mean():.4f}  median={df[c].median():.4f}  "
          f"std={df[c].std():.4f}  p5={df[c].quantile(0.05):.4f}  p95={df[c].quantile(0.95):.4f}")

# Correlations — all three pairs
r_fr, _ = stats.pearsonr(df["p_ffpm"], df["p_recur"])
r_fc, _ = stats.pearsonr(df["p_ffpm"], df["p_conf"])
r_rc, _ = stats.pearsonr(df["p_recur"], df["p_conf"])
print(f"\n[MEASURED] Pearson r(p_ffpm, p_recur) = {r_fr:.4f}")
print(f"[MEASURED] Pearson r(p_ffpm, p_conf)  = {r_fc:.4f}")
print(f"[MEASURED] Pearson r(p_recur, p_conf) = {r_rc:.4f}")

rs_fr, _ = stats.spearmanr(df["p_ffpm"], df["p_recur"])
rs_fc, _ = stats.spearmanr(df["p_ffpm"], df["p_conf"])
rs_rc, _ = stats.spearmanr(df["p_recur"], df["p_conf"])
print(f"[MEASURED] Spearman r(p_ffpm, p_recur) = {rs_fr:.4f}")
print(f"[MEASURED] Spearman r(p_ffpm, p_conf)  = {rs_fc:.4f}")
print(f"[MEASURED] Spearman r(p_recur, p_conf) = {rs_rc:.4f}")

# Cross-tab: high ffpm vs high recur
thr = 0.5
df["hi_ffpm"]  = df["p_ffpm"]  >= thr
df["hi_recur"] = df["p_recur"] >= thr
df["hi_conf"]  = df["p_conf"]  >= thr
ct = pd.crosstab(df["hi_ffpm"], df["hi_recur"], margins=True)
print(f"\n[MEASURED] Crosstab hi_ffpm × hi_recur (threshold={thr}):\n{ct}")

# Expected under independence
p_hf = df["hi_ffpm"].mean()
p_hr = df["hi_recur"].mean()
obs_joint = (df["hi_ffpm"] & df["hi_recur"]).mean()
exp_joint = p_hf * p_hr
print(f"\n[MEASURED] P(hi_ffpm) = {p_hf:.4f}")
print(f"[MEASURED] P(hi_recur) = {p_hr:.4f}")
print(f"[MEASURED] P(both) observed  = {obs_joint:.4f}")
print(f"[MEASURED] P(both) expected  = {exp_joint:.4f}  (independence)")
print(f"[MEASURED] Ratio obs/exp = {obs_joint/exp_joint:.4f}")

print("\n" + "=" * 72)
print("VERDICT AGAINST PRE-DECLARATION:")
print(f"  Pre-declared: |r(p_ffpm, p_recur)| < 0.50 (approximately independent)")
if abs(r_fr) < 0.50:
    print(f"  CONFIRMED — Pearson r={r_fr:.4f} (Spearman r={rs_fr:.4f})")
    print("  Noisy-OR independence assumption is approximately satisfied.")
else:
    print(f"  NOT CONFIRMED — |r|={abs(r_fr):.4f} >= 0.50")
    print(f"  WARNING: Noisy-OR inflation likely; components are correlated.")
print("=" * 72)
