"""
T11 supplement — corr(n_layers, study_frequency) — item 11 from remediation list.
Tests whether the n_layers=2 fame bias (L2 lines are 1.19× more likely top-20%)
is explained by well-studied/frequently-published cancer lines having both
protein coverage AND inherent oncogenic biology.
Run from repo root: python diagnostics/T11_study_freq.py
"""
import sys
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import stats
import duckdb

REPO = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(REPO / "architecture"))
from config import PREDICTIONS

print("=" * 72)
print("T11 SUPPLEMENT — corr(n_layers, study_frequency)")
print("=" * 72)

pred = pd.read_parquet(PREDICTIONS, columns=["model_id","n_layers","score_rank_pct"])
pred["top20"] = pred["score_rank_pct"] >= 0.80

# Per-cell-line: protein coverage flag
line_has_prot = pred.groupby("model_id")["n_layers"].max()  # 2 if any protein, else 1
print(f"[COUNT] Cell lines with protein (max_n_layers=2): "
      f"{(line_has_prot==2).sum():,} / {len(line_has_prot):,}")

# Study frequency proxy: mutational_burden from gdsc_models
# gdsc_models.model_id is already lowercase ach- format
DB = REPO / "architecture" / "outputs" / "celllineselector.db"
con = duckdb.connect(str(DB), read_only=True)

try:
    burden = con.execute(
        "SELECT lower(model_id) as model_id, mutational_burden "
        "FROM main.gdsc_models WHERE mutational_burden IS NOT NULL"
    ).df()
    print(f"[COUNT] gdsc_models with mutational_burden: {len(burden):,}")
    proxy_name = "mutational_burden (GDSC)"
except Exception as e:
    print(f"  gdsc_models failed: {e}")
    burden = pd.DataFrame(columns=["model_id","mutational_burden"])
    proxy_name = "none"

con.close()

if len(burden) == 0:
    print("[ERROR] No study proxy available.")
    sys.exit(1)

burden = burden.rename(columns={"mutational_burden":"study_proxy"})

# Merge
line_df = line_has_prot.reset_index()
line_df.columns = ["model_id","n_layers_max"]
line_df["has_protein"] = (line_df["n_layers_max"] == 2).astype(int)

merged = line_df.merge(burden, on="model_id", how="inner").dropna()
print(f"[COUNT] Cell lines in merged analysis: {len(merged):,}")

if len(merged) < 10:
    print("[ERROR] Insufficient overlap.")
    sys.exit(1)

r_pt, p_pt = stats.pointbiserialr(merged["has_protein"], merged["study_proxy"])
r_sp, p_sp = stats.spearmanr(merged["has_protein"], merged["study_proxy"])

print(f"\n[MEASURED] Point-biserial r(has_protein, {proxy_name}): {r_pt:.4f}  p={p_pt:.4e}")
print(f"[MEASURED] Spearman r(has_protein, {proxy_name}):        {r_sp:.4f}  p={p_sp:.4e}")

prot_lines    = merged.loc[merged["has_protein"]==1, "study_proxy"]
nonprot_lines = merged.loc[merged["has_protein"]==0, "study_proxy"]
print(f"\n[MEASURED] Protein lines    (n={len(prot_lines):,}):  "
      f"median={prot_lines.median():.2f}  mean={prot_lines.mean():.2f}  "
      f"IQR=[{prot_lines.quantile(.25):.2f},{prot_lines.quantile(.75):.2f}]")
print(f"[MEASURED] Non-protein lines (n={len(nonprot_lines):,}): "
      f"median={nonprot_lines.median():.2f}  mean={nonprot_lines.mean():.2f}  "
      f"IQR=[{nonprot_lines.quantile(.25):.2f},{nonprot_lines.quantile(.75):.2f}]")

u_stat, u_p = stats.mannwhitneyu(prot_lines, nonprot_lines, alternative="greater")
print(f"\n[MEASURED] Mann-Whitney U (protein > non-protein): p={u_p:.4e}")

# Also: corr(n_layers_max, study_proxy) — continuous version
r_cont, p_cont = stats.spearmanr(merged["n_layers_max"], merged["study_proxy"])
print(f"[MEASURED] Spearman r(n_layers_max, study_proxy) [continuous]: {r_cont:.4f}  p={p_cont:.4e}")

print("\n" + "=" * 72)
print("VERDICT:")
if abs(r_sp) > 0.20 and p_sp < 0.05:
    print(f"  CONFIRMED — protein coverage correlates with mutation burden")
    print(f"  (Spearman r={r_sp:.4f}, p={p_sp:.2e})")
    print("  The n_layers fame bias is at least partially explained by")
    print("  well-studied (high-burden) lines being proteomically profiled.")
elif p_sp < 0.05:
    print(f"  WEAK SIGNAL — Spearman r={r_sp:.4f} significant but small")
else:
    print(f"  NOT CONFIRMED — Spearman r={r_sp:.4f} (p={p_sp:.2e})")
    print("  The fame bias is not explained by mutational burden alone.")
print("=" * 72)
