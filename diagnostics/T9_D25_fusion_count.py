"""
T9 — D25: fusion_count semantics verification.
Does fusion_count vary within a gene (across lines) or is it constant per gene?
Run from repo root: python diagnostics/T9_D25_fusion_count.py
"""
import sys
from pathlib import Path
import pandas as pd
import numpy as np

REPO = Path(__file__).resolve().parent.parent

print("=" * 72)
print("T9 — D25: FUSION_COUNT SEMANTICS")
print("=" * 72)

raw = pd.read_parquet(REPO / "cleaned_track_data/fusions_gene_level.parquet")
print(f"[COUNT] fusions_gene_level rows: {len(raw):,}")
print(f"[COUNT] columns: {list(raw.columns)}")

# Is fusion_count constant per gene (ensg_id) across model_ids?
per_gene = raw.groupby("ensg_id")["fusion_count"].agg(["min","max","nunique"])
n_gene = len(per_gene)
n_constant = (per_gene["nunique"] == 1).sum()
n_varying  = (per_gene["nunique"] > 1).sum()
print(f"\n[COUNT] Distinct genes: {n_gene:,}")
print(f"[COUNT] Genes where fusion_count is constant across lines: {n_constant:,}  ({100*n_constant/n_gene:.1f}%)")
print(f"[COUNT] Genes where fusion_count varies across lines:      {n_varying:,}  ({100*n_varying/n_gene:.1f}%)")

# Distribution of fusion_count
vc = raw["fusion_count"].value_counts().sort_index().head(20)
print(f"\n[MEASURED] fusion_count value distribution (top 20 values):")
for k,v in vc.items():
    print(f"  {k}: {v:,}  ({100*v/len(raw):.2f}%)")

print(f"\n[MEASURED] fusion_count: mean={raw.fusion_count.mean():.2f}  "
      f"median={raw.fusion_count.median():.0f}  "
      f"p95={raw.fusion_count.quantile(0.95):.0f}  "
      f"max={raw.fusion_count.max():.0f}")

# What does fusion_count measure? Likely across-all-lines count
# vs best_ffpm (per-line value)
print(f"\n[MEASURED] best_ffpm: mean={raw.best_ffpm.mean():.4f}  "
      f"median={raw.best_ffpm.median():.4f}  "
      f"max={raw.best_ffpm.max():.4f}")

# Correlation between best_ffpm and fusion_count
corr = raw[["best_ffpm","fusion_count"]].corr().iloc[0,1]
print(f"[MEASURED] Pearson r(best_ffpm, fusion_count) = {corr:.4f}")

# Hill function impact
sys.path.insert(0, str(REPO / "final_pipeline" / "03_Altercations"))
from fusions_scoring import hill as fhill
raw["p_recur"] = fhill(raw["fusion_count"], p0=2.0, k=2.0)
print(f"\n[MEASURED] p_recur (hill, p0=2.0, k=2.0):")
for fc in [1,2,3,5,10,20,50]:
    print(f"  fusion_count={fc}: p_recur={fhill(fc, 2.0, 2.0):.4f}")

# Singletons (fusion_count==1)
n_sing = (raw["fusion_count"]==1).sum()
print(f"\n[MEASURED] Singletons (fusion_count=1): {n_sing:,}  ({100*n_sing/len(raw):.2f}%)")
print(f"[MEASURED] p_recur for singleton: {fhill(1, 2.0, 2.0):.4f}")

print("\n" + "=" * 72)
print("VERDICT AGAINST PRE-DECLARATION:")
print("  Pre-declared: fusion_count is constant within gene (across all lines)")
if n_constant / n_gene > 0.98:
    print(f"  CONFIRMED — {n_constant}/{n_gene} genes constant ({100*n_constant/n_gene:.1f}%)")
    print("  fusion_count = genome-wide recurrence count (not per-line).")
    print("  Implication: hill(fusion_count) gives same recurrence weight to ALL lines.")
elif n_varying / n_gene > 0.10:
    print(f"  NOT CONFIRMED — {n_varying}/{n_gene} genes vary ({100*n_varying/n_gene:.1f}%)")
    print("  fusion_count varies per line — it IS line-specific.")
else:
    print(f"  AMBIGUOUS — {100*n_constant/n_gene:.1f}% constant")
print("=" * 72)
