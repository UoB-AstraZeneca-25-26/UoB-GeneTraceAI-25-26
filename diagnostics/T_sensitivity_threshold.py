"""
T_sensitivity_threshold.py
Checks whether the global LN_IC50 ≤ 1.0 threshold in eval.py
encodes drug identity rather than line sensitivity.

Diagnostic: per drug, count % of lines called sensitive.
If distribution clusters near 0% and 100%, threshold is drug-specific
and hit@20 inherits drug-identity signal.

Run: python diagnostics/T_sensitivity_threshold.py
"""
import sys
from pathlib import Path
import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "final_pipeline"))
from config import VALIDATION

gdsc_path = VALIDATION / "gdsc_scored_ready.parquet"
gdsc = pd.read_parquet(gdsc_path)
ic50_col = "LN_IC50" if "LN_IC50" in gdsc.columns else "log_ic50"

print("=" * 60)
print("Stage 6 sensitivity threshold diagnostic")
print("=" * 60)
print(f"Columns: {list(gdsc.columns)}")
print(f"Rows: {len(gdsc):,}  |  Drugs: {gdsc.drug_name.nunique() if 'drug_name' in gdsc.columns else '?':,}")
print(f"Global {ic50_col}: mean={gdsc[ic50_col].mean():.3f}  SD={gdsc[ic50_col].std():.3f}")
print(f"Global sensitive (≤1.0): {(gdsc[ic50_col]<=1.0).mean()*100:.1f}%")

drug_col = "drug_name" if "drug_name" in gdsc.columns else (
    "DRUG_NAME" if "DRUG_NAME" in gdsc.columns else gdsc.columns[0])

per_drug = gdsc.groupby(drug_col)[ic50_col].agg(
    n="count",
    mean="mean",
    pct_sensitive=lambda x: (x <= 1.0).mean()
)
print(f"\nPer-drug sensitive% (global threshold LN_IC50 ≤ 1.0):")
print(per_drug["pct_sensitive"].describe().to_string())

bins = [0, 0.05, 0.20, 0.80, 0.95, 1.01]
labels = ["<5%", "5-20%", "20-80%", "80-95%", ">95%"]
per_drug["bucket"] = pd.cut(per_drug["pct_sensitive"], bins=bins, labels=labels, right=False)
print(f"\nDistribution across drugs:")
print(per_drug["bucket"].value_counts().sort_index().to_string())
print("\n--- Drugs with >80% lines sensitive (nearly always sensitive) ---")
print(per_drug[per_drug["pct_sensitive"] > 0.80].sort_values("pct_sensitive", ascending=False).head(10).to_string())
print("\n--- Drugs with <5% lines sensitive (nearly never sensitive) ---")
print(per_drug[per_drug["pct_sensitive"] < 0.05].sort_values("pct_sensitive").head(10).to_string())
