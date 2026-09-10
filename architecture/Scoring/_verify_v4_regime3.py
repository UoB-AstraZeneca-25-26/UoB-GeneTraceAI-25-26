import pandas as pd
import numpy as np

LIVE = "../outputs/core_score.parquet"
V4   = "../outputs/core_score_v4_regime3.parquet"

live = pd.read_parquet(LIVE)
v4   = pd.read_parquet(V4)

print("="*70)
print("STRUCTURAL CHECK")
print("="*70)
print(f"live: rows={len(live):,} genes={live.ensg_id.nunique():,} lines={live.model_id.nunique():,}")
print(f"v4:   rows={len(v4):,} genes={v4.ensg_id.nunique():,} lines={v4.model_id.nunique():,}")
assert len(live) == len(v4), "ROW COUNT MISMATCH"
assert live.ensg_id.nunique() == v4.ensg_id.nunique()
assert live.model_id.nunique() == v4.model_id.nunique()
print("Row/gene/line counts match.")

m = live[["ensg_id","model_id","core_score","n_layers"]].merge(
    v4[["ensg_id","model_id","core_score","n_layers","combination_variant"]],
    on=["ensg_id","model_id"], suffixes=("_live","_v4"), how="outer", indicator=True)
print(f"\nMerge indicator: {m._merge.value_counts().to_dict()}")
assert (m._merge == "both").all(), "ROWS DON'T ALIGN 1:1"
assert (m.n_layers_live == m.n_layers_v4).all(), "n_layers MISMATCH"

non_regime3 = m[m.combination_variant == "A"]
regime3     = m[m.combination_variant == "C_regime3"]
print(f"\nnon-regime3 rows: {len(non_regime3):,}  |  regime3 rows: {len(regime3):,}")

diff_nonregime3 = (non_regime3.core_score_live - non_regime3.core_score_v4).abs()
print(f"Max |diff| on non-regime3 rows: {diff_nonregime3.max():.2e}  "
      f"(bit-identity check: {'PASS' if diff_nonregime3.max() < 1e-9 else 'FAIL'})")

diff_regime3 = (regime3.core_score_live - regime3.core_score_v4).abs()
print(f"Regime3 rows: mean|diff|={diff_regime3.mean():.4f}  (expected to differ -- this IS the fix)")

print("\n" + "="*70)
print("SPIKE COMPARISON")
print("="*70)
def spike(df, col):
    s = df[col].dropna()
    return ((s<0.01)|(s>0.99)).mean()*100

n2_live = live[live.n_layers==2]
n2_v4   = v4[v4.n_layers==2]
print(f"live (A):  n_layers=2 spike={spike(n2_live,'core_score'):.4f}%   full-panel spike={spike(live,'core_score'):.4f}%")
print(f"v4 (A+regime3): n_layers=2 spike={spike(n2_v4,'core_score'):.4f}%   full-panel spike={spike(v4,'core_score'):.4f}%")

print("\n" + "="*70)
print("REGIME-3 WORKED EXAMPLE: ENSG00000164985 / ach-000679")
print("="*70)
for label, df in [("live (A)", live), ("v4 (A+regime3)", v4)]:
    row = df[(df.ensg_id=="ENSG00000164985") & (df.model_id.str.lower()=="ach-000679")]
    if len(row):
        r = row.iloc[0]
        extra = f"  variant={r.combination_variant}" if "combination_variant" in row.columns else ""
        print(f"{label}: core_score={r.core_score:.4f}  stratum_rank={r.stratum_rank}{extra}")
    else:
        print(f"{label}: NOT FOUND")
