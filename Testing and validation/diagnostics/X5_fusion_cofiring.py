"""
X5 — Fusion channel co-firing at the gate.
COWORK v4.0. No changes to p0. Read-only.
Run from repo root: python diagnostics/X5_fusion_cofiring.py
"""
import sys
from pathlib import Path
import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(REPO / "architecture"))

FUSIONS_SRC = REPO / "cleaned_track_data" / "fusions_gene_level.parquet"

P0_CONF  = 2.0; K_CONF  = 1.5
P0_FFPM  = 0.5; K_FFPM  = 2.0
P0_RECUR = 2.0; K_RECUR = 2.0
INFRAME_BOOST = 0.15
FUS_GATE      = 0.5

print("=" * 72)
print("X5 -- FUSION CHANNEL CO-FIRING AT THE GATE (COWORK v4.0)")
print("=" * 72)

def hill(x, p0, k):
    xk = np.power(np.clip(x, 0, None), k)
    return xk / (xk + p0 ** k)

CONF_MAP = {"low": 1, "medium": 2, "high": 3}

df = pd.read_parquet(FUSIONS_SRC)
print(f"[COUNT] fusions_gene_level rows: {len(df):,}  "
      f"genes: {df.ensg_id.nunique():,}  models: {df.model_id.nunique():,}")

df = df.copy()
df["conf_ord"] = df["max_confidence"].map(CONF_MAP).astype(float)
df["p_conf"]   = hill(df["conf_ord"],      P0_CONF,  K_CONF)
df["p_ffpm"]   = hill(df["best_ffpm"],     P0_FFPM,  K_FFPM)
df["p_recur"]  = hill(df["fusion_count"],  P0_RECUR, K_RECUR)

df["p_base"]  = 1.0 - (1.0 - df["p_conf"].fillna(0)) * \
                       (1.0 - df["p_ffpm"].fillna(0)) * \
                       (1.0 - df["p_recur"].fillna(0))
inframe = df["any_in_frame"].fillna(False).astype(int)
df["p_fusion"] = (df["p_base"] + INFRAME_BOOST * inframe).clip(upper=1.0)

passing = df[df["p_fusion"] >= FUS_GATE].copy()
n_pass  = len(passing)
n_total = len(df)
print(f"\n[COUNT] Fusions passing gate (p_fusion >= {FUS_GATE}): {n_pass:,} / {n_total:,} "
      f"({100*n_pass/n_total:.1f}%)")

# Channel analysis: which channels >= 0.4 among passing
THRESH = 0.4
passing["c_conf"]  = (passing["p_conf"]  >= THRESH).astype(int)
passing["c_ffpm"]  = (passing["p_ffpm"]  >= THRESH).astype(int)
passing["c_recur"] = (passing["p_recur"] >= THRESH).astype(int)
passing["n_chan"]  = passing["c_conf"] + passing["c_ffpm"] + passing["c_recur"]

print(f"\n--- CHANNEL COUNT >= {THRESH} AMONG PASSING FUSIONS ---")
for n, grp in passing.groupby("n_chan"):
    print(f"  {n} channels >= {THRESH}: {len(grp):,}  ({100*len(grp)/n_pass:.1f}%)")

# Which single channel dominates for 1-channel passes
one_chan = passing[passing["n_chan"] == 1].copy()
if len(one_chan) > 0:
    pass_conf  = (one_chan["c_conf"] == 1).sum()
    pass_ffpm  = (one_chan["c_ffpm"] == 1).sum()
    pass_recur = (one_chan["c_recur"]== 1).sum()
    print(f"\n--- SINGLE-CHANNEL PASSES (n={len(one_chan):,}) ---")
    print(f"  via p_conf alone:   {pass_conf:,}  ({100*pass_conf/len(one_chan):.1f}%)")
    print(f"  via p_ffpm alone:   {pass_ffpm:,}  ({100*pass_ffpm/len(one_chan):.1f}%)")
    print(f"  via p_recur alone:  {pass_recur:,}  ({100*pass_recur/len(one_chan):.1f}%)")

# Genuine Noisy-OR: no single channel >= 0.5 but p_fusion >= 0.5
genuine_nor = passing[(passing["p_conf"] < 0.5) & (passing["p_ffpm"] < 0.5) & (passing["p_recur"] < 0.5)]
# And not just from inframe boost
purely_boost = passing[(passing["p_base"] < 0.5) & (passing["p_fusion"] >= 0.5)]
print(f"\n--- GENUINE COMBINATION (no single channel >= 0.5) ---")
print(f"  {len(genuine_nor):,} / {n_pass:,}  ({100*len(genuine_nor)/n_pass:.1f}%)")

print(f"\n--- INFRAME BOOST ONLY (p_base<0.5, p_fusion>=0.5) ---")
print(f"  {len(purely_boost):,} / {n_pass:,}  ({100*len(purely_boost)/n_pass:.1f}%)  "
      f"(pre-declared: ~8%)")

# Channel distribution stats
print(f"\n--- CHANNEL SCORE DISTRIBUTION (all fusions) ---")
for ch, col in [("p_conf","p_conf"), ("p_ffpm","p_ffpm"), ("p_recur","p_recur")]:
    s = df[col].dropna()
    print(f"  {ch}: mean={s.mean():.3f}  median={s.median():.3f}  "
          f">=0.5: {(s>=0.5).mean()*100:.1f}%")

print(f"\n--- PRE-DECLARED CHECKS ---")
pct_1chan = 100 * len(one_chan) / n_pass
pct_genuine = 100 * len(genuine_nor) / n_pass
print(f"  pre-declared: >60% passes have exactly 1 channel >= 0.4")
print(f"  measured: {pct_1chan:.1f}%  {'CONFIRMED' if pct_1chan > 60 else 'WRONG'}")
print(f"  pre-declared: <20% passes through genuine combination (no channel >= 0.5)")
print(f"  measured: {pct_genuine:.1f}%  {'CONFIRMED' if pct_genuine < 20 else 'WRONG'}")

print("\n" + "=" * 72)
