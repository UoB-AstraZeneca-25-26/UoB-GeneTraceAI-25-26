"""
T6 — D23: Fusion in-frame boost dependency.
Run from repo root: python diagnostics/T6_D23_fusion_boost.py
"""
import sys
from pathlib import Path
import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(REPO / "architecture"))
from config import FUSIONS_SCR
sys.path.insert(0, str(REPO / "architecture" / "03_Altercations"))
from fusions_scoring import hill

print("=" * 72)
print("T6 — D23: FUSION IN-FRAME BOOST DEPENDENCY")
print("=" * 72)

# Verify constants in code
import fusions_scoring as fs
# Extract constants from the module
import ast, inspect
src = inspect.getsource(fs)
print("[DOCUMENTED] fusions_scoring.py constants (from source):")
print(f"  CONF_MAP      = {fs.CONF_MAP}")
print(f"  INFRAME_BOOST = {fs.INFRAME_BOOST}")

# Verify the hill function parameters actually used
# Parse from docstring/comments
print()
# Load the raw fusions data to recompute from scratch
raw = pd.read_parquet(REPO / "cleaned_track_data/fusions_gene_level.parquet")
print(f"[COUNT] fusions_gene_level: {len(raw):,} rows")
print(f"[COUNT] any_in_frame not-null: {raw.any_in_frame.notna().sum():,}")
print(f"[COUNT] any_in_frame TRUE: {raw.any_in_frame.fillna(False).sum():,}")

# Recompute with the FIXED parameters (p0_conf=2.0, p0_ffpm=0.5, p0_recur=2.0)
P0_CONF  = 2.0;  K_CONF  = 1.5
P0_FFPM  = 0.5;  K_FFPM  = 2.0
P0_RECUR = 2.0;  K_RECUR = 2.0
IN_FRAME_BOOST = 0.15
MUT_DRIVER_THRESHOLD = 0.50

# Assertions — check module source for inline hill call parameters
import re as _re
src = inspect.getsource(fs)
_p0c = _re.search(r'hill\(df\["conf_ord"\].*?p0=([\d.]+)', src)
_p0f = _re.search(r'hill\(df\["best_ffpm"\].*?p0=([\d.]+)', src)
_p0r = _re.search(r'hill\(df\["fusion_count"\].*?p0=([\d.]+)', src)
print(f"[DOCUMENTED] Source hill p0s: conf={_p0c.group(1) if _p0c else '?'}  "
      f"ffpm={_p0f.group(1) if _p0f else '?'}  recur={_p0r.group(1) if _p0r else '?'}")
assert _p0c and float(_p0c.group(1)) == 2.0, f"p0_conf != 2.0 in source: {_p0c}"
assert _p0f and float(_p0f.group(1)) == 0.5, f"p0_ffpm != 0.5 in source: {_p0f}"
assert _p0r and float(_p0r.group(1)) == 2.0, f"p0_recur != 2.0 in source: {_p0r}"
assert fs.INFRAME_BOOST == 0.15,             f"INFRAME_BOOST mismatch: got {fs.INFRAME_BOOST}"
IN_FRAME_BOOST = fs.INFRAME_BOOST
assert P0_CONF  == 2.0 and K_CONF  == 1.5,  f"Local P0_CONF/K_CONF mismatch"
assert P0_FFPM  == 0.5 and K_FFPM  == 2.0,  f"Local P0_FFPM/K_FFPM mismatch"
assert P0_RECUR == 2.0 and K_RECUR == 2.0,  f"Local P0_RECUR/K_RECUR mismatch"
print("[MEASURED] All constant assertions PASS")

df = raw.copy()
df["conf_ord"] = df["max_confidence"].map({"low":1,"medium":2,"high":3}).astype(float)
df["p_conf"]   = hill(df["conf_ord"],   p0=P0_CONF,  k=K_CONF)
df["p_ffpm"]   = hill(df["best_ffpm"],  p0=P0_FFPM,  k=K_FFPM)
df["p_recur"]  = hill(df["fusion_count"], p0=P0_RECUR, k=K_RECUR)
df["inframe"]  = df["any_in_frame"].fillna(False).astype(int)

p_base = 1 - (1-df["p_conf"].fillna(0))*(1-df["p_ffpm"].fillna(0))*(1-df["p_recur"].fillna(0))
df["p_base"]   = p_base
df["p_fusion"] = (p_base + IN_FRAME_BOOST * df["inframe"]).clip(upper=1.0)

N_total   = len(df)
N_pass    = (df["p_fusion"] >= MUT_DRIVER_THRESHOLD).sum()
N_pass_boostonly = ((df["p_base"] < MUT_DRIVER_THRESHOLD) & (df["p_fusion"] >= MUT_DRIVER_THRESHOLD)).sum()

print(f"\n[MEASURED] N_total:           {N_total:,}")
print(f"[MEASURED] N_pass (p>=0.50):  {N_pass:,}  ({100*N_pass/N_total:.2f}%)")
print(f"[MEASURED] N_pass_boostonly:  {N_pass_boostonly:,}  ({100*N_pass_boostonly/N_pass:.2f}% of passing)")
print(f"[MEASURED] fraction_boostonly: {N_pass_boostonly/N_pass:.4f}")

# Cross-check N_pass/N_total against spec's 0.648 figure
ratio = N_pass / N_total
print(f"\n[MEASURED] N_pass/N_total = {ratio:.4f}  (spec claims 0.648, diff={abs(ratio-0.648):.4f})")
if abs(ratio - 0.648) <= 0.01:
    print("  MATCHES spec within 0.01 tolerance")
else:
    print(f"  MISMATCH: spec 0.648 vs measured {ratio:.4f}  diff={abs(ratio-0.648):.4f}")

# p_base histogram in [0.30, 0.50) with in_frame breakdown
band = df[(df["p_base"] >= 0.30) & (df["p_base"] < 0.50)]
print(f"\n[MEASURED] p_base in [0.30, 0.50): {len(band):,} rows")
print(f"[MEASURED]   any_in_frame=TRUE:  {band.inframe.sum():,}  ({100*band.inframe.mean():.1f}%)")
print(f"[MEASURED]   any_in_frame=FALSE: {(band.inframe==0).sum():,}")
print(f"[MEASURED]   p_fusion >= 0.50:   {(band.p_fusion>=0.50).sum():,}  (all boosted-to-pass)")

# P(in_frame | p_base < 0.40) vs P(in_frame | p_base >= 0.70)
low  = df[df["p_base"] < 0.40]
high = df[df["p_base"] >= 0.70]
print(f"\n[MEASURED] P(in_frame | p_base < 0.40)  = {low.inframe.mean():.4f}  (n={len(low):,})")
print(f"[MEASURED] P(in_frame | p_base >= 0.70) = {high.inframe.mean():.4f}  (n={len(high):,})")

print("\n" + "=" * 72)
print("VERDICT AGAINST PRE-DECLARATION:")
print(f"  Pre-declared: fraction_boostonly > 0.05; N_pass/N_total ≈ 0.648")
if N_pass_boostonly/N_pass > 0.05:
    print(f"  CONFIRMED — fraction_boostonly = {N_pass_boostonly/N_pass:.4f} > 0.05")
    print("  Single boolean override is material.")
else:
    print(f"  NOT CONFIRMED — fraction_boostonly = {N_pass_boostonly/N_pass:.4f}")
print("=" * 72)
