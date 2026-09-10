"""
X4 — Stage 5 co-selection: tie inflation from min().
COWORK v4.0. No changes. Read-only.
Run from repo root: python diagnostics/X4_coselect_ties.py
"""
import sys
from pathlib import Path
import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(REPO / "architecture"))
from config import PREDICTIONS

N_PAIRS  = 100
SEED     = 99
FLOOR    = 0.50

print("=" * 72)
print("X4 -- CO-SELECTION TIE INFLATION FROM min() (COWORK v4.0)")
print("=" * 72)

pred = pd.read_parquet(PREDICTIONS, columns=["ensg_id","model_id","core_score"])
print(f"[COUNT] predictions rows: {len(pred):,}  genes: {pred.ensg_id.nunique():,}  "
      f"models: {pred.model_id.nunique():,}")

all_genes = pred.ensg_id.unique()
rng = np.random.default_rng(SEED)
idx  = rng.choice(len(all_genes), size=N_PAIRS * 2, replace=False)
pairs = [(all_genes[idx[i]], all_genes[idx[i+N_PAIRS]]) for i in range(N_PAIRS)]

single_tie_rates   = []
joint_tie_rates    = []
lines_above_floor  = []
spread_nontop      = []
tied_in_top20      = []

for gene_a, gene_b in pairs:
    a = pred[pred.ensg_id == gene_a][["model_id","core_score"]].rename(columns={"core_score":"sa"})
    b = pred[pred.ensg_id == gene_b][["model_id","core_score"]].rename(columns={"core_score":"sb"})
    merged = a.merge(b, on="model_id", how="inner").dropna()
    if len(merged) < 10:
        continue

    merged["joint"] = merged[["sa","sb"]].min(axis=1)
    above = merged[merged["joint"] >= FLOOR]
    n_above = len(above)
    lines_above_floor.append(n_above)
    if n_above == 0:
        continue

    # Distinct joint scores
    n_distinct_joint = above["joint"].nunique()
    joint_tie_rate   = 1 - n_distinct_joint / n_above
    joint_tie_rates.append(joint_tie_rate)

    # Single-gene tie rate (for gene_a)
    n_above_a = (merged["sa"] >= FLOOR).sum()
    n_dist_a  = merged.loc[merged["sa"] >= FLOOR, "sa"].nunique()
    single_tie_rate = 1 - n_dist_a / max(n_above_a, 1)
    single_tie_rates.append(single_tie_rate)

    # Top-20 analysis
    if n_above >= 20:
        top20 = above.nlargest(20, "joint")
        n_tied_in_top20 = (top20["joint"] == top20["joint"].min()).sum()
        tied_in_top20.append(n_tied_in_top20)
        # Spread of non-minimum score among tied lines
        tied = top20[top20["joint"] == top20["joint"].min()]
        if len(tied) > 1:
            # Spread of the OTHER score (not the minimum)
            # For each line: the non-min score is max(sa,sb)
            non_min = tied.apply(lambda r: max(r["sa"],r["sb"]), axis=1)
            spread_nontop.append(non_min.max() - non_min.min())

n_eligible = len(joint_tie_rates)
print(f"[COUNT] Eligible pairs (n_above_floor >= 1): {n_eligible}")

if n_eligible > 0:
    jtr = np.array(joint_tie_rates)
    str_ = np.array(single_tie_rates)
    laf = np.array(lines_above_floor)

    print(f"\n--- LINES ABOVE FLOOR (joint >= {FLOOR}) ---")
    print(f"  mean={laf.mean():.1f}  median={np.median(laf):.1f}  "
          f"min={laf.min()}  max={laf.max()}")

    print(f"\n--- TIE RATES ---")
    print(f"  Joint min() tie rate:   mean={jtr.mean()*100:.1f}%  "
          f"median={np.median(jtr)*100:.1f}%  "
          f"IQR=[{np.percentile(jtr,25)*100:.1f}%, {np.percentile(jtr,75)*100:.1f}%]")
    print(f"  Single-gene tie rate:   mean={str_.mean()*100:.1f}%  "
          f"median={np.median(str_)*100:.1f}%")
    diff = jtr - str_
    print(f"  EXCESS tie rate (joint - single):  "
          f"mean={diff.mean()*100:+.1f} pp  median={np.median(diff)*100:+.1f} pp")

    if tied_in_top20:
        t20 = np.array(tied_in_top20)
        print(f"\n--- TOP-20 TIE ANALYSIS ---")
        print(f"  n_pairs with >= 20 above-floor lines: {len(t20)}")
        print(f"  Lines sharing min(joint) in top-20:  "
              f"mean={t20.mean():.1f}  median={np.median(t20):.1f}  max={t20.max()}")

    if spread_nontop:
        sp = np.array(spread_nontop)
        print(f"\n--- NON-MINIMUM SCORE SPREAD AMONG TOP-20 TIES ---")
        print(f"  (large spread = substantively different lines scoring identically)")
        print(f"  mean={sp.mean():.4f}  median={np.median(sp):.4f}  max={sp.max():.4f}")

    print(f"\n--- PRE-DECLARED CHECK ---")
    print(f"  pre-declared: joint tie rate > single-gene tie rate by >=10 pp")
    excess = float(np.median(diff)) * 100
    print(f"  measured excess (median): {excess:+.1f} pp")
    if excess >= 10:
        print("  verdict: PRE-DECLARED CONFIRMED")
    else:
        print(f"  verdict: PRE-DECLARED WRONG -- excess only {excess:.1f} pp")

print("\n" + "=" * 72)
