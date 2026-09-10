"""
X3 — Stage 5 selectivity: gene-pair independence audit.
COWORK v4.0. No changes. Read-only.
Run from repo root: python diagnostics/X3_selectivity_corr.py
"""
import sys
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import stats

REPO = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(REPO / "architecture"))
from config import PREDICTIONS

N_PAIRS = 1000
SEED    = 99

print("=" * 72)
print("X3 -- SELECTIVITY GENE-PAIR INDEPENDENCE (COWORK v4.0)")
print("=" * 72)

pred = pd.read_parquet(PREDICTIONS, columns=["ensg_id","model_id","core_score","n_layers"])
print(f"[COUNT] predictions rows: {len(pred):,}  genes: {pred.ensg_id.nunique():,}  "
      f"models: {pred.model_id.nunique():,}")

all_genes = pred.ensg_id.unique()
rng = np.random.default_rng(SEED)
idx = rng.choice(len(all_genes), size=N_PAIRS * 2, replace=False)
gene_list = all_genes[idx]
pairs = [(gene_list[i], gene_list[i + N_PAIRS]) for i in range(N_PAIRS)]
print(f"[COUNT] Gene pairs sampled: {len(pairs):,}  (seed={SEED})")

pearson_r  = []
spearman_r = []
n_shared   = []

# Stratification bins
corr_L2_only   = []
corr_L1_only   = []
corr_same_gene = []

for gene_a, gene_b in pairs:
    a = pred[pred.ensg_id == gene_a][["model_id","core_score","n_layers"]].set_index("model_id")
    b = pred[pred.ensg_id == gene_b][["model_id","core_score","n_layers"]].set_index("model_id")
    merged = a.join(b, how="inner", lsuffix="_a", rsuffix="_b").dropna()
    n = len(merged)
    if n < 10:
        continue
    n_shared.append(n)
    pr = merged["core_score_a"].corr(merged["core_score_b"])
    sr, _ = stats.spearmanr(merged["core_score_a"], merged["core_score_b"])
    pearson_r.append(pr)
    spearman_r.append(sr)
    # Stratify by n_layers
    l2 = merged[(merged.n_layers_a == 2) & (merged.n_layers_b == 2)]
    l1 = merged[(merged.n_layers_a == 1) & (merged.n_layers_b == 1)]
    if len(l2) >= 10:
        corr_L2_only.append(l2["core_score_a"].corr(l2["core_score_b"]))
    if len(l1) >= 10:
        corr_L1_only.append(l1["core_score_a"].corr(l1["core_score_b"]))

pearson_r  = np.array(pearson_r)
spearman_r = np.array(spearman_r)
n_shared   = np.array(n_shared)

print(f"\n[COUNT] Eligible pairs (>=10 shared lines): {len(pearson_r):,}")
print(f"\n--- PEARSON r(score_A, score_B) ---")
print(f"  mean={pearson_r.mean():.4f}  median={np.median(pearson_r):.4f}  SD={pearson_r.std():.4f}")
print(f"  IQR=[{np.percentile(pearson_r,25):.4f}, {np.percentile(pearson_r,75):.4f}]")
print(f"  p5={np.percentile(pearson_r,5):.4f}  p95={np.percentile(pearson_r,95):.4f}")
print(f"  >0.1: {(pearson_r>0.1).mean()*100:.1f}%  >0.2: {(pearson_r>0.2).mean()*100:.1f}%  "
      f">0.3: {(pearson_r>0.3).mean()*100:.1f}%")

print(f"\n--- SPEARMAN r(score_A, score_B) ---")
print(f"  mean={spearman_r.mean():.4f}  median={np.median(spearman_r):.4f}  SD={spearman_r.std():.4f}")
print(f"  IQR=[{np.percentile(spearman_r,25):.4f}, {np.percentile(spearman_r,75):.4f}]")
print(f"  >0.1: {(spearman_r>0.1).mean()*100:.1f}%  >0.2: {(spearman_r>0.2).mean()*100:.1f}%  "
      f">0.3: {(spearman_r>0.3).mean()*100:.1f}%")

print(f"\n--- N_SHARED LINES per pair ---")
print(f"  mean={n_shared.mean():.0f}  median={np.median(n_shared):.0f}  "
      f"min={n_shared.min()}  max={n_shared.max()}")

print(f"\n--- STRATIFICATION BY n_layers ---")
if corr_L2_only:
    arr = np.array(corr_L2_only)
    print(f"  L2-only lines: n_pairs={len(arr)}  "
          f"median r = {np.median(arr):.4f}  mean = {arr.mean():.4f}")
if corr_L1_only:
    arr = np.array(corr_L1_only)
    print(f"  L1-only lines: n_pairs={len(arr)}  "
          f"median r = {np.median(arr):.4f}  mean = {arr.mean():.4f}")

# Rank displacement in top-20 selectivity
print(f"\n--- RANK DISPLACEMENT: selectivity score vs rank-based ---")
# Pick 100 pairs where both genes have many lines
sub_pairs = [(a,b) for a,b,n in zip(
    [pairs[i][0] for i in range(len(pearson_r))],
    [pairs[i][1] for i in range(len(pearson_r))],
    n_shared) if n >= 50][:100]

displacements = []
for gene_a, gene_b in sub_pairs:
    a = pred[pred.ensg_id == gene_a][["model_id","core_score"]].set_index("model_id")
    b = pred[pred.ensg_id == gene_b][["model_id","core_score"]].set_index("model_id")
    merged = a.join(b, how="inner", lsuffix="_a", rsuffix="_b").dropna()
    if len(merged) < 20:
        continue
    merged["sel_raw"]   = merged["core_score_a"] * (1 - merged["core_score_b"])
    # Correlation-adjusted: centre each score then take product
    merged["a_c"] = merged["core_score_a"] - merged["core_score_a"].mean()
    merged["b_c"] = merged["core_score_b"] - merged["core_score_b"].mean()
    merged["sel_adj"]   = merged["a_c"] * (1 - merged["core_score_b"])
    top_raw = set(merged.nlargest(20, "sel_raw").index)
    top_adj = set(merged.nlargest(20, "sel_adj").index)
    overlap = len(top_raw & top_adj)
    displacements.append(20 - overlap)

if displacements:
    d = np.array(displacements)
    print(f"  n_pairs evaluated: {len(d)}")
    print(f"  mean rank displacement in top-20: {d.mean():.2f}  (0=identical, 20=complete swap)")
    print(f"  median: {np.median(d):.1f}  max: {d.max()}")

print(f"\n--- PRE-DECLARED CHECK ---")
print(f"  pre-declared: median corr(score_A, score_B) > 0.2")
med_pearson = float(np.median(pearson_r))
print(f"  measured median Pearson: {med_pearson:.4f}")
if med_pearson > 0.2:
    print("  verdict: PRE-DECLARED CONFIRMED")
else:
    print("  verdict: PRE-DECLARED WRONG -- correlation < 0.2")

print("\n" + "=" * 72)
