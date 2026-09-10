"""
X3 — Stage 5 selectivity: gene-pair independence audit (FAST version).
COWORK v4.0. No changes. Read-only.
Pivots the predictions matrix first then correlates gene columns directly.
Run from repo root: python diagnostics/X3_selectivity_corr_fast.py
"""
import sys, time
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import stats

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "final_pipeline"))
from config import PREDICTIONS

N_PAIRS = 1000
SEED    = 99

print("=" * 72)
print("X3 -- SELECTIVITY GENE-PAIR INDEPENDENCE (FAST, COWORK v4.0)")
print("=" * 72)

t0 = time.time()
pred = pd.read_parquet(PREDICTIONS, columns=["ensg_id","model_id","core_score","n_layers"])
print(f"[COUNT] predictions rows: {len(pred):,}  genes: {pred.ensg_id.nunique():,}  "
      f"models: {pred.model_id.nunique():,}  [{time.time()-t0:.0f}s]")

# Sample genes
all_genes = pred.ensg_id.unique()
rng = np.random.default_rng(SEED)
idx = rng.choice(len(all_genes), size=N_PAIRS * 2, replace=False)
gene_list = all_genes[idx]
sampled_genes = set(gene_list)
print(f"[COUNT] Gene pairs sampled: {N_PAIRS:,}  (seed={SEED})")

# Restrict to sampled genes and pivot wide
print(f"[STATUS] Pivoting to wide matrix...  [{time.time()-t0:.0f}s]")
sub = pred[pred.ensg_id.isin(sampled_genes)].copy()
wide = sub.pivot_table(index="model_id", columns="ensg_id", values="core_score", aggfunc="first")
# n_layers pivot for stratification
nl_wide = sub.pivot_table(index="model_id", columns="ensg_id", values="n_layers", aggfunc="first")
print(f"[STATUS] Wide matrix: {wide.shape[0]} models × {wide.shape[1]} genes  [{time.time()-t0:.0f}s]")

pairs = [(gene_list[i], gene_list[i + N_PAIRS]) for i in range(N_PAIRS)]

pearson_r  = []
spearman_r = []
n_shared   = []
corr_L2_only = []
corr_L1_only = []

for gene_a, gene_b in pairs:
    if gene_a not in wide.columns or gene_b not in wide.columns:
        continue
    col_a = wide[gene_a]
    col_b = wide[gene_b]
    mask = col_a.notna() & col_b.notna()
    n = mask.sum()
    if n < 10:
        continue
    a = col_a[mask].values
    b = col_b[mask].values
    n_shared.append(n)
    pr = np.corrcoef(a, b)[0, 1]
    sr, _ = stats.spearmanr(a, b)
    pearson_r.append(pr)
    spearman_r.append(sr)

    # Stratify by n_layers
    nl_a = nl_wide[gene_a][mask]
    nl_b = nl_wide[gene_b][mask]
    l2_mask = (nl_a == 2) & (nl_b == 2)
    l1_mask = (nl_a == 1) & (nl_b == 1)
    if l2_mask.sum() >= 10:
        corr_L2_only.append(np.corrcoef(a[l2_mask.values], b[l2_mask.values])[0, 1])
    if l1_mask.sum() >= 10:
        corr_L1_only.append(np.corrcoef(a[l1_mask.values], b[l1_mask.values])[0, 1])

pearson_r  = np.array(pearson_r)
spearman_r = np.array(spearman_r)
n_shared   = np.array(n_shared)

print(f"\n[COUNT] Eligible pairs (>=10 shared lines): {len(pearson_r):,}  [{time.time()-t0:.0f}s]")

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

# Rank displacement in top-20 selectivity (subsample: 100 pairs with >=50 shared)
print(f"\n--- RANK DISPLACEMENT: selectivity score vs correlation-adjusted ---")
eligible_for_disp = [(gene_list[i], gene_list[i+N_PAIRS]) for i, n in enumerate(n_shared) if n >= 50][:100]
displacements = []
for gene_a, gene_b in eligible_for_disp:
    if gene_a not in wide.columns or gene_b not in wide.columns:
        continue
    col_a = wide[gene_a]
    col_b = wide[gene_b]
    mask = col_a.notna() & col_b.notna()
    if mask.sum() < 20:
        continue
    df2 = pd.DataFrame({"sa": col_a[mask].values, "sb": col_b[mask].values})
    df2["sel_raw"] = df2["sa"] * (1 - df2["sb"])
    df2["a_c"] = df2["sa"] - df2["sa"].mean()
    df2["sel_adj"] = df2["a_c"] * (1 - df2["sb"])
    top_raw = set(df2.nlargest(20, "sel_raw").index)
    top_adj = set(df2.nlargest(20, "sel_adj").index)
    displacements.append(20 - len(top_raw & top_adj))

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

print(f"\n[ELAPSED] {time.time()-t0:.0f}s")
print("\n" + "=" * 72)
