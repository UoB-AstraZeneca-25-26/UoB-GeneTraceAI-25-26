"""
T4 — D15: Does protein coverage determine top-20% membership?
Run from repo root: python diagnostics/T4_D15_nlayers_top20.py
"""
import hashlib, sys
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import stats

REPO = Path(__file__).resolve().parent.parent.parent
OUT  = REPO / "Testing and validation" / "diagnostics" / "out"
sys.path.insert(0, str(REPO / "architecture"))
from config import PREDICTIONS, DB
import duckdb

def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()

print("=" * 72)
print("T4 — D15: TOP-20% BY N_LAYERS BRANCH")
print("=" * 72)
print(f"[PROVENANCE] predictions sha256 = {sha256(PREDICTIONS)}")

pred = pd.read_parquet(PREDICTIONS)
print(f"[COUNT] predictions total rows: {len(pred):,}")
print(f"[COUNT] n_layers=1: {(pred.n_layers==1).sum():,}")
print(f"[COUNT] n_layers=2: {(pred.n_layers==2).sum():,}")
print()

# Join lineage
con = duckdb.connect(str(DB), read_only=True)
lin = con.execute(
    "SELECT model_id, lineage FROM main.sample_info WHERE lineage IS NOT NULL"
).df().drop_duplicates("model_id")
con.close()
pred = pred.merge(lin, on="model_id", how="left")
n_missing_lin = pred["lineage"].isna().sum()
print(f"[COUNT] Rows missing lineage after join: {n_missing_lin:,}")
print(f"[COUNT] Distinct lineages: {pred['lineage'].nunique():,}")
assert pred["lineage"].nunique() >= 20, \
    f"INSUFFICIENT CLUSTERS: only {pred['lineage'].nunique()} lineages"

pred["top20"] = pred["score_rank_pct"] >= 0.80
pred["core_z"] = stats.norm.ppf(np.clip(pred["core_score"], 1e-6, 1-1e-6))

# Basic counts
p1_rows = pred[pred.n_layers == 1]
p2_rows = pred[pred.n_layers == 2]
p1_top = p1_rows["top20"].mean()
p2_top = p2_rows["top20"].mean()
sd1 = p1_rows["core_z"].std()
sd2 = p2_rows["core_z"].std()

print(f"[MEASURED] P(top20% | n_layers=1) = {p1_top:.6f}")
print(f"[MEASURED] P(top20% | n_layers=2) = {p2_top:.6f}")
print(f"[MEASURED] Ratio (L1/L2):           {p1_top/p2_top:.4f}")
print(f"[MEASURED] Difference (L1-L2):      {(p1_top-p2_top)*100:.4f} pp")
print(f"[MEASURED] SD(core_z | n_layers=1) = {sd1:.6f}")
print(f"[MEASURED] SD(core_z | n_layers=2) = {sd2:.6f}")
print(f"[MEASURED] SD ratio (L2/L1):        {sd2/sd1:.4f}  (predicted ~0.95)")
print()

# Cluster bootstrap on lineage (1000 resamples)
print("[MEASURED] Cluster bootstrap CI on P(top20%) (1000 resamples, cluster=lineage) ...")
rng = np.random.default_rng(42)
lineages = pred["lineage"].dropna().unique()
n_lin    = len(lineages)
boot_p1, boot_p2 = [], []

for _ in range(1000):
    sampled_lin = rng.choice(lineages, size=n_lin, replace=True)
    lin_counts  = pd.Series(sampled_lin).value_counts()
    frames = []
    for lin_name, cnt in lin_counts.items():
        rows = pred[pred.lineage == lin_name]
        frames.extend([rows] * cnt)
    boot_df = pd.concat(frames, ignore_index=True)
    l1 = boot_df[boot_df.n_layers == 1]["top20"].mean()
    l2 = boot_df[boot_df.n_layers == 2]["top20"].mean()
    boot_p1.append(l1)
    boot_p2.append(l2)

boot_p1 = np.array(boot_p1)
boot_p2 = np.array(boot_p2)
boot_ratio = boot_p1 / np.where(boot_p2 > 0, boot_p2, np.nan)

print(f"[MEASURED] P(top20%|L1) 95% CI: [{np.percentile(boot_p1,2.5):.4f}, {np.percentile(boot_p1,97.5):.4f}]")
print(f"[MEASURED] P(top20%|L2) 95% CI: [{np.percentile(boot_p2,2.5):.4f}, {np.percentile(boot_p2,97.5):.4f}]")
print(f"[MEASURED] Ratio 95% CI:         [{np.nanpercentile(boot_ratio,2.5):.4f}, {np.nanpercentile(boot_ratio,97.5):.4f}]")
print()

# Naive binomial CIs for comparison
def binom_ci(p, n, z=1.96):
    se = np.sqrt(p*(1-p)/n)
    return p - z*se, p + z*se

n1 = len(p1_rows); n2 = len(p2_rows)
ci1_naive = binom_ci(p1_top, n1); ci2_naive = binom_ci(p2_top, n2)
print(f"[MEASURED] Naive CI P(top20%|L1): [{ci1_naive[0]:.4f}, {ci1_naive[1]:.4f}]")
print(f"[MEASURED] Naive CI P(top20%|L2): [{ci2_naive[0]:.4f}, {ci2_naive[1]:.4f}]")
cluster_se_1 = np.std(boot_p1)
naive_se_1   = np.sqrt(p1_top*(1-p1_top)/n1)
print(f"[MEASURED] Cluster SE / Naive SE (L1): {cluster_se_1/naive_se_1:.2f}  (design effect)")
print()

# Stratified by gene: median over genes of [P(top20%|L1) - P(top20%|L2)]
print("[MEASURED] Within-gene stratified difference ...")
per_gene = pred.groupby("ensg_id").apply(lambda g: pd.Series({
    "p1": g.loc[g.n_layers==1,"top20"].mean() if (g.n_layers==1).any() else np.nan,
    "p2": g.loc[g.n_layers==2,"top20"].mean() if (g.n_layers==2).any() else np.nan,
}), include_groups=False).dropna()
per_gene["diff"] = per_gene["p1"] - per_gene["p2"]
per_gene_valid = per_gene.dropna(subset=["diff"])
print(f"[MEASURED] Genes with both branches: {len(per_gene_valid):,}")
print(f"[MEASURED] Median within-gene [P(L1)-P(L2)]: {per_gene_valid['diff'].median()*100:.4f} pp")
print(f"[MEASURED] Mean within-gene [P(L1)-P(L2)]:   {per_gene_valid['diff'].mean()*100:.4f} pp")
print()

# Verdict
diff_pp = (p1_top - p2_top) * 100
print("=" * 72)
print("VERDICT AGAINST PRE-DECLARATION:")
print(f"  Pre-declared: P(top20%|L1) > P(top20%|L2), diff >= 1pp, SD(L2)/SD(L1) ~ 0.95")
if p1_top > p2_top and diff_pp >= 1.0:
    print(f"  CONFIRMED — diff = {diff_pp:.2f} pp (>= 1pp), SD ratio = {sd2/sd1:.4f}")
    print(f"  This is a fame-bias channel: protein-free lines over-represented in top-20%")
elif p1_top > p2_top and diff_pp < 1.0:
    print(f"  PARTIAL — direction confirmed but diff = {diff_pp:.4f} pp < 1pp threshold")
else:
    print(f"  NOT CONFIRMED — direction reversed: L1={p1_top:.4f}, L2={p2_top:.4f}")
print("=" * 72)
