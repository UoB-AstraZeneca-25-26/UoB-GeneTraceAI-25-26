"""
T4 — D15 (fast): cluster bootstrap on lineage-level rates (not raw rows).
Run from repo root: python diagnostics/T4_D15_fast.py
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
print("T4 — D15: TOP-20% BY N_LAYERS BRANCH (cluster bootstrap, fast)")
print("=" * 72)
print(f"[PROVENANCE] predictions sha256 = {sha256(PREDICTIONS)}")

pred = pd.read_parquet(PREDICTIONS, columns=["model_id","ensg_id","core_score","n_layers","score_rank_pct"])
print(f"[COUNT] predictions: {len(pred):,} rows")

# Join lineage
con = duckdb.connect(str(DB), read_only=True)
lin = con.execute("SELECT model_id, lineage FROM main.sample_info WHERE lineage IS NOT NULL").df()
lin = lin.drop_duplicates("model_id")
con.close()

pred = pred.merge(lin, on="model_id", how="left")
n_miss = pred["lineage"].isna().sum()
print(f"[COUNT] Rows missing lineage: {n_miss:,}")
print(f"[COUNT] Distinct lineages: {pred['lineage'].nunique():,}")
assert pred["lineage"].nunique() >= 20, f"INSUFFICIENT CLUSTERS: {pred['lineage'].nunique()}"

pred["top20"] = pred["score_rank_pct"] >= 0.80
pred["core_z"] = stats.norm.ppf(np.clip(pred["core_score"], 1e-9, 1-1e-9))

# ── Point estimates ─────────────────────────────────────────────────────────
p1 = pred[pred.n_layers == 1]
p2 = pred[pred.n_layers == 2]
P1_obs = p1["top20"].mean()
P2_obs = p2["top20"].mean()
sd1    = p1["core_z"].std()
sd2    = p2["core_z"].std()
print(f"\n[MEASURED] P(top20% | n_layers=1) = {P1_obs:.6f}  (n={len(p1):,})")
print(f"[MEASURED] P(top20% | n_layers=2) = {P2_obs:.6f}  (n={len(p2):,})")
print(f"[MEASURED] Ratio  (L1/L2):          {P1_obs/P2_obs:.4f}")
print(f"[MEASURED] Difference (L1-L2):      {(P1_obs-P2_obs)*100:.4f} pp")
print(f"[MEASURED] SD(core_z | n_layers=1) = {sd1:.6f}")
print(f"[MEASURED] SD(core_z | n_layers=2) = {sd2:.6f}")
print(f"[MEASURED] SD ratio SD(L2)/SD(L1):  {sd2/sd1:.4f}  (predicted ~0.95)")

# ── Fast cluster bootstrap: resample lineages, use pre-agg rates ────────────
print("\n[MEASURED] Cluster bootstrap (1000 resamples, cluster=lineage) ...")
agg = pred.groupby(["lineage","n_layers"]).agg(
    n_top=("top20","sum"), n_total=("top20","count")
).reset_index()
lineages = agg["lineage"].unique()
n_lin    = len(lineages)
print(f"  Lineages used: {n_lin}")

rng = np.random.default_rng(42)
boot_p1, boot_p2 = [], []
for _ in range(1000):
    sampled = rng.choice(lineages, size=n_lin, replace=True)
    counts  = pd.Series(sampled).value_counts().reset_index()
    counts.columns = ["lineage","multiplicity"]
    tmp = agg.merge(counts, on="lineage")
    tmp["n_top_w"]   = tmp["n_top"]   * tmp["multiplicity"]
    tmp["n_total_w"] = tmp["n_total"] * tmp["multiplicity"]
    by_layer = tmp.groupby("n_layers")[["n_top_w","n_total_w"]].sum()
    boot_p1.append(by_layer.loc[1,"n_top_w"] / by_layer.loc[1,"n_total_w"] if 1 in by_layer.index else np.nan)
    boot_p2.append(by_layer.loc[2,"n_top_w"] / by_layer.loc[2,"n_total_w"] if 2 in by_layer.index else np.nan)

boot_p1 = np.array(boot_p1); boot_p2 = np.array(boot_p2)
boot_ratio = np.where(boot_p2 > 0, boot_p1/boot_p2, np.nan)

print(f"\n[MEASURED] Cluster-bootstrap 95% CI P(top20%|L1): [{np.percentile(boot_p1,2.5):.6f}, {np.percentile(boot_p1,97.5):.6f}]")
print(f"[MEASURED] Cluster-bootstrap 95% CI P(top20%|L2): [{np.percentile(boot_p2,2.5):.6f}, {np.percentile(boot_p2,97.5):.6f}]")
print(f"[MEASURED] Cluster-bootstrap 95% CI Ratio:         [{np.nanpercentile(boot_ratio,2.5):.4f}, {np.nanpercentile(boot_ratio,97.5):.4f}]")

# Naive CIs for comparison
se1_naive = np.sqrt(P1_obs*(1-P1_obs)/len(p1))
se2_naive = np.sqrt(P2_obs*(1-P2_obs)/len(p2))
se1_clust = np.std(boot_p1)
se2_clust = np.std(boot_p2)
print(f"\n[MEASURED] Naive 95% CI P(L1): [{P1_obs-1.96*se1_naive:.6f}, {P1_obs+1.96*se1_naive:.6f}]")
print(f"[MEASURED] Naive 95% CI P(L2): [{P2_obs-1.96*se2_naive:.6f}, {P2_obs+1.96*se2_naive:.6f}]")
print(f"[MEASURED] Cluster SE / Naive SE (L1): {se1_clust/se1_naive:.2f}  (design effect)")
print(f"[MEASURED] Cluster SE / Naive SE (L2): {se2_clust/se2_naive:.2f}")

# ── Stratified within-gene ──────────────────────────────────────────────────
print("\n[MEASURED] Within-gene stratification ...")
per_gene = pred.groupby(["ensg_id","n_layers"])["top20"].mean().unstack("n_layers")
per_gene.columns = [f"p_L{int(c)}" for c in per_gene.columns]
per_gene_both = per_gene.dropna()
diff_col = per_gene_both.get("p_L1", per_gene_both.iloc[:,0]) - per_gene_both.get("p_L2", per_gene_both.iloc[:,1])
# re-compute properly:
if "p_L1" in per_gene_both.columns and "p_L2" in per_gene_both.columns:
    diff_col = per_gene_both["p_L1"] - per_gene_both["p_L2"]
    print(f"[MEASURED] Genes with both L1 and L2 rows: {len(per_gene_both):,}")
    print(f"[MEASURED] Median within-gene [P(L1)-P(L2)]: {diff_col.median()*100:.4f} pp")
    print(f"[MEASURED] Mean within-gene [P(L1)-P(L2)]:   {diff_col.mean()*100:.4f} pp")
else:
    print("[MEASURED] Columns available:", list(per_gene_both.columns), "— cannot compute diff")

# ── Verdict ──────────────────────────────────────────────────────────────────
diff_pp = (P1_obs - P2_obs) * 100
ci_diff_lo = np.percentile(boot_p1-boot_p2, 2.5)*100
ci_diff_hi = np.percentile(boot_p1-boot_p2, 97.5)*100
print("\n" + "=" * 72)
print("VERDICT AGAINST PRE-DECLARATION:")
print(f"  Pre-declared: P(L1) > P(L2), diff >= 1pp, SD(L2)/SD(L1) ~ 0.95")
print(f"  Measured: diff = {diff_pp:.2f} pp  CI=[{ci_diff_lo:.2f}, {ci_diff_hi:.2f}] pp")
print(f"  SD ratio = {sd2/sd1:.4f}")
if P1_obs > P2_obs and diff_pp >= 1.0 and ci_diff_lo > 0:
    print("  CONFIRMED — direction predicted, diff >= 1pp, CI excludes 0")
    print("  Fame-bias channel: protein-free lines over-represented in top-20%")
elif P1_obs > P2_obs and diff_pp < 1.0:
    print(f"  PARTIAL — correct direction but diff {diff_pp:.2f} pp < 1pp threshold")
elif P1_obs > P2_obs and ci_diff_lo <= 0:
    print(f"  INCONCLUSIVE — direction correct, diff {diff_pp:.2f} pp, but CI includes 0")
else:
    print(f"  NOT CONFIRMED — L1={P1_obs:.4f}, L2={P2_obs:.4f}")
print("=" * 72)
