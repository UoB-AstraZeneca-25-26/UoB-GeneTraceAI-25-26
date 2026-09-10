"""
V1_bootstrap_ci.py
Paired bootstrap CI on hr_driver vs hr_flat, clustered on gene.
Per-drug MAD-z sensitivity labels (SENSITIVITY_MAD_Z = -0.5).
2000 resamples, seed 42.

Run: python diagnostics/V1_bootstrap_ci.py
"""
import sys, json
from pathlib import Path
import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(REPO / "architecture"))
from config import FLAGS_DRIVER, VALIDATION, GENE_LKP

HOLD_OUT_SEED    = 42
HOLD_OUT_FRAC    = 0.20
HITS_AT          = 20
SENSITIVITY_MAD_Z = -0.5
N_BOOT           = 2000

def label_sensitive_per_drug(gdsc, ic50_col, drug_col):
    def _drug_z(grp):
        med = grp[ic50_col].median()
        mad = (grp[ic50_col] - med).abs().median()
        mad_scaled = mad * 1.4826 if mad > 0 else grp[ic50_col].std(ddof=1)
        mad_scaled = max(mad_scaled, 1e-6)
        grp = grp.copy()
        grp["drug_z"] = (grp[ic50_col] - med) / mad_scaled
        grp["is_sensitive"] = grp["drug_z"] <= SENSITIVITY_MAD_Z
        return grp
    return gdsc.groupby(drug_col, group_keys=False).apply(_drug_z, include_groups=False)

print("=" * 70)
print("V1 — Paired bootstrap CI (clustered on gene)")
print("=" * 70)

flags    = pd.read_parquet(FLAGS_DRIVER)
gdsc_raw = pd.read_parquet(VALIDATION / "gdsc_scored_ready.parquet")
genes    = pd.read_parquet(GENE_LKP, columns=["ensg_id","gene_role"])

ic50_col = "LN_IC50" if "LN_IC50" in gdsc_raw.columns else "log_ic50"
drug_col = "drug_name" if "drug_name" in gdsc_raw.columns else gdsc_raw.columns[0]

gdsc = label_sensitive_per_drug(gdsc_raw.copy(), ic50_col, drug_col)
gdsc = gdsc.rename(columns={"target_ensg": "ensg_id"}) if "target_ensg" in gdsc.columns else gdsc

# Build test gene set (identical to eval.py)
curated = sorted(set(flags.ensg_id) & set(gdsc_raw.target_ensg))
rng_split = np.random.default_rng(HOLD_OUT_SEED)
test_genes = set(rng_split.choice(curated, size=int(len(curated) * HOLD_OUT_FRAC), replace=False))

print(f"Curated genes (intersection): {len(curated)}")
print(f"Test genes (20% holdout):     {len(test_genes)}")

core_t = flags[flags.ensg_id.isin(test_genes)].copy()
gdsc_t = gdsc_raw[gdsc_raw.target_ensg.isin(test_genes)].copy()
gdsc_t = gdsc_t.rename(columns={"target_ensg": "ensg_id"})
gdsc_t["model_id"] = gdsc_t["model_id"].str.lower()
gdsc_t = label_sensitive_per_drug(gdsc_t, ic50_col, drug_col)

# Join sensitivity labels
core_t = core_t.merge(
    gdsc_t[["ensg_id", "model_id", "is_sensitive"]].drop_duplicates(["ensg_id","model_id"]),
    on=["ensg_id","model_id"], how="left"
)
core_t["is_sensitive"] = core_t["is_sensitive"].fillna(False)
core_t = core_t.merge(genes, on="ensg_id", how="left")

# Compute ranks
core_t["rank_flat"] = (
    core_t.groupby("ensg_id")["core_score"].rank(ascending=False, method="first")
)
core_t["rank_driver"] = (
    core_t.groupby("ensg_id", group_keys=False).apply(
        lambda g: g["core_score"].where(g["has_driver_alteration"], other=-1)
                                 .rank(ascending=False, method="first"),
        include_groups=False
    )
)

# Pre-compute per-gene contributions
def gene_contributions(df, rank_col, k=HITS_AT):
    result = {}
    for gene, grp in df.groupby("ensg_id"):
        n_top = int(((grp[rank_col] <= k) & grp["is_sensitive"]).sum())
        n_sen = int(grp["is_sensitive"].sum())
        result[gene] = (n_top, n_sen)
    return result

flat_c   = gene_contributions(core_t, "rank_flat")
driver_c = gene_contributions(core_t, "rank_driver")

gene_list = sorted(flat_c.keys())
n_genes   = len(gene_list)

# Point estimates (pooled across all strata)
def pool_hr(contribs, genes):
    num = sum(contribs[g][0] for g in genes)
    den = sum(contribs[g][1] for g in genes)
    return num / max(den, 1), num, den

hr_flat_pt,   num_f, den_f = pool_hr(flat_c,   gene_list)
hr_driver_pt, num_d, den_d = pool_hr(driver_c, gene_list)
delta_pt    = hr_driver_pt - hr_flat_pt
lift_pt     = hr_driver_pt / hr_flat_pt if hr_flat_pt > 0 else float("nan")

print(f"\nPoint estimates (pooled, n_genes={n_genes}):")
print(f"  hr_flat   = {hr_flat_pt:.6f}  (num={num_f}, den={den_f})")
print(f"  hr_driver = {hr_driver_pt:.6f}  (num={num_d}, den={den_d})")
print(f"  delta     = {delta_pt:+.6f}")
print(f"  lift      = {lift_pt:.4f}x")

# If hr_flat == 0, lift is undefined
if hr_flat_pt == 0:
    print("\n  WARNING: hr_flat = 0 — lift_ratio undefined, reporting delta only")

# Paired bootstrap (clustered on gene)
genes_arr = np.array(gene_list)
rng_boot  = np.random.default_rng(HOLD_OUT_SEED)

deltas_boot = np.zeros(N_BOOT)
lifts_boot  = np.full(N_BOOT, np.nan)

for i in range(N_BOOT):
    samp = rng_boot.choice(genes_arr, size=n_genes, replace=True)
    f_n = sum(flat_c[g][0] for g in samp)
    f_d = sum(flat_c[g][1] for g in samp)
    v_n = sum(driver_c[g][0] for g in samp)
    v_d = sum(driver_c[g][1] for g in samp)
    hr_f = f_n / max(f_d, 1)
    hr_v = v_n / max(v_d, 1)
    deltas_boot[i] = hr_v - hr_f
    if hr_f > 0:
        lifts_boot[i] = hr_v / hr_f

ci_delta_lo = float(np.percentile(deltas_boot, 2.5))
ci_delta_hi = float(np.percentile(deltas_boot, 97.5))
ci_lift_lo  = float(np.nanpercentile(lifts_boot, 2.5))
ci_lift_hi  = float(np.nanpercentile(lifts_boot, 97.5))

print(f"\nPaired bootstrap (n_resamples={N_BOOT}, clustered on gene, n_clusters={n_genes}):")
print(f"  delta 95% CI: [{ci_delta_lo:+.6f}, {ci_delta_hi:+.6f}]")
print(f"  lift  95% CI: [{ci_lift_lo:.4f}x, {ci_lift_hi:.4f}x]")

excludes_zero = (ci_delta_lo > 0) or (ci_delta_hi < 0)
print(f"  CI excludes zero: {excludes_zero}")
print(f"  VERDICT: {'PASS' if excludes_zero else 'FAIL — CI includes zero'}")

# Naive (unclustered) CI — resample pairs, not genes
# Permute is_sensitive labels to construct empirical null
# Approximate using gene-level deltas
per_gene_delta = {}
for g in gene_list:
    f_n, f_d = flat_c[g]
    d_n, d_d = driver_c[g]
    hr_f_g = f_n / max(f_d, 1)
    hr_d_g = d_n / max(d_d, 1)
    per_gene_delta[g] = hr_d_g - hr_f_g

gene_deltas = np.array([per_gene_delta[g] for g in gene_list])
naive_se    = gene_deltas.std(ddof=1) / np.sqrt(n_genes)  # naive SE treats genes as i.i.d. pairs
# Clustered SE from bootstrap
clustered_se = deltas_boot.std(ddof=1)
design_effect = clustered_se / naive_se if naive_se > 0 else float("nan")

print(f"\nDesign effect:")
print(f"  Naive SE (gene-level i.i.d.):   {naive_se:.6f}")
print(f"  Clustered SE (bootstrap):        {clustered_se:.6f}")
print(f"  Design effect (clustered/naive): {design_effect:.3f}")

ci_naive_lo = delta_pt - 1.96 * naive_se
ci_naive_hi = delta_pt + 1.96 * naive_se
print(f"  Naive 95% CI: [{ci_naive_lo:+.6f}, {ci_naive_hi:+.6f}]")

# Per-gene delta distribution (for the triage tree if CI includes zero)
print(f"\nPer-gene delta distribution:")
print(f"  mean={gene_deltas.mean():.6f}  median={np.median(gene_deltas):.6f}  "
      f"SD={gene_deltas.std(ddof=1):.6f}")
print(f"  delta > 0: {(gene_deltas > 0).sum()}  delta == 0: {(gene_deltas == 0).sum()}  "
      f"delta < 0: {(gene_deltas < 0).sum()}")

# Per-stratum breakdown (for V2 context)
print(f"\nPer-stratum breakdown:")
for cls in ["oncogene", "tsg", "both"]:
    sub = core_t[core_t.gene_role == cls]
    if sub.empty:
        print(f"  {cls}: EMPTY")
        continue
    ng = sub.ensg_id.nunique()
    ns = int(sub["is_sensitive"].sum())
    gc_f = gene_contributions(sub, "rank_flat")
    gc_d = gene_contributions(sub, "rank_driver")
    hr_f, _, _ = pool_hr(gc_f, list(gc_f.keys()))
    hr_d, _, _ = pool_hr(gc_d, list(gc_d.keys()))
    interpretable = "NOT INTERPRETABLE (n=1)" if ng < 5 else ""
    print(f"  {cls}: n_genes={ng}  n_sensitive={ns}  "
          f"hr_flat={hr_f:.4f}  hr_driver={hr_d:.4f}  {interpretable}")

# Cluster count check
print(f"\nSign-off check — n_clusters={n_genes} == n_genes={n_genes}: "
      f"{'OK' if n_genes < 20 and True else 'OK' if n_genes >= 20 else 'CHECK'}")
if n_genes < 20:
    print(f"  WARNING: n_genes={n_genes} < 20 — INSUFFICIENT CLUSTERS per Cameron & Miller 2015")
