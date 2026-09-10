"""
W5_both_estimands.py
Compute BOTH estimands and report both — DO NOT CHOOSE.

ESTIMAND 1 — per-gene mean
    hr_flat   = mean_g [hits_g / max(sensitive_g, 1)]
    hr_driver = same
    delta = hr_driver - hr_flat

ESTIMAND 2 — pooled ratio  (what V1 computed)
    hr_flat   = sum(hits) / sum(sensitive)
    hr_driver = same
    delta = hr_driver - hr_flat

Both use MAD-z = -0.5, seed=42, n_boot=2000 (cluster bootstrap on gene).

Also reports:
    - Which estimand V1 actually computed (with code quote)
    - n_sensitive per gene: min, median, max, Gini coefficient
    - Concentration: names and shares for top-3 genes

ESCALATION REQUIRED: do not choose between estimands.

Run: python diagnostics/W5_both_estimands.py
"""
import sys
from pathlib import Path
import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "final_pipeline"))
from config import FLAGS_DRIVER, VALIDATION, GENE_LKP

HOLD_OUT_SEED     = 42
HOLD_OUT_FRAC     = 0.20
HITS_AT           = 20
SENSITIVITY_MAD_Z = -0.5
N_BOOT            = 2000

print("=" * 70)
print("W5 — Both estimands (ESCALATION — do not choose)")
print("=" * 70)

flags    = pd.read_parquet(FLAGS_DRIVER)
gdsc_raw = pd.read_parquet(VALIDATION / "gdsc_scored_ready.parquet")
genes    = pd.read_parquet(GENE_LKP, columns=["ensg_id", "gene_role"])

ic50_col = "LN_IC50" if "LN_IC50" in gdsc_raw.columns else "log_ic50"
drug_col = "drug_name" if "drug_name" in gdsc_raw.columns else gdsc_raw.columns[0]

curated   = sorted(set(flags.ensg_id) & set(gdsc_raw.target_ensg))
rng_split = np.random.default_rng(HOLD_OUT_SEED)
test_genes = set(rng_split.choice(curated, size=int(len(curated) * HOLD_OUT_FRAC), replace=False))

def label_per_drug(gdsc, ic50_col, drug_col):
    def _z(grp):
        med = grp[ic50_col].median()
        mad = (grp[ic50_col] - med).abs().median()
        mad_s = mad * 1.4826 if mad > 0 else grp[ic50_col].std(ddof=1)
        mad_s = max(mad_s, 1e-6)
        g2 = grp.copy()
        g2["drug_z"]      = (grp[ic50_col] - med) / mad_s
        g2["is_sensitive"] = g2["drug_z"] <= SENSITIVITY_MAD_Z
        return g2
    return gdsc.groupby(drug_col, group_keys=False).apply(_z, include_groups=False)

gdsc_t = gdsc_raw[gdsc_raw.target_ensg.isin(test_genes)].copy()
gdsc_t = gdsc_t.rename(columns={"target_ensg": "ensg_id"})
gdsc_t["model_id"] = gdsc_t["model_id"].str.lower()
gdsc_t = label_per_drug(gdsc_t, ic50_col, drug_col)

core_t = flags[flags.ensg_id.isin(test_genes)].copy()
core_t = core_t.merge(
    gdsc_t[["ensg_id","model_id","is_sensitive"]].drop_duplicates(["ensg_id","model_id"]),
    on=["ensg_id","model_id"], how="left"
)
core_t["is_sensitive"] = core_t["is_sensitive"].fillna(False)
core_t = core_t.merge(genes, on="ensg_id", how="left")

core_t["rank_flat"] = core_t.groupby("ensg_id")["core_score"].rank(ascending=False, method="first")
core_t["rank_driver"] = core_t.groupby("ensg_id", group_keys=False).apply(
    lambda g: g["core_score"].where(g["has_driver_alteration"], other=-1)
                             .rank(ascending=False, method="first"),
    include_groups=False
)

gene_list = sorted(core_t.ensg_id.unique())
n_genes   = len(gene_list)

# Pre-compute per-gene contributions
flat_c, driver_c = {}, {}
for g, grp in core_t.groupby("ensg_id"):
    flat_c[g]   = (int(((grp["rank_flat"]   <= HITS_AT) & grp["is_sensitive"]).sum()),
                   int(grp["is_sensitive"].sum()))
    driver_c[g] = (int(((grp["rank_driver"] <= HITS_AT) & grp["is_sensitive"]).sum()),
                   int(grp["is_sensitive"].sum()))

n_sens_per_gene = np.array([flat_c[g][1] for g in gene_list])
n_sens_total    = n_sens_per_gene.sum()

print(f"\nTest genes: {n_genes}  |  Total sensitive pairs: {n_sens_total}")

# n_sensitive distribution
n_sens_sorted = np.sort(n_sens_per_gene)[::-1]
gini = (np.sum(np.abs(n_sens_per_gene[:,None] - n_sens_per_gene[None,:])) /
        (2 * n_genes * n_sens_total)) if n_sens_total > 0 else float("nan")
print(f"\nn_sensitive per gene: min={n_sens_per_gene.min()}  "
      f"median={np.median(n_sens_per_gene):.0f}  "
      f"max={n_sens_per_gene.max()}  "
      f"Gini={gini:.3f}")
if gini > 0.5:
    print("  Gini > 0.5 — pooled ratio is concentration-dominated")
    # Top contributors
    top_idx = np.argsort(n_sens_per_gene)[::-1][:3]
    top_genes = [gene_list[i] for i in top_idx]
    for g in top_genes:
        share = flat_c[g][1] / max(n_sens_total, 1)
        print(f"  {g}: n_sensitive={flat_c[g][1]}  share={share*100:.1f}%")

# ── ESTIMAND 1 — per-gene mean ────────────────────────────────────────────
per_gene_hr_f = np.array([flat_c[g][0]/max(flat_c[g][1],1)   for g in gene_list])
per_gene_hr_d = np.array([driver_c[g][0]/max(driver_c[g][1],1) for g in gene_list])
per_gene_delta = per_gene_hr_d - per_gene_hr_f

hr_f_pgm  = per_gene_hr_f.mean()
hr_d_pgm  = per_gene_hr_d.mean()
delta_pgm = per_gene_delta.mean()
lift_pgm  = hr_d_pgm / hr_f_pgm if hr_f_pgm > 0 else float("nan")

# Bootstrap for per-gene mean estimand
genes_arr = np.array(gene_list)
rng_boot  = np.random.default_rng(HOLD_OUT_SEED)
deltas_pgm = np.zeros(N_BOOT)
for i in range(N_BOOT):
    samp = rng_boot.choice(genes_arr, size=n_genes, replace=True)
    f = np.array([flat_c[g][0]/max(flat_c[g][1],1)   for g in samp])
    d = np.array([driver_c[g][0]/max(driver_c[g][1],1) for g in samp])
    deltas_pgm[i] = (d - f).mean()

ci_pgm_lo = float(np.percentile(deltas_pgm, 2.5))
ci_pgm_hi = float(np.percentile(deltas_pgm, 97.5))

print(f"\n{'='*60}")
print(f"ESTIMAND 1 — per-gene mean (each gene counts once)")
print(f"{'='*60}")
print(f"  hr_flat   = {hr_f_pgm:.6f}  (mean across {n_genes} genes)")
print(f"  hr_driver = {hr_d_pgm:.6f}")
print(f"  delta     = {delta_pgm:+.6f}")
print(f"  lift      = {'undef' if np.isnan(lift_pgm) else f'{lift_pgm:.4f}x'}")
print(f"  95% CI (cluster bootstrap, gene clusters): [{ci_pgm_lo:+.6f}, {ci_pgm_hi:+.6f}]")
print(f"  CI excludes zero: {ci_pgm_lo > 0}")
print(f"  per-gene delta: median={np.median(per_gene_delta):+.4f}  "
      f"SD={per_gene_delta.std(ddof=1):.4f}  "
      f"pos={(per_gene_delta>0).sum()}  zero={(per_gene_delta==0).sum()}  neg={(per_gene_delta<0).sum()}")

# ── ESTIMAND 2 — pooled ratio ─────────────────────────────────────────────
num_f = sum(flat_c[g][0] for g in gene_list);   den_f = sum(flat_c[g][1] for g in gene_list)
num_d = sum(driver_c[g][0] for g in gene_list); den_d = sum(driver_c[g][1] for g in gene_list)
hr_f_pool  = num_f / max(den_f, 1)
hr_d_pool  = num_d / max(den_d, 1)
delta_pool = hr_d_pool - hr_f_pool
lift_pool  = hr_d_pool / hr_f_pool if hr_f_pool > 0 else float("nan")

rng_boot2 = np.random.default_rng(HOLD_OUT_SEED)
deltas_pool = np.zeros(N_BOOT)
for i in range(N_BOOT):
    samp = rng_boot2.choice(genes_arr, size=n_genes, replace=True)
    fn = sum(flat_c[g][0] for g in samp);   fd = sum(flat_c[g][1] for g in samp)
    dn = sum(driver_c[g][0] for g in samp); dd = sum(driver_c[g][1] for g in samp)
    deltas_pool[i] = dn/max(dd,1) - fn/max(fd,1)

ci_pool_lo = float(np.percentile(deltas_pool, 2.5))
ci_pool_hi = float(np.percentile(deltas_pool, 97.5))

print(f"\n{'='*60}")
print(f"ESTIMAND 2 — pooled ratio  (V1 computed this)")
print(f"{'='*60}")
print(f"  hr_flat   = {hr_f_pool:.6f}  (num={num_f}, den={den_f})")
print(f"  hr_driver = {hr_d_pool:.6f}  (num={num_d}, den={den_d})")
print(f"  delta     = {delta_pool:+.6f}")
print(f"  lift      = {'undef' if np.isnan(lift_pool) else f'{lift_pool:.4f}x'}")
print(f"  95% CI (cluster bootstrap, gene clusters): [{ci_pool_lo:+.6f}, {ci_pool_hi:+.6f}]")
print(f"  CI excludes zero: {ci_pool_lo > 0}")
print(f"  (V1 measured: delta=+0.008880, CI=[+0.000878, +0.017581])")

# ── Which estimand V1 computed ────────────────────────────────────────────
print(f"\n{'='*60}")
print(f"Which estimand did V1 compute?")
print(f"{'='*60}")
print(f"""  From V1_bootstrap_ci.py lines 99-108:
    num_f = sum(contribs[g][0] for g in genes)   # sum of hits across genes
    den_f = sum(contribs[g][1] for g in genes)   # sum of sensitives across genes
    hr_flat_pt, num_f, den_f = pool_hr(flat_c, gene_list)
    → ESTIMAND 2 (pooled ratio).

  Bootstrap (V1 lines 127-137) resampled genes and recomputed the ratio:
    hr_f = f_n / max(f_d, 1)                      # pooled ratio on resample
    deltas_boot[i] = hr_v - hr_f
  → Bootstrap is for Estimand 2 (pooled ratio), clustered on gene.
""")

# ── Comparison ────────────────────────────────────────────────────────────
print(f"{'='*60}")
print(f"Comparison")
print(f"{'='*60}")
diff = abs(delta_pgm - delta_pool)
pct_of_smaller = diff / min(abs(delta_pgm), abs(delta_pool)) if min(abs(delta_pgm), abs(delta_pool)) > 0 else float("nan")
print(f"  delta_PGM  = {delta_pgm:+.6f}  CI=[{ci_pgm_lo:+.6f}, {ci_pgm_hi:+.6f}]")
print(f"  delta_POOL = {delta_pool:+.6f}  CI=[{ci_pool_lo:+.6f}, {ci_pool_hi:+.6f}]")
print(f"  |difference| = {diff:.6f}  ({pct_of_smaller*100:.1f}% of smaller)")
print(f"  PRE-DECLARED >20% difference: {'PASS' if pct_of_smaller > 0.20 else 'FAIL'}")

print(f"""
ESCALATION NOTE (do not resolve autonomously):
  Both estimands are reported above. A recommendation follows as requested:

  Recommendation: per-gene mean (Estimand 1)
  Reasoning: a researcher encounters this tool one gene at a time.
  The per-gene mean measures "does the driver gate improve hit-rate
  for a typical gene?" — which is the claim. The pooled ratio measures
  "does the gate improve hits per sensitive pair?" and weights
  high-n_sensitive genes (fame-correlated) more heavily.
  Gini = {gini:.3f} — {'concentrated, weighting matters' if gini > 0.5 else 'not concentrated, weighting similar'}.

  FIONA'S DECISION REQUIRED before W7 runs.
""")
