"""
W7_consolidated_bootstrap.py
Consolidated bootstrap matrix: threshold x estimand x method.

Thresholds : -0.25, -0.5, -1.0, -1.5, -2.0
Estimands  : per-gene mean (E1), pooled ratio (E2)
Methods    : A (percentile cluster bootstrap, 2000 resamples)
             B (paired permutation test / Rademacher sign test, 10000 perms)
             C (wild cluster bootstrap, Rademacher weights, 2000 resamples)

References:
  MacKinnon & Webb 2017 (DOI 10.1002/jae.2508) -- percentile bootstrap under-covers at few clusters
  Cameron & Miller 2015 (DOI 10.3368/jhr.50.2.317) -- cluster-robust inference
  Good 2005 (ISBN 0387988645) -- permutation tests as exact small-sample inference
  Iorio et al. 2016 (PMID 27397505) -- GDSC per-drug convention (~17% sensitive)

Run: python diagnostics/W7_consolidated_bootstrap.py
"""
import sys
from pathlib import Path
import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "final_pipeline"))
from config import FLAGS_DRIVER, VALIDATION, GENE_LKP

HOLD_OUT_SEED  = 42
HOLD_OUT_FRAC  = 0.20
HITS_AT        = 20
THRESHOLDS     = [-0.25, -0.5, -1.0, -1.5, -2.0]
N_BOOT         = 2000
N_PERM         = 10000
STABILITY_THRESH = 5    # flag if median n_sensitive/gene falls below this

print("=" * 70)
print("W7 -- Consolidated bootstrap matrix")
print("=" * 70)

flags    = pd.read_parquet(FLAGS_DRIVER)
gdsc_raw = pd.read_parquet(VALIDATION / "gdsc_scored_ready.parquet")
genes    = pd.read_parquet(GENE_LKP, columns=["ensg_id", "gene_role"])

ic50_col = "LN_IC50" if "LN_IC50" in gdsc_raw.columns else "log_ic50"
drug_col = "drug_name" if "drug_name" in gdsc_raw.columns else gdsc_raw.columns[0]

curated   = sorted(set(flags.ensg_id) & set(gdsc_raw.target_ensg))
rng_split = np.random.default_rng(HOLD_OUT_SEED)
test_genes = set(rng_split.choice(curated, size=int(len(curated) * HOLD_OUT_FRAC), replace=False))

# Pre-compute drug stats (for fast per-drug z)
drug_stats = gdsc_raw.groupby(drug_col)[ic50_col].agg(
    med="median",
    mad=lambda x: (x - x.median()).abs().median()
).reset_index()
drug_stats["mad_s"] = (drug_stats["mad"] * 1.4826).clip(lower=1e-6)
gdsc_z = gdsc_raw.merge(drug_stats[[drug_col, "med", "mad_s"]], on=drug_col, how="left")
gdsc_z["drug_z"] = (gdsc_z[ic50_col] - gdsc_z["med"]) / gdsc_z["mad_s"]

# Pre-filter to test genes and normalise model_id
gdsc_t_base = gdsc_z[gdsc_z.target_ensg.isin(test_genes)].copy()
gdsc_t_base = gdsc_t_base.rename(columns={"target_ensg": "ensg_id"})
gdsc_t_base["model_id"] = gdsc_t_base["model_id"].str.lower()

# Compute ranks once (invariant to sensitivity label)
core_base = flags[flags.ensg_id.isin(test_genes)].copy()
core_base = core_base.merge(genes, on="ensg_id", how="left")
core_base["rank_flat"] = core_base.groupby("ensg_id")["core_score"].rank(ascending=False, method="first")
core_base["rank_driver"] = core_base.groupby("ensg_id", group_keys=False).apply(
    lambda g: g["core_score"].where(g["has_driver_alteration"], other=-1)
                             .rank(ascending=False, method="first"),
    include_groups=False
)

gene_list = sorted(core_base.ensg_id.unique())
n_genes   = len(gene_list)

print(f"Test genes: {n_genes}  |  Curated intersection: {len(curated)}")
print(f"Methods: A=percentile cluster bootstrap ({N_BOOT} resamples)  "
      f"B=permutation sign test ({N_PERM} perms)  "
      f"C=wild cluster bootstrap ({N_BOOT} resamples)")

results = []

for thr in THRESHOLDS:
    # Label sensitivity
    gdsc_t = gdsc_t_base[["ensg_id","model_id","drug_z"]].copy()
    gdsc_t["is_sensitive"] = gdsc_t["drug_z"] <= thr

    core_t = core_base.merge(
        gdsc_t[["ensg_id","model_id","is_sensitive"]].drop_duplicates(["ensg_id","model_id"]),
        on=["ensg_id","model_id"], how="left"
    )
    core_t["is_sensitive"] = core_t["is_sensitive"].fillna(False)

    # Guard and stability
    per_drug_rate = gdsc_t_base.assign(is_s=gdsc_t_base["drug_z"] <= thr).groupby(drug_col)["is_s"].mean()
    guard = ((per_drug_rate >= 0.05) & (per_drug_rate <= 0.40)).mean()
    med_rate = per_drug_rate.median()
    n_sensitive_total = int(core_t["is_sensitive"].sum())

    # Per-gene contributions
    flat_n,   flat_d   = {}, {}
    driver_n, driver_d = {}, {}
    pg_delta_flat = {}   # per-gene hit@k by strategy (for E1)
    for g, grp in core_t.groupby("ensg_id"):
        fn = int(((grp["rank_flat"]   <= HITS_AT) & grp["is_sensitive"]).sum())
        dn = int(((grp["rank_driver"] <= HITS_AT) & grp["is_sensitive"]).sum())
        d  = int(grp["is_sensitive"].sum())
        flat_n[g] = fn; flat_d[g] = d
        driver_n[g] = dn; driver_d[g] = d
        pg_delta_flat[g] = dn/max(d,1) - fn/max(d,1)

    n_sens_per_gene = np.array([flat_d[g] for g in gene_list])
    med_n_sens = float(np.median(n_sens_per_gene))
    n_stable = int((n_sens_per_gene >= STABILITY_THRESH).sum())
    stable = med_n_sens >= STABILITY_THRESH

    # ── ESTIMAND 1: per-gene mean ─────────────────────────────────────────
    pg_deltas = np.array([pg_delta_flat[g] for g in gene_list])
    e1_flat   = float(np.mean([flat_n[g]/max(flat_d[g],1) for g in gene_list]))
    e1_driver = float(np.mean([driver_n[g]/max(driver_d[g],1) for g in gene_list]))
    e1_delta  = float(np.mean(pg_deltas))
    e1_lift   = e1_driver / e1_flat if e1_flat > 0 else float("nan")

    # Method A: percentile bootstrap (per-gene mean)
    rng_a = np.random.default_rng(HOLD_OUT_SEED)
    genes_arr = np.array(gene_list)
    bs_a_e1 = np.zeros(N_BOOT)
    for i in range(N_BOOT):
        samp = rng_a.choice(len(gene_list), size=n_genes, replace=True)
        bs_a_e1[i] = pg_deltas[samp].mean()
    ci_a_e1_lo = float(np.percentile(bs_a_e1, 2.5))
    ci_a_e1_hi = float(np.percentile(bs_a_e1, 97.5))

    # Method B: permutation sign test (per-gene mean)
    rng_b = np.random.default_rng(HOLD_OUT_SEED)
    perm_b_e1 = np.zeros(N_PERM)
    for i in range(N_PERM):
        w = rng_b.choice([-1, 1], size=n_genes)
        perm_b_e1[i] = np.mean(w * pg_deltas)
    p_b_e1 = float(np.mean(np.abs(perm_b_e1) >= abs(e1_delta)))

    # Method C: wild cluster bootstrap (per-gene mean)
    rng_c = np.random.default_rng(HOLD_OUT_SEED)
    eps_e1 = pg_deltas - e1_delta   # centered residuals
    bs_c_e1 = np.zeros(N_BOOT)
    for i in range(N_BOOT):
        w = rng_c.choice([-1, 1], size=n_genes)
        bs_c_e1[i] = e1_delta + np.mean(w * eps_e1)
    ci_c_e1_lo = float(np.percentile(bs_c_e1, 2.5))
    ci_c_e1_hi = float(np.percentile(bs_c_e1, 97.5))

    # ── ESTIMAND 2: pooled ratio ──────────────────────────────────────────
    sum_fn = sum(flat_n[g] for g in gene_list)
    sum_fd = sum(flat_d[g] for g in gene_list)
    sum_dn = sum(driver_n[g] for g in gene_list)
    sum_dd = sum(driver_d[g] for g in gene_list)
    e2_flat   = sum_fn / max(sum_fd, 1)
    e2_driver = sum_dn / max(sum_dd, 1)
    e2_delta  = e2_driver - e2_flat
    e2_lift   = e2_driver / e2_flat if e2_flat > 0 else float("nan")

    # Method A: percentile bootstrap (pooled ratio)
    rng_a2 = np.random.default_rng(HOLD_OUT_SEED)
    # Pre-build arrays
    fn_arr = np.array([flat_n[g]   for g in gene_list])
    fd_arr = np.array([flat_d[g]   for g in gene_list])
    dn_arr = np.array([driver_n[g] for g in gene_list])
    dd_arr = np.array([driver_d[g] for g in gene_list])
    bs_a_e2 = np.zeros(N_BOOT)
    for i in range(N_BOOT):
        idx = rng_a2.choice(n_genes, size=n_genes, replace=True)
        hf = fn_arr[idx].sum() / max(fd_arr[idx].sum(), 1)
        hd = dn_arr[idx].sum() / max(dd_arr[idx].sum(), 1)
        bs_a_e2[i] = hd - hf
    ci_a_e2_lo = float(np.percentile(bs_a_e2, 2.5))
    ci_a_e2_hi = float(np.percentile(bs_a_e2, 97.5))

    # Method B: permutation sign test (pooled ratio)
    # Under null: within each gene, swap hits_driver and hits_flat
    rng_b2 = np.random.default_rng(HOLD_OUT_SEED)
    perm_b_e2 = np.zeros(N_PERM)
    num_delta_arr = (dn_arr - fn_arr).astype(float)  # per-gene numerator contribution
    D = sum_fd
    for i in range(N_PERM):
        w = rng_b2.choice([-1, 1], size=n_genes)
        perm_b_e2[i] = (w * num_delta_arr).sum() / max(D, 1)
    p_b_e2 = float(np.mean(np.abs(perm_b_e2) >= abs(e2_delta)))

    # Method C: wild cluster bootstrap (pooled ratio)
    rng_c2 = np.random.default_rng(HOLD_OUT_SEED)
    eps_e2 = num_delta_arr - num_delta_arr.mean()
    bs_c_e2 = np.zeros(N_BOOT)
    for i in range(N_BOOT):
        w = rng_c2.choice([-1, 1], size=n_genes)
        bs_delta = e2_delta + (w * eps_e2).sum() / max(D, 1)
        bs_c_e2[i] = bs_delta
    ci_c_e2_lo = float(np.percentile(bs_c_e2, 2.5))
    ci_c_e2_hi = float(np.percentile(bs_c_e2, 97.5))

    # Width comparison (wild vs percentile, expected wild wider)
    width_a_e1 = ci_a_e1_hi - ci_a_e1_lo
    width_c_e1 = ci_c_e1_hi - ci_c_e1_lo
    width_a_e2 = ci_a_e2_hi - ci_a_e2_lo
    width_c_e2 = ci_c_e2_hi - ci_c_e2_lo

    row = {
        "thr": thr,
        "guard_pct": round(guard*100, 1),
        "med_sens_pct": round(med_rate*100, 1),
        "n_sensitive": n_sensitive_total,
        "med_n_sens_gene": round(med_n_sens, 1),
        "n_stable": n_stable,
        "stable": stable,
        # E1
        "e1_flat": round(e1_flat, 6),
        "e1_driver": round(e1_driver, 6),
        "e1_delta": round(e1_delta, 6),
        "e1_lift": round(e1_lift, 4) if not (isinstance(e1_lift, float) and np.isnan(e1_lift)) else "undef",
        "e1_ci_A_lo": round(ci_a_e1_lo, 6), "e1_ci_A_hi": round(ci_a_e1_hi, 6),
        "e1_A_excl0": ci_a_e1_lo > 0,
        "e1_p_B": round(p_b_e1, 4),
        "e1_B_sig": p_b_e1 < 0.05,
        "e1_ci_C_lo": round(ci_c_e1_lo, 6), "e1_ci_C_hi": round(ci_c_e1_hi, 6),
        "e1_C_excl0": ci_c_e1_lo > 0,
        "e1_wild_wider": width_c_e1 > width_a_e1,
        # E2
        "e2_flat": round(e2_flat, 6),
        "e2_driver": round(e2_driver, 6),
        "e2_delta": round(e2_delta, 6),
        "e2_lift": round(e2_lift, 4) if not (isinstance(e2_lift, float) and np.isnan(e2_lift)) else "undef",
        "e2_ci_A_lo": round(ci_a_e2_lo, 6), "e2_ci_A_hi": round(ci_a_e2_hi, 6),
        "e2_A_excl0": ci_a_e2_lo > 0,
        "e2_p_B": round(p_b_e2, 4),
        "e2_B_sig": p_b_e2 < 0.05,
        "e2_ci_C_lo": round(ci_c_e2_lo, 6), "e2_ci_C_hi": round(ci_c_e2_hi, 6),
        "e2_C_excl0": ci_c_e2_lo > 0,
        "e2_wild_wider": width_c_e2 > width_a_e2,
    }
    results.append(row)

    stab_flag = " [UNSTABLE: med_n_sens < 5]" if not stable else ""
    print(f"\n{'='*60}")
    print(f"Threshold {thr:5.2f}  guard={guard*100:.0f}%  med_sens={med_rate*100:.0f}%  "
          f"n_sensitive={n_sensitive_total}  med_n_sens/gene={med_n_sens:.0f}{stab_flag}")
    print(f"  ESTIMAND 1 (per-gene mean):")
    print(f"    delta={e1_delta:+.6f}  lift={e1_lift:.4f}x  (flat={e1_flat:.6f}  driver={e1_driver:.6f})")
    print(f"    Method A (percentile):  [{ci_a_e1_lo:+.6f}, {ci_a_e1_hi:+.6f}]  excl0={ci_a_e1_lo>0}")
    print(f"    Method B (permutation): p={p_b_e1:.4f}  sig={p_b_e1<0.05}  "
          f"(width_A={width_a_e1:.6f})")
    print(f"    Method C (wild boot):   [{ci_c_e1_lo:+.6f}, {ci_c_e1_hi:+.6f}]  excl0={ci_c_e1_lo>0}  "
          f"wider={width_c_e1>width_a_e1}  (width_C={width_c_e1:.6f})")
    print(f"  ESTIMAND 2 (pooled ratio):")
    print(f"    delta={e2_delta:+.6f}  lift={e2_lift:.4f}x  (flat={e2_flat:.6f}  driver={e2_driver:.6f})")
    print(f"    Method A (percentile):  [{ci_a_e2_lo:+.6f}, {ci_a_e2_hi:+.6f}]  excl0={ci_a_e2_lo>0}")
    print(f"    Method B (permutation): p={p_b_e2:.4f}  sig={p_b_e2<0.05}  "
          f"(width_A={width_a_e2:.6f})")
    print(f"    Method C (wild boot):   [{ci_c_e2_lo:+.6f}, {ci_c_e2_hi:+.6f}]  excl0={ci_c_e2_lo>0}  "
          f"wider={width_c_e2>width_a_e2}  (width_C={width_c_e2:.6f})")

print(f"\n{'='*70}")
print("FULL RESULTS TABLE")
print(f"{'='*70}")
df = pd.DataFrame(results)
print(df.to_string(index=False))

# Pre-declared checks
print(f"\n{'='*70}")
print("PRE-DECLARED CHECKS")
print(f"{'='*70}")

for thr in [-0.5, -1.0]:
    r = next(x for x in results if x["thr"] == thr)
    print(f"\nThreshold {thr}:")
    # E1
    print(f"  E1 Method A (percentile) excl 0: {r['e1_A_excl0']}  (pre-declared: True)")
    print(f"  E1 Method B (permutation) p<0.05: {r['e1_B_sig']}  p={r['e1_p_B']}  (pre-declared: True)")
    print(f"  E1 Method C (wild) excl 0: {r['e1_C_excl0']}  (pre-declared: True)")
    print(f"  E1 wild wider than percentile: {r['e1_wild_wider']}  (pre-declared: True)")
    # E2
    print(f"  E2 Method A (percentile) excl 0: {r['e2_A_excl0']}  (pre-declared: True)")
    print(f"  E2 Method B (permutation) p<0.05: {r['e2_B_sig']}  p={r['e2_p_B']}  (pre-declared: True)")
    print(f"  E2 Method C (wild) excl 0: {r['e2_C_excl0']}  (pre-declared: True)")
    print(f"  E2 wild wider than percentile: {r['e2_wild_wider']}  (pre-declared: True)")

# Stability check
print(f"\nStability (med_n_sens_gene >= {STABILITY_THRESH}):")
for r in results:
    flag = "[UNSTABLE]" if not r['stable'] else "OK"
    print(f"  thr={r['thr']:5.2f}  med_n_sens={r['med_n_sens_gene']}  {flag}")

print(f"\nThreshold selection:")
print(f"  Chosen: -1.0  (GDSC waterfall convention, ~17% sensitive per drug, Iorio et al. 2016 PMID 27397505)")
print(f"  Guard: -1.0 passes 5-40% window at {next(r['guard_pct'] for r in results if r['thr']==-1.0):.0f}%")
print(f"  Note: -0.5 also passes guard ({next(r['guard_pct'] for r in results if r['thr']==-0.5):.0f}%).")
print(f"  Selection justified by convention, NOT by delta. Full sweep is reported.")
print(f"  E1 delta(-0.5)={next(r['e1_delta'] for r in results if r['thr']==-0.5):+.4f}  "
      f"E1 delta(-1.0)={next(r['e1_delta'] for r in results if r['thr']==-1.0):+.4f}")
print(f"  E2 delta(-0.5)={next(r['e2_delta'] for r in results if r['thr']==-0.5):+.4f}  "
      f"E2 delta(-1.0)={next(r['e2_delta'] for r in results if r['thr']==-1.0):+.4f}")
