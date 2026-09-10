"""
V4_threshold_sweep.py
MAD-z threshold sweep: -0.25, -0.5, -1.0, -1.5, -2.0
1000 resamples per threshold (paired bootstrap clustered on gene).
Reports ALL rows; no threshold recommended.

Run: python diagnostics/V4_threshold_sweep.py
"""
import sys
from pathlib import Path
import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(REPO / "architecture"))
from config import FLAGS_DRIVER, VALIDATION, GENE_LKP

HOLD_OUT_SEED = 42
HOLD_OUT_FRAC = 0.20
HITS_AT       = 20
N_BOOT        = 1000
THRESHOLDS    = [-0.25, -0.5, -1.0, -1.5, -2.0]

print("=" * 70)
print("V4 — MAD-z threshold sweep")
print("=" * 70)

flags    = pd.read_parquet(FLAGS_DRIVER)
gdsc_raw = pd.read_parquet(VALIDATION / "gdsc_scored_ready.parquet")
genes    = pd.read_parquet(GENE_LKP, columns=["ensg_id","gene_role"])

ic50_col = "LN_IC50" if "LN_IC50" in gdsc_raw.columns else "log_ic50"
drug_col = "drug_name" if "drug_name" in gdsc_raw.columns else gdsc_raw.columns[0]

curated   = sorted(set(flags.ensg_id) & set(gdsc_raw.target_ensg))
rng_split = np.random.default_rng(HOLD_OUT_SEED)
test_genes = set(rng_split.choice(curated, size=int(len(curated) * HOLD_OUT_FRAC), replace=False))

# Pre-filter to test genes for efficiency
flags_t = flags[flags.ensg_id.isin(test_genes)].copy()
gdsc_t  = gdsc_raw[gdsc_raw.target_ensg.isin(test_genes)].copy()
gdsc_t  = gdsc_t.rename(columns={"target_ensg":"ensg_id"})
gdsc_t["model_id"] = gdsc_t["model_id"].str.lower()

# Pre-compute drug median and MAD for speed
drug_stats = gdsc_t.groupby(drug_col)[ic50_col].agg(
    med="median",
    mad=lambda x: (x - x.median()).abs().median()
).reset_index()
drug_stats["mad_s"] = (drug_stats["mad"] * 1.4826).clip(lower=1e-6)
gdsc_t = gdsc_t.merge(drug_stats[[drug_col,"med","mad_s"]], on=drug_col, how="left")
gdsc_t["drug_z"] = (gdsc_t[ic50_col] - gdsc_t["med"]) / gdsc_t["mad_s"]

# Compute ranks once (independent of sensitivity label)
core_t_base = flags_t.copy()
core_t_base["rank_flat"] = (
    core_t_base.groupby("ensg_id")["core_score"].rank(ascending=False, method="first")
)
core_t_base["rank_driver"] = (
    core_t_base.groupby("ensg_id", group_keys=False).apply(
        lambda g: g["core_score"].where(g["has_driver_alteration"], other=-1)
                                 .rank(ascending=False, method="first"),
        include_groups=False
    )
)

results = []

for thr in THRESHOLDS:
    # Label sensitivity at this threshold
    gdsc_labeled = gdsc_t[["ensg_id","model_id","drug_z"]].copy()
    gdsc_labeled["is_sensitive"] = gdsc_labeled["drug_z"] <= thr

    core_t = core_t_base.merge(
        gdsc_labeled[["ensg_id","model_id","is_sensitive"]].drop_duplicates(["ensg_id","model_id"]),
        on=["ensg_id","model_id"], how="left"
    )
    core_t["is_sensitive"] = core_t["is_sensitive"].fillna(False)

    # Guard: % drugs in 5-40% window
    per_drug_rate = gdsc_t.assign(is_s=gdsc_t["drug_z"] <= thr).groupby(drug_col)["is_s"].mean()
    in_window = ((per_drug_rate >= 0.05) & (per_drug_rate <= 0.40)).mean()
    median_rate = per_drug_rate.median()
    n_sensitive = int(core_t["is_sensitive"].sum())

    # Per-gene contributions
    gene_list = sorted(core_t.ensg_id.unique())
    flat_c   = {}
    driver_c = {}
    for g, grp in core_t.groupby("ensg_id"):
        flat_c[g]   = (int(((grp["rank_flat"]   <= HITS_AT) & grp["is_sensitive"]).sum()),
                       int(grp["is_sensitive"].sum()))
        driver_c[g] = (int(((grp["rank_driver"] <= HITS_AT) & grp["is_sensitive"]).sum()),
                       int(grp["is_sensitive"].sum()))

    def pool_hr(c, gs):
        num = sum(c[g][0] for g in gs); den = sum(c[g][1] for g in gs)
        return num / max(den, 1)

    hr_f = pool_hr(flat_c,   gene_list)
    hr_d = pool_hr(driver_c, gene_list)
    delta_pt  = hr_d - hr_f
    lift_pt   = hr_d / hr_f if hr_f > 0 else float("nan")

    # Bootstrap
    genes_arr = np.array(gene_list)
    n_genes   = len(gene_list)
    rng_b     = np.random.default_rng(HOLD_OUT_SEED)
    deltas    = np.zeros(N_BOOT)
    lifts     = np.full(N_BOOT, np.nan)
    for i in range(N_BOOT):
        samp = rng_b.choice(genes_arr, size=n_genes, replace=True)
        fn = sum(flat_c[g][0] for g in samp);   fd = sum(flat_c[g][1] for g in samp)
        dn = sum(driver_c[g][0] for g in samp); dd = sum(driver_c[g][1] for g in samp)
        hf = fn / max(fd, 1); hd = dn / max(dd, 1)
        deltas[i] = hd - hf
        if hf > 0: lifts[i] = hd / hf

    ci_lo = float(np.percentile(deltas, 2.5))
    ci_hi = float(np.percentile(deltas, 97.5))
    excludes_zero = (ci_lo > 0) or (ci_hi < 0)

    results.append({
        "threshold": thr,
        "guard_%_in_window": round(in_window*100, 1),
        "median_%_sensitive": round(median_rate*100, 1),
        "n_sensitive": n_sensitive,
        "hr_flat":  round(hr_f, 6),
        "hr_driver": round(hr_d, 6),
        "delta": round(delta_pt, 6),
        "lift": round(lift_pt, 4) if not np.isnan(lift_pt) else "undef",
        "ci_95_lo": round(ci_lo, 6),
        "ci_95_hi": round(ci_hi, 6),
        "excludes_zero": excludes_zero,
    })
    print(f"  thr={thr:5.2f}  guard={in_window*100:.0f}%  "
          f"med_sens={median_rate*100:.0f}%  n_sens={n_sensitive:,}  "
          f"hr_f={hr_f:.4f}  hr_d={hr_d:.4f}  delta={delta_pt:+.4f}  "
          f"CI=[{ci_lo:+.4f},{ci_hi:+.4f}]  excl_0={excludes_zero}")

print(f"\nFull table:")
df_res = pd.DataFrame(results)
print(df_res.to_string(index=False))

# Is the parameter load-bearing?
deltas_list = [r["delta"] for r in results]
delta_range = max(deltas_list) - min(deltas_list)
print(f"\nParameter load-bearing check:")
print(f"  delta range across sweep: {delta_range:.4f}  "
      f"({'LOAD-BEARING (>2x span)' if delta_range > abs(deltas_list[0])*2 else 'not load-bearing'} )")
# Peak at -0.5?
delta_at_chosen = [r["delta"] for r in results if r["threshold"] == -0.5][0]
peak_thr = results[deltas_list.index(max(deltas_list))]["threshold"]
print(f"  delta at chosen threshold (-0.5): {delta_at_chosen:+.4f}")
print(f"  peak delta at threshold: {peak_thr}")
if peak_thr == -0.5:
    print(f"  PARAMETER SELECTION CONCERN: delta peaks at the chosen threshold (-0.5)")
else:
    print(f"  No selection concern — peak is at {peak_thr}, not at -0.5")
