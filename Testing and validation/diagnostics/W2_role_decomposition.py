"""
W2_role_decomposition.py
Decompose V1 delta by role group:
  Group A: gene_role in {oncogene, tsg, both}
  Group B: gene_role NOT in {oncogene, tsg, both} and NOT NaN (the 13 excluded genes)

Uses the same V1 labels (MAD-z = -0.5, seed=42).

Run: python diagnostics/W2_role_decomposition.py
"""
import sys
from pathlib import Path
import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(REPO / "architecture"))
from config import FLAGS_DRIVER, VALIDATION, GENE_LKP

HOLD_OUT_SEED     = 42
HOLD_OUT_FRAC     = 0.20
HITS_AT           = 20
SENSITIVITY_MAD_Z = -0.5
N_BOOT            = 2000
KNOWN_ROLES       = {"oncogene", "tsg", "both"}

print("=" * 70)
print("W2 — Delta decomposition by gene role group")
print("=" * 70)

flags    = pd.read_parquet(FLAGS_DRIVER)
gdsc_raw = pd.read_parquet(VALIDATION / "gdsc_scored_ready.parquet")
genes    = pd.read_parquet(GENE_LKP, columns=["ensg_id", "gene_role"])

ic50_col = "LN_IC50" if "LN_IC50" in gdsc_raw.columns else "log_ic50"
drug_col = "drug_name" if "drug_name" in gdsc_raw.columns else gdsc_raw.columns[0]

curated   = sorted(set(flags.ensg_id) & set(gdsc_raw.target_ensg))
rng_split = np.random.default_rng(HOLD_OUT_SEED)
test_genes = set(rng_split.choice(curated, size=int(len(curated) * HOLD_OUT_FRAC), replace=False))

# Per-drug sensitivity labels (identical to V1)
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
    gdsc_t[["ensg_id", "model_id", "is_sensitive"]].drop_duplicates(["ensg_id","model_id"]),
    on=["ensg_id","model_id"], how="left"
)
core_t["is_sensitive"] = core_t["is_sensitive"].fillna(False)
core_t = core_t.merge(genes, on="ensg_id", how="left")

# Ranks
core_t["rank_flat"] = core_t.groupby("ensg_id")["core_score"].rank(ascending=False, method="first")
core_t["rank_driver"] = core_t.groupby("ensg_id", group_keys=False).apply(
    lambda g: g["core_score"].where(g["has_driver_alteration"], other=-1)
                             .rank(ascending=False, method="first"),
    include_groups=False
)

# Group assignment
core_t["group"] = core_t["gene_role"].apply(
    lambda r: "A" if r in KNOWN_ROLES else ("B" if pd.notna(r) else "NaN")
)

def bootstrap_ci(flat_c, driver_c, gene_list, n_boot=N_BOOT, seed=HOLD_OUT_SEED):
    genes_arr = np.array(gene_list)
    n = len(gene_list)
    rng = np.random.default_rng(seed)
    deltas = np.zeros(n_boot)
    for i in range(n_boot):
        samp = rng.choice(genes_arr, size=n, replace=True)
        fn = sum(flat_c[g][0] for g in samp);   fd = sum(flat_c[g][1] for g in samp)
        dn = sum(driver_c[g][0] for g in samp); dd = sum(driver_c[g][1] for g in samp)
        deltas[i] = dn/max(dd,1) - fn/max(fd,1)
    return float(np.percentile(deltas, 2.5)), float(np.percentile(deltas, 97.5))

def group_stats(df, label):
    gene_list = sorted(df.ensg_id.unique())
    n_genes = len(gene_list)
    flat_c, driver_c = {}, {}
    for g, grp in df.groupby("ensg_id"):
        flat_c[g]   = (int(((grp["rank_flat"]   <= HITS_AT) & grp["is_sensitive"]).sum()),
                       int(grp["is_sensitive"].sum()))
        driver_c[g] = (int(((grp["rank_driver"] <= HITS_AT) & grp["is_sensitive"]).sum()),
                       int(grp["is_sensitive"].sum()))

    num_f = sum(v[0] for v in flat_c.values());   den_f = sum(v[1] for v in flat_c.values())
    num_d = sum(v[0] for v in driver_c.values()); den_d = sum(v[1] for v in driver_c.values())
    hr_f = num_f / max(den_f, 1)
    hr_d = num_d / max(den_d, 1)
    delta = hr_d - hr_f
    lift  = hr_d / hr_f if hr_f > 0 else float("nan")

    ci_lo, ci_hi = bootstrap_ci(flat_c, driver_c, gene_list)
    excl_0 = ci_lo > 0

    per_gene_deltas = np.array([driver_c[g][0]/max(driver_c[g][1],1) - flat_c[g][0]/max(flat_c[g][1],1)
                                 for g in gene_list])

    n_sens_per_gene = np.array([flat_c[g][1] for g in gene_list])
    median_n_sens = float(np.median(n_sens_per_gene))

    print(f"\n  Group {label}: n_genes={n_genes}  n_sensitive_total={den_f}  "
          f"median_n_sensitive/gene={median_n_sens:.0f}")
    print(f"    hr_flat={hr_f:.6f}  hr_driver={hr_d:.6f}  delta={delta:+.6f}  "
          f"lift={'undef' if np.isnan(lift) else f'{lift:.4f}x'}")
    print(f"    95% CI (cluster bootstrap): [{ci_lo:+.6f}, {ci_hi:+.6f}]  "
          f"excludes_0={excl_0}")
    if n_genes < 5:
        print(f"    WARNING: n_genes={n_genes} < 5 — NOT INTERPRETABLE (CI unreliable)")
    print(f"    per-gene delta: mean={per_gene_deltas.mean():+.4f}  "
          f"median={np.median(per_gene_deltas):+.4f}  "
          f"pos={int((per_gene_deltas>0).sum())}  "
          f"zero={int((per_gene_deltas==0).sum())}  "
          f"neg={int((per_gene_deltas<0).sum())}")

    # B-group extra: which channels drive alterations?
    if label == "B":
        print(f"    Alteration channels (group B):")
        for col, thresh in [("has_cna_alteration", None), ("p_mutation", 0.5), ("p_fusion", 0.5)]:
            if col in df.columns:
                if thresh is None:
                    n = int(df[col].astype(bool).sum())
                else:
                    n = int((df[col] >= thresh).sum())
                print(f"      {col}: {n} rows fire")

    return {"n_genes": n_genes, "n_sensitive": den_f, "hr_flat": hr_f, "hr_driver": hr_d,
            "delta": delta, "ci_lo": ci_lo, "ci_hi": ci_hi, "flat_c": flat_c, "driver_c": driver_c,
            "gene_list": gene_list}

print(f"\nTest genes: {len(test_genes)}")

grp_a = core_t[core_t["group"] == "A"]
grp_b = core_t[core_t["group"] == "B"]
grp_nan = core_t[core_t["group"] == "NaN"]

stats_a = group_stats(grp_a, "A (oncogene/tsg/both)")
stats_b = group_stats(grp_b, "B (other role, n=?)")
if len(grp_nan) > 0:
    print(f"\n  NaN-role genes: n={grp_nan.ensg_id.nunique()} (these were counted separately in V2)")

# Pooled (all 26 genes)
all_genes = sorted(core_t.ensg_id.unique())
flat_all, driver_all = {}, {}
for g in all_genes:
    grp = core_t[core_t.ensg_id == g]
    flat_all[g]   = (int(((grp["rank_flat"]   <= HITS_AT) & grp["is_sensitive"]).sum()),
                     int(grp["is_sensitive"].sum()))
    driver_all[g] = (int(((grp["rank_driver"] <= HITS_AT) & grp["is_sensitive"]).sum()),
                     int(grp["is_sensitive"].sum()))

fn = sum(v[0] for v in flat_all.values());   fd = sum(v[1] for v in flat_all.values())
dn = sum(v[0] for v in driver_all.values()); dd = sum(v[1] for v in driver_all.values())
pooled_delta = dn/max(dd,1) - fn/max(fd,1)

print(f"\n  Pooled (all {len(all_genes)} genes): hr_flat={fn/max(fd,1):.6f}  "
      f"hr_driver={dn/max(dd,1):.6f}  delta={pooled_delta:+.6f}  "
      f"(V1 measured: +0.008880)")

print(f"\n  Decomposition arithmetic check:")
w_a = stats_a["n_sensitive"] / max(stats_a["n_sensitive"] + stats_b["n_sensitive"], 1)
w_b = stats_b["n_sensitive"] / max(stats_a["n_sensitive"] + stats_b["n_sensitive"], 1)
print(f"    n_sensitive: A={stats_a['n_sensitive']}  B={stats_b['n_sensitive']}  total={stats_a['n_sensitive']+stats_b['n_sensitive']}")
print(f"    Weight_A={w_a:.3f}  Weight_B={w_b:.3f}  (pooled ratio weights by n_sensitive)")
approx_pooled = (stats_a["hr_driver"]*w_a + stats_b["hr_driver"]*w_b) - \
                (stats_a["hr_flat"]*w_a + stats_b["hr_flat"]*w_b)
print(f"    Weighted delta (approx): {approx_pooled:+.6f}  vs actual pooled: {pooled_delta:+.6f}")

print(f"\nPRE-DECLARED CHECK: delta_A > delta_B?")
print(f"  delta_A={stats_a['delta']:+.6f}  delta_B={stats_b['delta']:+.6f}  "
      f"→ {'PASS' if stats_a['delta'] > stats_b['delta'] else 'FAIL (unexpected: B >= A)'}")
