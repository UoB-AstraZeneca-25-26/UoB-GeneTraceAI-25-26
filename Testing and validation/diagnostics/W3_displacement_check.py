"""
W3_displacement_check.py
Displacement audit for the oncogene stratum (n=11 genes):
  displaced  = sensitive lines in flat top-20 NOT in driver-gated top-20
  promoted   = sensitive lines in driver-gated top-20 NOT in flat top-20
  net        = promoted - displaced (> 0 means gate helped)

Also checks the boundary zone (flat ranks 15-25) and the alteration type of
promoted lines where net < 0.

Run: python diagnostics/W3_displacement_check.py
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
BOUNDARY_LO       = 15
BOUNDARY_HI       = 25

print("=" * 70)
print("W3 — Displacement check, oncogene stratum")
print("=" * 70)

flags    = pd.read_parquet(FLAGS_DRIVER)
gdsc_raw = pd.read_parquet(VALIDATION / "gdsc_scored_ready.parquet")
genes    = pd.read_parquet(GENE_LKP, columns=["ensg_id", "gene_role"])

ic50_col = "LN_IC50" if "LN_IC50" in gdsc_raw.columns else "log_ic50"
drug_col = "drug_name" if "drug_name" in gdsc_raw.columns else gdsc_raw.columns[0]

curated   = sorted(set(flags.ensg_id) & set(gdsc_raw.target_ensg))
rng_split = np.random.default_rng(HOLD_OUT_SEED)
test_genes = set(rng_split.choice(curated, size=int(len(curated) * HOLD_OUT_FRAC), replace=False))

# Per-drug sensitivity labels
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

core_t["rank_flat"] = core_t.groupby("ensg_id")["core_score"].rank(ascending=False, method="first")
core_t["rank_driver"] = core_t.groupby("ensg_id", group_keys=False).apply(
    lambda g: g["core_score"].where(g["has_driver_alteration"], other=-1)
                             .rank(ascending=False, method="first"),
    include_groups=False
)

# Restrict to oncogene stratum
onco = core_t[core_t["gene_role"] == "oncogene"].copy()
onco_genes = sorted(onco.ensg_id.unique())
print(f"Oncogene stratum: n_genes={len(onco_genes)}")

rows = []
for g in onco_genes:
    sub = onco[onco.ensg_id == g]
    in_flat_top20    = set(sub[sub["rank_flat"]   <= HITS_AT]["model_id"])
    in_driver_top20  = set(sub[sub["rank_driver"] <= HITS_AT]["model_id"])
    sensitive_ids    = set(sub[sub["is_sensitive"]]["model_id"])

    # Boundary zone sensitive lines (flat ranks 15-25, likely to move)
    boundary = sub[sub["rank_flat"].between(BOUNDARY_LO, BOUNDARY_HI)]
    n_boundary_sens = int(boundary["is_sensitive"].sum())

    # Displaced: in flat top 20 AND sensitive AND NOT in driver top 20
    displaced = len(sensitive_ids & in_flat_top20 - in_driver_top20)
    # Promoted: in driver top 20 AND sensitive AND NOT in flat top 20
    promoted  = len(sensitive_ids & in_driver_top20 - in_flat_top20)
    net       = promoted - displaced

    # For negative net genes: what alteration type do promoted non-sensitive lines carry?
    if net < 0:
        promoted_lines = sub[(sub["rank_driver"] <= HITS_AT) &
                             (~sub["model_id"].isin(in_flat_top20)) &
                             (~sub["is_sensitive"])]
        alt_types = {}
        for col in ["has_cna_alteration", "p_mutation", "p_fusion"]:
            if col in promoted_lines.columns:
                if col.startswith("has_"):
                    alt_types[col] = int(promoted_lines[col].astype(bool).sum())
                else:
                    alt_types[col] = int((promoted_lines[col] >= 0.5).sum())
    else:
        alt_types = {}

    gene_role = sub["gene_role"].iloc[0] if len(sub) > 0 else "?"
    print(f"  {g} ({gene_role}): "
          f"n_sens={len(sensitive_ids):3d}  boundary_sens={n_boundary_sens}  "
          f"displaced={displaced}  promoted={promoted}  net={net:+d}")
    if net < 0:
        print(f"    [net<0] promoted non-sensitive lines by channel: {alt_types}")

    rows.append({"ensg_id": g, "n_sensitive": len(sensitive_ids), "boundary_sens": n_boundary_sens,
                 "displaced": displaced, "promoted": promoted, "net": net})

df = pd.DataFrame(rows)

print(f"\n--- Aggregate (oncogene, n=11 genes) ---")
print(f"  Total displaced:  {df['displaced'].sum()}")
print(f"  Total promoted:   {df['promoted'].sum()}")
print(f"  Net:              {df['net'].sum():+d}")
print(f"  n_genes net > 0:  {(df['net'] > 0).sum()}")
print(f"  n_genes net = 0:  {(df['net'] == 0).sum()}")
print(f"  n_genes net < 0:  {(df['net'] < 0).sum()}")
print(f"  n_genes with any boundary_sens: {(df['boundary_sens'] > 0).sum()}")

print(f"\nPRE-DECLARED CHECK:")
print(f"  net > 0 overall: {df['net'].sum() > 0}  (expected: True, since V1 delta is positive)")
neg_count = (df['net'] < 0).sum()
print(f"  n_genes net < 0: {neg_count}  (pre-declared: ≥2 of 11)")
print(f"  PRE-DECLARED {'PASS' if neg_count >= 2 else 'FAIL (fewer than 2 genes with net<0)'}")

print(f"\nFull displacement table:")
print(df.to_string(index=False))
