"""
V2_V3_V6_checks.py
V2: Gene counts per stratum
V3: Reconcile BEFORE numbers against ALTERATION_DEFENCE.md
V6: both stratum — hr_driver = 0.000 investigation

Run: python diagnostics/V2_V3_V6_checks.py
"""
import sys
from pathlib import Path
import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(REPO / "architecture"))
from config import FLAGS_DRIVER, VALIDATION, GENE_LKP, PREDICTIONS

HOLD_OUT_SEED  = 42
HOLD_OUT_FRAC  = 0.20
HITS_AT        = 20
GLOBAL_IC50_THRESHOLD = 1.0
SENSITIVITY_MAD_Z = -0.5


def label_global(gdsc, ic50_col):
    gdsc = gdsc.copy()
    gdsc["is_sensitive"] = gdsc[ic50_col] <= GLOBAL_IC50_THRESHOLD
    return gdsc


def label_per_drug(gdsc, ic50_col, drug_col):
    def _z(grp):
        med = grp[ic50_col].median()
        mad = (grp[ic50_col] - med).abs().median()
        mad_s = mad * 1.4826 if mad > 0 else grp[ic50_col].std(ddof=1)
        mad_s = max(mad_s, 1e-6)
        g2 = grp.copy()
        g2["drug_z"] = (grp[ic50_col] - med) / mad_s
        g2["is_sensitive"] = g2["drug_z"] <= SENSITIVITY_MAD_Z
        return g2
    return gdsc.groupby(drug_col, group_keys=False).apply(_z, include_groups=False)


def compute_ranks(df):
    df = df.copy()
    df["rank_flat"] = df.groupby("ensg_id")["core_score"].rank(ascending=False, method="first")
    df["rank_driver"] = df.groupby("ensg_id", group_keys=False).apply(
        lambda g: g["core_score"].where(g["has_driver_alteration"], other=-1)
                                 .rank(ascending=False, method="first"),
        include_groups=False
    )
    return df


def hits_at_k(df, rank_col, k=HITS_AT):
    in_top = df[df[rank_col] <= k]
    return in_top["is_sensitive"].sum() / max(df["is_sensitive"].sum(), 1)


flags    = pd.read_parquet(FLAGS_DRIVER)
gdsc_raw = pd.read_parquet(VALIDATION / "gdsc_scored_ready.parquet")
genes    = pd.read_parquet(GENE_LKP, columns=["ensg_id","gene_role"])

ic50_col = "LN_IC50" if "LN_IC50" in gdsc_raw.columns else "log_ic50"
drug_col = "drug_name" if "drug_name" in gdsc_raw.columns else gdsc_raw.columns[0]

curated = sorted(set(flags.ensg_id) & set(gdsc_raw.target_ensg))
rng_split = np.random.default_rng(HOLD_OUT_SEED)
test_genes = set(rng_split.choice(curated, size=int(len(curated) * HOLD_OUT_FRAC), replace=False))

# ── V2 — gene counts per stratum ──────────────────────────────────────────────
print("=" * 70)
print("V2 — Gene counts per stratum")
print("=" * 70)

core_t = flags[flags.ensg_id.isin(test_genes)].copy()
gdsc_t = gdsc_raw[gdsc_raw.target_ensg.isin(test_genes)].copy().rename(columns={"target_ensg":"ensg_id"})
gdsc_t["model_id"] = gdsc_t["model_id"].str.lower()
gdsc_t_madz = label_per_drug(gdsc_t, ic50_col, drug_col)

core_t = core_t.merge(
    gdsc_t_madz[["ensg_id","model_id","is_sensitive"]].drop_duplicates(["ensg_id","model_id"]),
    on=["ensg_id","model_id"], how="left"
)
core_t["is_sensitive"] = core_t["is_sensitive"].fillna(False)
core_t = core_t.merge(genes, on="ensg_id", how="left")
core_t = compute_ranks(core_t)

print(f"Total test genes: {core_t.ensg_id.nunique()}")
print(f"Total test pairs: {len(core_t):,}")

for cls in ["oncogene", "tsg", "both", None]:
    if cls is None:
        sub = core_t[core_t.gene_role.isna()]
        label = "no_role (excluded from strata)"
    else:
        sub = core_t[core_t.gene_role == cls]
        label = cls

    n_genes   = sub.ensg_id.nunique()
    n_sens    = int(sub["is_sensitive"].sum())
    median_per_gene = sub.groupby("ensg_id")["is_sensitive"].sum().median() if n_genes > 0 else 0
    n_genes_with_sensitive = int((sub.groupby("ensg_id")["is_sensitive"].sum() >= 1).sum())
    n_genes_driver_top20 = int(
        sub[sub["rank_driver"] <= HITS_AT].groupby("ensg_id")["has_driver_alteration"]
        .any().sum()
    ) if n_genes > 0 else 0
    interpretable = "" if n_genes >= 5 else f"  ← NOT INTERPRETABLE (n={n_genes})"
    print(f"\n  {label}:")
    print(f"    n_genes                        = {n_genes}{interpretable}")
    print(f"    n_genes with ≥1 sensitive line = {n_genes_with_sensitive}")
    print(f"    n_sensitive_lines (total)       = {n_sens}")
    print(f"    n_sensitive_lines (median/gene) = {median_per_gene:.1f}")
    print(f"    n_genes with ≥1 driver in top20 = {n_genes_driver_top20}")

print(f"\nINVARIANT CHECK: sum of stratum counts = "
      f"{sum(core_t[core_t.gene_role==c].ensg_id.nunique() for c in ['oncogene','tsg','both'])} "
      f"+ {core_t[core_t.gene_role.isna()].ensg_id.nunique()} (no_role) "
      f"= {core_t.ensg_id.nunique()} total")

# ── V3 — reconcile BEFORE ─────────────────────────────────────────────────────
print("\n" + "=" * 70)
print("V3 — Reconcile BEFORE against ALTERATION_DEFENCE.md")
print("=" * 70)
print("Committed: hr_flat=0.014025  hr_driver=0.017735  (28 genes, global LN_IC50<=1.0)")
print("Recent BEFORE: hr_flat=0.006  hr_driver=0.014  (oncogene stratum, 11 genes)")

# Reproduce with global label
gdsc_t_global = label_global(gdsc_t, ic50_col)
core_v3 = flags[flags.ensg_id.isin(test_genes)].copy()
core_v3 = core_v3.merge(
    gdsc_t_global[["ensg_id","model_id","is_sensitive"]].drop_duplicates(["ensg_id","model_id"]),
    on=["ensg_id","model_id"], how="left"
)
core_v3["is_sensitive"] = core_v3["is_sensitive"].fillna(False)
core_v3 = core_v3.merge(genes, on="ensg_id", how="left")
core_v3 = compute_ranks(core_v3)

print(f"\nReproduced (global LN_IC50<=1.0, seed=42):")
print(f"  n_test_genes = {core_v3.ensg_id.nunique()}")
for cls in ["oncogene", "tsg", "both"]:
    sub = core_v3[core_v3.gene_role == cls]
    if sub.empty: continue
    hr_f = hits_at_k(sub, "rank_flat")
    hr_d = hits_at_k(sub, "rank_driver")
    print(f"  {cls} (n={sub.ensg_id.nunique()}): flat={hr_f:.6f}  driver={hr_d:.6f}")

# Pooled across all strata
hr_f_all = hits_at_k(core_v3, "rank_flat")
hr_d_all = hits_at_k(core_v3, "rank_driver")
print(f"  POOLED all roles: flat={hr_f_all:.6f}  driver={hr_d_all:.6f}")
print(f"  Committed:        flat=0.014025        driver=0.017735")
print(f"  Match (pooled):   flat={'YES' if abs(hr_f_all-0.014025)<0.001 else 'NO'}  "
      f"driver={'YES' if abs(hr_d_all-0.017735)<0.001 else 'NO'}")
print(f"\n  HYPOTHESIS D CHECK: oncogene flat=0.006, tsg flat=0.014 → pooled flat={hr_f_all:.4f}")
print(f"  0.014025 sits between 0.006 (oncogene, n=11) and 0.014 (tsg, n=1) — "
      f"pooled accounts for more genes per stratum")

# ── V6 — both stratum hr_driver = 0.000 ──────────────────────────────────────
print("\n" + "=" * 70)
print("V6 — 'both' stratum: hr_driver = 0.000")
print("=" * 70)

both_genes = core_t[core_t.gene_role == "both"].ensg_id.unique()
if len(both_genes) == 0:
    print("  No 'both' genes in test set — cannot investigate")
else:
    for gene_id in both_genes:
        sub = core_t[core_t.ensg_id == gene_id].copy()
        # Try to get gene name
        gname = genes[genes.ensg_id == gene_id]["gene_role"].values
        # Get from GENE_LKP for symbol
        try:
            glk = pd.read_parquet(REPO / "reference" / "gene_lookup.parquet",
                                  columns=["ensg_id","hgnc_symbol","gene_role"])
            sym = glk[glk.ensg_id == gene_id]["hgnc_symbol"].values
            sym = sym[0] if len(sym) > 0 else gene_id
        except Exception:
            sym = gene_id

        print(f"\n  Gene: {sym} ({gene_id})  gene_role=both  n_lines={len(sub)}")
        n_driver = int(sub["has_driver_alteration"].sum())
        n_sensitive = int(sub["is_sensitive"].sum())
        print(f"  has_driver_alteration=TRUE:  {n_driver}")
        print(f"  is_sensitive=TRUE:           {n_sensitive}")

        # Triage: step 1 — is n_driver == 0?
        if n_driver == 0:
            print(f"  TRIAGE: n_driver=0 → ARTEFACT (join/coverage issue)")
            # Check which alteration layers have rows for this gene
            mut_path = REPO / "architecture" / "outputs" / "mutations_scores.parquet"
            fus_path = REPO / "architecture" / "outputs" / "fusions_scores.parquet"
            cna_path = REPO / "architecture" / "outputs" / "cna_flags.parquet"
            for label, path in [("mutations", mut_path), ("fusions", fus_path), ("CNA", cna_path)]:
                if path.exists():
                    try:
                        df = pd.read_parquet(path)
                        col = "ensg_id" if "ensg_id" in df.columns else df.columns[0]
                        n = (df[col].str.upper() == gene_id.upper()).sum()
                        print(f"    {label}: {n} rows for this gene")
                    except Exception as e:
                        print(f"    {label}: could not read ({e})")
        else:
            # Top 20 lines by core_score
            top20 = sub.nlargest(20, "core_score")[
                ["model_id","core_score","has_driver_alteration","is_sensitive",
                 "p_mutation","p_fusion","has_cna_alteration","rank_flat","rank_driver"]
            ]
            print(f"\n  Top 20 by core_score:")
            print(top20.to_string(index=False))

            # Where do driver-flagged lines rank?
            driver_ranks = sub[sub["has_driver_alteration"]]["rank_flat"].sort_values()
            print(f"\n  Driver-flagged line ranks (rank_flat): {driver_ranks.tolist()[:20]}")
            print(f"  Driver-flagged lines in top 20 (rank_flat): {(driver_ranks <= 20).sum()}")

            # Where do sensitive lines rank?
            sens_ranks = sub[sub["is_sensitive"]]["rank_flat"].sort_values()
            print(f"  Sensitive line ranks (rank_flat): {sens_ranks.tolist()[:20]}")

            # What are driver-line scores relative to sensitive lines?
            print(f"\n  TRIAGE: n_driver={n_driver} > 0")
            driver_in_top20 = sub[sub["has_driver_alteration"] & (sub["rank_driver"] <= 20)]
            print(f"  Driver lines in driver-gated top 20: {len(driver_in_top20)}")
            if len(driver_in_top20) > 0:
                driver_sens = driver_in_top20["is_sensitive"].sum()
                print(f"  Sensitive among driver-gated top 20: {driver_sens}")
                if driver_sens == 0:
                    print(f"  → Gate selected correctly; drug label disagrees (non-discriminating compound?)")
            else:
                print(f"  → No driver lines in top 20 even with driver gate")
                # Check: are driver lines anti-correlated with sensitivity?
                driver_lines = sub[sub["has_driver_alteration"]]
                print(f"    Driver lines: n={len(driver_lines)}, "
                      f"mean core_score={driver_lines.core_score.mean():.4f}, "
                      f"mean rank_flat={driver_lines.rank_flat.mean():.1f}")
                non_driver = sub[~sub["has_driver_alteration"]]
                print(f"    Non-driver lines: n={len(non_driver)}, "
                      f"mean core_score={non_driver.core_score.mean():.4f}, "
                      f"mean rank_flat={non_driver.rank_flat.mean():.1f}")
