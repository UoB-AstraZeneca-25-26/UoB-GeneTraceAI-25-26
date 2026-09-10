"""
W1_abundance_tracking.py
Characterise the gene_role value(s) that caused 13/26 test genes to be
excluded from eval.py's stratum loop.

Run: python diagnostics/W1_abundance_tracking.py
"""
import sys
from pathlib import Path
import pandas as pd
import numpy as np

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "final_pipeline"))
from config import GENE_LKP, FLAGS_DRIVER, REPO as _

print("=" * 70)
print("W1 — Characterise the 13 excluded test genes")
print("=" * 70)

# 1. Distinct gene_role values in the full gene lookup
genes = pd.read_parquet(GENE_LKP)
print(f"\nColumns in gene_lookup.parquet: {list(genes.columns)}")
print(f"Total genes in gene_lookup: {len(genes):,}")
role_counts = genes["gene_role"].value_counts(dropna=False)
print(f"\n1. gene_role value counts (all genes):")
for val, n in role_counts.items():
    print(f"   {repr(val):30s}  {n:,}")

# 2. Identify the 13 excluded test genes (same holdout as eval.py)
flags = pd.read_parquet(FLAGS_DRIVER)
import duckdb
con = duckdb.connect(str(REPO / "src/pipeline/outputs/celllineselector.db"), read_only=True)
try:
    gdsc_raw = None
    gdsc_path = REPO / "validation" / "prepared" / "gdsc_scored_ready.parquet"
    if gdsc_path.exists():
        gdsc_raw = pd.read_parquet(gdsc_path)
    con.close()
except Exception as e:
    con.close()
    print(f"  [DB warning: {e}]")

if gdsc_raw is not None:
    curated = sorted(set(flags.ensg_id) & set(gdsc_raw.target_ensg))
    rng = np.random.default_rng(42)
    test_genes_set = set(rng.choice(curated, size=int(len(curated) * 0.20), replace=False))
    core_t = flags[flags.ensg_id.isin(test_genes_set)].copy()
    core_t = core_t.merge(genes[["ensg_id", "gene_role"]], on="ensg_id", how="left")

    print(f"\n2. Test gene set: {len(test_genes_set)} genes")
    test_gene_roles = core_t.drop_duplicates("ensg_id")[["ensg_id", "gene_role"]]
    role_in_test = test_gene_roles["gene_role"].value_counts(dropna=False)
    print(f"   gene_role distribution in test genes:")
    for val, n in role_in_test.items():
        print(f"   {repr(val):30s}  {n}")

    known_roles = {"oncogene", "tsg", "both"}
    excluded_test = test_gene_roles[~test_gene_roles["gene_role"].isin(known_roles)]
    print(f"\n3. Excluded genes (role not in {{oncogene, tsg, both}}):")
    print(f"   n = {len(excluded_test)}")
    print(f"   Excluded gene_role values: {excluded_test['gene_role'].value_counts(dropna=False).to_dict()}")
    print(f"   Excluded ensg_ids: {sorted(excluded_test['ensg_id'].tolist())}")

    # 4. COSMIC CGC check
    # gene_lookup might have a 'cgc' or 'is_cgc' or 'tier' column
    cgc_cols = [c for c in genes.columns if "cgc" in c.lower() or "tier" in c.lower() or "cosmic" in c.lower()]
    print(f"\n4. CGC-related columns in gene_lookup: {cgc_cols}")
    if cgc_cols:
        col = cgc_cols[0]
        excluded_full = genes[~genes["gene_role"].isin(known_roles) & genes["gene_role"].notna()]
        print(f"   Excluded genes (project-wide) with {col}:")
        print(f"   {excluded_full[col].value_counts(dropna=False).to_dict()}")
        cgc_overlap = excluded_test.merge(genes[["ensg_id", col]], on="ensg_id", how="left")
        print(f"   CGC status of 13 excluded test genes:")
        print(cgc_overlap[["ensg_id", col]].to_string(index=False))
    else:
        print("   No CGC column found in gene_lookup — checking for 'hgnc_symbol' to cross-ref")
        if "hgnc_symbol" in genes.columns:
            excluded_syms = genes[genes["ensg_id"].isin(excluded_test["ensg_id"])]["hgnc_symbol"].tolist()
            print(f"   Excluded gene symbols: {excluded_syms}")

    # 5. Which alteration channels can fire for excluded genes?
    print(f"\n5. Alteration channel audit for excluded test genes:")
    excl_ids = set(excluded_test["ensg_id"])
    core_excl = core_t[core_t["ensg_id"].isin(excl_ids)]

    for col in ["has_driver_alteration", "has_cna_alteration", "p_mutation", "p_fusion"]:
        if col in core_excl.columns:
            if core_excl[col].dtype == bool or col.startswith("has_"):
                n_true = int(core_excl[col].astype(bool).sum())
                print(f"   {col} = TRUE: {n_true:,}  ({100*n_true/max(len(core_excl),1):.1f}%)")
            else:
                n_fire = int((core_excl[col] >= 0.5).sum())
                print(f"   {col} >= 0.5: {n_fire:,}  ({100*n_fire/max(len(core_excl),1):.1f}%)")
        else:
            print(f"   {col}: column not found in flags_with_driver")

    # 5b. Check if excluded genes can receive CNA alteration calls
    cna_col = "has_cna_alteration"
    if cna_col in core_excl.columns:
        n_cna = int(core_excl[cna_col].astype(bool).sum())
        pre_declared = "PASS" if n_cna == 0 else "FAIL — pre-declared 0"
        print(f"\n   CNA PRE-DECLARED CHECK: n_cna_true={n_cna} (expect 0)  → {pre_declared}")
    else:
        print(f"\n   {cna_col} column not found")

else:
    print("  [GDSC file not found — cannot reproduce test gene split]")

# 6. Exact exclusion line in eval.py
print(f"\n6. Exclusion in eval.py:")
eval_path = REPO / "final_pipeline" / "Validation" / "eval.py"
if eval_path.exists():
    lines = eval_path.read_text(encoding="utf-8").splitlines()
    for i, line in enumerate(lines, 1):
        if "for cls in" in line and ("oncogene" in line or "tsg" in line or "both" in line):
            print(f"   Line {i}: {line.strip()}")
        if "gene_role" in line and "==" in line:
            print(f"   Line {i}: {line.strip()}")

print(f"\n7. Is the exclusion logged anywhere?")
print(f"   Searching for '[EXCLUSION]' or similar log lines in eval.py...")
excl_logged = False
for i, line in enumerate(lines, 1):
    if "exclusion" in line.lower() or "abundance" in line.lower() or "dropping" in line.lower():
        print(f"   Line {i}: {line.strip()}")
        excl_logged = True
if not excl_logged:
    print("   NOT LOGGED — confirmed silent exclusion")

print(f"\nPRE-DECLARED CHECK:")
if gdsc_raw is not None:
    excl_ids_check = set(excluded_test["ensg_id"])
    core_excl_check = core_t[core_t["ensg_id"].isin(excl_ids_check)]
    if "has_cna_alteration" in core_excl_check.columns:
        n_cna = int(core_excl_check["has_cna_alteration"].astype(bool).sum())
        print(f"  CNA true for excluded genes: {n_cna}  (pre-declared: 0)  → {'PASS' if n_cna==0 else 'FAIL'}")
print(f"  Exclusion logged: {excl_logged}  (pre-declared: False)  → {'FAIL' if excl_logged else 'PASS'}")
