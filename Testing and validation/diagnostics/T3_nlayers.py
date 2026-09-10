"""
T3 — D21: Is the protein arm live?
Run from repo root: python diagnostics/T3_nlayers.py
"""
import hashlib, sys
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import stats

REPO = Path(__file__).resolve().parent.parent.parent
OUT  = REPO / "Testing and validation" / "diagnostics" / "out"
sys.path.insert(0, str(REPO / "architecture"))
from config import CORE_SCORE, RNA_Z

def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()

print("=" * 72)
print("T3 — D21: IS THE PROTEIN ARM LIVE?")
print("=" * 72)

print(f"[PROVENANCE] core_score sha256 = {sha256(CORE_SCORE)}")
core = pd.read_parquet(CORE_SCORE)
print(f"[COUNT] core_score total rows: {len(core):,}")
print()

# 1. n_layers value counts
vc = core["n_layers"].value_counts().sort_index()
total = len(core)
print("[MEASURED] n_layers.value_counts():")
for v, n in vc.items():
    print(f"  n_layers={v}: {n:,}  ({100*n/total:.2f}%)")
print()

# 2. Distinct lines per branch
n1_lines = core[core.n_layers == 1]["model_id"].nunique()
n2_lines = core[core.n_layers == 2]["model_id"].nunique()
print(f"[MEASURED] Distinct lines in n_layers=1 branch: {n1_lines:,}")
print(f"[MEASURED] Distinct lines in n_layers=2 branch: {n2_lines:,}")
print()

# 3. Verify protein genuinely contributed to n_layers=2 rows
# For n_layers=2, core_score should NOT equal Phi(rna_pct_z)
# Strategy: recompute rna_pct_z for a sample and compare
print("[MEASURED] Verifying protein contribution to n_layers=2 ...")
rna = pd.read_parquet(RNA_Z).rename(columns={"z_t": "rna_z", "gene_id": "ensg_id"})

# Get lineage from DuckDB
import duckdb
from config import DB
con = duckdb.connect(str(DB), read_only=True)
lineage_map = con.execute(
    "SELECT model_id, lineage FROM main.sample_info WHERE lineage IS NOT NULL"
).df().drop_duplicates("model_id").set_index("model_id")["lineage"]
con.close()

# Sample 1000 n_layers=2 rows
sample_n2 = core[core.n_layers == 2].sample(1000, random_state=42)

# Recompute rna_pct_z for these (model, gene) pairs
rna_wide = rna.pivot(index="model_id", columns="ensg_id", values="rna_z")
rna_wide["lineage"] = lineage_map.reindex(rna_wide.index)
rna_pct = rna_wide.groupby("lineage", group_keys=False).apply(
    lambda g: g.drop(columns=["lineage"]).rank(pct=True, ascending=True, na_option="keep"),
    include_groups=False
)
rna_pct_z_wide = pd.DataFrame(
    stats.norm.ppf(np.clip(rna_pct.values, 0.001, 0.999)),
    index=rna_pct.index, columns=rna_pct.columns
)
rna_only_score_wide = pd.DataFrame(
    stats.norm.cdf(rna_pct_z_wide.values),
    index=rna_pct_z_wide.index, columns=rna_pct_z_wide.columns
)

matches = 0
checked = 0
for _, row in sample_n2.iterrows():
    m, g = row["model_id"], row["ensg_id"]
    if m in rna_only_score_wide.index and g in rna_only_score_wide.columns:
        rna_only = rna_only_score_wide.loc[m, g]
        actual   = row["core_score"]
        if np.isclose(actual, rna_only, atol=1e-5):
            matches += 1
        checked += 1

print(f"  Checked: {checked:,} n_layers=2 pairs")
print(f"  Pairs where core_score == Phi(rna_pct_z) exactly: {matches:,}")
print(f"  Pairs where protein genuinely contributed: {checked - matches:,}")
if checked > 0:
    print(f"  Fraction with protein contribution: {(checked-matches)/checked:.4f}")
print()

# 4. grep docs/ for protein arm exclusion claims
print("[DOCUMENTED] Searching docs/ for protein-exclusion statements:")
docs_dir = REPO / "architecture" / "docs"
if docs_dir.exists():
    for p in docs_dir.rglob("*.md"):
        for i, line in enumerate(p.read_text(encoding="utf-8", errors="ignore").splitlines(), 1):
            if any(kw in line.lower() for kw in ["protein arm excluded", "protein excluded",
                                                   "protein not used", "protein disabled",
                                                   "n_layers=1 only", "nlayers=1 only"]):
                print(f"  {p.relative_to(REPO)}:{i}: {line.strip()}")
else:
    print("  docs/ directory does not exist")
print()

# Verdict
print("=" * 72)
n2_count = vc.get(2, 0)
frac_n2  = n2_count / total
if n2_count > 0 and frac_n2 > 0.50:
    print(f"VERDICT: EXPECTED — n_layers=2 IS LIVE: {n2_count:,} pairs ({100*frac_n2:.1f}% of total)")
    print("  Any docs claiming protein arm is excluded are STALE.")
else:
    print(f"VERDICT: UNEXPECTED — n_layers=2 NOT dominant: {n2_count:,} pairs ({100*frac_n2:.1f}%)")
print("=" * 72)
