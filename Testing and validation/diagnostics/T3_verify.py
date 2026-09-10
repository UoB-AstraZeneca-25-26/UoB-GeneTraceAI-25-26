"""T3 step 2 — verify protein genuinely contributed to n_layers=2 rows."""
import sys
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import stats

REPO = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(REPO / "architecture"))
from config import CORE_SCORE, RNA_Z, DB
import duckdb

core = pd.read_parquet(CORE_SCORE)
rna  = pd.read_parquet(RNA_Z).rename(columns={"z_t": "rna_z", "gene_id": "ensg_id"})

con = duckdb.connect(str(DB), read_only=True)
lineage_map = con.execute(
    "SELECT model_id, lineage FROM main.sample_info WHERE lineage IS NOT NULL"
).df().drop_duplicates("model_id").set_index("model_id")["lineage"]
con.close()

# Build rna_pct_z the same way core_score.py does
rna_wide = rna.pivot(index="model_id", columns="ensg_id", values="rna_z")
rna_wide_with_lin = rna_wide.copy()
rna_wide_with_lin["lineage"] = lineage_map.reindex(rna_wide_with_lin.index)
rna_wide_with_lin = rna_wide_with_lin.dropna(subset=["lineage"])

# rank within lineage (include_groups=False means lineage is excluded from g)
rna_pct = rna_wide_with_lin.groupby("lineage", group_keys=False).apply(
    lambda g: g.rank(pct=True, ascending=True, na_option="keep"),
    include_groups=False
)
rna_pct_z_wide = pd.DataFrame(
    stats.norm.ppf(np.clip(rna_pct.values, 0.001, 0.999)),
    index=rna_pct.index, columns=rna_pct.columns
)
rna_only_core_wide = pd.DataFrame(
    stats.norm.cdf(rna_pct_z_wide.values),
    index=rna_pct_z_wide.index, columns=rna_pct_z_wide.columns
)

# Sample 1000 n_layers=2 rows
sample_n2 = core[core.n_layers == 2].sample(1000, random_state=42)

matches = 0
checked = 0
diffs   = []
for _, row in sample_n2.iterrows():
    m, g = row["model_id"], row["ensg_id"]
    if m in rna_only_core_wide.index and g in rna_only_core_wide.columns:
        rna_only = rna_only_core_wide.loc[m, g]
        actual   = row["core_score"]
        diff     = abs(actual - rna_only)
        diffs.append(diff)
        if diff < 1e-4:
            matches += 1
        checked += 1

print(f"[MEASURED] Checked: {checked:,} n_layers=2 pairs against pure-RNA score")
print(f"[MEASURED] Pairs where core_score ≈ Phi(rna_pct_z)  (diff < 1e-4): {matches:,}")
print(f"[MEASURED] Pairs where protein genuinely changed score:               {checked - matches:,}")
if checked > 0:
    print(f"[MEASURED] Fraction with protein contribution: {(checked-matches)/checked:.4f}")
    print(f"[MEASURED] Mean abs diff (core vs RNA-only):   {np.mean(diffs):.4f}")
    print(f"[MEASURED] Median abs diff:                    {np.median(diffs):.4f}")
    print(f"[MEASURED] Max abs diff:                       {np.max(diffs):.4f}")
