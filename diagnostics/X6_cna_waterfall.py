"""
X6 — CNA attrition waterfall from the raw cosmic_cna table.
COWORK v4.0. Read-only. Starts from ALL rows in main.cosmic_cna.
Run from repo root: python diagnostics/X6_cna_waterfall.py
"""
import sys
from pathlib import Path
import numpy as np
import pandas as pd
import duckdb
from scipy import stats

REPO = Path(__file__).resolve().parent.parent
DB   = REPO / "src" / "pipeline" / "outputs" / "celllineselector.db"
GENE_LKP = REPO / "reference" / "gene_lookup.parquet"
CNA_FLAGS = REPO / "final_pipeline" / "outputs" / "cna_flags.parquet"

print("=" * 72)
print("X6 -- CNA ATTRITION WATERFALL FROM SOURCE (COWORK v4.0)")
print("=" * 72)

gl = pd.read_parquet(GENE_LKP, columns=["ensg_id","gene_role"])
print(f"[COUNT] Gene lookup rows: {len(gl):,}")
print(f"[COUNT] Role breakdown: {gl.gene_role.value_counts().to_dict()}")
n_curated = gl[gl.gene_role.isin(["oncogene","tsg","both"])].ensg_id.nunique()
print(f"[COUNT] Genes with role in {{oncogene,tsg,both}}: {n_curated:,}")

con = duckdb.connect(str(DB), read_only=True)

# ── Step 1: ALL rows in cosmic_cna ────────────────────────────────────────────
print("\n--- STEP 1: ALL ROWS in main.cosmic_cna ---")
raw = con.execute("""
    SELECT model_id, lower(gene_id) AS ensg_id, total_cn, cna_call, is_ambiguous
    FROM main.cosmic_cna
    WHERE model_id IS NOT NULL AND gene_id IS NOT NULL
""").df()
print(f"[COUNT] Step 1 raw: {len(raw):,} rows  "
      f"genes: {raw.ensg_id.nunique():,}  models: {raw.model_id.nunique():,}")
print(f"  cna_call values: {raw.cna_call.value_counts().to_dict()}")

# ── Step 2: pipeline lines (from sample_info) ─────────────────────────────────
print("\n--- STEP 2: After sample-ID join to pipeline lines ---")
si = con.execute("SELECT lower(model_id) AS model_id FROM main.sample_info WHERE lineage IS NOT NULL").df()
pipeline_ids = set(si.model_id.unique())
raw["model_id_lower"] = raw["model_id"].str.lower()
step2 = raw[raw["model_id_lower"].isin(pipeline_ids)].copy()
print(f"[COUNT] Step 2 after join: {len(step2):,}  lost: {len(raw)-len(step2):,}")
print(f"  models matched: {step2.model_id.nunique():,} / {len(pipeline_ids):,} pipeline lines")
unmatched_examples = list(set(raw.model_id.unique()) - {m.upper() for m in pipeline_ids})[:10]
print(f"  unmatched sample_id examples (first 10): {unmatched_examples[:10]}")

# ── Step 3: is_ambiguous = FALSE ──────────────────────────────────────────────
print("\n--- STEP 3: After is_ambiguous = FALSE ---")
step3 = step2[step2.is_ambiguous == False].copy()
print(f"[COUNT] Step 3 after ambiguity filter: {len(step3):,}  lost: {len(step2)-len(step3):,}")
print(f"  cna_call distribution: {step3.cna_call.value_counts().to_dict()}")

# ── Step 4: cna_call value counts ─────────────────────────────────────────────
print("\n--- STEP 4: cna_call value_counts ---")
print(step3.cna_call.value_counts().to_string())

# ── Step 5: gene_role join ────────────────────────────────────────────────────
print("\n--- STEP 5: After gene_role join ---")
gl2 = gl.copy()
gl2["ensg_id"] = gl2["ensg_id"].str.lower()
step5 = step3.merge(gl2, on="ensg_id", how="left")
n_no_role = step5.gene_role.isna().sum()
step5["gene_role"] = step5["gene_role"].fillna("unknown")
print(f"[COUNT] Step 5 after gene_role join: {len(step5):,}  lost: {len(step3)-len(step5):,}")
print(f"  rows with no gene_role (not in lookup): {n_no_role:,}")

# ── Step 6: gene_role distribution ────────────────────────────────────────────
print("\n--- STEP 6: gene_role value_counts ---")
print(step5.gene_role.value_counts().to_string())

# ── Step 7: direction gate → has_cna_alteration ───────────────────────────────
print("\n--- STEP 7: Direction gate (oncogene→amp, tsg→del, both→either) ---")
step5["is_amp"] = step5["cna_call"] == "amplification"
step5["is_del"] = step5["cna_call"] == "deletion"
step5["has_cna_alteration"] = (
    (step5["is_amp"] & step5["gene_role"].isin(["oncogene","both"])) |
    (step5["is_del"] & step5["gene_role"].isin(["tsg","both"]))
)
n_pass = step5["has_cna_alteration"].sum()
print(f"[COUNT] Step 7 has_cna_alteration=TRUE: {n_pass:,}  lost: {len(step5)-n_pass:,}")

print(f"\n--- FULL WATERFALL SUMMARY ---")
print(f"  Step 1 (all cosmic_cna):          {len(raw):>10,}")
print(f"  Step 2 (after sample join):        {len(step2):>10,}  lost {len(raw)-len(step2):,}")
print(f"  Step 3 (is_ambiguous=FALSE):       {len(step3):>10,}  lost {len(step2)-len(step3):,}")
print(f"  Step 5 (gene_role join):           {len(step5):>10,}  [inner join, no loss expected]")
print(f"  Step 7 (direction gate passes):    {n_pass:>10,}  lost {len(step5)-n_pass:,}")
print(f"  Final saved cna_flags has_cna_alt: {pd.read_parquet(CNA_FLAGS)['has_cna_alteration'].sum():>10,}")

# ── Dominant loss ──────────────────────────────────────────────────────────────
roles_with_dir = step5[step5.gene_role.isin(["oncogene","tsg","both"])]
print(f"\n--- DOMINANT LOSS ANALYSIS ---")
print(f"  Rows for role-annotated genes (oncogene/tsg/both): {len(roles_with_dir):,}")
print(f"  Rows for gene_role=unknown (structurally unevaluable): "
      f"{(step5.gene_role=='unknown').sum():,}")
print(f"  Of role-annotated rows, wrong CNA direction: "
      f"{len(roles_with_dir)-roles_with_dir['has_cna_alteration'].sum():,}")

# ── Ploidy correlation ─────────────────────────────────────────────────────────
print("\n--- PLOIDY CORRELATION ---")
step5["total_cn"] = pd.to_numeric(step5["total_cn"], errors="coerce")
step5_num = step5.dropna(subset=["total_cn"])
print(f"  median total_cn: {step5_num.total_cn.median():.1f}")
print(f"  mean total_cn:   {step5_num.total_cn.mean():.2f}")
print(f"  lines with mean_CN > 4.5: "
      f"{(step5_num.groupby('model_id')['total_cn'].mean() > 4.5).mean()*100:.1f}%")

line_stats = step5_num.groupby("model_id").agg(
    mean_cn=("total_cn","mean"),
    alt_rate=("has_cna_alteration","mean")
).reset_index()
if len(line_stats) > 10:
    r_sp, p_sp = __import__("scipy.stats", fromlist=["spearmanr"]).spearmanr(
        line_stats["mean_cn"], line_stats["alt_rate"])
    print(f"  Spearman r(per_line_mean_CN, per_line_alt_rate): {r_sp:.4f}  p={p_sp:.3e}")
    print(f"  pre-declared: |corr| < 0.15  "
          f"{'CONFIRMED' if abs(r_sp) < 0.15 else 'WRONG -- ploidy artefact still present'}")

con.close()
print("\n" + "=" * 72)
