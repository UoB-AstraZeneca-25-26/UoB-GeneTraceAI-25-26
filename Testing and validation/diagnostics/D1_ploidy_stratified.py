"""Stratify ploidy correlation by cna_call to separate artefact from biology."""
import sys, duckdb
import pandas as pd
import numpy as np
from scipy import stats
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent.parent
DB = REPO / "architecture/outputs/celllineselector.db"
GENE_LKP = REPO / "reference/gene_lookup.parquet"

gl = pd.read_parquet(GENE_LKP, columns=["ensg_id","gene_role"])
gl["ensg_id"] = gl["ensg_id"].str.lower()

con = duckdb.connect(str(DB), read_only=True)
cna = con.execute("""
    SELECT lower(model_id) AS model_id, lower(gene_id) AS ensg_id,
           cna_call, total_cn
    FROM main.cosmic_cna
    WHERE is_ambiguous = FALSE AND model_id IS NOT NULL AND gene_id IS NOT NULL
""").df()
con.close()

cna = cna.merge(gl, on="ensg_id", how="left")
cna["gene_role"] = cna["gene_role"].fillna("unknown")
cna["total_cn"] = pd.to_numeric(cna["total_cn"], errors="coerce")

line_cn = cna.groupby("model_id")["total_cn"].mean().reset_index(name="mean_cn")
print(f"Lines with CN data: {len(line_cn)}")
print(f"median mean_CN: {line_cn.mean_cn.median():.2f}  mean: {line_cn.mean_cn.mean():.2f}")
print()

for call in ["amplification", "deletion", "neutral"]:
    sub = cna[cna.cna_call == call].copy()
    # For amplification: oncogene/both = "correct direction"
    # For deletion: tsg/both = "correct direction"
    if call == "amplification":
        sub["is_correct_dir"] = sub["gene_role"].isin(["oncogene","both"])
    elif call == "deletion":
        sub["is_correct_dir"] = sub["gene_role"].isin(["tsg","both"])
    else:
        sub["is_correct_dir"] = False

    per_line = sub.groupby("model_id").agg(
        n_correct=("is_correct_dir","sum"),
        n_total=("model_id","size")
    ).reset_index()
    per_line = per_line.merge(line_cn, on="model_id")
    per_line["alt_rate"] = per_line["n_correct"] / per_line["n_total"]

    r, p = stats.spearmanr(per_line["mean_cn"], per_line["alt_rate"])
    r_n, p_n = stats.spearmanr(per_line["mean_cn"], per_line["n_correct"])
    print(f"=== {call.upper()} ===")
    print(f"  n lines: {len(per_line)}")
    print(f"  r(mean_CN, correct-direction alt rate): {r:.4f}  p={p:.2e}")
    print(f"  r(mean_CN, n_correct-direction hits):   {r_n:.4f}  p={p_n:.2e}")
    # Overall rate
    r_all, p_all = stats.spearmanr(per_line["mean_cn"], per_line["n_total"])
    print(f"  r(mean_CN, total {call} calls): {r_all:.4f}  p={p_all:.2e}")
    print()

print("INTERPRETATION:")
print("  If AMP r >> 0, DEL r ~ 0 or negative -> threshold artefact (high CN lines look amplified)")
print("  If both AMP r >> 0 and DEL r >> 0 -> genome instability (biology)")
