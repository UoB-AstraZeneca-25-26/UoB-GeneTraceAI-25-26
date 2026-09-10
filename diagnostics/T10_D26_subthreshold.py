"""
T10 — D26: Sub-threshold expression — measured_not_expressed vs never_measured.
Verifies that RNA-only pairs scoring below threshold split into two categories:
  A) measured and expressed below threshold
  B) never measured (absent from RNA assay)
Run from repo root: python diagnostics/T10_D26_subthreshold.py
"""
import sys
from pathlib import Path
import pandas as pd
import numpy as np

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "final_pipeline"))
from config import PREDICTIONS

print("=" * 72)
print("T10 — D26: SUB-THRESHOLD EXPRESSION LANDING")
print("=" * 72)

pred = pd.read_parquet(PREDICTIONS)
print(f"[COUNT] predictions: {len(pred):,}")

# Load RNA output to check coverage
rna_files = list((REPO / "final_pipeline/output").glob("rna*.parquet"))
print(f"[COUNT] RNA output files found: {len(rna_files)}")
for f in rna_files:
    print(f"  {f.name}")

# Load the rna_scores output
rna_path = REPO / "final_pipeline/output/rna_scores.parquet"
if not rna_path.exists():
    rna_path = REPO / "final_pipeline/output"
    candidates = list(rna_path.glob("*.parquet"))
    print(f"  Available parquets: {[c.name for c in candidates]}")
    rna_path = None

if rna_path and rna_path.exists():
    rna = pd.read_parquet(rna_path)
    print(f"[COUNT] rna_scores rows: {len(rna):,}")
    print(f"[COUNT] rna_scores columns: {list(rna.columns)}")

# n_layers=1 pairs = RNA only
rna_only = pred[pred["n_layers"] == 1].copy()
print(f"\n[COUNT] n_layers=1 pairs (RNA-only): {len(rna_only):,}")

# Sub-threshold: score < 0 (not top performers)
sub = rna_only[rna_only["core_score"] < 0].copy()
above = rna_only[rna_only["core_score"] >= 0].copy()
print(f"[COUNT] RNA-only sub-zero score: {len(sub):,}  ({100*len(sub)/len(rna_only):.1f}%)")
print(f"[COUNT] RNA-only zero-or-above: {len(above):,}  ({100*len(above)/len(rna_only):.1f}%)")

# In n_layers=1 pairs, stratum_rank is the RNA robust z-score wrapper
# Check how many have extreme low scores (indicating never measured vs. measured low)
if "stratum_rank" in pred.columns:
    print(f"\n[MEASURED] RNA-only stratum_rank stats:")
    print(f"  min={rna_only.stratum_rank.min():.4f}  max={rna_only.stratum_rank.max():.4f}")
    print(f"  p1={rna_only.stratum_rank.quantile(0.01):.4f}  p99={rna_only.stratum_rank.quantile(0.99):.4f}")
    floor = (rna_only.stratum_rank == rna_only.stratum_rank.min()).sum()
    print(f"  Rows at minimum (potential never-measured floor): {floor:,}")

# Load the RNA scores output to check _measured_ flag
output_dir = REPO / "final_pipeline/output"
rna_score_candidates = ["rna_scores.parquet","rna_scored.parquet","stage1_rna.parquet",
                         "rna_output.parquet"]
rna_df = None
for name in rna_score_candidates:
    p = output_dir / name
    if p.exists():
        rna_df = pd.read_parquet(p)
        print(f"\n[COUNT] Loaded {name}: {len(rna_df):,} rows, cols: {list(rna_df.columns)[:15]}")
        break
if rna_df is None:
    all_parqs = sorted(output_dir.glob("*.parquet"))
    print(f"\n[COUNT] output dir parquets: {[p.name for p in all_parqs]}")
    # try loading stage_0 or similar
    for p in all_parqs:
        if "rna" in p.name.lower() or "stage1" in p.name.lower():
            rna_df = pd.read_parquet(p)
            print(f"  Loaded {p.name}: {len(rna_df):,} rows")
            break

if rna_df is not None:
    # Look for measured/expressed flag
    for col in ["is_expressed","measured","is_measured","silent","below_threshold","expressed"]:
        if col in rna_df.columns:
            print(f"[MEASURED] {col}: TRUE={rna_df[col].sum():,}  FALSE={(~rna_df[col]).sum():,}")

    # Check if never_measured exists (gene absent from sample)
    # Proxy: count (gene, line) pairs NOT in RNA output at all
    rna_pairs = set(zip(rna_df["ensg_id"].str.lower(), rna_df["model_id"].str.lower())) \
                if "ensg_id" in rna_df.columns and "model_id" in rna_df.columns else set()
    pred_pairs = set(zip(pred["ensg_id"], pred["model_id"]))
    n_never_measured = len(pred_pairs - rna_pairs)
    n_in_rna = len(pred_pairs & rna_pairs)
    print(f"\n[MEASURED] (gene,line) pairs in predictions: {len(pred_pairs):,}")
    print(f"[MEASURED] of which in rna_output: {n_in_rna:,}")
    print(f"[MEASURED] never in rna_output:    {n_never_measured:,}  ({100*n_never_measured/len(pred_pairs):.1f}%)")

print("\n" + "=" * 72)
print("VERDICT AGAINST PRE-DECLARATION:")
print("  Pre-declared: measured_not_expressed ≠ never_measured; both landing below threshold")
print("  Result: see counts above. Quantification depends on RNA output columns available.")
print("=" * 72)
