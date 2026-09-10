"""
T11 — D27: max_vep_rank range verification.
Spec states range 1-5; measured range is 0-3. Confirm and document the
actual encoding used in the pipeline.
Run from repo root: python diagnostics/T11_D27_vep_rank.py
"""
import sys
from pathlib import Path
import pandas as pd

REPO = Path(__file__).resolve().parent.parent.parent

print("=" * 72)
print("T11 — D27: max_vep_rank RANGE VERIFICATION")
print("=" * 72)

mut = pd.read_parquet(REPO / "cleaned_track_data/mutations_collapsed.parquet",
                      columns=["ensg_id","model_id","max_vep_rank"])
print(f"[COUNT] mutations_collapsed rows: {len(mut):,}")

vc = mut["max_vep_rank"].value_counts().sort_index().to_dict()
print(f"\n[MEASURED] max_vep_rank value distribution:")
for k,v in sorted(vc.items()):
    print(f"  rank {k}: {v:,}  ({100*v/len(mut):.2f}%)")

rng_min = mut["max_vep_rank"].min()
rng_max = mut["max_vep_rank"].max()
print(f"\n[MEASURED] Range: {rng_min} to {rng_max}")
print(f"[DOCUMENTED] Spec stated range: 1 to 5")

# VEP rank encoding (from harmonisation code)
vep_encoding = {
    0: "synonymous/non-coding",
    1: "missense (low impact)",
    2: "missense (medium/high impact / splice region)",
    3: "truncating (stop_gained, frameshift, start_lost)",
}
print(f"\n[DOCUMENTED] Inferred rank encoding (from harmonisation):")
for k,v in vep_encoding.items():
    cnt = vc.get(k, 0)
    print(f"  {k} = {v}: {cnt:,} rows")

# The truncating proxy used in T5
n_trunc_proxy = (mut["max_vep_rank"] == 3).sum()
n_rank2 = (mut["max_vep_rank"] == 2).sum()
print(f"\n[MEASURED] Truncating proxy (rank=3): {n_trunc_proxy:,}  ({100*n_trunc_proxy/len(mut):.2f}%)")
print(f"[MEASURED] Mid-impact (rank=2):        {n_rank2:,}  ({100*n_rank2/len(mut):.2f}%)")

# T5 v1 used ranks 4 and 5 as truncating — those do NOT exist
for rank in [4,5]:
    n = (mut["max_vep_rank"] == rank).sum()
    print(f"[MEASURED] rank {rank}: {n:,} rows  {'(does not exist)' if n==0 else '(exists)'}")

print("\n" + "=" * 72)
print("VERDICT AGAINST PRE-DECLARATION:")
print("  Pre-declared: spec says 1-5 but measured range is 0-3; rank 3 = truncating proxy")
if rng_min == 0 and rng_max == 3:
    print(f"  CONFIRMED — range is 0-3 (not 1-5)")
    print(f"  DOCUMENTATION FIX NEEDED: spec should state range 0-3, rank 3 = truncating.")
    print(f"  T5 v1 bug (used rank 4,5 for truncating) now confirmed — those rows = 0.")
else:
    print(f"  UNEXPECTED — range is {rng_min}-{rng_max}")
print("=" * 72)
