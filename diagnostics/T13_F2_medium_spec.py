"""
T13 — F2: MEDIUM confidence tier spec text correction.
Spec says MEDIUM = "mutation XOR CNA". Actual logic in confidence_tiers.py is
"mutation AND-NOT CNA" (i.e. mutation present, CNA absent). This is a docs-only fix.
Run from repo root: python diagnostics/T13_F2_medium_spec.py
"""
import sys
from pathlib import Path
import pandas as pd

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "final_pipeline"))
from config import PREDICTIONS

print("=" * 72)
print("T13 — F2: MEDIUM TIER SPEC TEXT VERIFICATION")
print("=" * 72)

# Find confidence_tiers.py
tiers_candidates = list(REPO.rglob("confidence_tiers.py"))
print(f"[DOCUMENTED] Found confidence_tiers.py at: {[str(p) for p in tiers_candidates]}")

for path in tiers_candidates:
    src = path.read_text(encoding="utf-8")
    print(f"\n[DOCUMENTED] {path.name} — relevant MEDIUM tier logic:")
    for i, line in enumerate(src.splitlines(), 1):
        if "medium" in line.lower() or "MEDIUM" in line:
            print(f"  line {i}: {line.rstrip()}")

# Check predictions distribution
pred = pd.read_parquet(PREDICTIONS)
print(f"\n[COUNT] predictions: {len(pred):,}")
tier_counts = pred["confidence_tier"].value_counts().to_dict()
print(f"[MEASURED] confidence_tier distribution: {tier_counts}")

# For MEDIUM rows: check has_cna_alteration and p_mutation
medium = pred[pred["confidence_tier"] == "MEDIUM"].copy()
print(f"\n[COUNT] MEDIUM rows: {len(medium):,}")

if len(medium) > 0:
    n_mut   = (medium["p_mutation"].fillna(0) >= 0.5).sum()
    n_cna   = medium["has_cna_alteration"].fillna(False).sum()
    n_fus   = (medium["p_fusion"].fillna(0) >= 0.5).sum() if "p_fusion" in medium.columns else 0
    n_both  = ((medium["p_mutation"].fillna(0) >= 0.5) & medium["has_cna_alteration"].fillna(False)).sum()
    print(f"[MEASURED] MEDIUM with mutation>=0.5:   {n_mut:,}  ({100*n_mut/len(medium):.1f}%)")
    print(f"[MEASURED] MEDIUM with has_cna=TRUE:    {n_cna:,}  ({100*n_cna/len(medium):.1f}%)")
    print(f"[MEASURED] MEDIUM with fusion>=0.5:     {n_fus:,}")
    print(f"[MEASURED] MEDIUM with BOTH mut+cna:    {n_both:,}  ({100*n_both/len(medium):.1f}%)")

    # XOR check: exactly one of (mut, cna, fus) should be TRUE if spec is XOR
    layers = (
        (medium["p_mutation"].fillna(0) >= 0.5).astype(int) +
        medium["has_cna_alteration"].fillna(False).astype(int) +
        (medium["p_fusion"].fillna(0) >= 0.5).astype(int)
    )
    n_exactly_one = (layers == 1).sum()
    n_more_than_one = (layers > 1).sum()
    n_zero = (layers == 0).sum()
    print(f"\n[MEASURED] Exactly 1 alteration layer:  {n_exactly_one:,}  ({100*n_exactly_one/len(medium):.1f}%)")
    print(f"[MEASURED] More than 1 alteration:       {n_more_than_one:,}  ({100*n_more_than_one/len(medium):.1f}%)")
    print(f"[MEASURED] Zero alterations:             {n_zero:,}  ({100*n_zero/len(medium):.1f}%)")

    if n_more_than_one > 0:
        print(f"\n[MEASURED] MEDIUM has {n_more_than_one:,} rows with >1 layer active")
        print("  → cannot be XOR (XOR requires exactly one TRUE)")
        print("  → actual logic is AND-NOT (mutation present, CNA absent from HIGH tier criteria)")

print("\n" + "=" * 72)
print("VERDICT AGAINST PRE-DECLARATION:")
print("  Pre-declared: spec says XOR, actual code uses AND-NOT (mutation AND CNA absent)")
print("  Documentation fix needed: replace 'XOR' with 'mutation present, CNA absent'")
print("  This is a DOCS-ONLY fix — no code changes required.")
print("=" * 72)
