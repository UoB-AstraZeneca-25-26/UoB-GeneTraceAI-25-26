"""
T12 delta — Chronos fix impact (remediation item 6).
Compares chronos_validation_OLD.parquet (built from Project Score negated)
vs chronos_validation.parquet (built from real Chronos HDF5).

Produces:
  - Tier counts before/after
  - Confusion matrix of chronos_check labels
  - 3 named genes: GAPDH (housekeeping), EGFR (oncogene), TP53 (TSG)

Run from repo root: python diagnostics/T12_delta.py
"""
import sys
from pathlib import Path
import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parent.parent.parent
OLD = REPO / "architecture/outputs/chronos_validation_OLD.parquet"
NEW = REPO / "architecture/outputs/chronos_validation.parquet"
DB  = REPO / "architecture/outputs/celllineselector.db"

print("=" * 72)
print("T12 DELTA — CHRONOS FIX IMPACT")
print("=" * 72)

for p, label in [(OLD, "OLD (Project Score negated)"), (NEW, "NEW (real Chronos HDF5)")]:
    if not p.exists():
        print(f"[ERROR] {label}: file not found at {p}")
        sys.exit(1)

old = pd.read_parquet(OLD)
new = pd.read_parquet(NEW)

print(f"\n[COUNT] OLD rows: {len(old):,}  NEW rows: {len(new):,}")

# ── 1. Tier counts ────────────────────────────────────────────────────────────
print("\n--- CHRONOS_CHECK TIER COUNTS ---")
cats = ["validated","inverted","weak_positive","weak_negative","none"]

old_counts = old["chronos_check"].value_counts().reindex(cats, fill_value=0)
new_counts = new["chronos_check"].value_counts().reindex(cats, fill_value=0)

print(f"\n{'Tier':<18} {'OLD':>10} {'NEW':>10} {'Delta':>10}")
print("-" * 50)
for cat in cats:
    delta = int(new_counts[cat]) - int(old_counts[cat])
    print(f"{cat:<18} {old_counts[cat]:>10,} {new_counts[cat]:>10,} {delta:>+10,}")

# ── 2. Confusion matrix ───────────────────────────────────────────────────────
merged = old[["ensg_id","chronos_check","rho"]].rename(
    columns={"chronos_check":"check_old","rho":"rho_old"}
).merge(
    new[["ensg_id","chronos_check","rho"]].rename(
        columns={"chronos_check":"check_new","rho":"rho_new"}
    ),
    on="ensg_id", how="inner"
)
print(f"\n[COUNT] Genes in both OLD and NEW: {len(merged):,}")

print("\n--- CONFUSION MATRIX (rows=OLD, cols=NEW) ---")
conf = pd.crosstab(merged["check_old"], merged["check_new"],
                   margins=True, margins_name="TOTAL")
print(conf.to_string())

changed = merged[merged["check_old"] != merged["check_new"]]
print(f"\n[MEASURED] Genes that changed label: {len(changed):,} / {len(merged):,} "
      f"({100*len(changed)/len(merged):.1f}%)")

# Direction-flip (validated→inverted or vice versa)
flips = changed[
    (changed.check_old.isin(["validated","inverted"])) &
    (changed.check_new.isin(["validated","inverted"])) &
    (changed.check_old != changed.check_new)
]
print(f"[MEASURED] Direction-flipped (validated↔inverted): {len(flips):,}")

# ── 3. Named genes ────────────────────────────────────────────────────────────
print("\n--- NAMED GENE COMPARISON ---")
NAMED = {
    "GAPDH":  "ensg00000111640",   # housekeeping — should be ~zero rho with both
    "EGFR":   "ensg00000146648",   # oncogene
    "TP53":   "ensg00000141510",   # TSG
}
# chronos_validation uses lowercase ensg_id
print(f"\n{'Gene':<8} {'OLD_rho':>9} {'OLD_check':<16} {'NEW_rho':>9} {'NEW_check':<16} {'Delta_rho':>10}")
print("-" * 70)
for gene, ensg in NAMED.items():
    ensg_lo = ensg.lower()
    o = old[old["ensg_id"] == ensg_lo]
    n = new[new["ensg_id"] == ensg_lo]
    if len(o) == 0 or len(n) == 0:
        print(f"{gene:<8}  not found in {'OLD' if len(o)==0 else 'NEW'}")
        continue
    o_rho = o.iloc[0]["rho"]
    n_rho = n.iloc[0]["rho"]
    o_chk = o.iloc[0]["chronos_check"]
    n_chk = n.iloc[0]["chronos_check"]
    print(f"{gene:<8} {o_rho:>+9.4f} {o_chk:<16} {n_rho:>+9.4f} {n_chk:<16} {n_rho-o_rho:>+10.4f}")

# ── 4. rho distribution shift ─────────────────────────────────────────────────
print("\n--- RHO DISTRIBUTION: OLD vs NEW ---")
for label, df in [("OLD", old), ("NEW", new)]:
    r = df["rho"]
    print(f"  {label}: mean={r.mean():.4f}  median={r.median():.4f}  "
          f"SD={r.std():.4f}  p5={r.quantile(0.05):.4f}  p95={r.quantile(0.95):.4f}")

print("\n" + "=" * 72)
print("VERDICT:")
n_chng = len(changed)
n_tot  = len(merged)
print(f"  {n_chng:,}/{n_tot:,} genes ({100*n_chng/n_tot:.1f}%) changed chronos_check label.")
if len(flips) > 0:
    print(f"  {len(flips)} genes FLIPPED direction (validated↔inverted).")
print("  Distribution shift above shows systematic bias from old Chronos data.")
print("=" * 72)
