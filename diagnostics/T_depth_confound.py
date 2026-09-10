"""
T_depth_confound.py — Is r(has_protein, mut_burden) = 0.308 driven by sequencing depth?

Hypothesis: well-studied lines are more deeply sequenced, so more variants are called,
inflating mut_burden independent of true biology. If sequencing depth explains the
has_protein ↔ mut_burden correlation, the 0.308 is an artefact of fame/attention bias,
not a biological signal.

Checks:
  (1) corr(mut_burden, seq_depth)   — does depth explain mutation count?
  (2) corr(has_protein, seq_depth)  — are protein-covered lines more deeply sequenced?
  (3) partial corr(has_protein, mut_burden | seq_depth) — does the 0.308 survive adjustment?

seq_depth proxy: median(dp) per cell line from the mutations table.
dp = depth at variant site (not total coverage, but a reasonable proxy when aggregated).

Source: Lawrence et al. (2013), Nature 499:214 — background mutation rate and coverage
as confounders in driver detection. PMID 23770567.
"""
import sys
from pathlib import Path
import numpy as np
import pandas as pd
import duckdb
from scipy import stats

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "final_pipeline"))
from config import DB, PROT_Z

print("=" * 70)
print("T_depth_confound — Sequencing depth as confounder in has_protein ↔ mut_burden")
print("=" * 70)

con = duckdb.connect(str(DB), read_only=True)

# Per-line mut_burden from gdsc_models
burden = con.execute(
    "SELECT lower(model_id) AS model_id, mutational_burden "
    "FROM main.gdsc_models WHERE mutational_burden IS NOT NULL"
).df()
print(f"[COUNT] Lines with mut_burden: {len(burden):,}")

# Per-line sequencing depth proxy: median dp from mutations table
depth = con.execute(
    "SELECT lower(model_id) AS model_id, "
    "       median(dp)   AS seq_depth_median, "
    "       count(*)     AS n_variants "
    "FROM main.mutations WHERE dp IS NOT NULL AND dp > 0 "
    "GROUP BY model_id "
    "HAVING count(*) >= 5"
).df()
con.close()
print(f"[COUNT] Lines with depth proxy: {len(depth):,}")

# Has-protein flag from prot_z output
prot = pd.read_parquet(PROT_Z, columns=["gene_id","model_id"])
has_prot = prot.groupby("model_id").size().reset_index(name="n_genes_prot")
has_prot["has_protein"] = 1

# Merge all three
df = burden.merge(depth, on="model_id", how="inner")
df = df.merge(has_prot[["model_id","has_protein"]], on="model_id", how="left")
df["has_protein"] = df["has_protein"].fillna(0).astype(int)
print(f"[COUNT] Lines in full merge: {len(df):,}  |  protein-covered: {df.has_protein.sum():,}")

# ── Correlations ────────────────────────────────────────────────────────────
r1, p1 = stats.spearmanr(df["mutational_burden"], df["seq_depth_median"])
r2, p2 = stats.spearmanr(df["has_protein"],       df["seq_depth_median"])
r3, p3 = stats.spearmanr(df["has_protein"],       df["mutational_burden"])

print(f"\n[MEASURED] r(mut_burden,  seq_depth)   = {r1:+.4f}  p={p1:.2e}")
print(f"[MEASURED] r(has_protein, seq_depth)   = {r2:+.4f}  p={p2:.2e}")
print(f"[MEASURED] r(has_protein, mut_burden)  = {r3:+.4f}  p={p3:.2e}  (ref: 0.308 from T11)")

# ── Partial correlation controlling for depth ────────────────────────────────
# r_partial(X,Y|Z) = (r_XY - r_XZ*r_YZ) / sqrt((1-r_XZ^2)(1-r_YZ^2))
def partial_r(rxy, rxz, ryz):
    num = rxy - rxz * ryz
    den = np.sqrt((1 - rxz**2) * (1 - ryz**2))
    return num / den if den > 0 else np.nan

r_partial = partial_r(r3, r1, r2)
print(f"\n[MEASURED] Partial r(has_protein, mut_burden | seq_depth) = {r_partial:+.4f}")
print(f"[MEASURED] Attenuation: {abs(r3):.4f} → {abs(r_partial):.4f}  "
      f"({100*(abs(r3)-abs(r_partial))/abs(r3):.1f}% reduction)")

# ── n_variants as a simpler depth proxy ─────────────────────────────────────
r4, p4 = stats.spearmanr(df["has_protein"],  df["n_variants"])
r5, p5 = stats.spearmanr(df["mutational_burden"], df["n_variants"])
r_partial2 = partial_r(r3, r5, r4)
print(f"\n[ALTERNATIVE] Using n_variants as depth proxy:")
print(f"  r(mut_burden, n_variants)   = {r5:+.4f}  p={p5:.2e}")
print(f"  r(has_protein, n_variants)  = {r4:+.4f}  p={p4:.2e}")
print(f"  Partial r(has_protein, mut_burden | n_variants) = {r_partial2:+.4f}")

# ── Verdict ─────────────────────────────────────────────────────────────────
print("\n" + "=" * 70)
attenuation_pct = 100 * (abs(r3) - abs(r_partial)) / abs(r3)
if attenuation_pct > 50:
    print(f"VERDICT: DEPTH CONFOUNDED — {attenuation_pct:.0f}% of r(has_protein, mut_burden)")
    print("  is explained by sequencing depth. The 0.308 is largely an artefact.")
elif attenuation_pct > 20:
    print(f"VERDICT: PARTIAL CONFOUND — {attenuation_pct:.0f}% of the correlation")
    print("  is explained by depth; the signal is real but inflated.")
else:
    print(f"VERDICT: NOT DEPTH-DRIVEN — only {attenuation_pct:.0f}% reduction after depth control.")
    print("  The r=0.308 reflects biology, not sequencing attention bias.")
print("=" * 70)
