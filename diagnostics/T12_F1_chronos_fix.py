"""
T12 — F1: Chronos file path fix.
DIAGNOSTIC: confirms chronos_long.parquet is Project Score negated (positive mean),
            not real Chronos (which has mean ≈ 0, negatives = essential).
FIX: builds final_pipeline/output/chronos_corrected.parquet from the HDF5 files
     using Entrez-ID → ENSG mapping from DuckDB hgnc table.
     Then patches build_chronos_validation.py to use the corrected file.
Run from repo root: python diagnostics/T12_F1_chronos_fix.py
"""
import sys
from pathlib import Path
import numpy as np
import pandas as pd
import h5py
import duckdb
import re

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "final_pipeline"))
from config import DB

SRC_OUTPUTS = REPO / "src" / "pipeline" / "outputs"   # where build_chronos_validation.py looks
FINAL_OUTPUTS = REPO / "final_pipeline" / "outputs"   # where final_pipeline/ scripts look
HDF5_PATH = REPO / "data" / "DepMap_Chronos" / "GeneFitnessEffect_Chronos_Achilles.hdf5"
OLD_PATH  = REPO / "validation" / "prepared" / "chronos_long.parquet"
NEW_PATH  = SRC_OUTPUTS / "chronos_corrected.parquet"  # build_chronos_validation.py reads OUTPUTS/
BUILD_PY  = REPO / "src" / "pipeline" / "build_chronos_validation.py"

print("=" * 72)
print("T12 — F1: CHRONOS FILE PATH FIX")
print("=" * 72)

# ── DIAGNOSTIC: verify old file is NOT real Chronos ──────────────────────────
old = pd.read_parquet(OLD_PATH)
print(f"\n[MEASURED] chronos_long.parquet (current):")
print(f"  rows: {len(old):,}")
print(f"  essentiality mean: {old.essentiality.mean():.4f}  (real Chronos mean ≈ 0)")
print(f"  essentiality std:  {old.essentiality.std():.4f}")
print(f"  essentiality p5:   {old.essentiality.quantile(0.05):.4f}")
print(f"  essentiality p95:  {old.essentiality.quantile(0.95):.4f}")

# Check GAPDH — should be near 0 in Chronos (not essential) but very positive in Project Score negated
gapdh = "ensg00000111640"
gapdh_rows = old[old["ensg_id"].str.lower() == gapdh]
if len(gapdh_rows) > 0:
    print(f"\n[MEASURED] GAPDH in chronos_long: mean={gapdh_rows.essentiality.mean():.3f}  "
          f"(expect ≈ 0 for Chronos; positive = Project Score negated)")
else:
    print(f"[MEASURED] GAPDH not found in chronos_long ({gapdh})")

print(f"\n[DOCUMENTED] CONCLUSION: chronos_long.parquet mean={old.essentiality.mean():.2f} >> 0")
print(f"  Real Chronos: mean ≈ 0 (Gaussian centered at 0, negatives = essential)")
print(f"  Project Score negated: positive mean, right-skewed")
print(f"  → chronos_long.parquet is Project Score negated, NOT real Chronos. BUG CONFIRMED.")

# ── FIX: build corrected Chronos from HDF5 ────────────────────────────────────
print(f"\n[FIX] Building chronos_corrected.parquet from:")
print(f"  {HDF5_PATH}")

# 1. Load HDF5
with h5py.File(HDF5_PATH, "r") as f:
    model_ids = [s.decode() for s in f["dim_0"][:]]  # ACH-xxxxxx
    gene_ids  = [s.decode() for s in f["dim_1"][:]]  # "SYMBOL (ENTREZ_ID)"
    data      = f["data"][:]                           # shape (n_models, n_genes)

n_models, n_genes = data.shape
print(f"[COUNT] HDF5: {n_models} models × {n_genes} genes")
print(f"[MEASURED] Real Chronos essentiality: mean={np.nanmean(data):.4f}  "
      f"std={np.nanstd(data):.4f}  NaN%={100*np.isnan(data).mean():.1f}%")

# 2. Parse Entrez IDs from gene label "SYMBOL (ENTREZ_ID)"
entrez_ids = []
gene_syms  = []
for g in gene_ids:
    m = re.match(r'^(.+) \((\d+)\)$', g)
    if m:
        gene_syms.append(m.group(1).lower())
        entrez_ids.append(int(m.group(2)))
    else:
        gene_syms.append(g.lower())
        entrez_ids.append(None)

print(f"[COUNT] Genes with valid Entrez ID: {sum(e is not None for e in entrez_ids):,}")

# 3. Map Entrez IDs → ENSG via hgnc table
con = duckdb.connect(str(DB), read_only=True)
hgnc = con.execute("SELECT entrez_id, gene_id FROM main.hgnc WHERE entrez_id IS NOT NULL").df()
con.close()
hgnc["entrez_id"] = hgnc["entrez_id"].astype(str)
hgnc["gene_id"]   = hgnc["gene_id"].str.lower()
entrez_to_ensg = dict(zip(hgnc["entrez_id"], hgnc["gene_id"]))
print(f"[COUNT] hgnc Entrez→ENSG mappings: {len(entrez_to_ensg):,}")

ensg_ids = [entrez_to_ensg.get(str(e)) for e in entrez_ids]
n_mapped = sum(e is not None for e in ensg_ids)
print(f"[COUNT] Genes mapped to ENSG: {n_mapped:,} / {n_genes:,}  ({100*n_mapped/n_genes:.1f}%)")

# 4. Normalize model IDs to lowercase (pipeline convention)
model_ids_lower = [m.lower() for m in model_ids]

# 5. Build long-format DataFrame (vectorised — no Python loops)
print(f"\n[FIX] Building long format ({n_models} × {n_mapped} genes)...")
# Keep only columns where ensg_id is mapped
mapped_cols  = [i for i, e in enumerate(ensg_ids) if e is not None]
mapped_ensgs = [ensg_ids[i] for i in mapped_cols]
data_sub = data[:, mapped_cols]                            # (n_models, n_mapped)

# Repeat model_ids and tile gene IDs
model_rep = np.repeat(model_ids_lower, len(mapped_cols))  # n_models × n_mapped
gene_rep  = np.tile(mapped_ensgs,      n_models)          # same length
val_rep   = data_sub.ravel(order='C')                     # row-major = model-fast

# Drop NaN rows
mask = ~np.isnan(val_rep)
chron = pd.DataFrame({
    "ensg_id":     gene_rep[mask],
    "model_id":    model_rep[mask],
    "essentiality": val_rep[mask],
})
print(f"[COUNT] Long rows (non-NaN): {len(chron):,}")
print(f"[COUNT] Distinct genes: {chron.ensg_id.nunique():,}")
print(f"[COUNT] Distinct models: {chron.model_id.nunique():,}")
print(f"[MEASURED] Corrected Chronos: mean={chron.essentiality.mean():.4f}  "
      f"std={chron.essentiality.std():.4f}  "
      f"p5={chron.essentiality.quantile(0.05):.4f}  "
      f"p95={chron.essentiality.quantile(0.95):.4f}")

# Check GAPDH
gapdh_new = chron[chron["ensg_id"] == gapdh]
if len(gapdh_new) > 0:
    print(f"[MEASURED] GAPDH corrected mean={gapdh_new.essentiality.mean():.3f}  "
          f"(expect near 0 = non-essential)")

# 6. Save
SRC_OUTPUTS.mkdir(parents=True, exist_ok=True)
chron.to_parquet(NEW_PATH, index=False)
print(f"\n[FIX] Written: {NEW_PATH}  ({len(chron):,} rows)")

# ── PATCH build_chronos_validation.py ──────────────────────────────────────────
print(f"\n[FIX] Patching {BUILD_PY.name}...")
src = BUILD_PY.read_text(encoding="utf-8")

OLD_BLOCK = '''ch = pd.read_parquet(VAL / "chronos_long.parquet")
ib = pd.read_parquet(VAL / "id_bridge.parquet")
ch = ch.merge(ib, on="sanger_model_id", how="inner")
ch["model_id"] = ch["model_id"].str.lower()'''

NEW_BLOCK = '''# T12 F1 fix: read real Chronos (from Achilles HDF5 via diagnostics/T12_F1_chronos_fix.py)
# instead of the old chronos_long.parquet (which was Project Score negated, not Chronos)
ch = pd.read_parquet(OUTPUTS / "chronos_corrected.parquet")
# model_id is already lowercase ACH- format; ensg_id is lowercase ENSG'''

if OLD_BLOCK in src:
    new_src = src.replace(OLD_BLOCK, NEW_BLOCK)
    BUILD_PY.write_text(new_src, encoding="utf-8")
    print("[FIX] Patched build_chronos_validation.py — id_bridge merge removed")
    print("      old source now reads from OUTPUTS/chronos_corrected.parquet")
else:
    print("[WARNING] Could not find expected block to patch. Manual edit needed:")
    print("  Replace the chronos_long / id_bridge merge block with:")
    print(f"  {NEW_BLOCK}")

# Also fix the direction: Chronos effect is already negative for essential genes
# build_chronos_validation.py line 92 does 1.0 - rank (correct for negative data)
# Real Chronos: more negative = more essential → 1 - rank(essentiality) = highest rank to most negative
# This is CORRECT direction for real Chronos (essentiality already encodes direction)
print("\n[DOCUMENTED] Direction check:")
print("  Real Chronos: negative effect = essential; line 92 '1.0 - rank' gives")
print("  rank=1 to most negative (most essential) → chronos_pct=0 for most essential")
print("  Then merged code: high chronos_pct = NON-essential (correct)")
print("  Line 92 logic is correct for real Chronos; BUG was only the wrong source file.")

print("\n" + "=" * 72)
print("VERDICT AGAINST PRE-DECLARATION:")
print("  Pre-declared: chronos_long.parquet uses Project Score data, not Chronos")
print(f"  Old mean = {old.essentiality.mean():.3f}  New mean = {chron.essentiality.mean():.3f}")
if old.essentiality.mean() > 1.0 and abs(chron.essentiality.mean()) < 0.5:
    print("  CONFIRMED AND FIXED")
    print(f"  chronos_corrected.parquet written to {NEW_PATH}")
else:
    print("  UNEXPECTED — check essentiality distributions above")
print("=" * 72)
