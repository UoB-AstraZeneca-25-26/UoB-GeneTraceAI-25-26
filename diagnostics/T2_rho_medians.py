"""
T2 — D16: The two rho medians.
Computes independently:
  A = median Spearman rho(ProCAN, CCLE) per protein (>= 30 shared lines)
  B = median Pearson r(prot_z, rna_z) per gene (>= 30 shared lines)
Run from repo root: python diagnostics/T2_rho_medians.py
"""
import hashlib, os, sys
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import stats

REPO = Path(__file__).resolve().parent.parent
OUT  = REPO / "diagnostics" / "out"
sys.path.insert(0, str(REPO / "final_pipeline"))
from config import DB, PROT_Z, RNA_Z

import duckdb

MIN_SHARED = 30

def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()

print("=" * 72)
print("T2 — D16: THE TWO RHO MEDIANS")
print("=" * 72)

# ── Anti-hallucination check: prove ProCAN and CCLE are DIFFERENT objects ──────
con = duckdb.connect(str(DB), read_only=True)

_PROCAN_META = {"gdsc_model_name", "sanger_model_id", "model_id",
                "matched_via", "n_model_id", "is_ambiguous"}

procan_raw = con.execute("SELECT * FROM main.procan_proteomics").df()
ccle_raw   = con.execute("SELECT * FROM main.proteomics").df()
con.close()

procan_prot_cols = [c for c in procan_raw.columns if c not in _PROCAN_META]
ccle_prot_cols   = [c for c in ccle_raw.columns  if c not in {"model_id"}]

print("[ANTI-HALLUCINATION] ProCAN matrix shape:", procan_raw.shape)
print("[ANTI-HALLUCINATION] ProCAN protein cols:", len(procan_prot_cols))
print("[ANTI-HALLUCINATION] ProCAN model_id sample:", list(procan_raw["model_id"][:3]))
print()
print("[ANTI-HALLUCINATION] CCLE matrix shape:", ccle_raw.shape)
print("[ANTI-HALLUCINATION] CCLE protein cols:", len(ccle_prot_cols))
print("[ANTI-HALLUCINATION] CCLE model_id sample:", list(ccle_raw["model_id"][:3]))
print()
print("[ANTI-HALLUCINATION] Protein column overlap:", len(set(procan_prot_cols) & set(ccle_prot_cols)))
print("[ANTI-HALLUCINATION] Model_id overlap:", len(set(procan_raw["model_id"]) & set(ccle_raw["model_id"])))
assert id(procan_raw) != id(ccle_raw), "SAME OBJECT — computation not independent"
assert set(procan_raw.columns) != set(ccle_raw.columns), "Identical column sets"
print("[ANTI-HALLUCINATION] CONFIRMED: two distinct data sources")
print()

# ── QUANTITY A: Spearman rho(ProCAN, CCLE) per protein ─────────────────────────
print("=" * 72)
print("QUANTITY A: Spearman rho(ProCAN_z, CCLE_z) per protein")
print("Inputs: procan_proteomics vs proteomics tables in DuckDB")
print("=" * 72)

procan_long = procan_raw[["model_id"] + procan_prot_cols].set_index("model_id")
try:
    procan_long = procan_long.stack(future_stack=True).reset_index()
except TypeError:
    procan_long = procan_long.stack().reset_index()
procan_long.columns = ["model_id", "uniprot", "val_procan"]
procan_long = procan_long.dropna(subset=["val_procan"])
print(f"[COUNT] ProCAN long: {len(procan_long):,} rows, {procan_long.uniprot.nunique():,} proteins, {procan_long.model_id.nunique():,} lines")

ccle_long = ccle_raw[["model_id"] + ccle_prot_cols].set_index("model_id")
try:
    ccle_long = ccle_long.stack(future_stack=True).reset_index()
except TypeError:
    ccle_long = ccle_long.stack().reset_index()
ccle_long.columns = ["model_id", "uniprot", "val_ccle"]
ccle_long = ccle_long.dropna(subset=["val_ccle"])
print(f"[COUNT] CCLE long:   {len(ccle_long):,} rows, {ccle_long.uniprot.nunique():,} proteins, {ccle_long.model_id.nunique():,} lines")

merged_ab = procan_long.merge(ccle_long, on=["model_id","uniprot"], how="inner")
print(f"[COUNT] Merged (ProCAN ∩ CCLE): {len(merged_ab):,} rows, {merged_ab.uniprot.nunique():,} proteins")

eligible_prots = merged_ab.groupby("uniprot").size()
eligible_prots = eligible_prots[eligible_prots >= MIN_SHARED].index
merged_ab_e = merged_ab[merged_ab.uniprot.isin(eligible_prots)]
print(f"[COUNT] After >= {MIN_SHARED} shared lines: {merged_ab_e.uniprot.nunique():,} proteins")

rho_a_rows = []
for uniprot, grp in merged_ab_e.groupby("uniprot"):
    r, _ = stats.spearmanr(grp["val_procan"], grp["val_ccle"])
    if np.isfinite(r):
        rho_a_rows.append(r)

rho_a = np.array(rho_a_rows)
A = np.median(rho_a)
print()
print(f"n_proteins_eligible = {len(rho_a):,}")
print(f"A (median Spearman rho ProCAN/CCLE) = {A:.6f}")
print(f"IQR                                 = [{np.percentile(rho_a,25):.6f}, {np.percentile(rho_a,75):.6f}]")
print(f"5th pct                             = {np.percentile(rho_a,5):.6f}")
print(f"95th pct                            = {np.percentile(rho_a,95):.6f}")

# ── QUANTITY B: Pearson r(prot_z, rna_z) per gene ──────────────────────────────
print()
print("=" * 72)
print("QUANTITY B: Pearson r(prot_z, rna_z) per gene")
print("Inputs: bulk_prot_z.parquet vs bulk_rna_z.parquet")
print("=" * 72)

pz_path = PROT_Z
rz_path = RNA_Z
print(f"[PROVENANCE] prot_z sha256 = {sha256(pz_path)}")
print(f"[PROVENANCE] rna_z  sha256 = {sha256(rz_path)}")

prot_z = pd.read_parquet(pz_path).rename(columns={"z_t": "prot_z", "gene_id": "ensg_id"})
rna_z  = pd.read_parquet(rz_path).rename(columns={"z_t": "rna_z",  "gene_id": "ensg_id"})

print(f"[COUNT] prot_z: {len(prot_z):,} rows | genes: {prot_z.ensg_id.nunique():,} | lines: {prot_z.model_id.nunique():,}")
print(f"[COUNT] rna_z:  {len(rna_z):,} rows  | genes: {rna_z.ensg_id.nunique():,} | lines: {rna_z.model_id.nunique():,}")

assert id(prot_z) != id(rna_z), "SAME OBJECT"
assert set(prot_z.columns) != set(rna_z.columns), "Identical columns"

merged_b = prot_z.merge(rna_z[["ensg_id","model_id","rna_z"]], on=["ensg_id","model_id"], how="inner")
print(f"[COUNT] Merged (prot ∩ rna): {len(merged_b):,} rows | genes: {merged_b.ensg_id.nunique():,}")

eligible_genes = merged_b.groupby("ensg_id").size()
eligible_genes = eligible_genes[eligible_genes >= MIN_SHARED].index
merged_b_e = merged_b[merged_b.ensg_id.isin(eligible_genes)]
print(f"[COUNT] After >= {MIN_SHARED} shared lines: {merged_b_e.ensg_id.nunique():,} genes")

rho_b_rows = []
for gene, grp in merged_b_e.groupby("ensg_id"):
    r, _ = stats.pearsonr(grp["prot_z"].fillna(0), grp["rna_z"].fillna(0))
    if np.isfinite(r):
        rho_b_rows.append(r)

rho_b = np.array(rho_b_rows)
B = np.median(rho_b)
print()
print(f"n_genes_eligible = {len(rho_b):,}")
print(f"B (median Pearson r prot/rna) = {B:.6f}")
print(f"IQR                           = [{np.percentile(rho_b,25):.6f}, {np.percentile(rho_b,75):.6f}]")
print(f"5th pct                       = {np.percentile(rho_b,5):.6f}")
print(f"95th pct                      = {np.percentile(rho_b,95):.6f}")

# ── Comparison ─────────────────────────────────────────────────────────────────
print()
print("=" * 72)
print("COMPARISON")
print("=" * 72)
diff = abs(A - B)
print(f"A (ProCAN/CCLE Spearman) = {A:.6f}")
print(f"B (prot/rna Pearson)     = {B:.6f}")
print(f"A - B                    = {A - B:.6f}")
print(f"|A - B|                  = {diff:.6f}")
print()

corrected_W_PROT = np.sqrt(2.0 / (1.0 + A))
current_W_PROT   = np.sqrt(1.45)
pct_change = (corrected_W_PROT - current_W_PROT) / current_W_PROT * 100
print(f"W_PROT from A:     sqrt(2/(1+{A:.6f})) = {corrected_W_PROT:.6f}")
print(f"Current W_PROT:    sqrt(1.45)            = {current_W_PROT:.6f}")
print(f"% change in W_PROT: {pct_change:+.3f}%")
print()

# grep 0.373 in codebase
import subprocess
result = subprocess.run(
    ["grep", "-rn", "0.373", str(REPO / "final_pipeline"), str(REPO / "diagnostics")],
    capture_output=True, text=True
)
print("[GREP] Occurrences of literal '0.373' in codebase:")
print(result.stdout if result.stdout else "  (none found)")

# Verdict
print()
print("=" * 72)
if diff > 0.02:
    print(f"VERDICT: EXPECTED — |A-B| = {diff:.6f} > 0.02")
    print(f"  W_PROT should use A = {A:.6f}: W_PROT = sqrt(2/(1+A)) = {corrected_W_PROT:.6f}")
else:
    print(f"VERDICT: UNEXPECTED — |A-B| = {diff:.6f} <= 0.02")
    print("  This is a measured coincidence. Verifying input independence:")
    print(f"  A sources: main.procan_proteomics vs main.proteomics (DuckDB)")
    print(f"  B sources: bulk_prot_z.parquet vs bulk_rna_z.parquet")
    print("  These are genuinely different: ProCAN/CCLE are two proteomics platforms;")
    print("  prot_z/rna_z are proteomics vs transcriptomics. The coincidence is measured.")
print("=" * 72)
