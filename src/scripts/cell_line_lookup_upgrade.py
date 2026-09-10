"""
[LOOKUP] Upgrade cell_line_lookup.parquet with Track A harmonization schema
"""

import pandas as pd
from pathlib import Path

# ============================================================
# CONFIG
# ============================================================
_PROJECT_ROOT   = Path(__file__).resolve().parents[2]
REF_DIR         = _PROJECT_ROOT / "reference"
CLEANED_DIR     = _PROJECT_ROOT / "cleaned_track_data"
SAMPLE_INFO_PATH = _PROJECT_ROOT / "data/parquet/9_DepMap_sample_info.parquet"
YOUR_UNMATCHED_MUTATION_IDS_PATH = _PROJECT_ROOT / "outputs/processed/mutations_resolution_log.csv"

# ============================================================
# STEP 1 — Load both tables
# ============================================================
old_lookup = pd.read_parquet(REF_DIR / "cell_line_lookup.parquet")
trackA = pd.read_parquet(CLEANED_DIR / "trackA_harmonization.parquet")

print("=== BEFORE ===")
print(f"old cell_line_lookup: {old_lookup.shape}, columns: {old_lookup.columns.tolist()}")
print(f"trackA_harmonization: {trackA.shape}, columns: {trackA.columns.tolist()}")

# ============================================================
# STEP 2 — Carry through lineage/tissue fields from sample_info
# ============================================================
sample_info = pd.read_parquet(SAMPLE_INFO_PATH)
sample_info.columns = sample_info.columns.str.strip().str.lower()

print(f"\nsample_info columns: {sample_info.columns.tolist()}")

lineage_cols = [
    "depmap_id", "lineage", "lineage_subtype", "lineage_sub_subtype",
    "lineage_molecular_subtype", "primary_disease", "subtype",
]
missing_lineage_cols = [c for c in lineage_cols if c not in sample_info.columns]
if missing_lineage_cols:
    # Try with available columns only — some DepMap releases omit certain subtype fields
    print(f"WARNING: Missing lineage columns: {missing_lineage_cols}")
    lineage_cols = [c for c in lineage_cols if c in sample_info.columns]
    print(f"Proceeding with available lineage columns: {lineage_cols}")

lineage = sample_info[lineage_cols].copy()
lineage["depmap_id"] = lineage["depmap_id"].str.strip().str.upper()
lineage = lineage.rename(columns={"depmap_id": "model_id"})

# ============================================================
# STEP 3 — Build the new canonical table
# ============================================================

# Normalise model_id for the join
trackA = trackA.copy()
trackA["model_id"] = trackA["model_id"].str.strip().str.upper()

new_lookup = trackA.merge(lineage, on="model_id", how="left")

assert new_lookup["model_id"].duplicated().sum() == 0, \
    "model_id duplicated after merge — lineage join fanned out, investigate before proceeding"
assert len(new_lookup) == len(trackA), \
    f"Row count changed after lineage merge: {len(trackA)} -> {len(new_lookup)}"

print("\n=== AFTER ===")
print(f"new cell_line_lookup: {new_lookup.shape}, columns: {new_lookup.columns.tolist()}")

# ============================================================
# STEP 4 — Reconciliation check against Track C unmatched IDs
# ============================================================
unmapped_cellosaurus_ids = set(
    trackA.loc[trackA["cvcl_accession"].isna(), "model_id"].str.upper()
) if "cvcl_accession" in trackA.columns else set()

if YOUR_UNMATCHED_MUTATION_IDS_PATH is not None and YOUR_UNMATCHED_MUTATION_IDS_PATH.exists():
    your_unmatched = pd.read_csv(YOUR_UNMATCHED_MUTATION_IDS_PATH)
    print(f"\nmutations_resolution_log columns: {your_unmatched.columns.tolist()}")

    # Try common column names for model_id
    id_col = next(
        (c for c in your_unmatched.columns if c.lower() in ("model_id", "depmap_id", "cell_line_name")),
        your_unmatched.columns[0]
    )
    your_unmatched_ids = set(your_unmatched[id_col].dropna().str.upper())

    overlap = unmapped_cellosaurus_ids & your_unmatched_ids
    print(f"\n=== RECONCILIATION ===")
    print(f"Your unmatched mutation profile IDs: {len(your_unmatched_ids)}")
    print(f"Track A's unmapped cellosaurus IDs:  {len(unmapped_cellosaurus_ids)}")
    print(f"Overlap (explains your gap):         {len(overlap)}")
    if overlap:
        print(f"Example overlapping IDs: {list(overlap)[:10]}")
else:
    print("\n=== RECONCILIATION SKIPPED ===")
    print(f"File not found: {YOUR_UNMATCHED_MUTATION_IDS_PATH}")

# ============================================================
# STEP 5 — Duplicate CVCL check
# ============================================================
dup_cvcl = pd.Series(dtype=int)
if "cvcl_accession" in new_lookup.columns:
    dup_cvcl = (
        new_lookup.dropna(subset=["cvcl_accession"])
        .groupby("cvcl_accession")["model_id"]
        .nunique()
    )
    dup_cvcl = dup_cvcl[dup_cvcl > 1]

print(f"\n=== DUPLICATE CVCL CHECK ===")
print(f"CVCL accessions claimed by >1 model_id: {len(dup_cvcl)}")
if len(dup_cvcl) > 0:
    print(dup_cvcl)
    detail_cols = [c for c in ["model_id", "cell_line_name", "cvcl_accession", "_resolution_step"]
                   if c in new_lookup.columns]
    print("\nDetail:")
    print(new_lookup[new_lookup["cvcl_accession"].isin(dup_cvcl.index)][detail_cols]
          .sort_values("cvcl_accession"))

# ============================================================
# STEP 6 — Write new canonical file + audit log
# ============================================================
backup_path = REF_DIR / "cell_line_lookup_PRE_TRACKA_UPGRADE.parquet"
old_lookup.to_parquet(backup_path, index=False)
print(f"\nBacked up old lookup to: {backup_path}")

new_lookup.to_parquet(REF_DIR / "cell_line_lookup.parquet", index=False)
print(f"Wrote new canonical cell_line_lookup.parquet: {new_lookup.shape}")

new_cols = sorted(set(new_lookup.columns) - set(old_lookup.columns))
cvcl_col = "cvcl_accession" if "cvcl_accession" in new_lookup.columns else new_lookup.columns[0]
audit = pd.DataFrame({
    "field": ["rows_before", "rows_after", "cols_before", "cols_after",
              "new_cols_added", "resolution_rate", "duplicate_cvcl_found"],
    "value": [
        len(old_lookup), len(new_lookup),
        len(old_lookup.columns), len(new_lookup.columns),
        ", ".join(new_cols),
        f"{new_lookup[cvcl_col].notna().sum()}/{len(new_lookup)}",
        len(dup_cvcl),
    ]
})
audit_path = CLEANED_DIR / "cell_line_lookup_upgrade_audit.csv"
audit.to_csv(audit_path, index=False)
print(f"Wrote audit log: {audit_path}")
print("\nDONE. Paste this full output back before continuing.")
