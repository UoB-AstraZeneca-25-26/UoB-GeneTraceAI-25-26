"""
Build cell_line_lookup.parquet and union_lookup_entity.parquet.
Run from project root: python src/cell_line_lookup_build/build_lookup.py
"""

import pandas as pd
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[2]
CLEAN        = PROJECT_ROOT / "data/parquet/data_clean"
RAW          = PROJECT_ROOT / "data/parquet"
REF          = PROJECT_ROOT / "reference"
REF.mkdir(exist_ok=True)

# ============================================================
# TASK 1 — Build accession resolver
# ============================================================
print("=" * 60)
print("TASK 1 — Cellosaurus accession resolver")
print("=" * 60)

cello = pd.read_parquet(CLEAN / "cellosaurus_clean.parquet")
print(f"Loaded cellosaurus_clean: {cello.shape}")

ACC_COL = "cellosaurus_accession"
SEC_COL = "secondary accession number(s)"

# Confirm lowercase
assert cello[ACC_COL].str.match(r"^cvcl_").all(), "Primary accessions not all lowercase cvcl_"
print("Accession case: confirmed lowercase")

accession_to_primary: dict[str, str] = {}
collisions: list[tuple] = []

for _, row in cello.iterrows():
    primary = row[ACC_COL].strip()

    # map primary → itself
    if primary in accession_to_primary and accession_to_primary[primary] != primary:
        collisions.append((primary, accession_to_primary[primary], primary))
    accession_to_primary[primary] = primary

    # parse and map secondaries
    raw_sec = row[SEC_COL]
    if pd.isna(raw_sec):
        continue
    secondaries = [s.strip() for s in str(raw_sec).split(";") if s.strip()]
    for sec in secondaries:
        if sec in accession_to_primary and accession_to_primary[sec] != primary:
            collisions.append((sec, accession_to_primary[sec], primary))
        else:
            accession_to_primary[sec] = primary

assert len(collisions) == 0, (
    f"Secondary accession collision(s) detected:\n" +
    "\n".join(f"  {s!r} maps to both {a!r} and {b!r}" for s, a, b in collisions)
)

print(f"Resolver entries: {len(accession_to_primary):,}  (expect > 152,231)")
assert len(accession_to_primary) > len(cello), "Resolver should have more entries than rows"
print("TASK 1 PASSED\n")

# ============================================================
# TASK 2 — Resolve sample_info
# ============================================================
print("=" * 60)
print("TASK 2 — Resolve sample_info")
print("=" * 60)

si = pd.read_parquet(CLEAN / "sample_info_clean.parquet")
print(f"Loaded sample_info_clean: {si.shape}")

si["model_id"]  = si["depmap_id"].str.strip().str.upper()
si["rrid_norm"] = si["rrid"].str.strip().str.lower()

assert si["model_id"].str.match(r"^ACH-").all(), "Not all model_ids match ^ACH-"
print("model_id format: confirmed ^ACH-")

si["primary_cvcl"] = si["rrid_norm"].map(accession_to_primary)

class_A = si[si["rrid_norm"].notna() & si["primary_cvcl"].notna()].copy()
class_B = si[si["rrid_norm"].notna() & si["primary_cvcl"].isna()].copy()
class_C = si[si["rrid_norm"].isna()].copy()

assert len(class_A) + len(class_B) + len(class_C) == 1840, "Classes don't sum to 1840"
assert si["model_id"].is_unique, "model_id not unique across sample_info"
if len(class_C) != 22:
    print(f"WARNING: expected 22 gap-class rows, got {len(class_C)}")

print(f"Class A (resolved)           : {len(class_A)}")
print(f"Class B (rrid unresolved)    : {len(class_B)}")
print(f"Class C (gap / null rrid)    : {len(class_C)}")

# Save logs
log_B = class_B[["model_id", "rrid_norm", "cell_line_name"]].rename(columns={"rrid_norm": "rrid"})
log_B.to_csv(REF / "resolution_log_rrid_unresolved.csv", index=False)

log_C = class_C[["model_id", "cell_line_name"]].copy()
log_C["reason"] = "null_rrid"
log_C.to_csv(REF / "resolution_log_gap_class.csv", index=False)

print(f"Written: resolution_log_rrid_unresolved.csv ({len(log_B)} rows)")
print(f"Written: resolution_log_gap_class.csv ({len(log_C)} rows)")
print("TASK 2 PASSED\n")

# ============================================================
# TASK 3 — Entity grain (one row per biological entity)
# ============================================================
print("=" * 60)
print("TASK 3 — Entity grain")
print("=" * 60)

entity = (
    class_A
    .groupby("primary_cvcl")
    .agg(
        model_id_set   = ("model_id",       lambda s: set(s)),
        cell_line_names= ("cell_line_name",  lambda s: set(s.dropna())),
    )
    .reset_index()
)

entity["n_model_ids"]    = entity["model_id_set"].apply(len)
entity["is_multi_model"] = entity["n_model_ids"] > 1
entity["rrid"]           = entity["primary_cvcl"]
entity["has_depmap"]     = True
entity["gap_class"]      = False

# Attach Cellosaurus annotation columns
ANNO_WANT = {
    "cellosaurus_cell_line_name": "canonical_name",
    "diseases"                  : "diseases",
    "species of origin"         : "species_of_origin",
    "sex of cell"               : "sex",
    "age of donor at sampling"  : "age_at_sampling",
    "category"                  : "category",
    "hierarchy"                 : "parent_hierarchy",
    "comments"                  : "comments",
    "secondary accession number(s)": "secondary_accessions",
}

cello_anno = cello[[ACC_COL] + [c for c in ANNO_WANT if c in cello.columns]].copy()
cello_anno = cello_anno.rename(columns={**{ACC_COL: "primary_cvcl"}, **{k: v for k, v in ANNO_WANT.items() if k in cello.columns}})

attached = [v for k, v in ANNO_WANT.items() if k in cello.columns]
skipped  = [k for k in ANNO_WANT if k not in cello.columns]
print(f"Annotation fields attached : {attached}")
if skipped:
    print(f"Skipped (not in file)      : {skipped}")

entity = entity.merge(cello_anno, on="primary_cvcl", how="left")

assert entity["primary_cvcl"].is_unique, "primary_cvcl not unique at entity grain"
assert (entity["n_model_ids"] >= 1).all(), "Entity with zero model_ids"

assert entity.loc[entity["primary_cvcl"] == "cvcl_1122", "n_model_ids"].values[0] == 2, \
    "cvcl_1122 (chl1/chl1dm) should have n_model_ids == 2"
assert entity.loc[entity["primary_cvcl"] == "cvcl_1150", "n_model_ids"].values[0] == 2, \
    "cvcl_1150 (ctv1/ctv1dm) should have n_model_ids == 2"

entity.to_parquet(REF / "union_lookup_entity.parquet", engine="pyarrow", index=False)
print(f"Saved union_lookup_entity.parquet: {entity.shape}")
print(f"  Multi-model entities: {entity['is_multi_model'].sum()}")
print("TASK 3 PASSED\n")

# ============================================================
# TASK 4 — Scoring anchor (one row per model_id)
# ============================================================
print("=" * 60)
print("TASK 4 — Scoring anchor (cell_line_lookup)")
print("=" * 60)

# Explode entity → one row per model_id
exploded_rows = []
for _, row in entity.iterrows():
    for mid in row["model_id_set"]:
        r = row.drop(labels=["model_id_set", "cell_line_names"]).to_dict()
        r["model_id"] = mid
        exploded_rows.append(r)

anchor = pd.DataFrame(exploded_rows)

# Add cvcl_dup_group_size
dup_size = entity.set_index("primary_cvcl")["n_model_ids"].to_dict()
anchor["cvcl_dup_group_size"] = anchor["primary_cvcl"].map(dup_size)

# Carry cell_line_name from sample_info for display
si_name_map = si.set_index("model_id")["cell_line_name"].to_dict()
anchor["cell_line_name"] = anchor["model_id"].map(si_name_map)

# Add gap-class rows (class_C)
gap_rows = class_C[["model_id", "cell_line_name"]].copy()
gap_rows["primary_cvcl"]        = None
gap_rows["canonical_name"]      = None
gap_rows["rrid"]                = None
gap_rows["n_model_ids"]         = None
gap_rows["is_multi_model"]      = False
gap_rows["cvcl_dup_group_size"] = None
gap_rows["has_depmap"]          = True
gap_rows["gap_class"]           = True
for col in attached:
    gap_rows[col] = None

# Gap status classification — three-way: clean / confirmed_gone / disputed
try:
    dm_gap = pd.read_parquet(RAW / "8_DepMap_OmicsProfiles.parquet")
    dm_gap["ModelID"] = dm_gap["ModelID"].str.strip().str.upper()
    models_with_omics = set(dm_gap["ModelID"].unique())
except Exception as e:
    print(f"OmicsProfiles not available for gap classification: {e}")
    models_with_omics = set()

def classify_gap_full(row):
    name = str(row["cell_line_name"]).upper() if pd.notna(row["cell_line_name"]) else ""
    if any(tag in name for tag in ["KO", " DM", "GR", "OE", "KD", "STAG2"]):
        return "clean"
    if row["model_id"] not in models_with_omics:
        return "confirmed_gone"
    return "disputed"

gap_rows["gap_status"] = gap_rows.apply(classify_gap_full, axis=1)
anchor["gap_status"]   = None

# Union
anchor = pd.concat([anchor, gap_rows], ignore_index=True)

# Canonical column order
BASE_COLS = [
    "model_id", "primary_cvcl", "canonical_name", "rrid",
    "is_multi_model", "cvcl_dup_group_size",
    "has_depmap", "gap_class", "gap_status", "cell_line_name",
]
extra_cols = [c for c in anchor.columns if c not in BASE_COLS]
anchor = anchor[BASE_COLS + extra_cols]

# ASSERTIONS
assert anchor["model_id"].is_unique,                            "model_id not unique"
assert anchor["model_id"].str.match(r"^ACH-").all(),           "model_id not all ^ACH-"
assert len(anchor) == 1840,                                     f"Expected 1840 rows, got {len(anchor)}"
assert anchor["model_id"].str[0].str.isupper().all(),          "model_id case invariant broken"
resolved_cvcl = anchor.loc[anchor["primary_cvcl"].notna(), "primary_cvcl"]
assert resolved_cvcl.str.match(r"^cvcl_").all(),               "primary_cvcl not all ^cvcl_"
assert anchor.loc[~anchor["gap_class"], "primary_cvcl"].notna().all(), \
    "Non-gap rows have null primary_cvcl"

anchor.to_parquet(REF / "cell_line_lookup.parquet", engine="pyarrow", index=False)
print(f"Saved cell_line_lookup.parquet: {anchor.shape}")
print(f"  gap_class rows     : {anchor['gap_class'].sum()}")
print(f"  is_multi_model rows: {anchor['is_multi_model'].sum()}")
gap_counts = anchor.loc[anchor["gap_class"], "gap_status"].value_counts(dropna=False)
print(f"  gap_status counts (clean / confirmed_gone / disputed):\n{gap_counts.to_string()}")
print("TASK 4 PASSED\n")

# ============================================================
# TASK 5 — optional OmicsProfiles flags


# Optional — OmicsProfiles scoreability flags
try:
    dm = pd.read_parquet(RAW / "8_DepMap_OmicsProfiles.parquet")
    dm["ModelID"] = dm["ModelID"].str.strip().str.upper()
    dm["Datatype_low"] = dm["Datatype"].str.lower()
    flags = (
        dm.groupby("ModelID")["Datatype_low"]
        .apply(lambda s: set(s))
        .reset_index()
    )
    flags["has_rna"] = flags["Datatype_low"].apply(lambda s: "rna" in s)
    flags["has_wes"] = flags["Datatype_low"].apply(lambda s: "wes" in s)
    flags["has_wgs"] = flags["Datatype_low"].apply(lambda s: "wgs" in s)
    flags = flags.rename(columns={"ModelID": "model_id"})[["model_id","has_rna","has_wes","has_wgs"]]
    anchor = anchor.merge(flags, on="model_id", how="left")
    anchor.to_parquet(REF / "cell_line_lookup.parquet", engine="pyarrow", index=False)
    print(f"OmicsProfiles flags joined: has_rna={anchor['has_rna'].sum()}, has_wes={anchor['has_wes'].sum()}, has_wgs={anchor['has_wgs'].sum()}")
except Exception as e:
    print(f"OmicsProfiles flags skipped: {e}")

print("\nFiles written:")
for f in ["union_lookup_entity.parquet", "cell_line_lookup.parquet",
          "resolution_log_rrid_unresolved.csv", "resolution_log_gap_class.csv"]:
    p = REF / f
    print(f"  {p}  exists={p.exists()}")

