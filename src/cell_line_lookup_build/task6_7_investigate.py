"""
Task 6 — Identify and log the 4 unknown multi-model entities
Task 7 — Investigate the 17 disputed gap rows
Run from project root: python src/cell_line_lookup_build/task6_7_investigate.py
"""

import pandas as pd
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[2]
REF          = PROJECT_ROOT / "reference"

entity = pd.read_parquet(REF / "union_lookup_entity.parquet")
anchor = pd.read_parquet(REF / "cell_line_lookup.parquet")
gap_log = pd.read_csv(REF / "resolution_log_gap_class.csv")

# ============================================================
# TASK 6 — Multi-model entities
# ============================================================
print("=" * 65)
print("TASK 6 — Multi-model entities")
print("=" * 65)

multi = entity[entity["is_multi_model"]].copy()
print(f"Total multi-model entities: {len(multi)}")
print()

# Available columns to display
anno_cols = ["primary_cvcl", "canonical_name", "n_model_ids",
             "diseases", "species_of_origin", "parent_hierarchy"]
display_cols = [c for c in anno_cols if c in multi.columns]
print("All 6 multi-model rows:")
pd.set_option("display.max_colwidth", 80)
pd.set_option("display.width", 200)
print(multi[display_cols + ["model_id_set"]].to_string(index=False))

# ── Join cell_line_names per model_id from scoring anchor ───────────
anchor_map = anchor.set_index("model_id")[["cell_line_name", "gap_class",
                                           "has_rna", "has_wes", "has_wgs"]].to_dict("index")

rows = []
KNOWN = {"cvcl_1122", "cvcl_1150"}

for _, row in multi.iterrows():
    pcvcl = row["primary_cvcl"]
    model_ids = sorted(row["model_id_set"])
    cell_names = []
    all_non_gap = []
    for mid in model_ids:
        info = anchor_map.get(mid, {})
        cell_names.append(info.get("cell_line_name") or "?")
        all_non_gap.append(not info.get("gap_class", True))

    rows.append({
        "primary_cvcl"   : pcvcl,
        "canonical_name" : row.get("canonical_name", ""),
        "n_model_ids"    : row["n_model_ids"],
        "model_ids"      : " | ".join(model_ids),
        "cell_line_names": " | ".join(cell_names),
        "tissue"         : row.get("species_of_origin", ""),
        "diseases"       : str(row.get("diseases", ""))[:80],
        "all_gap_false"  : all(all_non_gap),
        "known"          : "known (spec)" if pcvcl in KNOWN else "discovered",
    })

summary = pd.DataFrame(rows)

print()
print("=" * 65)
print("CLEAN TABLE — all multi-model entities")
print("=" * 65)
print(summary[["primary_cvcl","canonical_name","model_ids","cell_line_names",
               "diseases","all_gap_false","known"]].to_string(index=False))

summary.to_csv(REF / "multi_model_entities.csv", index=False)
print(f"\nSaved: reference/multi_model_entities.csv")

# ============================================================
# TASK 7 — Disputed gap rows
# ============================================================
print()
print("=" * 65)
print("TASK 7 — Disputed gap rows")
print("=" * 65)

# gap_status lives on the scoring anchor — join it onto the gap log
gap_anchor_cols = ["model_id", "gap_status", "has_rna", "has_wes", "has_wgs"]
gap_log = gap_log.merge(anchor[gap_anchor_cols], on="model_id", how="left")

disputed = gap_log[gap_log["gap_status"] == "disputed"].copy()
print(f"Disputed gap rows: {len(disputed)}  (expected 17)")

# Fill missing omics flags as False (not profiled at all)
for col in ["has_rna", "has_wes", "has_wgs"]:
    disputed[col] = disputed[col].fillna(False)

print("All 17 disputed rows:")
print(disputed[["model_id", "cell_line_name", "gap_status",
                "has_rna", "has_wes", "has_wgs"]].to_string(index=False))

# Classify
disputed["has_any_omics"]       = disputed[["has_rna","has_wes","has_wgs"]].any(axis=1)
disputed["exclude_candidate"]   = ~disputed["has_any_omics"]

with_omics   = disputed[disputed["has_any_omics"]]
without_omics= disputed[~disputed["has_any_omics"]]

print()
print(f"Disputed WITH omics  (scoreable, identity-uncertain) : {len(with_omics)}")
print(f"Disputed WITH NO omics (dead weight, exclude candidate): {len(without_omics)}")

print()
print("=" * 65)
print("NO-OMICS DISPUTED ROWS — raise with Daniel / cross-check Track A")
print("=" * 65)
if len(without_omics):
    print(without_omics[["model_id", "cell_line_name",
                          "has_rna", "has_wes", "has_wgs",
                          "exclude_candidate"]].to_string(index=False))
else:
    print("None — all disputed rows have at least one omics modality.")

# Save reviewed gap log — omics flags already on gap_log from the earlier join
gap_log_full = gap_log.copy()
for col in ["has_rna", "has_wes", "has_wgs"]:
    gap_log_full[col] = gap_log_full[col].fillna(False)
gap_log_full["has_any_omics"]     = gap_log_full[["has_rna","has_wes","has_wgs"]].any(axis=1)
gap_log_full["exclude_candidate"] = (gap_log_full["gap_status"] == "disputed") & ~gap_log_full["has_any_omics"]

gap_log_full.to_csv(REF / "resolution_log_gap_class_reviewed.csv", index=False)
print(f"\nSaved: reference/resolution_log_gap_class_reviewed.csv ({len(gap_log_full)} rows)")

print()
print("ALL TASKS COMPLETE")
