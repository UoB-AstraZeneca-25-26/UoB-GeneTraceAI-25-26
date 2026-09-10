"""
final_pipeline/reference/export_sample_info.py
------------------------------------------------
Exports the `sample_info` table from celllineselector.db to a small standalone
parquet, columns limited to what the live API actually reads (lineage/subtype/
disease/metadata fields for /gene/detail, the lineage filter, and the
lineage-distribution breakdown).

The live API only ever queries this one table out of celllineselector.db's
100+ tables -- everything else (20M+-row transcriptomics/proteomics tables,
pipeline test/scratch tables) is dead weight for serving. Dropping the 3.7GB
.db from the Lambda image in favor of this ~KB-sized file is the single
largest lever on cold-start latency.

Reads:  final_pipeline/outputs/celllineselector.db (sample_info table, read-only)
Writes: final_pipeline/reference/sample_info.parquet
"""
from pathlib import Path

import duckdb

DB_PATH = Path(__file__).resolve().parent.parent / "outputs" / "celllineselector.db"
OUT_PATH = Path(__file__).resolve().parent / "sample_info.parquet"

# Exactly the columns read across ranking.py's _lineage_model_ids / _lineage_meta
# / _full_metadata / the startup _lineage_map+_meta_map load.
COLUMNS = [
    "model_id",
    "lineage",
    "lineage_subtype",
    "primary_disease",
    "subtype",
    "cellosaurus_ncit_disease",
    "primary_or_metastasis",
    "sample_collection_site",
    "default_growth_pattern",
    "sex",
    "age",
]


def run():
    print("=" * 70)
    print("Reference — exporting sample_info from celllineselector.db")
    print("=" * 70)
    if not DB_PATH.exists():
        raise SystemExit(f"Missing input: {DB_PATH}")

    con = duckdb.connect(str(DB_PATH), read_only=True)
    cols = ", ".join(COLUMNS)
    df = con.execute(f"SELECT {cols} FROM main.sample_info").df()
    con.close()

    df.to_parquet(OUT_PATH, index=False)
    print(f"Wrote {len(df):,} rows x {len(COLUMNS)} cols -> {OUT_PATH}")
    print(f"Size: {OUT_PATH.stat().st_size / 1024:.1f} KB "
          f"(vs {DB_PATH.stat().st_size / 1024**3:.2f} GB for the full .db)")


if __name__ == "__main__":
    run()
