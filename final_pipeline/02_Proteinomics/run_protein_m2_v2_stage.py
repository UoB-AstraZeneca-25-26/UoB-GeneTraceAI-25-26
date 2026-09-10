"""
02_Proteinomics/run_protein_m2_v2_stage.py  (TEST HARNESS -- not wired into
run_all.py, not the live scorer. Mirrors run_protein_m2_stage.py exactly,
except it imports 06_protein_score_M2_v2 and writes to separate DuckDB
tables and a separate output parquet, so it cannot collide with the live
protein_z_m2_run_all table or outputs/bulk_prot_z.parquet.)

Runs the Stage-2-shrinkage fix (per-(gene,lineage) MAD shrunk toward each
gene's panel-wide, own-platform MAD) end to end on the real warehouse, using
the SAME finalized thresholds (results/m2_thresholds.json) as the live
pipeline, so isoform collapse and platform-agreement confidence are held
constant -- only the z-scoring step under investigation differs.
"""
import json
import sys
from pathlib import Path

import duckdb

PROJECT_ROOT = Path(__file__).resolve().parent.parent
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

HERE = Path(__file__).resolve().parent
DB_PATH = PROJECT_ROOT / "outputs" / "celllineselector.db"
THRESHOLDS_PATH = HERE / "results" / "m2_thresholds.json"
OUTPUT_TABLE = "protein_z_m2_v2"
OUT_PARQUET = PROJECT_ROOT / "outputs" / "bulk_prot_z_v2.parquet"


def run():
    if str(HERE) not in sys.path:
        sys.path.insert(0, str(HERE))
    import importlib
    m2 = importlib.import_module("06_protein_score_M2_v2")

    thresholds = json.loads(THRESHOLDS_PATH.read_text())

    con = duckdb.connect(str(DB_PATH), read_only=False)
    try:
        result, isoform_qc, thr = m2.compute_protein_z(con, thresholds=thresholds)
        m2.write_protein_z(con, result, OUTPUT_TABLE)
        m2._write(con, isoform_qc, "protein_isoform_qc_v2")
    finally:
        con.close()

    _export_to_bulk_prot_z_v2()


def _export_to_bulk_prot_z_v2():
    con = duckdb.connect(str(DB_PATH), read_only=True)
    try:
        df = con.execute(f"""
            SELECT gene_id, LOWER(model_id) AS model_id, z_score AS z_t, n_sources
            FROM "{OUTPUT_TABLE}"
        """).df()
    finally:
        con.close()
    df.to_parquet(OUT_PARQUET, index=False)
    print(f"\nExported {len(df):,} rows -> {OUT_PARQUET}")
    print(f"  genes: {df.gene_id.nunique():,}  |  lines: {df.model_id.nunique():,}")


if __name__ == "__main__":
    run()
