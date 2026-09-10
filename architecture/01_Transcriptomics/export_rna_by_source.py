"""
01_Transcriptomics/export_rna_by_source.py
--------------------------------------------
Standalone export, NOT part of run_all.py: pulls the per-source RNA z-scores
(DepMap / HPA / GEO) that 05_transcriptomics_stats_layer.py already computed
and left sitting in celllineselector.db's transcriptomics_z_run_all table --
the combined z_lineage is all that gets exported to bulk_rna_z.parquet today;
this recovers the three inputs behind it. No recompute: same source table
run_transcriptomics_stage.py already reads, same status='ok' filter, just
three extra columns.

Writes to: final_pipeline/outputs/bulk_rna_z_by_source.parquet
Schema:    gene_id, model_id, z_depmap, z_hpa_rna, z_geo (any may be null --
           that source didn't measure this gene/line)
"""
import sys
from pathlib import Path

import duckdb

PROJECT_ROOT = Path(__file__).resolve().parent.parent
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

DB_PATH = PROJECT_ROOT / "outputs" / "celllineselector.db"
OUT_PATH = PROJECT_ROOT / "outputs" / "bulk_rna_z_by_source.parquet"
READ_TABLE = "transcriptomics_z_run_all"  # same table run_transcriptomics_stage.py exports from


def run():
    con = duckdb.connect(str(DB_PATH), read_only=True)
    try:
        df = con.execute(f"""
            SELECT ensg AS gene_id, model_id, z_depmap, z_hpa_rna, z_geo
            FROM "{READ_TABLE}"
            WHERE status = 'ok'
        """).df()
    finally:
        con.close()
    df.to_parquet(OUT_PATH, index=False)
    print(f"Exported {len(df):,} rows -> {OUT_PATH}")
    print(f"  genes: {df.gene_id.nunique():,}  |  lines: {df.model_id.nunique():,}")
    for col in ["z_depmap", "z_hpa_rna", "z_geo"]:
        print(f"  {col}: {df[col].notna().sum():,} non-null")


if __name__ == "__main__":
    run()
