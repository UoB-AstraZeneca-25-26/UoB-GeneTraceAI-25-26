"""
02_Proteinomics/run_protein_m2_stage.py
------------------------------------------
Thin run_all.py wrapper for Stage 2's protein scorer. Replaces
protein_scorer.py with 06_protein_score_M2.py, using the finalized,
previously-derived thresholds (results/m2_thresholds.json) so this stays
reproducible across reruns rather than re-deriving thresholds from scratch
each time. Exports the result to outputs/bulk_prot_z.parquet in the schema
core_score.py expects.

Field-mapping note (found while wiring this up, not previously checked):
protein_scorer.py's own schema docstring specifies "gene_id (ENSG uppercase),
model_id (ACH lowercase), z_t (f32), n_sources (int8)" and writes straight to
that schema. 06_protein_score_M2.py's compute_protein_z()/write_protein_z()
instead uppercases model_id (`.str.upper()`) and names its z-score column
"z_score", not "z_t". Left uncorrected, the model_id case mismatch would have
silently zeroed out core_score.py's RNA-protein inner join on
["gene_id","model_id"] (case-sensitive string comparison, no error raised) --
fixed here at the export boundary, not inside core_score.py itself.

06_protein_score_M2.py does not read protein_platform_tier.parquet or
reference "tier" anywhere (checked directly) -- it is not an input dependency
for M2's own scoring. platform_tier.py must still run before this stage
because core_score.py has a hard existence check on PROT_TIER before it will
run at all, not because M2 consumes its contents.
"""
import json
import sys
from pathlib import Path

import duckdb

PROJECT_ROOT = Path(__file__).resolve().parent.parent
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))
from config import PROT_Z

HERE = Path(__file__).resolve().parent
DB_PATH = PROJECT_ROOT / "outputs" / "celllineselector.db"
THRESHOLDS_PATH = HERE / "results" / "m2_thresholds.json"
OUTPUT_TABLE = "protein_z_m2_run_all"


def run():
    if str(HERE) not in sys.path:
        sys.path.insert(0, str(HERE))
    import importlib
    m2 = importlib.import_module("06_protein_score_M2")

    thresholds = json.loads(THRESHOLDS_PATH.read_text())

    con = duckdb.connect(str(DB_PATH), read_only=False)
    try:
        result, isoform_qc, thr = m2.compute_protein_z(con, thresholds=thresholds)
        m2.write_protein_z(con, result, OUTPUT_TABLE)
        m2._write(con, isoform_qc, "protein_isoform_qc_run_all")
    finally:
        con.close()

    _export_to_bulk_prot_z()


def _export_to_bulk_prot_z():
    con = duckdb.connect(str(DB_PATH), read_only=True)
    try:
        df = con.execute(f"""
            SELECT gene_id, LOWER(model_id) AS model_id, z_score AS z_t, n_sources
            FROM "{OUTPUT_TABLE}"
        """).df()
    finally:
        con.close()
    df.to_parquet(PROT_Z, index=False)
    print(f"\nExported {len(df):,} rows -> {PROT_Z}")
    print(f"  genes: {df.gene_id.nunique():,}  |  lines: {df.model_id.nunique():,}")


if __name__ == "__main__":
    run()
