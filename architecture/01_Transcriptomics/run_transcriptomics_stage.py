"""
01_Transcriptomics/run_transcriptomics_stage.py
--------------------------------------------------
Thin run_all.py wrapper for Stage 1 (RNA). Replaces rna_scorer.py with
transcriptomics.py's full-genome scoring, driven through
05_transcriptomics_stats_layer.py with the finalized parameters
(results/transcriptomics_cached/chosen_params_v2_widened.json), then exports
the result to outputs/bulk_rna_z.parquet in the schema core_score.py expects.

05_transcriptomics_stats_layer.py has no run() -- it is an argparse CLI whose
main() writes to a DuckDB table, not a parquet file. This wrapper drives that
CLI (via a temporary sys.argv swap, so the CLI stays the single source of
truth for the scoring sequence) and then bridges its output to the
parquet schema core_score.py's RNA_Z path expects: gene_id, model_id, z_t,
k_src (core_score.py auto-detects k_src vs the old pipeline's n_sources).

Only status="ok" rows are exported. "censored_ranked" rows (~7.8% of the
full run: 1,589,909 / 20,286,012) use a different, rank-based scoring regime
for below-detection genes and are not on the same z-scale as "ok" rows --
mixing them into core_score.py's per-(gene, k_src) SD=1 standardization
without review would silently conflate two different score types with one
one another. Excluded here as a deliberate, flagged choice, not a silent
default. Revisit as its own scoped decision if RNA coverage below the
detection floor is wanted later.
"""
import sys
from pathlib import Path

import duckdb

PROJECT_ROOT = Path(__file__).resolve().parent.parent
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))
from config import RNA_Z

HERE = Path(__file__).resolve().parent
DB_PATH = PROJECT_ROOT / "outputs" / "celllineselector.db"
PARAMS_PATH = HERE / "results" / "transcriptomics_cached" / "chosen_params_v2_widened.json"
WRITE_TABLE = "transcriptomics_z_run_all"


def run():
    if str(HERE) not in sys.path:
        sys.path.insert(0, str(HERE))
    import importlib
    stats_layer = importlib.import_module("05_transcriptomics_stats_layer")

    argv_backup = sys.argv
    try:
        sys.argv = [
            "05_transcriptomics_stats_layer.py",
            "--db", str(DB_PATH),
            "--params", str(PARAMS_PATH),
            "--write-table", WRITE_TABLE,
        ]
        stats_layer.main()
    finally:
        sys.argv = argv_backup

    _export_to_bulk_rna_z()


def _export_to_bulk_rna_z():
    con = duckdb.connect(str(DB_PATH), read_only=True)
    try:
        df = con.execute(f"""
            SELECT ensg AS gene_id, model_id, z_lineage AS z_t, k_src
            FROM "{WRITE_TABLE}"
            WHERE status = 'ok'
        """).df()
    finally:
        con.close()
    df.to_parquet(RNA_Z, index=False)
    print(f"\nExported {len(df):,} rows (status='ok' only) -> {RNA_Z}")
    print(f"  genes: {df.gene_id.nunique():,}  |  lines: {df.model_id.nunique():,}")


if __name__ == "__main__":
    run()
