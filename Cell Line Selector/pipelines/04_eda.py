"""
04_eda.py

Exploratory data analysis over every cleaned table + roster (from
PARQUET_CLEAN), backed by DuckDB rather than pandas: each cleaned
parquet file is registered as a DuckDB view (read directly off disk,
never fully loaded into memory), so this scales to very wide tables
(depmap_expr, geo_expr, proteomics can have thousands of columns).

depmap_expr and geo_expr are restricted to gene_roster's
protein-coding gene_ids (see GENE_MATRIX_TABLES in eda_functions.py).

Output is saved per-table under reports/eda/<table_name>/, plus a
cross-dataset overview and two cross-table summary plots
(_summary_gaussian_per_column.png, _summary_missingness_by_table.png).

Results are cached on disk (missing_value_report.csv /
numeric_summary.csv under each table's folder) -- a rerun reuses
them instead of recomputing. Pass --force to ignore the cache and
recompute everything from scratch (e.g. after the cleaned data has
changed).

Run with:
    python pipelines/04_eda.py
    python pipelines/04_eda.py --force

Run this AFTER 01_data_cleaning.py and 02_build_rosters.py, since it
analyzes whatever those two have written to PARQUET_CLEAN.
"""

import sys
from pathlib import Path

project_root = Path(__file__).resolve().parent
for candidate in [project_root, *project_root.parents]:
    if (candidate / "src").exists():
        project_root = candidate
        break
sys.path.insert(0, str(project_root))

try:
    import duckdb  # type: ignore
except ImportError as e:
    print(f"Error: duckdb is not installed. Please install it with: pip install duckdb")
    raise

from src.scripts.data_utils import PARQUET_CLEAN, EDA_DIR
from src.scripts.eda_functions import run_eda_pipeline
from src.scripts.logging_utils import get_logger, log_error

logger = get_logger("04_eda")


def load_duckdb_views():
    """
    Register every cleaned parquet file in PARQUET_CLEAN as a DuckDB
    view. DuckDB reads directly off disk on query rather than
    materializing the file in memory up front, which is what lets
    the EDA suite run over thousand-column tables without loading
    them as pandas DataFrames.

    Returns (con, table_names).
    """
    print(f"Registering cleaned parquet files from {PARQUET_CLEAN} as DuckDB views ...")
    con = duckdb.connect(database=":memory:")
    table_names = []
    try:
        for path in sorted(PARQUET_CLEAN.glob("*.parquet")):
            name = path.stem
            con.execute(f'CREATE VIEW "{name}" AS SELECT * FROM read_parquet(\'{path.as_posix()}\')')
            table_names.append(name)
    except Exception as e:
        log_error(logger, step="load_duckdb_views", error=e, source_dir=str(PARQUET_CLEAN))
        raise
    print(f"Registered {len(table_names)} tables: {', '.join(table_names)}")
    return con, table_names


def main() -> None:
    force = "--force" in sys.argv
    if force:
        print("--force: ignoring cached stats, recomputing everything")

    con, table_names = load_duckdb_views()
    overview = run_eda_pipeline(con, table_names, EDA_DIR, force=force)
    print("\nCross-dataset overview:")
    print(overview.to_string(index=False))
    con.close()


if __name__ == "__main__":
    main()