"""
04_eda.py

Exploratory data analysis stage of the pipeline. Runs over every cleaned
table and roster under ``PARQUET_CLEAN``, backed by DuckDB rather than
pandas: each cleaned Parquet file is registered as a DuckDB view, read
directly off disk and never fully loaded into memory, so the suite scales
to very wide tables — ``depmap_expr``, ``geo_expr`` and ``proteomics`` can
each carry thousands of columns.

Scope
-----
``depmap_expr`` and ``geo_expr`` are restricted to ``gene_roster``'s
protein-coding ``gene_id`` values — see ``GENE_MATRIX_TABLES`` in
``eda_functions.py``. Every other registered table is analysed in full.

Outputs
-------
Per-table results under ``reports/eda/<table_name>/``, plus a
cross-dataset overview and two cross-table summary plots
(``_summary_gaussian_per_column.png``,
``_summary_missingness_by_table.png``).

Caching
-------
Results are cached on disk as ``missing_value_report.csv`` and
``numeric_summary.csv`` inside each table's folder. A rerun reuses them
instead of recomputing. Pass ``--force`` to ignore the cache and recompute
everything from scratch — for example after the cleaned data has changed.

Inputs
------
Cleaned Parquet files under ``PARQUET_CLEAN``.

Run with
--------
    python pipelines/04_eda.py
    python pipelines/04_eda.py --force

Run this AFTER ``01_data_cleaning.py`` and ``02_data_harmonisation.py``,
since it analyses whatever those two have written to ``PARQUET_CLEAN``.
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
    Register every cleaned Parquet file in ``PARQUET_CLEAN`` as a DuckDB view.

    Opens an in-memory DuckDB connection and creates one view per
    ``*.parquet`` file, named after the file's stem, wrapping
    ``read_parquet`` over the path. DuckDB reads directly off disk when a
    query runs rather than materialising the file up front, which is what
    lets the EDA suite work over thousand-column tables without ever
    holding them as pandas DataFrames.

    Returns
    -------
    con : duckdb.DuckDBPyConnection
        Open in-memory connection with all views registered. The caller
        owns it and is responsible for closing it.
    table_names : list of str
        View names, in sorted filename order — the same keys used for the
        per-table output folders under ``EDA_DIR``.

    Raises
    ------
    Exception
        Any failure registering a view is logged with the source
        directory and re-raised. This is fatal by design: a partially
        registered set of views would produce a silently incomplete EDA
        report.

    Notes
    -----
    The database itself is ``:memory:`` — only view definitions live
    there, so nothing is written to disk by this function and no state
    carries between runs.
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
    """
    Run the full EDA stage end to end.

    Parses ``--force`` from the command line, registers every cleaned
    Parquet file as a DuckDB view, runs the EDA suite over all of them via
    :func:`run_eda_pipeline`, prints the cross-dataset overview, and
    closes the connection.

    Returns
    -------
    None
        Per-table reports, plots and the cross-table summaries are written
        under ``EDA_DIR`` as side effects; the overview is printed to
        stdout.

    Other Parameters
    ----------------
    --force : flag
        Read from ``sys.argv``. Ignores the on-disk cache and recomputes
        every statistic from scratch. Without it, existing
        ``missing_value_report.csv`` and ``numeric_summary.csv`` files are
        reused.

    Raises
    ------
    Exception
        Propagated from :func:`load_duckdb_views` if the cleaned Parquet
        files cannot be registered.
    """
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