"""
03_load_duckdb.py

Database loading stage of the pipeline. Loads every cleaned Parquet file
under ``PARQUET_CLEAN`` — both the cleaned tables from
``01_data_cleaning.py`` and the rosters from ``02_data_harmonisation.py``
— into a DuckDB database at ``db/celllineselector.duckdb``, one table per
file, then verifies what actually landed.

Design notes
------------
* The stage is deliberately thin: discovery, loading and row counting all
  live in ``src/scripts/duckdb_utils.py``, so this file stays a readable
  two-step entry point.
* Nothing is filtered or selected here. Whatever the two preceding stages
  wrote to ``PARQUET_CLEAN`` is what gets loaded, table name taken from
  the filename stem.
* Import failures are caught up front and answered with install
  instructions rather than a bare traceback, since a missing ``duckdb``
  or ``pyarrow`` is the most likely first-run problem.

Inputs
------
Cleaned Parquet files under ``PARQUET_CLEAN``.

Outputs
-------
A DuckDB database at ``DUCKDB_PATH``, plus log entries recording how many
tables loaded successfully.

Run with
--------
    python pipelines/03_load_duckdb.py

Run this AFTER ``01_data_cleaning.py`` and ``02_data_harmonisation.py``,
since it loads whatever those two have written to ``PARQUET_CLEAN``.
"""

import sys
from pathlib import Path

# Locate project root (the repository top containing `src`) so imports work
project_root = Path(__file__).resolve().parent
for candidate in [project_root, *project_root.parents]:
    if (candidate / "src").exists():
        project_root = candidate
        break
sys.path.insert(0, str(project_root))

# Import project modules; if a dependency is missing, show actionable instructions
try:
    from src.scripts.data_utils import PARQUET_CLEAN, DUCKDB_PATH
    from src.scripts.duckdb_utils import load_all_parquets_to_duckdb, list_duckdb_tables
    from src.scripts.logging_utils import get_logger
except ModuleNotFoundError as e:
    missing = e.name if hasattr(e, "name") else str(e)
    msg = (
        f"Missing Python dependency: {missing}\n"
        "Install the runtime dependencies and try again. Recommended steps:\n"
        "  1) Create and activate a virtual environment: `python3 -m venv .venv && source .venv/bin/activate`\n"
        "  2) Install requirements: `python -m pip install --upgrade pip && pip install -r requirements.txt`\n"
        "Alternately, install the package system-wide: `pip install pandas duckdb pyarrow`\n"
    )
    print(msg)
    sys.exit(2)

logger = get_logger("03_load_duckdb")


def load_tables() -> None:
    """
    Load every cleaned Parquet file into the DuckDB database.

    Delegates to :func:`load_all_parquets_to_duckdb`, which reads each
    ``*.parquet`` under ``PARQUET_CLEAN`` and writes it to ``DUCKDB_PATH``
    as a table named after the file's stem. The per-file results table it
    returns is printed, then reduced to a pass/fail count for the log.

    Returns
    -------
    None
        The database is written and log entries emitted as side effects;
        the results table is printed to stdout.

    Notes
    -----
    Rows in the results table carry a ``Status`` column; anything other
    than ``OK`` counts as a failure. Failures are counted and warned about
    rather than raised, so one bad Parquet file does not stop the rest
    from loading. Check the log for which ones and why.
    """
    print(f"Loading parquet files from {PARQUET_CLEAN} into {DUCKDB_PATH} ...")
    results = load_all_parquets_to_duckdb()
    print("\nLoad results:")
    print(results.to_string(index=False))

    n_ok = (results["Status"] == "OK").sum()
    logger.info(f"Loaded {n_ok} / {len(results)} tables into DuckDB.")
    if n_ok < len(results):
        logger.warning(f"{len(results) - n_ok} table(s) failed to load — see log for details")


def verify_tables() -> None:
    """
    Print every table currently in the DuckDB database, with row counts.

    Reads back from ``DUCKDB_PATH`` rather than reporting what
    :func:`load_tables` believed it wrote, so the output reflects the
    database's actual state — including any tables left over from an
    earlier run.

    Returns
    -------
    None
        Output is printed to stdout.
    """
    print(f"\nVerifying tables in {DUCKDB_PATH} ...")
    tables = list_duckdb_tables()
    print(tables.to_string(index=False))


def main() -> None:
    """
    Run the full DuckDB loading stage end to end.

    Loads every cleaned Parquet file into the database, then reads the
    database back to report the tables and row counts that landed.

    Returns
    -------
    None
        The DuckDB database and log entries are written as side effects;
        both the load results and the verification listing are printed to
        stdout.
    """
    load_tables()
    verify_tables()


if __name__ == "__main__":
    main()