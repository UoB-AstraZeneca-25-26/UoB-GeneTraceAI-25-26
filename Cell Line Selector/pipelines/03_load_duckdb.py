"""
03_load_duckdb.py

Loads every cleaned parquet file (cleaned tables + rosters, from
PARQUET_CLEAN) into a DuckDB database at db/celllineselector.duckdb,
one table per file, then verifies what landed.

Run with:
    python pipelines/03_load_duckdb.py

Run this AFTER 01_data_cleaning.py and 02_build_rosters.py, since it
loads whatever those two have written to PARQUET_CLEAN.
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
    """Load every cleaned parquet file into the DuckDB database."""
    print(f"Loading parquet files from {PARQUET_CLEAN} into {DUCKDB_PATH} ...")
    results = load_all_parquets_to_duckdb()
    print("\nLoad results:")
    print(results.to_string(index=False))

    n_ok = (results["Status"] == "OK").sum()
    logger.info(f"Loaded {n_ok} / {len(results)} tables into DuckDB.")
    if n_ok < len(results):
        logger.warning(f"{len(results) - n_ok} table(s) failed to load — see log for details")


def verify_tables() -> None:
    """Print every table currently in the DuckDB database, with row counts."""
    print(f"\nVerifying tables in {DUCKDB_PATH} ...")
    tables = list_duckdb_tables()
    print(tables.to_string(index=False))


def main() -> None:
    load_tables()
    verify_tables()


if __name__ == "__main__":
    main()