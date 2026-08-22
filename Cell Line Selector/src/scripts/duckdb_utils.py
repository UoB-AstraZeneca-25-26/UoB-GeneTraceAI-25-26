"""
duckdb_utils.py

Load cleaned parquet files (from PARQUET_CLEAN — cleaned tables +
rosters) into a DuckDB database file (db/celllineselector.duckdb) so
the pipeline output can be queried directly with SQL.

Each parquet file becomes one table, named after the file's stem
(e.g. hpa_rna.parquet -> table 'hpa_rna'). Tables are created with
CREATE OR REPLACE, so re-running the load after upstream data changes
overwrites the old table content rather than duplicating it.
"""

from pathlib import Path
import importlib
from typing import TYPE_CHECKING

import pandas as pd

# Import duckdb for runtime; keep static type imports only for type checkers
if TYPE_CHECKING:
    import duckdb  # type: ignore
else:
    try:
        duckdb = importlib.import_module("duckdb")
    except Exception:
        duckdb = None

from src.scripts.data_utils import PARQUET_CLEAN, DUCKDB_PATH
from src.scripts.logging_utils import get_logger, log_file_error

logger = get_logger("duckdb_utils")


def get_connection(db_path: Path = DUCKDB_PATH) -> "duckdb.DuckDBPyConnection":
    """
    Open (creating if needed) a connection to the DuckDB database file
    at `db_path`. Creates the parent directory if it doesn't exist.
    """
    db_path = Path(db_path)
    db_path.parent.mkdir(parents=True, exist_ok=True)
    return duckdb.connect(str(db_path))


def load_parquet_to_table(con: "duckdb.DuckDBPyConnection", parquet_path: Path, table_name: str) -> int:
    """
    Load a single parquet file into a DuckDB table via
    CREATE OR REPLACE TABLE ... AS SELECT * FROM read_parquet(...).
    Returns the row count of the resulting table.
    """
    parquet_path = Path(parquet_path)
    con.execute(
        f'CREATE OR REPLACE TABLE "{table_name}" AS SELECT * FROM read_parquet(?)',
        [str(parquet_path)],
    )
    return con.execute(f'SELECT COUNT(*) FROM "{table_name}"').fetchone()[0]


def load_all_parquets_to_duckdb(parquet_dir: Path = PARQUET_CLEAN, db_path: Path = DUCKDB_PATH) -> pd.DataFrame:
    """
    Load every .parquet file in `parquet_dir` into its own table in
    the DuckDB database at `db_path`. A failure on one file is logged
    and skipped — loading continues for the rest.

    Returns a summary DataFrame of what was loaded (table, rows, status).
    """
    parquet_dir = Path(parquet_dir)
    if not parquet_dir.exists():
        raise FileNotFoundError(f"Cleaned parquet directory not found: {parquet_dir}")

    con = get_connection(db_path)
    results = []

    try:
        for path in sorted(parquet_dir.glob("*.parquet")):
            table_name = path.stem
            try:
                row_count = load_parquet_to_table(con, path, table_name)
                print(f"  loaded {table_name} ({row_count:,} rows)")
                results.append({"Table": table_name, "Rows": f"{row_count:,}", "Status": "OK"})
            except Exception as e:
                log_file_error(
                    logger, file_path=path, step="load_all_parquets_to_duckdb",
                    error=e, table=table_name,
                )
                print(f"  [FAILED] {table_name}: load into DuckDB failed — see log")
                results.append({"Table": table_name, "Rows": "—", "Status": "ERROR — see log"})
    finally:
        con.close()

    return pd.DataFrame(results)


def list_duckdb_tables(db_path: Path = DUCKDB_PATH) -> pd.DataFrame:
    """
    Return a DataFrame of every table currently in the DuckDB database
    at `db_path`, with its row count.
    """
    con = get_connection(db_path)
    try:
        tables = con.execute("SHOW TABLES").fetchdf()
        rows = []
        for table_name in tables["name"]:
            row_count = con.execute(f'SELECT COUNT(*) FROM "{table_name}"').fetchone()[0]
            rows.append({"Table": table_name, "Rows": f"{row_count:,}"})
    finally:
        con.close()

    return pd.DataFrame(rows)