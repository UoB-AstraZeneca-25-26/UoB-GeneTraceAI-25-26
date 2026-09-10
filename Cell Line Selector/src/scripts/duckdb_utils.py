"""
duckdb_utils.py

Load cleaned parquet files (from PARQUET_CLEAN — cleaned tables +
rosters) into a DuckDB database file (db/celllineselector.duckdb) so
the pipeline output can be queried directly with SQL.

Each parquet file becomes one table, named after the file's stem
(e.g. hpa_rna.parquet -> table 'hpa_rna'). Tables are created with
CREATE OR REPLACE, so re-running the load after upstream data changes
overwrites the old table content rather than duplicating it.

Design notes
------------
* DuckDB reads the Parquet itself. Nothing passes through pandas on the
  way in — ``read_parquet`` inside a ``CREATE TABLE AS`` means a wide
  expression matrix is never materialised in Python memory.
* Failures are per-file, never fatal. A file that will not load is
  logged with its path and intended table name, reported on stdout, and
  skipped; the rest still load, and the returned summary says which
  failed.
* Every connection is closed in a ``finally`` block, so a failure
  mid-load does not leave the database file locked.
* Idempotent by construction. ``CREATE OR REPLACE`` means a re-run
  refreshes each table rather than appending, so the loader can be run
  as often as the upstream data changes.

One thing the replace semantics do NOT do: remove tables whose Parquet
file has since disappeared. A dataset dropped upstream keeps its stale
table until the database file is deleted — :func:`list_duckdb_tables`
reads back from the database, so it will show it.
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

    Parameters
    ----------
    db_path : Path or str, optional
        Database file. Defaults to ``DUCKDB_PATH``. The file itself is
        created by DuckDB on first connect; only its parent directory is
        created here.

    Returns
    -------
    duckdb.DuckDBPyConnection
        Open read-write connection. The caller owns it and must close
        it — every caller in this module does so in a ``finally`` block.

    Notes
    -----
    DuckDB permits only one write connection to a database file at a
    time, so a connection left open elsewhere will make this raise.
    """
    db_path = Path(db_path)
    db_path.parent.mkdir(parents=True, exist_ok=True)
    return duckdb.connect(str(db_path))


def load_parquet_to_table(con: "duckdb.DuckDBPyConnection", parquet_path: Path, table_name: str) -> int:
    """
    Load a single parquet file into a DuckDB table via
    CREATE OR REPLACE TABLE ... AS SELECT * FROM read_parquet(...).
    Returns the row count of the resulting table.

    DuckDB reads the file directly, so nothing passes through pandas and
    the table's types come from the Parquet schema rather than being
    re-inferred.

    Parameters
    ----------
    con : duckdb.DuckDBPyConnection
        Open, writable connection.
    parquet_path : Path or str
        Source file. Passed as a bound parameter, not interpolated.
    table_name : str
        Destination table. Replaced outright if it already exists.

    Returns
    -------
    int
        Rows in the resulting table, read back from the database rather
        than reported from the write — so it reflects what actually
        landed.

    Notes
    -----
    The table name is interpolated into the SQL, so it must come from a
    trusted source; here it is always a filename stem from the pipeline's
    own output directory.
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

    Whatever the preceding stages wrote is what gets loaded — there is no
    filtering or selection here, and the table name is always the
    filename stem.

    Parameters
    ----------
    parquet_dir : Path or str, optional
        Directory of cleaned Parquet files. Defaults to
        ``PARQUET_CLEAN``.
    db_path : Path or str, optional
        Destination database. Defaults to ``DUCKDB_PATH``.

    Returns
    -------
    pandas.DataFrame
        One row per file, with columns ``Table``, ``Rows`` and
        ``Status``. Successful loads carry ``"OK"``; failures carry
        ``"—"`` for rows and ``"ERROR — see log"`` for status. The caller
        counts the non-``OK`` rows to decide whether to warn.

    Raises
    ------
    FileNotFoundError
        If ``parquet_dir`` does not exist. This one is fatal by design:
        an absent directory means the cleaning stages have not run, so
        there is nothing to load rather than something that failed to
        load.

    Notes
    -----
    Files are processed in sorted order, so the console output and the
    returned summary are stable across runs. The connection is closed in
    a ``finally`` block, so a failure part-way through still releases the
    database file.
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

    Reads back from the database rather than reporting what a load
    believed it wrote, so the result reflects the actual state —
    including tables left over from an earlier run whose Parquet file no
    longer exists.

    Parameters
    ----------
    db_path : Path or str, optional
        Database to inspect. Defaults to ``DUCKDB_PATH``.

    Returns
    -------
    pandas.DataFrame
        Columns ``Table`` and ``Rows``, one row per table. Empty if the
        database has no tables.

    Notes
    -----
    Runs one ``COUNT(*)`` per table. DuckDB answers these from metadata,
    so it stays fast even over the wide expression tables.

    Opens a connection through :func:`get_connection`, which will create
    an empty database file if none exists — so calling this against a
    wrong path returns an empty frame rather than raising.
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