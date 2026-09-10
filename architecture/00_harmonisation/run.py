"""
00_harmonisation/run.py
------------------------
Stage 0 — checks that the DuckDB warehouse exists and summarises its contents.
The warehouse is built once by the harmonisation notebooks in src/pipeline/:
  - 00_harmonisation.ipynb   (DepMap expression, GEO, cell-line metadata)
  - 00b_harmonisation_hpa.ipynb  (HPA RNA)
  - 00c_geo_rebuild.ipynb    (GEO scale correction)

If the DB is missing, this script prints the notebooks to run and exits with
a non-zero code so run_all.py stops early.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from config import DB

import duckdb


def run():
    """
    Verify the Stage-0 warehouse exists and print a row-count summary.

    Parameters
    ----------
    None

    Returns
    -------
    None
        Prints either the missing-warehouse rebuild instructions (see
        Notes) or, if the warehouse exists, its filename and a
        per-table row count for every base table (views excluded).

    Notes
    -----
    Raises ``SystemExit(1)`` if ``DB`` does not exist, so
    ``run_all.py`` stops before attempting any later stage. The
    printed rebuild instructions still point at
    ``src/pipeline/00_harmonisation.ipynb`` and
    ``00b_enriched_harmonisation.ipynb`` — those notebooks lived under
    the old ``src/`` layout, which was removed in the September 2026
    reorganisation into ``architecture/``. If you actually hit this
    error, those exact notebooks no longer exist in this repo; you'll
    need to either recover them from an earlier commit/branch or
    rebuild the warehouse a different way.
    """
    if not DB.exists():
        print(f"[00_harmonisation] MISSING: {DB}")
        print("The pre-built DB should be at architecture/outputs/celllineselector.db")
        print("If absent, rebuild it by running these notebooks (CWD = src/pipeline/):")
        print("  1. src/pipeline/00_harmonisation.ipynb")
        print("  2. src/pipeline/00b_enriched_harmonisation.ipynb")
        print("Then move the output DB here:")
        print(f"  src/pipeline/outputs/celllineselector.db  ->  {DB}")
        sys.exit(1)

    con = duckdb.connect(str(DB), read_only=True)
    # Only count base tables — views that reference the old src/pipeline/outputs/
    # "out" schema would fail if that directory no longer exists.
    base_tables = con.execute(
        "SELECT table_name FROM information_schema.tables "
        "WHERE table_schema='main' AND table_type='BASE TABLE' "
        "ORDER BY table_name"
    ).fetchdf()
    row_counts = {}
    for t in base_tables["table_name"].tolist():
        n = con.execute(f"SELECT COUNT(*) FROM {t}").fetchone()[0]
        row_counts[t] = n
    con.close()

    print("[00_harmonisation] Warehouse OK:", DB.name)
    for t, n in sorted(row_counts.items()):
        print(f"  {t:<30s} {n:>10,} rows")


if __name__ == "__main__":
    run()
