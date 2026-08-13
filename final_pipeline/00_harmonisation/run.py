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
    if not DB.exists():
        print(f"[00_harmonisation] MISSING: {DB}")
        print("Run these notebooks in order from src/pipeline/:")
        print("  1. 00_harmonisation.ipynb")
        print("  2. 00b_harmonisation_hpa.ipynb")
        print("  3. 00c_geo_rebuild.ipynb")
        sys.exit(1)

    con = duckdb.connect(str(DB), read_only=True)
    tables = con.execute("SHOW TABLES").fetchdf()
    row_counts = {}
    for t in tables["name"].tolist():
        n = con.execute(f"SELECT COUNT(*) FROM {t}").fetchone()[0]
        row_counts[t] = n
    con.close()

    print("[00_harmonisation] Warehouse OK:", DB.name)
    for t, n in sorted(row_counts.items()):
        print(f"  {t:<30s} {n:>10,} rows")


if __name__ == "__main__":
    run()
