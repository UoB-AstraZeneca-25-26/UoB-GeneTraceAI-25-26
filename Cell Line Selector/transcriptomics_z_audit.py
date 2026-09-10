#!/usr/bin/env python3
"""Audit the generated transcriptomics long-format score table.

Usage:
    python transcriptomics_z_audit.py
    python transcriptomics_z_audit.py --table transcriptomics_z
"""

from __future__ import annotations

import argparse
from pathlib import Path

import duckdb
import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parent


def parse_args():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--db", default=str(PROJECT_ROOT / "db" / "celllineselector.duckdb"))
    ap.add_argument("--table", default="transcriptomics_z")
    ap.add_argument("--limit", type=int, default=10)
    return ap.parse_args()


def main():
    args = parse_args()
    con = duckdb.connect(args.db, read_only=True)

    table = args.table
    try:
        cols = con.execute(f'PRAGMA table_info("{table}")').fetchdf()
    except Exception as e:
        raise RuntimeError(f"Table '{table}' not found in database: {args.db}\n{e}") from e

    print(f"Database: {args.db}")
    print(f"Table: {table}")
    print(f"Columns: {cols['name'].tolist()}")

    summary = con.execute(f"""
        SELECT
            count(*) AS rows,
            count(DISTINCT ensg) AS genes,
            count(DISTINCT model_id) AS cell_lines,
            round(avg(k_src), 2) AS mean_k_src,
            round(avg(z_lineage), 4) AS mean_z_lineage,
            round(stddev_samp(z_lineage), 4) AS sd_z_lineage
        FROM "{table}"
    """).fetchdf()
    print("\nSummary:")
    print(summary.to_string(index=False))

    sample = con.execute(f"SELECT * FROM \"{table}\" LIMIT {args.limit}").fetchdf()
    print(f"\nSample rows (first {len(sample)}):")
    print(sample.to_string(index=False))

    # Also show k_src distribution and status counts
    src_counts = con.execute(f"""
        SELECT k_src, count(*) AS n
        FROM "{table}"
        GROUP BY 1
        ORDER BY 1
    """).fetchdf()
    print("\nSource coverage by k_src:")
    print(src_counts.to_string(index=False))

    status_counts = con.execute(f"""
        SELECT status, count(*) AS n
        FROM "{table}"
        GROUP BY 1
        ORDER BY 1
    """).fetchdf()
    print("\nStatus counts:")
    print(status_counts.to_string(index=False))


if __name__ == "__main__":
    raise SystemExit(main())
