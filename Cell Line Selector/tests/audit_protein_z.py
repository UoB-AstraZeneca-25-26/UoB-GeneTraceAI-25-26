#!/usr/bin/env python3
"""Audit the generated protein long-format z-score table.

Usage:
    python audit_protein_z.py
    python audit_protein_z.py --table protein_z
"""

from __future__ import annotations

import argparse
from pathlib import Path

import duckdb

PROJECT_ROOT = Path(__file__).resolve().parent


def parse_args():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--db", default=str(PROJECT_ROOT / "db" / "celllineselector.duckdb"))
    ap.add_argument("--table", default="protein_z")
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
            count(DISTINCT gene_id) AS genes,
            count(DISTINCT model_id) AS cell_lines,
            count(DISTINCT lineage) AS lineages,
            round(avg(z_score), 4) AS mean_z_score,
            round(stddev_samp(z_score), 4) AS sd_z_score,
            round(avg(n_sources), 2) AS mean_n_sources,
            min(z_score) AS min_z_score,
            max(z_score) AS max_z_score
        FROM "{table}"
    """).fetchdf()
    print("\nSummary:")
    print(summary.to_string(index=False))

    sample = con.execute(f'SELECT * FROM "{table}" LIMIT {args.limit}').fetchdf()
    print(f"\nSample rows (first {len(sample)}):")
    print(sample.to_string(index=False))

    src_counts = con.execute(f"""
        SELECT n_sources, count(*) AS n_rows
        FROM "{table}"
        GROUP BY 1
        ORDER BY 1
    """).fetchdf()
    print("\nCoverage by source count:")
    print(src_counts.to_string(index=False))

    lineage_counts = con.execute(f"""
        SELECT lineage, count(*) AS n_rows
        FROM "{table}"
        GROUP BY 1
        ORDER BY 2 DESC, 1
    """).fetchdf()
    print("\nLineage counts:")
    print(lineage_counts.to_string(index=False))

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
