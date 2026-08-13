"""
duckdb_catalog.py
Registers every parquet/CSV dataset in the repo as a DuckDB view, so any table
can be queried with SQL without loading it into memory first.

The catalog holds *views over the parquet files on disk* — not copies. The
.duckdb file is a few hundred KB of view definitions; the parquet files stay the
single source of truth, and a rebuilt parquet is picked up on the next query
with no re-import step.

Five schemas, one per zone of the repo:

    raw     data/parquet/raw_data        14 source extracts, uncleaned
    clean   data/parquet/data_clean      cleaned per-dataset tables
    track   cleaned_track_data           Track A-D outputs
    ref     reference                    lookup tables (join anchors)
    out     src/pipeline/outputs         Stage 0-7 scoring outputs

Usage from Python:

    from duckdb_catalog import connect, q, tables

    con = connect()                       # builds/refreshes the catalog
    q("SELECT count(*) FROM ref.gene_lookup")
    q("SELECT * FROM out.core_score WHERE ensg_id = 'ENSG00000157764' "
      "ORDER BY core_score DESC LIMIT 10")

Usage from the command line:

    python src/scripts/duckdb_catalog.py --list
    python src/scripts/duckdb_catalog.py --schema out.core_score
    python src/scripts/duckdb_catalog.py --sql "SELECT * FROM ref.gene_lookup LIMIT 5"
"""

from __future__ import annotations

import argparse
import os
import re

import duckdb
import pandas as pd

# ── Paths ────────────────────────────────────────────────────────
ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
DB_PATH = os.path.join(ROOT, "genetrace.duckdb")

# schema name -> directory holding its parquet/CSV files
ZONES = {
    "raw":   os.path.join(ROOT, "data", "parquet", "raw_data"),
    "clean": os.path.join(ROOT, "data", "parquet", "data_clean"),
    "track": os.path.join(ROOT, "cleaned_track_data"),
    "ref":   os.path.join(ROOT, "reference"),
    "out":   os.path.join(ROOT, "src", "pipeline", "outputs"),
}

# Readable view names for the numbered raw extracts and the _clean suffixes.
# Anything not listed here is auto-named from its filename (see _view_name).
NAME_OVERRIDES = {
    "raw": {
        "1_4_hpa_rna_celline":                                "hpa_rna",
        "2_DepMap_OmicsExpressionAllGenesTPMLogp1Profile":    "depmap_expr",
        "3_GEOexpression":                                    "geo_expr",
        "4_Harmonized_MS_CCLE_Gygi_subsetted":                "proteomics",
        "5_OmicsFusionFilteredSupplementary":                 "fusions",
        "6_OmicsSomaticMutationsProfile":                     "mutations",
        "7_cellosaurus":                                      "cellosaurus",
        "8_DepMap_OmicsProfiles":                             "depmap_profiles",
        "9_DepMap_sample_info":                               "sample_info",
        "10_GEOInfo":                                         "geo_info",
        "11_hpa_rna_celline_description":                     "hpa_desc",
        "12_CCLE_metabolomics_20190502":                      "metabolomics",
        "13_CCLE_miRNA_20181103":                             "mirna",
        "14_OmicsGlobalSignatures":                           "signatures",
    },
}

READABLE_EXT = {".parquet": "read_parquet", ".csv": "read_csv"}

# Subdirectories inside a zone that should also be registered, prefixed.
NESTED = {"out": ["variants"]}


def _view_name(zone: str, stem: str) -> str:
    """Map a filename stem to a legal, readable DuckDB view name."""
    override = NAME_OVERRIDES.get(zone, {}).get(stem)
    if override:
        return override
    name = re.sub(r"[^0-9a-zA-Z]+", "_", stem).strip("_").lower()
    # the clean zone suffixes every file with _clean; the schema already says so
    if zone == "clean" and name.endswith("_clean"):
        name = name[: -len("_clean")]
    if name and name[0].isdigit():
        name = "t_" + name
    return name


def _files_in(directory: str) -> list[str]:
    if not os.path.isdir(directory):
        return []
    out = []
    for fname in sorted(os.listdir(directory)):
        path = os.path.join(directory, fname)
        if os.path.isfile(path) and os.path.splitext(fname)[1].lower() in READABLE_EXT:
            out.append(path)
    return out


def build_catalog(con: duckdb.DuckDBPyConnection, verbose: bool = False) -> dict:
    """Create one schema per zone and a view per data file. Idempotent.

    Views are CREATE OR REPLACE, and any view whose backing file has since been
    deleted is dropped — otherwise it would linger as a broken view and keep the
    catalog permanently out of step with the disk. Files that DuckDB cannot read
    (e.g. a malformed CSV) are skipped with a warning rather than aborting.
    """
    registered: dict[str, list[str]] = {}
    before = _registered_views(con)

    for zone, directory in ZONES.items():
        con.execute(f'CREATE SCHEMA IF NOT EXISTS "{zone}"')
        registered[zone] = []

        targets = [("", directory)]
        for sub in NESTED.get(zone, []):
            targets.append((f"{sub}_", os.path.join(directory, sub)))

        for prefix, target_dir in targets:
            for path in _files_in(target_dir):
                stem, ext = os.path.splitext(os.path.basename(path))
                view = prefix + _view_name(zone, stem)
                reader = READABLE_EXT[ext.lower()]
                sql_path = path.replace("\\", "/").replace("'", "''")
                try:
                    con.execute(
                        f'CREATE OR REPLACE VIEW "{zone}"."{view}" AS '
                        f"SELECT * FROM {reader}('{sql_path}')"
                    )
                except duckdb.Error as exc:
                    if verbose:
                        print(f"  skip {zone}.{view}: {str(exc).splitlines()[0]}")
                    continue
                registered[zone].append(view)
                if verbose:
                    print(f"  {zone}.{view:<40} <- {os.path.relpath(path, ROOT)}")

    # drop views left behind by deleted files
    live = {(z, v) for z, vs in registered.items() for v in vs}
    for zone, view in sorted(before - live):
        con.execute(f'DROP VIEW IF EXISTS "{zone}"."{view}"')
        if verbose:
            print(f"  dropped {zone}.{view} (source file gone)")

    return registered


def _expected_views() -> set:
    """The (zone, view) pairs implied by what is currently on disk. Filesystem
    only — no parquet is opened, so this is milliseconds."""
    expected = set()
    for zone, directory in ZONES.items():
        targets = [("", directory)]
        for sub in NESTED.get(zone, []):
            targets.append((f"{sub}_", os.path.join(directory, sub)))
        for prefix, target_dir in targets:
            for path in _files_in(target_dir):
                stem = os.path.splitext(os.path.basename(path))[0]
                expected.add((zone, prefix + _view_name(zone, stem)))
    return expected


def _registered_views(con: duckdb.DuckDBPyConnection) -> set:
    rows = con.execute(
        "SELECT table_schema, table_name FROM information_schema.tables "
        "WHERE table_schema IN ('" + "', '".join(ZONES) + "')"
    ).fetchall()
    return {(z, v) for z, v in rows}


def connect(db_path: str = DB_PATH, refresh: str | bool = "auto",
            read_only: bool = False) -> duckdb.DuckDBPyConnection:
    """Open the catalog, rebuilding its views only when necessary.

    A view is `SELECT * FROM read_parquet(...)`, re-bound by DuckDB on every
    query, so new *rows* and new *columns* in an existing file show up with no
    rebuild at all. Only adding or deleting a file changes the catalog.

    refresh="auto" (default) compares the files on disk against the registered
    views and rebuilds only on a mismatch — a filesystem listing, ~10 ms.
    refresh=True forces a full rebuild (~12 s, since every parquet schema is
    re-read); refresh=False skips the check entirely.

    read_only=True opens without writing — use it when another process holds the
    database, since DuckDB allows only one writer at a time. It implies no
    refresh, so build the catalog first if the file does not exist yet.
    """
    if read_only:
        return duckdb.connect(db_path, read_only=True)
    con = duckdb.connect(db_path)
    if refresh is True:
        build_catalog(con)
    elif refresh == "auto" and _expected_views() != _registered_views(con):
        build_catalog(con)
    return con


_default_con: duckdb.DuckDBPyConnection | None = None


def _shared() -> duckdb.DuckDBPyConnection:
    global _default_con
    if _default_con is None:
        _default_con = connect()
    return _default_con


def q(sql: str, **params) -> pd.DataFrame:
    """Run SQL against the catalog and return a pandas DataFrame.

    Pass values as named parameters rather than formatting them into the string:

        q("SELECT * FROM out.core_score WHERE ensg_id = $gene", gene="ENSG00000157764")
    """
    con = _shared()
    if params:
        return con.execute(sql, params).df()
    return con.execute(sql).df()


def tables(zone: str | None = None) -> pd.DataFrame:
    """List every registered view with its column count and row count."""
    con = _shared()
    where = "WHERE table_schema = $zone" if zone else \
            "WHERE table_schema IN ('" + "', '".join(ZONES) + "')"
    listing = con.execute(
        f"SELECT table_schema AS zone, table_name AS view "
        f"FROM information_schema.tables {where} ORDER BY 1, 2",
        {"zone": zone} if zone else None,
    ).df()

    rows, cols = [], []
    for _, r in listing.iterrows():
        # index by key, not attribute: Series.view is a pandas method
        zone_name, view_name = r["zone"], r["view"]
        try:
            rows.append(con.execute(
                f'SELECT count(*) FROM "{zone_name}"."{view_name}"').fetchone()[0])
            cols.append(len(con.execute(
                f'SELECT * FROM "{zone_name}"."{view_name}" LIMIT 0').description))
        except duckdb.Error:
            rows.append(None)
            cols.append(None)
    listing["rows"] = rows
    listing["cols"] = cols
    return listing


def schema_of(qualified: str) -> pd.DataFrame:
    """Column names and types for a `zone.view` reference."""
    zone, view = qualified.split(".", 1)
    return _shared().execute(
        "SELECT column_name, data_type FROM information_schema.columns "
        "WHERE table_schema = $z AND table_name = $v ORDER BY ordinal_position",
        {"z": zone, "v": view},
    ).df()


def main() -> None:
    ap = argparse.ArgumentParser(description="GeneTraceAI DuckDB catalog")
    ap.add_argument("--build", action="store_true",
                    help="(re)build the catalog and print every view registered")
    ap.add_argument("--list", action="store_true",
                    help="list views with row and column counts")
    ap.add_argument("--zone", help="restrict --list to one zone")
    ap.add_argument("--schema", metavar="ZONE.VIEW",
                    help="print the columns and types of one view")
    ap.add_argument("--sql", help="run a query and print the result")
    args = ap.parse_args()

    if args.build or not any([args.list, args.schema, args.sql]):
        con = duckdb.connect(DB_PATH)
        print(f"catalog: {DB_PATH}")
        registered = build_catalog(con, verbose=args.build)
        total = sum(len(v) for v in registered.values())
        print(f"{total} views across {len(registered)} schemas: " +
              ", ".join(f"{z}={len(v)}" for z, v in registered.items()))
        con.close()
        if not any([args.list, args.schema, args.sql]):
            return

    pd.set_option("display.max_rows", 200, "display.width", 160)
    if args.list:
        print(tables(args.zone).to_string(index=False))
    if args.schema:
        print(schema_of(args.schema).to_string(index=False))
    if args.sql:
        print(q(args.sql).to_string(index=False))


if __name__ == "__main__":
    main()
