"""
build_warehouse_views.py
------------------------
Exposes every pipeline output inside the harmonisation warehouse, so
`outputs/celllineselector.db` is the single SQL entry point for the whole
project rather than one of two half-databases.

    python src/pipeline/build_warehouse_views.py

WHY VIEWS AND NOT TABLES
    `00_harmonisation` and `00b_enriched_harmonisation` write *materialised*
    tables into the warehouse, because harmonisation needs SQL while it runs.
    The Stage 2-7 outputs are different: they are already parquet on disk, and
    the largest of them (`core_score.parquet`, ~28M rows) would add over a
    gigabyte to the database if copied in.

    A view is `SELECT * FROM read_parquet(...)`, re-bound by DuckDB on every
    query. So a rebuilt parquet is picked up on the next query with no re-import
    step, and the parquet files stay the single source of truth. This is the
    same design `src/scripts/duckdb_catalog.py` uses for `genetrace.duckdb`;
    the difference is that these views live in the warehouse, alongside the
    harmonisation tables, so a query can join a score to the identity hub
    without attaching a second database.

SCHEMAS
    main    materialised    the 28 harmonisation tables (written by 00 / 00b)
    out     views           src/pipeline/outputs      Stage 0-7 scoring outputs
    ref     views           reference                 lookup tables (join anchors)
    track   views           cleaned_track_data        Track A-D outputs
    clean   views           data/parquet/data_clean   cleaned per-dataset tables

Idempotent: CREATE OR REPLACE throughout, and views whose backing file has been
deleted are dropped rather than left behind as broken references.
"""

from __future__ import annotations

import argparse
import os
import re

import duckdb
import pandas as pd

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
DB_PATH = os.path.join(ROOT, "src", "pipeline", "outputs", "celllineselector.db")

# schema name -> directory of parquet/CSV files to expose as views
ZONES = {
    "out":   os.path.join(ROOT, "src", "pipeline", "outputs"),
    "ref":   os.path.join(ROOT, "reference"),
    "track": os.path.join(ROOT, "cleaned_track_data"),
    "clean": os.path.join(ROOT, "data", "parquet", "data_clean"),
}

# Subdirectories inside a zone that are also registered, with a name prefix.
NESTED = {"out": ["variants"]}

READABLE_EXT = {".parquet": "read_parquet", ".csv": "read_csv"}


def _view_name(zone: str, stem: str) -> str:
    """Map a filename stem to a legal, readable DuckDB view name."""
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
    return [
        os.path.join(directory, f)
        for f in sorted(os.listdir(directory))
        if os.path.isfile(os.path.join(directory, f))
        and os.path.splitext(f)[1].lower() in READABLE_EXT
    ]


def _registered(con: duckdb.DuckDBPyConnection) -> set:
    rows = con.execute(
        "SELECT table_schema, table_name FROM information_schema.tables "
        "WHERE table_schema IN ('" + "', '".join(ZONES) + "')"
    ).fetchall()
    return {(z, v) for z, v in rows}


def build(con: duckdb.DuckDBPyConnection, verbose: bool = False) -> dict:
    registered: dict[str, list[str]] = {}
    before = _registered(con)

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
                    # a malformed file should not abort the whole catalog
                    print(f"  skip {zone}.{view}: {str(exc).splitlines()[0]}")
                    continue
                registered[zone].append(view)
                if verbose:
                    print(f"  {zone}.{view:<42} <- {os.path.relpath(path, ROOT)}")

    live = {(z, v) for z, vs in registered.items() for v in vs}
    for zone, view in sorted(before - live):
        con.execute(f'DROP VIEW IF EXISTS "{zone}"."{view}"')
        print(f"  dropped {zone}.{view} (source file gone)")

    return registered


# ── Bridge views ─────────────────────────────────────────────────
# THE CASING TRAP: the Stage 2-7 outputs key on UPPERCASE ACH- / ENSG (rule 1 and
# rule 2 in 02_core_score.ipynb), while the harmonisation warehouse stores both
# lowercase. A cross-schema join written the obvious way returns ZERO rows and
# looks like missing data rather than a key mismatch.
#
# These views do the normalisation once, so nobody has to remember. They are
# views, so they cost nothing until queried.
BRIDGE_VIEWS = {
    # ---- v_cell_line: every fact the warehouse holds about one cell line -----
    # Identity, clinical annotation, and modality coverage in one row. Assembled
    # once here so the CLI, the web export and any ad-hoc query all describe a
    # cell line the same way instead of each re-joining four tables.
    "v_cell_line": ("""
        SELECT c.model_id,
               h.canonical_name,
               s.cell_line_name,
               s.stripped_cell_line_name        AS stripped_name,
               s.ccle_name,
               -- clinical / biological annotation (DepMap)
               s.lineage, s.lineage_subtype, s.primary_disease, s.subtype,
               s.sex, s.age, s.primary_or_metastasis, s.sample_collection_site,
               s.source,
               -- Sanger/GDSC annotation, independent of DepMap's
               g.tissue, g.cancer_type, g.cancer_type_detail, g.tissue_status,
               g.model_type, g.is_organoid, g.growth_properties,
               g.msi_status, g.ploidy_wes, g.mutational_burden,
               g.gender AS gdsc_gender, g.ethnicity, g.age_at_sampling,
               -- identity axes
               c.rrids, c.sanger_ids, c.cosmic_ids, c.ccle_ids, c.profile_ids,
               c.geo_accessions, c.in_roster, c.source_wave, c.sanger_id_conflict,
               h.gap_class,
               -- coverage: which of the 11 measured layers reach this line
               m.n_modalities, m.complete_cycle,
               m.depmap_expr, m.geo_expr, m.hpa_rna, m.mutations, m.fusions,
               m.proteomics, m.procan_proteomics, m.metabolomics, m.mirna,
               m.signatures, m.cosmic_cna
        FROM cell_line_connection_enriched c
        LEFT JOIN sample_info s              USING (model_id)
        LEFT JOIN gdsc_models g              USING (model_id)
        LEFT JOIN coverage_matrix_enriched m USING (model_id)
        LEFT JOIN out.harmonised_enriched h  USING (model_id)
    """, ["cell_line_connection_enriched", "sample_info", "gdsc_models",
          "coverage_matrix_enriched", "out.harmonised_enriched"]),

    # ---- v_gene: every fact about one gene ----------------------------------
    # gene_dispersion is joined in because signal_spread is what tells a reader
    # whether a ranking for this gene means anything at all (see
    # evidence_state.py) -- it belongs with the gene's identity, not buried in a
    # separate table.
    "v_gene": ("""
        SELECT g.gene_id,
               g.hugo_symbol, g.gene_full_name, g.gene_names, g.hugo_symbol_all,
               g.hgnc_id, g.locus_group, g.locus_type, g.hgnc_status,
               g.is_protein_coding, g.gene_role, g.cgc_tier,
               g.uniprot_ids, g.n_uniprot_ids, g.name_conflict,
               d.n_lines            AS n_lines_expressed,
               d.median_log2tpm, d.dynamic_range, d.top_vs_median_fold,
               d.signal_spread
        FROM gene_enriched g
        LEFT JOIN out.gene_dispersion d ON d.ensg_id = g.gene_id
    """, ["gene_enriched", "out.gene_dispersion"]),

    # core_score with normalised keys alongside the originals
    "v_core_score": ("""
        SELECT s.*,
               lower(s.model_id) AS model_id_lc,
               lower(s.ensg_id)  AS gene_id
        FROM out.core_score s
    """, ["out.core_score"]),

    # the join this whole warehouse exists to make one query: a score, the gene
    # it belongs to, and everything harmonisation knows about the cell line
    "v_score_annotated": ("""
        SELECT s.model_id,
               s.ensg_id,
               s.core_score,
               s.n_layers,
               s.stratum_rank,
               g.hugo_symbol,
               g.gene_role,
               g.is_protein_coding,
               g.signal_spread,
               g.top_vs_median_fold,
               c.canonical_name,
               c.lineage,
               c.primary_disease,
               c.tissue,
               c.cancer_type,
               c.msi_status,
               c.n_modalities,
               c.procan_proteomics,
               c.cosmic_cna
        FROM v_core_score s
        LEFT JOIN v_gene g      ON g.gene_id  = s.gene_id
        LEFT JOIN v_cell_line c ON c.model_id = s.model_id_lc
    """, ["v_core_score", "v_gene", "v_cell_line"]),
}


def _exists(con: duckdb.DuckDBPyConnection, qualified: str) -> bool:
    if "." in qualified:
        schema, name = qualified.split(".", 1)
    else:
        schema, name = "main", qualified
    return con.execute(
        "SELECT count(*) FROM information_schema.tables "
        "WHERE table_schema = ? AND lower(table_name) = ?",
        [schema, name.lower()]).fetchone()[0] > 0


def build_bridges(con: duckdb.DuckDBPyConnection, verbose: bool = False) -> list[str]:
    """Create the key-normalising convenience views, skipping any whose inputs
    are absent (e.g. before 00b has run, or before Stage 2 has produced a score)."""
    made = []
    for name, (sql, deps) in BRIDGE_VIEWS.items():
        missing = [d for d in deps if not _exists(con, d)]
        if missing:
            print(f"  skip {name}: needs {', '.join(missing)}")
            continue
        con.execute(f'CREATE OR REPLACE VIEW "{name}" AS {sql}')
        made.append(name)
        if verbose:
            print(f"  main.{name:<40} <- {', '.join(deps)}")
    return made


def verify(con: duckdb.DuckDBPyConnection) -> pd.DataFrame:
    """Row/column counts for every object in the warehouse, materialised or not."""
    rows = []
    listing = con.execute(
        "SELECT table_schema, table_name, table_type FROM information_schema.tables "
        "WHERE table_schema IN ('main', '" + "', '".join(ZONES) + "') "
        "ORDER BY table_schema, table_name"
    ).fetchall()
    for schema, name, ttype in listing:
        try:
            n = con.execute(f'SELECT count(*) FROM "{schema}"."{name}"').fetchone()[0]
            c = len(con.execute(f'SELECT * FROM "{schema}"."{name}" LIMIT 0').description)
        except duckdb.Error:
            n, c = None, None
        rows.append({
            "schema": schema,
            "object": name,
            "kind": "table" if ttype == "BASE TABLE" else "view",
            "rows": n,
            "cols": c,
        })
    return pd.DataFrame(rows)


def main() -> None:
    ap = argparse.ArgumentParser(
        description="Expose the pipeline outputs as views inside the warehouse")
    ap.add_argument("--verbose", action="store_true",
                    help="print every view as it is registered")
    ap.add_argument("--verify", action="store_true",
                    help="print row/column counts for the whole warehouse")
    args = ap.parse_args()

    if not os.path.isfile(DB_PATH):
        raise SystemExit(
            f"{DB_PATH} not found.\n"
            "Run 00_harmonisation.ipynb (then 00b_enriched_harmonisation.ipynb) first — "
            "this script adds views to that warehouse, it does not create one.")

    con = duckdb.connect(DB_PATH)
    print(f"warehouse: {DB_PATH}")
    registered = build(con, verbose=args.verbose)
    bridges = build_bridges(con, verbose=args.verbose)
    total = sum(len(v) for v in registered.values())
    n_tables = con.execute(
        "SELECT count(*) FROM information_schema.tables "
        "WHERE table_schema = 'main' AND table_type = 'BASE TABLE'").fetchone()[0]
    print(f"\n{n_tables} materialised tables (main) + {total} views across "
          f"{len(registered)} schemas: "
          + ", ".join(f"{z}={len(v)}" for z, v in registered.items()))
    if bridges:
        print(f"bridge views (keys normalised): {', '.join(bridges)}")

    if args.verify:
        pd.set_option("display.max_rows", 300, "display.width", 160)
        print("\n" + verify(con).to_string(index=False))

    con.close()


if __name__ == "__main__":
    main()
