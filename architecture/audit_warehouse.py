"""
audit_warehouse.py
-------------------
READ-ONLY diagnostic pass over celllineselector.db.

Two independent jobs, run in sequence:
  PART A — species / non-human contamination check
  PART B — enrichment discovery: what columns already sit in the warehouse
           unused by any current pipeline output, that could be surfaced
           in the frontend (cell line cards, gene detail pages, etc.)

Nothing here writes to the DB. Every connection is opened read_only=True.
Run this from anywhere; it only needs the DB path.

Usage:
    python audit_warehouse.py --db /path/to/celllineselector.db
    python audit_warehouse.py            # uses config.DB if importable
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import duckdb
import pandas as pd

pd.set_option("display.max_rows", 200)
pd.set_option("display.width", 140)


def resolve_db_path(cli_arg: str | None) -> Path:
    if cli_arg:
        return Path(cli_arg)
    try:
        sys.path.insert(0, str(Path(__file__).resolve().parent))
        from config import DB  # type: ignore
        return DB
    except Exception:
        raise SystemExit(
            "No --db given and config.DB not importable. "
            "Run with: python audit_warehouse.py --db /full/path/to/celllineselector.db"
        )


def section(title: str):
    print("\n" + "=" * 78)
    print(title)
    print("=" * 78)


def safe_query(con, sql: str, label: str) -> pd.DataFrame | None:
    """Run a query; on failure (missing table/column), report and continue
    rather than crashing the whole audit."""
    try:
        return con.execute(sql).fetchdf()
    except Exception as e:
        print(f"  [skip] {label}: {e}")
        return None


def list_tables(con) -> list[str]:
    return con.execute(
        "SELECT table_name FROM information_schema.tables "
        "WHERE table_schema='main' AND table_type='BASE TABLE' "
        "ORDER BY table_name"
    ).fetchdf()["table_name"].tolist()


def list_columns(con, table: str) -> pd.DataFrame:
    return con.execute(f"PRAGMA table_info('{table}')").fetchdf()


# ======================================================================
# PART A — species / contamination audit
# ======================================================================

def part_a_species_audit(con, tables: list[str]):
    section("PART A — SPECIES / NON-HUMAN CONTAMINATION AUDIT")

    # A1. Find every column across the warehouse that might carry species info.
    print("\n--- A1. Scanning all tables for species/organism-like columns ---")
    candidate_cols = []
    for t in tables:
        cols = list_columns(con, t)
        hits = cols[cols["name"].str.lower().str.contains(
            "species|organism|taxon|host", regex=True, na=False)]
        for _, row in hits.iterrows():
            candidate_cols.append((t, row["name"], row["type"]))
    if candidate_cols:
        print(pd.DataFrame(candidate_cols, columns=["table", "column", "type"])
              .to_string(index=False))
    else:
        print("  No column across any base table matches "
              "species|organism|taxon|host by name.")
        print("  -> species was likely never captured explicitly during "
              "harmonisation. See A5 for an identity-based workaround.")

    # A2. If sample_info exists, profile it directly regardless of A1 hits —
    # some pipelines store species as a free-text field with a different name
    # (e.g. 'source', 'lineage' subtype, a Cellosaurus cross-ref column).
    if "sample_info" in tables:
        print("\n--- A2. sample_info: full column list (manual eyeball for species-like fields) ---")
        print(list_columns(con, "sample_info")[["name", "type"]].to_string(index=False))

        df = safe_query(con,
            "SELECT * FROM sample_info LIMIT 5", "sample_info sample rows")
        if df is not None:
            print("\n  sample_info — first 5 rows (all columns):")
            print(df.to_string(index=False))

    # A3. Direct species tally, if a species-like column was found in A1.
    for t, c, _ in candidate_cols:
        df = safe_query(con,
            f'SELECT "{c}" AS value, count(*) AS n FROM "{t}" '
            f'GROUP BY 1 ORDER BY 2 DESC', f"{t}.{c} tally")
        if df is not None:
            print(f"\n--- A3. {t}.{c} — value counts ---")
            print(df.to_string(index=False))

    # A4. GEO-specific check — GEO is the most likely entry point for
    # non-human material since deposits routinely mix organisms.
    if "geo_info" in tables:
        print("\n--- A4. geo_info: full column list ---")
        print(list_columns(con, "geo_info")[["name", "type"]].to_string(index=False))
        df = safe_query(con, "SELECT * FROM geo_info LIMIT 5", "geo_info sample")
        if df is not None:
            print("\n  geo_info — first 5 rows:")
            print(df.to_string(index=False))

    if "geo_platform" in tables:
        # characteristics / source_name / molecule fields from SOFT metadata
        # often contain the organism as free text even with no dedicated column.
        df = safe_query(con,
            "SELECT source_name, molecule, characteristics FROM geo_platform LIMIT 20",
            "geo_platform free-text fields")
        if df is not None:
            print("\n--- A4b. geo_platform — free-text fields that may reveal organism ---")
            print(df.to_string(index=False))
        df2 = safe_query(con, """
            SELECT
                CASE
                    WHEN lower(source_name) LIKE '%mouse%' OR lower(characteristics) LIKE '%mus musculus%'
                        THEN 'possible_mouse'
                    WHEN lower(source_name) LIKE '%rat%' OR lower(characteristics) LIKE '%rattus%'
                        THEN 'possible_rat'
                    WHEN lower(source_name) LIKE '%human%' OR lower(characteristics) LIKE '%homo sapiens%'
                        THEN 'looks_human'
                    ELSE 'unclear'
                END AS guess,
                count(*) AS n
            FROM geo_platform
            GROUP BY 1 ORDER BY 2 DESC
        """, "geo_platform organism keyword scan")
        if df2 is not None:
            print("\n--- A4c. geo_platform — keyword-based organism guess (free text, not authoritative) ---")
            print(df2.to_string(index=False))
            print("  NOTE: this is a crude keyword scan, not a real species field.")
            print("  Use it only to flag candidates for the Cellosaurus join in A5,")
            print("  never as the actual filter — free-text guesses will have false positives")
            print("  (e.g. a human line's characteristics mentioning a mouse xenograft origin).")

    # A5. RRID -> Cellosaurus is the correct authoritative check, per the
    # pipeline's own confirmed architecture (RRID-based identity resolution).
    # This just checks whether the join key exists and is populated — it does
    # NOT hit the network. Cross-referencing against a downloaded Cellosaurus
    # export is a separate offline step once RRID coverage is confirmed here.
    print("\n--- A5. RRID coverage check (for an offline Cellosaurus species join) ---")
    for t in tables:
        cols = list_columns(con, t)
        rrid_cols = cols[cols["name"].str.lower().str.contains("rrid", na=False)]
        for _, row in rrid_cols.iterrows():
            c = row["name"]
            df = safe_query(con, f"""
                SELECT
                    count(*) AS total_rows,
                    count("{c}") AS non_null_rrid,
                    count(DISTINCT "{c}") AS distinct_rrid
                FROM "{t}"
            """, f"{t}.{c} RRID coverage")
            if df is not None:
                print(f"  {t}.{c}:")
                print(f"    {df.to_string(index=False)}")
    print("\n  If RRID coverage above is high, the next offline step is:")
    print("    1. Download the current Cellosaurus release (cellosaurus.txt or .obo)")
    print("    2. Parse AC (accession/RRID) -> OX (NCBI TaxID) / OS (organism species)")
    print("    3. Left-join your RRID column to that mapping")
    print("    4. Anything not resolving to Homo sapiens (TaxID 9606) is a candidate")
    print("       for removal — but confirm each hit manually before deleting, since")
    print("       misassigned RRIDs are also a known failure mode distinct from true")
    print("       non-human contamination.")

    # A6. Cross-check cell line COUNT against the "original 14 datasets"
    # assumption — a blunt but useful sanity check on scale creep.
    if "sample_info" in tables:
        df = safe_query(con, "SELECT count(DISTINCT model_id) AS n_models FROM sample_info",
                        "distinct model_id count")
        if df is not None:
            print("\n--- A6. Total distinct model_id count in sample_info ---")
            print(df.to_string(index=False))
            print("  Compare this against known human cancer cell line panel sizes")
            print("  (e.g. DepMap ~1800, CCLE ~1000, GDSC ~1000) to sanity-check scale.")


# ======================================================================
# PART B — enrichment discovery for the frontend
# ======================================================================

def part_b_enrichment_scan(con, tables: list[str]):
    section("PART B — ENRICHMENT DISCOVERY (unused warehouse columns for frontend)")

    print("\n--- B1. Full schema inventory: every table x every column ---")
    all_cols = []
    for t in tables:
        cols = list_columns(con, t)
        n_rows = con.execute(f'SELECT count(*) FROM "{t}"').fetchone()[0]
        for _, row in cols.iterrows():
            all_cols.append((t, n_rows, row["name"], row["type"]))
    schema_df = pd.DataFrame(all_cols, columns=["table", "table_rows", "column", "type"])
    print(schema_df.to_string(index=False))

    # B2. Flag columns that are plainly metadata/descriptive rather than
    # numeric scoring data — these are the best frontend-enrichment candidates
    # (display strings, categorical tags, external IDs, links) as opposed to
    # z-scores/statistics which are already surfaced via the *_Z parquet outputs.
    print("\n--- B2. Likely display/enrichment columns (heuristic: text/categorical, "
          "not already a *_z / *_score output) ---")
    display_like = schema_df[
        schema_df["column"].str.lower().str.contains(
            "name|symbol|description|synonym|alias|url|link|image|disease|"
            "tissue|subtype|histology|sex|age|source|reference|pubmed|doi|"
            "vendor|catalog|supplier|comment|note", regex=True, na=False)
    ]
    if not display_like.empty:
        print(display_like.to_string(index=False))
    else:
        print("  No obviously descriptive columns found by name pattern — "
              "inspect schema_df above manually.")

    # B3. Sample a few rows from every table with likely-enrichment columns,
    # so it's obvious at a glance what's actually populated vs. mostly NULL.
    print("\n--- B3. Sample rows + null-rate for candidate enrichment tables ---")
    candidate_tables = display_like["table"].unique().tolist()
    for t in candidate_tables:
        cols = list_columns(con, t)["name"].tolist()
        print(f"\n  >> {t}")
        df = safe_query(con, f'SELECT * FROM "{t}" LIMIT 3', f"{t} sample")
        if df is not None:
            print(df.to_string(index=False))
        # null-rate per column, useful to know which "enrichment" fields are
        # actually usable versus sparsely populated and not worth surfacing
        null_rates = []
        for c in cols:
            r = safe_query(con, f"""
                SELECT
                    round(100.0 * count(*) FILTER (WHERE "{c}" IS NULL) / count(*), 1) AS pct_null
                FROM "{t}"
            """, f"{t}.{c} null rate")
            if r is not None:
                null_rates.append((c, r.iloc[0]["pct_null"]))
        if null_rates:
            nr_df = pd.DataFrame(null_rates, columns=["column", "pct_null"]).sort_values("pct_null")
            print(f"  null rates for {t}:")
            print(nr_df.to_string(index=False))

    # B4. Specifically check known reference tables that are commonly rich
    # in frontend-useful metadata but easy to forget once the numeric
    # pipeline stages are built: gene, sample_info, cell line lookups.
    print("\n--- B4. Targeted checks on likely-rich reference tables ---")
    for t in ["gene", "gene_roster", "sample_info", "geo_info"]:
        if t in tables:
            print(f"\n  >> {t} — all columns:")
            print(list_columns(con, t)[["name", "type"]].to_string(index=False))

    # B5. Compare warehouse schema against known pipeline OUTPUT columns
    # (from config.py) to make explicit what's raw-only vs. already exposed.
    print("\n--- B5. Cross-reference: pipeline output parquet files (from config.py) ---")
    print("  These already reach the frontend via existing outputs:")
    known_outputs = [
        "bulk_rna_z.parquet", "bulk_prot_z.parquet", "protein_platform_tier.parquet",
        "cna_flags.parquet", "mutations_collapsed.parquet", "mutations_scores.parquet",
        "fusions_gene_level.parquet", "fusions_scores.parquet", "core_score.parquet",
        "gene_dispersion.parquet", "flags_with_driver.parquet",
        "predictions_with_confidence.parquet", "evidence_ledger.parquet",
    ]
    for f in known_outputs:
        print(f"    - {f}")
    print("\n  Anything found in B1-B4 that is NOT one of the above outputs, and is")
    print("  not purely a join key or an intermediate numeric column already folded")
    print("  into one of the above, is a genuine enrichment candidate: metadata that")
    print("  exists in the warehouse today but never reaches any current output file.")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--db", default=None, help="Path to celllineselector.db")
    args = ap.parse_args()

    db_path = resolve_db_path(args.db)
    if not Path(db_path).exists():
        raise SystemExit(f"DB not found at {db_path}")

    con = duckdb.connect(str(db_path), read_only=True)
    try:
        tables = list_tables(con)
        section("WAREHOUSE OVERVIEW")
        print(f"DB path: {db_path}")
        print(f"Base tables ({len(tables)}): {tables}")

        part_a_species_audit(con, tables)
        part_b_enrichment_scan(con, tables)

        section("DONE")
        print("Nothing was written. This connection was read-only throughout.")
    finally:
        con.close()


if __name__ == "__main__":
    main()
