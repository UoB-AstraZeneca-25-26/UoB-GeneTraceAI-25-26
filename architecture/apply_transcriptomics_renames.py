"""
apply_transcriptomics_renames.py
----------------------------------
Applies the schema/data fixes needed to make the warehouse compatible
with the new transcriptomics.py (Triveni's branch). Idempotent -- every
step checks state before acting, skips if already done. Safe to rerun.

PART A -- 4 column renames:
    depmap_profiles.profileid -> profile_id
    depmap_expr.profileid     -> profile_id   (separate column, same old name)
    geo_expr.sample           -> gsm_id
    hpa_rna.gene              -> gene_id

Confirmed safe: the only other code touching these old names was
utils/common.py's gene_universe() and fetch_depmap(), both used
exclusively by rna_scorer.py (itself being retired in this same
migration).

PART B -- 2026-08-24: geo_expr duplicate-GSM row removal.

Two GSMs in geo_expr each had 2 rows sharing the same gsm_id, differing
ONLY in model_id (identical gsm_id and all 16,063 gene expression
values in both cases) -- a many-to-one join fan-out during geo_expr's
own model_id assignment, not corrupted expression data. Confirmed this
predates the transcriptomics migration; it just wasn't reachable until
this pipeline's align_axis_rows() started reindexing on gsm_id.

The new transcriptomics.py pipeline cannot ingest a raw id (gsm_id)
mapped to more than one model_id -- make_axis() assumes a 1:1 raw_id ->
model_id mapping, so align_axis_rows()'s reindex crashes on the
duplicate index. Each pair is resolved here, not dropped outright,
because the evidence for the surviving row is specific and one-sided in
both cases, not a coin flip -- see per-row justification below.

  GSM219757:
    kept   ach-000833 (rh-30 / rh30)
    dropped ach-001189 (sjrh30)
    why: geo_info resolves this GSM to ach-000833 confidently
         (n_model_id=1, is_ambiguous=False) -- the warehouse's own
         authoritative per-GSM resolver already disambiguated this.
         ach-001189 is not corroborated by geo_info at all.

  GSM1374458:
    kept   ach-002222 (ctv-1 / ctv1)
    dropped ach-001737 (ctv-1-dm / ctv1dm)
    why: geo_info could not auto-resolve this GSM to a model_id
         (n_model_id=0), but its own cell_line_trimmed value is an
         EXACT string match to ach-002222's stripped_cell_line_name
         ("ctv1" == "ctv1"). ach-001737 is a distinct, named subline
         ("-dm" suffix), not a name collision with the trimmed value.

Usage:
    python apply_transcriptomics_renames.py --db outputs/celllineselector.db
"""
import argparse

import duckdb

RENAMES = [
    ("depmap_profiles", "profileid", "profile_id"),
    ("depmap_expr", "profileid", "profile_id"),
    ("geo_expr", "sample", "gsm_id"),
    ("hpa_rna", "gene", "gene_id"),
]

# (gsm_id, model_id_to_drop, justification)
DUPLICATE_GSM_ROWS_TO_DROP = [
    ("gsm219757", "ach-001189",
     "geo_info confidently resolves this GSM to ach-000833 (rh-30) instead "
     "(n_model_id=1, is_ambiguous=False); ach-001189 (sjrh30) is uncorroborated."),
    ("gsm1374458", "ach-001737",
     "geo_info's cell_line_trimmed='ctv1' is an exact match to ach-002222's "
     "stripped name; ach-001737 (ctv-1-dm) is a distinct named subline."),
]


def apply_renames(con):
    for table, old, new in RENAMES:
        cols = con.execute(f"PRAGMA table_info('{table}')").fetchdf()["name"].tolist()
        if new in cols:
            print(f"  [skip] {table}.{new} already present")
            continue
        if old not in cols:
            print(f"  [WARN] {table}.{old} not found -- nothing to rename")
            continue
        con.execute(f'ALTER TABLE {table} RENAME COLUMN "{old}" TO "{new}"')
        print(f"  [done] {table}.{old} -> {new}")


def drop_duplicate_geo_rows(con):
    for gsm_id, model_id, why in DUPLICATE_GSM_ROWS_TO_DROP:
        n = con.execute(
            "SELECT count(*) FROM geo_expr WHERE gsm_id = ? AND model_id = ?",
            [gsm_id, model_id],
        ).fetchone()[0]
        if n == 0:
            print(f"  [skip] {gsm_id}/{model_id} already absent")
            continue
        con.execute(
            "DELETE FROM geo_expr WHERE gsm_id = ? AND model_id = ?",
            [gsm_id, model_id],
        )
        print(f"  [done] dropped {gsm_id}/{model_id} -- {why}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--db", required=True)
    args = ap.parse_args()

    con = duckdb.connect(args.db, read_only=False)
    try:
        print("-- renames --")
        apply_renames(con)
        print("-- duplicate GEO row removal --")
        drop_duplicate_geo_rows(con)
    finally:
        con.close()


if __name__ == "__main__":
    main()
