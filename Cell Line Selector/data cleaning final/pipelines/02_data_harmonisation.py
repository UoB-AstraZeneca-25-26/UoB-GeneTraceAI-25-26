"""
02_data_harmonisation.py

Harmonisation stage of the pipeline. Builds the two cross-dataset roster
tables from the cleaned Parquet files produced by ``01_data_cleaning.py``,
then uses one of them to attach ``model_id`` back onto every cleaned table.

Rosters built
-------------
cell_line_roster
    One row per ``model_id``, aggregating every ``cell_line_name`` /
    ``cvcl_id`` / ``ccle_id`` seen for that model across all cleaned
    tables, then widened in both directions against Cellosaurus
    (accession -> synonyms, then synonyms -> further accessions).
gene_roster
    One row per ``gene_id``, restricted to protein-coding genes per the
    mutations table, aggregating every ``gene_name`` / ``geo_gsm_id`` /
    ``profile_id`` seen for that gene, plus UniProt and COSMIC bridges.

After the rosters are built, ``model_id`` is attached to every cleaned
table that lacks one — using ``cell_line_roster`` as the lookup, with
extra name- and RRID-based paths for the protein matrix and model list —
and the updated tables are re-saved as both Parquet and CSV.

Design notes
------------
* Order matters in several places and is called out in comments at each
  site: Cellosaurus synonyms before CVCL expansion; ``build_gene_roster``
  before the proteomics backfill; the backfill before the COSMIC bridges.
* Each sub-step inside :func:`build_rosters` carries its own
  ``try``/``except``, so a failure names the specific step rather than
  silently discarding a whole roster.
* Errors are logged to ``logs/pipeline_<date>.log`` on the same pattern
  as ``01_data_cleaning.py``. Diagnostic CSVs (gene name mismatches,
  ambiguous ``model_id`` rows, attachment summary) are written to
  ``LOG_DIR`` alongside it.

Inputs
------
Cleaned Parquet files under ``PARQUET_CLEAN``, plus the raw ProCan TSV
under ``data/raw data/gene expression/`` if present.

Outputs
-------
``cell_line_roster`` and ``gene_roster`` written to ``PARQUET_CLEAN`` and
``CSV_CLEAN``, the re-saved cleaned tables with ``model_id`` attached, and
diagnostic CSVs plus a dated log file under ``LOG_DIR``.

Run with
--------
    python pipelines/02_data_harmonisation.py
"""

import sys
from datetime import datetime
from pathlib import Path

project_root = Path(__file__).resolve().parent
for candidate in [project_root, *project_root.parents]:
    if (candidate / "src").exists():
        project_root = candidate
        break
sys.path.insert(0, str(project_root))

import pandas as pd

from src.scripts.data_utils import PARQUET_CLEAN, CSV_CLEAN
from src.scripts.harmonisation_functions import (
    build_cell_line_roster,
    add_cellosaurus_synonyms_to_cell_line_roster,
    expand_cvcl_ids_from_names,
    build_gene_roster,
    build_combined_proteomics_gene_map,
    add_missing_proteomics_genes_to_gene_roster,
    filter_protein_coding_expression_tables,
    check_gene_name_consistency,
    attach_model_id_everywhere,
    attach_model_id_to_protein_matrix,
    attach_model_id_to_model_list,
    summarize_model_id_attachment,
    add_cosmic_cnv_id_to_gene_roster,
    add_gene_id_to_cosmic_cna,
    add_cosmic_gene_id_to_gene_roster,
    add_gene_id_to_cosmic_gene_census,
)
from src.scripts.export_data import save_cleaned_parquets, save_cleaned_csvs
from src.scripts.logging_utils import get_logger, log_error, LOG_DIR

logger = get_logger("02_build_rosters")


def load_cleaned_tables() -> dict:
    """
    Load every cleaned Parquet file in ``PARQUET_CLEAN`` into a dict.

    Globs ``*.parquet`` in sorted order and keys each table by its
    filename stem, so the dict keys match the dataset keys used
    throughout the cleaning stage.

    Returns
    -------
    dict
        Mapping of dataset key (str) to :class:`pandas.DataFrame`.

    Raises
    ------
    Exception
        Any read failure is logged with the source directory and
        re-raised. This one is fatal: with no cleaned tables there is
        nothing to harmonise.
    """
    print(f"Loading cleaned parquet files from {PARQUET_CLEAN} ...")
    tables = {}
    try:
        for path in sorted(PARQUET_CLEAN.glob("*.parquet")):
            tables[path.stem] = pd.read_parquet(path)
    except Exception as e:
        log_error(logger, step="load_cleaned_tables", error=e, source_dir=str(PARQUET_CLEAN))
        raise
    print(f"Loaded {len(tables)} cleaned tables: {', '.join(sorted(tables))}")
    return tables


def check_cell_line_consistency(direct: pd.DataFrame, max_expected_models: int = 2) -> pd.DataFrame:
    """
    Flag model_id/identifier combinations that look like genuine mapping problems.

    Distinguishes real mapping errors from the normal multi-value cases
    that ``cell_line_roster``'s aggregation exists to handle. Two classes
    of issue are reported:

    * a ``cvcl_id`` or ``ccle_id`` attached to more than
      ``max_expected_models`` distinct ``model_id`` values. Up to 2 is
      treated as normal — a cell line re-derived or re-registered under a
      second ``model_id`` — and is not flagged. Beyond that is the
      signature of a likely bad mapping in one of the source tables.
    * a ``model_id`` with no ``cell_line_name``, ``cvcl_id`` *or*
      ``ccle_id`` found in any table (thin or absent coverage).

    A ``model_id`` carrying several ``cell_line_name`` or ``cvcl_id``
    values is *not* flagged — aliases and multiple accessions are
    expected.

    Parameters
    ----------
    direct : pandas.DataFrame
        Long-form identifier table, one row per identifier occurrence,
        with columns ``model_id``, ``cell_line_name``, ``cvcl_id``,
        ``ccle_id`` and ``_source_table``.
    max_expected_models : int, optional
        Number of distinct ``model_id`` values an accession may map to
        before being flagged. Default 2.

    Returns
    -------
    pandas.DataFrame
        One row per issue with columns ``issue``, ``value``,
        ``n_model_ids``, ``model_ids`` and ``source_tables``. Empty if
        nothing was flagged.

    Warnings
    --------
    **Not currently called from** :func:`main`. It requires ``direct`` to
    carry a ``_source_table`` column identifying which raw table each row
    came from, and :func:`build_cell_line_roster` does not currently
    produce one. Wire this in once that tracking is added to
    ``roster_utils.py`` / ``harmonisation_functions.py``; calling it as-is
    against the current roster output will raise ``KeyError`` on
    ``_source_table``.
    """
    issues = []

    for id_col in ["cvcl_id", "ccle_id"]:
        if id_col not in direct.columns:
            continue
        sub = direct.dropna(subset=[id_col, "model_id"])
        counts = sub.groupby(id_col)["model_id"].nunique()
        conflicted = counts[counts > max_expected_models]
        for val, n in conflicted.items():
            model_ids = sorted(sub.loc[sub[id_col] == val, "model_id"].unique())
            issues.append({
                "issue": f"{id_col}_maps_to_multiple_model_ids",
                "value": val,
                "n_model_ids": n,
                "model_ids": "; ".join(model_ids),
                "source_tables": "; ".join(sorted(sub.loc[sub[id_col] == val, "_source_table"].unique())),
            })

    all_model_ids = set(direct["model_id"].dropna().unique())
    has_any_id = direct.dropna(subset=["cell_line_name", "cvcl_id", "ccle_id"], how="all")
    covered = set(has_any_id["model_id"].dropna().unique())
    for model_id in sorted(all_model_ids - covered):
        issues.append({
            "issue": "no_identifiers_found",
            "value": model_id,
            "n_model_ids": None,
            "model_ids": model_id,
            "source_tables": "",
        })

    return pd.DataFrame(issues)


def run_diagnostics(tables: dict) -> None:
    """
    Cross-check gene_id/gene_name consistency between sources before building rosters.

    Compares the gene names attached to each ``gene_id`` in the mutations
    and fusions tables against those in ``hpa_rna``, and reports any
    ``gene_id`` whose name sets are entirely disjoint — a sign that two
    sources disagree about what a given Ensembl ID refers to, which would
    otherwise be silently merged by ``gene_roster``'s aggregation.

    Parameters
    ----------
    tables : dict
        Mapping of dataset key to cleaned :class:`pandas.DataFrame`, as
        returned by :func:`load_cleaned_tables`.

    Returns
    -------
    None
        Findings are printed to stdout; any mismatches are written to
        ``LOG_DIR/gene_name_mismatches.csv``.

    Notes
    -----
    Advisory only — nothing downstream branches on the result. A failure
    is logged and skipped so it cannot block roster building.
    """
    print("\nChecking gene_id/gene_name consistency (mutations+fusions vs hpa_rna)...")
    try:
        mismatches = check_gene_name_consistency(tables)
        if mismatches.empty:
            print("  No mismatches found.")
        else:
            print(f"  {len(mismatches)} gene_id(s) with disjoint names — see logs/gene_name_mismatches.csv")
            mismatches.to_csv(LOG_DIR / "gene_name_mismatches.csv", index=False)
    except Exception as e:
        log_error(logger, step="check_gene_name_consistency", error=e)
        print("  [SKIPPED] consistency check failed — see log")

def build_rosters(tables: dict) -> dict:
    """
    Build both roster tables, in the order their dependencies require.

    Runs five phases:

    1. **cell_line_roster.** Aggregates identifiers per ``model_id``, then
       widens them against Cellosaurus in both directions — accession to
       synonyms first, then the widened name set back out to further
       accessions. Running the second step first would miss every
       synonym-derived accession.
    2. **Protein map check.** ``01_data_cleaning`` already writes
       ``proteomics_gene_map`` with ``procan_uniprot_ids`` merged in, so
       the map is used as-is when that column is present, and a warning is
       printed if the column exists but is entirely null. Only when the
       column is missing is the map rebuilt — see the inline comment for
       why an unconditional rebuild would be destructive.
    3. **Raw ProCan header rows.** The first two rows of the raw ProCan
       TSV (symbol and UniProt ID) are loaded as ``raw_procan`` to give
       ``build_gene_roster`` a direct bridge. If the file is absent,
       ``procan_uniprot_id`` falls back to ``proteomics_gene_map``.
    4. **gene_roster.** Built, then used to filter the expression tables
       down to protein-coding ENSG IDs in place, followed by a per-column
       coverage report — an all-null column here means a bridge matched
       nothing, which is otherwise easy to miss.
    5. **Proteomics backfill and COSMIC bridges.** The backfill must run
       after ``build_gene_roster`` but before the COSMIC bridges, so
       proteomics-only genes also get a chance to match a COSMIC symbol.
       Each COSMIC bridge is two-way: an ID onto ``gene_roster``, and
       ``gene_id`` back onto the COSMIC table.

    Parameters
    ----------
    tables : dict
        Mapping of dataset key to cleaned :class:`pandas.DataFrame`.
        **Mutated in place:** the expression tables are replaced with
        their protein-coding-only versions, ``raw_procan`` is added, the
        COSMIC tables gain a ``gene_id`` column, and
        ``proteomics_gene_map`` may be rebuilt.

    Returns
    -------
    dict
        Mapping with keys ``cell_line_roster`` and ``gene_roster``. A
        roster that failed to build is absent, and the steps that depend
        on it are skipped rather than raising.

    Notes
    -----
    Each sub-step carries its own ``try``/``except`` deliberately. Under a
    single shared block, a failure anywhere in the ``gene_roster`` chain
    discarded ``gene_roster`` entirely with no indication of which step
    broke.
    """
    rosters = {}

    # ── cell_line_roster ─────────────────────────────────────────────
    print("\nBuilding cell_line_roster...")
    try:
        rosters["cell_line_roster"] = build_cell_line_roster(tables)

        if "cellosaurus" in tables:
            # CVCL -> names: pull Cellosaurus synonyms in via the accession.
            rosters["cell_line_roster"] = add_cellosaurus_synonyms_to_cell_line_roster(
                rosters["cell_line_roster"],
                tables["cellosaurus"],
            )
            print("  Added Cellosaurus synonyms to cell_line_roster via CVCL match.")

            # names -> CVCL: reverse direction, using the widened name set
            # from the step above. Order matters; running this first would
            # miss every synonym-derived accession.
            rosters["cell_line_roster"] = expand_cvcl_ids_from_names(
                rosters["cell_line_roster"],
                tables["cellosaurus"],
            )

        print(f"  cell_line_roster: {rosters['cell_line_roster'].shape[0]:,} rows")
    except Exception as e:
        log_error(logger, step="build_cell_line_roster", error=e)
        print("  [SKIPPED] cell_line_roster: build failed — see log")

    # ── protein map ──────────────────────────────────────────────────
    # 01_data_cleaning writes proteomics_gene_map with procan_uniprot_ids
    # already merged in. build_combined_proteomics_gene_map() expects a
    # SEPARATE tables['procan_gene_map'], which 01 never creates -- calling
    # it here would overwrite the working map with an all-null procan
    # column. Only rebuild if the existing map is missing that column.
    print("\nChecking protein map...")
    try:
        existing = tables.get("proteomics_gene_map")
        if existing is not None and "procan_uniprot_ids" in existing.columns:
            n_procan = existing["procan_uniprot_ids"].notna().sum()
            print(f"  proteomics_gene_map already has procan_uniprot_ids "
                  f"({n_procan:,} / {len(existing):,} non-null) — using as-is.")
            if n_procan == 0:
                print("  [WARNING] procan_uniprot_ids is entirely null — the ProCan "
                      "bridge in add_procan_uniprot_ids_to_proteomics_map found no "
                      "matches. procan_uniprot_id will be empty in gene_roster.")
        else:
            tables["proteomics_gene_map"] = build_combined_proteomics_gene_map(tables)
            print(f"  protein map rebuilt: {tables['proteomics_gene_map'].shape[0]:,} rows, "
                  f"columns: {list(tables['proteomics_gene_map'].columns)}")
    except Exception as e:
        log_error(logger, step="build_combined_proteomics_gene_map", error=e)
        print("  [SKIPPED] protein map check failed — see log. gene_roster's "
              "uniprot_id(s)/procan_uniprot_id and the backfill below will be "
              "incomplete as a result.")

    # ── raw ProCan header rows (symbol + uniprot_id) for direct bridge ──
    _procan_path = project_root / "data" / "raw data" / "gene expression" / "Protein_matrix_averaged_20250211.tsv"
    if _procan_path.exists():
        tables["raw_procan"] = pd.read_csv(_procan_path, sep="\t", header=None, low_memory=False, nrows=2)
        print(f"  raw_procan header rows loaded ({tables['raw_procan'].shape[1]:,} columns)")
    else:
        print("  [WARNING] raw ProCan TSV not found — procan_uniprot_id will fall back to proteomics_gene_map")

    # ── gene_roster ──────────────────────────────────────────────────
    print("\nBuilding gene_roster...")
    try:
        rosters["gene_roster"] = build_gene_roster(tables)
        print(f"  gene_roster: {rosters['gene_roster'].shape[0]:,} rows")

        # Keep only protein-coding ENSG entries in the expression tables,
        # so non-coding genes are removed before downstream model building.
        filtered = filter_protein_coding_expression_tables(tables, rosters["gene_roster"])
        for key, value in filtered.items():
            tables[key] = value
        print("  expression tables filtered to protein-coding ENSG IDs")

        # Coverage report -- an all-null column here means a bridge
        # silently matched nothing, which is easy to miss otherwise.
        for col in ["gene_name(s)", "geo_gsm_id(s)", "profile_id(s)",
                    "uniprot_id(s)", "procan_uniprot_id"]:
            if col in rosters["gene_roster"].columns:
                n = rosters["gene_roster"][col].notna().sum()
                flag = "  <-- EMPTY" if n == 0 else ""
                print(f"    {col:<22} {n:>8,} non-null{flag}")
            else:
                print(f"    {col:<22} {'MISSING':>8}  <-- column not created")
    except Exception as e:
        log_error(logger, step="build_gene_roster", error=e)
        print("  [SKIPPED] gene_roster: build failed — see log")

    # ── proteomics backfill ──────────────────────────────────────────
    # Must run AFTER build_gene_roster but BEFORE the COSMIC bridges, so
    # newly-added proteomics-only genes also get a chance to match a
    # COSMIC symbol.
    if "gene_roster" in rosters and "proteomics_gene_map" in tables:
        print("\nBackfilling proteomics-only genes into gene_roster...")
        try:
            before = rosters["gene_roster"].shape[0]
            rosters["gene_roster"] = add_missing_proteomics_genes_to_gene_roster(
                rosters["gene_roster"],
                tables["proteomics_gene_map"],
            )
            after = rosters["gene_roster"].shape[0]
            print(f"  gene_roster after backfill: {after:,} rows (+{after - before:,})")
        except Exception as e:
            log_error(logger, step="add_missing_proteomics_genes_to_gene_roster", error=e)
            print("  [SKIPPED] proteomics backfill failed — see log")
    elif "gene_roster" in rosters:
        print("\nSkipping proteomics backfill — no protein map available.")

    # ── COSMIC bridges ───────────────────────────────────────────────
    if "gene_roster" in rosters:
        print("\nBridging gene_roster to COSMIC tables...")

        cna_key = "celllinesproject_completecna_v104_grch37"
        if cna_key in tables:
            try:
                rosters["gene_roster"] = add_cosmic_cnv_id_to_gene_roster(
                    rosters["gene_roster"], tables[cna_key],
                )
                tables[cna_key] = add_gene_id_to_cosmic_cna(
                    tables[cna_key], rosters["gene_roster"],
                )
                n = rosters["gene_roster"]["cosmic_cnv_id"].notna().sum()
                print(f"  cosmic_cnv_id -> gene_roster ({n:,} non-null), "
                      f"gene_id -> {cna_key}")
            except Exception as e:
                log_error(logger, step="cosmic_cna_bridge", error=e)
                print(f"  [SKIPPED] {cna_key} bridge failed — see log")

        census_key = "cosmic_cancergenecensus_v104_grch37"
        if census_key in tables:
            try:
                rosters["gene_roster"] = add_cosmic_gene_id_to_gene_roster(
                    rosters["gene_roster"], tables[census_key],
                )
                tables[census_key] = add_gene_id_to_cosmic_gene_census(
                    tables[census_key], rosters["gene_roster"],
                )
                n = rosters["gene_roster"]["cosmic_gene_id(s)"].notna().sum()
                print(f"  cosmic_gene_id(s) -> gene_roster ({n:,} non-null), "
                      f"gene_id -> {census_key}")
            except Exception as e:
                log_error(logger, step="cosmic_census_bridge", error=e)
                print(f"  [SKIPPED] {census_key} bridge failed — see log")

    return rosters


def spot_check_rosters(rosters: dict) -> None:
    """
    Print sample rows and basic coverage stats for each built roster.

    Shows the head of each roster, the ``model_id`` count for
    ``cell_line_roster``, and the number of ``gene_id`` entries left with
    an empty ``gene_name(s)`` for ``gene_roster``.

    Parameters
    ----------
    rosters : dict
        Mapping of roster name to :class:`pandas.DataFrame`, as returned
        by :func:`build_rosters`. Rosters that failed to build are absent
        and are silently skipped.

    Returns
    -------
    None
        Display only — output goes to stdout and nothing is written or
        mutated.
    """
    print("\nSpot checks:")

    if "cell_line_roster" in rosters:
        df = rosters["cell_line_roster"]
        print("\n-- cell_line_roster --")
        print(df.head())
        print(f"  model_ids: {df.shape[0]:,}")

    if "gene_roster" in rosters:
        df = rosters["gene_roster"]
        print("\n-- gene_roster --")
        print(df.head())
        empty_names = (df["gene_name(s)"] == "").sum()
        print(f"  gene_ids with no gene_name: {empty_names:,} / {df.shape[0]:,}")


def attach_ids_and_resave(tables: dict, rosters: dict) -> dict:
    """
    Attach model_id onto every cleaned table and re-save the results.

    Uses ``cell_line_roster`` as the lookup, via
    :func:`attach_model_id_everywhere`, then adds two table-specific
    paths for sources that carry no direct identifier match:

    * the protein matrix, matched on ``ccle_name`` against the roster's
      ``cell_line_name(s)``;
    * the model list, matched on ``rrid`` against ``cvcl_id(s)``.

    A before/after comparison is written as an attachment summary, any
    rows whose identifiers matched more than one ``model_id`` are logged,
    and every updated table is re-saved to both Parquet and CSV.

    Parameters
    ----------
    tables : dict
        Mapping of dataset key to cleaned :class:`pandas.DataFrame`.
    rosters : dict
        Mapping of roster name to :class:`pandas.DataFrame`; must contain
        ``cell_line_roster`` for anything to happen.

    Returns
    -------
    dict
        The updated tables, with ``model_id`` attached where a match was
        found. Returned unchanged if ``cell_line_roster`` was not built,
        or if attachment failed.

    Notes
    -----
    Writes ``LOG_DIR/model_id_attachment_summary.csv`` always, and
    ``LOG_DIR/model_id_ambiguous.csv`` when ambiguous rows exist. Save
    failures are collected across both formats and logged as one warning
    rather than raised.
    """
    if "cell_line_roster" not in rosters:
        print("\nSkipping model_id attachment — cell_line_roster wasn't built.")
        return tables

    print("\nAttaching model_id back to all tables...")
    try:
        tables_before = {name: df.copy() for name, df in tables.items()}
        tables, model_id_flags = attach_model_id_everywhere(tables, rosters["cell_line_roster"])

        protein_key = next((name for name in tables if "protein_matrix" in name.lower()), None)
        if protein_key is not None:
            tables[protein_key] = attach_model_id_to_protein_matrix(
                tables[protein_key],
                rosters["cell_line_roster"],
                protein_name_col="ccle_name",
                roster_name_col="cell_line_name(s)",
                output_col="model_id",
            )
            print(f"  Added model_id to {protein_key} via ccle_name -> cell_line_name(s) match.")

        model_list_key = next((name for name in tables if "model_list" in name.lower()), None)
        if model_list_key is not None:
            tables[model_list_key] = attach_model_id_to_model_list(
                tables[model_list_key],
                rosters["cell_line_roster"],
                rrid_col="rrid",
                roster_cvcl_col="cvcl_id(s)",
                output_col="model_id",
            )
            print(f"  Added model_id to {model_list_key} via rrid -> cvcl_id(s) match.")

        summary = summarize_model_id_attachment(tables_before, tables)
        print("\nmodel_id attachment summary:")
        print(summary.to_string(index=False))
        summary.to_csv(LOG_DIR / "model_id_attachment_summary.csv", index=False)

        if model_id_flags.empty:
            print("\n  No ambiguous identifier rows found.")
        else:
            print(f"\n  {len(model_id_flags)} row(s) with ambiguous model_id — see logs/model_id_ambiguous.csv")
            model_id_flags.to_csv(LOG_DIR / "model_id_ambiguous.csv", index=False)

        print(f"\nRe-saving updated tables (with model_id attached) to {PARQUET_CLEAN} ...")
        parquet_failed = save_cleaned_parquets(tables, out_dir=PARQUET_CLEAN)
        print(f"Re-saving updated tables (with model_id attached) to {CSV_CLEAN} ...")
        csv_failed = save_cleaned_csvs(tables, out_dir=CSV_CLEAN)

        failed = sorted(set(parquet_failed) | set(csv_failed))
        if failed:
            log_file = LOG_DIR / f"pipeline_{datetime.now():%Y-%m-%d}.log"
            logger.warning(
                f"{len(failed)} table(s) had a save failure after attach_model_id: "
                f"{', '.join(failed)} — see {log_file}"
            )
    except Exception as e:
        log_error(logger, step="attach_model_id_everywhere", error=e)
        print("  [SKIPPED] attach_model_id_everywhere failed — see log")

    return tables


def save_rosters(rosters: dict) -> None:
    """
    Save the roster DataFrames to both Parquet and CSV.

    Writes to ``PARQUET_CLEAN`` and ``CSV_CLEAN``, the same directories as
    the cleaned tables, so downstream stages pick the rosters up from a
    single location.

    Parameters
    ----------
    rosters : dict
        Mapping of roster name to :class:`pandas.DataFrame`, as returned
        by :func:`build_rosters`.

    Returns
    -------
    None
        Files and log entries are written as side effects. Returns early
        if no rosters were built.

    Notes
    -----
    Failures across the two formats are unioned and logged as a single
    warning naming the dated log file, so a roster that failed in only one
    format is still surfaced. Nothing raises.
    """
    if not rosters:
        print("\nNo rosters to save.")
        return

    print(f"\nSaving roster parquet files to {PARQUET_CLEAN} ...")
    parquet_failed = save_cleaned_parquets(rosters, out_dir=PARQUET_CLEAN)

    print(f"\nSaving roster CSV files to {CSV_CLEAN} ...")
    csv_failed = save_cleaned_csvs(rosters, out_dir=CSV_CLEAN)

    failed = sorted(set(parquet_failed) | set(csv_failed))
    if failed:
        log_file = LOG_DIR / f"pipeline_{datetime.now():%Y-%m-%d}.log"
        logger.warning(f"{len(failed)} roster(s) had a save failure: {', '.join(failed)} — see {log_file}")
    print("Done.")


def main() -> None:
    """
    Run the full harmonisation stage end to end.

    Loads the cleaned tables, runs the gene-name consistency diagnostic,
    builds both rosters, spot-checks them, attaches ``model_id`` back onto
    the cleaned tables and re-saves them, then saves the rosters
    themselves.

    Returns
    -------
    None
        Rosters, updated tables, diagnostic CSVs and log entries are
        written as side effects; summaries are printed to stdout.

    Raises
    ------
    Exception
        Propagated from :func:`load_cleaned_tables` if the cleaned Parquet
        files cannot be read. Every later failure is per-step: logged and
        skipped rather than raised.

    Notes
    -----
    Order is load -> diagnose -> build -> attach -> save. The tables are
    re-saved inside :func:`attach_ids_and_resave`, which also picks up the
    in-place changes :func:`build_rosters` made to the expression and
    COSMIC tables.
    """
    tables = load_cleaned_tables()
    run_diagnostics(tables)
    rosters = build_rosters(tables)
    spot_check_rosters(rosters)
    tables = attach_ids_and_resave(tables, rosters)
    save_rosters(rosters)


if __name__ == "__main__":
    main()