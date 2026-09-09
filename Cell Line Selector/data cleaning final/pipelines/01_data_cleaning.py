"""
01_data_cleaning.py

Data cleaning stage of the pipeline (converted from ``01_data_cleaning.ipynb``).

Loads the raw Parquet files written by ``00_data_loading.py``, applies a
per-dataset cleaning function to each, purges non-human cell lines, builds
the proteomics gene<->UniProt lookup, inspects and spot-checks the cleaned
tables, summarises missing values, and writes the cleaned tables back out
as both Parquet and CSV.

Design notes
------------
* Failures are per-dataset, not fatal. Anything that raises during
  cleaning, inspection or saving is logged with the dataset name, source
  file, step and full traceback, then skipped, so one bad file does not
  stop the rest of the run.
* Errors go to ``logs/pipeline_<date>.log``, giving enough context to fix
  the specific offending file without re-running the whole pipeline.
* Non-human filtering matches on CVCL accession only, never on cell-line
  name — see :func:`clean_data`.

Inputs
------
Raw Parquet files under ``PARQUET_RAW``, plus the raw ProCan TSV under
``data/raw data/gene expression/`` if present.

Outputs
-------
Cleaned tables under ``PARQUET_CLEAN`` and ``CSV_CLEAN``, plus a dated log
file under ``LOG_DIR``.

Run with
--------
    python pipelines/01_data_cleaning.py
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

from src.scripts.data_utils import load_raw_parquets, PARQUET_RAW, PARQUET_CLEAN, CSV_CLEAN
from src.scripts.cleaning_functions import (
    CLEANERS,
    base_clean,
    clean_cellosaurus,
    filter_non_human_rows,
    filter_geo_expr_by_kept_gsms,
    extract_proteomics_gene_map,
    add_procan_uniprot_ids_to_proteomics_map,
)
from src.scripts.export_data import save_cleaned_parquets, save_cleaned_csvs
from src.scripts.logging_utils import get_logger, log_file_error, log_error, LOG_DIR

logger = get_logger("01_data_cleaning")


def load_data() -> dict:
    """
    Load every raw Parquet file into a dict of DataFrames.

    Delegates to :func:`load_raw_parquets`, which reads each
    ``<key>.parquet`` under ``PARQUET_RAW`` produced by the data loading
    stage.

    Returns
    -------
    dict
        Mapping of dataset key (str) to :class:`pandas.DataFrame`.

    Raises
    ------
    Exception
        Any load failure is logged with the source directory and
        re-raised. Unlike the per-dataset failures elsewhere in this
        stage, this one is fatal: with no raw tables there is nothing
        to clean.
    """
    print("Loading raw parquet files...")
    try:
        tables = load_raw_parquets()
    except Exception as e:
        log_error(logger, step="load_data", error=e, source_dir=str(PARQUET_RAW))
        raise
    print(f"Loaded {len(tables)} raw tables.")
    return tables


def clean_data(tables: dict) -> dict:
    """
    Clean every dataset, purge non-human cell lines, and build the proteomics map.

    Runs in four phases:

    1. **Cellosaurus first.** :func:`clean_cellosaurus` must run before
       anything else because it populates ``NON_HUMAN_CVCLS`` as a
       side-effect — the blocklist every other table is filtered against.
       If cellosaurus is missing, a warning is printed and no blocklist
       is built.
    2. **Per-dataset cleaning.** Each remaining table is cleaned by its
       own cleaner from ``CLEANERS``, falling back to
       :func:`base_clean`, then passed through
       :func:`filter_non_human_rows`. Matching is on CVCL accession only,
       never on cell-line name — see the ``NON_HUMAN_CVCLS`` comment in
       ``cleaning_functions.py`` for the measured reason. Tables with no
       CVCL column pass through unfiltered; their non-human rows are
       excluded downstream by the ``cell_line_roster`` ``model_id`` join.
    3. **GEO alignment.** ``geo_expr`` has no species or CVCL column of
       its own, so it is restricted to the GSMs still present in cleaned
       ``geo_info``.
    4. **Proteomics gene map.** Builds the gene<->UniProt lookup from the
       proteomics table and enriches it with ProCan UniProt IDs, read
       from the raw ProCan TSV on disk if present, otherwise from a
       matching key already in ``tables``. If neither is available the
       map is still produced, without ``procan_uniprot_ids``.

    Parameters
    ----------
    tables : dict
        Mapping of dataset key to raw :class:`pandas.DataFrame`, as
        returned by :func:`load_data`.

    Returns
    -------
    dict
        Mapping of dataset key to cleaned :class:`pandas.DataFrame`.
        Datasets that failed to clean are absent. Includes the extra key
        ``proteomics_gene_map`` when the proteomics table was present and
        the map built successfully.

    Notes
    -----
    Per-dataset failures are logged with the source Parquet path, dataset
    name, cleaner name and input shape, reported on stdout as
    ``[SKIPPED]``, and the loop continues. The running count of dropped
    non-human rows is printed as a total.
    """
    cleaned = {}

    print("Applying cleaning functions...")
    if "cellosaurus" in tables:
        try:
            cleaned["cellosaurus"] = clean_cellosaurus(tables["cellosaurus"])
        except Exception as e:
            log_file_error(
                logger, file_path=PARQUET_RAW / "cellosaurus.parquet",
                step="clean_data", error=e, dataset="cellosaurus",
                cleaner="clean_cellosaurus", shape=str(tables["cellosaurus"].shape),
            )
            print("  [SKIPPED] cellosaurus: cleaning failed - see log")
    else:
        print("  [WARNING] 'cellosaurus' not in tables - no blocklist built.")

    print("\nChecking for non-human cell lines in each table...")
    total_dropped = 0

    for name, df in tables.items():
        if name == "cellosaurus":
            continue
        cleaner = CLEANERS.get(name, base_clean)
        try:
            out = cleaner(df)
            out, n_dropped = filter_non_human_rows(name, out)
            total_dropped += n_dropped
            cleaned[name] = out
        except Exception as e:
            log_file_error(
                logger, file_path=PARQUET_RAW / f"{name}.parquet",
                step="clean_data", error=e, dataset=name,
                cleaner=cleaner.__name__, shape=str(df.shape),
            )
            print(f"  [SKIPPED] {name}: failed - see log")

    print(f"\n  TOTAL non-human rows dropped: {total_dropped:,}")

    # geo_expr has no species or CVCL column of its own -- restrict it to
    # the GSMs still present in cleaned geo_info.
    if "geo_expr" in cleaned and "geo_info" in cleaned:
        print("\nAligning geo_expr to cleaned geo_info...")
        cleaned["geo_expr"], n_geo = filter_geo_expr_by_kept_gsms(
            cleaned["geo_expr"], cleaned["geo_info"]
        )
        total_dropped += n_geo

    if "proteomics" in tables:
        try:
            proteomics_map = extract_proteomics_gene_map(tables["proteomics"])
            raw_procan_path = (
                project_root / "data" / "raw data" / "gene expression"
                / "Protein_matrix_averaged_20250211.tsv"
            )
            if raw_procan_path.exists():
                raw_procan_df = pd.read_csv(
                    raw_procan_path, sep="\t", header=None, low_memory=False
                )
            else:
                procan_key = next(
                    (k for k in tables
                     if "protein_matrix_averaged_20250211" in k
                     or k in {"procan", "procan_tsv"}),
                    None,
                )
                raw_procan_df = tables[procan_key] if procan_key is not None else None

            if raw_procan_df is not None:
                proteomics_map = add_procan_uniprot_ids_to_proteomics_map(
                    proteomics_map, raw_procan_df
                )
            else:
                print("  [WARNING] ProCan source not found - "
                      "proteomics_gene_map will have no procan_uniprot_ids.")

            cleaned["proteomics_gene_map"] = proteomics_map
        except Exception as e:
            log_file_error(
                logger, file_path=PARQUET_RAW / "proteomics.parquet",
                step="clean_data", error=e, dataset="proteomics_gene_map",
                cleaner="extract_proteomics_gene_map",
                shape=str(tables["proteomics"].shape),
            )
            print("  [SKIPPED] proteomics_gene_map: build failed - see log")

    print(f"\nCleaned {len(cleaned)} / {len(tables)} datasets.")
    return cleaned


def inspect_shapes(cleaned: dict) -> None:
    """
    Print row and column counts for every cleaned table.

    A quick eyeball check that cleaning removed what was expected and
    nothing collapsed to zero rows.

    Parameters
    ----------
    cleaned : dict
        Mapping of dataset key to cleaned :class:`pandas.DataFrame`, as
        returned by :func:`clean_data`.

    Returns
    -------
    None
        Output is printed to stdout.
    """
    print("\nDataset shapes:")
    for name, df in cleaned.items():
            print(f"  {name:40s}: {df.shape[0]:>10,} rows x {df.shape[1]:>6,} cols")


def spot_check(cleaned: dict) -> None:
    """
    Print sample rows from key datasets as a sanity check on the cleaning.

    Shows the head of each dataset named in the internal ``checks`` map.
    Where a column list is given (currently only ``fusions``), only those
    columns are shown, which verifies that the expected renamed columns
    survived cleaning.

    Parameters
    ----------
    cleaned : dict
        Mapping of dataset key to cleaned :class:`pandas.DataFrame`, as
        returned by :func:`clean_data`.

    Returns
    -------
    None
        Output is printed to stdout.

    Notes
    -----
    Datasets absent from ``cleaned`` are silently skipped. A failure on
    one dataset's check — for example an expected column missing after a
    schema change upstream — is logged and skipped rather than crashing
    the run.
    """
    print("\nSpot checks:")

    checks = {
        "hpa_rna": None,
        "geo_expr": None,
        "fusions": ["gene1_name", "gene1_ens_id", "gene2_name", "gene2_ens_id"],
        "proteomics": None,
        "proteomics_gene_map": None,
    }

    for name, cols in checks.items():
        if name not in cleaned:
            continue
        try:
            df = cleaned[name]
            print(f"\n-- {name} --")
            print(df[cols].head() if cols else df.head())
        except Exception as e:
            log_error(logger, step="spot_check", error=e, dataset=name)
            print(f"  [SKIPPED] {name}: spot check failed - see log")


def missing_value_summary(cleaned: dict) -> pd.DataFrame:
    """
    Build and print a missing-value summary across all cleaned tables.

    Reports each table's shape and the percentage of missing cells over
    the whole frame, giving a post-cleaning completeness check comparable
    to the one produced in the data loading stage.

    Parameters
    ----------
    cleaned : dict
        Mapping of dataset key to cleaned :class:`pandas.DataFrame`, as
        returned by :func:`clean_data`.

    Returns
    -------
    pandas.DataFrame
        One row per dataset with columns ``Dataset``, ``Rows``, ``Cols``
        and ``Missing %``. Datasets whose summary failed are absent.

    Notes
    -----
    The table is also printed to stdout. Empty DataFrames are handled
    without raising — their missing percentage is reported as 0.
    """
    rows = []
    for name, df in cleaned.items():
        try:
            missing_pct = df.isnull().sum().sum() / df.size * 100 if df.size else 0
            rows.append({
                "Dataset": name,
                "Rows": f"{df.shape[0]:,}",
                "Cols": df.shape[1],
                "Missing %": f"{missing_pct:.1f}%",
            })
        except Exception as e:
            log_error(logger, step="missing_value_summary", error=e, dataset=name)
            print(f"  [SKIPPED] {name}: missing value summary failed - see log")

    summary = pd.DataFrame(rows)
    print("\nMissing value summary:")
    print(summary.to_string(index=False))
    return summary


def save_data(cleaned: dict) -> None:
    """
    Save the cleaned DataFrames to both Parquet and CSV.

    Writes Parquet to ``PARQUET_CLEAN`` and CSV to ``CSV_CLEAN``. Both
    formats are written so downstream stages can read the fast binary
    version while the CSVs stay inspectable by hand and by collaborators
    without a Parquet reader.

    Parameters
    ----------
    cleaned : dict
        Mapping of dataset key to cleaned :class:`pandas.DataFrame`, as
        returned by :func:`clean_data`.

    Returns
    -------
    None
        Files and log entries are written as side effects.

    Notes
    -----
    Both save helpers return the names of datasets that failed. The union
    of the two lists is logged as a single warning naming the dated log
    file, so a dataset that failed in only one format is still surfaced.
    A save failure does not raise — the run completes and reports.
    """
    print(f"\nSaving cleaned parquet files to {PARQUET_CLEAN} ...")
    parquet_failed = save_cleaned_parquets(cleaned, out_dir=PARQUET_CLEAN)

    print(f"\nSaving cleaned CSV files to {CSV_CLEAN} ...")
    csv_failed = save_cleaned_csvs(cleaned, out_dir=CSV_CLEAN)

    failed = sorted(set(parquet_failed) | set(csv_failed))
    if failed:
        log_file = LOG_DIR / f"pipeline_{datetime.now():%Y-%m-%d}.log"
        logger.warning(
            f"{len(failed)} dataset(s) had a save failure: "
            f"{', '.join(failed)} - see {log_file}"
        )
    print("Done.")


def main() -> None:
    """
    Run the full data cleaning stage end to end.

    Loads the raw Parquet tables, cleans them and drops non-human cell
    lines, prints shapes and spot checks, summarises missing values, then
    saves the cleaned tables as Parquet and CSV.

    Returns
    -------
    None
        Cleaned files and log entries are written as side effects;
        summaries are printed to stdout.

    Raises
    ------
    Exception
        Propagated from :func:`load_data` if the raw Parquet files cannot
        be read. All later failures are per-dataset: they are logged and
        skipped rather than raised.
    """
    tables = load_data()
    cleaned = clean_data(tables)
    inspect_shapes(cleaned)
    spot_check(cleaned)
    missing_value_summary(cleaned)
    save_data(cleaned)


if __name__ == "__main__":
    main()