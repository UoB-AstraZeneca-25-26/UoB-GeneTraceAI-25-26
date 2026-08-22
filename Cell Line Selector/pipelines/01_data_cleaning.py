"""
01_data_cleaning.py

Data Cleaning Pipeline (converted from 01_data_cleaning.ipynb)

Loads raw parquet files, applies per-dataset cleaning functions, purges
non-human cell lines, builds the proteomics gene<->uniprot lookup,
inspects the cleaned tables, summarizes missing values, and saves the
cleaned tables back out as both parquet and CSV.

Errors are logged to logs/pipeline_<date>.log with the dataset name,
source file, step, and full traceback, so an admin can find and fix
the specific bad file without re-running the whole pipeline.

Run with:
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
    """Load all raw parquet files into a dict of DataFrames."""
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
    Clean every dataset, then purge non-human cell lines.

    clean_cellosaurus() runs first because it populates NON_HUMAN_CVCLS as
    a side-effect. Every other table is cleaned by its own cleaner, then
    passed through filter_non_human_rows(), which matches on CVCL
    accession only -- never on cell-line name (see the NON_HUMAN_CVCLS
    comment in cleaning_functions.py for the measured reason).

    Tables with no CVCL column pass through unfiltered; their non-human
    rows are excluded downstream by the cell_line_roster model_id join.
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
    """Print row/column counts for every cleaned table."""
    print("\nDataset shapes:")
    for name, df in cleaned.items():
            print(f"  {name:40s}: {df.shape[0]:>10,} rows x {df.shape[1]:>6,} cols")


def spot_check(cleaned: dict) -> None:
    """
    Print sample rows from key datasets as a sanity check on the cleaning.
    A failure on one dataset's spot check (e.g. an expected column is
    missing) is logged and skipped rather than crashing the run.
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
    """Build and print a missing-value % summary across all cleaned tables."""
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
    """Save cleaned DataFrames to both parquet (PARQUET_CLEAN) and CSV (CSV_CLEAN)."""
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
    tables = load_data()
    cleaned = clean_data(tables)
    inspect_shapes(cleaned)
    spot_check(cleaned)
    missing_value_summary(cleaned)
    save_data(cleaned)


if __name__ == "__main__":
    main()