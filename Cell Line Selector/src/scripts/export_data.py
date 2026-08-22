"""
export_data.py

Save cleaned DataFrames to parquet and CSV.
"""

from pathlib import Path

from src.scripts.logging_utils import get_logger, log_file_error

logger = get_logger("export_data")


def save_cleaned_parquets(cleaned: dict, out_dir: Path) -> list:
    """
    Write each DataFrame in `cleaned` to `out_dir/<name>.parquet`.
    A failure on one dataset is logged with full context and that
    dataset is skipped — saving continues for the rest.
    Returns a list of dataset names that failed.
    """
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    failed = []
    for name, df in cleaned.items():
        out_path = out_dir / f"{name}.parquet"
        try:
            df.to_parquet(out_path, index=False)
            print(f"  saved {name} -> {out_path}")
        except Exception as e:
            log_file_error(
                logger, file_path=out_path, step="save_cleaned_parquets",
                error=e, dataset=name, shape=str(df.shape),
            )
            print(f"  [FAILED] {name}: save failed — see log")
            failed.append(name)

    return failed


def save_cleaned_csvs(cleaned: dict, out_dir: Path) -> list:
    """
    Write each DataFrame in `cleaned` to `out_dir/<name>.csv`.
    Same per-dataset error handling as save_cleaned_parquets: a
    failure on one dataset is logged and skipped, the rest continue.
    Returns a list of dataset names that failed.
    """
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    failed = []
    for name, df in cleaned.items():
        out_path = out_dir / f"{name}.csv"
        try:
            df.to_csv(out_path, index=False)
            print(f"  saved {name} -> {out_path}")
        except Exception as e:
            log_file_error(
                logger, file_path=out_path, step="save_cleaned_csvs",
                error=e, dataset=name, shape=str(df.shape),
            )
            print(f"  [FAILED] {name}: CSV save failed — see log")
            failed.append(name)

    return failed