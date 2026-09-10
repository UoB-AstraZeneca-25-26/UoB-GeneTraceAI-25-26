"""
export_data.py

Save cleaned DataFrames to parquet and CSV.

Two writers with identical semantics, one per format. Both stages that
persist cleaned data (``01_data_cleaning.py`` and
``02_data_harmonisation.py``) call both, so every cleaned table exists in
each format under its own directory.

Design notes
------------
* Both formats are written on purpose. Parquet is what downstream stages
  read — typed, compressed, and fast over the wide expression tables —
  while the CSVs stay inspectable by hand and by collaborators without a
  Parquet reader.
* Failures are per-dataset, never fatal. One dataset that cannot be
  written is logged with its name and shape, reported on stdout, and
  skipped; the rest still save. Each writer returns the names that
  failed so the caller can report them together.
* The index is never written. Cleaned tables carry their identifiers as
  columns, so a written index would appear downstream as a spurious
  ``Unnamed: 0``.
* Output directories are created on demand, so a first run needs no
  setup.

Both functions are called for their side effects and their failure list;
neither modifies the frames it is given.
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

    This is the format downstream stages read from, so a dataset absent
    from the returned failure list is one they can expect to find.

    Parameters
    ----------
    cleaned : dict
        Mapping of dataset name to :class:`pandas.DataFrame`. The key
        becomes the filename stem, which is also how later stages address
        the table.
    out_dir : Path or str
        Destination directory, created with any missing parents if it
        does not exist.

    Returns
    -------
    list of str
        Names of datasets that failed to save. Empty on full success.

    Notes
    -----
    Written with ``index=False``. A table whose row labels matter must
    carry them as a column before reaching here, or they are lost.

    Existing files are overwritten without warning, so a rerun replaces
    the previous output rather than appending to it. Note that a dataset
    dropped between runs leaves its stale file behind — nothing here
    removes files for keys that are no longer present.
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

    These CSVs are the human-readable mirror of the Parquet output —
    nothing in the pipeline reads them back, so a CSV failure does not
    block downstream stages the way a Parquet failure would.

    Parameters
    ----------
    cleaned : dict
        Mapping of dataset name to :class:`pandas.DataFrame`. The key
        becomes the filename stem.
    out_dir : Path or str
        Destination directory, created with any missing parents if it
        does not exist.

    Returns
    -------
    list of str
        Names of datasets that failed to save. Empty on full success.

    Notes
    -----
    Written with ``index=False``, and existing files are overwritten, as
    with the Parquet writer.

    CSV loses the typing Parquet preserves, so a value read back from
    here will not necessarily match its Parquet counterpart in dtype.
    That is acceptable for inspection but not for re-entering the
    pipeline.

    The wide expression tables produce very large CSVs — this is usually
    the slowest step in a cleaning run, and the likeliest place to hit a
    disk-space failure.
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