"""
00_data_loading.py

Data loading stage of the pipeline (converted from ``00_data_loading.ipynb``).

Walks every subfolder of ``BASE``, discovers each raw data file it finds
(rather than working from a hard-coded list of 14 — see
:func:`discover_raw_files`), converts two tab-delimited GEO text files to
CSV, previews every dataset, builds a missing-value/shape summary table,
and writes each dataset out as Parquet under ``PARQUET_RAW``.

Design notes
------------
* Each file is read from disk exactly once. :func:`load_all_datasets`
  loads and previews it, then :func:`save_all_as_parquet` writes out the
  same already-loaded DataFrame rather than re-reading it.
* Failures are per-file, not fatal. A file that fails to load or save is
  logged with full context (path, separator, skiprows, index_col, key)
  and skipped, so one malformed dataset does not stop the rest.
* Newly added files need no manual configuration: files absent from
  :data:`KNOWN_FILES` get their separator guessed from their extension
  and their Parquet key derived from their filename.

Inputs
------
Raw ``.csv`` / ``.tsv`` / ``.txt`` / ``.gct`` files under ``BASE``.

Outputs
-------
One Parquet file per dataset at ``PARQUET_RAW/<key>.parquet``, plus a
dated log file under ``LOG_DIR``.

Run with
--------
    python pipelines/00_data_loading.py
"""

import os
import re
import sys
from datetime import datetime
from pathlib import Path

# --- Make sure the project root (the folder containing 'src') is importable ---
project_root = Path(__file__).resolve().parent
for candidate in [project_root, *project_root.parents]:
    if (candidate / "src").exists():
        project_root = candidate
        break
sys.path.insert(0, str(project_root))

import pandas as pd

from src.scripts.cleaning_functions import clean_procan_raw_tsv
from src.scripts.data_utils import preview, BASE, PARQUET_RAW
from src.scripts.logging_utils import get_logger, log_file_error, LOG_DIR

logger = get_logger("00_data_loading")


def configure_pandas_display() -> None:
    """
    Set the pandas display options used throughout this stage.

    Widens the console output so that :func:`preview` prints readable
    tables: up to 20 columns, a 120-character line width, and a
    40-character cap on individual cell values.

    Returns
    -------
    None
        The options are set globally as a side effect.

    Notes
    -----
    This affects display only — no data is altered.
    """
    pd.set_option("display.max_columns", 20)
    pd.set_option("display.width", 120)
    pd.set_option("display.max_colwidth", 40)


def convert_text_files_to_csv() -> None:
    """
    Convert the raw tab-delimited GEO text files to CSV.

    Two GEO files (``3_GEOexpression.txt`` and ``10_GEOInfo.txt``) ship as
    tab-delimited ``.txt``; both are rewritten as ``.csv`` alongside the
    source so the rest of the pipeline sees a consistent format.

    Conversion is idempotent-friendly: if the ``.txt`` source is missing
    but the ``.csv`` output already exists, that pair is skipped and the
    run continues.

    Returns
    -------
    None
        Files are written to disk as a side effect.

    Raises
    ------
    FileNotFoundError
        If neither the ``.txt`` source nor the ``.csv`` output exists for
        a configured pair. Logged with file context before being raised.
    Exception
        Any read/write failure is logged with the file path, step name,
        and separator, then re-raised unchanged.

    See Also
    --------
    discover_raw_files : Picks up the converted CSVs on the next walk.
    """
    print(f"BASE = {BASE}")

    conversions = [
        ("gene expression/3_GEOexpression.txt", "gene expression/3_GEOexpression.csv"),
        ("nomenclature/10_GEOInfo.txt", "nomenclature/10_GEOInfo.csv"),
    ]

    for src_rel, dst_rel in conversions:
        src = BASE / src_rel
        dst = BASE / dst_rel

        if not src.exists():
            if dst.exists():
                print(f"  {src.name} not found, but {dst.name} already exists — skipping conversion.")
                continue
            error = FileNotFoundError(f"Neither raw .txt nor converted .csv found for {src.name}")
            log_file_error(logger, file_path=src, step="convert_text_files_to_csv", error=error)
            raise error

        try:
            df = pd.read_csv(src, sep="\t", low_memory=False)
            df.to_csv(dst, index=False)
            print(f"  converted {src.name} -> {dst.name}")
        except Exception as e:
            log_file_error(logger, file_path=src, step="convert_text_files_to_csv", error=e, sep=repr("\t"))
            raise

    print("Text files converted to CSV (where needed).")


#: Curated read settings for the 14 known raw files, keyed by filename.
#:
#: Each value is a dict of ``sep``, ``skiprows``, ``index_col`` and ``key``,
#: where ``key`` becomes the Parquet filename stem. These settings were
#: established by inspecting each file (e.g. the GCT miRNA matrix carries
#: two header lines; the DepMap expression matrix has cell lines in the
#: index). Any file *not* listed here is treated as newly added and has
#: its settings auto-detected in :func:`discover_raw_files`.
KNOWN_FILES = {
    "1_4_hpa_rna_celline.tsv":                             dict(sep="\t", skiprows=0, index_col=None, key="hpa_rna"),
    "2_DepMap_OmicsExpressionAllGenesTPMLogp1Profile.csv": dict(sep=",",  skiprows=0, index_col=0,    key="depmap_expr"),
    "3_GEOexpression.csv":                                 dict(sep="\t", skiprows=0, index_col=None, key="geo_expr"),
    "4_Harmonized_MS_CCLE_Gygi_subsetted.csv":              dict(sep=",",  skiprows=0, index_col=None, key="proteomics"),
    "5_OmicsFusionFilteredSupplementary.csv":               dict(sep=",",  skiprows=0, index_col=None, key="fusions"),
    "6_OmicsSomaticMutationsProfile.csv":                   dict(sep=",",  skiprows=0, index_col=None, key="mutations"),
    "7_cellosaurus.csv":                                    dict(sep=",",  skiprows=0, index_col=None, key="cellosaurus"),
    "8_DepMap_OmicsProfiles.csv":                           dict(sep=",",  skiprows=0, index_col=None, key="depmap_profiles"),
    "9_DepMap_sample_info.csv":                             dict(sep=",",  skiprows=0, index_col=None, key="sample_info"),
    "10_GEOInfo.txt":                                       dict(sep="\t", skiprows=0, index_col=None, key="geo_info"),
    "11_hpa_rna_celline_description.tsv":                   dict(sep="\t", skiprows=0, index_col=None, key="hpa_desc"),
    "12_CCLE_metabolomics_20190502.csv":                    dict(sep=",",  skiprows=0, index_col=None, key="metabolomics"),
    "13_CCLE_miRNA_20181103.gct":                           dict(sep="\t", skiprows=2, index_col=None, key="mirna"),
    "14_OmicsGlobalSignatures.csv":                         dict(sep=",",  skiprows=0, index_col=None, key="signatures"),
    "10_GEOInfo.csv":                                       dict(sep=",",  skiprows=0, index_col=None, key="geo_info"),
}

#: Fallback separator per file extension.
#:
#: Used only for files absent from :data:`KNOWN_FILES` — i.e. newly added
#: raw files — so they can be ingested without manual configuration.
EXT_SEP_DEFAULTS = {".csv": ",", ".tsv": "\t", ".txt": "\t", ".gct": "\t"}


def slugify_filename(filename: str) -> str:
    """
    Derive a Parquet key from a raw filename.

    Strips a leading numeric prefix (e.g. ``15_``), drops the extension,
    collapses runs of non-alphanumeric characters to single underscores,
    trims leading/trailing underscores, and lowercases the result.

    Parameters
    ----------
    filename : str
        Raw filename or path; only the stem is used.

    Returns
    -------
    str
        A lowercase, underscore-separated key safe to use as a Parquet
        filename stem.

    Examples
    --------
    >>> slugify_filename("15_NewProteomicsPanel.csv")
    'newproteomicspanel'
    >>> slugify_filename("CCLE_miRNA_20181103.gct")
    'ccle_mirna_20181103'
    """
    stem = Path(filename).stem
    stem = re.sub(r"^\d+_", "", stem)
    stem = re.sub(r"[^0-9a-zA-Z]+", "_", stem).strip("_")
    return stem.lower()


def discover_raw_files(base: Path = BASE) -> list:
    """
    Walk every subfolder of ``base`` and return the read settings for each raw file.

    Recurses through ``base``, keeping files with a ``.csv``, ``.tsv``,
    ``.txt`` or ``.gct`` extension and skipping anything under a
    ``parquet`` folder (which holds this stage's own output).

    Files listed in :data:`KNOWN_FILES` use their curated read settings.
    Anything else is treated as newly added: its separator is guessed
    from its extension via :data:`EXT_SEP_DEFAULTS` and its key is derived
    by :func:`slugify_filename`, so it is ingested without manual setup.

    Parameters
    ----------
    base : Path, optional
        Root directory to walk. Defaults to ``BASE``.

    Returns
    -------
    list of tuple
        One ``(path, sep, skiprows, index_col, key)`` tuple per discovered
        file, sorted by path, where:

        * ``path`` (:class:`~pathlib.Path`) — absolute file path;
        * ``sep`` (str) — field delimiter passed to :func:`pandas.read_csv`;
        * ``skiprows`` (int) — leading lines to skip (0 for most files);
        * ``index_col`` (int or None) — column to use as the index;
        * ``key`` (str) — Parquet filename stem for this dataset.

    Notes
    -----
    Newly detected files are announced on stdout so an unexpected addition
    to the raw data folder is visible in the run output.
    """
    discovered = []
    for path in sorted(base.rglob("*")):
        if not path.is_file() or path.suffix.lower() not in EXT_SEP_DEFAULTS:
            continue
        if "parquet" in path.parts:
            continue

        if path.name in KNOWN_FILES:
            cfg = KNOWN_FILES[path.name]
            discovered.append((path, cfg["sep"], cfg["skiprows"], cfg["index_col"], cfg["key"]))
        else:
            sep = EXT_SEP_DEFAULTS[path.suffix.lower()]
            key = slugify_filename(path.name)
            print(f"  [new file detected] {path.name} -> key '{key}' (sep={sep!r})")
            discovered.append((path, sep, 0, None, key))

    return discovered


def load_all_datasets() -> dict:
    """
    Load and preview every raw file returned by :func:`discover_raw_files`.

    Each file is read once with its configured separator, skiprows and
    index column, then printed via :func:`preview`. The ``index_col`` is
    carried alongside the DataFrame so :func:`save_all_as_parquet` can
    write the index back out correctly for files that were read with one
    (e.g. ``depmap_expr``) without re-reading from disk.

    Returns
    -------
    dict
        Mapping of ``key`` (str) to ``(df, index_col)``, where ``df`` is a
        :class:`pandas.DataFrame` and ``index_col`` is the int or None
        used when reading it. Files that failed to load are absent.

    Notes
    -----
    * Failures are non-fatal: a file that raises during read is logged
      with full context, reported on stdout as ``[SKIPPED]``, and loading
      continues with the remaining files.
    * ``Protein_matrix_averaged_20250211.tsv`` is read with
      ``header=None`` because it carries no header row.
    * Retained from the original notebook: the GEO expression file is
      sometimes comma-separated despite its extension, so if the
      tab-separated read yields a single column it is re-read with pandas'
      default separator.
    """
    loaded = {}
    for path, sep, skiprows, index_col, key in discover_raw_files():
        try:
            read_kwargs = dict(
                sep=sep,
                skiprows=skiprows if skiprows > 0 else None,
                index_col=index_col,
                low_memory=False,
            )
            if path.name == "Protein_matrix_averaged_20250211.tsv":
                read_kwargs["header"] = None

            df = pd.read_csv(path, **read_kwargs)

            # Retained from the original notebook: GEO expression sometimes
            # isn't actually tab-separated despite the file extension.
            if key == "geo_expr" and df.shape[1] == 1:
                df = pd.read_csv(path, low_memory=False)
                index_col = None
        except Exception as e:
            log_file_error(
                logger, file_path=path, step="load_all_datasets", error=e,
                sep=repr(sep), skiprows=skiprows, index_col=index_col, key=key,
            )
            print(f"  [SKIPPED] {path.name}: load failed — see log")
            continue

        preview(df, f"{key} — {path.name}")
        loaded[key] = (df, index_col)

    return loaded


def build_summary_table(loaded: dict) -> pd.DataFrame:
    """
    Build a missing-value and shape summary across all loaded datasets.

    For each dataset, reports its row and column counts, the percentage of
    missing cells overall, the column with the most missing values, and
    what percentage of that column is missing — a quick check on data
    completeness before any downstream cleaning.

    Parameters
    ----------
    loaded : dict
        Mapping of ``key`` to ``(df, index_col)``, as returned by
        :func:`load_all_datasets`. The ``index_col`` element is ignored.

    Returns
    -------
    pandas.DataFrame
        One row per dataset with columns ``Dataset``, ``Rows``, ``Cols``,
        ``Missing %``, ``Worst col`` and ``Worst col %``. Counts and
        percentages are pre-formatted as display strings.

    Notes
    -----
    The table is also printed to stdout. Empty DataFrames are handled
    without raising: their missing percentage is reported as 0 and the
    worst column as ``None``.
    """
    rows = []
    for key, (df, _) in loaded.items():
        missing_pct = df.isnull().sum().sum() / df.size * 100 if df.size else 0
        worst_col = df.isnull().sum().idxmax() if df.size else None
        worst_pct = df.isnull().sum().max() / len(df) * 100 if len(df) else 0
        rows.append({
            "Dataset": key,
            "Rows": f"{df.shape[0]:,}",
            "Cols": f"{df.shape[1]:,}",
            "Missing %": f"{missing_pct:.1f}%",
            "Worst col": worst_col,
            "Worst col %": f"{worst_pct:.1f}%",
        })
    summary = pd.DataFrame(rows)
    print("\nDataset summary:")
    print(summary.to_string(index=False))
    return summary


def save_all_as_parquet(loaded: dict) -> pd.DataFrame:
    """
    Write the already-loaded DataFrames out as Parquet files.

    Each dataset is saved to ``PARQUET_RAW/<key>.parquet``, creating the
    output directory if needed. The index is written only for datasets
    that were read with an ``index_col``, preserving row labels (e.g. cell
    line identifiers) without adding a spurious integer column elsewhere.

    Parameters
    ----------
    loaded : dict
        Mapping of ``key`` to ``(df, index_col)``, as returned by
        :func:`load_all_datasets`.

    Returns
    -------
    pandas.DataFrame
        One row per dataset with columns ``Dataset``, ``Saved as``,
        ``Shape`` and ``Parquet MB``. Failed saves carry ``—`` in the
        first three fields and ``ERROR — see log`` in the last.

    Notes
    -----
    Per-dataset failures are logged with the output path, dataset key and
    shape, then skipped — one bad dataset does not stop the rest from
    saving. A summary count is written to the log, and a warning naming
    the dated log file is emitted if any dataset failed.
    """
    PARQUET_RAW.mkdir(parents=True, exist_ok=True)

    results = []
    for key, (df, index_col) in loaded.items():
        out_path = PARQUET_RAW / f"{key}.parquet"
        try:
            df.to_parquet(out_path, index=(index_col is not None))
            size_mb = os.path.getsize(out_path) / 1_048_576
            results.append({
                "Dataset": key, "Saved as": out_path.name,
                "Shape": f"{df.shape[0]:,} x {df.shape[1]:,}", "Parquet MB": f"{size_mb:.1f}",
            })
        except Exception as e:
            log_file_error(
                logger, file_path=out_path, step="save_all_as_parquet", error=e,
                dataset=key, shape=str(df.shape),
            )
            results.append({"Dataset": key, "Saved as": "—", "Shape": "—", "Parquet MB": "ERROR — see log"})

    results_df = pd.DataFrame(results)
    print("\nParquet save results:")
    print(results_df.to_string(index=False))

    n_ok = (results_df["Saved as"] != "—").sum()
    logger.info(f"Saved {n_ok} / {len(results_df)} datasets to parquet.")
    if n_ok < len(results_df):
        log_file = LOG_DIR / f"pipeline_{datetime.now():%Y-%m-%d}.log"
        n_failed = len(results_df) - n_ok
        logger.warning(f"{n_failed} dataset(s) failed to save — see {log_file} for details")

    return results_df


def main() -> None:
    """
    Run the full data loading stage end to end.

    Sets pandas display options, converts the GEO text files to CSV,
    loads and previews every discovered raw file, prints the
    missing-value summary, and writes each dataset out as Parquet.

    Returns
    -------
    None
        Parquet files and log entries are written as side effects; the
        summary and save tables are printed to stdout.

    Raises
    ------
    FileNotFoundError
        Propagated from :func:`convert_text_files_to_csv` if a configured
        GEO file is missing in both ``.txt`` and ``.csv`` form. Load and
        save failures for individual datasets are logged and skipped
        rather than raised.
    """
    configure_pandas_display()
    convert_text_files_to_csv()
    loaded = load_all_datasets()
    build_summary_table(loaded)
    save_all_as_parquet(loaded)


if __name__ == "__main__":
    main()