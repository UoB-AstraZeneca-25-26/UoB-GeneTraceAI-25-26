"""
data_utils.py

Shared paths and utilities used across the data pipeline scripts.

Every stage resolves its input and output locations from the constants
here rather than building paths of its own, so the directory layout is
defined in exactly one place and moving a directory is a one-line change.

Directory flow through the pipeline
-----------------------------------
    BASE            raw source files                (input to 00)
    PARQUET_RAW     raw, converted to parquet       (00 -> 01)
    PARQUET_CLEAN   cleaned tables + rosters        (01, 02 -> 03, 04)
    CSV_CLEAN       same, human-readable            (01, 02; not read back)
    DUCKDB_PATH     queryable database              (03)
    EDA_DIR         reports and figures             (04)

Note that PARQUET_RAW sits *inside* BASE, so anything walking BASE for
source files has to skip it — see ``discover_raw_files`` in
``00_data_loading.py``.

No side effects at import: the constants are computed but no directory is
created here. Each stage creates what it writes to.
"""

from pathlib import Path
import pandas as pd

# --- Project-relative data directories ---
PROJECT_ROOT = Path(__file__).resolve().parents[2]  # .../project

BASE = PROJECT_ROOT / "data" / "raw data"      # raw source files (csv/tsv/txt/gct), in subfolders:
#   gene expression/, gene properties/,
#   nomenclature/, non gene expression/
PARQUET_RAW = BASE / "parquet"     # raw datasets saved as parquet (output of 00_data_loading)
PARQUET_CLEAN = PROJECT_ROOT / "data" / "clean data"  # cleaned datasets saved as parquet (output of 01_data_cleaning)
CSV_CLEAN = PROJECT_ROOT / "data" / "clean_csv"  # cleaned datasets saved as CSV
DUCKDB_PATH = PROJECT_ROOT / "db" / "celllineselector.duckdb"
EDA_DIR = PROJECT_ROOT / "reports" / "eda"


def preview(df: pd.DataFrame, title: str = None, n: int = 5) -> None:
    """Print a quick shape + head preview of a DataFrame.

    Used during loading and cleaning to eyeball each dataset as it goes
    past — the shape catches a file that parsed into one column, the head
    catches a header row read as data.

    Parameters
    ----------
    df : pandas.DataFrame
        Frame to preview.
    title : str, optional
        Heading, printed inside a rule so datasets are easy to find when
        scrolling a long run. Omitted entirely when None.
    n : int, optional
        Rows to show. Default 5.

    Returns
    -------
    None
        Output goes to stdout.

    Notes
    -----
    How much of the head is shown depends on the pandas display options
    the caller has set — see ``configure_pandas_display`` in
    ``00_data_loading.py``.
    """
    if title:
        print(f"\n{'=' * 80}\n{title}\n{'=' * 80}")
    print(f"Shape: {df.shape[0]:,} rows x {df.shape[1]:,} cols")
    print(df.head(n))


def load_raw_parquets(raw_dir: Path = PARQUET_RAW) -> dict:
    """
    Load every .parquet file in `raw_dir` into a dict of DataFrames,
    keyed by filename (without extension).

    e.g. data/raw data/parquet/hpa_rna.parquet -> tables["hpa_rna"]

    The stem-as-key convention is what carries a dataset's identity
    through the pipeline: the same key names the table in ``CLEANERS``,
    the cleaned output file, and eventually the DuckDB table.

    Parameters
    ----------
    raw_dir : Path or str, optional
        Directory of raw Parquet files. Defaults to ``PARQUET_RAW``.

    Returns
    -------
    dict
        Mapping of filename stem to :class:`pandas.DataFrame`, built in
        sorted filename order.

    Raises
    ------
    FileNotFoundError
        If the directory does not exist, or exists but holds no Parquet
        files. Both are fatal: they mean ``00_data_loading.py`` has not
        run, so there is nothing to clean rather than something that
        failed to load.

    Notes
    -----
    Every file is read fully into memory, so peak usage is the whole raw
    dataset at once — the wide expression matrices dominate this.
    """
    raw_dir = Path(raw_dir)
    if not raw_dir.exists():
        raise FileNotFoundError(f"Raw parquet directory not found: {raw_dir}")

    tables = {}
    for file in sorted(raw_dir.glob("*.parquet")):
        tables[file.stem] = pd.read_parquet(file)

    if not tables:
        raise FileNotFoundError(f"No .parquet files found in {raw_dir}")

    return tables