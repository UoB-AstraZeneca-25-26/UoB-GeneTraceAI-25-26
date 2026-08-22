"""
data_utils.py

Shared paths and utilities used across the data pipeline scripts.
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
    """Print a quick shape + head preview of a DataFrame."""
    if title:
        print(f"\n{'=' * 80}\n{title}\n{'=' * 80}")
    print(f"Shape: {df.shape[0]:,} rows x {df.shape[1]:,} cols")
    print(df.head(n))


def load_raw_parquets(raw_dir: Path = PARQUET_RAW) -> dict:
    """
    Load every .parquet file in `raw_dir` into a dict of DataFrames,
    keyed by filename (without extension).

    e.g. data/raw data/parquet/hpa_rna.parquet -> tables["hpa_rna"]
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