"""
00_data_loading.py

Data Loading Pipeline (converted from 00_data_loading.ipynb)

Loads every raw data file found under BASE's subfolders (not a fixed
list of 14 -- see discover_raw_files), converts a couple of
tab-delimited text files to CSV, previews every dataset, builds a
summary table, and saves everything as Parquet files.

Each file is read from disk exactly once: load_all_datasets() loads
and previews it, then save_all_as_parquet() writes out the same
already-loaded DataFrame rather than re-reading it.

Run with:
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
    """Set pandas display options used throughout the notebook."""
    pd.set_option("display.max_columns", 20)
    pd.set_option("display.width", 120)
    pd.set_option("display.max_colwidth", 40)


def convert_text_files_to_csv() -> None:
    """
    Convert raw tab-delimited GEO text files to CSV.
    If a .txt source doesn't exist but a .csv version already does,
    conversion is skipped for it. Any other failure is logged with
    full file context before being raised.
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


# Curated read settings for the 14 known files (folder, filename, sep, skiprows, index_col, key)
# Anything NOT in this dict is a newly added file — its settings are auto-detected below.
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

# Extension -> default separator, used only for files NOT in KNOWN_FILES (i.e. newly added ones)
EXT_SEP_DEFAULTS = {".csv": ",", ".tsv": "\t", ".txt": "\t", ".gct": "\t"}



def slugify_filename(filename: str) -> str:
    """
    Derive a parquet key from a raw filename: strip a leading numeric
    prefix (e.g. '15_'), drop the extension,and collapse
    non-alphanumeric runs to underscores.
    e.g. '15_NewProteomicsPanel.csv' -> 'new_proteomics_panel'
    """
    stem = Path(filename).stem
    stem = re.sub(r"^\d+_", "", stem)
    stem = re.sub(r"[^0-9a-zA-Z]+", "_", stem).strip("_")
    return stem.lower()


def discover_raw_files(base: Path = BASE) -> list:
    """
    Walk every subfolder of `base` and return
    (path, sep, skiprows, index_col, key) for every data file found
    (.csv, .tsv, .txt, .gct), skipping the parquet output folder.

    Known files (the original 14) use their curated read settings.
    Any file not in KNOWN_FILES is treated as newly added: its
    separator is guessed from its extension and its key is derived
    from its filename, so it gets converted without any manual setup.
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
    Load every discovered raw file (see discover_raw_files) and
    preview it. Returns {key: (df, index_col)} — index_col is carried
    alongside each DataFrame so save_all_as_parquet() can write the
    index back out correctly for files (like depmap_expr) that were
    read with one, without needing to re-read the file from disk.

    A failure loading one file is logged with full context and that
    file is skipped — loading continues for the rest.
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
    """Build the missing-value / shape summary table across all loaded datasets."""
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
    Save already-loaded DataFrames (from load_all_datasets) to
    parquet at PARQUET_RAW/<key>.parquet. Per-dataset failures are
    logged with full context and the loop continues — one bad
    dataset doesn't stop the rest from saving.
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
    configure_pandas_display()
    convert_text_files_to_csv()
    loaded = load_all_datasets()
    build_summary_table(loaded)
    save_all_as_parquet(loaded)


if __name__ == "__main__":
    main()