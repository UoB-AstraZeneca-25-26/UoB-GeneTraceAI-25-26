"""
data_utils.py
Reusable functions for loading, previewing, and saving datasets.
"""

import os
import pandas as pd

BASE          = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', 'data'))
REF           = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', 'reference'))
PARQUET_RAW   = os.path.join(BASE, "parquet/raw_data")
PARQUET_CLEAN = os.path.join(BASE, "parquet/data_clean")


def preview(df: pd.DataFrame, name: str) -> None:
    """Print shape, dtypes, missing % and first 5 rows of a DataFrame."""
    print(f"{name}")
    print(f"Shape : {df.shape[0]:,} rows x {df.shape[1]:,} cols")
    print(f"Columns ({len(df.columns)}): {list(df.columns[:10])}" + (" ..." if len(df.columns) > 10 else ""))
    print(f"Dtypes :\n{df.dtypes.value_counts().to_string()}")
    null_pct = df.isnull().sum().sum() / df.size * 100
    print(f"Missing: {null_pct:.1f}% of all values")
    print("\n Top 5 rows")
    try:
        from IPython.display import display
        display(df.head(5))
    except Exception:
        print(df.head(5))


def load_raw_parquets(parquet_dir: str = PARQUET_RAW) -> dict:
    """Load all 14 raw parquet files and return a dict keyed by variable name."""
    files = {
        "hpa_rna":         "1_4_hpa_rna_celline.parquet",
        "depmap_expr":     "2_DepMap_OmicsExpressionAllGenesTPMLogp1Profile.parquet",
        "geo_expr":        "3_GEOexpression.parquet",
        "proteomics":      "4_Harmonized_MS_CCLE_Gygi_subsetted.parquet",
        "fusions":         "5_OmicsFusionFilteredSupplementary.parquet",
        "mutations":       "6_OmicsSomaticMutationsProfile.parquet",
        "cellosaurus":     "7_cellosaurus.parquet",
        "depmap_profiles": "8_DepMap_OmicsProfiles.parquet",
        "sample_info":     "9_DepMap_sample_info.parquet",
        "geo_info":        "10_GEOInfo.parquet",
        "hpa_desc":        "11_hpa_rna_celline_description.parquet",
        "metabolomics":    "12_CCLE_metabolomics_20190502.parquet",
        "mirna":           "13_CCLE_miRNA_20181103.parquet",
        "signatures":      "14_OmicsGlobalSignatures.parquet",
    }
    tables = {}
    for name, fname in files.items():
        path = os.path.join(parquet_dir, fname)
        tables[name] = pd.read_parquet(path)
        print(f"Loaded {name}: {tables[name].shape}")
    return tables


def load_clean_parquets(parquet_dir: str = PARQUET_CLEAN) -> dict:
    """Load all cleaned parquet files and return a dict keyed by variable name."""
    files = {
        "hpa_rna":         "hpa_rna_clean.parquet",
        "depmap_expr":     "depmap_expr_clean.parquet",
        "geo_expr":        "geo_expr_clean.parquet",
        "proteomics":      "proteomics_clean.parquet",
        "protein_map":     "protein_map.parquet",
        "fusions":         "fusions_clean.parquet",
        "mutations":       "mutations_clean.parquet",
        "cellosaurus":     "cellosaurus_clean.parquet",
        "depmap_profiles": "depmap_profiles_clean.parquet",
        "sample_info":     "sample_info_clean.parquet",
        "geo_info":        "geo_info_clean.parquet",
        "hpa_desc":        "hpa_desc_clean.parquet",
        "metabolomics":    "metabolomics_clean.parquet",
        "mirna":           "mirna_clean.parquet",
        "signatures":      "signatures_clean.parquet",
    }
    tables = {}
    for name, fname in files.items():
        path = os.path.join(parquet_dir, fname)
        tables[name] = pd.read_parquet(path)
        print(f"Loaded {name}: {tables[name].shape}")
    return tables


def save_dataframes_to_csv(tables: dict, out_dir: str) -> None:
    """Save a dict of DataFrames to CSV files in out_dir."""
    os.makedirs(out_dir, exist_ok=True)
    for name, df in tables.items():
        path = os.path.join(out_dir, f"{name}.csv")
        df.to_csv(path, index=False)
        print(f"Saved {name} -> {path}")