"""
export_data.py
Export processed/cleaned data to parquet and generate summary reports.
"""

import os
import pandas as pd


_ROOT         = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
PARQUET_CLEAN = os.path.join(_ROOT, "data", "parquet", "data_clean")
PROCESSED_DIR = os.path.join(_ROOT, "outputs", "processed")
REPORT_PATH   = os.path.join(_ROOT, "outputs", "clean_data_report.xlsx")

CLEAN_FILE_MAP = {
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


def save_cleaned_parquets(cleaned: dict, out_dir: str = PARQUET_CLEAN) -> None:
    """Save a dict of cleaned DataFrames to parquet files."""
    os.makedirs(out_dir, exist_ok=True)
    for name, fname in CLEAN_FILE_MAP.items():
        if name not in cleaned:
            continue
        path = os.path.join(out_dir, fname)
        index = name == "depmap_expr"
        cleaned[name].to_parquet(path, index=index)
        size_mb = os.path.getsize(path) / 1_048_576
        print(f"Saved {name} -> {fname} ({size_mb:.1f} MB)")


def build_table_summary(df: pd.DataFrame, n_samples: int = 3) -> pd.DataFrame:
    """Return a per-column quality summary for a DataFrame."""
    rows = []
    n_rows = len(df)
    for col in df.columns:
        s = df[col]
        n_missing = s.isna().sum()
        n_present = n_rows - n_missing
        n_unique = s.nunique(dropna=True)
        n_duplicates = n_present - n_unique
        missing_pct = round((n_missing / n_rows) * 100, 2) if n_rows else 0
        is_pk = (n_missing == 0) and (n_duplicates == 0) and (n_rows > 0)
        sample = s.dropna().astype(str).head(n_samples).tolist()
        rows.append({
            "column": col,
            "dtype": str(s.dtype),
            "sample_values": ", ".join(sample),
            "missing_%": missing_pct,
            "unique_values": n_unique,
            "primary_key": is_pk,
            "duplicates": n_duplicates,
        })
    return pd.DataFrame(rows)


def export_clean_report(
    tables: dict,
    out_path: str = REPORT_PATH,
    n_samples: int = 3,
) -> None:
    """Write one Excel sheet per table with a per-column quality summary."""
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    with pd.ExcelWriter(out_path, engine="openpyxl") as writer:
        for name, df in tables.items():
            summary = build_table_summary(df, n_samples=n_samples)
            sheet_name = name[:31]
            header = pd.DataFrame({
                "info": [f"TABLE: {name}", f"shape: {df.shape[0]} rows x {df.shape[1]} cols"]
            })
            header.to_excel(writer, sheet_name=sheet_name, index=False, header=False, startrow=0)
            summary.to_excel(writer, sheet_name=sheet_name, index=False, startrow=3)
    print(f"Written {len(tables)} sheets to: {out_path}")


def export_processed_csvs(tables: dict, out_dir: str = PROCESSED_DIR) -> None:
    """Export a dict of DataFrames as CSV files to the processed directory."""
    os.makedirs(out_dir, exist_ok=True)
    for name, df in tables.items():
        path = os.path.join(out_dir, f"{name}.csv")
        df.to_csv(path, index=False)
        print(f"Exported {name} -> {path}")


if __name__ == "__main__":
    from data_utils import load_clean_parquets
    tables = load_clean_parquets()
    export_processed_csvs(tables)
    export_clean_report(tables)