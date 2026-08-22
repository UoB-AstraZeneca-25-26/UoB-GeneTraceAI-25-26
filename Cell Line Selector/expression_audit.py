from pathlib import Path
import re

import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parent
CLEAN_DIR = PROJECT_ROOT / "data" / "clean data"
PATTERN = re.compile(r"^ENSG\d+", re.IGNORECASE)


def canonical_roster_ids() -> set[str]:
    roster_path = CLEAN_DIR / "gene_roster.parquet"
    if not roster_path.exists():
        raise FileNotFoundError(f"Missing roster file: {roster_path}")

    roster = pd.read_parquet(roster_path)
    gene_ids = roster["gene_id"].dropna().astype(str)
    return {str(v).upper() for v in gene_ids if PATTERN.match(str(v))}


def count_invalid_ensg_columns(columns: list[str], valid_ids: set[str]) -> tuple[int, int]:
    ensg_cols = [str(c).upper() for c in columns if isinstance(c, str) and PATTERN.match(c)]
    valid = sum(1 for c in ensg_cols if c in valid_ids)
    invalid = sum(1 for c in ensg_cols if c not in valid_ids)
    return valid, invalid


def count_invalid_gene_ids(values: pd.Series, valid_ids: set[str]) -> tuple[int, int]:
    clean = values.dropna().astype(str)
    ensg_values = [str(v).upper() for v in clean if PATTERN.match(str(v))]
    valid = sum(1 for v in ensg_values if v in valid_ids)
    invalid = sum(1 for v in ensg_values if v not in valid_ids)
    return valid, invalid


def main() -> None:
    valid_ids = canonical_roster_ids()
    print(f"Roster gene count: {len(valid_ids):,}")
    print(f"Clean data dir: {CLEAN_DIR}")
    print("-" * 80)

    for table_name in ["depmap_expr", "geo_expr", "hpa_rna", "gene_roster"]:
        path = CLEAN_DIR / f"{table_name}.parquet"
        if not path.exists():
            print(f"{table_name}: MISSING ({path})")
            continue

        df = pd.read_parquet(path)
        print(f"{table_name}: shape={df.shape}")

        if "gene_id" in df.columns:
            valid, invalid = count_invalid_gene_ids(df["gene_id"], valid_ids)
            print(f"  gene_id rows: {len(df):,}")
            print(f"  valid roster ENSG rows: {valid:,}")
            print(f"  invalid non-roster ENSG rows: {invalid:,}")
        else:
            valid, invalid = count_invalid_ensg_columns(list(df.columns), valid_ids)
            print(f"  gene columns: {len(df.columns):,}")
            print(f"  valid roster ENSG columns: {valid:,}")
            print(f"  invalid non-roster ENSG columns: {invalid:,}")

        print()

    print("Done.")


if __name__ == "__main__":
    main()
