"""
audit_non_human.py

Post-hoc audit: verify that no non-human cell lines survived the
cleaning pipeline.

This runs INDEPENDENTLY of clean_data(). It rebuilds the non-human
blocklist from the RAW cellosaurus file, then scans every CLEANED
parquet for any surviving match. If the filter worked, every table
reports CLEAN. Anything else is a leak and is reported with the
offending column, count, and example values.

The audit deliberately re-derives the blocklist from raw rather than
importing NON_HUMAN_NAMES from a previous run — importing the same
set that did the filtering would make the check circular and it
would pass even if the filter never ran.

Run with:
    python pipelines/audit_non_human.py
    python pipelines/audit_non_human.py --csv          # audit CSV outputs instead
    python pipelines/audit_non_human.py --save-report  # write leaked rows to disk
"""

import argparse
import sys
from pathlib import Path

project_root = Path(__file__).resolve().parent
for candidate in [project_root, *project_root.parents]:
    if (candidate / "src").exists():
        project_root = candidate
        break
sys.path.insert(0, str(project_root))

import pandas as pd

from src.scripts.data_utils import PARQUET_RAW, PARQUET_CLEAN, CSV_CLEAN
from src.scripts.cleaning_functions import base_clean, rename_by_prefix


# Columns that can carry a cell-line identity. Kept in sync with
# filter_non_human_rows() in cleaning_functions.py — if you add a
# column there, add it here too or the audit will under-report.
NAME_COLS = [
    "cell_line_name",
    "cell_line",
    "cellline",
    "stripped_cell_line_name",
    "ccle_name",
    "ccle_id",
    "model_name",
]
CVCL_COLS = ["cvcl_id", "cellosaurus_id", "rrid"]

HUMAN_RE = r"9606|homo\s*sapiens|\bhuman\b"


# ============================================================
# Blocklist rebuild (from RAW cellosaurus)
# ============================================================

def build_blocklist_from_raw() -> tuple:
    """
    Rebuild the non-human name + cvcl blocklist from the RAW cellosaurus
    parquet. Mirrors clean_cellosaurus()'s logic exactly, but reads raw
    so the audit is independent of whatever the pipeline actually did.

    Returns (names, cvcls, n_non_human_rows) — both sets lowercased.
    """
    raw_path = PARQUET_RAW / "cellosaurus.parquet"
    if not raw_path.exists():
        raise FileNotFoundError(
            f"Raw cellosaurus not found at {raw_path} — cannot build blocklist."
        )

    df = pd.read_parquet(raw_path)
    df = base_clean(df)
    df = rename_by_prefix(df, {
        "identifier": "cell_line_name",
        "accession": "cvcl_id",
    })

    species_col = next(
        (c for c in df.columns if c.startswith("species_of_origin")), None
    )
    if species_col is None:
        raise KeyError(
            f"No 'species_of_origin' column in raw cellosaurus. "
            f"Columns: {list(df.columns)}"
        )

    is_human = df[species_col].fillna("").astype(str).str.contains(
        HUMAN_RE, regex=True, na=False
    )
    non_human = df[~is_human]

    names = set(
        non_human["cell_line_name"].dropna().astype(str).str.strip().str.lower()
    )

    syn_col = next((c for c in df.columns if c.startswith("synonym")), None)
    if syn_col is not None:
        for cell in non_human[syn_col].dropna().astype(str):
            for part in cell.split(";"):
                part = part.strip().lower()
                if part:
                    names.add(part)
    names -= {"", "nan"}

    cvcls = set(
        non_human["cvcl_id"].dropna().astype(str).str.strip().str.lower()
    ) - {"", "nan"}

    return names, cvcls, len(non_human)


# ============================================================
# Per-table audit
# ============================================================

def audit_table(name: str, df: pd.DataFrame, names: set, cvcls: set) -> dict:
    """
    Scan one cleaned table for surviving non-human rows.

    Returns a dict with:
      status    — CLEAN / LEAK / NO_ID_COL
      n_rows    — total rows in the table
      n_leaked  — rows matching the blocklist
      hits      — {column: count} for columns with matches
      examples  — up to 5 offending values, in original capitalization
    """
    name_cols = [c for c in NAME_COLS if c in df.columns]
    cvcl_cols = [c for c in CVCL_COLS if c in df.columns]

    if not name_cols and not cvcl_cols:
        return {
            "status": "NO_ID_COL",
            "n_rows": len(df),
            "n_leaked": 0,
            "hits": {},
            "examples": [],
            "checked": [],
        }

    mask = pd.Series(False, index=df.index)
    hits = {}

    for col in name_cols:
        m = df[col].fillna("").astype(str).str.strip().str.lower().isin(names)
        if m.any():
            hits[col] = int(m.sum())
        mask |= m

    for col in cvcl_cols:
        m = df[col].fillna("").astype(str).str.strip().str.lower().isin(cvcls)
        if m.any():
            hits[col] = int(m.sum())
        mask |= m

    n_leaked = int(mask.sum())

    examples = []
    if n_leaked:
        seen = set()
        for col in list(hits.keys()):
            for v in df.loc[mask, col].dropna().astype(str).str.strip().head(5):
                if v not in seen:
                    seen.add(v)
                    examples.append(f"{col}={v}")
                if len(examples) >= 5:
                    break
            if len(examples) >= 5:
                break

    return {
        "status": "LEAK" if n_leaked else "CLEAN",
        "n_rows": len(df),
        "n_leaked": n_leaked,
        "hits": hits,
        "examples": examples,
        "checked": name_cols + cvcl_cols,
        "mask": mask if n_leaked else None,
    }


# ============================================================
# Main
# ============================================================

def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--csv", action="store_true",
        help="Audit the cleaned CSV outputs instead of parquet.",
    )
    parser.add_argument(
        "--save-report", action="store_true",
        help="Write the leaked rows of each failing table to logs/audit_leaks/.",
    )
    args = parser.parse_args()

    src_dir = CSV_CLEAN if args.csv else PARQUET_CLEAN
    suffix = "*.csv" if args.csv else "*.parquet"
    reader = pd.read_csv if args.csv else pd.read_parquet

    print("=" * 78)
    print("NON-HUMAN CELL LINE AUDIT")
    print("=" * 78)

    # ── Rebuild blocklist from RAW ───────────────────────────────────────
    print("\nRebuilding blocklist from raw cellosaurus...")
    try:
        names, cvcls, n_non_human = build_blocklist_from_raw()
    except (FileNotFoundError, KeyError) as e:
        print(f"  [FATAL] {e}")
        return 2

    print(f"  {n_non_human:,} non-human cellosaurus rows")
    print(f"  {len(names):,} blocked names/synonyms")
    print(f"  {len(cvcls):,} blocked cvcl_ids")

    if not names and not cvcls:
        print("  [FATAL] Blocklist is empty — audit cannot verify anything.")
        return 2

    # ── Scan cleaned tables ──────────────────────────────────────────────
    files = sorted(src_dir.glob(suffix))
    if not files:
        print(f"\n  [FATAL] No {suffix} files found in {src_dir}")
        return 2

    print(f"\nScanning {len(files)} cleaned tables in {src_dir}\n")
    print(f"  {'Table':<42} {'Status':<10} {'Rows':>10} {'Leaked':>8}")
    print("  " + "-" * 74)

    results = {}
    for path in files:
        table = path.stem
        try:
            df = reader(path)
        except Exception as e:
            print(f"  {table:<42} {'READ_FAIL':<10} {'—':>10} {'—':>8}")
            print(f"      {type(e).__name__}: {e}")
            continue

        r = audit_table(table, df, names, cvcls)
        results[table] = (r, df)

        leaked_str = f"{r['n_leaked']:,}" if r["status"] == "LEAK" else "—"
        marker = "  ***" if r["status"] == "LEAK" else ""
        print(
            f"  {table:<42} {r['status']:<10} {r['n_rows']:>10,} "
            f"{leaked_str:>8}{marker}"
        )

        if r["status"] == "LEAK":
            print(f"      hits: {r['hits']}")
            print(f"      e.g. {r['examples']}")

    print("  " + "-" * 74)

    # ── Summary ──────────────────────────────────────────────────────────
    leaks = {k: v for k, (v, _) in results.items() if v["status"] == "LEAK"}
    no_id = [k for k, (v, _) in results.items() if v["status"] == "NO_ID_COL"]
    clean = [k for k, (v, _) in results.items() if v["status"] == "CLEAN"]

    print(f"\n  CLEAN      : {len(clean)} tables")
    print(f"  NO_ID_COL  : {len(no_id)} tables  {no_id if no_id else ''}")
    print(f"  LEAK       : {len(leaks)} tables")

    if no_id:
        print(
            "\n  Note: NO_ID_COL tables have no cell-line name or CVCL column,"
            "\n  so they cannot be audited directly. Their cell-line identity"
            "\n  is resolved downstream via model_id/profile_id joins — verify"
            "\n  those separately after harmonisation."
        )

    if leaks:
        total_leaked = sum(v["n_leaked"] for v in leaks.values())
        print(f"\n  FAILED: {total_leaked:,} non-human rows survived cleaning.")

        if args.save_report:
            out_dir = project_root / "logs" / "audit_leaks"
            out_dir.mkdir(parents=True, exist_ok=True)
            for table in leaks:
                r, df = results[table]
                leaked_rows = df[r["mask"]]
                out_path = out_dir / f"{table}_leaked.csv"
                leaked_rows.to_csv(out_path, index=False)
                print(f"    wrote {len(leaked_rows):,} rows -> {out_path}")

        return 1

    print("\n  PASSED: no non-human cell lines found in any auditable table.")
    return 0


if __name__ == "__main__":
    sys.exit(main())