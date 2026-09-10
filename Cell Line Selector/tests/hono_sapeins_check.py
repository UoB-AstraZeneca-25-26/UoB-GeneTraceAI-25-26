"""
audit_name_overlap.py

Check whether any cell-line NAME or SYNONYM is used by BOTH a human and a
non-human Cellosaurus entry.

Why this matters
----------------
filter_non_human_rows() drops rows from every downstream table whose
cell_line_name matches the non-human blocklist. That blocklist is built
from names + semicolon-separated synonyms of non-human Cellosaurus rows.

If a name is shared — e.g. 'L2' is both a human line and a rat line, or a
non-human line's synonym happens to be a human line's primary name — then
filtering by name will silently delete legitimate HUMAN rows from hpa_rna,
sample_info, etc. That's a false-positive purge, and it's invisible unless
you check for it explicitly.

This script does NOT modify anything. It reports:
  1. Primary-name collisions   (human cell_line_name == non-human cell_line_name)
  2. Synonym-mediated collisions (human name appears in a non-human row's synonyms,
     or vice versa)
  3. Short/generic blocklist entries (<= 3 chars) that are the highest-risk
     candidates for accidental collision
  4. Which downstream tables would actually lose human rows to each collision

Run with:
    python pipelines/audit_name_overlap.py
    python pipelines/audit_name_overlap.py --min-len 4   # simulate a length guard
    python pipelines/audit_name_overlap.py --save        # write CSVs to logs/audit_overlap/
    python pipelines/audit_name_overlap.py --check-tables  # scan cleaned tables for impact
"""

import argparse
import re
import sys
from pathlib import Path

project_root = Path(__file__).resolve().parent
for candidate in [project_root, *project_root.parents]:
    if (candidate / "src").exists():
        project_root = candidate
        break
sys.path.insert(0, str(project_root))

import pandas as pd

from src.scripts.data_utils import PARQUET_RAW, PARQUET_CLEAN
from src.scripts.cleaning_functions import base_clean, rename_by_prefix


# Must stay in sync with clean_cellosaurus(). Duplicated deliberately so the
# audit is independent of whatever the pipeline actually did.
HUMAN_ORIGIN_RE = re.compile(
    r"""
      \b9606\b
    | homo\s*sapiens
    | \bh\.\s*sapiens\b
    | \bhuman\b
    | \bhumans\b
    """,
    re.IGNORECASE | re.VERBOSE,
)

# Must stay in sync with filter_non_human_rows().
NAME_COLS = ["cell_line_name", "cell_line", "cellline",
             "stripped_cell_line_name", "ccle_name", "ccle_id", "model_name"]


# ============================================================
# Load + split cellosaurus
# ============================================================

def load_cellosaurus() -> tuple:
    """
    Read raw cellosaurus, apply base_clean + the standard renames, and split
    into (human_df, non_human_df, species_col, synonym_col).
    """
    path = PARQUET_RAW / "cellosaurus.parquet"
    if not path.exists():
        raise FileNotFoundError(f"Raw cellosaurus not found at {path}")

    df = pd.read_parquet(path)
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
            f"No 'species_of_origin' column. Columns: {list(df.columns)}"
        )

    syn_col = next((c for c in df.columns if c.startswith("synonym")), None)

    species = df[species_col].fillna("").astype(str)
    is_human = species.str.contains(HUMAN_ORIGIN_RE, na=False)

    return df[is_human].copy(), df[~is_human].copy(), species_col, syn_col


def name_map(sub: pd.DataFrame, syn_col, include_synonyms: bool = True) -> dict:
    """
    Build {lowercased_name: [(original_name, cvcl_id, source), ...]} for a
    subset of cellosaurus. `source` is 'primary' or 'synonym' so a collision
    report can say WHERE the shared name came from.
    """
    out = {}

    for _, row in sub.iterrows():
        cvcl = str(row.get("cvcl_id", "")).strip()

        primary = str(row.get("cell_line_name", "")).strip()
        if primary and primary.lower() != "nan":
            out.setdefault(primary.lower(), []).append((primary, cvcl, "primary"))

        if include_synonyms and syn_col is not None:
            raw_syn = row.get(syn_col)
            if pd.notna(raw_syn):
                for part in str(raw_syn).split(";"):
                    part = part.strip()
                    if part and part.lower() != "nan":
                        out.setdefault(part.lower(), []).append((part, cvcl, "synonym"))

    return out


# ============================================================
# Collision analysis
# ============================================================

def find_collisions(human_map: dict, non_human_map: dict) -> pd.DataFrame:
    """
    Every lowercased name present in BOTH maps, with where it came from on
    each side and the CVCL IDs involved.
    """
    shared = set(human_map) & set(non_human_map)
    if not shared:
        return pd.DataFrame(columns=[
            "name_lc", "n_chars", "human_sources", "human_cvcls",
            "non_human_sources", "non_human_cvcls", "collision_type",
        ])

    rows = []
    for key in sorted(shared):
        h = human_map[key]
        nh = non_human_map[key]

        h_src = sorted({s for _, _, s in h})
        nh_src = sorted({s for _, _, s in nh})

        if "primary" in h_src and "primary" in nh_src:
            ctype = "primary<->primary"
        elif "primary" in h_src:
            ctype = "human primary <-> non-human synonym"
        elif "primary" in nh_src:
            ctype = "human synonym <-> non-human primary"
        else:
            ctype = "synonym<->synonym"

        rows.append({
            "name_lc": key,
            "n_chars": len(key),
            "human_sources": "; ".join(h_src),
            "human_cvcls": "; ".join(sorted({c for _, c, _ in h if c})[:5]),
            "non_human_sources": "; ".join(nh_src),
            "non_human_cvcls": "; ".join(sorted({c for _, c, _ in nh if c})[:5]),
            "collision_type": ctype,
        })

    return pd.DataFrame(rows).sort_values(["n_chars", "name_lc"]).reset_index(drop=True)


def impact_on_tables(collisions: pd.DataFrame) -> pd.DataFrame:
    """
    For each cleaned table, count how many rows would be dropped by a
    colliding name — i.e. rows that are plausibly HUMAN but match the
    non-human blocklist purely because the name is shared.
    """
    if collisions.empty:
        return pd.DataFrame()

    colliding = set(collisions["name_lc"])
    rows = []

    for path in sorted(PARQUET_CLEAN.glob("*.parquet")):
        table = path.stem
        if table == "cellosaurus":
            continue
        try:
            df = pd.read_parquet(path)
        except Exception as e:
            print(f"  [skip] {table}: {type(e).__name__}: {e}")
            continue

        cols = [c for c in NAME_COLS if c in df.columns]
        if not cols:
            continue

        mask = pd.Series(False, index=df.index)
        for col in cols:
            mask |= df[col].fillna("").astype(str).str.strip().str.lower().isin(colliding)

        n = int(mask.sum())
        if n:
            hit_names = sorted(set(
                df.loc[mask, cols[0]].dropna().astype(str).str.strip()
            ))[:5]
            rows.append({
                "table": table,
                "rows_at_risk": n,
                "total_rows": len(df),
                "pct": round(n / len(df) * 100, 3) if len(df) else 0,
                "example_names": "; ".join(hit_names),
            })

    return pd.DataFrame(rows).sort_values("rows_at_risk", ascending=False).reset_index(drop=True)


# ============================================================
# Main
# ============================================================

def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--min-len", type=int, default=None,
                        help="Simulate a blocklist length guard: report how many "
                             "collisions a min-length filter would eliminate.")
    parser.add_argument("--save", action="store_true",
                        help="Write collision tables to logs/audit_overlap/")
    parser.add_argument("--check-tables", action="store_true",
                        help="Also scan cleaned parquet tables for rows at risk.")
    args = parser.parse_args()

    print("=" * 78)
    print("CELLOSAURUS HUMAN / NON-HUMAN NAME OVERLAP AUDIT")
    print("=" * 78)

    try:
        human, non_human, species_col, syn_col = load_cellosaurus()
    except (FileNotFoundError, KeyError) as e:
        print(f"\n  [FATAL] {e}")
        return 2

    print(f"\n  species column : {species_col}")
    print(f"  synonym column : {syn_col if syn_col else '(none found)'}")
    print(f"  human rows     : {len(human):,}")
    print(f"  non-human rows : {len(non_human):,}")

    # ── Primary names only ───────────────────────────────────────────────
    h_primary = name_map(human, syn_col, include_synonyms=False)
    nh_primary = name_map(non_human, syn_col, include_synonyms=False)
    primary_only = find_collisions(h_primary, nh_primary)

    print(f"\n{'-' * 78}")
    print("1. PRIMARY NAME COLLISIONS (cell_line_name only, synonyms ignored)")
    print(f"{'-' * 78}")
    print(f"  {len(primary_only):,} names used by both a human and a non-human line.")
    if not primary_only.empty:
        print(primary_only[["name_lc", "n_chars", "human_cvcls", "non_human_cvcls"]]
              .head(25).to_string(index=False))
        if len(primary_only) > 25:
            print(f"  ... and {len(primary_only) - 25:,} more")

    # ── Including synonyms ───────────────────────────────────────────────
    h_all = name_map(human, syn_col, include_synonyms=True)
    nh_all = name_map(non_human, syn_col, include_synonyms=True)
    all_collisions = find_collisions(h_all, nh_all)

    print(f"\n{'-' * 78}")
    print("2. ALL COLLISIONS (names + synonyms, both sides)")
    print(f"{'-' * 78}")
    print(f"  human names+synonyms     : {len(h_all):,}")
    print(f"  non-human names+synonyms : {len(nh_all):,}")
    print(f"  COLLIDING                : {len(all_collisions):,}")

    if not all_collisions.empty:
        print("\n  By collision type:")
        for ctype, n in all_collisions["collision_type"].value_counts().items():
            print(f"    {ctype:<40} {n:>7,}")

        print("\n  Shortest colliding names (highest false-positive risk):")
        print(all_collisions[["name_lc", "n_chars", "collision_type",
                              "human_cvcls", "non_human_cvcls"]]
              .head(25).to_string(index=False))

    # ── Length guard simulation ──────────────────────────────────────────
    if args.min_len is not None and not all_collisions.empty:
        kept = all_collisions[all_collisions["n_chars"] >= args.min_len]
        removed = len(all_collisions) - len(kept)
        print(f"\n{'-' * 78}")
        print(f"3. LENGTH GUARD SIMULATION (min_len={args.min_len})")
        print(f"{'-' * 78}")
        print(f"  {removed:,} collisions eliminated by dropping blocklist entries "
              f"shorter than {args.min_len} chars.")
        print(f"  {len(kept):,} collisions would REMAIN.")
        if not kept.empty:
            print("\n  Remaining (these need name-based filtering disabled, or "
                  "CVCL-only matching):")
            print(kept[["name_lc", "n_chars", "collision_type"]]
                  .head(20).to_string(index=False))

    # ── Downstream impact ────────────────────────────────────────────────
    if args.check_tables:
        print(f"\n{'-' * 78}")
        print("4. DOWNSTREAM IMPACT (rows in cleaned tables matching a colliding name)")
        print(f"{'-' * 78}")
        if not PARQUET_CLEAN.exists():
            print(f"  [skip] {PARQUET_CLEAN} does not exist — run 01_data_cleaning first.")
        else:
            impact = impact_on_tables(all_collisions)
            if impact.empty:
                print("  No cleaned-table rows match any colliding name. "
                      "Name-based filtering is safe on this data.")
            else:
                print(impact.to_string(index=False))
                print(f"\n  TOTAL rows at risk: {impact['rows_at_risk'].sum():,}")
                print("  These rows may be HUMAN but would be dropped by name matching.")

    # ── Save ─────────────────────────────────────────────────────────────
    if args.save:
        out_dir = project_root / "logs" / "audit_overlap"
        out_dir.mkdir(parents=True, exist_ok=True)
        primary_only.to_csv(out_dir / "collisions_primary_only.csv", index=False)
        all_collisions.to_csv(out_dir / "collisions_all.csv", index=False)
        print(f"\n  Wrote collision tables to {out_dir}")
        if args.check_tables and not all_collisions.empty and PARQUET_CLEAN.exists():
            impact_on_tables(all_collisions).to_csv(
                out_dir / "downstream_impact.csv", index=False)

    # ── Verdict ──────────────────────────────────────────────────────────
    print(f"\n{'=' * 78}")
    if all_collisions.empty:
        print("VERDICT: No overlap. Name-based non-human filtering is safe.")
        return 0

    n_primary = len(primary_only)
    n_short = int((all_collisions["n_chars"] <= 3).sum())
    print(f"VERDICT: {len(all_collisions):,} colliding names "
          f"({n_primary:,} primary-to-primary, {n_short:,} are <=3 chars).")
    print("Name-based filtering WILL drop some human rows. Options:")
    print("  (a) Filter on cvcl_id only — accessions are unique, names are not.")
    print("  (b) Keep name filtering but subtract human names from the blocklist.")
    print("  (c) Add a min-length guard (see --min-len) to cut the worst cases.")
    print("=" * 78)
    return 1

if __name__ == "__main__":
    sys.exit(main())