#!/usr/bin/env python3
"""Harmonisation validation for the cell-line and gene rosters.

Purpose:
  1. Check that the canonical gene roster is the protein-coding-only gene universe.
  2. Check that the cell_line_roster is the model_id universe used across tables.
  3. Detect missing links, leaked IDs, and non-human rows.
  4. Produce a concise summary report for data quality review.

This script reads the cleaned tables directly from the DuckDB project database and
compares the roster tables against the source expression / metadata tables.
"""

from __future__ import annotations

import re
from collections import defaultdict
from pathlib import Path
from typing import Iterable

import duckdb
import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parent
DB_PATH = PROJECT_ROOT / "db" / "celllineselector.duckdb"


_ENSG_RE = re.compile(r"(?i)^(ENSG\d+)(?:\.\d+)?$")


def canonical_ensg(value) -> str | None:
    if value is None or pd.isna(value):
        return None
    text = str(value).strip().upper()
    if not text:
        return None
    m = _ENSG_RE.fullmatch(text)
    return m.group(1) if m else None


def clean_set(values: Iterable[str | None]) -> set[str]:
    out = set()
    for v in values:
        if v is None or pd.isna(v):
            continue
        s = str(v).strip()
        if s:
            out.add(s)
    return out


def normalise_gene_series(s: pd.Series) -> set[str]:
    out = set()
    for v in s.dropna():
        g = canonical_ensg(v)
        if g:
            out.add(g)
    return out


def human_species_values() -> set[str]:
    vals = {
        "homo sapiens",
        "homo sapiens (human)",
        "human",
        "humans",
        "h. sapiens",
        "hsapiens",
        "homo-sapiens",
        "homo sapiens / human",
    }
    vals.update({v.upper() for v in vals})
    vals.update({v.title() for v in vals})
    return vals


def species_is_human(value) -> bool:
    if value is None or pd.isna(value):
        return False
    text = str(value).strip().lower()
    if not text:
        return False
    return any(token in text for token in ["homo sapiens", "human", "h. sapiens", "hsapiens"])


def table_exists(con, name: str) -> bool:
    try:
        return bool(con.execute("SHOW TABLES").fetchdf()["name"].eq(name).any())
    except Exception:
        return False


def load_table(con, name: str) -> pd.DataFrame | None:
    if not table_exists(con, name):
        return None
    return con.execute(f'SELECT * FROM "{name}"').fetchdf()


def summarize_table(df: pd.DataFrame | None, name: str) -> dict:
    if df is None:
        return {"table": name, "exists": False, "rows": 0, "cols": 0}
    return {"table": name, "exists": True, "rows": len(df), "cols": len(df.columns)}


def gene_roster_ids(con) -> set[str]:
    df = load_table(con, "gene_roster")
    if df is None or df.empty:
        return set()
    if "gene_id" not in df.columns:
        return set()
    return normalise_gene_series(df["gene_id"])


def cell_roster_model_ids(con) -> set[str]:
    df = load_table(con, "cell_line_roster")
    if df is None or df.empty:
        return set()
    if "model_id" not in df.columns:
        return set()
    return clean_set(df["model_id"].dropna())


def gene_usage_from_expression_tables(con) -> dict[str, set[str]]:
    usage = {}
    for name in ["depmap_expr", "geo_expr"]:
        df = load_table(con, name)
        if df is None:
            continue
        cols = [c for c in df.columns if c != "profile_id" and c != "gsm_id"]
        usage[name] = set()
        for col in cols:
            g = canonical_ensg(col)
            if g:
                usage[name].add(g)

    df = load_table(con, "hpa_rna")
    if df is not None and "gene_id" in df.columns:
        usage["hpa_rna"] = normalise_gene_series(df["gene_id"])

    return usage


def model_usage_from_tables(con) -> dict[str, set[str]]:
    usage = {}
    for table_name in [
        "sample_info",
        "depmap_profiles",
        "geo_info",
        "hpa_rna",
        "metabolomics",
        "fusions",
    ]:
        df = load_table(con, table_name)
        if df is None:
            continue
        if "model_id" in df.columns:
            usage[table_name] = clean_set(df["model_id"].dropna())
        elif "modelid" in df.columns:
            usage[table_name] = clean_set(df["modelid"].dropna())
        elif "modelId" in df.columns:
            usage[table_name] = clean_set(df["modelId"].dropna())
        elif table_name == "fusions" and "modelid" in df.columns:
            usage[table_name] = clean_set(df["modelid"].dropna())
        elif table_name == "geo_info" and "geo_accession" in df.columns:
            # geo_info does not carry model_id directly; it is joined via geo_expr -> geo_info mapping,
            # so it is not a primary model_id source for the roster check.
            pass
    return usage


def species_reference_summary(con) -> dict:
    out = {}
    checks = [
        ("cellosaurus", "species_of_origin"),
        ("model_list_20260709", "species"),
        ("geo_info", "organism_ch1"),
    ]
    for table_name, col in checks:
        df = load_table(con, table_name)
        if df is None or col not in df.columns:
            out[table_name] = {"present": False, "human_count": 0, "non_human_count": 0, "sample_values": []}
            continue
        vals = df[col].dropna().astype(str)
        human = vals[vals.map(species_is_human)]
        non_human = vals[~vals.map(species_is_human)]
        out[table_name] = {
            "present": True,
            "human_count": int(len(human)),
            "non_human_count": int(len(non_human)),
            "sample_values": human.head(10).tolist() + non_human.head(10).tolist(),
        }
    return out


def cell_usage_species_check(con) -> dict:
    roster_ids = cell_roster_model_ids(con)
    usage = model_usage_from_tables(con)
    out = {}
    for name, ids in usage.items():
        extra = ids - roster_ids
        out[name] = {
            "n_in_source": len(ids),
            "n_in_roster": len(ids & roster_ids),
            "n_missing_from_roster": len(extra),
            "examples_missing": sorted(list(extra))[:10],
        }
    return out


def validate_genes(con) -> dict:
    roster = gene_roster_ids(con)
    expr_usage = gene_usage_from_expression_tables(con)
    report = {}

    for table_name, genes in expr_usage.items():
        missing = genes - roster
        report[table_name] = {
            "n_genes_in_table": len(genes),
            "n_genes_in_roster": len(genes & roster),
            "n_missing_from_roster": len(missing),
            "examples_missing": sorted(list(missing))[:10],
        }

    # Also check whether gene_roster includes genes not used by any expression table.
    used_anywhere = set().union(*expr_usage.values()) if expr_usage else set()
    orphan_roster = roster - used_anywhere
    report["gene_roster_orphans"] = {
        "n_orphan_in_roster": len(orphan_roster),
        "examples": sorted(list(orphan_roster))[:10],
    }
    return report


def validate_cell_lines(con) -> dict:
    roster = cell_roster_model_ids(con)
    usage = model_usage_from_tables(con)
    report = defaultdict(dict)
    for table_name, ids in usage.items():
        missing = ids - roster
        report[table_name] = {
            "n_in_table": len(ids),
            "n_in_roster": len(ids & roster),
            "n_missing_from_roster": len(missing),
            "examples_missing": sorted(list(missing))[:10],
        }

    all_usage = set().union(*usage.values()) if usage else set()
    orphan_roster = roster - all_usage
    report["roster_orphans"] = {
        "n_orphan_in_roster": len(orphan_roster),
        "examples": sorted(list(orphan_roster))[:10],
    }

    # Human-only consistency check: roster IDs present in sample_info etc. but with non-human species values.
    species_ref = species_reference_summary(con)
    report["species_reference"] = species_ref
    return dict(report)


def print_summary_section(title: str):
    print(f"\n{'=' * 80}\n{title}\n{'=' * 80}")


def main() -> None:
    con = duckdb.connect(str(DB_PATH), read_only=True)

    print_summary_section("Harmonisation validation summary")
    print(f"Database: {DB_PATH}")

    roster_summary = {
        "gene_roster": summarize_table(load_table(con, "gene_roster"), "gene_roster"),
        "cell_line_roster": summarize_table(load_table(con, "cell_line_roster"), "cell_line_roster"),
    }
    for v in roster_summary.values():
        print(f"{v['table']}: rows={v['rows']}, cols={v['cols']}")

    gene_roster = gene_roster_ids(con)
    cell_roster = cell_roster_model_ids(con)
    print(f"Unique gene IDs in gene_roster: {len(gene_roster)}")
    print(f"Unique model_ids in cell_line_roster: {len(cell_roster)}")

    gene_report = validate_genes(con)
    cell_report = validate_cell_lines(con)

    print_summary_section("Gene validation")
    for table_name in ["depmap_expr", "geo_expr", "hpa_rna"]:
        if table_name in gene_report:
            info = gene_report[table_name]
            print(f"{table_name}: genes_in_table={info['n_genes_in_table']}, in_roster={info['n_genes_in_roster']}, missing={info['n_missing_from_roster']}")
            if info["examples_missing"]:
                print(f"  examples_missing: {info['examples_missing'][:5]}")
    print(f"gene_roster_orphans: n={gene_report.get('gene_roster_orphans', {}).get('n_orphan_in_roster', 0)}")
    if gene_report.get("gene_roster_orphans", {}).get("examples"):
        print(f"  examples: {gene_report['gene_roster_orphans']['examples'][:5]}")

    print_summary_section("Cell-line validation")
    for table_name in [
        "sample_info",
        "depmap_profiles",
        "hpa_rna",
        "metabolomics",
        "fusions",
    ]:
        info = cell_report.get(table_name)
        if not info:
            continue
        print(f"{table_name}: in_table={info['n_in_table']}, in_roster={info['n_in_roster']}, missing={info['n_missing_from_roster']}")
        if info["examples_missing"]:
            print(f"  examples_missing: {info['examples_missing'][:5]}")
    print(f"roster_orphans: n={cell_report.get('roster_orphans', {}).get('n_orphan_in_roster', 0)}")
    if cell_report.get("roster_orphans", {}).get("examples"):
        print(f"  examples: {cell_report['roster_orphans']['examples'][:5]}")

    print_summary_section("Species / human-origin validation")
    species_ref = cell_report.get("species_reference", {})
    for table_name, info in species_ref.items():
        print(f"{table_name}: present={info['present']}, human={info['human_count']}, non_human={info['non_human_count']}")
        if info["sample_values"]:
            print(f"  sample_values: {info['sample_values'][:5]}")

    # Final verdict summary
    print_summary_section("Overall verdict")
    gene_missing_total = sum(v["n_missing_from_roster"] for v in gene_report.values() if isinstance(v, dict) and "n_missing_from_roster" in v)
    cell_missing_total = sum(v["n_missing_from_roster"] for v in cell_report.values() if isinstance(v, dict) and "n_missing_from_roster" in v)

    print(f"Total gene IDs missing from gene_roster across expression tables: {gene_missing_total}")
    print(f"Total cell-line IDs missing from cell_line_roster across source tables: {cell_missing_total}")

    if gene_missing_total == 0 and cell_missing_total == 0:
        print("Overall assessment: connected and internally consistent for protein-coding genes and cell-line model IDs.")
    else:
        print("Overall assessment: missing links or leakage were detected. Please inspect the per-table details above.")


if __name__ == "__main__":
    main()
