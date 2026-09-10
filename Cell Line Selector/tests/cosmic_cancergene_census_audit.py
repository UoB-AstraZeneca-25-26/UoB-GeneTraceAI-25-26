#!/usr/bin/env python3
"""
Audit raw COSMIC Cancer Gene Census TSV for nulls and duplicate mappings.

This checks the key identity columns used in the raw file:
- GENE_SYMBOL
- COSMIC_GENE_ID

It also checks the reverse relationship:
- GENE_SYMBOL -> number of distinct COSMIC_GENE_ID values
- COSMIC_GENE_ID -> number of distinct GENE_SYMBOL values

Usage:
    python cosmic_cancergene_census_audit.py
    python cosmic_cancergene_census_audit.py "data/raw data/cosmic/Cosmic_CancerGeneCensus_v104_GRCh37.tsv"
"""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


def default_raw_path() -> Path:
    """Return the canonical raw Cancer Gene Census TSV path."""
    here = Path(__file__).resolve().parent
    candidates = [
        here / "data" / "raw data" / "cosmic" / "Cosmic_CancerGeneCensus_v104_GRCh37.tsv",
        here / "data" / "raw data" / "Cosmic_CancerGeneCensus_v104_GRCh37.tsv",
    ]
    for candidate in candidates:
        if candidate.exists():
            return candidate
    return candidates[0]


def summarize_column(series: pd.Series) -> dict:
    s = series.fillna("").astype(str).str.strip()
    duplicated = s[s.duplicated(keep=False)]
    return {
        "total_rows": int(len(s)),
        "nulls": int(series.isna().sum()),
        "unique_values": int(s.nunique()),
        "duplicate_rows": int(duplicated.shape[0]),
        "duplicate_unique_values": int(duplicated.nunique()),
        "is_unique": bool(s.nunique() == len(s)),
    }


def audit_cancer_gene_census(raw_path: str | Path) -> dict:
    """Read and audit the Cancer Gene Census raw TSV."""
    path = Path(raw_path)
    if not path.exists():
        raise FileNotFoundError(f"Cancer Gene Census raw file not found: {path}")

    raw = pd.read_csv(path, sep="\t", low_memory=False)

    for col in ["GENE_SYMBOL", "COSMIC_GENE_ID"]:
        if col not in raw.columns:
            raise ValueError(f"Missing expected column: {col}")

    gene_to_geneids = (
        raw.groupby("GENE_SYMBOL")["COSMIC_GENE_ID"].nunique().reset_index(name="n_cosmic_gene_ids")
    )
    geneid_to_genes = (
        raw.groupby("COSMIC_GENE_ID")["GENE_SYMBOL"].nunique().reset_index(name="n_gene_symbols")
    )

    return {
        "file": str(path),
        "shape": raw.shape,
        "columns": {
            "GENE_SYMBOL": summarize_column(raw["GENE_SYMBOL"]),
            "COSMIC_GENE_ID": summarize_column(raw["COSMIC_GENE_ID"]),
        },
        "gene_symbol_to_cosmic_gene_id": gene_to_geneids,
        "cosmic_gene_id_to_gene_symbol": geneid_to_genes,
    }


def print_summary(summary: dict) -> None:
    print(f"File: {summary['file']}")
    print(f"Shape: {summary['shape'][0]} rows x {summary['shape'][1]} cols")

    for col in ["GENE_SYMBOL", "COSMIC_GENE_ID"]:
        info = summary["columns"][col]
        print(f"\n{col}:")
        print(f"  nulls: {info['nulls']}")
        print(f"  unique values: {info['unique_values']}")
        print(f"  duplicate rows: {info['duplicate_rows']}")
        print(f"  is_unique: {info['is_unique']}")

    gene_to_geneids = summary["gene_symbol_to_cosmic_gene_id"]
    multi_gene = gene_to_geneids[gene_to_geneids["n_cosmic_gene_ids"] > 1]
    geneid_to_genes = summary["cosmic_gene_id_to_gene_symbol"]
    multi_id = geneid_to_genes[geneid_to_genes["n_gene_symbols"] > 1]

    print("\nGENE_SYMBOL -> distinct COSMIC_GENE_ID counts:")
    print(gene_to_geneids.sort_values("n_cosmic_gene_ids", ascending=False).head(10).to_string(index=False))
    print(f"\nGenes with >1 COSMIC_GENE_ID: {len(multi_gene)}")
    if len(multi_gene) == 0:
        print("Result: no gene_symbol maps to more than one COSMIC_GENE_ID.")
    else:
        print(multi_gene.head(20).to_string(index=False))

    print("\nCOSMIC_GENE_ID -> distinct GENE_SYMBOL counts:")
    print(geneid_to_genes.sort_values("n_gene_symbols", ascending=False).head(10).to_string(index=False))
    print(f"\nGene IDs with >1 GENE_SYMBOL: {len(multi_id)}")
    if len(multi_id) == 0:
        print("Result: no COSMIC_GENE_ID maps to more than one GENE_SYMBOL.")
    else:
        print(multi_id.head(20).to_string(index=False))


def main() -> None:
    parser = argparse.ArgumentParser(description="Audit nulls and duplicate mappings in the Cancer Gene Census TSV.")
    parser.add_argument(
        "path",
        nargs="?",
        default=None,
        help="Path to the raw Cosmic_CancerGeneCensus_v104_GRCh37.tsv file.",
    )
    args = parser.parse_args()

    raw_path = args.path if args.path else default_raw_path()
    summary = audit_cancer_gene_census(raw_path)
    print_summary(summary)


if __name__ == "__main__":
    main()
