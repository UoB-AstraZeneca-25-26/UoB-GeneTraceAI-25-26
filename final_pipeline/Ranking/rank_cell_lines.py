"""
Ranking/rank_cell_lines.py
---------------------------
Rank cell lines for a gene. Output: ACH ID | Cell Line Name | Score.
Use cli.py for the full interactive interface (detail view, co-selection, exclusion).

Usage:
  python rank_cell_lines.py BRAF
  python rank_cell_lines.py ENSG00000157764
"""
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from config import PREDICTIONS, CELL_LKP, GENE_LKP

TOP_N = 30


def resolve_gene(query: str, gene_lkp: pd.DataFrame) -> tuple[str, str]:
    q = query.strip().upper()
    match = gene_lkp[gene_lkp["hgnc_symbol"].str.upper() == q]
    if not match.empty:
        row = match.iloc[0]
        return row["ensg_id"], row["hgnc_symbol"]
    q2 = query.strip().lower()
    match2 = gene_lkp[gene_lkp["ensg_id"].str.lower() == q2]
    if not match2.empty:
        row = match2.iloc[0]
        return row["ensg_id"], row["hgnc_symbol"]
    raise ValueError(f"Gene not found: {query}")


def rank(gene_query: str) -> pd.DataFrame:
    for p in [PREDICTIONS, CELL_LKP, GENE_LKP]:
        if not p.exists():
            raise SystemExit(f"Missing: {p} — run scoring stages first.")

    genes = pd.read_parquet(GENE_LKP, columns=["ensg_id", "hgnc_symbol"])
    lines = pd.read_parquet(CELL_LKP, columns=["model_id", "cell_line_name"])
    lines["model_id"] = lines["model_id"].str.lower()
    ensg, gene_sym = resolve_gene(gene_query, genes)

    pred = pd.read_parquet(PREDICTIONS)
    sub  = pred[pred.ensg_id == ensg].copy()
    if sub.empty:
        print(f"No predictions for {gene_sym} ({ensg})")
        return pd.DataFrame()

    sub = sub.sort_values("core_score", ascending=False).reset_index(drop=True)
    sub = sub.merge(lines, on="model_id", how="left")

    print(f"\nGene: {gene_sym}  ({ensg})")
    print("-" * 60)
    print(f"  {'ACH ID':<13} {'Cell Line Name':<30} {'Score':>6}  Tier")
    print("-" * 60)
    for _, row in sub.head(TOP_N).iterrows():
        name = str(row.get("cell_line_name", "?"))[:29]
        print(f"  {row['model_id']:<13} {name:<30} {row['core_score']:>6.3f}  {row.get('confidence_tier', '?')}")
    print(f"\n  Showing top {min(TOP_N, len(sub))} of {len(sub):,} lines")
    return sub


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python rank_cell_lines.py <GENE_SYMBOL_OR_ENSG>")
        sys.exit(1)
    rank(sys.argv[1])
