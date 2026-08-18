"""
api/queries.py
--------------
JSON-returning query functions for the Lambda handler.
Uses DuckDB to query predictions_with_confidence.parquet directly —
reads only the ~1,500 rows for the requested gene instead of loading
all 28M rows into pandas memory.
"""
from __future__ import annotations
import os
from pathlib import Path

import duckdb
import pandas as pd

_BASE = Path(os.environ.get("DATA_DIR", str(Path(__file__).resolve().parents[1] / "outputs")))
_REF  = Path(os.environ.get("REF_DIR",  str(Path(__file__).resolve().parents[1] / "reference")))
_SIM  = Path(os.environ.get("SIM_DIR",  str(Path(__file__).resolve().parents[2] / "cell_similarity" / "outputs")))

PREDICTIONS = _BASE / "predictions_with_confidence.parquet"
GENE_LKP    = _REF  / "gene_lookup.parquet"
CELL_LKP    = _REF  / "cell_line_lookup.parquet"
RNA_NBRS    = _SIM  / "rna_neighbours.parquet"

TOP_N    = 30
SIM_TOPK = 5

# Small lookup tables loaded once into memory (~2 MB total)
_genes = None
_lines = None
_nbrs  = None


def _load_lookups():
    global _genes, _lines, _nbrs
    if _genes is not None:
        return
    _genes = pd.read_parquet(GENE_LKP, columns=["ensg_id", "hgnc_symbol", "chromosomal_location"])
    _lines = pd.read_parquet(CELL_LKP, columns=["model_id", "cell_line_name"])
    _lines["model_id"] = _lines["model_id"].str.lower()
    _nbrs  = pd.read_parquet(RNA_NBRS) if RNA_NBRS.exists() else pd.DataFrame()


def _con():
    """Fresh read-only DuckDB connection — cheap to open, no persistent state needed."""
    return duckdb.connect()


def _fetch_gene(ensg: str) -> pd.DataFrame:
    """Query only the rows for one gene from the parquet — ~1,500 rows, not 28M."""
    con = _con()
    df = con.execute(
        f"""
        SELECT model_id, ensg_id, core_score, n_layers, confidence_tier,
               has_driver_alteration, p_mutation, p_fusion, has_cna_alteration
        FROM read_parquet('{PREDICTIONS}')
        WHERE ensg_id = ?
        """,
        [ensg]
    ).df()
    con.close()
    df["model_id"] = df["model_id"].str.lower()
    return df


def _resolve(query: str) -> tuple[str, str]:
    q = query.strip().upper()
    m = _genes[_genes["hgnc_symbol"].str.upper() == q]
    if not m.empty:
        return m.iloc[0]["ensg_id"], m.iloc[0]["hgnc_symbol"]
    m2 = _genes[_genes["ensg_id"].str.upper() == q]
    if not m2.empty:
        return m2.iloc[0]["ensg_id"], m2.iloc[0]["hgnc_symbol"]
    raise ValueError(f"Gene not found: {query!r}")


def _y_linked_set() -> set:
    loc = _genes["chromosomal_location"].fillna("")
    return set(_genes.loc[loc.str.startswith("Y"), "ensg_id"])


def _similar(model_id: str) -> list[dict]:
    if _nbrs.empty:
        return []
    sub = (_nbrs[_nbrs["model_id"].str.lower() == model_id.lower()]
           .sort_values("rank")
           .head(SIM_TOPK)
           .merge(_lines[["model_id", "cell_line_name"]].rename(
               columns={"model_id": "neighbour"}), on="neighbour", how="left"))
    return [
        {"model_id": r["neighbour"].upper(),
         "name": str(r.get("cell_line_name", "?")),
         "similarity": round(float(r.get("similarity", 0)), 3)}
        for _, r in sub.iterrows()
    ]


# ── public query functions ─────────────────────────────────────────────────────

def query_gene(gene: str, cell_line: str | None = None, top_n: int = TOP_N) -> dict:
    _load_lookups()
    ensg, sym = _resolve(gene)
    sub = _fetch_gene(ensg)
    if sub.empty:
        return {"error": f"No predictions for {gene}"}

    sub = sub.sort_values("core_score", ascending=False).reset_index(drop=True)
    sub = sub.merge(_lines, on="model_id", how="left")

    if cell_line:
        row = sub[sub.model_id.str.lower() == cell_line.lower()]
        if row.empty:
            return {"error": f"No prediction for {gene} + {cell_line}"}
        row = row.iloc[0]
        rank_pos = int(sub[sub.model_id.str.lower() == cell_line.lower()].index[0]) + 1
        return {
            "gene": sym, "ensg": ensg,
            "cell_line": {"model_id": row["model_id"].upper(),
                          "name": str(row.get("cell_line_name", cell_line))},
            "rank": rank_pos, "total": len(sub),
            "score": round(float(row["core_score"]), 4),
            "tier": row.get("confidence_tier"),
            "n_layers": int(row.get("n_layers", 0)),
            "driver_alteration": bool(row.get("has_driver_alteration", False)),
            "rna_alternatives": _similar(cell_line),
        }

    rows = [
        {
            "rank": i + 1,
            "model_id": r["model_id"].upper(),
            "name": str(r.get("cell_line_name", "?")),
            "score": round(float(r["core_score"]), 4),
            "tier": r.get("confidence_tier"),
            "n_layers": int(r.get("n_layers", 0)),
            "driver_alteration": bool(r.get("has_driver_alteration", False)),
        }
        for i, (_, r) in enumerate(sub.head(top_n).iterrows())
    ]

    return {
        "gene": sym, "ensg": ensg,
        "total": len(sub), "showing": len(rows),
        "lines": rows,
        "rna_alternatives_for_rank1": _similar(sub.iloc[0]["model_id"]),
    }


def query_genes(gene_list: list[str], top_n: int = TOP_N) -> dict:
    _load_lookups()
    FLOOR = 0.50
    resolved = [_resolve(g) for g in gene_list]
    ensgs = [e for e, _ in resolved]
    syms  = [s for _, s in resolved]

    wide = None
    for ensg, sym in zip(ensgs, syms):
        s = _fetch_gene(ensg)[["model_id", "core_score"]].rename(columns={"core_score": sym})
        wide = s if wide is None else wide.merge(s, on="model_id", how="inner")

    if wide is None or wide.empty:
        return {"error": "No shared cell lines found across all genes"}

    wide["joint_score"] = wide[syms].min(axis=1)
    wide = wide[wide["joint_score"] >= FLOOR].sort_values("joint_score", ascending=False)
    wide = wide.merge(_lines, on="model_id", how="left").reset_index(drop=True)

    rows = [
        {
            "model_id": r["model_id"].upper(),
            "name": str(r.get("cell_line_name", "?")),
            "joint_score": round(float(r["joint_score"]), 4),
            "scores": {sym: round(float(r[sym]), 4) for sym in syms},
        }
        for _, r in wide.head(top_n).iterrows()
    ]
    return {"genes": syms, "floor": FLOOR, "total_passing": len(wide), "lines": rows}


def query_exclude(gene_a: str, gene_b: str, top_n: int = TOP_N) -> dict:
    _load_lookups()
    ensg_a, sym_a = _resolve(gene_a)
    ensg_b, sym_b = _resolve(gene_b)

    a = _fetch_gene(ensg_a)[["model_id", "core_score"]].rename(columns={"core_score": "score_a"})
    b = _fetch_gene(ensg_b)[["model_id", "core_score"]].rename(columns={"core_score": "score_b"})

    merged = a.merge(b, on="model_id", how="inner")
    merged["selectivity"] = merged["score_a"] * (1.0 - merged["score_b"])
    merged = merged.sort_values("selectivity", ascending=False).reset_index(drop=True)
    merged = merged.merge(_lines, on="model_id", how="left")

    rows = [
        {
            "model_id": r["model_id"].upper(),
            "name": str(r.get("cell_line_name", "?")),
            f"score_{sym_a}": round(float(r["score_a"]), 4),
            f"score_{sym_b}": round(float(r["score_b"]), 4),
            "selectivity": round(float(r["selectivity"]), 4),
        }
        for _, r in merged.head(top_n).iterrows()
    ]
    return {
        "gene_high": sym_a, "gene_low": sym_b,
        "formula": f"score_{sym_a} × (1 − score_{sym_b})",
        "total_ranked": len(merged),
        "lines": rows,
    }
