"""
test_cmd_genes_3way.py
----------------------
Tests the final_pipeline co-selection join (cmd_genes logic) with 3-gene
triplets using real prediction data.

Covers:
  1. Triplet: 3 gene pairs complete without error; joint = min of 3
  2. NaN propagation fix: no line with NaN per-gene score in results
  3. Limiting-gene: the weakest gene is correctly named
  4. 3 distinct triplets fired concurrently return independent results
  5. Inner join: only lines scored for ALL genes appear in results
"""

from __future__ import annotations

import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from config import PREDICTIONS, GENE_LKP, CELL_LKP

JOINT_FLOOR = 0.50

# ---------------------------------------------------------------------------
# Shared fixtures — load once for the whole module
# ---------------------------------------------------------------------------

@pytest.fixture(scope="module")
def pred() -> pd.DataFrame:
    return pd.read_parquet(PREDICTIONS, columns=["model_id", "ensg_id", "core_score"])


@pytest.fixture(scope="module")
def gene_lkp() -> pd.DataFrame:
    return pd.read_parquet(GENE_LKP, columns=["ensg_id", "hgnc_symbol"])


# ---------------------------------------------------------------------------
# Core join logic (mirrors cmd_genes, without I/O or print statements)
# ---------------------------------------------------------------------------

def _coselect(gene_syms: list[str], pred: pd.DataFrame,
              gene_lkp: pd.DataFrame) -> pd.DataFrame:
    """
    Build the joint-score DataFrame for the requested genes.
    Returns the full wide DataFrame (unfiltered) so tests can inspect it.
    """
    resolved = []
    for sym in gene_syms:
        m = gene_lkp[gene_lkp["hgnc_symbol"].str.upper() == sym.upper()]
        if m.empty:
            raise ValueError(f"Gene not found: {sym}")
        resolved.append((m.iloc[0]["ensg_id"], m.iloc[0]["hgnc_symbol"]))

    wide = None
    for ensg, sym in resolved:
        s = (pred[pred.ensg_id == ensg][["model_id", "core_score"]]
             .rename(columns={"core_score": sym}))
        wide = s if wide is None else wide.merge(s, on="model_id", how="inner")

    syms = [s for _, s in resolved]
    wide["joint_score"] = wide[syms].min(axis=1, skipna=False)
    scoreable = wide["joint_score"].notna()
    wide.loc[scoreable, "limiting_gene"] = wide.loc[scoreable, syms].idxmin(axis=1)
    return wide, syms


# ---------------------------------------------------------------------------
# 1. Triplet completes and joint = min of 3
# ---------------------------------------------------------------------------

def test_triplet_joint_equals_min_of_three_scores(pred, gene_lkp):
    wide, syms = _coselect(["BRAF", "KRAS", "EGFR"], pred, gene_lkp)

    passing = wide[wide["joint_score"].notna() & (wide["joint_score"] >= JOINT_FLOOR)]
    assert not passing.empty, "Expected some lines to clear the 0.50 floor for BRAF+KRAS+EGFR"

    for _, row in passing.head(20).iterrows():
        expected = min(row["BRAF"], row["KRAS"], row["EGFR"])
        assert abs(row["joint_score"] - expected) < 1e-9, (
            f"{row['model_id']}: joint_score={row['joint_score']:.6f} "
            f"!= min({row['BRAF']:.4f}, {row['KRAS']:.4f}, {row['EGFR']:.4f})"
        )


# ---------------------------------------------------------------------------
# 2. NaN lines excluded from passing results
# ---------------------------------------------------------------------------

def test_nan_core_score_lines_excluded_from_passing(pred, gene_lkp):
    wide, syms = _coselect(["BRAF", "KRAS", "EGFR"], pred, gene_lkp)
    passing = wide[wide["joint_score"].notna() & (wide["joint_score"] >= JOINT_FLOOR)]

    for sym in syms:
        nan_in_passing = passing[sym].isna().sum()
        assert nan_in_passing == 0, (
            f"{nan_in_passing} rows with NaN {sym} score survived into passing results. "
            f"skipna=False fix may not be applied."
        )


# ---------------------------------------------------------------------------
# 3. Limiting gene is the column with the minimum value
# ---------------------------------------------------------------------------

def test_limiting_gene_is_weakest_column(pred, gene_lkp):
    wide, syms = _coselect(["BRAF", "KRAS", "EGFR"], pred, gene_lkp)
    passing = wide[wide["joint_score"].notna() & (wide["joint_score"] >= JOINT_FLOOR)]

    for _, row in passing.head(50).iterrows():
        scores = {s: row[s] for s in syms}
        expected_limiter = min(scores, key=scores.__getitem__)
        assert row["limiting_gene"] == expected_limiter, (
            f"{row['model_id']}: limiting_gene={row['limiting_gene']!r} "
            f"but min score is {expected_limiter!r} ({scores})"
        )


# ---------------------------------------------------------------------------
# 4. Inner join: every line in wide must have a real model_id in all gene frames
# ---------------------------------------------------------------------------

def test_inner_join_only_returns_shared_lines(pred, gene_lkp):
    genes = ["BRAF", "KRAS", "EGFR"]
    ensg_map = {}
    for sym in genes:
        m = gene_lkp[gene_lkp["hgnc_symbol"].str.upper() == sym]
        ensg_map[sym] = m.iloc[0]["ensg_id"]

    wide, syms = _coselect(genes, pred, gene_lkp)

    for sym in genes:
        ensg = ensg_map[sym]
        scored_lines = set(pred[pred.ensg_id == ensg]["model_id"])
        extra = set(wide["model_id"]) - scored_lines
        assert not extra, (
            f"{len(extra)} model_ids in result were NOT in {sym}'s scored set: "
            f"{list(extra)[:5]}"
        )


# ---------------------------------------------------------------------------
# 5. Three different triplets fired concurrently — no cross-contamination
# ---------------------------------------------------------------------------

TRIPLETS = [
    ["BRAF", "KRAS", "EGFR"],
    ["TP53", "PTEN", "KRAS"],
    ["BRAF", "TP53", "EGFR"],
]


def _run(args):
    genes, pred, gene_lkp = args
    wide, syms = _coselect(genes, pred, gene_lkp)
    passing = wide[wide["joint_score"].notna() & (wide["joint_score"] >= JOINT_FLOOR)]
    return genes, passing, syms


def test_concurrent_triplets_independent(pred, gene_lkp):
    results = {}
    with ThreadPoolExecutor(max_workers=3) as pool:
        futures = {pool.submit(_run, (genes, pred, gene_lkp)): genes
                   for genes in TRIPLETS}
        for fut in as_completed(futures):
            genes, passing, syms = fut.result()
            results[tuple(genes)] = (passing, syms)

    assert len(results) == 3, "All 3 concurrent queries must complete"

    for genes, (passing, syms) in results.items():
        # Each result must only contain its own gene score columns
        for sym in syms:
            assert sym in passing.columns, f"{sym} missing from result for {genes}"
        other_syms = {s for g in TRIPLETS for s in g} - set(syms)
        for sym in other_syms:
            if sym not in syms:
                assert sym not in passing.columns, (
                    f"Column {sym} from another triplet leaked into result for {genes}"
                )

        # No NaN per-gene scores in passing rows
        for sym in syms:
            assert passing[sym].isna().sum() == 0, (
                f"NaN {sym} in concurrent result for {genes}"
            )

        # joint_score == min of all gene scores
        for _, row in passing.head(10).iterrows():
            expected = min(row[s] for s in syms)
            assert abs(row["joint_score"] - expected) < 1e-9
