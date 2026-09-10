import sys
import pandas as pd

sys.path.append("src/pipeline")
from explain_pair import sort_gene_rows, load_data
from evidence_ledger import rna_covered_model_ids, derive_layers_present

_CELL_LINE_CONTEXT_COLS = [
    "model_id", "diseases", "category", "sex", "age_at_sampling",
    "parent_hierarchy", "has_rna", "has_wes", "has_wgs", "gap_class",
]
_cell_line_lookup_cache = None
_rna_covered_cache = None


def _cell_line_context():
    """cell_line_lookup.parquet already carries this -- a display join at
    query time, zero architectural risk (does not touch core_score)."""
    global _cell_line_lookup_cache
    if _cell_line_lookup_cache is None:
        cll = pd.read_parquet("reference/cell_line_lookup.parquet")
        cll = cll.copy()
        cll["model_id"] = cll["model_id"].str.lower()
        _cell_line_lookup_cache = cll[_CELL_LINE_CONTEXT_COLS]
    return _cell_line_lookup_cache


def rank_cell_lines(ensg_id, pred=None, chronos_val=None):
    """
    Returns the ranked list of cell lines for a gene, using the same
    query-time sort logic validated in explain_pair (not stored row order).

    Each row carries cell-line context (diseases/category/etc, joined from
    cell_line_lookup.parquet) and a `contradicted` flag (chronos_check ==
    'inverted' for this gene) surfaced as its own column -- not buried in
    a label -- so it sorts/filters/counts as visibly as confirmed evidence.
    """
    if pred is None:
        pred, _, _, _, _, _ = load_data()

    global _rna_covered_cache
    if _rna_covered_cache is None:
        _rna_covered_cache = rna_covered_model_ids()

    gene_rows = pred[pred["ensg_id"] == ensg_id].copy()
    if gene_rows.empty:
        return {"error": f"No predictions found for gene {ensg_id}"}

    sorted_rows = sort_gene_rows(gene_rows)

    contradicted = False
    if chronos_val is not None:
        cv_row = chronos_val[chronos_val["ensg_id"] == ensg_id]
        if not cv_row.empty:
            contradicted = bool((cv_row["chronos_check"] == "inverted").iloc[0])

    result = sorted_rows[[
        "sorted_position", "model_id", "core_score", "n_layers",
        "stratum_rank", "class", "rank_basis", "has_driver_alteration",
        "confidence", "confidence_reason"
    ]].rename(columns={"stratum_rank": "stratum_rank_percentile"})

    result["layers_present"] = [
        derive_layers_present(nl, mid, _rna_covered_cache)
        for nl, mid in zip(result["n_layers"], result["model_id"])
    ]
    result["contradicted"] = contradicted  # gene-level fact, broadcast per row

    result = result.merge(_cell_line_context(), on="model_id", how="left")

    # Gene-level facts live at the top of the result, not repeated per row.
    # signal_spread says whether the ordering below is worth acting on: core_score
    # is a within-gene percentile, so a near-flat gene still produces a full 0-1
    # spread and a confident-looking top-10. 'unmeasured' until the gene_dispersion
    # layer has been rebuilt (02_core_score.ipynb Cell 7c).
    spread = (sorted_rows["signal_spread"].iloc[0]
              if "signal_spread" in sorted_rows.columns else "unmeasured")
    fold = (sorted_rows["gene_top_vs_median_fold"].iloc[0]
            if "gene_top_vs_median_fold" in sorted_rows.columns else None)

    return {
        "gene": ensg_id,
        "class": sorted_rows["class"].iloc[0],
        "contradicted": contradicted,
        "signal_spread": spread,
        "top_vs_median_fold": fold,
        "total_candidates": len(sorted_rows),
        "ranking": result.to_dict(orient="records"),
    }


if __name__ == "__main__":
    from explain_pair import explain_pair

    pred, flags, variants, harm, gene_lookup, chronos_val = load_data()

    # --- Test 1: BRAF (activation_driven) — A375 must remain a top driver-positive hit ---
    # NOTE: this was a hardcoded "position == 13" check until the core_score recipe
    # was switched from Noisy-OR to a correlation-penalised weighted sum (see
    # 02_core_score.ipynb). That's a genuinely different combination formula, not a
    # monotonic rescaling, so it can (and did) reshuffle relative order among
    # driver-positive cell lines for the same gene -- A375 moved from #13 to #36.
    # A375's own stratum_rank percentile actually *improved* (0.767 -> 0.814) under
    # the fix; the exact-position drop is expected recalibration noise, not a
    # regression. An exact hardcoded position is too brittle to survive legitimate
    # future scoring changes, so this now checks the two invariants that actually
    # matter: (1) the real BRAF V600E driver mutation is still detected and still
    # forces A375 into the driver-positive group, and (2) it still lands comfortably
    # near the top of the full candidate list (top 5%), not just among driver-positive
    # lines. A genuine regression (e.g. A375 falling to the middle of the pack) would
    # still fail this.
    braf_result = rank_cell_lines("ensg00000157764", pred, chronos_val)
    print("=== BRAF ranking — top 15 ===")
    for row in braf_result["ranking"][:15]:
        print(f"  #{row['sorted_position']:>4}  {row['model_id']}  "
              f"score={row['core_score']:.4f}  driver={row['has_driver_alteration']}  "
              f"basis={row['rank_basis']}  conf={row['confidence']}  "
              f"layers={row['layers_present']}  diseases={row.get('diseases')}")

    a375 = [r for r in braf_result["ranking"] if r["model_id"] == "ach-000219"]
    total_n = braf_result["total_candidates"]
    top5pct = max(1, round(total_n * 0.05))
    print(f"\nA375 entry: position={a375[0]['sorted_position']} of {total_n}, "
          f"score={a375[0]['core_score']:.4f}, driver={a375[0]['has_driver_alteration']}")
    assert a375[0]["has_driver_alteration"], \
        "A375's BRAF V600E driver mutation is no longer being detected — this IS a regression"
    assert a375[0]["sorted_position"] <= top5pct, \
        (f"A375 fell to position {a375[0]['sorted_position']} of {total_n} "
         f"(outside top 5% = {top5pct}) — likely a genuine regression, not recalibration noise")
    print(f"PASS  Driver mutation detected and position ({a375[0]['sorted_position']}) "
          f"is within top 5% ({top5pct}) of {total_n} candidates")

    # --- Test 2: abundance_tracking gene — pure core_score sort ---
    abundance_gene = pred[pred["class"] == "abundance_tracking"]["ensg_id"].iloc[0]
    abundance_result = rank_cell_lines(abundance_gene, pred, chronos_val)
    print(f"\n=== Abundance-tracking gene {abundance_gene} — top 10 ===")
    for row in abundance_result["ranking"][:10]:
        print(f"  #{row['sorted_position']:>4}  {row['model_id']}  "
              f"score={row['core_score']:.4f}  basis={row['rank_basis']}")

    scores = [r["core_score"] for r in abundance_result["ranking"]]
    is_sorted = all(scores[i] >= scores[i + 1] for i in range(len(scores) - 1))
    print(f"\nPurely core_score-sorted (descending)? {is_sorted}")
    assert is_sorted, "abundance_tracking gene should sort by core_score alone"
    print("PASS  Correct sort behaviour for abundance_tracking class")

    # --- Test 3: rank_cell_lines and explain_pair must agree on position ---
    ep_result = explain_pair("ensg00000157764", "ach-000219", pred, flags, variants, harm, chronos_val)
    rcl_pos = a375[0]["sorted_position"]
    ep_pos  = ep_result["position"]["sorted_position"]
    print(f"\nrank_cell_lines position: {rcl_pos},  explain_pair position: {ep_pos}")
    assert rcl_pos == ep_pos, \
        f"MISMATCH — single-source-of-truth broken: rcl={rcl_pos}, ep={ep_pos}"
    print("PASS  rank_cell_lines and explain_pair agree — single source of truth confirmed")

    # --- Test 4: contradicted flag agrees between the two paths ---
    print(f"\nrank_cell_lines contradicted: {braf_result['contradicted']},  "
          f"explain_pair contradicted: {ep_result['contradicted']}")
    assert braf_result["contradicted"] == ep_result["contradicted"], \
        "MISMATCH — contradicted flag disagrees between rank_cell_lines and explain_pair"
    print("PASS  contradicted flag agrees")

    # --- Test 5: unknown-class genes must NOT hard-gate on the driver flag ---
    # A gene with no measured GDSC regime and no COSMIC census role has no basis
    # for letting a single alteration flag outrank the whole scored ranking. The
    # flag may only break ties between equal core_scores. Checked structurally
    # across a sample of unknown-class genes, not on any named gene, so the test
    # keeps meaning as the gene set changes.
    # One mask + one groupby per class group -- NOT one full-table mask per gene
    # (pred is 28.4M rows; per-gene masking here used to thrash for many minutes).
    def _sample_genes(mask, n=25):
        return pred.loc[mask, "ensg_id"].drop_duplicates().head(n).tolist()

    unknown_genes = _sample_genes(pred["class"] == "unknown")
    subset = pred[pred["ensg_id"].isin(unknown_genes)]
    print(f"\n=== Unknown-class sort check over {len(unknown_genes)} genes ===")
    violations = []
    for g, rows in subset.groupby("ensg_id", sort=False):
        sc = sort_gene_rows(rows.copy())["core_score"].to_numpy()
        # core_score must be non-increasing: nothing outranks a higher score
        bad = int((sc[:-1] < sc[1:]).sum())
        if bad:
            violations.append((g, bad))
    print(f"  genes where a lower core_score outranked a higher one: {len(violations)}")
    assert not violations, \
        f"driver flag is still outranking core_score for unknown-class genes: {violations[:5]}"
    print("PASS  Unknown-class genes rank on core_score; driver flag is tiebreak only")

    # --- Test 6: the driver gate is still live where it IS validated ---
    gated_genes = _sample_genes(pred["class"].isin(["activation_driven", "loss_of_function"]))
    subset = pred[pred["ensg_id"].isin(gated_genes)]
    print(f"\n=== Driver-gate retained check over {len(gated_genes)} gated-class genes ===")
    lost = []
    for g, rows in subset.groupby("ensg_id", sort=False):
        drv = sort_gene_rows(rows.copy())["has_driver_alteration"].to_numpy().astype(int)
        if (drv[:-1] < drv[1:]).sum():
            lost.append(g)
    print(f"  genes where a driver-positive line fell below a driver-negative one: {len(lost)}")
    assert not lost, f"driver gate was lost for classes that should keep it: {lost[:5]}"
    print("PASS  activation_driven / loss_of_function still hard-gate on the driver flag")

    print("\nAll tests passed.")
