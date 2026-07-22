import sys
import pandas as pd

sys.path.append("src/pipeline")
from explain_pair import sort_gene_rows, load_data


def rank_cell_lines(ensg_id, pred=None):
    """
    Returns the ranked list of cell lines for a gene, using the same
    query-time sort logic validated in explain_pair (not stored row order).
    """
    if pred is None:
        pred, _, _, _, _ = load_data()

    gene_rows = pred[pred["ensg_id"] == ensg_id].copy()
    if gene_rows.empty:
        return {"error": f"No predictions found for gene {ensg_id}"}

    sorted_rows = sort_gene_rows(gene_rows)

    result = sorted_rows[[
        "sorted_position", "model_id", "core_score", "n_layers",
        "stratum_rank", "class", "rank_basis", "has_driver_alteration",
        "confidence", "confidence_reason"
    ]].rename(columns={"stratum_rank": "stratum_rank_percentile"})

    return {
        "gene": ensg_id,
        "class": sorted_rows["class"].iloc[0],
        "total_candidates": len(sorted_rows),
        "ranking": result.to_dict(orient="records"),
    }


if __name__ == "__main__":
    from explain_pair import explain_pair

    pred, flags, variants, harm, gene_lookup = load_data()

    # --- Test 1: BRAF (activation_driven) — A375 must land at position 13 ---
    braf_result = rank_cell_lines("ENSG00000157764", pred)
    print("=== BRAF ranking — top 15 ===")
    for row in braf_result["ranking"][:15]:
        print(f"  #{row['sorted_position']:>4}  {row['model_id']}  "
              f"score={row['core_score']:.4f}  driver={row['has_driver_alteration']}  "
              f"basis={row['rank_basis']}  conf={row['confidence']}")

    a375 = [r for r in braf_result["ranking"] if r["model_id"] == "ACH-000219"]
    print(f"\nA375 entry: position={a375[0]['sorted_position']}, "
          f"score={a375[0]['core_score']:.4f}, driver={a375[0]['has_driver_alteration']}")
    assert a375[0]["sorted_position"] == 13, \
        f"Expected position 13, got {a375[0]['sorted_position']}"
    print("PASS  Position matches explain_pair's validated result")

    # --- Test 2: abundance_tracking gene — pure core_score sort ---
    abundance_gene = pred[pred["class"] == "abundance_tracking"]["ensg_id"].iloc[0]
    abundance_result = rank_cell_lines(abundance_gene, pred)
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
    ep_result = explain_pair("ENSG00000157764", "ACH-000219", pred, flags, variants, harm)
    rcl_pos = a375[0]["sorted_position"]
    ep_pos  = ep_result["position"]["sorted_position"]
    print(f"\nrank_cell_lines position: {rcl_pos},  explain_pair position: {ep_pos}")
    assert rcl_pos == ep_pos, \
        f"MISMATCH — single-source-of-truth broken: rcl={rcl_pos}, ep={ep_pos}"
    print("PASS  rank_cell_lines and explain_pair agree — single source of truth confirmed")

    print("\nAll tests passed.")
