"""
Minimal CLI for GeneTraceAI agentic layer v0.
Two modes only, per spec — no dashboard creep:
  1. gene lookup  -> rank_cell_lines
  2. pair lookup  -> explain_pair

Accepts ENSG IDs or gene symbols (e.g. BRAF, SHANK3) for the gene argument.
"""
import sys
import json
sys.path.append("src/pipeline")
from explain_pair import explain_pair, load_data, resolve_gene
from rank_cell_lines import rank_cell_lines


def print_ranking(gene_input, pred, gene_lookup):
    try:
        ensg_id, symbol = resolve_gene(gene_input, gene_lookup)
    except ValueError as e:
        print(e)
        return

    result = rank_cell_lines(ensg_id, pred)
    if "error" in result:
        # Gene known to lookup but not in scored set — give a useful message
        print(f"{symbol} ({ensg_id}): {result['error']}")
        print("Note: only the 141 curated GDSC-targeted genes have predictions in v0.")
        return

    print(f"\nGene: {symbol} ({ensg_id}) | class: {result['class']} | {result['total_candidates']} candidates\n")
    print(f"{'Pos':<5}{'model_id':<14}{'core_score':<12}{'n_layers':<10}{'stratum_pct':<13}{'confidence':<12}{'rank_basis'}")
    for row in result["ranking"][:20]:
        print(f"{row['sorted_position']:<5}{row['model_id']:<14}{row['core_score']:<12.4f}"
              f"{row['n_layers']:<10}{row['stratum_rank_percentile']:<13.4f}"
              f"{row['confidence']:<12}{row['rank_basis']}")
    if result["total_candidates"] > 20:
        print(f"... ({result['total_candidates'] - 20} more)")


def print_explanation(gene_input, model_id, pred, flags, variants, harm, gene_lookup):
    try:
        ensg_id, symbol = resolve_gene(gene_input, gene_lookup)
    except ValueError as e:
        print(e)
        return

    result = explain_pair(ensg_id, model_id, pred, flags, variants, harm)
    if "error" in result:
        print(result["error"])
        return

    print(f"\nGene: {symbol} ({result['gene']}) | Cell line: {result['cell_line']['model_id']} "
          f"({result['cell_line']['canonical_name']})")
    print(f"Position: #{result['position']['sorted_position']} of {result['position']['total_candidates']}"
          f"  |  Stratum percentile: {result['position']['stratum_rank_percentile']:.4f}")
    print(f"Confidence: {result['confidence']}  |  rank_basis: {result['rank_basis']}")
    print(f"Reason: {result['confidence_reason']}\n")

    print("RANKING layers (drove score):")
    rt = result["ranking_tier"]
    print(f"  - core_score: {rt['core_score']:.4f}  n_layers: {rt['n_layers']}")
    print(f"  - {rt['note']}")

    print("\nCONFIDENCE MODIFIER layers (do not affect stratum_rank):")
    cmt = result["confidence_modifier_tier"]
    if cmt:
        for layer, detail in cmt.items():
            print(f"  - {layer}: {json.dumps(detail, default=str)}")
    else:
        print("  (none present)")

    print("\nABSENT layers:")
    for a in result["absent_layers"]:
        print(f"  - {a}")

    if result["inconsistencies_flagged"]:
        print("\nINCONSISTENCIES FLAGGED:")
        for i in result["inconsistencies_flagged"]:
            print(f"  - {i}")
    else:
        print("\nNo inconsistencies flagged.")


def main():
    pred, flags, variants, harm, gene_lookup = load_data()

    print("GeneTraceAI v0 — CLI")
    print("Commands: 'gene <ENSG_ID|symbol>' | 'pair <ENSG_ID|symbol> <MODEL_ID>' | 'quit'\n")

    while True:
        try:
            raw = input("> ").strip()
        except EOFError:
            break
        if not raw or raw.lower() == "quit":
            break

        parts = raw.split()
        if parts[0] == "gene" and len(parts) == 2:
            print_ranking(parts[1], pred, gene_lookup)
        elif parts[0] == "pair" and len(parts) == 3:
            print_explanation(parts[1], parts[2], pred, flags, variants, harm, gene_lookup)
        else:
            print("Usage: gene <ENSG_ID|symbol>  |  pair <ENSG_ID|symbol> <MODEL_ID>")


if __name__ == "__main__":
    main()
