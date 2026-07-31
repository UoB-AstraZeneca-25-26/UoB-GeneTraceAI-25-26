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


def print_ranking(gene_input, pred, gene_lookup, chronos_val):
    try:
        ensg_id, symbol = resolve_gene(gene_input, gene_lookup)
    except ValueError as e:
        print(e)
        return

    result = rank_cell_lines(ensg_id, pred, chronos_val)
    if "error" in result:
        # Gene known to lookup but not in scored set — give a useful message
        print(f"{symbol} ({ensg_id}): {result['error']}")
        print("Note: only the 141 curated GDSC-targeted genes have predictions in v0.")
        return

    warn = "  *** CONTRADICTED: abundance ranking correlates NEGATIVELY with real CRISPR essentiality for this gene ***" \
        if result["contradicted"] else ""
    print(f"\nGene: {symbol} ({ensg_id}) | class: {result['class']} | {result['total_candidates']} candidates{warn}\n")
    print(f"{'Pos':<5}{'model_id':<14}{'core_score':<12}{'n_layers':<10}{'stratum_pct':<13}"
          f"{'confidence':<12}{'layers_present':<16}{'diseases'}")
    for row in result["ranking"][:10]:
        print(f"{row['sorted_position']:<5}{row['model_id']:<14}{row['core_score']:<12.4f}"
              f"{row['n_layers']:<10}{row['stratum_rank_percentile']:<13.4f}"
              f"{row['confidence']:<12}{row['layers_present']:<16}{row.get('diseases') or ''}")
    if result["total_candidates"] > 10:
        print(f"... ({result['total_candidates'] - 10} more)")


def print_explanation(gene_input, model_id, pred, flags, variants, harm, gene_lookup, chronos_val):
    try:
        ensg_id, symbol = resolve_gene(gene_input, gene_lookup)
    except ValueError as e:
        print(e)
        return

    result = explain_pair(ensg_id, model_id, pred, flags, variants, harm, chronos_val)
    if "error" in result:
        print(result["error"])
        return

    print(f"\nGene: {symbol} ({result['gene']}) | Cell line: {result['cell_line']['model_id']} "
          f"({result['cell_line']['canonical_name']})")
    print(f"Position: #{result['position']['sorted_position']} of {result['position']['total_candidates']}"
          f"  |  Stratum percentile: {result['position']['stratum_rank_percentile']:.4f}")
    print(f"Confidence: {result['confidence']}  |  rank_basis: {result['rank_basis']}")
    print(f"Reason: {result['confidence_reason']}")
    if result["contradicted"]:
        print("*** CONTRADICTED: abundance ranking correlates NEGATIVELY with real CRISPR "
              "essentiality for this gene — treat with particular caution ***")
    print()

    # Evidence ledger — mechanical, sourced from evidence_ledger.py, same fields
    # and same tiers as the bulk 28.4M-row table. Ranking-tier fields first
    # (what drove the score), then confidence-modifier fields (what only
    # touched the label) — not additive, each field is an independent fact.
    led = result["evidence_ledger"]
    tiers = led["field_tiers"]
    print("EVIDENCE LEDGER")
    print("  -- RANKING tier (drove core_score / stratum_rank) --")
    for field, tier in tiers.items():
        if tier == "RANKING":
            print(f"    {field}: {led[field]}")
    print("  -- CONFIDENCE MODIFIER tier (affects rank_basis/confidence label only) --")
    for field, tier in tiers.items():
        if tier == "CONFIDENCE_MODIFIER":
            marker = "  [CONTRADICTS]" if field == "essentiality_check" and result["contradicted"] else ""
            print(f"    {field}: {led[field]}{marker}")

    print("\nDETAIL (mutation/fusion, richer than the ledger's terse driver_alteration string):")
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
    pred, flags, variants, harm, gene_lookup, chronos_val = load_data()

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
            print_ranking(parts[1], pred, gene_lookup, chronos_val)
        elif parts[0] == "pair" and len(parts) == 3:
            print_explanation(parts[1], parts[2], pred, flags, variants, harm, gene_lookup, chronos_val)
        else:
            print("Usage: gene <ENSG_ID|symbol>  |  pair <ENSG_ID|symbol> <MODEL_ID>")


if __name__ == "__main__":
    main()
