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
from evidence_state import gene_verdict, pair_verdict, load_dispersion, render
from multi_gene import co_select, format_result
from metadata import (load_cell_lines, load_genes, cell_line_meta, gene_meta,
                      format_block, modality_line)

_DISPERSION = None


def print_ranking(gene_input, pred, gene_lookup, chronos_val):
    try:
        ensg_id, symbol = resolve_gene(gene_input, gene_lookup)
    except ValueError as e:
        print(e)
        return

    result = rank_cell_lines(ensg_id, pred, chronos_val)
    if "error" in result:
        # Gene known to the lookup but never scored. Say so as a state, in the
        # same voice as every other "cannot answer" case, rather than as an
        # incidental error string.
        print()
        print(render(gene_verdict(ensg_id, symbol, 0)))
        print("  Predictions cover the 19,176 genes with at least one scored omics "
              "layer.\n  This gene is in the reference lookup but no layer measures it.")
        return

    print(f"\nGene: {symbol} ({ensg_id}) | class: {result['class']} | "
          f"{result['total_candidates']} candidates")

    # The verdict comes BEFORE the table and can suppress it entirely. core_score
    # is a within-gene percentile, so a gene that barely varies still yields a
    # full 0-1 spread and a top-10 that reads as confidently as a truly
    # differential gene's. A footnote under a confident-looking table does not
    # correct that impression; a headline that refuses to rank does.
    # States and thresholds: evidence_state.py.
    global _DISPERSION
    if _DISPERSION is None:
        _DISPERSION = load_dispersion()
    disp = (_DISPERSION.loc[[ensg_id]]
            if _DISPERSION is not None and ensg_id in _DISPERSION.index else None)

    v = gene_verdict(ensg_id, symbol, result["total_candidates"], disp,
                     contradicted=result["contradicted"],
                     gene_class=result.get("class"))
    text = render(v)
    if text:
        print(text)
    if not v.show_ranking:
        return
    if v.qualify_ranking:
        print("  (unranked)")
    print()
    cells = load_cell_lines()
    print(f"{'Pos':<4}{'model_id':<13}{'name':<14}{'score':<8}{'lay':<5}{'conf':<10}"
          f"{'mods':<6}{'lineage':<14}{'disease'}")
    for row in result["ranking"][:10]:
        cm = cell_line_meta(row["model_id"], cells)
        print(f"{row['sorted_position']:<4}{row['model_id']:<13}"
              f"{str(cm.get('name') or '-')[:13]:<14}{row['core_score']:<8.4f}"
              f"{row['n_layers']:<5}{row['confidence'][:9]:<10}"
              f"{str(cm.get('_n_modalities') or '-'):<6}"
              f"{str(cm.get('lineage') or '-')[:13]:<14}"
              f"{str(cm.get('disease') or cm.get('cancer type (Sanger)') or '-')[:34]}")
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
        # Absence is reported as absence, never as a weak result (Design Record §5:
        # "missing != against").
        print()
        print(render(pair_verdict(symbol, model_id, row=None)))
        return

    print(f"\nGene: {symbol} ({result['gene']}) | Cell line: {result['cell_line']['model_id']} "
          f"({result['cell_line']['canonical_name']})")
    print(f"Position: #{result['position']['sorted_position']} of {result['position']['total_candidates']}"
          f"  |  Stratum percentile: {result['position']['stratum_rank_percentile']:.4f}")
    print(f"Confidence: {result['confidence']}  |  rank_basis: {result['rank_basis']}")
    print(f"Reason: {result['confidence_reason']}")

    gmeta = gene_meta(ensg_id, load_genes())
    if gmeta:
        print("\nGENE METADATA")
        print(format_block(gmeta))
    cmeta = cell_line_meta(result["cell_line"]["model_id"], load_cell_lines())
    if cmeta:
        print("\nCELL LINE METADATA")
        print(format_block(cmeta))
        print(modality_line(cmeta))
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


def print_co_selection(gene_inputs, pred, gene_lookup, chronos_val):
    """Cell lines that satisfy SEVERAL genes jointly — weakest-link bound."""
    global _DISPERSION
    if _DISPERSION is None:
        _DISPERSION = load_dispersion()
    res = co_select(gene_inputs, pred, gene_lookup, resolve_gene,
                    dispersion=_DISPERSION, chronos_val=chronos_val)
    print(f"\nJOINT QUERY: {' + '.join(g.upper() for g in gene_inputs)}")
    print(format_result(res))


def main():
    pred, flags, variants, harm, gene_lookup, chronos_val = load_data()

    print("GeneTraceAI v0 — CLI")
    print("Commands: 'gene <ENSG_ID|symbol>' | 'pair <ENSG_ID|symbol> <MODEL_ID>'")
    print("          'genes <A> <B> [C ...]'  — cell lines suiting ALL of them | 'quit'\n")

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
        elif parts[0] == "genes" and len(parts) >= 3:
            print_co_selection(parts[1:], pred, gene_lookup, chronos_val)
        else:
            print("Usage: gene <ENSG_ID|symbol>  |  pair <ENSG_ID|symbol> <MODEL_ID>"
                  "  |  genes <A> <B> [C ...]")


if __name__ == "__main__":
    main()
