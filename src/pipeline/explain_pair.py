import pandas as pd
import numpy as np
import json


def load_data(base="src/pipeline/outputs"):
    pred     = pd.read_parquet(f"{base}/predictions_with_confidence.parquet")
    flags    = pd.read_parquet(f"{base}/flags_with_driver.parquet")
    variants = pd.read_parquet("cleaned_track_data/mutations_variant_detail.parquet")
    harm     = pd.read_parquet(f"{base}/harmonised_enriched.parquet")
    gene_lookup = pd.read_parquet("reference/gene_lookup.parquet")
    return pred, flags, variants, harm, gene_lookup


def _load_data_compat(base="src/pipeline/outputs"):
    """Back-compat wrapper that returns the original 4-tuple."""
    pred, flags, variants, harm, _ = load_data(base)
    return pred, flags, variants, harm


def resolve_gene(gene_input, gene_lookup):
    """
    Accept an ENSG ID or a gene symbol (case-insensitive).
    Returns (ensg_id, hgnc_symbol) or raises ValueError with a clear message.
    """
    g = gene_input.strip()
    if g.upper().startswith("ENSG"):
        row = gene_lookup[gene_lookup["ensg_id"] == g.upper()]
        symbol = row["hgnc_symbol"].iloc[0] if not row.empty else g
        return g.upper(), symbol

    # Try exact symbol match first, then alias
    row = gene_lookup[gene_lookup["hgnc_symbol"].str.upper() == g.upper()]
    if row.empty:
        # Search in alias_symbols (pipe-separated strings); avoid regex groups warning
        pattern = r"(?:^|[|])" + g.upper() + r"(?:$|[|])"
        alias_mask = gene_lookup["alias_symbols"].fillna("").str.upper().str.contains(
            pattern, regex=True
        )
        row = gene_lookup[alias_mask]
    if row.empty:
        raise ValueError(f"Gene '{gene_input}' not found in gene_lookup (tried symbol and aliases).")
    ensg = row["ensg_id"].iloc[0]
    symbol = row["hgnc_symbol"].iloc[0]
    return ensg, symbol
    # normalise model_id case
    for df in (pred, flags, variants):
        df["model_id"] = df["model_id"].str.upper()
    return pred, flags, variants, harm


def sort_gene_rows(gene_rows):
    """Re-derive the correct ranked order at query time — never trust stored row order."""
    gene_class = gene_rows["class"].iloc[0]
    if gene_class == "abundance_tracking":
        sorted_rows = gene_rows.sort_values("core_score", ascending=False).reset_index(drop=True)
    else:
        sorted_rows = gene_rows.sort_values(
            ["has_driver_alteration", "core_score"], ascending=[False, False]
        ).reset_index(drop=True)
    sorted_rows["sorted_position"] = sorted_rows.index + 1
    return sorted_rows


def explain_pair(ensg_id, model_id, pred, flags, variants, harm):
    model_id_upper = model_id.upper()

    gene_rows = pred[pred["ensg_id"] == ensg_id].copy()
    if gene_rows.empty:
        return {"error": f"No predictions found for gene {ensg_id}"}

    sorted_rows = sort_gene_rows(gene_rows)
    total_n = len(sorted_rows)

    row = sorted_rows[sorted_rows["model_id"] == model_id_upper]
    if row.empty:
        return {"error": f"No prediction found for pair ({ensg_id}, {model_id_upper})"}
    row = row.iloc[0]

    # Cell-line metadata (harmonised_enriched uses lowercase model_id)
    meta_row = harm[harm["model_id"] == model_id.lower()]
    meta = {
        "canonical_name": meta_row["canonical_name"].iloc[0] if not meta_row.empty else None,
        "gap_class":      meta_row["gap_class"].iloc[0]      if not meta_row.empty else None,
    }

    # RANKING tier (expression + proteomics via core_score/n_layers)
    ranking_tier = {
        "core_score":              float(row["core_score"]),
        "n_layers":                int(row["n_layers"]),
        "stratum_rank_percentile": float(row["stratum_rank"]),
        "tag": "RANKING",
        "note": ("core_score reflects expression and/or proteomics presence per n_layers; "
                 "individual expression/proteomics values not yet wired into explain_pair (TODO v1)."),
    }

    # CONFIDENCE MODIFIER tier (mutation, fusion)
    flag_row = flags[(flags["ensg_id"] == ensg_id) & (flags["model_id"] == model_id_upper)]
    modifier_tier = {}
    absent = []

    if not flag_row.empty:
        fr = flag_row.iloc[0]

        if pd.notna(fr.get("p_mutation")) and fr["p_mutation"] > 0:
            mut_entry = {
                "p_mutation":  float(fr["p_mutation"]),
                "driver_flag": bool(fr["mut_driver"]),
                "tag": "CONFIDENCE MODIFIER — does not affect stratum_rank",
            }
            vrows = variants[
                (variants["ensg_id"] == ensg_id) & (variants["model_id"] == model_id_upper)
            ]
            if not vrows.empty:
                v = vrows.iloc[0]
                mut_entry["variant_detail"] = {
                    "protein_change": v.get("proteinchange"),
                    "hotspot":        bool(v["hotspot"]) if pd.notna(v.get("hotspot")) else None,
                    "vep_impact":     v.get("vepimpact"),
                }
            modifier_tier["mutation"] = mut_entry
        else:
            absent.append("mutation")

        if pd.notna(fr.get("p_fusion")) and fr["p_fusion"] > 0:
            modifier_tier["fusion"] = {
                "p_fusion":    float(fr["p_fusion"]),
                "driver_flag": bool(fr["fusion_driver"]),
                "tag": "CONFIDENCE MODIFIER — does not affect stratum_rank",
            }
        else:
            absent.append("fusion")
    else:
        absent.extend(["mutation", "fusion"])

    absent.append("methylation [not wired in v0]")

    # Self-check
    inconsistencies = []
    expected_basis = (
        "core_score" if row["class"] == "abundance_tracking"
        else ("driver_flag+score" if row["has_driver_alteration"] else "score_only")
    )
    if row["rank_basis"] != expected_basis:
        inconsistencies.append(
            f"rank_basis mismatch: stored='{row['rank_basis']}', expected='{expected_basis}' "
            f"given class={row['class']} and has_driver_alteration={row['has_driver_alteration']}"
        )

    # Check that driver pairs actually rank above non-driver pairs of the same gene
    if row["rank_basis"] == "driver_flag+score":
        pos = int(row["sorted_position"])
        non_driver_above = sorted_rows[
            (sorted_rows["sorted_position"] < pos) & (~sorted_rows["has_driver_alteration"])
        ]
        if not non_driver_above.empty:
            inconsistencies.append(
                f"sort inconsistency: {len(non_driver_above)} non-driver row(s) rank above "
                f"this driver pair at position {pos}"
            )

    return {
        "gene":       ensg_id,
        "cell_line":  {"model_id": model_id_upper, **meta},
        "position": {
            "sorted_position":         int(row["sorted_position"]),
            "total_candidates":        total_n,
            "stratum_rank_percentile": float(row["stratum_rank"]),
            "note": (
                "sorted_position is query-time re-derived rank "
                "([has_driver_alteration, core_score] desc for activation_driven; "
                "core_score desc for abundance_tracking). "
                "stratum_rank_percentile is a separate fact: percentile of core_score "
                "within the (gene, n_layers) stratum. These are NOT the same measurement."
            ),
        },
        "class":              row["class"],
        "rank_basis":         row["rank_basis"],
        "confidence":         row["confidence"],
        "confidence_reason":  row["confidence_reason"],
        "ranking_tier":       ranking_tier,
        "confidence_modifier_tier": modifier_tier,
        "absent_layers":      absent,
        "inconsistencies_flagged": inconsistencies,
    }


if __name__ == "__main__":
    pred, flags, variants, harm = load_data()
    result = explain_pair("ENSG00000157764", "ACH-000219", pred, flags, variants, harm)
    print(json.dumps(result, indent=2, default=str))
