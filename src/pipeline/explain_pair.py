import pandas as pd
import numpy as np
import json
import sys
sys.path.insert(0, "src/pipeline")
from evidence_ledger import (
    rna_covered_model_ids, derive_layers_present, derive_driver_alteration,
    derive_cna_role_consistent, derive_essentiality_check, derive_census_role,
    derive_validation_source, derive_signal_spread, FIELD_TIERS,
)

_RNA_COVERED = None  # lazily built, cached module-level (cheap: ~1495 model_ids)


def load_data(base="final_pipeline/outputs"):
    pred     = pd.read_parquet(f"{base}/predictions_with_confidence.parquet")
    flags    = pd.read_parquet(f"{base}/flags_with_driver.parquet")
    variants = pd.read_parquet("cleaned_track_data/mutations_variant_detail.parquet")
    harm     = pd.read_parquet(f"{base}/harmonised_enriched.parquet")
    gene_lookup = pd.read_parquet("reference/gene_lookup.parquet")
    try:
        chronos_val = pd.read_parquet(f"{base}/chronos_validation.parquet")
    except FileNotFoundError:
        chronos_val = None

    # Canonical key case is LOWERCASE, matching the harmonisation warehouse.
    # Normalising here rather than at each use makes this idempotent: it is a
    # no-op on files already written lowercase, and it silently upgrades the
    # legacy UPPERCASE artefacts (predictions_with_confidence, chronos_validation)
    # that predate the convention and cannot currently be rebuilt.
    for _df, _cols in [(pred, ("model_id", "ensg_id")), (flags, ("model_id", "ensg_id")),
                       (variants, ("model_id", "ensg_id")), (harm, ("model_id",)),
                       (chronos_val, ("model_id", "ensg_id")),
                       (gene_lookup, ("ensg_id", "hgnc_symbol", "alias_symbols"))]:
        for _c in _cols:
            if _c in _df.columns:
                _df[_c] = _df[_c].astype("string").str.lower()

    return pred, flags, variants, harm, gene_lookup, chronos_val


def _load_data_compat(base="src/pipeline/outputs"):
    """Back-compat wrapper that returns the original 4-tuple."""
    pred, flags, variants, harm, _, _ = load_data(base)
    return pred, flags, variants, harm


def resolve_gene(gene_input, gene_lookup):
    """
    Accept an ENSG ID or a gene symbol (case-insensitive).
    Returns (ensg_id, hgnc_symbol) or raises ValueError with a clear message.
    """
    g = gene_input.strip()
    if g.lower().startswith("ensg"):
        row = gene_lookup[gene_lookup["ensg_id"] == g.lower()]
        symbol = row["hgnc_symbol"].iloc[0] if not row.empty else g
        return g.lower(), symbol

    # Try exact symbol match first, then alias
    row = gene_lookup[gene_lookup["hgnc_symbol"].str.lower() == g.lower()]
    if row.empty:
        # Search in alias_symbols (pipe-separated strings); avoid regex groups warning
        pattern = r"(?:^|[|])" + g.lower() + r"(?:$|[|])"
        alias_mask = gene_lookup["alias_symbols"].fillna("").str.lower().str.contains(
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
        df["model_id"] = df["model_id"].str.lower()
    return pred, flags, variants, harm


# Classes for which a driver alteration is allowed to OUTRANK core_score
# (a hard gate: every driver-positive line sorts above every driver-negative one).
#
# The gate is only defensible where the gene has an established role that makes
# the flag interpretable — a mutation in a census oncogene/TSG is mechanistic
# evidence, a mutation in a gene with no known role is just a mutation. Genes
# with class 'unknown' have neither a measured GDSC regime nor a COSMIC census
# role, so for them the flag drops to a tiebreaker *after* core_score.
#
# Why this matters: 18,514 of 19,176 scored genes are class 'unknown'. Under the
# previous rule they were hard-gated too, which promoted single flagged lines
# above the entire scored ranking — e.g. MUC1 ranked a uveal melanoma line at
# core_score 0.284 above pancreatic/cholangio/breast lines at 1.000-0.995, and
# ASGR1 ranked a bladder line at 0.574 above hepatoblastoma/HCC lines at 1.000.
# In both cases the score itself already had the biology right.
DRIVER_GATED_CLASSES = frozenset({"activation_driven", "loss_of_function"})


def driver_is_primary_sort_key(gene_class):
    """True where has_driver_alteration may outrank core_score for this class."""
    return gene_class in DRIVER_GATED_CLASSES


def sort_gene_rows(gene_rows):
    """Re-derive the correct ranked order at query time — never trust stored row order.

    core_score is already direction-adjusted for TSGs (inverted at build time),
    so descending sort is always correct here.
    """
    gene_class = gene_rows["class"].iloc[0]
    if gene_class == "abundance_tracking":
        # Empirically validated to rank on abundance alone; flag plays no part.
        sort_keys = ["core_score"]
    elif driver_is_primary_sort_key(gene_class):
        # Established role: driver alteration is mechanistic evidence and gates.
        # (core_score already inverted at build time for loss_of_function.)
        sort_keys = ["has_driver_alteration", "core_score"]
    else:
        # No established role: the flag cannot outrank measured abundance, but it
        # is still evidence, so it breaks ties between equal scores.
        sort_keys = ["core_score", "has_driver_alteration"]

    sorted_rows = gene_rows.sort_values(
        sort_keys, ascending=[False] * len(sort_keys)
    ).reset_index(drop=True)
    sorted_rows["sorted_position"] = sorted_rows.index + 1
    return sorted_rows


def explain_pair(ensg_id, model_id, pred, flags, variants, harm, chronos_val=None):
    global _RNA_COVERED
    model_id_lc = model_id.lower()

    gene_rows = pred[pred["ensg_id"] == ensg_id].copy()
    if gene_rows.empty:
        return {"error": f"No predictions found for gene {ensg_id}"}

    sorted_rows = sort_gene_rows(gene_rows)
    total_n = len(sorted_rows)

    row = sorted_rows[sorted_rows["model_id"] == model_id_lc]
    if row.empty:
        return {"error": f"No prediction found for pair ({ensg_id}, {model_id_lc})"}
    row = row.iloc[0]

    # Cell-line metadata (every key in the chain is lowercase)
    meta_row = harm[harm["model_id"] == model_id_lc]
    meta = {
        "canonical_name": meta_row["canonical_name"].iloc[0] if not meta_row.empty else None,
        "gap_class":      meta_row["gap_class"].iloc[0]      if not meta_row.empty else None,
    }

    gene_role = row.get("gene_role", "unknown") if "gene_role" in row.index else "unknown"
    if gene_role == "tsg":
        direction = "inverted (low abundance = high relevance) — tumour suppressor"
    elif gene_role == "oncogene":
        direction = "normal (high abundance = high relevance) — oncogene"
    elif gene_role == "both":
        direction = "normal (high abundance = high relevance) — oncogene/TSG dual role"
    else:
        direction = "normal (high abundance = high relevance) — role unknown"

    # RANKING tier (expression + proteomics via core_score/n_layers)
    ranking_tier = {
        "core_score":              float(row["core_score"]),
        "n_layers":                int(row["n_layers"]),
        "stratum_rank_percentile": float(row["stratum_rank"]),
        "gene_role":               gene_role,
        "direction":               direction,
        "tag": "RANKING",
        "note": ("core_score reflects expression and/or proteomics presence per n_layers; "
                 "individual expression/proteomics values not yet wired into explain_pair (TODO v1)."),
    }

    # CONFIDENCE MODIFIER tier (mutation, fusion)
    flag_row = flags[(flags["ensg_id"] == ensg_id) & (flags["model_id"] == model_id_lc)]
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
                (variants["ensg_id"] == ensg_id) & (variants["model_id"] == model_id_lc)
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

    # ── Evidence ledger — mechanical, not hand-written (see evidence_ledger.py) ──
    # Same derivation functions as build_evidence_ledger.py's bulk pass, so a
    # single pair's explanation and the 28.4M-row bulk table can never drift
    # apart. Fields are NOT additive -- each is an independent fact, tagged
    # with the tier that says whether it drove the ranking or only the label.
    if _RNA_COVERED is None:
        _RNA_COVERED = rna_covered_model_ids()

    fr = flag_row.iloc[0] if not flag_row.empty else pd.Series(dtype=object)
    chronos_check = row.get("chronos_check", "n/a")
    rho = pval = None
    if chronos_val is not None:
        cv_row = chronos_val[chronos_val["ensg_id"] == ensg_id]
        if not cv_row.empty:
            rho, pval = float(cv_row.iloc[0]["rho"]), float(cv_row.iloc[0]["pval"])
    contradicted = chronos_check == "inverted"

    evidence_ledger = {
        "layers_present":      derive_layers_present(int(row["n_layers"]), model_id_lc, _RNA_COVERED),
        "n_layers":             int(row["n_layers"]),
        "signal_spread":        derive_signal_spread(
            row.get("signal_spread"), row.get("gene_top_vs_median_fold")),
        "driver_alteration":    derive_driver_alteration(fr.get("mut_driver"), fr.get("fusion_driver")),
        "cna_role_consistent":  derive_cna_role_consistent(
            gene_role, row["class"], bool(row.get("has_cna_alteration", False)),
            bool(fr.get("has_cna_context", False)) if not flag_row.empty else False,
        ),
        "essentiality_check":   derive_essentiality_check(chronos_check, rho, pval),
        "census_role":          derive_census_role(gene_role),
        "validation_source":    derive_validation_source(row.get("regime_source")),
        "field_tiers":          FIELD_TIERS,
    }

    # Self-check
    inconsistencies = []
    if row["class"] == "loss_of_function":
        expected_basis = "lof_flag+inverted" if row["has_driver_alteration"] else "lof_inverted"
    elif row["class"] == "unknown":
        expected_basis = "score+driver_tiebreak"
    elif row["class"] == "abundance_tracking":
        expected_basis = "core_score"
    elif row["has_driver_alteration"]:
        expected_basis = "driver_flag+score"
    else:
        expected_basis = "score_only"
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
        "gene_role":  gene_role,
        "cell_line":  {"model_id": model_id_lc, **meta},
        "position": {
            "sorted_position":         int(row["sorted_position"]),
            "total_candidates":        total_n,
            "stratum_rank_percentile": float(row["stratum_rank"]),
            "note": (
                "sorted_position is query-time re-derived rank "
                "([has_driver_alteration, core_score] desc for activation_driven and "
                "loss_of_function, where the census/GDSC role makes the driver flag "
                "mechanistic evidence; core_score desc for abundance_tracking; "
                "[core_score, has_driver_alteration] desc for unknown-class genes, "
                "where the flag breaks ties but cannot outrank measured abundance). "
                "stratum_rank_percentile is a separate fact: percentile of core_score "
                "within the (gene, n_layers) stratum. These are NOT the same measurement."
            ),
        },
        "class":              row["class"],
        "rank_basis":         row["rank_basis"],
        "confidence":         row["confidence"],
        "confidence_reason":  row["confidence_reason"],
        "contradicted":       contradicted,
        "evidence_ledger":    evidence_ledger,
        "ranking_tier":       ranking_tier,
        "confidence_modifier_tier": modifier_tier,
        "absent_layers":      absent,
        "inconsistencies_flagged": inconsistencies,
    }


if __name__ == "__main__":
    pred, flags, variants, harm, gene_lookup, chronos_val = load_data()
    result = explain_pair("ensg00000157764", "ach-000219", pred, flags, variants, harm, chronos_val)
    print(json.dumps(result, indent=2, default=str))
