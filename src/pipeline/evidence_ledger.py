"""
evidence_ledger.py
-------------------
Shared derivation logic for the evidence ledger — a projection of columns
that already exist (Stage 1-6 outputs), not a new model. Used by both
explain_pair.py (single-pair, live query) and build_evidence_ledger.py
(bulk, all 28.4M rows) so the two never drift apart.

Design rules (per review):
  - Fields are NOT additive. There is no combined "evidence score" anywhere
    in this module — each field is an independent fact, tagged with the
    tier that says whether it drove the ranking or only touched the label.
  - The 5 existing confidence labels (high/moderate/low/prior/unknown) are
    kept as-is, as a sort key. This module explains them, it doesn't
    replace them.
  - `contradicted` (chronos_check == 'inverted', 1,258 genes) is surfaced
    as its own boolean field, not buried inside a text string, so it can
    be filtered/sorted/counted as prominently as "confirmed" evidence.

Tiers (matches the RANKING / CONFIDENCE MODIFIER split explain_pair.py
already used):
  RANKING            -- drove core_score / stratum_rank
  CONFIDENCE_MODIFIER -- affects rank_basis / sorted_position / confidence
                         label only, never stratum_rank
"""
import pandas as pd
from pathlib import Path

REF = Path("reference")

FIELD_TIERS = {
    "layers_present":      "RANKING",
    "n_layers":             "RANKING",
    "signal_spread":        "RANKING",
    "driver_alteration":    "CONFIDENCE_MODIFIER",
    "cna_role_consistent":  "CONFIDENCE_MODIFIER",
    "essentiality_check":   "CONFIDENCE_MODIFIER",
    "census_role":          "CONFIDENCE_MODIFIER",
    "validation_source":    "CONFIDENCE_MODIFIER",
}
# Column display order: ranking tier first, then confidence-modifier tier —
# what drove the ranking, then what only touched the label.
LEDGER_FIELD_ORDER = [f for f, t in FIELD_TIERS.items() if t == "RANKING"] + \
                     [f for f, t in FIELD_TIERS.items() if t == "CONFIDENCE_MODIFIER"]


def rna_covered_model_ids():
    """Cell lines with an RNA profile at all. The universe-filtered
    expression matrix has 0% missing cells within this set (verified),
    so for n_layers==1 rows, RNA-profile membership alone disambiguates
    layers_present without an expensive per-pair presence join."""
    prof = pd.read_parquet(REF / "depmap_profiles.parquet")
    rna = prof[prof["datatype"] == "rna"]
    return set(rna["modelid"].str.lower())


def derive_layers_present(n_layers, model_id, rna_covered):
    if n_layers == 2:
        return "RNA + Protein"
    if n_layers == 1:
        return "RNA only" if model_id in rna_covered else "Protein only"
    return "none"  # shouldn't occur -- rows with 0 layers are dropped upstream


def derive_signal_spread(signal_spread, top_vs_median_fold=None):
    """How far the top of this gene's abundance distribution sits above a typical
    cell line, measured on the RAW pre-percentile values (02_core_score Cell 7c).

    core_score is a within-gene percentile, so every gene is rescaled onto the
    same uniform 0-1 range. That makes a near-flat gene's top-10 look exactly as
    confident as a sharply tissue-restricted gene's -- GAPDH's best line scores
    0.9993 against ERBB2's 1.0000. This field is the missing context: the
    ordering is still real, but 'narrow' says the differences behind it are
    small enough that acting on the top of the list may not be meaningful.

    Gene-level fact, identical for every cell line of the gene. It does not
    enter core_score and does not reorder anything.
    """
    if signal_spread is None or (isinstance(signal_spread, float) and pd.isna(signal_spread)):
        signal_spread = "unmeasured"
    if signal_spread in ("unmeasured", ""):
        return "unmeasured (no RNA layer for this gene)"

    fold = ("n/a" if top_vs_median_fold is None or pd.isna(top_vs_median_fold)
            else f"{float(top_vs_median_fold):.1f}x")
    if signal_spread == "narrow":
        return (f"NARROW ({fold} top-vs-median) — abundance is close to flat across "
                "cell lines; ranking is valid but the differences behind it are small")
    if signal_spread == "wide":
        return (f"wide ({fold} top-vs-median) — abundance strongly differentiates "
                "cell lines")
    return f"moderate ({fold} top-vs-median)"


def derive_driver_alteration(mut_driver, fusion_driver):
    mut = bool(mut_driver) if pd.notna(mut_driver) else False
    fus = bool(fusion_driver) if pd.notna(fusion_driver) else False
    if mut and fus:
        return "mutation + fusion driver detected"
    if mut:
        return "mutation driver detected"
    if fus:
        return "fusion driver detected"
    return "none detected"


def derive_cna_role_consistent(gene_role, regime_class, has_cna_alteration, has_cna_context):
    no_role = pd.isna(gene_role) or gene_role in (None, "unknown", "")
    no_regime_role = regime_class not in ("loss_of_function", "activation_driven", "abundance_tracking")
    if no_role and no_regime_role:
        return "not applicable (no census role)"
    if has_cna_alteration:
        if regime_class == "loss_of_function" or gene_role in ("tsg",):
            return "consistent (deletion in TSG)"
        return "consistent (amplification in oncogene)"
    if has_cna_context:
        return "CNA present but inconsistent with gene role"
    return "not applicable (no CNA event in this line)"


def derive_essentiality_check(chronos_check, rho=None, pval=None):
    if chronos_check == "n/a":
        return "not applicable (class already empirically validated)"
    if chronos_check == "untested":
        return "untested (no coverage)"
    stat = (f" (rho={rho:.3f}, p={pval:.3g})"
            if rho is not None and pd.notna(rho) else "")
    if chronos_check == "validated":
        return f"validated{stat}"
    if chronos_check == "inverted":
        return f"CONTRADICTED{stat} — ranking may be backwards for this gene"
    # Weak bands: direction agrees but the effect is below the floor that would
    # justify a confidence upgrade or a contradiction warning (see
    # build_chronos_validation.py for why the floor is |rho| >= 0.30).
    if chronos_check == "weak_positive":
        return (f"weak agreement{stat} — direction matches real CRISPR essentiality "
                "but the effect is too small to treat as validation")
    if chronos_check == "weak_negative":
        return (f"weak disagreement{stat} — direction is opposite to real CRISPR "
                "essentiality but the effect is too small to treat as a contradiction")
    if chronos_check == "none":
        return "tested, no significant relationship"
    return "unknown"


def derive_census_role(gene_role):
    if pd.isna(gene_role) or gene_role in (None, "unknown", ""):
        return "none recorded"
    return gene_role


def derive_validation_source(regime_source):
    if regime_source == "measured":
        return "curated GDSC-validated (141-gene set)"
    if regime_source == "prior":
        return "COSMIC-prior (not empirically validated)"
    return "not in curated set"
