"""
Scoring/confidence_tiers.py
-----------------------------
Assigns confidence tiers to (gene, cell-line) predictions using:
  - core_score percentile within gene
  - has_driver_alteration flag
  - gene regime (activation_driven vs loss_of_function vs abundance_tracking)

Tier definitions:
  HIGH    — top 20% score AND has driver alteration
  MEDIUM  — top 20% score AND NOT has driver alteration  (high score, no alteration support)
  CONTEXT — NOT top 20% AND has driver alteration        (alteration present, score low)
  LOW     — NOT top 20% AND NOT has driver alteration

Reads from:  final_pipeline/outputs/{core_score, flags_with_driver}.parquet
             reference/gene_lookup.parquet  (gene regime via gene_role)
Writes to:   final_pipeline/outputs/predictions_with_confidence.parquet
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from config import CORE_SCORE, FLAGS_DRIVER, GENE_LKP, PREDICTIONS

TOP_SCORE_PERCENTILE = 0.80   # within-gene cutoff for HIGH/MEDIUM tier


def run():
    print("=" * 70)
    print("Scoring — confidence tiers")
    print("=" * 70)

    for p in [FLAGS_DRIVER, CORE_SCORE]:
        if not p.exists():
            raise SystemExit(f"Missing input: {p}")

    flags = pd.read_parquet(FLAGS_DRIVER)
    glk   = pd.read_parquet(GENE_LKP, columns=["ensg_id","gene_role","hgnc_symbol"])

    flags = flags.merge(glk, on="ensg_id", how="left")

    # Within-gene score rank (for percentile cutoff)
    flags["score_rank_pct"] = (
        flags.groupby("ensg_id")["core_score"]
        .rank(pct=True, ascending=True, method="average")
    )

    top = flags["score_rank_pct"] >= TOP_SCORE_PERCENTILE
    drv = flags["has_driver_alteration"]
    is_lof = flags["is_lof_alteration"]

    # directional_top: for driver tiers (HIGH/CONTEXT) only.
    # LOF alterations (truncating mutation or deletion in TSG) rank by
    # LOW expression (bottom 20% = top rank); all other alterations rank by
    # HIGH expression (top 20%). MEDIUM/LOW tiers remain role-blind (use top).
    # gene_role=="both" is treated as GOF regardless of alteration type —
    # the "both" stratum has no empirical basis for directional separation.
    is_lof_for_direction = is_lof & (flags["gene_role"] != "both")
    directional_top = is_lof_for_direction & (flags["score_rank_pct"] <= (1.0 - TOP_SCORE_PERCENTILE))
    directional_top = directional_top | (~is_lof_for_direction & top)

    flags["confidence_tier"] = "LOW"
    flags.loc[top & ~drv,             "confidence_tier"] = "MEDIUM"
    flags.loc[~directional_top & drv, "confidence_tier"] = "CONTEXT"
    flags.loc[directional_top & drv,  "confidence_tier"] = "HIGH"

    # Partition assertions:
    #   HIGH + CONTEXT == has_driver_alteration   (directional_top covers all driver pairs)
    #   HIGH + MEDIUM + CONTEXT + LOW == total    (complete partition)
    #   MEDIUM + LOW == ~has_driver_alteration    (MEDIUM/LOW are role-blind, use top)
    assert (flags.loc[flags.confidence_tier.isin(["HIGH","CONTEXT"]), "has_driver_alteration"] == True).all(), \
        "HIGH/CONTEXT contains non-driver rows"
    assert not flags.loc[flags.confidence_tier == "MEDIUM", "has_driver_alteration"].any(), \
        "MEDIUM predicate regression: MEDIUM contains driver-alteration rows"
    assert flags.confidence_tier.value_counts().sum() == len(flags), \
        "Tier counts don't sum to total — tier assignment has gaps or overlaps"

    flags.to_parquet(PREDICTIONS, index=False)
    print(f"Predictions: {len(flags):,} rows  ->  {PREDICTIONS}")
    tier_counts = flags["confidence_tier"].value_counts()
    for t, n in tier_counts.items():
        print(f"  {t:<10s} {n:>8,}  ({n/len(flags)*100:.1f}%)")


if __name__ == "__main__":
    run()
