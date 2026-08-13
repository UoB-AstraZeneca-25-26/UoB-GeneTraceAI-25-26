"""
Scoring/confidence_tiers.py
-----------------------------
Assigns confidence tiers to (gene, cell-line) predictions using:
  - core_score percentile within gene
  - has_driver_alteration flag
  - gene regime (activation_driven vs loss_of_function vs abundance_tracking)

Tier definitions:
  HIGH    — top 20% score AND has driver alteration
  MEDIUM  — top 20% score OR has driver alteration (not both)
  LOW     — bottom 80%, no driver alteration
  CONTEXT — gene has known alteration but score is low

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

    flags["confidence_tier"] = "LOW"
    flags.loc[top & ~drv,  "confidence_tier"] = "MEDIUM"
    flags.loc[~top & drv,  "confidence_tier"] = "CONTEXT"
    flags.loc[top & drv,   "confidence_tier"] = "HIGH"

    flags.to_parquet(PREDICTIONS, index=False)
    print(f"Predictions: {len(flags):,} rows  ->  {PREDICTIONS}")
    tier_counts = flags["confidence_tier"].value_counts()
    for t, n in tier_counts.items():
        print(f"  {t:<10s} {n:>8,}  ({n/len(flags)*100:.1f}%)")


if __name__ == "__main__":
    run()
