"""
Scoring/confidence_tiers_v2_no_evidence.py  (CANDIDATE -- builds on the live
confidence_tiers.py, unchanged except for one targeted addition: a distinct
NO_EVIDENCE tier for core_score.isna() rows.)
-----------------------------------------------------------------------------
Assigns confidence tiers to (gene, cell-line) predictions using:
  - core_score percentile within gene
  - has_driver_alteration flag
  - gene regime (activation_driven vs loss_of_function vs abundance_tracking)

Tier definitions:
  HIGH        — top 20% score AND has driver alteration
  MEDIUM      — top 20% score AND NOT has driver alteration
  CONTEXT     — NOT top 20% AND has driver alteration
  LOW         — NOT top 20% AND NOT has driver alteration
  NO_EVIDENCE — core_score is NaN (no interpretable expression evidence at
                all -- distinct from "measured and found unremarkable")

NO_EVIDENCE fix (this task, not yet promoted):
  The live confidence_tiers.py sets `confidence_tier = "LOW"` as an
  unconditional default before any score-based override. For a row where
  `core_score` is NaN (no RNA/protein score exists for this gene x line at
  all), `score_rank_pct` is also NaN, so `top` (`score_rank_pct >= 0.80`)
  evaluates False for that row -- comparisons against NaN are always False,
  never an error. This means a NaN row is indistinguishable, downstream,
  from a row that was genuinely measured and found NOT in the top 20% --
  "no evidence" silently collapses into "measured and unremarkable".

  Two independent sources of core_score=NaN rows exist in this pipeline:
    1. The core_score.py `.stack()` NaN-row-survival bug (see
       core_score_v9_nan_fix.py, this task's Part A) -- fixed there, at the
       source, not here.
    2. driver_routing.py's OWN, separate, INTENTIONAL recovery step: it adds
       back alteration-driver pairs with no expression score at all (the
       gene was silent in that lineage, or the line was never profiled),
       deliberately setting core_score=NaN, n_layers=0 -- documented in
       driver_routing.py's own comment as recovering ~32% of mutation
       drivers, ~37% of fusion drivers, ~56% of CNA alterations that would
       otherwise be silently lost. This source is NOT affected by Part A's
       fix at all -- it is a legitimate, separate NaN-producing path.
       Confirmed directly (this task): re-running driver_routing.py's logic
       against the Part-A-fixed core_score_v9_nan_fix.parquet still
       produces core_score=NaN rows via this recovery step -- Part B is
       independently necessary regardless of Part A.

  Because `has_driver_alteration` is True by construction for every
  recovery-step row (it only exists because a mutation/fusion/CNA driver
  was detected), such rows currently land in CONTEXT under the live logic
  (not LOW) -- CONTEXT's definition ("alteration present, score low")
  is still a misrepresentation for these rows: there IS no score to be
  "low", only n_layers=0. Before Part A is applied, phantom-bug NaN rows
  WITHOUT independent driver evidence do collapse into LOW (confirmed
  directly against the live pipeline: CD86 had 1,397 phantom rows tiered
  LOW, 36 tiered CONTEXT). Either way, `core_score.isna()` rows should never
  be tiered as if they were a real, measured outcome.

  Fix: after the existing HIGH/MEDIUM/CONTEXT/LOW assignment (order
  preserved exactly, so every score-based/driver-based row is tiered
  identically to before), a FINAL override sets confidence_tier="NO_EVIDENCE"
  wherever core_score.isna() -- this must run LAST, not as an early default,
  because the existing CONTEXT assignment's `~directional_top & drv` term
  evaluates True for every driver-flagged NaN row (directional_top is False
  for NaN inputs, same NaN-comparison-is-False mechanism as `top` above) and
  would silently overwrite an early NO_EVIDENCE default.

  This changes behaviour ONLY for core_score.isna() rows. Every row with a
  real, non-NaN score is tiered by exactly the same boolean logic as the
  live file, unchanged.

Reads from:  final_pipeline/outputs/{core_score, flags_with_driver}.parquet
             reference/gene_lookup.parquet  (gene regime via gene_role)
Writes to:   final_pipeline/outputs/predictions_with_confidence_v2_no_evidence.parquet
             (candidate output -- NOT the live predictions_with_confidence.parquet)
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from config import CORE_SCORE, FLAGS_DRIVER, GENE_LKP, PREDICTIONS

TOP_SCORE_PERCENTILE = 0.80

OUT_PREDICTIONS = PREDICTIONS.with_name("predictions_with_confidence_v2_no_evidence.parquet")


def run():
    print("=" * 70)
    print("Scoring — confidence tiers v2 (NO_EVIDENCE fix, candidate)")
    print("=" * 70)

    for p in [FLAGS_DRIVER, CORE_SCORE]:
        if not p.exists():
            raise SystemExit(f"Missing input: {p}")

    flags = pd.read_parquet(FLAGS_DRIVER)
    glk   = pd.read_parquet(GENE_LKP, columns=["ensg_id","gene_role","hgnc_symbol"])

    flags = flags.merge(glk, on="ensg_id", how="left")

    flags["score_rank_pct"] = (
        flags.groupby("ensg_id")["core_score"]
        .rank(pct=True, ascending=True, method="average")
    )

    top = flags["score_rank_pct"] >= TOP_SCORE_PERCENTILE
    drv = flags["has_driver_alteration"]
    is_lof = flags["is_lof_alteration"]

    is_lof_for_direction = is_lof & (flags["gene_role"] != "both")
    directional_top = is_lof_for_direction & (flags["score_rank_pct"] <= (1.0 - TOP_SCORE_PERCENTILE))
    directional_top = directional_top | (~is_lof_for_direction & top)

    flags["confidence_tier"] = "LOW"
    flags.loc[top & ~drv,             "confidence_tier"] = "MEDIUM"
    flags.loc[~directional_top & drv, "confidence_tier"] = "CONTEXT"
    flags.loc[directional_top & drv,  "confidence_tier"] = "HIGH"

    # NO_EVIDENCE fix: must run LAST -- see module docstring for why an early
    # default would be silently overwritten by the CONTEXT assignment above
    # for every driver-flagged NaN row.
    no_evidence = flags["core_score"].isna()
    flags.loc[no_evidence, "confidence_tier"] = "NO_EVIDENCE"

    # Partition assertions, updated to account for NO_EVIDENCE:
    #   Among SCORED rows only: HIGH + CONTEXT == has_driver_alteration
    #   Among SCORED rows only: MEDIUM + LOW == ~has_driver_alteration
    #   HIGH + MEDIUM + CONTEXT + LOW + NO_EVIDENCE == total (complete partition)
    #   NO_EVIDENCE rows are EXACTLY the core_score.isna() rows (nothing more, nothing less)
    scored = flags.loc[~no_evidence]
    assert (scored.loc[scored.confidence_tier.isin(["HIGH","CONTEXT"]), "has_driver_alteration"] == True).all(), \
        "HIGH/CONTEXT (scored rows) contains non-driver rows"
    assert not scored.loc[scored.confidence_tier == "MEDIUM", "has_driver_alteration"].any(), \
        "MEDIUM predicate regression: MEDIUM contains driver-alteration rows"
    assert flags.confidence_tier.value_counts().sum() == len(flags), \
        "Tier counts don't sum to total — tier assignment has gaps or overlaps"
    assert (flags.loc[flags.confidence_tier == "NO_EVIDENCE", "core_score"].isna()).all(), \
        "NO_EVIDENCE contains a row with a real core_score"
    assert (flags.loc[no_evidence, "confidence_tier"] == "NO_EVIDENCE").all(), \
        "A core_score.isna() row escaped the NO_EVIDENCE override"

    flags.to_parquet(OUT_PREDICTIONS, index=False)
    print(f"Predictions: {len(flags):,} rows  ->  {OUT_PREDICTIONS}")
    tier_counts = flags["confidence_tier"].value_counts()
    for t, n in tier_counts.items():
        print(f"  {t:<12s} {n:>8,}  ({n/len(flags)*100:.1f}%)")


if __name__ == "__main__":
    run()
