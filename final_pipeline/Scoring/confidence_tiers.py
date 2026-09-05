"""
Scoring/confidence_tiers.py
-----------------------------
Assigns confidence tiers to (gene, cell-line) predictions using:
  - core_score percentile within gene
  - has_driver_alteration flag
  - gene regime (activation_driven vs loss_of_function vs abundance_tracking)

Tier definitions:
  HIGH        — top 20% score AND has driver alteration
  MEDIUM      — top 20% score AND NOT has driver alteration  (high score, no alteration support)
  CONTEXT     — NOT top 20% AND has driver alteration        (alteration present, score low)
  LOW         — NOT top 20% AND NOT has driver alteration
  NO_EVIDENCE — core_score is NaN (no interpretable expression evidence at all)

NO_EVIDENCE fix (promoted 2026-08-31): a distinct tier for core_score.isna() rows.

Previously, `confidence_tier = "LOW"` was an unconditional default before any
score-based override. For a core_score.isna() row, score_rank_pct is also
NaN, so `top` (score_rank_pct >= 0.80) evaluates False -- comparisons
against NaN are always False, never an error -- so a NaN row was
indistinguishable, downstream, from a row genuinely measured and found NOT
in the top 20%: "no evidence" silently collapsed into "measured and
unremarkable".

Two independent sources of core_score=NaN rows feed this file:
  1. core_score.py's `.stack()` NaN-row-survival bug -- fixed at the source
     (see core_score.py's own docstring, promoted the same day); no longer
     produces NaN rows as of that fix.
  2. driver_routing.py's own, separate, INTENTIONAL recovery step: it adds
     back alteration-driver pairs with no expression score at all,
     deliberately setting core_score=NaN, n_layers=0. This source is NOT
     affected by the core_score.py fix -- confirmed directly (this
     promotion's own verification): re-running driver_routing.py against
     the NaN-fixed core_score.parquet still produces 333,417 such rows,
     100% correctly tagged n_layers=0. NO_EVIDENCE is necessary regardless
     of the core_score.py fix, not made redundant by it.

Because has_driver_alteration is True by construction for every recovery-
step row, such rows land in CONTEXT under the pre-fix logic (not LOW) --
still a misrepresentation, since there is no score to be "low", only
n_layers=0. core_score.isna() rows should never be tiered as if they were a
real, measured outcome.

Implementation note: the NO_EVIDENCE override runs LAST, after the existing
HIGH/MEDIUM/CONTEXT/LOW assignment, not as an early default -- an early
default would be silently overwritten by the CONTEXT assignment's
`~directional_top & drv` term, which evaluates True for every driver-
flagged NaN row (directional_top is False for NaN inputs, the same
NaN-comparison-is-False mechanism as `top` above). This changes behaviour
ONLY for core_score.isna() rows; every row with a real score is tiered by
exactly the same boolean logic as before.

Verified before promotion (confidence_tiers_v2_no_evidence.py, this fix's
candidate file): in isolation and combined with the core_score.py fix,
LOW/MEDIUM/CONTEXT/HIGH counts for real rows were byte-identical
(13,630,183 / 3,406,451 / 255,900 / 77,815) regardless of which NaN source
fed NO_EVIDENCE; combined with the core_score.py fix, NO_EVIDENCE = 333,417,
exactly the legitimate orphan-recovery count. See docs/WHY_REFERENCE.md for
the full verification record and
confidence_tiers_v1_pre_no_evidence_fix_20260831_112601.py.bak /
predictions_with_confidence_v1_pre_no_evidence_fix_20260831_112601.parquet.bak
for the pre-promotion snapshots.

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

    flags.to_parquet(PREDICTIONS, index=False)
    print(f"Predictions: {len(flags):,} rows  ->  {PREDICTIONS}")
    tier_counts = flags["confidence_tier"].value_counts()
    for t, n in tier_counts.items():
        print(f"  {t:<12s} {n:>8,}  ({n/len(flags)*100:.1f}%)")


if __name__ == "__main__":
    run()
