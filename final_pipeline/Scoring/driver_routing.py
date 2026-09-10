"""
Scoring/driver_routing.py
--------------------------
Driver-gated routing: merges core_score with alteration evidence (mutations,
fusions, CNA) to produce the has_driver_alteration flag used by confidence tiers.

Gate logic:
  has_driver_alteration = mut_driver OR fusion_driver OR cna_alteration
  is_lof_alteration     = alteration is loss-of-function (truncating mutation or
                          deletion in TSG/both gene) — used for directional ranking

Reads from:  final_pipeline/outputs/{core_score, mutations_scores, fusions_scores, cna_flags}.parquet
             reference/gene_lookup.parquet  (gene_role for is_lof_alteration)
Writes to:   final_pipeline/outputs/flags_with_driver.parquet

Schema: model_id, ensg_id, core_score, n_layers, stratum_rank,
        p_mutation, max_vep_rank, p_fusion, has_cna_alteration,
        has_driver_alteration, has_alteration, is_lof_alteration
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from config import CORE_SCORE, MUTATIONS_SCR, FUSIONS_SCR, CNA_FLAGS, FLAGS_DRIVER, GENE_LKP

MUT_DRIVER_THRESHOLD  = 0.5   # p_mutation above this = driver mutation
FUS_DRIVER_THRESHOLD  = 0.5   # p_fusion above this = driver fusion


def run():
    print("=" * 70)
    print("Scoring — driver-gated routing")
    print("=" * 70)

    for p in [CORE_SCORE, MUTATIONS_SCR, FUSIONS_SCR, CNA_FLAGS]:
        if not p.exists():
            raise SystemExit(f"Missing input: {p}\nRun prior stages first.")

    core     = pd.read_parquet(CORE_SCORE)
    muts     = pd.read_parquet(MUTATIONS_SCR)  # now includes max_vep_rank
    fus      = pd.read_parquet(FUSIONS_SCR, columns=["ensg_id","model_id","p_fusion"])
    cna      = pd.read_parquet(CNA_FLAGS,   columns=["ensg_id","model_id","has_cna_alteration"])
    gene_lkp = pd.read_parquet(GENE_LKP,   columns=["ensg_id","gene_role"])
    # cna_layer.py writes lowercase ensg_id; upstream files use uppercase → normalise here
    cna["ensg_id"] = cna["ensg_id"].str.upper()

    print(f"core_score: {len(core):,}  mutations: {len(muts):,}  "
          f"fusions: {len(fus):,}  cna: {len(cna):,}")

    flags = core.merge(
        muts[["ensg_id","model_id","p_mutation","max_vep_rank"]],
        on=["ensg_id","model_id"], how="left")
    flags = flags.merge(fus,  on=["ensg_id","model_id"], how="left")
    flags = flags.merge(cna,  on=["ensg_id","model_id"], how="left")

    # Recover alteration-driver pairs absent from core_score.
    # These are (gene, cell-line) pairs where a mutation/fusion/CNA driver was
    # detected but the line has no RNA/protein score for that gene — either the
    # gene is silent in the line's lineage (silence guard masked it) or the line
    # was never profiled for this gene. Without recovery, ~32% of mutation drivers,
    # ~37% of fusion drivers, and ~56% of CNA alterations are silently lost.
    # Recovered pairs receive core_score=NaN (no interpretable expression evidence);
    # eval.py excludes NaN pairs from ranking.
    core_key = frozenset(core["ensg_id"] + "|||" + core["model_id"])

    muts_drv = muts.loc[muts["p_mutation"] >= MUT_DRIVER_THRESHOLD,
                         ["ensg_id", "model_id", "p_mutation", "max_vep_rank"]].copy()
    fus_drv  = fus.loc[fus["p_fusion"] >= FUS_DRIVER_THRESHOLD,
                        ["ensg_id", "model_id", "p_fusion"]].copy()
    cna_drv  = cna.loc[cna["has_cna_alteration"],
                        ["ensg_id", "model_id", "has_cna_alteration"]].copy()

    def _orphans(df):
        key = df["ensg_id"] + "|||" + df["model_id"]
        return df[~key.isin(core_key)]

    orphan_muts = _orphans(muts_drv)
    orphan_fus  = _orphans(fus_drv)
    orphan_cna  = _orphans(cna_drv)

    orphan_pairs = pd.concat([
        orphan_muts[["ensg_id", "model_id"]],
        orphan_fus[["ensg_id", "model_id"]],
        orphan_cna[["ensg_id", "model_id"]],
    ]).drop_duplicates()

    if len(orphan_pairs) > 0:
        orphan_pairs = orphan_pairs.merge(
            orphan_muts[["ensg_id", "model_id", "p_mutation", "max_vep_rank"]],
            on=["ensg_id", "model_id"], how="left")
        orphan_pairs = orphan_pairs.merge(
            orphan_fus[["ensg_id", "model_id", "p_fusion"]],
            on=["ensg_id", "model_id"], how="left")
        orphan_pairs = orphan_pairs.merge(
            orphan_cna[["ensg_id", "model_id", "has_cna_alteration"]],
            on=["ensg_id", "model_id"], how="left")
        orphan_pairs["core_score"]         = np.nan   # no interpretable expression evidence
        orphan_pairs["n_layers"]           = 0
        orphan_pairs["stratum_rank"]       = np.nan   # undefined for unscored pairs
        orphan_pairs["p_mutation"]         = orphan_pairs["p_mutation"].fillna(0.0)
        orphan_pairs["max_vep_rank"]       = orphan_pairs["max_vep_rank"].fillna(0).astype(int)
        orphan_pairs["p_fusion"]           = orphan_pairs["p_fusion"].fillna(0.0)
        orphan_pairs["has_cna_alteration"] = orphan_pairs["has_cna_alteration"].fillna(False)
        flags = pd.concat([flags, orphan_pairs], ignore_index=True)
        print(f"  Recovered {len(orphan_pairs):,} alteration-driver pairs absent from "
              f"core_score (mut: {len(orphan_muts):,}  fus: {len(orphan_fus):,}  "
              f"cna: {len(orphan_cna):,})")

    flags["p_mutation"]         = flags["p_mutation"].fillna(0.0)
    flags["max_vep_rank"]       = flags["max_vep_rank"].fillna(0).astype(int)
    flags["p_fusion"]           = flags["p_fusion"].fillna(0.0)
    flags["has_cna_alteration"] = flags["has_cna_alteration"].fillna(False)

    flags["mut_driver"]    = flags["p_mutation"] >= MUT_DRIVER_THRESHOLD
    flags["fusion_driver"] = flags["p_fusion"]   >= FUS_DRIVER_THRESHOLD

    flags["has_driver_alteration"] = (
        flags["mut_driver"] | flags["fusion_driver"] | flags["has_cna_alteration"]
    )
    flags["has_alteration"] = flags["p_mutation"] > 0

    # is_lof_alteration: alteration direction, not gene role alone.
    # A TSG gene with a truncating mutation OR a CNA deletion (which cna_layer.py
    # only assigns to tsg/both by construction) is classified LOF.
    # max_vep_rank >= 3 is the truncating proxy (HIGH/frameshift consequence rank).
    flags = flags.merge(gene_lkp, on="ensg_id", how="left")
    flags["gene_role"] = flags["gene_role"].fillna("unknown")
    is_tsg_or_both = flags["gene_role"].isin(["tsg", "both"])
    mut_truncating = flags["mut_driver"] & (flags["max_vep_rank"] >= 3)
    flags["is_lof_alteration"] = (
        (flags["has_cna_alteration"] & is_tsg_or_both) |
        (mut_truncating & is_tsg_or_both)
    )

    out_cols = ["model_id","ensg_id","core_score","n_layers","stratum_rank",
                "p_mutation","max_vep_rank","p_fusion","has_cna_alteration",
                "has_driver_alteration","has_alteration","is_lof_alteration"]
    flags[out_cols].to_parquet(FLAGS_DRIVER, index=False)

    n_driver = flags["has_driver_alteration"].sum()
    print(f"\nhas_driver_alteration: {n_driver:,} / {len(flags):,} pairs "
          f"({n_driver/len(flags)*100:.1f}%)")
    print(f"Written: {FLAGS_DRIVER}")


if __name__ == "__main__":
    run()
