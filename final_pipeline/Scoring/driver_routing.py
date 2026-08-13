"""
Scoring/driver_routing.py
--------------------------
Driver-gated routing: merges core_score with alteration evidence (mutations,
fusions, CNA) to produce the has_driver_alteration flag used by confidence tiers.

Gate logic:
  has_driver_alteration = mut_driver OR fusion_driver OR cna_alteration

Reads from:  final_pipeline/outputs/{core_score, mutations_scores, fusions_scores, cna_flags}.parquet
Writes to:   final_pipeline/outputs/flags_with_driver.parquet

Schema: model_id, ensg_id, core_score, n_layers, stratum_rank,
        p_mutation, p_fusion, has_cna_alteration,
        has_driver_alteration, has_alteration
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from config import CORE_SCORE, MUTATIONS_SCR, FUSIONS_SCR, CNA_FLAGS, FLAGS_DRIVER

MUT_DRIVER_THRESHOLD  = 0.5   # p_mutation above this = driver mutation
FUS_DRIVER_THRESHOLD  = 0.5   # p_fusion above this = driver fusion


def run():
    print("=" * 70)
    print("Scoring — driver-gated routing")
    print("=" * 70)

    for p in [CORE_SCORE, MUTATIONS_SCR, FUSIONS_SCR, CNA_FLAGS]:
        if not p.exists():
            raise SystemExit(f"Missing input: {p}\nRun prior stages first.")

    core  = pd.read_parquet(CORE_SCORE)
    muts  = pd.read_parquet(MUTATIONS_SCR)
    fus   = pd.read_parquet(FUSIONS_SCR, columns=["ensg_id","model_id","p_fusion"])
    cna   = pd.read_parquet(CNA_FLAGS,   columns=["ensg_id","model_id","has_cna_alteration"])

    print(f"core_score: {len(core):,}  mutations: {len(muts):,}  "
          f"fusions: {len(fus):,}  cna: {len(cna):,}")

    flags = core.merge(muts[["ensg_id","model_id","p_mutation"]], on=["ensg_id","model_id"], how="left")
    flags = flags.merge(fus,  on=["ensg_id","model_id"], how="left")
    flags = flags.merge(cna,  on=["ensg_id","model_id"], how="left")

    flags["p_mutation"]         = flags["p_mutation"].fillna(0.0)
    flags["p_fusion"]           = flags["p_fusion"].fillna(0.0)
    flags["has_cna_alteration"] = flags["has_cna_alteration"].fillna(False)

    flags["mut_driver"]    = flags["p_mutation"] >= MUT_DRIVER_THRESHOLD
    flags["fusion_driver"] = flags["p_fusion"]   >= FUS_DRIVER_THRESHOLD

    flags["has_driver_alteration"] = (
        flags["mut_driver"] | flags["fusion_driver"] | flags["has_cna_alteration"]
    )
    flags["has_alteration"] = flags["p_mutation"] > 0

    out_cols = ["model_id","ensg_id","core_score","n_layers","stratum_rank",
                "p_mutation","p_fusion","has_cna_alteration",
                "has_driver_alteration","has_alteration"]
    flags[out_cols].to_parquet(FLAGS_DRIVER, index=False)

    n_driver = flags["has_driver_alteration"].sum()
    print(f"\nhas_driver_alteration: {n_driver:,} / {len(flags):,} pairs "
          f"({n_driver/len(flags)*100:.1f}%)")
    print(f"Written: {FLAGS_DRIVER}")


if __name__ == "__main__":
    run()
