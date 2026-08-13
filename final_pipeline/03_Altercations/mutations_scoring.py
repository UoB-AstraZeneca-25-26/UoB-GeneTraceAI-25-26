"""
03_Altercations/mutations_scoring.py
--------------------------------------
Scores somatic mutations per (gene, cell-line) as a probability of functional
impact using a Noisy-OR combination of three evidence channels:

  p_vep    — VEP consequence rank (SIFT/PolyPhen: 1=low, 5=high)
  p_path   — AlphaMissense / CADD pathogenicity score (0-1)
  p_burden — variant count per (gene, cell-line)

Combination: p_base = 1 - (1-p_vep)(1-p_path)(1-p_burden)
             p_mutation = min(p_base + 0.15 * driver_flag, 1.0)

Hill function:  hill(x, p0, k) = x^k / (x^k + p0^k)
Hyperparameters:
  VEP     p0=3.0  k=2.0
  Path    p0=0.5  k=2.0
  Burden  p0=3.0  k=1.5
  Driver boost  b=0.15  (oncogene or TSG hit)

Reads from:  cleaned_track_data/mutations_collapsed.parquet
             reference/gene_lookup.parquet
Writes to:   final_pipeline/outputs/mutations_scores.parquet

Schema: ensg_id, model_id, p_mutation
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from config import CLEANED, MUTATIONS_SCR

P0_VEP    = 3.0;  K_VEP    = 2.0
P0_PATH   = 0.5;  K_PATH   = 2.0
P0_BURDEN = 3.0;  K_BURDEN = 1.5
DRIVER_BOOST = 0.15


def hill(x: pd.Series, p0: float, k: float) -> pd.Series:
    xk = np.power(x.clip(lower=0), k)
    return xk / (xk + p0 ** k)


def run():
    print("=" * 70)
    print("03_Altercations — mutations scoring")
    print("=" * 70)

    mut = pd.read_parquet(CLEANED / "mutations_collapsed.parquet")
    print(f"mutations_collapsed: {len(mut):,} rows x {mut.shape[1]} cols")

    df = mut[["ensg_id", "model_id", "max_vep_rank", "max_pathogenicity",
              "variant_count", "multi_hit_high_impact", "any_driver"]].copy()

    # variant_burden = variant_count, capped at 5 to avoid saturation
    df["variant_burden"] = df["variant_count"].clip(upper=5)

    df["p_vep"]    = hill(df["max_vep_rank"].astype(float),    P0_VEP,    K_VEP)
    df["p_path"]   = hill(df["max_pathogenicity"].astype(float), P0_PATH,  K_PATH)
    df["p_burden"] = hill(df["variant_burden"].astype(float),  P0_BURDEN, K_BURDEN)

    p_vep    = df["p_vep"].fillna(0.0)
    p_path   = df["p_path"].fillna(0.0)
    p_burden = df["p_burden"].fillna(0.0)
    p_base   = 1.0 - (1.0 - p_vep) * (1.0 - p_path) * (1.0 - p_burden)

    # any_driver already annotated in mutations_collapsed (oncogene hit or TSG hit)
    driver_flag = df["any_driver"].fillna(False).astype(int)

    df["p_mutation"] = (p_base + DRIVER_BOOST * driver_flag).clip(upper=1.0)

    scores = (
        df.groupby(["ensg_id", "model_id"], sort=False)["p_mutation"]
        .max()
        .reset_index()
    )
    scores["model_id"] = scores["model_id"].str.lower()
    scores.to_parquet(MUTATIONS_SCR, index=False)
    print(f"Output rows: {len(scores):,}  ->  {MUTATIONS_SCR}")
    print(f"p_mutation: mean={scores.p_mutation.mean():.3f}  median={scores.p_mutation.median():.3f}")


if __name__ == "__main__":
    run()
