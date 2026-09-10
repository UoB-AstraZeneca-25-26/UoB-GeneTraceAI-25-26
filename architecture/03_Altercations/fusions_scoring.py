"""
03_Altercations/fusions_scoring.py
------------------------------------
Scores gene fusions per (gene, cell-line) as p_fusion, a probability of
functional impact using three evidence channels:

  p_conf  — fusion caller confidence (low=1, medium=2, high=3 ordinal)
  p_ffpm  — FFPM (fragments per fused-gene million) — ABSOLUTE scale
  p_recur — fusion count across samples — ABSOLUTE count scale

Combination: p_base = 1 - (1-p_conf)(1-p_ffpm)(1-p_recur)
             p_fusion = min(p_base + 0.15 * in_frame_flag, 1.0)

Hill function:  hill(x, p0, k) = x^k / (x^k + p0^k)
Hyperparameters:
  Conf ordinal   p0=2.0 (half-max at medium); k=1.5
                 low≈0.26, medium=0.50, high≈0.65
                 (p0=1.0 set floor at 0.50, making all fusions drivers)
  FFPM absolute  p0=0.5  k=2.0  (half-max at 0.5 FFPM; median raw=0.08)
                 Replaces within-gene percentile rank, which inflated
                 median score to 0.774 regardless of absolute evidence.
  Recur count    p0=2.0  k=2.0  (half-max at 2 occurrences)
                 singleton→0.20, seen in 2+→0.50, 5+→0.86
                 (87% of fusions are singletons; percentile rank made
                  them all score 0.50+, violating Noisy-OR independence)
  In-frame boost b=0.15

Reads from:  cleaned_track_data/fusions_gene_level.parquet
Writes to:   final_pipeline/outputs/fusions_scores.parquet

Schema: ensg_id, model_id, p_fusion, conf_scored, ffpm_scored, recur_scored
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from config import CLEANED, FUSIONS_SCR

CONF_MAP      = {"low": 1, "medium": 2, "high": 3}
INFRAME_BOOST = 0.15


def hill(x: pd.Series, p0: float, k: float) -> pd.Series:
    xk = np.power(np.clip(x, 0, None), k)
    return xk / (xk + p0 ** k)


def run():
    print("=" * 70)
    print("03_Altercations — fusions scoring")
    print("=" * 70)

    df = pd.read_parquet(CLEANED / "fusions_gene_level.parquet")
    print(f"fusions_gene_level: {len(df):,} rows x {df.shape[1]} cols")

    df = df.copy()
    df["conf_ord"] = df["max_confidence"].map(CONF_MAP).astype(float)

    # Absolute scales — percentile ranks inflated all fusions to ~0.75+ regardless
    # of evidence, making Noisy-OR always ≥0.95 and fusion_driver always TRUE.
    df["p_conf"]  = hill(df["conf_ord"],   p0=2.0, k=1.5)  # half-max at medium; low≈0.26
    df["p_ffpm"]  = hill(df["best_ffpm"],  p0=0.5, k=2.0)  # half-max at 0.5 FFPM
    df["p_recur"] = hill(df["fusion_count"], p0=2.0, k=2.0)  # half-max at 2 occurrences

    df["conf_scored"]  = df["conf_ord"].notna()
    df["ffpm_scored"]  = df["best_ffpm"].notna()
    df["recur_scored"] = df["fusion_count"].notna()

    p_base  = 1.0 - (1.0 - df["p_conf"].fillna(0)) * \
                    (1.0 - df["p_ffpm"].fillna(0)) * \
                    (1.0 - df["p_recur"].fillna(0))
    inframe = df["any_in_frame"].fillna(False).astype(int)
    df["p_fusion"] = (p_base + INFRAME_BOOST * inframe).clip(upper=1.0)

    scores = (
        df.groupby(["ensg_id", "model_id"], sort=False)
        .agg(
            p_fusion    =("p_fusion",    "max"),
            conf_scored =("conf_scored", "any"),
            ffpm_scored =("ffpm_scored", "any"),
            recur_scored=("recur_scored","any"),
        )
        .reset_index()
    )
    scores["model_id"] = scores["model_id"].str.lower()
    scores.to_parquet(FUSIONS_SCR, index=False)
    print(f"Output rows: {len(scores):,}  ->  {FUSIONS_SCR}")
    print(f"p_fusion: mean={scores.p_fusion.mean():.3f}  median={scores.p_fusion.median():.3f}")


if __name__ == "__main__":
    run()
