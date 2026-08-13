"""
Validation/eval.py
-------------------
Held-out evaluation: hit@20 on a 20% gene hold-out (seed=42).

Three ranking strategies compared:
  flat          — core_score only, no alteration gating
  any-alt       — has_alteration flag boosts rank within gene class
  driver-gated  — has_driver_alteration flag applied per gene class

Metric: hit@20 — fraction of known sensitive lines in top 20 ranked lines.
Evaluated separately for activation_driven and abundance_tracking gene classes.

Reads from:  final_pipeline/outputs/{predictions_with_confidence, flags_with_driver}.parquet
             validation/prepared/gdsc_scored_ready.parquet
Writes to:   final_pipeline/outputs/stage4_eval_summary.json
"""
import json, sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from config import PREDICTIONS, FLAGS_DRIVER, VALIDATION, EVAL_SUMMARY, GENE_LKP

HOLD_OUT_SEED = 42
HOLD_OUT_FRAC = 0.20
HITS_AT       = 20
SENSITIVITY_IC50_THRESHOLD = 1.0   # log IC50 below this = sensitive


def run():
    print("=" * 70)
    print("Validation — held-out evaluation (hit@20)")
    print("=" * 70)

    gdsc_path = VALIDATION / "gdsc_scored_ready.parquet"
    for p in [PREDICTIONS, FLAGS_DRIVER, gdsc_path]:
        if not p.exists():
            raise SystemExit(f"Missing: {p}")

    flags  = pd.read_parquet(FLAGS_DRIVER)
    gdsc   = pd.read_parquet(gdsc_path)
    genes  = pd.read_parquet(GENE_LKP, columns=["ensg_id","gene_role"])

    # build gene set (intersection of scored genes and GDSC targets)
    curated = sorted(set(flags.ensg_id) & set(gdsc.target_ensg))
    rng = np.random.default_rng(HOLD_OUT_SEED)
    test_genes = set(rng.choice(curated, size=int(len(curated) * HOLD_OUT_FRAC), replace=False))

    core_t  = flags[flags.ensg_id.isin(test_genes)].copy()
    gdsc_t  = gdsc[gdsc.target_ensg.isin(test_genes)].copy()
    gdsc_t  = gdsc_t.rename(columns={"target_ensg":"ensg_id"})
    gdsc_t["model_id"] = gdsc_t["model_id"].str.lower()
    ic50_col = "LN_IC50" if "LN_IC50" in gdsc_t.columns else "log_ic50"
    gdsc_t["is_sensitive"] = gdsc_t[ic50_col] <= SENSITIVITY_IC50_THRESHOLD

    core_t = core_t.merge(gdsc_t[["ensg_id","model_id","is_sensitive"]], on=["ensg_id","model_id"], how="left")
    core_t["is_sensitive"] = core_t["is_sensitive"].fillna(False)
    core_t = core_t.merge(genes, on="ensg_id", how="left")

    core_t["rank_flat"] = core_t.groupby("ensg_id")["core_score"].rank(ascending=False, method="first")
    core_t["rank_driver"] = core_t.groupby("ensg_id").apply(
        lambda g: g["core_score"].where(g["has_driver_alteration"], other=-1)
                                 .rank(ascending=False, method="first")
    ).reset_index(level=0, drop=True)

    def hits_at_k(df, rank_col, k=HITS_AT):
        in_top = df[df[rank_col] <= k]
        return in_top["is_sensitive"].sum() / max(df["is_sensitive"].sum(), 1)

    summary = {}
    for cls in ["activation_driven", "abundance_tracking"]:
        sub = core_t[core_t.gene_role == cls]
        if sub.empty:
            continue
        summary[cls] = {
            "n_genes":       int(sub.ensg_id.nunique()),
            "n_pairs":       int(len(sub)),
            "hits20_flat":   round(float(hits_at_k(sub, "rank_flat")), 4),
            "hits20_driver": round(float(hits_at_k(sub, "rank_driver")), 4),
        }
        print(f"  {cls}:  flat={summary[cls]['hits20_flat']:.3f}  "
              f"driver={summary[cls]['hits20_driver']:.3f}  "
              f"({summary[cls]['n_genes']} genes)")

    EVAL_SUMMARY.write_text(json.dumps(summary, indent=2))
    print(f"\nSummary written: {EVAL_SUMMARY}")


if __name__ == "__main__":
    run()
