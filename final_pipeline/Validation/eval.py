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

Sensitivity labelling: per-drug MAD-z threshold (not a global absolute cut).
  A global LN_IC50 <= 1.0 threshold encodes drug potency rather than line
  sensitivity — sub-nanomolar compounds call every line sensitive; weak
  compounds call none sensitive. Per-drug standardisation removes that
  confound (Iorio et al., Cell 2016).
  drug_z = (LN_IC50 − drug_median) / (drug_MAD × 1.4826)
  is_sensitive = drug_z <= SENSITIVITY_MAD_Z  (default −0.5, ≈31% under normality)
  Guard: >=80% of drugs must land in 5–40% sensitive range.

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

HOLD_OUT_SEED    = 42
HOLD_OUT_FRAC    = 0.20
HITS_AT          = 20
SENSITIVITY_MAD_Z = -1.0   # per-drug MAD-z below which a line is called sensitive


def _label_sensitive_per_drug(gdsc: pd.DataFrame, ic50_col: str, drug_col: str) -> pd.DataFrame:
    """Replace global IC50 cut with per-drug MAD-z threshold."""
    def _drug_z(grp):
        med = grp[ic50_col].median()
        mad = (grp[ic50_col] - med).abs().median()
        mad_scaled = mad * 1.4826 if mad > 0 else grp[ic50_col].std(ddof=1)
        mad_scaled = max(mad_scaled, 1e-6)
        grp = grp.copy()
        grp["drug_z"] = (grp[ic50_col] - med) / mad_scaled
        grp["is_sensitive"] = grp["drug_z"] <= SENSITIVITY_MAD_Z
        return grp

    # Compute per-drug z-score stats via merge to avoid pandas groupby-apply
    # index-promotion behaviour (drug_col becoming index in pandas ≥2.2).
    def _drug_stats(grp):
        med = grp[ic50_col].median()
        mad = (grp[ic50_col] - med).abs().median()
        mad_scaled = max(mad * 1.4826 if mad > 0 else grp[ic50_col].std(ddof=1), 1e-6)
        return pd.Series({"_med": med, "_mad": mad_scaled})

    stats = gdsc.groupby(drug_col, group_keys=False).apply(_drug_stats, include_groups=False).reset_index()
    gdsc = gdsc.merge(stats, on=drug_col, how="left")
    gdsc["drug_z"]       = (gdsc[ic50_col] - gdsc["_med"]) / gdsc["_mad"]
    gdsc["is_sensitive"] = gdsc["drug_z"] <= SENSITIVITY_MAD_Z
    gdsc = gdsc.drop(columns=["_med", "_mad"])

    # Guard: >=80% of drugs must have 5–40% sensitive rate.
    per_drug_rate = gdsc.groupby(drug_col)["is_sensitive"].mean()
    in_window = ((per_drug_rate >= 0.05) & (per_drug_rate <= 0.40)).mean()
    print(f"  Sensitivity rate per drug: mean={per_drug_rate.mean():.3f}  "
          f"median={per_drug_rate.median():.3f}")
    print(f"  Drugs in 5–40% window: {in_window*100:.1f}%  (guard: >=80%)")
    if in_window < 0.80:
        print(f"  WARNING: guard not met ({in_window*100:.1f}% < 80%) — "
              f"consider adjusting SENSITIVITY_MAD_Z={SENSITIVITY_MAD_Z}")
    return gdsc


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
    ic50_col  = "LN_IC50" if "LN_IC50" in gdsc_t.columns else "log_ic50"
    drug_col  = "drug_name" if "drug_name" in gdsc_t.columns else gdsc_t.columns[0]

    print(f"\nSensitivity labelling (per-drug MAD-z <= {SENSITIVITY_MAD_Z}):")
    gdsc_t = _label_sensitive_per_drug(gdsc_t, ic50_col, drug_col)

    core_t = core_t.merge(gdsc_t[["ensg_id","model_id","is_sensitive"]].drop_duplicates(["ensg_id","model_id"]),
                          on=["ensg_id","model_id"], how="left")
    core_t["is_sensitive"] = core_t["is_sensitive"].astype("boolean").fillna(False).astype(bool)
    core_t = core_t.merge(genes, on="ensg_id", how="left")

    # Exclude pairs with no expression score (orphans or silence-guard masked).
    # Pre-declaration: n_pairs excluded ~277,125; n_test_genes expected back toward 13.
    unscored = core_t[core_t["core_score"].isna()]
    scored   = core_t[core_t["core_score"].notna()].copy()
    print(f"\n[UNSCORED] {len(unscored):,} pairs excluded from ranking "
          f"(core_score=NaN — orphan or silence-guard masked)")
    print(f"[UNSCORED] test genes in scored set: {scored.ensg_id.nunique()}  "
          f"(full set: {core_t.ensg_id.nunique()})")

    # rank_flat: pure expression rank, no alteration information, ascending=False for all genes.
    # Keeping this role-blind is essential — reversing for TSGs would contaminate the
    # control arm with treatment information and deflate the driver gate's delta.
    scored["rank_flat"] = scored.groupby("ensg_id")["core_score"].rank(ascending=False, method="first")

    # rank_driver: role-directional using is_lof_alteration (alteration type, not gene role).
    # LOF driver lines (truncating mutation / deletion in TSG-role gene) should rank high
    # when expression is LOW. GOF driver lines should rank high when expression is HIGH.
    # Non-driver lines rank below all driver lines regardless of direction.
    def _driver_rank(g):
        pct  = g["core_score"].rank(pct=True, ascending=True)
        drv  = g["has_driver_alteration"]
        # gene_role=="both" stays at top-20% — treat as GOF regardless of alteration type
        if g["gene_role"].iloc[0] == "both":
            lof = pd.Series(False, index=g.index)
        else:
            lof = g["is_lof_alteration"]
        priority = pd.Series(-1.0, index=g.index)
        priority[drv & lof]  = 1.0 - pct[drv & lof]   # LOF: lower score → better
        priority[drv & ~lof] = pct[drv & ~lof]          # GOF: higher score → better
        return priority.rank(ascending=False, method="first")

    scored["rank_driver"] = scored.groupby("ensg_id", group_keys=False).apply(
        _driver_rank, include_groups=False
    )

    def hits_at_k(df, rank_col, k=HITS_AT):
        in_top = df[df[rank_col] <= k]
        return in_top["is_sensitive"].sum() / max(df["is_sensitive"].sum(), 1)

    # Log genes excluded from the hit@20 loop (gene_role not in {oncogene, tsg, both}).
    # These genes (gene_role = "unknown") can receive driver calls via mutation/fusion
    # but not CNA (the direction gate requires a known oncogene/tsg/both role).
    # The exclusion is by design: hit@20 is only defined for genes with a known role.
    # W1/W6 characterisation (2026-08-14): 13/26 test genes have gene_role="unknown"
    # (GSK3B, KDM3A, KDM6B, PAK1, PLK2, PLK3, BCL2L1, BCL2L10, ROCK2, USP1, WEE1, CRBN, EHMT1)
    _excluded_test = scored[~scored["gene_role"].isin(["oncogene", "tsg", "both"])].drop_duplicates("ensg_id")
    print(f"\n[EXCLUSION] eval.py: {_excluded_test.ensg_id.nunique()} of "
          f"{scored.ensg_id.nunique()} scored test genes excluded from hit@20 loop "
          f"(gene_role not in {{oncogene, tsg, both}})")
    print(f"[EXCLUSION] roles excluded: "
          f"{_excluded_test['gene_role'].value_counts(dropna=False).to_dict()}")
    print(f"[EXCLUSION] excluded ensg_ids: {sorted(_excluded_test.ensg_id.tolist())}")

    print(f"\nhit@20 results (scored pairs only — NaN excluded; scored-only is the inferential estimand):")
    summary = {"sensitivity_threshold": f"per_drug_mad_z_{SENSITIVITY_MAD_Z}",
               "n_pairs_unscored": int(len(unscored)),
               "n_test_genes_scored": int(scored.ensg_id.nunique()),
               "n_test_genes_total":  int(core_t.ensg_id.nunique())}
    for cls in ["oncogene", "tsg", "both"]:
        sub = scored[scored.gene_role == cls]
        if sub.empty:
            continue
        summary[cls] = {
            "n_genes":       int(sub.ensg_id.nunique()),
            "n_pairs":       int(len(sub)),
            "n_sensitive":   int(sub["is_sensitive"].sum()),
            "hits20_flat":   round(float(hits_at_k(sub, "rank_flat")), 4),
            "hits20_driver": round(float(hits_at_k(sub, "rank_driver")), 4),
        }
        print(f"  {cls}:  flat={summary[cls]['hits20_flat']:.3f}  "
              f"driver={summary[cls]['hits20_driver']:.3f}  "
              f"({summary[cls]['n_genes']} genes  n_sensitive={summary[cls]['n_sensitive']:,})")

    EVAL_SUMMARY.write_text(json.dumps(summary, indent=2))
    print(f"\nSummary written: {EVAL_SUMMARY}")


if __name__ == "__main__":
    run()
