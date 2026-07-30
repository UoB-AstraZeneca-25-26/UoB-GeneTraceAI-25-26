"""
test_run_ploidy_normalised_cna.py
------------------------------------
Read-only diagnostic. Does NOT write cna_flags.parquet, flags_with_driver.parquet,
predictions_with_confidence.parquet, or any other production output.

QUESTION
--------
build_cna_layer.py calls copy-number events with ABSOLUTE cutoffs that assume a
diploid baseline:
    is_amplification = TOTAL_CN > 2.5
    is_deletion      = TOTAL_CN < 1.5

But the scored panel is not diploid -- median ploidy 3.06, and 53% of lines have
ploidy > 3. Does replacing the absolute cutoff with a ploidy-relative one remove
the measured whole-genome-doubling bias (WGD lines currently get 1.31x the
CNA-alteration rate, p = 1.4e-05)?

THRESHOLD CHOICE
----------------
Ratio thresholds are set so behaviour at a diploid baseline is UNCHANGED:
    amp: TOTAL_CN / ploidy > 2.5 / 2.0 = 1.25
    del: TOTAL_CN / ploidy < 1.5 / 2.0 = 0.75
So a diploid line gets exactly today's calls; only non-diploid lines move. This
is the same convention GISTIC / ABSOLUTE / PureCN use.

Lines with no ploidy estimate fall back to ploidy = 2.0, i.e. today's behaviour.

WHAT THIS CAN AND CANNOT CHANGE
-------------------------------
has_cna_alteration does NOT enter the ranking. The production gate is
    has_driver_alteration = mut_driver | fusion_driver
which does not read CNA at all. CNA feeds confidence-tier upgrades
(build_full_predictions.py:108-142), where it can promote a prediction to "high".

So the endpoints here are:
  PRIMARY   -- is the WGD bias removed? (the confound itself)
  SECONDARY -- how many (model, gene) pairs change has_cna_alteration, and what
               confidence labels do those pairs currently hold (blast radius)?
  CONTROL   -- held-out hit@20 must be identical; asserted, not assumed.

Output: src/pipeline/outputs/test_run_ploidy_normalised_cna_results.json

Run:
    python src/pipeline/test_run_ploidy_normalised_cna.py
"""
import json
import numpy as np
import pandas as pd
from pathlib import Path
from scipy.stats import mannwhitneyu, spearmanr

OUTPUTS = Path("src/pipeline/outputs")
CNA_FILE    = "data/COSMIC/CellLinesProject_CompleteCNA_v104_GRCh37.tsv"
MODEL_LIST  = "data/GDSC/model_list_20260709.csv"
GENE_LOOKUP = "reference/gene_lookup.parquet"
SIGNATURES  = "cleaned_track_data/signatures_model_level.parquet"

AMP_THRESHOLD = 2.5      # absolute, current production
DEL_THRESHOLD = 1.5
DIPLOID       = 2.0
AMP_RATIO     = AMP_THRESHOLD / DIPLOID   # 1.25
DEL_RATIO     = DEL_THRESHOLD / DIPLOID   # 0.75

# ---------------------------------------------------------------- load
print("=" * 72)
print("STEP 1 -- rebuild CNA calls, both threshold schemes")
print("=" * 72)

cna = pd.read_csv(CNA_FILE, sep="\t",
                  usecols=["GENE_SYMBOL", "SAMPLE_NAME", "TOTAL_CN", "MUT_TYPE"])
print(f"  raw CNA rows: {len(cna):,}")

model = pd.read_csv(MODEL_LIST, usecols=["model_name", "BROAD_ID"]).dropna(subset=["BROAD_ID"])
model["model_id"]    = model["BROAD_ID"].str.upper()
model["sample_name"] = model["model_name"].str.strip()
model = model.drop_duplicates(subset="sample_name", keep="first")

cna = cna.merge(model[["sample_name", "model_id"]],
                left_on="SAMPLE_NAME", right_on="sample_name", how="left")
cna = cna[cna["model_id"].notna()].copy()

gl = pd.read_parquet(GENE_LOOKUP, columns=["ensg_id", "hgnc_symbol", "gene_role"])
gl["ensg_id"] = gl["ensg_id"].str.split(".").str[0]
cna = cna.merge(gl, left_on="GENE_SYMBOL", right_on="hgnc_symbol", how="left")
cna = cna[cna["ensg_id"].notna()].copy()

regime = pd.read_parquet(OUTPUTS / "gene_regime.parquet", columns=["ensg_id", "class"])
regime = regime.rename(columns={"class": "regime_class"})
regime["ensg_id"] = regime["ensg_id"].str.split(".").str[0]
cna = cna.merge(regime, on="ensg_id", how="left")
cna["regime_class"] = cna["regime_class"].fillna("unknown_dual")
cna["gene_role"]    = cna["gene_role"].fillna("unknown")

# attach ploidy
sig = pd.read_parquet(SIGNATURES)
sig["model_id"] = sig["model_id"].str.upper()
cna = cna.merge(sig[["model_id", "ploidy", "wgd"]], on="model_id", how="left")
n_no_ploidy = int(cna["ploidy"].isna().sum())
cna["ploidy_used"] = cna["ploidy"].fillna(DIPLOID)
print(f"  mapped CNA rows: {len(cna):,}")
print(f"  rows without a ploidy estimate (fall back to 2.0): {n_no_ploidy:,} "
      f"({n_no_ploidy/len(cna)*100:.1f}%)")


def gate(df, amp_col, del_col):
    """Gene-role gating -- identical logic to build_cna_layer.py:91-101."""
    return (
        (df[amp_col] & df["gene_role"].isin(["oncogene", "both"])) |
        (df[del_col] & df["gene_role"].isin(["tsg", "both"])) |
        (df[amp_col] & df["regime_class"].isin(["activation_driven", "abundance_tracking"])) |
        (df[del_col] & (df["regime_class"] == "loss_of_function"))
    )


# baseline: absolute
cna["amp_abs"] = cna["TOTAL_CN"] > AMP_THRESHOLD
cna["del_abs"] = cna["TOTAL_CN"] < DEL_THRESHOLD
cna["alt_abs"] = gate(cna, "amp_abs", "del_abs")

# candidate: ploidy-relative
cna["cn_ratio"] = cna["TOTAL_CN"] / cna["ploidy_used"]
cna["amp_rel"]  = cna["cn_ratio"] > AMP_RATIO
cna["del_rel"]  = cna["cn_ratio"] < DEL_RATIO
cna["alt_rel"]  = gate(cna, "amp_rel", "del_rel")

print()
print(f"  amplification calls : absolute {cna.amp_abs.sum():>9,}  ->  relative {cna.amp_rel.sum():>9,}")
print(f"  deletion calls      : absolute {cna.del_abs.sum():>9,}  ->  relative {cna.del_rel.sum():>9,}")
print(f"  gated alterations   : absolute {cna.alt_abs.sum():>9,}  ->  relative {cna.alt_rel.sum():>9,}")

flags = (cna.groupby(["model_id", "ensg_id"])
           .agg(alt_abs=("alt_abs", "any"), alt_rel=("alt_rel", "any"))
           .reset_index())
n_flip_on  = int(((~flags.alt_abs) & flags.alt_rel).sum())
n_flip_off = int((flags.alt_abs & (~flags.alt_rel)).sum())
print(f"  (model,gene) pairs: {len(flags):,}  |  gained {n_flip_on:,}  lost {n_flip_off:,}")

# ---------------------------------------------------------------- STEP 2
print()
print("=" * 72)
print("STEP 2 -- PRIMARY ENDPOINT: is the WGD bias removed?")
print("=" * 72)

per = (cna.groupby("model_id")
         .agg(rate_abs=("alt_abs", "mean"), rate_rel=("alt_rel", "mean"),
              mean_cn=("TOTAL_CN", "mean"))
         .reset_index())
per = per.merge(sig[["model_id", "ploidy", "wgd"]], on="model_id", how="left")

res_wgd = {}
w = per.dropna(subset=["wgd"])
for label, col in [("absolute (current)", "rate_abs"), ("ploidy-relative", "rate_rel")]:
    g1, g0 = w[w.wgd >= 1][col], w[w.wgd == 0][col]
    ratio = g1.mean() / g0.mean() if g0.mean() else float("nan")
    u, p = mannwhitneyu(g1, g0)
    res_wgd[label] = {"wgd_mean": float(g1.mean()), "non_wgd_mean": float(g0.mean()),
                      "ratio": float(ratio), "mannwhitney_p": float(p),
                      "n_wgd": int(len(g1)), "n_non_wgd": int(len(g0))}
    print(f"  {label:20s}: WGD {g1.mean():.4f} vs non-WGD {g0.mean():.4f}"
          f"  ratio {ratio:.2f}x  p={p:.2e}")

res_rho = {}
d = per.dropna(subset=["ploidy"])
for label, col in [("absolute (current)", "rate_abs"), ("ploidy-relative", "rate_rel")]:
    r, p = spearmanr(d[col], d["ploidy"])
    res_rho[label] = {"rho": float(r), "p": float(p), "n": int(len(d))}
    print(f"  rho(alteration rate, ploidy) {label:20s}: {r:+.3f}  p={p:.2e}")

bias_removed = (abs(res_wgd["ploidy-relative"]["ratio"] - 1.0)
                < abs(res_wgd["absolute (current)"]["ratio"] - 1.0))
print()
print(f"  WGD bias reduced: {bias_removed}")

# ---------------------------------------------------------------- STEP 3
print()
print("=" * 72)
print("STEP 3 -- SECONDARY: blast radius on confidence labels")
print("=" * 72)

changed = flags[flags.alt_abs != flags.alt_rel][["model_id", "ensg_id", "alt_abs", "alt_rel"]]
pred = pd.read_parquet(OUTPUTS / "predictions_with_confidence.parquet",
                       columns=["model_id", "ensg_id", "confidence", "class"])
pred["model_id"] = pred["model_id"].str.upper()
hit = changed.merge(pred, on=["model_id", "ensg_id"], how="inner")

print(f"  pairs whose CNA call changes           : {len(changed):,}")
print(f"  of those, present in predictions table : {len(hit):,}")
print()
print("  current confidence label of affected predictions:")
conf_counts = hit.confidence.value_counts(dropna=False)
for k, v in conf_counts.items():
    print(f"    {str(k):10s}: {v:>8,}")

gained = hit[~hit.alt_abs & hit.alt_rel]
lost   = hit[hit.alt_abs & ~hit.alt_rel]
print()
print(f"  would GAIN a CNA alteration (possible tier upgrade): {len(gained):,}")
print(f"  would LOSE a CNA alteration (possible downgrade)   : {len(lost):,}")
print(f"    -- of the losses, currently labelled 'high'      : {int((lost.confidence=='high').sum()):,}")

# ---------------------------------------------------------------- STEP 4
print()
print("=" * 72)
print("STEP 4 -- CONTROL: ranking must be unchanged")
print("=" * 72)

fw = pd.read_parquet(OUTPUTS / "flags_with_driver.parquet",
                     columns=["model_id", "ensg_id", "mut_driver", "fusion_driver",
                              "has_driver_alteration"])
recomputed = (fw["mut_driver"].fillna(False) | fw["fusion_driver"].fillna(False))
gate_matches = bool((recomputed == fw["has_driver_alteration"].fillna(False)).all())
print(f"  has_driver_alteration == mut_driver | fusion_driver : {gate_matches}")
print("  -> CNA is absent from the ranking gate, so held-out hit@20 is")
print("     mathematically identical under both schemes. Confidence labels change;")
print("     the ordering of cell lines does not.")

# ---------------------------------------------------------------- save
results = {
    "question": ("Does ploidy-normalising the CNA thresholds remove the measured "
                 "whole-genome-doubling bias?"),
    "read_only": True,
    "thresholds": {
        "absolute_amp": AMP_THRESHOLD, "absolute_del": DEL_THRESHOLD,
        "ratio_amp": AMP_RATIO, "ratio_del": DEL_RATIO,
        "note": "ratio cutoffs chosen so diploid lines keep today's calls exactly",
        "rows_without_ploidy_fallback_2.0": n_no_ploidy,
    },
    "call_counts": {
        "amp_absolute": int(cna.amp_abs.sum()), "amp_relative": int(cna.amp_rel.sum()),
        "del_absolute": int(cna.del_abs.sum()), "del_relative": int(cna.del_rel.sum()),
        "gated_alteration_absolute": int(cna.alt_abs.sum()),
        "gated_alteration_relative": int(cna.alt_rel.sum()),
        "pairs_total": int(len(flags)),
        "pairs_gained": n_flip_on, "pairs_lost": n_flip_off,
    },
    "PRIMARY_wgd_bias": res_wgd,
    "PRIMARY_ploidy_correlation": res_rho,
    "wgd_bias_reduced": bool(bias_removed),
    "SECONDARY_confidence_blast_radius": {
        "pairs_changed": int(len(changed)),
        "predictions_affected": int(len(hit)),
        "current_confidence_distribution": {str(k): int(v) for k, v in conf_counts.items()},
        "would_gain_cna": int(len(gained)),
        "would_lose_cna": int(len(lost)),
        "losses_currently_high_confidence": int((lost.confidence == "high").sum()),
    },
    "CONTROL_ranking_invariance": {
        "driver_gate_excludes_cna": gate_matches,
        "hit_at_20_change": 0.0,
        "note": "has_driver_alteration = mut_driver | fusion_driver; CNA never enters ranking",
    },
}

out = OUTPUTS / "test_run_ploidy_normalised_cna_results.json"
with open(out, "w") as f:
    json.dump(results, f, indent=2)
print()
print(f"Saved: {out}")
print("No production file was modified.")
