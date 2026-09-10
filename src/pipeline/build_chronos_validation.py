"""
build_chronos_validation.py
----------------------------
Independent validation check for "unknown"-class genes (no curated GDSC
regime, no COSMIC oncogene/TSG role) using genome-wide CRISPR essentiality
(Chronos) as ground truth for whether abundance (core_score) actually
tracks functional dependency for that gene.

This does NOT feed Chronos into the ranking score itself (that was tried
and rejected earlier -- BCL2's hit@20 collapsed to 0 when Chronos was
blended into core_score as a 3rd layer). This is a read-only diagnostic:
for each gene, does the cell-line ordering implied by core_score agree
with the cell-line ordering implied by real essentiality data?

Per-gene classification (Spearman rho between core_score and Chronos
essentiality percentile, across cell lines both cover):
    validated      : rho >=  0.30, p < 0.05 -- abundance ranking tracks real
                     dependency strongly enough to act on
    inverted       : rho <= -0.30, p < 0.05 -- abundance ranking is likely
                     backwards for this gene (e.g. low abundance = high
                     dependency, a threshold/buffering effect)
    weak_positive  : 0.10 <= rho < 0.30, p < 0.05 -- direction agrees but the
    weak_negative  : -0.30 < rho <= -0.10, p < 0.05 -- effect is too small to
                     carry a confidence upgrade or a contradiction warning
    none           : no significant relationship -- stays honestly unknown,
                     but a CHECKED unknown rather than an assumed one

WHY THE EFFECT FLOOR MOVED FROM 0.10 TO 0.30
The old floor was |rho| > 0.10, i.e. ~1% of shared variance. At the median
overlap of n=823 cell lines, p < 0.05 is reached by |rho| >= 0.068 -- BELOW the
rho floor -- so the significance test was entirely redundant and 0.10 was the
only active criterion. The consequence was visible in the outputs:

  - 769 genes were labelled "validated" with median rho 0.128 (IQR 0.111-0.165).
    Almost the whole set sat between 0.10 and 0.17.
  - 1,258 genes carried the user-facing "*** CONTRADICTED ***" banner on rho
    between -0.10 and -0.17.

Both claims are far stronger than ~2% of variance supports, and the failure mode
is systematic rather than random: any gene essential in nearly every cell line
(housekeeping and core-essential genes) will show a small but reliably non-zero
abundance-essentiality correlation purely because n is large. GAPDH passed at
rho=0.138, p=7e-5, which upgraded its confidence to "moderate" across all 1,485
of its rows -- leaving the pipeline more confident about a housekeeping gene than
about MUC1 (rho=-0.057, p=0.10, "unknown"). That is exactly backwards.

0.30 is Cohen's conventional floor for a medium correlation (~9% of variance).
It is a deliberately conservative line for a claim as strong as "safe to trust"
or "the ranking is backwards": validated 769 -> 29 genes, inverted 1,258 -> 14.
The 0.10-0.30 band is not discarded -- it is retained as weak_positive /
weak_negative so the direction is still recorded, it simply no longer moves a
confidence label or raises a contradiction warning.

Role in pipeline: read by build_full_predictions.py to adjust confidence
tier and signal_quality for unknown-class genes only. Curated (measured)
and COSMIC-prior genes are untouched by this file.

Output: src/pipeline/outputs/chronos_validation.parquet
    (ensg_id, n, rho, pval, chronos_check)
"""
import time
import numpy as np
import pandas as pd
from pathlib import Path
from scipy.stats import t as tdist

OUTPUTS = Path("src/pipeline/outputs")
VAL     = Path("validation/prepared")

# Effect-size floor for a claim strong enough to move a confidence label or raise
# a contradiction warning. Cohen's medium correlation (~9% shared variance).
RHO_THRESHOLD      = 0.30
# Floor for recording a direction without acting on it (Cohen's small).
RHO_WEAK_THRESHOLD = 0.10
P_THRESHOLD        = 0.05
MIN_N              = 30

t0 = time.time()

pred = pd.read_parquet(OUTPUTS / "predictions_with_confidence.parquet",
                        columns=["model_id", "ensg_id", "core_score", "class"])
pred["model_id"] = pred["model_id"].str.lower()
unknown = pred[pred["class"] == "unknown"][["model_id", "ensg_id", "core_score"]]
print(f"[{time.time()-t0:.1f}s] unknown-class rows: {len(unknown):,}  genes: {unknown.ensg_id.nunique():,}")

# T12 F1 fix: read real Chronos (from Achilles HDF5 via diagnostics/T12_F1_chronos_fix.py)
# instead of the old chronos_long.parquet (which was Project Score negated, not Chronos)
ch = pd.read_parquet(OUTPUTS / "chronos_corrected.parquet")
# model_id is already lowercase ACH- format; ensg_id is lowercase ENSG
ch["ensg_id"]  = ch["ensg_id"].str.lower()
# More negative raw essentiality = more essential -> invert so high chronos_pct = essential
ch["chronos_pct"] = 1.0 - ch.groupby("ensg_id")["essentiality"].rank(pct=True, method="average")
print(f"[{time.time()-t0:.1f}s] chronos_pct built for {ch.ensg_id.nunique():,} genes")

merged = unknown.merge(ch[["model_id", "ensg_id", "chronos_pct"]], on=["model_id", "ensg_id"], how="inner")
print(f"[{time.time()-t0:.1f}s] merged: {len(merged):,} rows, {merged.ensg_id.nunique():,} genes with any overlap")

# Vectorised Spearman (Pearson-on-ranks via groupby aggregation -- no per-gene python loop)
merged["r1"] = merged.groupby("ensg_id")["core_score"].rank(method="average")
merged["r2"] = merged.groupby("ensg_id")["chronos_pct"].rank(method="average")

g = merged.groupby("ensg_id")
merged["dx"]   = merged["r1"] - g["r1"].transform("mean")
merged["dy"]   = merged["r2"] - g["r2"].transform("mean")
merged["dxdy"] = merged["dx"] * merged["dy"]
merged["dx2"]  = merged["dx"] ** 2
merged["dy2"]  = merged["dy"] ** 2

agg = merged.groupby("ensg_id").agg(
    n=("r1", "size"),
    sum_dxdy=("dxdy", "sum"),
    sum_dx2=("dx2", "sum"),
    sum_dy2=("dy2", "sum"),
)
agg["rho"] = agg["sum_dxdy"] / np.sqrt(agg["sum_dx2"] * agg["sum_dy2"])
agg = agg[agg["n"] >= MIN_N].copy()

with np.errstate(invalid="ignore", divide="ignore"):
    tstat = agg["rho"] * np.sqrt((agg["n"] - 2) / (1 - agg["rho"] ** 2))
agg["pval"] = 2 * tdist.sf(np.abs(tstat), df=agg["n"] - 2)

# Fisher z 95% CI on rho -- reported, not gated on. Makes the precision of each
# per-gene estimate visible instead of collapsing it to a pass/fail label.
with np.errstate(invalid="ignore", divide="ignore"):
    _z  = np.arctanh(agg["rho"].clip(-0.999999, 0.999999))
    _se = 1.0 / np.sqrt(agg["n"] - 3)
    agg["rho_ci_lo"] = np.tanh(_z - 1.96 * _se)
    agg["rho_ci_hi"] = np.tanh(_z + 1.96 * _se)

_sig = agg.pval < P_THRESHOLD
agg["chronos_check"] = np.select(
    [
        _sig & (agg.rho >=  RHO_THRESHOLD),
        _sig & (agg.rho <= -RHO_THRESHOLD),
        _sig & (agg.rho >=  RHO_WEAK_THRESHOLD),
        _sig & (agg.rho <= -RHO_WEAK_THRESHOLD),
    ],
    ["validated", "inverted", "weak_positive", "weak_negative"],
    default="none",
)
agg = agg.reset_index()[["ensg_id", "n", "rho", "pval",
                         "rho_ci_lo", "rho_ci_hi", "chronos_check"]]

print(f"[{time.time()-t0:.1f}s] classified {len(agg):,} genes "
      f"(effect floor |rho| >= {RHO_THRESHOLD}, weak band >= {RHO_WEAK_THRESHOLD})")
print(agg["chronos_check"].value_counts().to_string())
print("\nrho by band:")
print(agg.groupby("chronos_check")["rho"].agg(["count", "min", "median", "max"]).round(3).to_string())

agg.to_parquet(OUTPUTS / "chronos_validation.parquet", index=False)
print(f"\n[{time.time()-t0:.1f}s] Written: {OUTPUTS / 'chronos_validation.parquet'}  ({len(agg):,} rows)")
