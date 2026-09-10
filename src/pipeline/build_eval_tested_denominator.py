"""
build_eval_tested_denominator.py
--------------------------------
Corrected evaluation denominator: score only (gene, line) pairs that were
actually TESTED in GDSC2.

THE DEFECT
The shipped harness left-joins GDSC sensitivity onto every scored pair and then
applies `.fillna(False)`:

    stage4_eval_save.py:41            core_t["is_sensitive"].fillna(False)
    06_held_out_eval.ipynb cell 4     same
    build_gene_regime.py:150          is_sensitive.astype("boolean").fillna(False)

Every (gene, line) pair GDSC2 never tested therefore enters the denominator as a
CONFIRMED NEGATIVE. Measured: the median gene was tested on 931 cell lines but
scored against ~1,563, so roughly 630 untested lines per gene were counted as
true negatives. Apparent prevalence collapses from 41.2% (among tested) to
~14.6% (among all scored), and every @k figure and its chance baseline is
computed on that inflated negative class.

This is not a small correction. It changes the hypergeometric baseline for
any-hit@20 from ~0.91 to ~1.00, which is what makes the previously reported
"93.6% top-20" unambiguously a chance result rather than merely a marginal one.

THE FIX
Inner-join instead of left-join. A pair with no GDSC2 row is `not_tested` and is
excluded from both numerator and denominator -- it is not a negative. This is the
three-state rule (docs/BUILD_SPEC.md C3) applied to the label column.

WHAT THIS DOES NOT AFFECT
The gate result (gate-region pAUC 0.596-0.608) is measured against CRISPR
dependency labels (Chronos / Project Score), not GDSC2, and does not pass through
this denominator. It is unaffected.

OUTPUT
    src/pipeline/outputs/eval_tested_denominator.parquet   per-gene, both denominators
    src/pipeline/outputs/eval_tested_denominator.json      headline comparison
"""
import json
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.special import gammaln

OUTPUTS = Path("src/pipeline/outputs")
VAL = Path("validation/prepared")
SEED = 42
KS = (1, 5, 20)

t0 = time.time()

core = pd.read_parquet(OUTPUTS / "core_score.parquet",
                       columns=["model_id", "ensg_id", "core_score"])
core["model_id"] = core["model_id"].str.lower()
regime = pd.read_parquet(OUTPUTS / "gene_regime.parquet", columns=["ensg_id"])
gdsc = pd.read_parquet(VAL / "gdsc_scored_ready.parquet",
                       columns=["model_id", "target_ensg", "sensitive"])
gdsc["model_id"] = gdsc["model_id"].str.lower()
gdsc["ensg_id"] = gdsc["target_ensg"].astype("string").str.lower()
print(f"[{time.time()-t0:.1f}s] core {len(core):,} pairs; gdsc {len(gdsc):,} rows, "
      f"{gdsc.ensg_id.nunique()} genes, {gdsc.model_id.nunique()} lines")

# one row per (gene, line): tested, and whether sensitive in any assay
tested = (gdsc.groupby(["ensg_id", "model_id"], as_index=False)["sensitive"]
          .max().rename(columns={"sensitive": "is_sensitive"}))
tested["tested"] = True
print(f"[{time.time()-t0:.1f}s] distinct tested (gene, line) pairs: {len(tested):,}")

# the same held-out split the shipped harness uses
curated = sorted(set(regime.ensg_id) & set(gdsc.ensg_id))
rng = np.random.default_rng(SEED)
test_genes = set(rng.choice(curated, size=max(1, len(curated) // 5), replace=False))
print(f"[{time.time()-t0:.1f}s] curated {len(curated)}, held out {len(test_genes)}")

c = core[core.ensg_id.isin(curated)].merge(tested, on=["ensg_id", "model_id"], how="left")
c["tested"] = c["tested"].fillna(False).infer_objects(copy=False).astype(bool)
c["is_sensitive_broken"] = c["is_sensitive"].fillna(False).infer_objects(copy=False).astype(bool)


def logC(n, k):
    return gammaln(n + 1) - gammaln(k + 1) - gammaln(n - k + 1)


def per_gene(df, label):
    """df already restricted to the rows that form the denominator."""
    out = []
    for g, d in df.groupby("ensg_id"):
        N = len(d)
        K = int(d["y"].sum())
        if K == 0 or N - K == 0:
            continue
        r = d["core_score"].rank(ascending=False, method="first")
        rec = {"ensg_id": g, "denominator": label, "N": N, "K": K,
               "prevalence": K / N}
        for k in KS:
            hits = int(d.loc[r <= k, "y"].sum())
            rec[f"recall@{k}"] = hits / K
            rec[f"anyhit@{k}"] = float(hits >= 1)
            rec[f"base_recall@{k}"] = min(k, N) / N
            rec[f"base_anyhit@{k}"] = 1.0 - (np.exp(logC(N - K, k) - logC(N, k))
                                             if N - K >= k else 0.0)
        out.append(rec)
    return pd.DataFrame(out)


broken = c.assign(y=c["is_sensitive_broken"])
fixed = c[c["tested"]].assign(y=c.loc[c["tested"], "is_sensitive"].astype(bool))
PB = per_gene(broken, "broken (untested = negative)")
PF = per_gene(fixed, "fixed (tested pairs only)")
P = pd.concat([PB, PF], ignore_index=True)
P["split"] = np.where(P.ensg_id.isin(test_genes), "held-out", "tuning")


def boot(vals, base, n=10000, seed=SEED):
    rg = np.random.default_rng(seed)
    d = np.asarray(vals) - np.asarray(base)
    if len(d) == 0:
        return (np.nan, np.nan)
    m = d[rg.integers(0, len(d), size=(n, len(d)))].mean(axis=1)
    return tuple(np.percentile(m, [2.5, 97.5]))


print(f"\n[{time.time()-t0:.1f}s] denominator sizes, per gene (median):")
for lab, sub in P.groupby("denominator"):
    print(f"  {lab:32s} N={sub.N.median():7.0f}  K={sub.K.median():6.0f}  "
          f"prevalence={sub.prevalence.median():.3f}")

rows = []
print("\n" + "=" * 104)
print(f"{'denominator':32s} {'split':9s} {'n':>4s} {'metric':11s} {'obs':>8s} "
      f"{'chance':>8s} {'excess':>9s}  95% CI")
print("=" * 104)
for lab in ["broken (untested = negative)", "fixed (tested pairs only)"]:
    for split in ["tuning", "held-out"]:
        s = P[(P.denominator == lab) & (P.split == split)]
        if not len(s):
            continue
        for k in KS:
            for nm in ("recall", "anyhit"):
                o = s[f"{nm}@{k}"].mean()
                b = s[f"base_{nm}@{k}"].mean()
                lo, hi = boot(s[f"{nm}@{k}"], s[f"base_{nm}@{k}"])
                star = "" if (lo <= 0 <= hi) else "  *"
                print(f"{lab:32s} {split:9s} {len(s):4d} {nm+'@'+str(k):11s} "
                      f"{o:8.4f} {b:8.4f} {o-b:+9.4f}  [{lo:+.4f}, {hi:+.4f}]{star}")
                rows.append({"denominator": lab, "split": split, "n_genes": len(s),
                             "metric": f"{nm}@{k}", "observed": float(o),
                             "chance": float(b), "excess": float(o - b),
                             "ci_lo": float(lo), "ci_hi": float(hi),
                             "significant": not (lo <= 0 <= hi)})
print("* = 95% paired-bootstrap CI on (observed - chance) excludes zero")

P.to_parquet(OUTPUTS / "eval_tested_denominator.parquet", index=False)
summary = {
    "defect": ("untested (gene, line) pairs entered the denominator as confirmed "
               "negatives via .fillna(False)"),
    "defect_sites": ["src/pipeline/stage4_eval_save.py:41",
                     "src/pipeline/06_held_out_eval.ipynb cell 4",
                     "src/pipeline/build_gene_regime.py:150"],
    "seed": SEED, "n_curated": len(curated), "n_held_out": len(test_genes),
    "median_N_broken": float(PB.N.median()), "median_N_fixed": float(PF.N.median()),
    "median_prevalence_broken": float(PB.prevalence.median()),
    "median_prevalence_fixed": float(PF.prevalence.median()),
    "results": rows,
    "unaffected": ("gate-region pAUC is measured against CRISPR dependency labels, "
                   "not GDSC2, and does not pass through this denominator"),
}
with open(OUTPUTS / "eval_tested_denominator.json", "w") as fh:
    json.dump(summary, fh, indent=2)
print(f"\n[{time.time()-t0:.1f}s] Written: {OUTPUTS/'eval_tested_denominator.parquet'}")
print(f"[{time.time()-t0:.1f}s] Written: {OUTPUTS/'eval_tested_denominator.json'}")
