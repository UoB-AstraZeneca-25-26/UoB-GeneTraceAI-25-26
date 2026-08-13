"""
test_run_driver_gate_ci.py -- confidence intervals on the driver-gate lift, and
the abundance_tracking bypass diagnostic.

Two questions the shipped stage4_eval_summary.json leaves open:

  1. The headline "26.45% driver lift" is a ratio of two means over n=28 genes
     with no interval attached. Does the lift survive a paired bootstrap?
  2. `sort_driver` bypasses the gate entirely for abundance_tracking genes
     (stage4_eval_save.py:45-47), so their delta is 0.0 by construction. Are
     those genes genuinely abundance-driven, or is the bypass hiding a gate
     that would have helped?

Read-only. Bootstraps the per-gene table the eval already wrote rather than
re-deriving it, so the numbers are the shipped ones by construction.

Bootstrap unit is the GENE, matching how the hit rate is aggregated (mean over
per-gene hit rates). Resampling pairs within a gene would understate the
variance that actually matters, which is gene-to-gene.

N_BOOT / SEED match the rest of the test_run_* suite.

Output: src/pipeline/outputs/test_run_driver_gate_ci_results.json

Usage:
    python src/pipeline/test_run_driver_gate_ci.py
"""
from __future__ import annotations

import json
from pathlib import Path

import duckdb
import numpy as np
import pandas as pd

OUTPUTS = Path("src/pipeline/outputs")
DB = OUTPUTS / "celllineselector.db"

N_BOOT = 2000
SEED = 42
CI = (2.5, 97.5)

res: dict = {}

per = pd.read_parquet(OUTPUTS / "stage4_eval_consistent_denom.parquet")
print(f"per-gene eval table: {per.shape[0]} genes "
      f"({(per['class'] == 'activation_driven').sum()} activation_driven, "
      f"{(per['class'] == 'abundance_tracking').sum()} abundance_tracking)")

# ---------------------------------------------------------------- Step 1
print()
print("=" * 72)
print("STEP 1 -- paired bootstrap over genes")
print("=" * 72)

rng = np.random.default_rng(SEED)
n = len(per)
flat = per.hr_flat.to_numpy()
driver = per.hr_driver.to_numpy()
anyalt = per.hr_any.to_numpy()

boot = {"delta_driver": [], "delta_any": [], "driver_vs_any": [], "lift_pct": []}
for _ in range(N_BOOT):
    idx = rng.integers(0, n, n)
    f, d, a = flat[idx].mean(), driver[idx].mean(), anyalt[idx].mean()
    boot["delta_driver"].append(d - f)
    boot["delta_any"].append(a - f)
    boot["driver_vs_any"].append(d - a)
    boot["lift_pct"].append((d - f) / f * 100 if f > 0 else np.nan)

for k, v in boot.items():
    v = np.asarray(v, dtype=float)
    v = v[np.isfinite(v)]
    lo, hi = np.percentile(v, CI)
    point = {
        "delta_driver": driver.mean() - flat.mean(),
        "delta_any": anyalt.mean() - flat.mean(),
        "driver_vs_any": driver.mean() - anyalt.mean(),
        "lift_pct": (driver.mean() - flat.mean()) / flat.mean() * 100,
    }[k]
    straddles = bool(lo <= 0 <= hi)
    res.setdefault("bootstrap", {})[k] = {
        "point": round(float(point), 6),
        "ci_lo": round(float(lo), 6),
        "ci_hi": round(float(hi), 6),
        "straddles_zero": straddles,
        "n_boot": N_BOOT,
    }
    flag = "  <-- STRADDLES ZERO" if straddles else ""
    print(f"  {k:15s} point={point:+.6f}  95% CI [{lo:+.6f}, {hi:+.6f}]{flag}")

# sign test on the per-gene deltas -- distribution-free companion
d = per.delta_driver.to_numpy()
res["sign_test"] = {
    "genes_improved": int((d > 0).sum()),
    "genes_unchanged": int((d == 0).sum()),
    "genes_worsened": int((d < 0).sum()),
    "n": int(len(d)),
}
print(f"\n  per-gene sign: {int((d > 0).sum())} improved, "
      f"{int((d == 0).sum())} unchanged, {int((d < 0).sum())} worsened (n={len(d)})")

# ---------------------------------------------------------------- Step 2
print()
print("=" * 72)
print("STEP 2 -- abundance_tracking bypass diagnostic")
print("=" * 72)

abund = per[per["class"] == "abundance_tracking"]
con = duckdb.connect(str(DB), read_only=True)
sym = con.execute("""
    select lower(gene_id) as ensg_id, hugo_symbol, gene_role, cgc_tier
    from main.gene_enriched
""").df().set_index("ensg_id")
con.close()

core = pd.read_parquet(OUTPUTS / "core_score.parquet",
                       columns=["model_id", "ensg_id", "core_score"])
core["model_id"] = core["model_id"].str.lower()
core["ensg_id"] = core["ensg_id"].astype("string").str.lower()
flags = pd.read_parquet(OUTPUTS / "flags_with_driver.parquet",
                        columns=["model_id", "ensg_id", "has_driver_alteration"])
flags["model_id"] = flags["model_id"].str.lower()
flags["ensg_id"] = flags["ensg_id"].astype("string").str.lower()
gdsc = pd.read_parquet("validation/prepared/gdsc_scored_ready.parquet",
                       columns=["model_id", "target_ensg", "sensitive"])
gdsc["model_id"] = gdsc["model_id"].str.lower()
gdsc["target_ensg"] = gdsc["target_ensg"].astype("string").str.lower()

rows = []
for ensg in abund.ensg_id:
    g = sym.loc[ensg] if ensg in sym.index else None
    cg = core[core.ensg_id == ensg]
    fg = flags[flags.ensg_id == ensg]
    sens = set(gdsc[(gdsc.target_ensg == ensg) & gdsc.sensitive].model_id)
    known = sens & set(cg.model_id)

    m = cg.merge(fg[["model_id", "has_driver_alteration"]], on="model_id", how="left")
    m["has_driver_alteration"] = m["has_driver_alteration"].astype("boolean").fillna(False)

    # shipped behaviour: gate bypassed, score only
    shipped = set(m.sort_values("core_score", ascending=False).head(20).model_id)
    # counterfactual: gate applied, as activation_driven genes get
    cf = set(m.sort_values(["has_driver_alteration", "core_score"],
                           ascending=[False, False]).head(20).model_id)

    hr_shipped = len(shipped & sens) / len(known) if known else float("nan")
    hr_cf = len(cf & sens) / len(known) if known else float("nan")
    n_flagged = int(m.has_driver_alteration.sum())

    rows.append({
        "ensg_id": ensg,
        "symbol": None if g is None else str(g.hugo_symbol),
        "gene_role": None if g is None else str(g.gene_role),
        "cgc_tier": None if g is None or pd.isna(g.cgc_tier) else str(g.cgc_tier),
        "known": int(len(known)),
        "lines_scored": int(len(m)),
        "lines_driver_flagged": n_flagged,
        "pct_driver_flagged": round(n_flagged / len(m) * 100, 2) if len(m) else None,
        "hr_shipped_bypassed": round(hr_shipped, 6),
        "hr_counterfactual_gated": round(hr_cf, 6),
        "delta_if_gate_applied": round(hr_cf - hr_shipped, 6),
    })
    print(f"  {rows[-1]['symbol']} ({ensg})")
    print(f"    role={rows[-1]['gene_role']} tier={rows[-1]['cgc_tier']} known={rows[-1]['known']}")
    print(f"    driver-flagged lines: {n_flagged:,} / {len(m):,} "
          f"({rows[-1]['pct_driver_flagged']}%)")
    print(f"    hr shipped (bypassed) = {hr_shipped:.6f}")
    print(f"    hr counterfactual (gated) = {hr_cf:.6f}   "
          f"delta = {hr_cf - hr_shipped:+.6f}")

res["abundance_bypass"] = rows

p = OUTPUTS / "test_run_driver_gate_ci_results.json"
p.write_text(json.dumps(res, indent=2), encoding="utf-8")
print(f"\nSaved: {p}")
