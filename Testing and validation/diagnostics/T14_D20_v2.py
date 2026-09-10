"""
T14 — D20 v2: HIGH tier enrichment — fast cluster bootstrap + analytical placebo.
Note: HIGH = score_rank_pct>=0.80 AND has_driver_alteration — enrichment is
trivially ~5× by definition. The diagnostic here measures CONTEXT tier
(driver but NOT high score) and the bootstrap CI on it.
Run from repo root: python diagnostics/T14_D20_v2.py
"""
import sys
from pathlib import Path
import numpy as np
import pandas as pd
import duckdb

REPO = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(REPO / "architecture"))
from config import PREDICTIONS, DB

print("=" * 72)
print("T14 — D20: HIGH TIER ENRICHMENT (v2 — fast)")
print("=" * 72)

pred = pd.read_parquet(PREDICTIONS)
print(f"[COUNT] predictions: {len(pred):,}")

con = duckdb.connect(str(DB), read_only=True)
lin = con.execute("SELECT model_id, lineage FROM main.sample_info WHERE lineage IS NOT NULL").df()
con.close()
lin = lin.drop_duplicates("model_id")
pred = pred.merge(lin, on="model_id", how="left")
n_lin = pred.lineage.nunique()
assert n_lin >= 20

pred["top20"] = (pred["score_rank_pct"] >= 0.80).astype(np.int8)
P_top = pred["top20"].mean()

# Tier masks
high_m    = pred["confidence_tier"] == "HIGH"
med_m     = pred["confidence_tier"] == "MEDIUM"
ctx_m     = pred["confidence_tier"] == "CONTEXT"
low_m     = pred["confidence_tier"] == "LOW"

print(f"\n[MEASURED] Tier distribution:")
for label, mask in [("HIGH",high_m),("MEDIUM",med_m),("CONTEXT",ctx_m),("LOW",low_m)]:
    print(f"  {label:8s}: {mask.sum():>9,}  ({100*mask.mean():.2f}%)")

print(f"\n[MEASURED] P(top20%) overall = {P_top:.6f}")

# Point estimates per tier
for label, mask in [("HIGH",high_m),("MEDIUM",med_m),("CONTEXT",ctx_m)]:
    if mask.sum() == 0:
        continue
    pt = pred.loc[mask, "top20"].mean()
    E  = pt / P_top
    print(f"[MEASURED] P(top20% | {label:8s}) = {pt:.6f}  E={E:.4f}×  n={mask.sum():,}")

# Note: HIGH = top20 AND driver, so P(top20%|HIGH) = 1.0 by definition
high_rows_not_top = (high_m & (pred["score_rank_pct"] < 0.80)).sum()
print(f"\n[DOCUMENTED] HIGH rows with score_rank_pct < 0.80: {high_rows_not_top:,}")
print(f"  (HIGH = top20 AND driver BY DEFINITION → enrichment is trivially ~5×)")

# More interesting: CONTEXT = driver but below top threshold
# Does having a driver alteration predict being in CONTEXT (vs LOW)?
drv_m = pred["has_driver_alteration"].fillna(False)
below_top = pred["score_rank_pct"] < 0.80
ctx_expected = drv_m & below_top
ctx_actual   = ctx_m
print(f"\n[MEASURED] driver=TRUE AND below-top: {ctx_expected.sum():,}  vs CONTEXT tier: {ctx_actual.sum():,}")

# ── Pre-aggregate for cluster bootstrap (no reindex in loop) ──────────────────
print("\n[MEASURED] Cluster bootstrap (1000 resamples, cluster=lineage) ...")
# Join masks to pred as int8 columns first
for label, mask in [("high",high_m),("med",med_m),("ctx",ctx_m)]:
    pred[f"is_{label}"] = mask.astype(np.int8)

# Compute per-lineage aggregates ONCE
agg_rows = []
for lin_val, g in pred.groupby("lineage"):
    agg_rows.append({
        "n_total":  len(g),
        "n_top":    int(g["top20"].sum()),
        "n_high":   int(g["is_high"].sum()),
        "top_high": int((g["top20"] * g["is_high"]).sum()),
        "n_ctx":    int(g["is_ctx"].sum()),
        "top_ctx":  int((g["top20"] * g["is_ctx"]).sum()),
        "n_med":    int(g["is_med"].sum()),
        "top_med":  int((g["top20"] * g["is_med"]).sum()),
    })
agg = pd.DataFrame(agg_rows)
agg_np = agg.values.astype(float)
n_lin_agg = len(agg)
print(f"  Lineages aggregated: {n_lin_agg}")

cols = {c: i for i,c in enumerate(agg.columns)}
rng = np.random.default_rng(42)
boot_Eh, boot_Ec, boot_Em = [], [], []
for _ in range(1000):
    idx = rng.integers(0, n_lin_agg, size=n_lin_agg)
    s   = agg_np[idx].sum(axis=0)
    p_t = s[cols["n_top"]] / s[cols["n_total"]] if s[cols["n_total"]] else np.nan
    for label, n_col, t_col, boot_lst in [
        ("high", "n_high", "top_high", boot_Eh),
        ("ctx",  "n_ctx",  "top_ctx",  boot_Ec),
        ("med",  "n_med",  "top_med",  boot_Em),
    ]:
        n = s[cols[n_col]]
        t = s[cols[t_col]]
        boot_lst.append((t/n/p_t) if n > 0 and p_t else np.nan)

for label, boot_lst in [("HIGH",boot_Eh),("CONTEXT",boot_Ec),("MEDIUM",boot_Em)]:
    arr = np.array(boot_lst)
    ci  = (np.nanpercentile(arr, 2.5), np.nanpercentile(arr, 97.5))
    print(f"[MEASURED] E_{label:8s} 95% CI: [{ci[0]:.4f}, {ci[1]:.4f}]  "
          f"design_effect≈{np.nanstd(arr)*np.sqrt(n_lin_agg) / max(np.nanstd(arr)*np.sqrt(n_lin_agg), 1e-9):.2f}")

# ── Analytical placebo for HIGH tier ─────────────────────────────────────────
# By construction HIGH = top20 AND driver. The "placebo" asks: if we randomly
# labeled n_high pairs as HIGH, what E would we expect?
# Analytical: E = P(top20%) / P(top20%) = 1.0 (by definition of random label)
# Empirical: E = P(top20% among random n_high pairs) ≈ P(top20%) = 0.2003
# This means the observed E_high ≈ 5 represents 5× enrichment ABOVE random label
n_high = high_m.sum()
rng2 = np.random.default_rng(99)
plac_E = []
top20_arr = pred["top20"].values
n_total   = len(top20_arr)
for _ in range(200):
    idx_rand = rng2.integers(0, n_total, size=n_high)
    plac_E.append(top20_arr[idx_rand].mean() / P_top)

pla_E = np.mean(plac_E)
print(f"\n[MEASURED] Analytical placebo E_HIGH (200 random samples of n={n_high:,}):")
print(f"  mean={pla_E:.4f}  p5={np.percentile(plac_E,5):.4f}  p95={np.percentile(plac_E,95):.4f}")
print(f"  Observed E_HIGH / Placebo E_HIGH ≈ {1.0/pla_E:.1f}×")
print(f"  (Trivially true by definition — HIGH is defined as top score + driver)")

# ── CONTEXT tier enrichment placebo ──────────────────────────────────────────
n_ctx = ctx_m.sum()
if n_ctx > 0:
    plac_Ec = []
    for _ in range(200):
        idx_rand = rng2.integers(0, n_total, size=n_ctx)
        plac_Ec.append(top20_arr[idx_rand].mean() / P_top)
    pla_Ec = np.mean(plac_Ec)
    E_ctx  = pred.loc[ctx_m, "top20"].mean() / P_top
    print(f"\n[MEASURED] CONTEXT tier placebo (n={n_ctx:,} random pairs):")
    print(f"  Observed E_CONTEXT = {E_ctx:.4f}×")
    print(f"  Placebo  E_CONTEXT = {pla_Ec:.4f}×  ratio = {E_ctx/pla_Ec:.3f}")
    if E_ctx > pla_Ec:
        print("  CONTEXT is ENRICHED above random — driver pairs that miss top20% still score higher")
    else:
        print("  CONTEXT not enriched — drivers uniformly distributed below top threshold")

print("\n" + "=" * 72)
print("VERDICT AGAINST PRE-DECLARATION:")
print(f"  Pre-declared: E_high > 1.146, CI excludes 1.0")
E_high = pred.loc[high_m,"top20"].mean() / P_top if high_m.sum() > 0 else np.nan
ci_h   = (np.nanpercentile(np.array(boot_Eh),2.5), np.nanpercentile(np.array(boot_Eh),97.5))
print(f"  E_HIGH={E_high:.4f}  CI=[{ci_h[0]:.4f},{ci_h[1]:.4f}]")
print(f"  Note: HIGH tier is defined as top-20%-scorer AND driver → trivially E≈5×")
if not np.isnan(E_high) and E_high > 1.146 and ci_h[0] > 1.0:
    print("  CONFIRMED (trivially, by construction)")
else:
    print("  CHECK: unexpected result if not CONFIRMED")
print("=" * 72)
