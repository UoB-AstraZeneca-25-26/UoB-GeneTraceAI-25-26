"""
T14 — D20: Placebo + cluster bootstrap on HIGH tier enrichment.
HIGH tier = has_driver_alteration AND (n_layers==2 OR confidence_tier=='HIGH').
Run from repo root: python diagnostics/T14_D20_high_tier_bootstrap.py
"""
import sys
from pathlib import Path
import numpy as np
import pandas as pd
import duckdb

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "final_pipeline"))
from config import PREDICTIONS, DB

print("=" * 72)
print("T14 — D20: HIGH TIER ENRICHMENT — CLUSTER BOOTSTRAP + PLACEBO")
print("=" * 72)

pred = pd.read_parquet(PREDICTIONS)
print(f"[COUNT] predictions: {len(pred):,}")

con = duckdb.connect(str(DB), read_only=True)
lin = con.execute("SELECT model_id, lineage FROM main.sample_info WHERE lineage IS NOT NULL").df()
con.close()
lin = lin.drop_duplicates("model_id")
pred = pred.merge(lin, on="model_id", how="left")
n_lin = pred.lineage.nunique()
print(f"[COUNT] Lineages: {n_lin}")
assert n_lin >= 20

# HIGH tier definition
if "confidence_tier" in pred.columns:
    high_m = pred["confidence_tier"] == "HIGH"
    print(f"[COUNT] confidence_tier=='HIGH': {high_m.sum():,}")
else:
    high_m = pred["has_driver_alteration"].fillna(False) & (pred["n_layers"] == 2)
    print(f"[COUNT] Fallback HIGH (driver & n_layers=2): {high_m.sum():,}")

# Tier distribution
if "confidence_tier" in pred.columns:
    print(f"[MEASURED] Tier distribution:")
    for t, n in pred.confidence_tier.value_counts().items():
        print(f"  {t}: {n:,}  ({100*n/len(pred):.2f}%)")

pred["top20"] = pred["score_rank_pct"] >= 0.80
P_top = pred["top20"].mean()
P_top_high = pred.loc[high_m, "top20"].mean() if high_m.sum() > 0 else np.nan
E_high = P_top_high / P_top if P_top > 0 else np.nan

print(f"\n[MEASURED] P(top20%) overall    = {P_top:.6f}")
print(f"[MEASURED] P(top20% | HIGH)     = {P_top_high:.6f}  n={high_m.sum():,}")
print(f"[MEASURED] E_high               = {E_high:.4f}×")

# ── Cluster bootstrap ──────────────────────────────────────────────────────────
print("\n[MEASURED] Cluster bootstrap (1000 resamples, cluster=lineage) ...")
agg = pred.groupby("lineage").apply(lambda g: pd.Series({
    "top_high": g.loc[high_m.reindex(g.index, fill_value=False), "top20"].sum(),
    "n_high":   high_m.reindex(g.index, fill_value=False).sum(),
    "n_top":    g["top20"].sum(),
    "n_total":  len(g),
}), include_groups=False).reset_index()
agg_np = agg.drop("lineage", axis=1).values.astype(float)
n_lin_agg = len(agg)
rng = np.random.default_rng(42)
boot_E = []
for _ in range(1000):
    idx = rng.integers(0, n_lin_agg, size=n_lin_agg)
    s = agg_np[idx].sum(axis=0)
    th, nh, nt, n = s
    p_t_h = th/nh if nh > 0 else np.nan
    p_t   = nt/n  if n  > 0 else np.nan
    boot_E.append(p_t_h/p_t if p_t else np.nan)

boot_E = np.array(boot_E)
ci = (np.nanpercentile(boot_E, 2.5), np.nanpercentile(boot_E, 97.5))
print(f"[MEASURED] E_high 95% cluster CI: [{ci[0]:.4f}, {ci[1]:.4f}]")
design_effect = (np.nanstd(boot_E) * np.sqrt(n_lin_agg)) / (np.sqrt(P_top_high*(1-P_top_high)/high_m.sum())) \
                if high_m.sum() > 0 and not np.isnan(P_top_high) else np.nan
print(f"[MEASURED] Design effect (SE_cluster / SE_naive): {design_effect:.2f}")

# ── Placebo: permute confidence_tier across pairs ─────────────────────────────
print("\n[MEASURED] PLACEBO: permuting confidence_tier (1000 permutations) ...")
if "confidence_tier" in pred.columns:
    tiers = pred["confidence_tier"].values.copy()
    rng2  = np.random.default_rng(99)
    plac_E = []
    for _ in range(1000):
        perm = rng2.permuted(tiers)
        high_p = perm == "HIGH"
        pt_h = pred.loc[high_p, "top20"].mean() if high_p.sum() > 0 else np.nan
        plac_E.append(pt_h/P_top if not np.isnan(pt_h) else np.nan)
    pla_E = np.nanmean(plac_E)
    print(f"[MEASURED] Placebo E_high mean: {pla_E:.4f}  p95: {np.nanpercentile(plac_E,95):.4f}")
    print(f"[MEASURED] Placebo/observed: {pla_E/E_high:.3f}" if E_high > 0 else "")
else:
    pla_E = np.nan
    print("[MEASURED] No confidence_tier column — skipping placebo")

# ── By tier (if available) ────────────────────────────────────────────────────
if "confidence_tier" in pred.columns:
    print(f"\n[MEASURED] Enrichment by tier:")
    for tier in ["HIGH","MEDIUM","LOW"]:
        m = pred["confidence_tier"] == tier
        if m.sum() == 0: continue
        pt = pred.loc[m,"top20"].mean()
        E  = pt / P_top
        print(f"  {tier:8s}: n={m.sum():>8,}  E_top20={E:.4f}×")

print("\n" + "=" * 72)
print("VERDICT AGAINST PRE-DECLARATION:")
print(f"  Pre-declared: E_high > 1.146 (prior pooled enrichment), CI excludes 1.0")
print(f"  E_high={E_high:.4f}  CI=[{ci[0]:.4f},{ci[1]:.4f}]")
if not np.isnan(E_high) and E_high > 1.146 and ci[0] > 1.0:
    print("  CONFIRMED — HIGH tier enriched above pooled baseline, CI excludes 1.0")
elif not np.isnan(E_high) and ci[0] > 1.0:
    print(f"  PARTIAL — E_high ({E_high:.4f}) < 1.146 but CI excludes 1.0")
elif not np.isnan(E_high) and E_high > 1.0:
    print(f"  WEAK — E_high > 1.0 but CI includes 1.0 or E_high < 1.146")
else:
    print("  NOT CONFIRMED")
if not np.isnan(pla_E) and E_high > 0 and pla_E/E_high >= 0.5:
    print(f"  WARNING: PLACEBO CONTAMINATED (placebo {pla_E:.4f} / observed {E_high:.4f} = {pla_E/E_high:.2f})")
print("=" * 72)
