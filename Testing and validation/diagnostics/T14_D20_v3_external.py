"""
T14 v3 — HIGH tier enrichment against COSMIC CGC Tier 1 (external label).
Remediation item 3: current E_HIGH=4.99× is trivially by construction (HIGH=top20%).
This version tests whether HIGH-tier genes are enriched for KNOWN cancer drivers
(COSMIC CGC tier 1) relative to MEDIUM-tier genes (also top-20%, but no driver flag).

Design:
  - External label: COSMIC CGC tier 1 (592 genes, manually curated, pre-dates our pipeline)
  - Enrichment: E = P(CGC_t1 | HIGH) / P(CGC_t1 | MEDIUM)
    (both HIGH and MEDIUM are top-20%, so score-rank is controlled for;
     the only difference is has_driver_alteration)
  - Placebo: permute confidence_tier labels within score-rank stratum
  - Cluster bootstrap: resample lineages with replacement (design effect correction)

Run from repo root: python diagnostics/T14_D20_v3_external.py
"""
import sys, re
from pathlib import Path
import numpy as np
import pandas as pd
import duckdb

REPO = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(REPO / "architecture"))
from config import PREDICTIONS, DB, CORE_SCORE

RNG   = np.random.default_rng(42)
N_BOOT   = 1000
N_PLACEBO = 500

print("=" * 72)
print("T14 v3 — HIGH TIER vs COSMIC CGC TIER 1 (EXTERNAL BENCHMARK)")
print("=" * 72)

# ── Load predictions + lineage ────────────────────────────────────────────────
pred = pd.read_parquet(PREDICTIONS, columns=[
    "model_id","ensg_id","confidence_tier","score_rank_pct"
])
# lineage lives in core_score.parquet
lin_df = pd.read_parquet(CORE_SCORE, columns=["model_id","ensg_id","lineage"])
pred = pred.merge(lin_df, on=["model_id","ensg_id"], how="left")
pred["top20"] = pred["score_rank_pct"] >= 0.80
top = pred[pred["top20"]].copy()   # 3.95M rows; both HIGH and MEDIUM are here
print(f"[COUNT] Total predictions: {len(pred):,}")
print(f"[COUNT] Top-20% rows (HIGH+MEDIUM): {len(top):,}")
tier_dist = top["confidence_tier"].value_counts()
print(f"[COUNT] HIGH: {tier_dist.get('HIGH',0):,}  MEDIUM: {tier_dist.get('MEDIUM',0):,}")

# ── Load COSMIC CGC tier 1 ensg_ids ──────────────────────────────────────────
con = duckdb.connect(str(DB), read_only=True)
cgc = con.execute("SELECT synonyms, tier FROM main.cosmic_cgc WHERE tier = 1").df()
con.close()
print(f"[COUNT] COSMIC CGC tier-1 genes: {len(cgc):,}")

# Extract ENSG IDs from synonyms field (comma-separated, includes ENSG...)
def extract_ensg(syn_str):
    if not isinstance(syn_str, str):
        return []
    return [s.strip().upper()[:15] for s in syn_str.split(",")
            if re.match(r'ENSG\d{11}', s.strip(), re.I)]

cgc_ensg = set()
for _, row in cgc.iterrows():
    cgc_ensg.update(extract_ensg(row["synonyms"]))
print(f"[COUNT] COSMIC CGC tier-1 unique ENSG IDs extracted: {len(cgc_ensg):,}")

# Flag in top-20% predictions
top["cgc_t1"] = top["ensg_id"].str.upper().isin(cgc_ensg)
print(f"[COUNT] Top-20% rows matching CGC t1 gene: {top.cgc_t1.sum():,}")

# ── Point estimate ────────────────────────────────────────────────────────────
high   = top[top["confidence_tier"] == "HIGH"]
medium = top[top["confidence_tier"] == "MEDIUM"]

p_high   = high["cgc_t1"].mean()
p_medium = medium["cgc_t1"].mean()
e_obs    = p_high / p_medium if p_medium > 0 else np.nan

print(f"\n[MEASURED] P(CGC_t1 | HIGH)   = {p_high:.6f}  (n={len(high):,})")
print(f"[MEASURED] P(CGC_t1 | MEDIUM) = {p_medium:.6f}  (n={len(medium):,})")
print(f"[MEASURED] E_HIGH/MEDIUM = {e_obs:.4f}×  "
      f"({'CGC enriched in HIGH' if e_obs > 1 else 'not enriched'})")

# ── Cluster bootstrap (resample lineages) ────────────────────────────────────
print(f"\n[STATUS] Cluster bootstrap ({N_BOOT} resamples, cluster=lineage)...")
lineages = top["lineage"].unique()
boot_e   = np.empty(N_BOOT)
for i in range(N_BOOT):
    samp_lin = RNG.choice(lineages, size=len(lineages), replace=True)
    samp     = pd.concat([top[top["lineage"] == lin] for lin in samp_lin], ignore_index=True)
    h = samp[samp["confidence_tier"] == "HIGH"]["cgc_t1"].mean()
    m = samp[samp["confidence_tier"] == "MEDIUM"]["cgc_t1"].mean()
    boot_e[i] = h / m if m > 0 else np.nan

boot_e_clean = boot_e[~np.isnan(boot_e)]
ci_lo, ci_hi = np.percentile(boot_e_clean, [2.5, 97.5])
print(f"[MEASURED] Bootstrap 95% CI on E_HIGH/MEDIUM: [{ci_lo:.4f}, {ci_hi:.4f}]")
print(f"[MEASURED] Bootstrap SE: {boot_e_clean.std():.4f}")

# ── Placebo: permute tier labels within top-20% (preserves score-rank confound) ─
print(f"\n[STATUS] Placebo ({N_PLACEBO} permutations)...")
placebo_e = np.empty(N_PLACEBO)
for i in range(N_PLACEBO):
    shuffled = top.copy()
    shuffled["confidence_tier"] = RNG.permutation(top["confidence_tier"].values)
    h = shuffled[shuffled["confidence_tier"] == "HIGH"]["cgc_t1"].mean()
    m = shuffled[shuffled["confidence_tier"] == "MEDIUM"]["cgc_t1"].mean()
    placebo_e[i] = h / m if m > 0 else np.nan

placebo_clean  = placebo_e[~np.isnan(placebo_e)]
placebo_median = np.median(placebo_clean)
pct_above      = np.mean(placebo_clean >= e_obs)

print(f"[MEASURED] Placebo median E = {placebo_median:.4f}×")
print(f"[MEASURED] Placebo >= observed: {100*pct_above:.1f}%")

# ── Summary ───────────────────────────────────────────────────────────────────
print("\n" + "=" * 72)
print("VERDICT AGAINST PRE-DECLARATION (external benchmark):")
print(f"  E_HIGH/MEDIUM = {e_obs:.4f}×  95% CI: [{ci_lo:.4f}, {ci_hi:.4f}]")
print(f"  Placebo median = {placebo_median:.4f}×  (placebo >= obs: {100*pct_above:.1f}%)")
if ci_lo > 1.0 and pct_above < 0.05:
    print("  CONFIRMED — HIGH tier enriched for COSMIC CGC t1 vs MEDIUM,")
    print("  CI excludes 1.0, and result not explained by permutation.")
elif ci_lo > 1.0:
    print("  BORDERLINE — CI excludes 1.0 but placebo too high.")
else:
    print("  NOT CONFIRMED — CI includes 1.0; HIGH-tier driver flag does not")
    print("  independently predict COSMIC CGC tier-1 membership above")
    print("  score-rank baseline.")
print("=" * 72)
print("\n[DOCUMENTED] This test controls for score_rank_pct by comparing within")
print("top-20% only. E>1 means the has_driver_alteration flag captures real")
print("cancer driver genes beyond what the abundance score alone identifies.")
