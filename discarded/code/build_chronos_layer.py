"""
Improvement 2 — Chronos as a third scoring layer.

Steps:
  1. Load Chronos, bridge SIDM -> ACH model IDs
  2. Percentile-normalise per gene; invert so high = essential
  3. 50/50 cell-line train/test split (anti-circularity)
  4. Measure pairwise correlations (expr-chronos, prot-chronos)
  5. Build 3-layer core_score on TRAIN lines, validate on TEST lines
  6. Compare 2-layer vs 3-layer hit@20 for BRAF and BCL2

Run:
    python src/pipeline/build_chronos_layer.py
"""
import numpy as np
import pandas as pd
from pathlib import Path
from scipy.stats import spearmanr
import sys
sys.path.insert(0, "src/pipeline")
from scoring_variants import combine_three

OUTPUTS   = Path("src/pipeline/outputs")
VAL       = Path("validation/prepared")

# ── 1. Load Chronos and bridge to ACH IDs ────────────────────────────────────
print("=== Step 1: Load Chronos ===")
ch = pd.read_parquet(VAL / "chronos_long.parquet")
ib = pd.read_parquet(VAL / "id_bridge.parquet")
ch = ch.merge(ib, on="sanger_model_id", how="inner")
ch = ch.rename(columns={"model_id": "model_id_ach", "essentiality": "chronos_raw"})
ch["model_id"] = ch["model_id_ach"].str.upper()
ch["ensg_id"]  = ch["ensg_id"].str.upper()

print(f"Chronos raw: {ch.shape[0]:,} rows")
print(f"  value range: min={ch.chronos_raw.min():.4f}  max={ch.chronos_raw.max():.4f}  mean={ch.chronos_raw.mean():.4f}")
print(f"  unique cell lines (ACH): {ch.model_id.nunique()}")
print(f"  unique genes: {ch.ensg_id.nunique()}")

# ── 2. Percentile-normalise per gene, invert so high = essential ──────────────
print("\n=== Step 2: Normalise ===")
# More negative chronos = more essential → rank descending then percentile
ch["chronos_pct"] = (
    ch.groupby("ensg_id")["chronos_raw"]
      .rank(ascending=True, pct=True, method="average")   # low raw → low rank
)
# Invert: essential (most negative) gets chronos_pct → 1.0
ch["chronos_pct"] = 1.0 - ch["chronos_pct"]
print(f"chronos_pct range: {ch.chronos_pct.min():.4f} – {ch.chronos_pct.max():.4f}")

# ── 3. Train/test split on cell lines ────────────────────────────────────────
print("\n=== Step 3: Train/test split ===")
all_lines = np.array(sorted(ch.model_id.unique()))
rng = np.random.default_rng(42)
rng.shuffle(all_lines)
cut = len(all_lines) // 2
train_lines = set(all_lines[:cut])
test_lines  = set(all_lines[cut:])
print(f"Train lines: {len(train_lines)}, Test lines: {len(test_lines)}")

ch_train = ch[ch.model_id.isin(train_lines)][["model_id","ensg_id","chronos_pct"]]
ch_test  = ch[ch.model_id.isin(test_lines)][["model_id","ensg_id","chronos_pct"]]

# ── 4. Measure pairwise correlations ─────────────────────────────────────────
print("\n=== Step 4: Pairwise correlations ===")
core = pd.read_parquet(OUTPUTS / "core_score.parquet")
core["model_id"] = core["model_id"].str.upper()
core["ensg_id"]  = core["ensg_id"].str.upper()

# Load expr and prot percentiles from harmonised parquet
# We only have core_score (combination), so proxy with n_layers==1 as expr-only rows
# and reconstruct from layer3_continuous if available
try:
    l3 = pd.read_parquet(OUTPUTS / "layer3_continuous.parquet")
    l3["model_id"] = l3["model_id"].str.upper()
    l3["ensg_id"]  = l3["ensg_id"].str.upper()
    has_l3 = True
    print("  layer3_continuous available, using for expr/prot signals")
    print(f"  layer3_continuous cols: {l3.columns.tolist()}")
except Exception:
    has_l3 = False
    print("  layer3_continuous not available — using core_score as expr proxy for correlation only")

if has_l3:
    # Use expr_pct and prot_pct from layer3
    expr_col = [c for c in l3.columns if "expr" in c.lower() or "mrna" in c.lower() or "tpm" in c.lower()]
    prot_col = [c for c in l3.columns if "prot" in c.lower() or "ms" in c.lower()]
    print(f"  expr cols: {expr_col}, prot cols: {prot_col}")

# Merge chronos_train with core_score (2-layer) to get shared rows for correlation
merged = core.merge(ch_train, on=["model_id","ensg_id"], how="inner")
merged_2layer = merged[merged.n_layers == 2].copy()  # rows with both expr and prot
print(f"  Shared (expr+prot, chronos) rows on train: {len(merged_2layer):,}")

if len(merged_2layer) > 1000:
    # core_score is combination of expr+prot; correlate each proxy
    # For expr-chronos and prot-chronos: use n_layers==1 rows as proxies
    expr_only = core[core.n_layers == 1].merge(ch_train, on=["model_id","ensg_id"], how="inner")

    rho_core_c, _ = spearmanr(merged_2layer["core_score"], merged_2layer["chronos_pct"])
    rho_expr_c, _ = spearmanr(expr_only["core_score"],    expr_only["chronos_pct"])
    print(f"\n  core_score(2-layer) vs chronos rho:  {rho_core_c:.4f}")
    print(f"  core_score(1-layer/expr) vs chronos rho: {rho_expr_c:.4f}")
    print(f"  (True expr-chronos and prot-chronos require un-combined signals)")
    print(f"  Interpretation: {'WEAK' if abs(rho_core_c) < 0.3 else 'MODERATE' if abs(rho_core_c) < 0.6 else 'STRONG'} "
          f"correlation — Chronos adds {'mostly independent' if abs(rho_core_c) < 0.3 else 'partially redundant'} information")
else:
    print("  Not enough shared rows for correlation")

# ── 5. Build 3-layer scores on TRAIN lines ───────────────────────────────────
print("\n=== Step 5: Build 3-layer scores (train lines) ===")

# Weights from combine_three: inverse-correlation penalised weighted sum
# rho_ep=0.46, rho_ec=0.005, rho_pc=0.005
rho_ep, rho_ec, rho_pc = 0.46, 0.005, 0.005
mean_rho_e = (abs(rho_ep) + abs(rho_ec)) / 2
mean_rho_p = (abs(rho_ep) + abs(rho_pc)) / 2
mean_rho_c = (abs(rho_ec) + abs(rho_pc)) / 2
w_e_raw = 1 / (1 + mean_rho_e)
w_p_raw = 1 / (1 + mean_rho_p)
w_c_raw = 1 / (1 + mean_rho_c)
total   = w_e_raw + w_p_raw + w_c_raw
w_e, w_p, w_c = w_e_raw / total, w_p_raw / total, w_c_raw / total
# core_score = combined E+P signal; merge as w_ep*core + w_c*chronos, re-normalised
w_ep     = (w_e + w_p) / (w_e + w_p + w_c)
w_c_norm = w_c / (w_e + w_p + w_c)
print(f"  Weighted sum: core_score(E+P) weight={w_ep:.3f}, chronos weight={w_c_norm:.3f}")

core_train = core[core.model_id.isin(train_lines)].copy()
core_test  = core[core.model_id.isin(test_lines)].copy()

c3 = core_train.merge(ch_train, on=["model_id","ensg_id"], how="left")
has_c = c3["chronos_pct"].notna()
c3["core_score_3"] = c3["core_score"].copy()
c3.loc[has_c, "core_score_3"] = (
    w_ep * c3.loc[has_c, "core_score"] + w_c_norm * c3.loc[has_c, "chronos_pct"]
)
c3.loc[has_c, "n_layers"] = c3.loc[has_c, "n_layers"] + 1

# Recompute stratum_rank with new n_layers
c3["stratum_rank_3"] = (
    c3.groupby(["ensg_id","n_layers"])["core_score_3"]
      .rank(pct=True, method="average")
)

n_three = has_c.sum()
print(f"  Pairs upgraded to 3-layer: {n_three:,} ({100*n_three/len(c3):.1f}%)")
print(f"  3-layer n_layers dist: {c3.n_layers.value_counts().sort_index().to_dict()}")

# ── 6. Validate on TEST lines: 2-layer vs 3-layer hit@20 ─────────────────────
print("\n=== Step 6: Validation on TEST lines ===")

# Load GDSC validation pairs
gdsc = pd.read_parquet(VAL / "gdsc_scored_ready.parquet")
print(f"  GDSC scored_ready cols: {gdsc.columns.tolist()}")
print(f"  GDSC shape: {gdsc.shape}")

# Use the curated pairs for BRAF and BCL2
try:
    curated = pd.read_parquet(VAL / "curated_validation_pairs.parquet")
    print(f"  Curated pairs shape: {curated.shape}")
    print(f"  Curated cols: {curated.columns.tolist()}")
    print(curated.head(5).to_string())
except Exception as e:
    print(f"  curated_validation_pairs error: {e}")

# Check what columns gdsc_scored_ready has for hit@20 reconstruction
print()
print(gdsc.head(3).to_string())
