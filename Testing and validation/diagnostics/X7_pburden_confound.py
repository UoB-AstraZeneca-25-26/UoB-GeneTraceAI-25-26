"""
X7 — p_burden confound quantification.
COWORK v4.0. No code changes. Read-only.
Run from repo root: python diagnostics/X7_pburden_confound.py
"""
import sys
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import stats

REPO = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(REPO / "architecture"))

MUT_SRC  = REPO / "cleaned_track_data" / "mutations_collapsed.parquet"
MUT_SCR  = REPO / "architecture" / "outputs" / "mutations_scores.parquet"
GENE_LKP = REPO / "reference" / "gene_lookup.parquet"

P0_VEP    = 2.5; K_VEP    = 2.0
P0_PATH   = 0.5; K_PATH   = 2.0
P0_BURDEN = 3.0; K_BURDEN = 1.5
DRIVER_BOOST   = 0.15
MUT_GATE_THRESH = 0.5

print("=" * 72)
print("X7 -- p_burden CONFOUND QUANTIFICATION (COWORK v4.0)")
print("=" * 72)

def hill(x, p0, k):
    xk = np.power(np.clip(x, 0, None), k)
    return xk / (xk + p0 ** k)

mut = pd.read_parquet(MUT_SRC)
print(f"[COUNT] mutations_collapsed rows: {len(mut):,}  "
      f"genes: {mut.ensg_id.nunique():,}  models: {mut.model_id.nunique():,}")
print(f"  columns: {list(mut.columns)}")

mut = mut.copy()
mut["variant_burden"] = mut["variant_count"].clip(upper=5)
mut["p_vep"]    = hill(mut["max_vep_rank"].astype(float), P0_VEP, K_VEP)
mut["p_path"]   = hill(mut["max_pathogenicity"].astype(float), P0_PATH, K_PATH)
mut["p_burden"] = hill(mut["variant_burden"].astype(float), P0_BURDEN, K_BURDEN)
p_base = 1.0 - (1.0-mut["p_vep"].fillna(0)) * (1.0-mut["p_path"].fillna(0)) * (1.0-mut["p_burden"].fillna(0))
driver_flag = mut["any_driver"].fillna(False).astype(int)
mut["p_mutation"] = (p_base + DRIVER_BOOST * driver_flag).clip(upper=1.0)
mut["mut_driver"] = mut["p_mutation"] >= MUT_GATE_THRESH

# Per-line total variant count (TMB proxy)
tmb = mut.groupby("model_id")["variant_count"].sum().reset_index(name="tmb")
mut2 = mut.merge(tmb, on="model_id", how="left")
print(f"\n[COUNT] TMB distribution:")
print(f"  mean={tmb.tmb.mean():.1f}  median={tmb.tmb.median():.1f}  "
      f"p95={tmb.tmb.quantile(0.95):.1f}  max={tmb.tmb.max():.1f}")

# Hypermutators: top 5% by TMB
tmb_thresh = tmb.tmb.quantile(0.95)
hyper_models = set(tmb[tmb.tmb >= tmb_thresh].model_id)
print(f"\n[COUNT] Hypermutators (top 5%, TMB >= {tmb_thresh:.0f}): "
      f"{len(hyper_models):,} / {tmb.model_id.nunique():,} lines")

# corr(variant_count, p_mutation)
r1, p1 = stats.pearsonr(mut2["variant_count"].clip(upper=100), mut2["p_mutation"])
print(f"\n[MEASURED] Pearson r(variant_count, p_mutation): {r1:.4f}  p={p1:.2e}")
r1s, p1s = stats.spearmanr(mut2["variant_count"], mut2["p_mutation"])
print(f"[MEASURED] Spearman r(variant_count, p_mutation): {r1s:.4f}  p={p1s:.2e}")

# corr(TMB, p_mutation) per pair
r2, p2 = stats.pearsonr(mut2["tmb"], mut2["p_mutation"])
r2s, p2s = stats.spearmanr(mut2["tmb"], mut2["p_mutation"])
print(f"\n[MEASURED] Pearson  r(per_line_tmb, p_mutation): {r2:.4f}  p={p2:.2e}")
print(f"[MEASURED] Spearman r(per_line_tmb, p_mutation): {r2s:.4f}  p={p2s:.2e}")

# mut_driver rate: hypermutators vs rest
hyper = mut2[mut2.model_id.isin(hyper_models)]
rest  = mut2[~mut2.model_id.isin(hyper_models)]
hyper_rate = hyper["mut_driver"].mean()
rest_rate  = rest["mut_driver"].mean()
ratio = hyper_rate / max(rest_rate, 1e-6)
print(f"\n[MEASURED] mut_driver rate in hypermutators:  {hyper_rate*100:.2f}%")
print(f"[MEASURED] mut_driver rate in rest:            {rest_rate*100:.2f}%")
print(f"[MEASURED] ratio: {ratio:.2f}×  (pre-declared: >= 1.5×)")

# Lineage distribution of hypermutators
if "model_id" in mut.columns:
    try:
        DB = REPO / "architecture" / "outputs" / "celllineselector.db"
        import duckdb
        con = duckdb.connect(str(DB), read_only=True)
        si = con.execute("SELECT lower(model_id) AS model_id, lineage FROM main.sample_info").df()
        con.close()
        hyper_si = si[si.model_id.isin({m.lower() for m in hyper_models})]
        print(f"\n[MEASURED] Lineage distribution of hypermutators (n={len(hyper_si)}):")
        print(hyper_si.lineage.value_counts().head(10).to_string())
    except Exception as e:
        print(f"[WARN] Could not get lineage: {e}")

# Ablation: mut_driver without p_burden
mut2["p_mutation_noburd"] = (
    1.0 - (1.0-mut2["p_vep"].fillna(0)) * (1.0-mut2["p_path"].fillna(0))
    + DRIVER_BOOST * mut2["any_driver"].fillna(False).astype(int)
).clip(upper=1.0)
mut2["mut_driver_noburd"] = mut2["p_mutation_noburd"] >= MUT_GATE_THRESH
n_changed = (mut2["mut_driver"] != mut2["mut_driver_noburd"]).sum()
n_lost     = (mut2["mut_driver"] & ~mut2["mut_driver_noburd"]).sum()
n_gained   = (~mut2["mut_driver"] & mut2["mut_driver_noburd"]).sum()
print(f"\n--- ABLATION: remove p_burden entirely ---")
print(f"  Changed mut_driver calls: {n_changed:,} / {len(mut2):,} "
      f"({100*n_changed/len(mut2):.2f}%)")
print(f"  Lost (was driver, now not): {n_lost:,}")
print(f"  Gained (was not, now is):   {n_gained:,}")

print(f"\n--- PRE-DECLARED CHECKS ---")
print(f"  pre-declared: r(p_mutation, per_line_tmb) > 0.2")
print(f"  measured Spearman: {r2s:.4f}  {'CONFIRMED' if r2s > 0.2 else 'WRONG'}")
print(f"  pre-declared: hypermutator mut_driver rate >= 1.5x rest")
print(f"  measured ratio: {ratio:.2f}x  {'CONFIRMED' if ratio >= 1.5 else 'WRONG'}")

print("\n" + "=" * 72)
