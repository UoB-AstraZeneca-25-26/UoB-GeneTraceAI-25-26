"""
T5 — D22 v2: Stratified directional enrichment (fixed CNA columns, fast placebo).
GAIN: oncogene/both AND (has_cna_alteration OR non-truncating mut_driver OR fusion_driver)
LOSS: tsg/both AND (has_cna_alteration OR truncating mut_driver)
Note: has_cna_alteration already encodes direction — oncogene rows are amplification,
      tsg rows are deletion — because cna_layer.py gates on gene_role.
Truncating proxy: max_vep_rank == 3 (frameshift/stop_gained/start_lost).
Run from repo root: python diagnostics/T5_D22_v2.py
"""
import sys
from pathlib import Path
import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "final_pipeline"))
from config import PREDICTIONS, CNA_FLAGS, MUTATIONS_SCR, FUSIONS_SCR, DB
import duckdb

print("=" * 72)
print("T5 — D22: STRATIFIED DIRECTIONAL ENRICHMENT (v2)")
print("=" * 72)

# ── Load predictions + lineage ─────────────────────────────────────────────────
pred = pd.read_parquet(PREDICTIONS)
print(f"[COUNT] predictions: {len(pred):,}")

con = duckdb.connect(str(DB), read_only=True)
lin = con.execute("SELECT model_id, lineage FROM main.sample_info WHERE lineage IS NOT NULL").df()
lin = lin.drop_duplicates("model_id")
con.close()
pred = pred.merge(lin, on="model_id", how="left")
print(f"[COUNT] Rows missing lineage: {pred.lineage.isna().sum():,}")
print(f"[COUNT] Distinct lineages: {pred.lineage.nunique():,}")
assert pred.lineage.nunique() >= 20

# ── Load mutation truncating flag ──────────────────────────────────────────────
mut_raw = pd.read_parquet(REPO / "cleaned_track_data/mutations_collapsed.parquet",
                           columns=["ensg_id","model_id","max_vep_rank"])
mut_raw["model_id"] = mut_raw["model_id"].str.lower()
# max_vep_rank==3 → frameshift/stop_gained/start_lost (truncating proxy)
mut_raw["is_truncating"] = mut_raw["max_vep_rank"] == 3
trunc_flags = mut_raw.groupby(["ensg_id","model_id"])["is_truncating"].any().reset_index()
print(f"[COUNT] Truncating flags: {trunc_flags.is_truncating.sum():,} pairs")

pred = pred.merge(trunc_flags, on=["ensg_id","model_id"], how="left")
pred["is_truncating"] = pred["is_truncating"].fillna(False)

# ── Load CNA + fusion flags ────────────────────────────────────────────────────
cna = pd.read_parquet(CNA_FLAGS, columns=["model_id","ensg_id","has_cna_alteration"])
fus = pd.read_parquet(FUSIONS_SCR, columns=["model_id","ensg_id","p_fusion"])
pred = pred.merge(cna, on=["model_id","ensg_id"], how="left")
pred = pred.merge(fus, on=["model_id","ensg_id"], how="left")
pred["cna_driver"]    = pred["has_cna_alteration"].fillna(False)
pred["fusion_driver"] = pred["p_fusion"].fillna(0) >= 0.50
pred["mut_driver"]    = pred["p_mutation"].fillna(0) >= 0.50

# ── Role flags ─────────────────────────────────────────────────────────────────
pred["is_onc"] = pred["gene_role"].isin(["oncogene","both"])
pred["is_tsg"] = pred["gene_role"].isin(["tsg","both"])

# ── GAIN / LOSS partition ──────────────────────────────────────────────────────
# has_cna_alteration is ALREADY direction-gated by cna_layer.py:
#   oncogene→amplification, tsg→deletion. So cna_driver AND is_onc = amplification of oncogene.
gain_mask = (
    (pred["cna_driver"] & pred["is_onc"]) |
    (pred["mut_driver"] & pred["is_onc"] & ~pred["is_truncating"]) |
    (pred["fusion_driver"] & pred["is_onc"])
)
loss_mask = (
    (pred["cna_driver"] & pred["is_tsg"]) |
    (pred["mut_driver"] & pred["is_tsg"] & pred["is_truncating"])
)
driver_mask = pred["has_driver_alteration"]
other_mask  = driver_mask & ~gain_mask & ~loss_mask

print(f"\n[MEASURED] Driver partition:")
print(f"  GAIN  (oncogene activation): {gain_mask.sum():,}")
print(f"  LOSS  (TSG inactivation):    {loss_mask.sum():,}")
print(f"  OTHER (unclassified):        {other_mask.sum():,}")
print(f"  Total (GAIN+LOSS+OTHER):     {gain_mask.sum()+loss_mask.sum()+other_mask.sum():,}  vs drivers={driver_mask.sum():,}")

# ── Point estimates ────────────────────────────────────────────────────────────
pred["top20"] = pred["score_rank_pct"] >= 0.80
pred["bot20"] = pred["score_rank_pct"] <  0.20
P_top = pred["top20"].mean()
P_bot = pred["bot20"].mean()
P_top_gain = pred.loc[gain_mask,"top20"].mean() if gain_mask.sum() > 0 else np.nan
P_bot_loss = pred.loc[loss_mask,"bot20"].mean() if loss_mask.sum() > 0 else np.nan
E_gain = P_top_gain / P_top
E_loss = P_bot_loss / P_bot if not np.isnan(P_bot_loss) else np.nan
pooled = pred.loc[driver_mask,"top20"].mean() / P_top

print(f"\n[MEASURED] P(top20%) overall = {P_top:.6f}")
print(f"[MEASURED] P(bot20%) overall = {P_bot:.6f}")
print(f"[MEASURED] P(top20% | GAIN)  = {P_top_gain:.6f}  (n={gain_mask.sum():,})")
print(f"[MEASURED] E_gain = {E_gain:.4f}")
print(f"[MEASURED] P(bot20% | LOSS)  = {P_bot_loss:.6f}  (n={loss_mask.sum():,})")
print(f"[MEASURED] E_loss = {E_loss:.4f}")
print(f"[MEASURED] Pooled enrichment (recomputed) = {pooled:.4f}×")

# ── Fast cluster bootstrap: per-lineage aggregates ────────────────────────────
print("\n[MEASURED] Cluster bootstrap (1000 resamples, cluster=lineage) ...")

# Pre-aggregate per lineage
agg = pred.groupby("lineage").apply(lambda g: pd.Series({
    "top_gain": int(g.loc[gain_mask.reindex(g.index,fill_value=False),"top20"].sum()),
    "n_gain":   int(gain_mask.reindex(g.index,fill_value=False).sum()),
    "bot_loss": int(g.loc[loss_mask.reindex(g.index,fill_value=False),"bot20"].sum()),
    "n_loss":   int(loss_mask.reindex(g.index,fill_value=False).sum()),
    "n_top":    int(g["top20"].sum()),
    "n_bot":    int(g["bot20"].sum()),
    "n_total":  len(g),
}), include_groups=False).reset_index()

lineages = agg["lineage"].values
n_lin    = len(lineages)
print(f"  Lineages: {n_lin}")
rng = np.random.default_rng(42)
boot_Eg, boot_El, boot_pool = [], [], []
for _ in range(1000):
    idx = rng.integers(0, n_lin, size=n_lin)
    s   = agg.iloc[idx].sum(numeric_only=True)
    p_t_g = s["top_gain"] / s["n_gain"] if s["n_gain"] > 0 else np.nan
    p_b_l = s["bot_loss"] / s["n_loss"] if s["n_loss"] > 0 else np.nan
    p_t   = s["n_top"]    / s["n_total"] if s["n_total"] > 0 else np.nan
    p_b   = s["n_bot"]    / s["n_total"] if s["n_total"] > 0 else np.nan
    boot_Eg.append(p_t_g/p_t if p_t else np.nan)
    boot_El.append(p_b_l/p_b if p_b else np.nan)

boot_Eg = np.array(boot_Eg); boot_El = np.array(boot_El)
ci_Eg = (np.nanpercentile(boot_Eg,2.5), np.nanpercentile(boot_Eg,97.5))
ci_El = (np.nanpercentile(boot_El,2.5), np.nanpercentile(boot_El,97.5))
print(f"[MEASURED] E_gain 95% CI: [{ci_Eg[0]:.4f}, {ci_Eg[1]:.4f}]")
print(f"[MEASURED] E_loss 95% CI: [{ci_El[0]:.4f}, {ci_El[1]:.4f}]")

# ── MANDATORY PLACEBO: permute gene_role ACROSS GENES ─────────────────────────
print("\n[MEASURED] PLACEBO: permuting gene_role across genes ...")
gene_roles = pred[["ensg_id","gene_role"]].drop_duplicates("ensg_id")
n_genes = len(gene_roles)
rng2 = np.random.default_rng(99)
plac_Eg, plac_El = [], []

for _ in range(1000):
    shuffled_roles = gene_roles["gene_role"].values.copy()
    rng2.shuffle(shuffled_roles)
    role_map = dict(zip(gene_roles["ensg_id"].values, shuffled_roles))
    role_perm = pred["ensg_id"].map(role_map)
    is_onc_p = role_perm.isin(["oncogene","both"])
    is_tsg_p = role_perm.isin(["tsg","both"])
    gain_p = (
        (pred["cna_driver"] & is_onc_p) |
        (pred["mut_driver"] & is_onc_p & ~pred["is_truncating"]) |
        (pred["fusion_driver"] & is_onc_p)
    )
    loss_p = (
        (pred["cna_driver"] & is_tsg_p) |
        (pred["mut_driver"] & is_tsg_p & pred["is_truncating"])
    )
    pt_g = pred.loc[gain_p,"top20"].mean() if gain_p.sum()>0 else np.nan
    pb_l = pred.loc[loss_p,"bot20"].mean() if loss_p.sum()>0 else np.nan
    plac_Eg.append(pt_g/P_top if not np.isnan(pt_g) else np.nan)
    plac_El.append(pb_l/P_bot if not np.isnan(pb_l) else np.nan)

pla_Eg = np.nanmean(plac_Eg); pla_El = np.nanmean(plac_El)
print(f"[MEASURED] Placebo E_gain mean:  {pla_Eg:.4f}  (observed {E_gain:.4f})")
print(f"[MEASURED] Placebo E_loss mean:  {pla_El:.4f}  (observed {E_loss:.4f})")
print(f"[MEASURED] Placebo p95 E_gain:   {np.nanpercentile(plac_Eg,95):.4f}")
print(f"[MEASURED] Placebo fraction/observed (gain): {pla_Eg/E_gain:.3f}" if E_gain>0 else "E_gain=0")
print(f"[MEASURED] Placebo fraction/observed (loss): {pla_El/E_loss:.3f}" if not np.isnan(E_loss) and E_loss>0 else "")

# ── Verdict ────────────────────────────────────────────────────────────────────
print("\n" + "=" * 72)
gain_conf = not np.isnan(E_gain) and E_gain > 1.3 and ci_Eg[0] > 1.0
loss_conf = not np.isnan(E_loss) and E_loss > 1.3 and ci_El[0] > 1.0
gain_plac = not np.isnan(pla_Eg) and E_gain > 0 and (pla_Eg/E_gain) >= 0.50
print(f"VERDICT AGAINST PRE-DECLARATION:")
print(f"  Pre-declared: E_gain > 1.3 AND E_loss > 1.3, CIs exclude 1.0")
print(f"  E_gain = {E_gain:.4f}  [{ci_Eg[0]:.4f},{ci_Eg[1]:.4f}]  placebo={pla_Eg:.4f}")
print(f"  E_loss = {E_loss:.4f}  [{ci_El[0]:.4f},{ci_El[1]:.4f}]  placebo={pla_El:.4f}")
if gain_conf and loss_conf:
    print("  CONFIRMED")
elif gain_conf:
    print("  PARTIAL — E_gain confirmed; E_loss NOT confirmed")
elif loss_conf:
    print("  PARTIAL — E_loss confirmed; E_gain NOT confirmed")
else:
    print("  NOT CONFIRMED")
if gain_plac:
    print("  WARNING: PLACEBO CONTAMINATED (gain placebo >=50% of observed)")
print("=" * 72)
