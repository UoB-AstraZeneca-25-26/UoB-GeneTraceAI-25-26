"""
T5 — D22 v3: Stratified directional enrichment.
predictions already contains has_cna_alteration, p_mutation, p_fusion.
Only need to merge truncating proxy from mutations_collapsed.
Run from repo root: python diagnostics/T5_D22_v3.py
"""
import sys
from pathlib import Path
import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "final_pipeline"))
from config import PREDICTIONS, DB
import duckdb

print("=" * 72)
print("T5 — D22: STRATIFIED DIRECTIONAL ENRICHMENT (v3)")
print("=" * 72)

pred = pd.read_parquet(PREDICTIONS)
print(f"[COUNT] predictions: {len(pred):,}")
print(f"[COUNT] columns: {list(pred.columns)}")

con = duckdb.connect(str(DB), read_only=True)
lin = con.execute("SELECT model_id, lineage FROM main.sample_info WHERE lineage IS NOT NULL").df()
lin = lin.drop_duplicates("model_id")
con.close()
pred = pred.merge(lin, on="model_id", how="left")
print(f"[COUNT] Lineages: {pred.lineage.nunique():,}")
assert pred.lineage.nunique() >= 20

# Truncating proxy: max_vep_rank==3 from mutations_collapsed
mut = pd.read_parquet(REPO / "cleaned_track_data/mutations_collapsed.parquet",
                      columns=["ensg_id","model_id","max_vep_rank"])
mut["model_id"] = mut["model_id"].str.lower()
mut["is_trunc"] = (mut["max_vep_rank"] == 3)
trunc = mut.groupby(["ensg_id","model_id"])["is_trunc"].any().reset_index()
before = len(pred)
pred = pred.merge(trunc, on=["ensg_id","model_id"], how="left")
print(f"[COUNT] pred after truncating merge: {len(pred):,}  (lost: {before-len(pred):,})")
pred["is_trunc"] = pred["is_trunc"].fillna(False).infer_objects(copy=False)

# Driver flags (already in pred)
pred["cna_drv"] = pred["has_cna_alteration"].fillna(False).infer_objects(copy=False)
pred["mut_drv"] = pred["p_mutation"].fillna(0) >= 0.50
pred["fus_drv"] = pred["p_fusion"].fillna(0)   >= 0.50
pred["is_onc"]  = pred["gene_role"].isin(["oncogene","both"])
pred["is_tsg"]  = pred["gene_role"].isin(["tsg","both"])

# GAIN / LOSS
gain_m = (
    (pred["cna_drv"] & pred["is_onc"]) |
    (pred["mut_drv"] & pred["is_onc"] & ~pred["is_trunc"]) |
    (pred["fus_drv"] & pred["is_onc"])
)
loss_m = (
    (pred["cna_drv"] & pred["is_tsg"]) |
    (pred["mut_drv"] & pred["is_tsg"] & pred["is_trunc"])
)
drv_m  = pred["has_driver_alteration"]
other_m = drv_m & ~gain_m & ~loss_m

print(f"\n[MEASURED] GAIN  (oncogene activation): {gain_m.sum():,}")
print(f"[MEASURED] LOSS  (TSG inactivation):    {loss_m.sum():,}")
print(f"[MEASURED] OTHER (unclassified):        {other_m.sum():,}")
print(f"[MEASURED] Total drivers:               {drv_m.sum():,}")

pred["top20"] = pred["score_rank_pct"] >= 0.80
pred["bot20"] = pred["score_rank_pct"] <  0.20

P_top = pred["top20"].mean()
P_bot = pred["bot20"].mean()
P_tg  = pred.loc[gain_m,"top20"].mean() if gain_m.sum()>0 else np.nan
P_bl  = pred.loc[loss_m,"bot20"].mean() if loss_m.sum()>0 else np.nan
E_gain = P_tg / P_top
E_loss = P_bl / P_bot if not np.isnan(P_bl) else np.nan
pooled = pred.loc[drv_m,"top20"].mean() / P_top

print(f"\n[MEASURED] P(top20%) overall    = {P_top:.6f}")
print(f"[MEASURED] P(bot20%) overall    = {P_bot:.6f}")
print(f"[MEASURED] P(top20% | GAIN)     = {P_tg:.6f}  n={gain_m.sum():,}")
print(f"[MEASURED] E_gain               = {E_gain:.4f}×")
print(f"[MEASURED] P(bot20% | LOSS)     = {P_bl:.6f}  n={loss_m.sum():,}")
print(f"[MEASURED] E_loss               = {E_loss:.4f}×")
print(f"[MEASURED] Pooled enrichment    = {pooled:.4f}×  (check vs spec 1.146×)")

# ── Fast cluster bootstrap (per-lineage sums) ──────────────────────────────────
print("\n[MEASURED] Cluster bootstrap (1000 resamples, cluster=lineage) ...")
agg = pred.groupby("lineage").apply(lambda g: pd.Series({
    "tg": g.loc[gain_m.reindex(g.index,fill_value=False),"top20"].sum(),
    "ng": gain_m.reindex(g.index,fill_value=False).sum(),
    "bl": g.loc[loss_m.reindex(g.index,fill_value=False),"bot20"].sum(),
    "nl": loss_m.reindex(g.index,fill_value=False).sum(),
    "nt": g["top20"].sum(), "nb": g["bot20"].sum(), "n":len(g),
}), include_groups=False).reset_index()
agg_np = agg.drop("lineage",axis=1).values.astype(float)
lineage_vals = agg["lineage"].values
n_lin = len(lineage_vals)
print(f"  Lineages: {n_lin}")

rng = np.random.default_rng(42)
boot_Eg, boot_El = [], []
for _ in range(1000):
    idx = rng.integers(0, n_lin, size=n_lin)
    s   = agg_np[idx].sum(axis=0)
    tg,ng,bl,nl,nt,nb,n = s
    pte = nt/n if n>0 else np.nan
    pbe = nb/n if n>0 else np.nan
    boot_Eg.append((tg/ng/pte) if ng>0 and pte>0 else np.nan)
    boot_El.append((bl/nl/pbe) if nl>0 and pbe>0 else np.nan)

boot_Eg = np.array(boot_Eg); boot_El = np.array(boot_El)
ci_Eg = (np.nanpercentile(boot_Eg,2.5), np.nanpercentile(boot_Eg,97.5))
ci_El = (np.nanpercentile(boot_El,2.5), np.nanpercentile(boot_El,97.5))
print(f"[MEASURED] E_gain 95% CI: [{ci_Eg[0]:.4f}, {ci_Eg[1]:.4f}]")
print(f"[MEASURED] E_loss 95% CI: [{ci_El[0]:.4f}, {ci_El[1]:.4f}]")

# ── PLACEBO: permute gene_role across genes (not within lineage) ───────────────
print("\n[MEASURED] PLACEBO: permuting gene_role across genes (1000 permutations) ...")
genes = pred[["ensg_id","gene_role"]].drop_duplicates("ensg_id").copy()
rng2 = np.random.default_rng(99)
plac_Eg, plac_El = [], []
roles_arr = genes["gene_role"].values.copy()
gene_ids  = genes["ensg_id"].values
for _ in range(1000):
    perm_roles = rng2.permuted(roles_arr)
    rm = dict(zip(gene_ids, perm_roles))
    rp = pred["ensg_id"].map(rm)
    onc_p = rp.isin(["oncogene","both"])
    tsg_p = rp.isin(["tsg","both"])
    gp = (pred["cna_drv"]&onc_p)|(pred["mut_drv"]&onc_p&~pred["is_trunc"])|(pred["fus_drv"]&onc_p)
    lp = (pred["cna_drv"]&tsg_p)|(pred["mut_drv"]&tsg_p&pred["is_trunc"])
    plac_Eg.append(pred.loc[gp,"top20"].mean()/P_top if gp.sum()>0 else np.nan)
    plac_El.append(pred.loc[lp,"bot20"].mean()/P_bot if lp.sum()>0 else np.nan)

pla_Eg = np.nanmean(plac_Eg); pla_El = np.nanmean(plac_El)
print(f"[MEASURED] Placebo E_gain: mean={pla_Eg:.4f}  p95={np.nanpercentile(plac_Eg,95):.4f}")
print(f"[MEASURED] Placebo E_loss: mean={pla_El:.4f}  p95={np.nanpercentile(plac_El,95):.4f}")
print(f"[MEASURED] Placebo/observed ratio (gain): {pla_Eg/E_gain:.3f}" if E_gain>0 else "")
print(f"[MEASURED] Placebo/observed ratio (loss): {pla_El/E_loss:.3f}" if not np.isnan(E_loss) and E_loss>0 else "")

# ── Verdict ────────────────────────────────────────────────────────────────────
gain_c = E_gain > 1.3 and ci_Eg[0] > 1.0
loss_c = not np.isnan(E_loss) and E_loss > 1.3 and ci_El[0] > 1.0
gain_pc = E_gain > 0 and pla_Eg/E_gain >= 0.50
print("\n" + "=" * 72)
print(f"VERDICT AGAINST PRE-DECLARATION:")
print(f"  Pre-declared: E_gain>1.3 AND E_loss>1.3, CIs exclude 1.0")
print(f"  E_gain={E_gain:.4f} [{ci_Eg[0]:.4f},{ci_Eg[1]:.4f}]  placebo={pla_Eg:.4f}")
print(f"  E_loss={E_loss:.4f} [{ci_El[0]:.4f},{ci_El[1]:.4f}]  placebo={pla_El:.4f}")
if gain_c and loss_c:   print("  CONFIRMED")
elif gain_c:             print("  PARTIAL — E_gain confirmed; E_loss NOT confirmed")
elif loss_c:             print("  PARTIAL — E_loss confirmed; E_gain NOT confirmed")
else:                    print("  NOT CONFIRMED")
if gain_pc:             print("  WARNING: PLACEBO CONTAMINATED (gain >= 50%)")
print("=" * 72)
