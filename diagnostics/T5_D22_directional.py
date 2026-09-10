"""
T5 — D22: Stratified directional enrichment.
Run from repo root: python diagnostics/T5_D22_directional.py
"""
import hashlib, sys
from pathlib import Path
import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parent.parent
OUT  = REPO / "diagnostics" / "out"
sys.path.insert(0, str(REPO / "final_pipeline"))
from config import PREDICTIONS, CNA_FLAGS, MUTATIONS_SCR, FUSIONS_SCR, DB
import duckdb

def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()

print("=" * 72)
print("T5 — D22: STRATIFIED DIRECTIONAL ENRICHMENT")
print("=" * 72)

pred = pd.read_parquet(PREDICTIONS)
print(f"[COUNT] predictions: {len(pred):,}")
print(f"[COUNT] has_driver_alteration TRUE: {pred.has_driver_alteration.sum():,}")
print(f"[COUNT] gene_role values: {pred.gene_role.value_counts().to_dict()}")

# Join lineage
con = duckdb.connect(str(DB), read_only=True)
lin = con.execute("SELECT model_id, lineage FROM main.sample_info WHERE lineage IS NOT NULL").df()
lin = lin.drop_duplicates("model_id")
con.close()
pred = pred.merge(lin, on="model_id", how="left")
print(f"[COUNT] After lineage join, missing: {pred.lineage.isna().sum():,}")

pred["top20"]   = pred["score_rank_pct"] >= 0.80
pred["bot20"]   = pred["score_rank_pct"] <  0.20

# Load mutation flags to get truncating flag
mut_raw = pd.read_parquet(REPO / "cleaned_track_data/mutations_collapsed.parquet",
                           columns=["ensg_id","model_id","any_driver","multi_hit_high_impact"])
# "truncating" proxy: check if max_vep_rank is 3 (stop_gained, frameshift etc.)
# We'll use the raw data
mut_raw2 = pd.read_parquet(REPO / "cleaned_track_data/mutations_collapsed.parquet")
print(f"[COUNT] mutations_collapsed cols: {list(mut_raw2.columns)}")

mut_raw2["model_id"] = mut_raw2["model_id"].str.lower()
# Check for truncating column
has_trunc = "is_truncating" in mut_raw2.columns or "truncating" in mut_raw2.columns
print(f"[DOCUMENTED] truncating column present: {has_trunc}")
if "max_vep_rank" in mut_raw2.columns:
    print(f"[MEASURED] max_vep_rank values: {mut_raw2.max_vep_rank.value_counts().sort_index().to_dict()}")

# Load CNA data
cna = pd.read_parquet(CNA_FLAGS, columns=["model_id","ensg_id","has_cna_alteration",
                                           "is_amplification","is_deletion"]
                      if all(c in pd.read_parquet(CNA_FLAGS, columns=["model_id"]).columns or True
                             for c in ["is_amplification","is_deletion"]) else
                      ["model_id","ensg_id","has_cna_alteration"])
# Check available CNA columns
cna_cols = pd.read_parquet(CNA_FLAGS).columns.tolist()
print(f"[DOCUMENTED] CNA flags columns: {cna_cols}")
cna = pd.read_parquet(CNA_FLAGS)

# Join CNA call into predictions
pred = pred.merge(cna[["model_id","ensg_id","has_cna_alteration"]].rename(
    columns={"has_cna_alteration":"cna_driver"}), on=["model_id","ensg_id"], how="left")
pred["cna_driver"] = pred["cna_driver"].fillna(False)

# Fusions
fus = pd.read_parquet(FUSIONS_SCR, columns=["model_id","ensg_id","p_fusion"])
pred = pred.merge(fus, on=["model_id","ensg_id"], how="left")
pred["fusion_driver"] = pred["p_fusion"].fillna(0) >= 0.50

# Mutation drivers
muts = pd.read_parquet(MUTATIONS_SCR)
pred = pred.merge(muts[["model_id","ensg_id","p_mutation"]], on=["model_id","ensg_id"], how="left")
pred["mut_driver_flag"] = pred["p_mutation"].fillna(0) >= 0.50

# Check if we have truncating info in mutations
# max_vep_rank >= 3 could proxy truncating (rank 3 = splice acceptor/donor, 4/5 = stop_gained/frameshift)
# Use max_vep_rank from mutations_collapsed
mut_join = mut_raw2[["ensg_id","model_id","max_vep_rank"]].copy()
mut_join["model_id"] = mut_join["model_id"].str.lower()
mut_join["is_truncating"] = mut_join["max_vep_rank"].isin([4, 5]) if "max_vep_rank" in mut_join.columns else False
# Aggregate to (gene, line)
mut_trunc = mut_join.groupby(["ensg_id","model_id"])["is_truncating"].any().reset_index()
pred = pred.merge(mut_trunc, on=["ensg_id","model_id"], how="left")
pred["is_truncating"] = pred["is_truncating"].fillna(False)

# ── Partition into GAIN / LOSS / OTHER ────────────────────────────────────────
# GAIN: amplification of oncogene OR non-truncating mut of oncogene OR fusion of oncogene
pred["is_oncogene"] = pred["gene_role"].isin(["oncogene","both"])
pred["is_tsg"]      = pred["gene_role"].isin(["tsg","both"])

gain_mask = (
    (pred["cna_driver"] & pred["is_oncogene"] & pred["has_cna_alteration"]) |
    (pred["mut_driver_flag"] & pred["is_oncogene"] & ~pred["is_truncating"]) |
    (pred["fusion_driver"] & pred["is_oncogene"])
)
loss_mask = (
    (pred["cna_driver"] & pred["is_tsg"]) |
    (pred["mut_driver_flag"] & pred["is_tsg"] & pred["is_truncating"])
)
driver_mask = pred["has_driver_alteration"]
other_mask  = driver_mask & ~gain_mask & ~loss_mask

N_gain  = gain_mask.sum()
N_loss  = loss_mask.sum()
N_other = other_mask.sum()
N_driver = driver_mask.sum()
print(f"\n[MEASURED] Driver partition:")
print(f"  GAIN  (oncogene activation): {N_gain:,}")
print(f"  LOSS  (TSG inactivation):    {N_loss:,}")
print(f"  OTHER (unclassified):        {N_other:,}")
print(f"  Total drivers:               {N_driver:,}")
print(f"  GAIN+LOSS+OTHER:             {N_gain+N_loss+N_other:,}  (vs drivers={N_driver:,})")

P_top20_overall = pred["top20"].mean()
P_bot20_overall = pred["bot20"].mean()
print(f"\n[MEASURED] P(top20%) overall: {P_top20_overall:.6f}")
print(f"[MEASURED] P(bot20%) overall: {P_bot20_overall:.6f}")

P_top20_gain = pred.loc[gain_mask,"top20"].mean() if N_gain > 0 else np.nan
P_bot20_loss = pred.loc[loss_mask,"bot20"].mean() if N_loss > 0 else np.nan
E_gain = P_top20_gain / P_top20_overall if P_top20_overall > 0 else np.nan
E_loss = P_bot20_loss / P_bot20_overall if P_bot20_overall > 0 else np.nan

print(f"\n[MEASURED] P(top20% | GAIN) = {P_top20_gain:.6f}")
print(f"[MEASURED] E_gain = P(top20%|GAIN)/P(top20%) = {E_gain:.4f}")
print(f"[MEASURED] P(bot20% | LOSS) = {P_bot20_loss:.6f}")
print(f"[MEASURED] E_loss = P(bot20%|LOSS)/P(bot20%) = {E_loss:.4f}")

# Pooled recompute for consistency
P_top20_driver = pred.loc[driver_mask,"top20"].mean()
pooled_enrich  = P_top20_driver / P_top20_overall
print(f"\n[MEASURED] Pooled enrichment (all drivers, top20%): {pooled_enrich:.4f}×")

# ── Cluster bootstrap on lineage (fast, per-lineage aggregation) ──────────────
lineages = pred["lineage"].dropna().unique()
n_lin    = len(lineages)
assert n_lin >= 20, f"INSUFFICIENT CLUSTERS: {n_lin}"
print(f"\n[MEASURED] Lineages: {n_lin}")
print("[MEASURED] Cluster bootstrap (1000 resamples) ...")

agg = pred.groupby("lineage").apply(lambda g: pd.Series({
    "top_gain": g.loc[gain_mask.reindex(g.index,fill_value=False),"top20"].sum(),
    "n_gain":   gain_mask.reindex(g.index,fill_value=False).sum(),
    "bot_loss": g.loc[loss_mask.reindex(g.index,fill_value=False),"bot20"].sum(),
    "n_loss":   loss_mask.reindex(g.index,fill_value=False).sum(),
    "n_top20":  g["top20"].sum(),
    "n_total":  len(g),
    "n_bot20":  g["bot20"].sum(),
}), include_groups=False).reset_index()

rng = np.random.default_rng(42)
boot_Eg, boot_El = [], []
for _ in range(1000):
    sampled = rng.choice(lineages, size=n_lin, replace=True)
    cnt = pd.Series(sampled).value_counts().reset_index()
    cnt.columns = ["lineage","mult"]
    tmp = agg.merge(cnt, on="lineage")
    t = (tmp * tmp["mult"].values.reshape(-1,1)).sum()
    p_top_g = t["top_gain"] / t["n_gain"] if t["n_gain"] > 0 else np.nan
    p_top_o = t["n_top20"] / t["n_total"] if t["n_total"] > 0 else np.nan
    p_bot_l = t["bot_loss"] / t["n_loss"] if t["n_loss"] > 0 else np.nan
    p_bot_o = t["n_bot20"] / t["n_total"] if t["n_total"] > 0 else np.nan
    boot_Eg.append(p_top_g / p_top_o if p_top_o else np.nan)
    boot_El.append(p_bot_l / p_bot_o if p_bot_o else np.nan)

boot_Eg = np.array(boot_Eg); boot_El = np.array(boot_El)
print(f"[MEASURED] E_gain 95% CI: [{np.nanpercentile(boot_Eg,2.5):.4f}, {np.nanpercentile(boot_Eg,97.5):.4f}]")
print(f"[MEASURED] E_loss 95% CI: [{np.nanpercentile(boot_El,2.5):.4f}, {np.nanpercentile(boot_El,97.5):.4f}]")

# ── MANDATORY PLACEBO: permute gene_role within lineage ──────────────────────
print("\n[MEASURED] PLACEBO: permuting gene_role within lineage ...")
rng2 = np.random.default_rng(99)
placebo_Eg, placebo_El = [], []
for _ in range(1000):
    # shuffle gene_role labels within each lineage
    pred_p = pred.copy()
    for lin_val in lineages:
        mask_l = pred_p["lineage"] == lin_val
        pred_p.loc[mask_l,"gene_role"] = rng2.permuted(pred_p.loc[mask_l,"gene_role"].values)
    pred_p["is_onc_p"] = pred_p["gene_role"].isin(["oncogene","both"])
    pred_p["is_tsg_p"] = pred_p["gene_role"].isin(["tsg","both"])
    gain_p = (pred_p["mut_driver_flag"] & pred_p["is_onc_p"] & ~pred_p["is_truncating"]) | \
             (pred_p["cna_driver"] & pred_p["is_onc_p"]) | (pred_p["fusion_driver"] & pred_p["is_onc_p"])
    loss_p = (pred_p["cna_driver"] & pred_p["is_tsg_p"]) | \
             (pred_p["mut_driver_flag"] & pred_p["is_tsg_p"] & pred_p["is_truncating"])
    pt_g = pred_p.loc[gain_p,"top20"].mean() if gain_p.sum()>0 else np.nan
    pb_l = pred_p.loc[loss_p,"bot20"].mean() if loss_p.sum()>0 else np.nan
    placebo_Eg.append(pt_g / P_top20_overall if P_top20_overall>0 else np.nan)
    placebo_El.append(pb_l / P_bot20_overall if P_bot20_overall>0 else np.nan)

pla_Eg = np.nanmean(placebo_Eg); pla_El = np.nanmean(placebo_El)
print(f"[MEASURED] Placebo E_gain mean:  {pla_Eg:.4f}  (observed={E_gain:.4f})")
print(f"[MEASURED] Placebo E_loss mean:  {pla_El:.4f}  (observed={E_loss:.4f})")
print(f"[MEASURED] Placebo E_gain p95:   {np.nanpercentile(placebo_Eg,95):.4f}")
print(f"[MEASURED] Placebo fraction of observed (gain): {pla_Eg/E_gain:.3f}" if E_gain>0 else "E_gain=0")
print(f"[MEASURED] Placebo fraction of observed (loss): {pla_El/E_loss:.3f}" if not np.isnan(E_loss) and E_loss>0 else "E_loss=0 or NaN")

# ── Verdict ──────────────────────────────────────────────────────────────────
ci_Eg_lo = np.nanpercentile(boot_Eg,2.5); ci_Eg_hi = np.nanpercentile(boot_Eg,97.5)
ci_El_lo = np.nanpercentile(boot_El,2.5); ci_El_hi = np.nanpercentile(boot_El,97.5)
print("\n" + "=" * 72)
print("VERDICT AGAINST PRE-DECLARATION:")
print(f"  Pre-declared: E_gain > 1.3 AND E_loss > 1.3, both CIs exclude 1.0")
print(f"  E_gain = {E_gain:.4f}  CI=[{ci_Eg_lo:.4f},{ci_Eg_hi:.4f}]  placebo={pla_Eg:.4f}")
print(f"  E_loss = {E_loss:.4f}  CI=[{ci_El_lo:.4f},{ci_El_hi:.4f}]  placebo={pla_El:.4f}")
gain_ok = not np.isnan(E_gain) and E_gain > 1.3 and ci_Eg_lo > 1.0
loss_ok = not np.isnan(E_loss) and E_loss > 1.3 and ci_El_lo > 1.0
if gain_ok and loss_ok:
    print("  CONFIRMED — both > 1.3 and CIs exclude 1.0")
elif gain_ok:
    print("  PARTIAL — E_gain confirmed; E_loss not confirmed")
elif loss_ok:
    print("  PARTIAL — E_loss confirmed; E_gain not confirmed")
else:
    print("  NOT CONFIRMED — neither E_gain nor E_loss meets threshold")
if not np.isnan(pla_Eg) and not np.isnan(E_gain) and E_gain > 0 and pla_Eg/E_gain >= 0.50:
    print("  WARNING: PLACEBO CONTAMINATED — placebo reproduces >=50% of observed E_gain")
print("=" * 72)
