"""
Scoring/core_score_variant_hybrid.py  (TEST HARNESS -- not the live scorer,
not wired into run_all.py.)
-----------------------------------------------------------------------------
Row-selective hybrid: Variant A (protein residualized against RNA) is the
default for every row. Variant C's score (no residualization) is substituted
ONLY for rows matching a data-derived version of the KRT86 failure pattern --
protein's own direct reading is near-neutral, yet residualization has
manufactured a large deviation from it.

Selection rule (derived in Step 1 of the cowork task, not guessed):
    |prot_std_z| < 0.15   AND   |prot_resid_z| > 1.5
  - prot_std_z:   protein's OWN direct per-stratum standardization (what
                  Variant C would use) -- near 0 means the raw, unadjusted
                  protein reading is unremarkable for this gene/cell line.
  - prot_resid_z: protein residualized against RNA, then standardized (what
                  Variant A actually uses) -- large means residualization has
                  turned that unremarkable reading into an extreme one.
  Empirically: n=8,500 rows (0.126% of the n_layers=2 population), 59.4%
  precision against the top-1%-most-divergent-from-Variant-C rows (vs 1%
  base rate -- ~59x enrichment), 42.0% of flagged rows show a genuinely large
  score change (score_diff>0.3), only 8.5% barely move (score_diff<0.05).
  See docs/SCORING_METHODS_FULL.md / the cowork task report for the full
  derivation, including the honest caveat that swapping to C does NOT
  uniformly make these rows "safer" -- C often defers entirely to RNA's own
  magnitude (unchecked by protein) and can be MORE extreme than A, not less,
  for a meaningful share of flagged rows. The aggregate spike test (this
  script's actual output) is the real arbiter, not this rule's precision
  alone.

Per Step 2.4: uses the ORIGINAL (non-Stage-2-shrunk) protein input for both
the A and C portions -- isolated deliberately from the Stage 2 question,
which was already shown to hurt both designs when combined pipeline-wide.

Adds `combination_variant` ("A" or "C") to every row for audit/transparency.

Writes to:  final_pipeline/outputs/core_score_hybrid.parquet
"""
import sys, time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from config import RNA_Z, PROT_Z, CORE_SCORE
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "utils"))
import common as C

OUT = CORE_SCORE.with_name("core_score_hybrid.parquet")

W_RNA  = np.sqrt(1.431)
W_PROT = np.sqrt(1.45)
W_NORM = np.sqrt(W_RNA**2 + W_PROT**2)

RHO_PRIOR = 0.353
N_MIN_RHO = 30
K_RHO     = 30

STRATUM_MIN_N = 5
STRATUM_N0    = 25
RESID_N0      = 25
WINSOR_BOUND  = 5.0

# ---- Step 1's derived selection rule -------------------------------------
PROT_STD_NEUTRAL_BAND = 0.15   # |prot_std_z| below this = "protein says nothing"
PROT_RESID_MANUFACTURED_THRESH = 1.5   # |prot_resid_z| above this = "A disagrees anyway"


def standardize_by_stratum(long_df: pd.DataFrame, value_col: str, nsrc_col: str) -> pd.DataFrame:
    grp_cols = ["gene_id", nsrc_col]
    stratum = long_df.groupby(grp_cols)[value_col].agg(
        n="size", med="median", sd=lambda x: x.std(ddof=1)).reset_index()
    gene_global = long_df.groupby("gene_id")[value_col].agg(
        med_g="median", sd_g=lambda x: x.std(ddof=1)).reset_index()
    stratum = stratum.merge(gene_global, on="gene_id", how="left")
    stratum["sd_g"] = stratum["sd_g"].fillna(1e-6).clip(lower=1e-6)
    stratum["sd"]   = stratum["sd"].fillna(stratum["sd_g"]).clip(lower=1e-6)
    stratum["med_g"] = stratum["med_g"].fillna(stratum["med"])
    w    = stratum["n"] / (stratum["n"] + STRATUM_N0)
    thin = stratum["n"] < STRATUM_MIN_N
    blended_var = w * stratum["sd"] ** 2 + (1 - w) * stratum["sd_g"] ** 2
    stratum["sd_final"]  = np.where(thin, np.sqrt(blended_var), stratum["sd"])
    stratum["med_final"] = np.where(thin, w*stratum["med"] + (1-w)*stratum["med_g"], stratum["med"])
    out = long_df.merge(stratum[grp_cols + ["med_final", "sd_final"]], on=grp_cols, how="left")
    out["std_z"] = ((out[value_col] - out["med_final"]) / out["sd_final"]).clip(-WINSOR_BOUND, WINSOR_BOUND)
    return out


def standardize_resid_by_gene(resid_wide: pd.DataFrame, n_by_gene: pd.Series) -> pd.DataFrame:
    med_g = resid_wide.median(axis=0)
    sd_g  = resid_wide.std(ddof=1, axis=0).clip(lower=1e-6)
    sd_global = max(resid_wide.stack().std(ddof=1), 1e-6)
    n = n_by_gene.reindex(resid_wide.columns).fillna(0)
    w = n / (n + RESID_N0)
    sd_final = np.sqrt(w * sd_g**2 + (1 - w) * sd_global**2)
    z = resid_wide.subtract(med_g, axis=1).div(sd_final, axis=1)
    return z.clip(lower=-WINSOR_BOUND, upper=WINSOR_BOUND)


def run():
    t0 = time.time()
    print("=" * 70)
    print("Hybrid -- Variant A default, Variant C for the KRT86 failure pattern")
    print("=" * 70)

    con = C.connect()
    lineage_map = C.load_lineage(con)
    con.close()

    rna  = pd.read_parquet(RNA_Z).rename(columns={"z_t": "rna_z"})
    prot = pd.read_parquet(PROT_Z).rename(columns={"z_t": "prot_z"})
    rna_nsrc_col = "n_sources" if "n_sources" in rna.columns else "k_src"
    print(f"RNA rows: {len(rna):,}  |  Prot rows: {len(prot):,}")

    # ---- rna_std_z: identical in every variant, needed for combination and
    #      for the n_layers=1 fallback -----------------------------------
    rna_std = standardize_by_stratum(rna, "rna_z", rna_nsrc_col)
    rna_std_z_wide = rna_std.pivot(index="model_id", columns="gene_id", values="std_z")
    print(f"  rna_std_z SD by {rna_nsrc_col} (target 1.00): "
          f"{rna_std.groupby(rna_nsrc_col)['std_z'].std(ddof=1).to_dict()}")

    # ---- prot_std_z: Variant C's direct protein standardization ---------
    prot_std = standardize_by_stratum(prot, "prot_z", "n_sources")
    prot_std_z_wide = prot_std.pivot(index="model_id", columns="gene_id", values="std_z")
    print(f"  prot_std_z SD by n_sources (target 1.00): "
          f"{prot_std.groupby('n_sources')['std_z'].std(ddof=1).to_dict()}")

    # ---- prot_resid_z: Variant A's residualized protein signal -----------
    genes_both = list(set(rna.gene_id) & set(prot.gene_id))
    merged = prot[prot.gene_id.isin(genes_both)].merge(
        rna[rna.gene_id.isin(genes_both)][["gene_id","model_id","rna_z"]],
        on=["gene_id","model_id"], how="inner"
    )
    rho_rows = []
    for g, grp in merged.groupby("gene_id"):
        if len(grp) < N_MIN_RHO:
            continue
        rho_raw, _ = stats.pearsonr(grp["prot_z"].fillna(0), grp["rna_z"].fillna(0))
        n = len(grp)
        lam = K_RHO / (n + K_RHO)
        rho_g = lam*RHO_PRIOR + (1-lam)*rho_raw
        sd_rna, sd_prot = grp["rna_z"].std(ddof=1), grp["prot_z"].std(ddof=1)
        beta_g = rho_g * (sd_prot/sd_rna) if sd_rna > 0 else rho_g
        rho_rows.append({"gene_id": g, "beta_g": beta_g, "n_shared": n})
    rho_df = pd.DataFrame(rho_rows)

    sd_rna_global, sd_prot_global = merged["rna_z"].std(ddof=1), merged["prot_z"].std(ddof=1)
    beta_fallback = RHO_PRIOR*(sd_prot_global/sd_rna_global) if sd_rna_global > 0 else RHO_PRIOR

    merged = merged.merge(rho_df, on="gene_id", how="left")
    merged["beta_g"] = merged["beta_g"].fillna(beta_fallback)
    merged["prot_resid"] = merged["prot_z"] - merged["beta_g"]*merged["rna_z"]

    prot_resid_wide = merged.pivot(index="model_id", columns="gene_id", values="prot_resid")
    n_by_gene = merged.groupby("gene_id").size()
    prot_resid_z_wide = standardize_resid_by_gene(prot_resid_wide, n_by_gene)
    print(f"  prot_resid_z SD (post-standardise, target 1.00): "
          f"mean={prot_resid_z_wide.std(ddof=1).mean():.4f}")

    # ---- both combinations, n_layers=2 -----------------------------------
    shared_genes  = list(set(rna_std_z_wide.columns) & set(prot_std_z_wide.columns) & set(prot_resid_z_wide.columns))
    shared_models = list(set(rna_std_z_wide.index)   & set(prot_std_z_wide.index)   & set(prot_resid_z_wide.index))
    print(f"n_layers=2 shared genes: {len(shared_genes):,}  |  shared models: {len(shared_models):,}")

    rna_m       = rna_std_z_wide.loc[shared_models, shared_genes].values
    prot_std_m  = prot_std_z_wide.loc[shared_models, shared_genes].values
    prot_resid_m = prot_resid_z_wide.loc[shared_models, shared_genes].values

    core_z_A = (W_RNA*rna_m + W_PROT*prot_resid_m) / W_NORM
    core_z_C = (W_RNA*rna_m + W_PROT*prot_std_m)   / W_NORM
    score_A = stats.norm.cdf(core_z_A)
    score_C = stats.norm.cdf(core_z_C)

    # ---- Step 1's selection rule, applied per (gene, model) -------------
    select_C = (np.abs(prot_std_m) < PROT_STD_NEUTRAL_BAND) & (np.abs(prot_resid_m) > PROT_RESID_MANUFACTURED_THRESH)
    n_selected = select_C.sum()
    print(f"Rows matching the hybrid selection rule (Variant C substituted): "
          f"{n_selected:,} / {select_C.size:,} ({n_selected/select_C.size*100:.4f}%)")

    hybrid_score = np.where(select_C, score_C, score_A)
    core_pct = pd.DataFrame(hybrid_score, index=shared_models, columns=shared_genes)
    variant_used = pd.DataFrame(np.where(select_C, "C", "A"), index=shared_models, columns=shared_genes)

    core_long = core_pct.stack().reset_index()
    core_long.columns = ["model_id", "ensg_id", "core_score"]
    variant_long = variant_used.stack().reset_index()
    variant_long.columns = ["model_id", "ensg_id", "combination_variant"]
    core_long = core_long.merge(variant_long, on=["model_id","ensg_id"], how="left")
    core_long["lineage"] = lineage_map.reindex(core_long["model_id"]).values
    core_long["n_layers"] = 2

    # ---- n_layers=1 fallback / RNA-only genes: always Variant A's plain
    #      rna_std_z path -- identical across every variant tested so far --
    shared_models_set = set(shared_models)
    rna_models_no_prot = [m for m in rna_std_z_wide.index if m not in shared_models_set]
    if rna_models_no_prot and shared_genes:
        fb = pd.DataFrame(stats.norm.cdf(rna_std_z_wide.loc[rna_models_no_prot, shared_genes].values),
                          index=rna_models_no_prot, columns=shared_genes)
        fb_long = fb.stack().reset_index()
        fb_long.columns = ["model_id", "ensg_id", "core_score"]
        fb_long["combination_variant"] = "A"   # RNA-only, no protein arm to select between
        fb_long["lineage"] = lineage_map.reindex(fb_long["model_id"]).values
        fb_long["n_layers"] = 1
        core_long = pd.concat([core_long, fb_long], ignore_index=True)
        print(f"n_layers=1 fallback: {len(rna_models_no_prot):,} lines, {len(fb_long):,} pairs")

    shared_genes_set = set(shared_genes)
    rna_only_genes = [g for g in rna_std_z_wide.columns if g not in shared_genes_set]
    if rna_only_genes:
        ro = pd.DataFrame(stats.norm.cdf(rna_std_z_wide[rna_only_genes].values),
                          index=rna_std_z_wide.index, columns=rna_only_genes)
        ro_long = ro.stack().reset_index()
        ro_long.columns = ["model_id", "ensg_id", "core_score"]
        ro_long["combination_variant"] = "A"
        ro_long["lineage"] = lineage_map.reindex(ro_long["model_id"]).values
        ro_long["n_layers"] = 1
        core_long = pd.concat([core_long, ro_long], ignore_index=True)
        print(f"RNA-only genes added: {len(rna_only_genes):,}")

    core_long = core_long.sort_values(
        ["ensg_id", "lineage", "core_score", "model_id"], ascending=[True, True, False, True]
    ).reset_index(drop=True)
    core_long["stratum_rank"] = core_long.groupby(["ensg_id","lineage"])["core_score"].rank(ascending=False, method="first")

    core_long.to_parquet(OUT, index=False)
    print(f"\nDone in {time.time()-t0:.0f}s -> {OUT}")
    print(f"rows: {len(core_long):,}  genes: {core_long.ensg_id.nunique():,}  lines: {core_long.model_id.nunique():,}")
    print(f"combination_variant counts: {core_long.combination_variant.value_counts().to_dict()}")


if __name__ == "__main__":
    run()
