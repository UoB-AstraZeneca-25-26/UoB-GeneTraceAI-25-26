"""
Scoring/core_score_variant_B_flipped.py  (TEST HARNESS -- not the live
scorer, not wired into run_all.py. Part of the three-way residualization
comparison: see core_score_variant_A_current.py and
core_score_variant_C_no_residualization.py.)
-----------------------------------------------------------------------------
Flips the residualization direction relative to the live formula: RNA is
regressed against protein and only the RNA-independent-of-protein residual
enters the combination; protein enters as its own direct standardization.

    prot_std_z  = standardize_by_stratum(prot_z, grouped by protein's own
                  n_sources)  -- mirrors rna_std_z's treatment exactly,
                  same function, protein's data instead of RNA's.
    rna_resid   = rna_z - beta_g' * prot_z
    rna_resid_z = standardize_with_shrinkage_and_clip(rna_resid)  -- same
                  per-gene shrinkage+clip function used for prot_resid_z in
                  the live formula, applied here to the RNA residual instead.
    core_z      = (W_rna*rna_resid_z + W_prot*prot_std_z) / norm

n_layers=1 (protein-free lines) and RNA-only genes still fall back to plain
rna_std_z (RNA's own direct stratum standardization) -- residualizing RNA
against protein is meaningless when protein is absent for that line, so this
fallback is IDENTICAL across all three variants by construction, matching
what Step 2.2 requires.

Weights (W_rna, W_prot) are UNCHANGED from the live values -- not the
question this task investigates.

Writes to:  final_pipeline/outputs/core_score_variant_B.parquet
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

OUT = CORE_SCORE.with_name("core_score_variant_B.parquet")

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


def standardize_by_stratum(long_df: pd.DataFrame, value_col: str, nsrc_col: str) -> pd.DataFrame:
    """Generic version of core_score.py's _standardize_rna_by_stratum --
    identical algorithm, parameterized on which column/stratum to use, so it
    can standardize either RNA or protein directly, mirroring each other
    exactly. Per-(gene_id, nsrc_col) median/SD, shrunk toward gene-global for
    strata with <STRATUM_MIN_N lines, clipped at +/-WINSOR_BOUND."""
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
    """Generic version of core_score.py's _standardize_prot_resid_by_gene --
    identical algorithm, works on any residual matrix (RNA-on-protein here,
    protein-on-RNA in the live formula). Per-gene SD shrunk toward the
    panel-wide (all genes pooled) residual SD, w=n/(n+RESID_N0), clipped."""
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
    print("Variant B -- FLIPPED: RNA residualized against protein")
    print("=" * 70)

    con = C.connect()
    lineage_map = C.load_lineage(con)
    con.close()

    rna  = pd.read_parquet(RNA_Z).rename(columns={"z_t": "rna_z"})
    prot = pd.read_parquet(PROT_Z).rename(columns={"z_t": "prot_z"})
    rna_nsrc_col = "n_sources" if "n_sources" in rna.columns else "k_src"

    print(f"RNA rows: {len(rna):,}  |  Prot rows: {len(prot):,}")

    # ---- plain rna_std_z (needed for n_layers=1 fallback, IDENTICAL across
    #      all three variants -- residualization direction only touches
    #      n_layers=2 rows) --------------------------------------------------
    rna_std = standardize_by_stratum(rna, "rna_z", rna_nsrc_col)
    rna_std_z_wide = rna_std.pivot(index="model_id", columns="gene_id", values="std_z")
    print(f"  rna_std_z SD by {rna_nsrc_col} (target 1.00): "
          f"{rna_std.groupby(rna_nsrc_col)['std_z'].std(ddof=1).to_dict()}")

    # ---- prot_std_z (protein's own direct standardization, mirrors rna_std_z)
    prot_std = standardize_by_stratum(prot, "prot_z", "n_sources")
    prot_std_z_wide = prot_std.pivot(index="model_id", columns="gene_id", values="std_z")
    print(f"  prot_std_z SD by n_sources (target 1.00): "
          f"{prot_std.groupby('n_sources')['std_z'].std(ddof=1).to_dict()}")

    # ---- RNA residualized against protein --------------------------------
    genes_both = list(set(rna.gene_id) & set(prot.gene_id))
    merged = rna[rna.gene_id.isin(genes_both)].merge(
        prot[prot.gene_id.isin(genes_both)][["gene_id","model_id","prot_z"]],
        on=["gene_id","model_id"], how="inner"
    )
    rho_rows = []
    for g, grp in merged.groupby("gene_id"):
        if len(grp) < N_MIN_RHO:
            continue
        rho_raw, _ = stats.pearsonr(grp["rna_z"].fillna(0), grp["prot_z"].fillna(0))
        n = len(grp)
        lam = K_RHO / (n + K_RHO)
        rho_g = lam*RHO_PRIOR + (1-lam)*rho_raw
        sd_prot, sd_rna = grp["prot_z"].std(ddof=1), grp["rna_z"].std(ddof=1)
        beta_g = rho_g * (sd_rna/sd_prot) if sd_prot > 0 else rho_g
        rho_rows.append({"gene_id": g, "beta_g": beta_g, "n_shared": n})
    rho_df = pd.DataFrame(rho_rows)

    sd_rna_global, sd_prot_global = merged["rna_z"].std(ddof=1), merged["prot_z"].std(ddof=1)
    beta_fallback = RHO_PRIOR*(sd_rna_global/sd_prot_global) if sd_prot_global > 0 else RHO_PRIOR

    merged = merged.merge(rho_df, on="gene_id", how="left")
    merged["beta_g"] = merged["beta_g"].fillna(beta_fallback)
    merged["rna_resid"] = merged["rna_z"] - merged["beta_g"]*merged["prot_z"]

    rna_resid_wide = merged.pivot(index="model_id", columns="gene_id", values="rna_resid")
    n_by_gene = merged.groupby("gene_id").size()
    rna_resid_z_wide = standardize_resid_by_gene(rna_resid_wide, n_by_gene)
    print(f"  rna_resid_z SD (post-standardise, target 1.00): "
          f"mean={rna_resid_z_wide.std(ddof=1).mean():.4f}")

    # ---- combine: n_layers=2 uses (rna_resid_z, prot_std_z) ---------------
    shared_genes  = list(set(rna_resid_z_wide.columns) & set(prot_std_z_wide.columns))
    shared_models = list(set(rna_resid_z_wide.index)   & set(prot_std_z_wide.index))
    rna_m  = rna_resid_z_wide.loc[shared_models, shared_genes].values
    prot_m = prot_std_z_wide.loc[shared_models, shared_genes].values
    core_z = (W_RNA*rna_m + W_PROT*prot_m) / W_NORM
    core_pct = pd.DataFrame(stats.norm.cdf(core_z), index=shared_models, columns=shared_genes)

    core_long = core_pct.stack().reset_index()
    core_long.columns = ["model_id", "ensg_id", "core_score"]
    core_long["lineage"] = lineage_map.reindex(core_long["model_id"]).values
    core_long["n_layers"] = 2

    shared_models_set = set(shared_models)
    rna_models_no_prot = [m for m in rna_std_z_wide.index if m not in shared_models_set]
    if rna_models_no_prot and shared_genes:
        fb = pd.DataFrame(stats.norm.cdf(rna_std_z_wide.loc[rna_models_no_prot, shared_genes].values),
                          index=rna_models_no_prot, columns=shared_genes)
        fb_long = fb.stack().reset_index()
        fb_long.columns = ["model_id", "ensg_id", "core_score"]
        fb_long["lineage"] = lineage_map.reindex(fb_long["model_id"]).values
        fb_long["n_layers"] = 1
        core_long = pd.concat([core_long, fb_long], ignore_index=True)
        print(f"n_layers=1 fallback: {len(rna_models_no_prot):,} lines, {len(fb_long):,} pairs")

    # NOTE: checked against shared_genes (not prot_std_z_wide.columns directly)
    # -- a small number of genes are nominally in both gene universes but have
    # ZERO actual overlapping (gene, model_id) rows after the inner join (13
    # such genes found empirically), so they never make it into shared_genes.
    # The live formula catches this automatically because its equivalent
    # check is against a column set DERIVED from the same merge; checking
    # against prot_std_z_wide.columns directly (built independently of any
    # merge) would silently drop these genes from every variant's output --
    # verified this was happening (missing 20,319 = 13*1563 rows) before
    # this fix.
    shared_genes_set = set(shared_genes)
    rna_only_genes = [g for g in rna_std_z_wide.columns if g not in shared_genes_set]
    if rna_only_genes:
        ro = pd.DataFrame(stats.norm.cdf(rna_std_z_wide[rna_only_genes].values),
                          index=rna_std_z_wide.index, columns=rna_only_genes)
        ro_long = ro.stack().reset_index()
        ro_long.columns = ["model_id", "ensg_id", "core_score"]
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


if __name__ == "__main__":
    run()
