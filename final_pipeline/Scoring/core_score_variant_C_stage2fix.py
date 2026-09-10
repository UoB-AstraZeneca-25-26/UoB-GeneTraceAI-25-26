"""
Scoring/core_score_variant_C_stage2fix.py  (TEST HARNESS -- not the live
scorer, not wired into run_all.py.)
-----------------------------------------------------------------------------
Variant C (no residualization) COMBINED with the already-validated Stage 2
shrinkage fix (06_protein_score_M2_v2.py -> bulk_prot_z_v2.parquet), instead
of the live bulk_prot_z.parquet.

Motivation: combining Stage 2's shrinkage with Variant A (the live formula,
protein-residualized-against-RNA) made the final core_score spike WORSE
(3.098%->3.256%), traced to core_score.py's own per-gene rho_g/beta_g
regression being refit fresh from whatever protein input it receives --
taming Stage 2's outliers shifted the regression fit enough to newly-extreme
OTHER cell lines of the same gene. Variant C has no regression step at all,
so that specific interaction cannot occur here by construction. This tests
whether Stage 2's raw-z-score improvement (|z|>=5 fraction: -68% on the full
panel) shows through cleanly when there's no regression step to absorb or
reverse it.

Otherwise identical to core_score_variant_C_no_residualization.py:
    rna_std_z  = standardize_by_stratum(rna_z,  grouped by RNA's n_sources)
    prot_std_z = standardize_by_stratum(prot_z_v2, grouped by protein's n_sources)
    core_z     = (W_rna*rna_std_z + W_prot*prot_std_z) / norm

    rna_std_z  = standardize_by_stratum(rna_z,  grouped by RNA's n_sources)
    prot_std_z = standardize_by_stratum(prot_z, grouped by protein's n_sources)
    core_z     = (W_rna*rna_std_z + W_prot*prot_std_z) / norm

No regression step at all -- by construction, immune to the OLS
outlier-sensitivity issue found in the Stage-2 investigation (there is no
beta_g/beta_g' to be distorted by outliers here).

n_layers=1 / RNA-only-gene fallback: plain rna_std_z, identical to the other
two variants.

Weights (W_rna, W_prot) unchanged from the live values.

Writes to:  final_pipeline/outputs/core_score_variant_C.parquet
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

PROT_Z_V2 = PROT_Z.with_name("bulk_prot_z_v2.parquet")   # Stage-2-shrunk input
OUT = CORE_SCORE.with_name("core_score_variant_C_stage2fix.parquet")

W_RNA  = np.sqrt(1.431)
W_PROT = np.sqrt(1.45)
W_NORM = np.sqrt(W_RNA**2 + W_PROT**2)

STRATUM_MIN_N = 5
STRATUM_N0    = 25
WINSOR_BOUND  = 5.0


def standardize_by_stratum(long_df: pd.DataFrame, value_col: str, nsrc_col: str) -> pd.DataFrame:
    """Identical algorithm to core_score.py's _standardize_rna_by_stratum and
    to variant B's copy of the same function -- kept as an independent copy
    here (rather than a shared import) so this file is fully self-contained
    and inspectable on its own, matching the other two variant scripts."""
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


def run():
    t0 = time.time()
    print("=" * 70)
    print("Variant C -- NO RESIDUALIZATION: both arms independently standardized")
    print("=" * 70)

    con = C.connect()
    lineage_map = C.load_lineage(con)
    con.close()

    rna  = pd.read_parquet(RNA_Z).rename(columns={"z_t": "rna_z"})
    prot = pd.read_parquet(PROT_Z_V2).rename(columns={"z_t": "prot_z"})
    rna_nsrc_col = "n_sources" if "n_sources" in rna.columns else "k_src"

    print(f"RNA rows: {len(rna):,}  |  Prot rows: {len(prot):,}")

    rna_std = standardize_by_stratum(rna, "rna_z", rna_nsrc_col)
    rna_std_z_wide = rna_std.pivot(index="model_id", columns="gene_id", values="std_z")
    print(f"  rna_std_z SD by {rna_nsrc_col} (target 1.00): "
          f"{rna_std.groupby(rna_nsrc_col)['std_z'].std(ddof=1).to_dict()}")

    prot_std = standardize_by_stratum(prot, "prot_z", "n_sources")
    prot_std_z_wide = prot_std.pivot(index="model_id", columns="gene_id", values="std_z")
    print(f"  prot_std_z SD by n_sources (target 1.00): "
          f"{prot_std.groupby('n_sources')['std_z'].std(ddof=1).to_dict()}")

    shared_genes  = list(set(rna_std_z_wide.columns) & set(prot_std_z_wide.columns))
    shared_models = list(set(rna_std_z_wide.index)   & set(prot_std_z_wide.index))
    rna_m  = rna_std_z_wide.loc[shared_models, shared_genes].values
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

    # See core_score_variant_B_flipped.py's comment: check against
    # shared_genes, not prot_std_z_wide.columns directly, to correctly catch
    # genes with zero actual overlapping (gene, model_id) rows after the
    # inner join (13 such genes found empirically) -- matches the live
    # formula's behaviour, which derives its equivalent check from the same
    # merged data rather than an independently-built column set.
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
