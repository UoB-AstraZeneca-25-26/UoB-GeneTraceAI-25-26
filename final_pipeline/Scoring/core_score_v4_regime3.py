"""
Scoring/core_score_v4_regime3.py  (PROMOTION CANDIDATE -- builds on the
currently-live core_score.py, which already has the protein-residual
shrinkage+clip fix. This adds ONE additional, independently-verified fix on
top: regime 3 ("genuine disagreement") substitution.)
-----------------------------------------------------------------------------
Regime 3 detection (matches the regime-discovery task's bucket-5 criteria
exactly, re-confirmed against current live inputs):
    |rna_std_z| > 1.0  AND  |prot_std_z| > 1.0  AND  sign(rna_std_z) != sign(prot_std_z)
  i.e. both arms are confidently signaled, and they disagree in direction.

For these rows only, core_z is computed as Variant C's direct average
(rna_std_z and prot_std_z both independently standardized, no
residualization) instead of A's default (rna_std_z + prot_resid_z). This is
mechanistically safe specifically because this regime requires protein to
ALSO be confidently signaled -- the "C lets RNA run unchecked" failure mode
(seen when C was tried broadly, or as a whole-bucket protein-neutral fix)
cannot occur here, since protein is never neutral in this regime by
construction.

All other rows (regimes 1/2/4, both-quiet, n_layers=1, RNA-only genes) are
computed IDENTICALLY to live core_score.py -- verified bit-identical in the
promotion check.

Adds `combination_variant` audit column: "A" (default) or "C_regime3"
(substituted). n_layers=1 and RNA-only-gene rows are always "A".

Verified (scratchpad regime_ab_test.py, re-confirmed here against live
inputs unchanged since):
  population: 75,117 / 6,762,902 n_layers=2 rows (1.111%)
  extremity:  A=0.2045 -> C_regime3=0.1594  (reduced)
  spike (n_layers=2): 3.0389% -> 3.0313%  (-0.0076 pp)

Reads from:  final_pipeline/outputs/{bulk_rna_z, bulk_prot_z, protein_platform_tier}.parquet
             celllineselector.db  (lineage, harmonised_enriched)
Writes to:   final_pipeline/outputs/{core_score_v4_regime3, gene_dispersion_v4_regime3}.parquet
"""
import sys, time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from config import RNA_Z, PROT_Z, PROT_TIER, CORE_SCORE, GENE_DISP
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "utils"))
import common as C

OUT_SCORE = CORE_SCORE.with_name("core_score_v4_regime3.parquet")
OUT_DISP  = GENE_DISP.with_name("gene_dispersion_v4_regime3.parquet")

W_RNA  = np.sqrt(1.431)
W_PROT = np.sqrt(1.45)
W_NORM = np.sqrt(W_RNA**2 + W_PROT**2)

RHO_PRIOR = 0.353
N_MIN_RHO = 30
K_RHO     = 30

RNA_STRATUM_MIN_N = 5
RNA_STRATUM_N0     = 25
RNA_WINSOR_BOUND   = 5.0

PROT_RESID_N0            = 25
PROT_RESID_WINSOR_BOUND  = 5.0

# ── Regime 3 ("genuine disagreement") criteria -- matches regime-discovery ──
REGIME3_SIG_THRESH = 1.0   # |rna_std_z| and |prot_std_z| both > this = "signaled"

# Protein's own independent standardization (Variant C's algorithm, needed
# only to DETECT regime 3 and to compute its substitute score -- prot_std_z
# is not used anywhere else, live core_score.py's prot_resid_z is untouched).
PROT_STD_STRATUM_MIN_N = 5
PROT_STD_STRATUM_N0    = 25
PROT_STD_WINSOR_BOUND  = 5.0


def _standardize_prot_resid_by_gene(prot_resid_wide: pd.DataFrame,
                                    n_by_gene: pd.Series) -> pd.DataFrame:
    med_g = prot_resid_wide.median(axis=0)
    sd_g  = prot_resid_wide.std(ddof=1, axis=0).clip(lower=1e-6)
    all_resid = prot_resid_wide.stack()
    sd_global = max(all_resid.std(ddof=1), 1e-6)
    n = n_by_gene.reindex(prot_resid_wide.columns).fillna(0)
    w = n / (n + PROT_RESID_N0)
    sd_final = np.sqrt(w * sd_g**2 + (1 - w) * sd_global**2)
    z = prot_resid_wide.subtract(med_g, axis=1).div(sd_final, axis=1)
    z = z.clip(lower=-PROT_RESID_WINSOR_BOUND, upper=PROT_RESID_WINSOR_BOUND)
    return z


def _standardize_rna_by_stratum(rna_long: pd.DataFrame, nsrc_col: str) -> pd.DataFrame:
    grp_cols = ["gene_id", nsrc_col]
    stratum = rna_long.groupby(grp_cols)["rna_z"].agg(
        n="size", med="median", sd=lambda x: x.std(ddof=1)).reset_index()
    gene_global = rna_long.groupby("gene_id")["rna_z"].agg(
        med_g="median", sd_g=lambda x: x.std(ddof=1)).reset_index()
    stratum = stratum.merge(gene_global, on="gene_id", how="left")
    stratum["sd_g"] = stratum["sd_g"].fillna(1e-6).clip(lower=1e-6)
    stratum["sd"]   = stratum["sd"].fillna(stratum["sd_g"]).clip(lower=1e-6)
    stratum["med_g"] = stratum["med_g"].fillna(stratum["med"])
    w    = stratum["n"] / (stratum["n"] + RNA_STRATUM_N0)
    thin = stratum["n"] < RNA_STRATUM_MIN_N
    blended_var = w * stratum["sd"] ** 2 + (1 - w) * stratum["sd_g"] ** 2
    stratum["sd_final"]  = np.where(thin, np.sqrt(blended_var), stratum["sd"])
    stratum["med_final"] = np.where(
        thin, w * stratum["med"] + (1 - w) * stratum["med_g"], stratum["med"])
    out = rna_long.merge(stratum[grp_cols + ["med_final", "sd_final"]],
                         on=grp_cols, how="left")
    out["rna_std_z"] = (out["rna_z"] - out["med_final"]) / out["sd_final"]
    out["rna_std_z"] = out["rna_std_z"].clip(-RNA_WINSOR_BOUND, RNA_WINSOR_BOUND)
    return out


def _standardize_by_stratum_generic(long_df: pd.DataFrame, value_col: str, nsrc_col: str) -> pd.DataFrame:
    """Same algorithm as _standardize_rna_by_stratum, generic value/stratum
    column -- used only to compute prot_std_z for regime-3 detection, matching
    core_score_variant_C_no_residualization.py's standardize_by_stratum()."""
    grp_cols = ["gene_id", nsrc_col]
    stratum = long_df.groupby(grp_cols)[value_col].agg(
        n="size", med="median", sd=lambda x: x.std(ddof=1)).reset_index()
    gene_global = long_df.groupby("gene_id")[value_col].agg(
        med_g="median", sd_g=lambda x: x.std(ddof=1)).reset_index()
    stratum = stratum.merge(gene_global, on="gene_id", how="left")
    stratum["sd_g"] = stratum["sd_g"].fillna(1e-6).clip(lower=1e-6)
    stratum["sd"]   = stratum["sd"].fillna(stratum["sd_g"]).clip(lower=1e-6)
    stratum["med_g"] = stratum["med_g"].fillna(stratum["med"])
    w    = stratum["n"] / (stratum["n"] + PROT_STD_STRATUM_N0)
    thin = stratum["n"] < PROT_STD_STRATUM_MIN_N
    blended_var = w * stratum["sd"] ** 2 + (1 - w) * stratum["sd_g"] ** 2
    stratum["sd_final"]  = np.where(thin, np.sqrt(blended_var), stratum["sd"])
    stratum["med_final"] = np.where(thin, w*stratum["med"] + (1-w)*stratum["med_g"], stratum["med"])
    out = long_df.merge(stratum[grp_cols + ["med_final", "sd_final"]], on=grp_cols, how="left")
    out["std_z"] = ((out[value_col] - out["med_final"]) / out["sd_final"]).clip(
        -PROT_STD_WINSOR_BOUND, PROT_STD_WINSOR_BOUND)
    return out


def run():
    t0 = time.time()
    print("=" * 70)
    print("Scoring — core_score v4 (regime-3 direct-average substitution)")
    print("=" * 70)

    for p in [RNA_Z, PROT_Z, PROT_TIER]:
        if not p.exists():
            raise SystemExit(f"Missing input: {p}\nRun prior stages first.")

    con         = C.connect()
    lineage_map = C.load_lineage(con)
    con.close()

    rna  = pd.read_parquet(RNA_Z).rename(columns={"z_t": "rna_z"})
    prot = pd.read_parquet(PROT_Z).rename(columns={"z_t": "prot_z"})
    tier = pd.read_parquet(PROT_TIER)

    print(f"RNA  rows: {len(rna):,}  |  genes: {rna.gene_id.nunique():,}")
    print(f"Prot rows: {len(prot):,}  |  genes: {prot.gene_id.nunique():,}")

    nsrc_col = "n_sources" if "n_sources" in rna.columns else "k_src"
    rna_std  = _standardize_rna_by_stratum(rna, nsrc_col)
    rna_std_z_wide = rna_std.pivot(index="model_id", columns="gene_id", values="rna_std_z")
    print(f"  rna_std_z SD by {nsrc_col} (post-standardise, target 1.00): "
          f"{rna_std.groupby(nsrc_col)['rna_std_z'].std(ddof=1).to_dict()}")

    # ── prot_std_z: Variant C's independent standardization, for regime-3
    # detection and substitution ONLY -- does not replace the residual arm ──
    prot_nsrc_col = "n_sources" if "n_sources" in prot.columns else "k_src"
    prot_std = _standardize_by_stratum_generic(prot, "prot_z", prot_nsrc_col)
    prot_std_z_wide = prot_std.pivot(index="model_id", columns="gene_id", values="std_z")
    print(f"  prot_std_z SD by {prot_nsrc_col} (post-standardise, target 1.00): "
          f"{prot_std.groupby(prot_nsrc_col)['std_z'].std(ddof=1).to_dict()}")

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
        n   = len(grp)
        lam = K_RHO / (n + K_RHO)
        rho_g = lam * RHO_PRIOR + (1 - lam) * rho_raw
        sd_rna  = grp["rna_z"].std(ddof=1)
        sd_prot = grp["prot_z"].std(ddof=1)
        beta_g  = rho_g * (sd_prot / sd_rna) if sd_rna > 0 else rho_g
        rho_rows.append({"gene_id": g, "rho_g": rho_g, "beta_g": beta_g, "n_shared": n})
    rho_df = pd.DataFrame(rho_rows)

    _sd_rna_global  = merged["rna_z"].std(ddof=1)
    _sd_prot_global = merged["prot_z"].std(ddof=1)
    _beta_fallback  = (RHO_PRIOR * (_sd_prot_global / _sd_rna_global)
                       if _sd_rna_global > 0 else RHO_PRIOR)

    merged = merged.merge(rho_df, on="gene_id", how="left")
    merged["rho_g"]  = merged["rho_g"].fillna(RHO_PRIOR)
    merged["beta_g"] = merged["beta_g"].fillna(_beta_fallback)
    merged["prot_resid"] = merged["prot_z"] - merged["beta_g"] * merged["rna_z"]

    prot_resid_wide = merged.pivot(index="model_id", columns="gene_id", values="prot_resid")

    n_by_gene = merged.groupby("gene_id").size()
    prot_resid_z_wide = _standardize_prot_resid_by_gene(prot_resid_wide, n_by_gene)
    print(f"  prot_resid_z SD (post-standardise): "
          f"mean={prot_resid_z_wide.std(ddof=1).mean():.4f}  (target 1.00)")

    # ── Orthogonal combination (default, Variant A) ─────────────────────────
    shared_genes   = list(set(rna_std_z_wide.columns) & set(prot_resid_z_wide.columns))
    shared_models  = list(set(rna_std_z_wide.index)   & set(prot_resid_z_wide.index))

    rna_m  = rna_std_z_wide.loc[shared_models, shared_genes].values
    prot_resid_m = prot_resid_z_wide.loc[shared_models, shared_genes].values
    core_z_A = (W_RNA * rna_m + W_PROT * prot_resid_m) / W_NORM

    # prot_std_z aligned to the same (shared_models, shared_genes) grid --
    # not every shared gene/model necessarily has a prot_std_z value (thin
    # strata), reindex with NaN fill so regime-3 detection safely evaluates
    # to False wherever prot_std_z is unavailable.
    prot_std_m = prot_std_z_wide.reindex(index=shared_models, columns=shared_genes).values

    # ── Regime 3 detection + substitution ───────────────────────────────────
    with np.errstate(invalid="ignore"):
        rna_sig  = np.abs(rna_m) > REGIME3_SIG_THRESH
        prot_sig = np.abs(prot_std_m) > REGIME3_SIG_THRESH
        opp_sign = np.sign(rna_m) != np.sign(prot_std_m)
    is_regime3 = rna_sig & prot_sig & opp_sign & ~np.isnan(prot_std_m)

    core_z_C = (W_RNA * rna_m + W_PROT * np.nan_to_num(prot_std_m)) / W_NORM
    core_z = np.where(is_regime3, core_z_C, core_z_A)

    core_pct = pd.DataFrame(stats.norm.cdf(core_z), index=shared_models, columns=shared_genes)
    variant_grid = pd.DataFrame(
        np.where(is_regime3, "C_regime3", "A"), index=shared_models, columns=shared_genes)

    n_regime3 = int(is_regime3.sum())
    print(f"n_layers=2 shared genes: {len(shared_genes):,}  |  shared models: {len(shared_models):,}")
    print(f"Regime 3 rows substituted (direct average, opposite-sign confident disagreement): "
          f"{n_regime3:,} / {is_regime3.size:,} ({n_regime3/is_regime3.size*100:.4f}%)")

    core_long = core_pct.stack().reset_index()
    core_long.columns = ["model_id", "ensg_id", "core_score"]
    core_long["lineage"] = lineage_map.reindex(core_long["model_id"]).values
    core_long["n_layers"] = 2
    variant_long = variant_grid.stack().reset_index()
    variant_long.columns = ["model_id", "ensg_id", "combination_variant"]
    core_long = core_long.merge(variant_long, on=["model_id","ensg_id"], how="left")

    # ── n_layers=1 fallback (identical to live core_score.py) ───────────────
    shared_models_set = set(shared_models)
    rna_models_no_prot = [m for m in rna_std_z_wide.index if m not in shared_models_set]
    if rna_models_no_prot and shared_genes:
        fallback_pct = pd.DataFrame(
            stats.norm.cdf(rna_std_z_wide.loc[rna_models_no_prot, shared_genes].values),
            index=rna_models_no_prot, columns=shared_genes,
        )
        fallback_long = fallback_pct.stack().reset_index()
        fallback_long.columns = ["model_id", "ensg_id", "core_score"]
        fallback_long["lineage"] = lineage_map.reindex(fallback_long["model_id"]).values
        fallback_long["n_layers"] = 1
        fallback_long["combination_variant"] = "A"
        core_long = pd.concat([core_long, fallback_long], ignore_index=True)
        print(f"n_layers=1 fallback (protein-free lines): {len(rna_models_no_prot):,} lines, "
              f"{len(fallback_long):,} pairs")

    rna_only_genes = [g for g in rna_std_z_wide.columns if g not in prot_resid_z_wide.columns]
    if rna_only_genes:
        rna_only_pct = pd.DataFrame(
            stats.norm.cdf(rna_std_z_wide[rna_only_genes].values),
            index=rna_std_z_wide.index, columns=rna_only_genes,
        )
        rna_only_long = rna_only_pct.stack().reset_index()
        rna_only_long.columns = ["model_id", "ensg_id", "core_score"]
        rna_only_long["lineage"] = lineage_map.reindex(rna_only_long["model_id"]).values
        rna_only_long["n_layers"] = 1
        rna_only_long["combination_variant"] = "A"
        core_long = pd.concat([core_long, rna_only_long], ignore_index=True)
        print(f"RNA-only genes added: {len(rna_only_genes):,}")

    core_long = core_long.sort_values(
        ["ensg_id", "lineage", "core_score", "model_id"],
        ascending=[True, True, False, True],
    ).reset_index(drop=True)
    core_long["stratum_rank"] = (
        core_long.groupby(["ensg_id", "lineage"])["core_score"]
        .rank(ascending=False, method="first")
    )

    disp = core_long.groupby("ensg_id")["core_score"].agg(
        mean="mean", std="std", cv=lambda x: x.std() / x.mean() if x.mean() > 0 else 0
    ).reset_index()

    core_long.to_parquet(OUT_SCORE, index=False)
    disp.to_parquet(OUT_DISP, index=False)

    elapsed = time.time() - t0
    print(f"\nDone in {elapsed:.0f}s")
    print(f"core_score: {len(core_long):,} rows  ->  {OUT_SCORE}")
    print(f"  genes: {core_long.ensg_id.nunique():,}  |  lines: {core_long.model_id.nunique():,}")
    print(f"  core_score mean: {core_long.core_score.mean():.3f}")
    print(f"  combination_variant counts: {core_long.combination_variant.value_counts().to_dict()}")


if __name__ == "__main__":
    run()
