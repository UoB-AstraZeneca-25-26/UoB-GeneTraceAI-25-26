"""
Scoring/core_score_v9_nan_fix.py  (CANDIDATE -- builds on the live, C3/C7/
regime-3-fixed core_score.py, unchanged except for one targeted addition:
dropping phantom NaN rows that currently survive `.stack()`.)
-----------------------------------------------------------------------------
NaN-row-survival fix (this task, not yet promoted):

`shared_models`/`shared_genes` (and, identically, the two fallback blocks'
own index/column choices) are each built from the FULL row/column universe
of a wide-pivoted table, not the true per-(model, gene) overlap of real
measurements. This is by itself harmless IF a cell lacking real data
produces NaN and that NaN then drops out before the row is written --which
is exactly what `.stack()` used to do by default. In the pandas version this
project is currently running (3.0.1), `.stack()`'s default implementation no
longer drops NaN rows (confirmed empirically: `.stack(dropna=True)` now
raises `ValueError: dropna must be unspecified as the new implementation
does not introduce rows of NA values` under the current default
`future_stack=True` -- ironic phrasing, since it refers only to
alignment-introduced NaNs, not pre-existing NaN cell values, which this
grid genuinely has). The result: rows with NO real underlying RNA and/or
protein measurement survive into the final output with `core_score=NaN`,
while `n_layers`/`combination_variant` are stamped unconditionally as if the
row were fully evidenced.

Traced and quantified in a prior diagnostic task: 9,188,147 / 26,558,496
rows (34.60%) of the live `core_score.parquet` are such phantom rows,
concentrated almost entirely in low-RNA-coverage genes (up to 97%+ phantom
for the least-covered decile, <0.3% for the best-covered). `core_score`
itself was never fabricated -- it was correctly NaN throughout the traced
example (ach-001333 x CD86) -- the bug is row *survival*, not score
*computation*.

Fix: `.dropna()` immediately on the stacked Series, before `.reset_index()`,
at all three sites that build a row-per-(model,gene) frame from a dense
grid:
  1. the main RNA+protein combination block
  2. the n_layers=1 protein-free-line fallback
  3. the n_layers=1 RNA-only-gene fallback
Chosen over `.stack(future_stack=False, dropna=True)` (the old
implementation) because pandas explicitly flags that old code path itself
for removal in a future version ("Pandas4Warning: The previous
implementation of stack is deprecated") -- using it would just trade one
version-fragile assumption for another. `.stack().dropna()` works
identically regardless of which `.stack()` implementation pandas defaults
to, since it operates on the already-produced Series, not on a `.stack()`
constructor argument.

This changes ONLY which rows survive -- every row with a real, non-NaN
`core_score` is completely unaffected; the arithmetic that produces
`core_score` is untouched. Verification (this task): row count now
17,370,349 (== the previously-confirmed "real" count); zero remaining NaN
core_score rows; bit-identical values on every row that was already real.

Writes to:  final_pipeline/outputs/core_score_v9_nan_fix.parquet,
            final_pipeline/outputs/gene_dispersion_v9_nan_fix.parquet
(candidate outputs -- NOT the live core_score.parquet/gene_dispersion.parquet)
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

OUT_CORE = CORE_SCORE.with_name("core_score_v9_nan_fix.parquet")
OUT_DISP = GENE_DISP.with_name("gene_dispersion_v9_nan_fix.parquet")

W_RNA  = np.sqrt(1.431)
W_PROT = np.sqrt(1.45)
W_NORM = np.sqrt(W_RNA**2 + W_PROT**2)

RHO_PRIOR = 0.353
N_MIN_RHO = 30
K_RHO     = 30

SD_SHRINK_N0 = 25

RNA_STRATUM_MIN_N = 5
RNA_STRATUM_N0     = 25
RNA_WINSOR_BOUND   = 5.0

PROT_RESID_N0            = 25
PROT_RESID_WINSOR_BOUND  = 5.0

REGIME3_SIG_THRESH = 1.0

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
    print("Scoring — core_score v9 (NaN-row-survival fix, candidate)")
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

    prot_nsrc_col = "n_sources" if "n_sources" in prot.columns else "k_src"
    prot_std = _standardize_by_stratum_generic(prot, "prot_z", prot_nsrc_col)
    prot_std_z_wide = prot_std.pivot(index="model_id", columns="gene_id", values="std_z")
    print(f"  prot_std_z SD by {prot_nsrc_col} (post-standardise, target 1.00, "
          f"regime-3 detection only): {prot_std.groupby(prot_nsrc_col)['std_z'].std(ddof=1).to_dict()}")

    genes_both = list(set(rna.gene_id) & set(prot.gene_id))
    merged = prot[prot.gene_id.isin(genes_both)].merge(
        rna[rna.gene_id.isin(genes_both)][["gene_id","model_id","rna_z"]],
        on=["gene_id","model_id"], how="inner"
    )

    _sd_rna_global  = merged["rna_z"].std(ddof=1)
    _sd_prot_global = merged["prot_z"].std(ddof=1)

    rho_rows = []
    for g, grp in merged.groupby("gene_id"):
        if len(grp) < N_MIN_RHO:
            continue
        rho_raw, _ = stats.pearsonr(grp["prot_z"].fillna(0), grp["rna_z"].fillna(0))
        n   = len(grp)
        lam = K_RHO / (n + K_RHO)
        rho_g = lam * RHO_PRIOR + (1 - lam) * rho_raw
        sd_rna_raw  = grp["rna_z"].std(ddof=1)
        sd_prot_raw = grp["prot_z"].std(ddof=1)
        w_sd = n / (n + SD_SHRINK_N0)
        sd_rna  = np.sqrt(w_sd * sd_rna_raw**2  + (1 - w_sd) * _sd_rna_global**2)
        sd_prot = np.sqrt(w_sd * sd_prot_raw**2 + (1 - w_sd) * _sd_prot_global**2)
        beta_g  = rho_g * (sd_prot / sd_rna) if sd_rna > 0 else rho_g
        rho_rows.append({"gene_id": g, "rho_g": rho_g, "beta_g": beta_g, "n_shared": n})
    rho_df = pd.DataFrame(rho_rows)

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

    shared_genes   = list(set(rna_std_z_wide.columns) & set(prot_resid_z_wide.columns))
    shared_models  = list(set(rna_std_z_wide.index)   & set(prot_resid_z_wide.index))

    rna_m  = rna_std_z_wide.loc[shared_models, shared_genes].values
    prot_resid_m = prot_resid_z_wide.loc[shared_models, shared_genes].values
    core_z_A = (W_RNA * rna_m + W_PROT * prot_resid_m) / W_NORM

    prot_std_m = prot_std_z_wide.reindex(index=shared_models, columns=shared_genes).values

    with np.errstate(invalid="ignore"):
        rna_sig  = np.abs(rna_m) > REGIME3_SIG_THRESH
        prot_sig = np.abs(prot_std_m) > REGIME3_SIG_THRESH
        opp_sign = np.sign(rna_m) != np.sign(prot_std_m)
    is_regime3 = rna_sig & prot_sig & opp_sign & ~np.isnan(prot_std_m)

    core_z_C = (W_RNA * rna_m + W_PROT * np.nan_to_num(prot_std_m)) / W_NORM
    core_z = np.where(is_regime3, core_z_C, core_z_A)

    core_pct = pd.DataFrame(
        stats.norm.cdf(core_z),
        index=shared_models, columns=shared_genes,
    )
    variant_grid = pd.DataFrame(
        np.where(is_regime3, "C_regime3", "A"), index=shared_models, columns=shared_genes)

    n_regime3 = int(is_regime3.sum())
    print(f"Regime 3 rows substituted (direct average, opposite-sign confident disagreement): "
          f"{n_regime3:,} / {is_regime3.size:,} ({n_regime3/is_regime3.size*100:.4f}%)")

    # ── NaN-row-survival fix, site 1/3: main combination ────────────────────
    # .dropna() on the stacked Series itself, before reset_index() -- works
    # regardless of which .stack() implementation pandas defaults to.
    n_before = core_pct.size
    core_long = core_pct.stack().dropna().reset_index()
    core_long.columns = ["model_id", "ensg_id", "core_score"]
    n_dropped = n_before - len(core_long)
    print(f"[NaN-fix site 1/main] {n_dropped:,} / {n_before:,} phantom NaN rows dropped "
          f"({n_dropped/n_before*100:.2f}%)")
    core_long["lineage"] = lineage_map.reindex(core_long["model_id"]).values
    core_long["n_layers"] = 2
    variant_long = variant_grid.stack().reset_index()
    variant_long.columns = ["model_id", "ensg_id", "combination_variant"]
    core_long = core_long.merge(variant_long, on=["model_id","ensg_id"], how="left")

    # ── n_layers=1 FALLBACK for shared genes at protein-free lines ──────────
    shared_models_set = set(shared_models)
    rna_models_no_prot = [m for m in rna_std_z_wide.index if m not in shared_models_set]
    if rna_models_no_prot and shared_genes:
        fallback_core_z = (W_RNA / W_NORM) * rna_std_z_wide.loc[rna_models_no_prot, shared_genes].values
        fallback_pct = pd.DataFrame(
            stats.norm.cdf(fallback_core_z),
            index=rna_models_no_prot, columns=shared_genes,
        )
        # ── NaN-row-survival fix, site 2/3: protein-free-line fallback ──────
        n_before_fb = fallback_pct.size
        fallback_long = fallback_pct.stack().dropna().reset_index()
        fallback_long.columns = ["model_id", "ensg_id", "core_score"]
        n_dropped_fb = n_before_fb - len(fallback_long)
        print(f"[NaN-fix site 2/fallback] {n_dropped_fb:,} / {n_before_fb:,} phantom NaN rows dropped "
              f"({n_dropped_fb/n_before_fb*100:.2f}%)")
        fallback_long["lineage"] = lineage_map.reindex(fallback_long["model_id"]).values
        fallback_long["n_layers"] = 1
        fallback_long["combination_variant"] = "A_fallback_scaled"
        core_long = pd.concat([core_long, fallback_long], ignore_index=True)
        print(f"n_layers=1 fallback (protein-free lines, C7-scaled): {len(rna_models_no_prot):,} lines, "
              f"{len(fallback_long):,} real pairs")

    # ── RNA-only genes (not in any proteomics platform) ─────────────────────
    rna_only_genes = [g for g in rna_std_z_wide.columns if g not in prot_resid_z_wide.columns]
    if rna_only_genes:
        rna_only_core_z = (W_RNA / W_NORM) * rna_std_z_wide[rna_only_genes].values
        rna_only_pct = pd.DataFrame(
            stats.norm.cdf(rna_only_core_z),
            index=rna_std_z_wide.index, columns=rna_only_genes,
        )
        # ── NaN-row-survival fix, site 3/3: RNA-only-gene fallback ──────────
        n_before_ro = rna_only_pct.size
        rna_only_long = rna_only_pct.stack().dropna().reset_index()
        rna_only_long.columns = ["model_id", "ensg_id", "core_score"]
        n_dropped_ro = n_before_ro - len(rna_only_long)
        print(f"[NaN-fix site 3/rna-only-genes] {n_dropped_ro:,} / {n_before_ro:,} phantom NaN rows dropped "
              f"({n_dropped_ro/n_before_ro*100:.2f}%)")
        rna_only_long["lineage"] = lineage_map.reindex(rna_only_long["model_id"]).values
        rna_only_long["n_layers"] = 1
        rna_only_long["combination_variant"] = "A_fallback_scaled"
        core_long = pd.concat([core_long, rna_only_long], ignore_index=True)
        print(f"RNA-only genes added (C7-scaled): {len(rna_only_genes):,} genes, "
              f"{len(rna_only_long):,} real pairs")

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

    core_long.to_parquet(OUT_CORE, index=False)
    disp.to_parquet(OUT_DISP, index=False)

    elapsed = time.time() - t0
    print(f"\nDone in {elapsed:.0f}s")
    print(f"core_score: {len(core_long):,} rows  ->  {OUT_CORE}")
    print(f"  genes: {core_long.ensg_id.nunique():,}  |  lines: {core_long.model_id.nunique():,}")
    print(f"  core_score mean: {core_long.core_score.mean():.3f}")
    print(f"  remaining NaN core_score rows: {core_long.core_score.isna().sum():,}  (expect 0)")
    print(f"  combination_variant counts: {core_long.combination_variant.value_counts().to_dict()}")


if __name__ == "__main__":
    run()
