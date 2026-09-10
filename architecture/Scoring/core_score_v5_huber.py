"""
Scoring/core_score_v5_huber.py  (TEST HARNESS -- not live, not wired into
run_all.py, does not modify or write to any live file.)
-----------------------------------------------------------------------------
Tests whether replacing the per-gene rho_g/beta_g OLS fit with a Huber
M-estimator (robust regression) makes core_score.py insensitive to the kind
of upstream protein-scale shift that broke three independently-tested
upstream fixes (Stage-2 shrinkage x Variant A, Stage-2 shrinkage x Variant C,
lineage-level partial-pooling shrinkage in zscore_platform() -- all three
traced to the same root cause: core_score.py refits rho_g/beta_g fresh via
plain Pearson correlation / OLS from whatever protein scale arrives, so a
well-motivated upstream rescaling reshapes the reference frame this file
measures every gene against).

Design choice (see run()'s modified block for the exact mechanism): a Huber
M-estimator (statsmodels.robust.robust_linear_model.RLM, HuberT norm --
already a project dependency, no new package added) regresses prot_z on
rna_z per gene, replacing BOTH halves of the live OLS-derived beta_g:
  1. the correlation (rho_raw, previously plain Pearson) -- now a robust
     correlation derived from the Huber slope via robust (MAD-based) scale
     ratios, fed into the EXISTING EB-shrinkage blend unchanged
     (rho_g = lam*RHO_PRIOR + (1-lam)*rho_robust, lam/K_RHO/N_MIN_RHO
     untouched) -- preserving that logic rather than redesigning it, per
     the task's explicit preference.
  2. the SD ratio in beta_g = rho_g * (sd_prot/sd_rna) -- the live code uses
     RAW std(ddof=1) here, which is exactly the OTHER half of the corruption
     mechanism already traced for test 3 (a single outlier inflated
     sd_prot from ~0.64 to 16.09 for one gene). Leaving this raw would let
     contamination back in even with a robust correlation, defeating the
     point of this test -- so this SD ratio is ALSO made robust (MAD-based)
     here, for both the per-gene path and the global fallback path (for
     genes below N_MIN_RHO), for internal consistency.

Everything else in core_score.py is unchanged: RNA standardization, the
protein-residual shrinkage+clip mechanism itself, regime-3 detection and
substitution, final weighted combination and CDF.

Reads from:  final_pipeline/outputs/{bulk_rna_z, bulk_prot_z, protein_platform_tier}.parquet
             celllineselector.db  (lineage, harmonised_enriched)
Writes to:   final_pipeline/outputs/{core_score_v5_huber, gene_dispersion_v5_huber}.parquet
             (overridable via module-level OUT_SCORE/OUT_DISP, same
             monkeypatch pattern used for the lineage-shrinkage downstream test)
"""
import sys, time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
import statsmodels.api as sm

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from config import RNA_Z, PROT_Z, PROT_TIER, CORE_SCORE, GENE_DISP
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "utils"))
import common as C

OUT_SCORE = CORE_SCORE.with_name("core_score_v5_huber.parquet")
OUT_DISP  = GENE_DISP.with_name("gene_dispersion_v5_huber.parquet")

MAD_CONST = 1.4826  # MAD -> SD scaling constant for a normal, matches project convention


def _mad_scale(x: np.ndarray) -> float:
    x = x[np.isfinite(x)]
    if x.size == 0:
        return np.nan
    med = np.median(x)
    mad = np.median(np.abs(x - med))
    return MAD_CONST * mad if mad > 1e-9 else (np.std(x, ddof=1) if x.size > 1 else np.nan)


def _huber_beta_raw(prot: np.ndarray, rna: np.ndarray) -> tuple[float, float]:
    """Robust regression of prot on rna via Huber's M-estimator.
    Returns (beta_huber, rho_robust) -- rho_robust derived from the Huber
    slope via MAD-based (robust) scale ratios, NOT raw std, so contamination
    from a small number of extreme points can't leak back in through the
    scale conversion the way it did in the live OLS path."""
    mask = np.isfinite(prot) & np.isfinite(rna)
    prot, rna = prot[mask], rna[mask]
    if len(prot) < 3 or np.std(rna) < 1e-9:
        return 0.0, 0.0
    X = sm.add_constant(rna)
    try:
        # update_scale=False: hold the initial (OLS-residual-based MAD) scale
        # estimate fixed through IRLS rather than letting it update each
        # iteration. DISCOVERED, NOT ASSUMED: with the default update_scale=
        # True, the scale estimate collapsed toward 0 for the MAJORITY of
        # genes tested in this project's data (confirmed empirically: 39/40
        # sampled well-behaved genes on a spot check), triggering
        # statsmodels' own "Estimated scale is 0.0 ... perfect fit" warning
        # and producing nonsense slopes (e.g. 1e-61) as IRLS weights
        # collapsed to ~0 for most points, keeping only a handful with
        # near-exact-zero residuals. This is a known failure mode of Huber's
        # default MAD-of-residuals scale estimator under data with many
        # tied/near-tied residuals (plausible here: fillna(0) and thin-
        # stratum fallbacks produce real ties in these z-scores). Holding
        # the scale fixed at its initial (pre-iteration) estimate avoids
        # this collapse and was verified to restore sane slopes on both the
        # well-behaved sample and the two known-contaminated genes.
        res = sm.RLM(prot, X, M=sm.robust.norms.HuberT()).fit(update_scale=False)
        beta_huber = float(res.params[1])
    except Exception:
        return 0.0, 0.0
    sd_rna_r  = _mad_scale(rna)
    sd_prot_r = _mad_scale(prot)
    if not np.isfinite(sd_rna_r) or not np.isfinite(sd_prot_r) or sd_prot_r < 1e-9:
        return beta_huber, 0.0
    rho_robust = float(np.clip(beta_huber * (sd_rna_r / sd_prot_r), -1.0, 1.0))
    return beta_huber, rho_robust


W_RNA  = np.sqrt(1.431)   # Kish n_eff: N=3 sources, ρ_avg=0.548 → n_eff=3/(1+2×0.548)=1.431; replaces provisional sqrt(2.00)
W_PROT = np.sqrt(1.45)    # 2 sources, ProCAN-CCLE rho=0.373 -> n_eff = 2/(1+0.373)
W_NORM = np.sqrt(W_RNA**2 + W_PROT**2)

RHO_PRIOR = 0.353   # median per-gene r(prot_z, rna_z) — RNA-protein EB prior
N_MIN_RHO = 30      # minimum shared lines to estimate per-gene rho
K_RHO     = 30      # EB pseudo-count: prior/data get equal weight at n = K_RHO

# ── RNA standardization (per-gene × n_sources, thin-stratum shrinkage) ──────
# Replaces the prior within-lineage percentile + clipped norm.ppf treatment.
# That approach double-conditioned RNA (already lineage-conditioned upstream
# in rna_scorer.py) and structurally suppressed RNA magnitude relative to
# protein's uncapped (x-median)/SD treatment.
#
# Per-n_sources SD is genuinely heterogeneous on the live RNA output (measured:
# 0.88/1.54/1.86 raw for n_sources=1/2/3), and a naive per-gene-global
# standardization does NOT fix this -- per-n_sources SD stays heterogeneous
# (0.65/1.08/1.35) even after "standardizing", because the heterogeneity is a
# stratum-level effect, not a gene-level scale issue. Standardizing within
# each (gene, n_sources) stratum directly targets this; confirmed via a
# worked example (ENSG00000178828) to correct real rank misordering rather
# than just relabelling the same distortion.
RNA_STRATUM_MIN_N = 5    # below this, a raw per-stratum SD is unreliably estimated
RNA_STRATUM_N0     = 25  # shrinkage pseudo-count toward the gene-global SD for thin
                         # strata -- James-Stein-style blend, same shape as
                         # transcriptomics.py's method_z() (scale = sqrt(w*stratum_var
                         # + (1-w)*global_var), w = n/(n+n0)). n0=25 mirrors that
                         # module's calibrated value as a FIRST-PASS choice only --
                         # not independently validated for this context.
RNA_WINSOR_BOUND   = 5.0  # safety cap replacing the old percentile step's implicit
                         # ±3.09 bound (norm.ppf(0.999)); observed to bind on <0.2%
                         # of rows post-standardization, so this rarely activates.

# ── Protein residual shrinkage + clip ───────────────────────────────────────
# core_score_v1_pre_rna_fix_*.py.bak diffed byte-for-byte against this file at
# the time: the RNA fix above added shrinkage AND a clip to the RNA arm; the
# protein residual standardization block was not touched at all. No doc or
# code comment anywhere stated this was a deliberate asymmetry -- confirmed
# via investigation to be an inherited gap, the same "fixed one arm, never
# revisited the other" pattern already found once in this project for the RNA
# percentile step's own history. This section closes that gap using the same
# w=n/(n+n0) shrinkage shape as RNA, continuously applied (not gated behind a
# hard "thin" cutoff the way RNA's is) -- a standalone test on known-outlier
# genes showed RNA's own MIN_N=5 hard cutoff would not even fire for the
# concrete outlier case found (n=19 shared lines, above 5), so a continuous
# blend is used here deliberately, not copied verbatim from RNA's exact
# on/off shape.
#
# "Wider reference" for protein mirrors RNA's own structure one level up:
# RNA blends a (gene, n_sources)-stratum estimate toward a gene-global one;
# protein has no stratum dimension at this stage (residuals are already
# per-gene, panel-wide), so it blends the gene-level SD toward a single
# panel-wide SD (pooling every gene's residuals together) instead.
PROT_RESID_N0            = 25   # SAME first-pass value as RNA_STRATUM_N0, reused
                                # by convention, NOT independently re-derived or
                                # calibrated for this context. Flagged as a
                                # first-pass choice needing its own validation.
PROT_RESID_WINSOR_BOUND  = 5.0  # matches RNA_WINSOR_BOUND. Chosen because
                                # prot_resid_z is standardized to the same
                                # target (SD=1.00, verified by the printed
                                # post-standardise check) as rna_std_z, so a
                                # 5-sigma bound is an equally rare event on
                                # either arm under that shared calibration --
                                # not blindly copied without checking the
                                # underlying scale matched.

# ── Regime 3 ("genuine disagreement") detection + substitution ─────────────
REGIME3_SIG_THRESH = 1.0   # |rna_std_z| and |prot_std_z| both > this = "signaled";
                          # matches the regime-discovery task's bucket-5
                          # threshold exactly (not independently re-tuned here).

# prot_std_z is computed ONLY to detect/substitute regime 3 -- it does not
# replace prot_resid_z anywhere else. Same per-(gene, n_sources) stratum
# shrinkage shape as RNA's, matching Variant C's standardize_by_stratum().
PROT_STD_STRATUM_MIN_N = 5
PROT_STD_STRATUM_N0    = 25
PROT_STD_WINSOR_BOUND  = 5.0


def _standardize_prot_resid_by_gene(prot_resid_wide: pd.DataFrame,
                                    n_by_gene: pd.Series) -> pd.DataFrame:
    """(model_id x gene_id) prot_resid -> prot_resid_z, with the gene's own
    residual SD shrunk toward the panel-wide residual SD, weighted by how
    many shared RNA/protein lines that gene had (w = n/(n+PROT_RESID_N0)),
    continuously for every gene (see module docstring for why this is
    continuous rather than RNA's hard-cutoff shape). Centre (median) is left
    as the gene's own -- only the scale, the actual source of the blow-up
    identified in the investigation, is shrunk. Result is winsorized at
    +/-PROT_RESID_WINSOR_BOUND.
    """
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
    """(gene_id, model_id, rna_z, <nsrc_col>) -> same + rna_std_z.

    Per-(gene_id, nsrc_col) median/SD standardization, with a shrinkage blend
    toward the gene-global median/SD for strata with fewer than
    RNA_STRATUM_MIN_N contributing lines (not a full fallback to the
    gene-global standardization -- that would reintroduce the exact
    suppression/amplification problem this replaces, just for a smaller,
    ~1.2-1.3% population of genes). Result is winsorized at ±RNA_WINSOR_BOUND.
    """
    grp_cols = ["gene_id", nsrc_col]

    stratum = rna_long.groupby(grp_cols)["rna_z"].agg(
        n="size", med="median", sd=lambda x: x.std(ddof=1)).reset_index()
    gene_global = rna_long.groupby("gene_id")["rna_z"].agg(
        med_g="median", sd_g=lambda x: x.std(ddof=1)).reset_index()
    stratum = stratum.merge(gene_global, on="gene_id", how="left")

    # sd_g can itself be NaN/degenerate for single-row genes; floor both scales
    # so division is always well-defined (mirrors MAD_FLOOR-style flooring
    # used elsewhere in this project, at a much smaller, purely-defensive scale
    # since this is a last-resort floor, not a calibrated minimum spread).
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
    column. Used only to compute prot_std_z for regime-3 detection and
    substitution -- matches core_score_variant_C_no_residualization.py's
    standardize_by_stratum() exactly, kept as an independent copy here (not
    a shared import) per this project's self-contained-script convention."""
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
    print("Scoring — core_score (RNA + Protein orthogonal combination)")
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

    # ── RNA: per-(gene, n_sources) standardization, winsorized ─────────────
    nsrc_col = "n_sources" if "n_sources" in rna.columns else "k_src"
    rna_std  = _standardize_rna_by_stratum(rna, nsrc_col)
    rna_std_z_wide = rna_std.pivot(index="model_id", columns="gene_id", values="rna_std_z")
    print(f"  rna_std_z SD by {nsrc_col} (post-standardise, target 1.00): "
          f"{rna_std.groupby(nsrc_col)['rna_std_z'].std(ddof=1).to_dict()}")

    # ── Protein: per-gene independent standardization, for regime-3 only ───
    prot_nsrc_col = "n_sources" if "n_sources" in prot.columns else "k_src"
    prot_std = _standardize_by_stratum_generic(prot, "prot_z", prot_nsrc_col)
    prot_std_z_wide = prot_std.pivot(index="model_id", columns="gene_id", values="std_z")
    print(f"  prot_std_z SD by {prot_nsrc_col} (post-standardise, target 1.00, "
          f"regime-3 detection only): {prot_std.groupby(prot_nsrc_col)['std_z'].std(ddof=1).to_dict()}")

    # ── Protein: per-gene empirical Bayes ρ blend → residualise against RNA ─
    genes_both = list(set(rna.gene_id) & set(prot.gene_id))
    merged = prot[prot.gene_id.isin(genes_both)].merge(
        rna[rna.gene_id.isin(genes_both)][["gene_id","model_id","rna_z"]],
        on=["gene_id","model_id"], how="inner"
    )

    rho_rows = []
    for g, grp in merged.groupby("gene_id"):
        if len(grp) < N_MIN_RHO:
            continue
        # HUBER (robust) fit in place of OLS: regress prot_z on rna_z via a
        # Huber M-estimator, then derive a robust correlation from the
        # robust slope using MAD-based (not raw-std) scale ratios -- see
        # _huber_beta_raw()'s docstring for why the raw-std conversion is
        # deliberately avoided here. rho_robust is then fed into the SAME
        # EB-shrinkage blend the live code uses, unchanged.
        _, rho_robust = _huber_beta_raw(
            grp["prot_z"].fillna(0).to_numpy(dtype=float),
            grp["rna_z"].fillna(0).to_numpy(dtype=float),
        )
        n   = len(grp)
        lam = K_RHO / (n + K_RHO)   # decreases with n: more data → less shrinkage
        rho_g = lam * RHO_PRIOR + (1 - lam) * rho_robust
        # β = ρ × SD(Y)/SD(X). Live code uses RAW std(ddof=1) here -- exactly
        # the other half of the corruption mechanism traced for test 3 (a
        # single outlier inflated sd_prot to 16x its robust value for one
        # gene). Using raw std here would let contamination back in even
        # with a robust correlation, so this ratio is ALSO MAD-based.
        sd_rna_r  = _mad_scale(grp["rna_z"].to_numpy(dtype=float))
        sd_prot_r = _mad_scale(grp["prot_z"].to_numpy(dtype=float))
        beta_g  = rho_g * (sd_prot_r / sd_rna_r) if sd_rna_r > 1e-9 else rho_g
        rho_rows.append({"gene_id": g, "rho_g": rho_g, "beta_g": beta_g, "n_shared": n})
    rho_df = pd.DataFrame(rho_rows)

    # Fallback β for genes with too few shared lines -- MAD-based global
    # scale ratio, for consistency with the per-gene path above (the live
    # code's raw-std version is less exposed to single-gene contamination
    # since it pools all 6.7M rows, but kept consistent here regardless).
    _sd_rna_global  = _mad_scale(merged["rna_z"].to_numpy(dtype=float))
    _sd_prot_global = _mad_scale(merged["prot_z"].to_numpy(dtype=float))
    _beta_fallback  = (RHO_PRIOR * (_sd_prot_global / _sd_rna_global)
                       if _sd_rna_global > 1e-9 else RHO_PRIOR)

    merged = merged.merge(rho_df, on="gene_id", how="left")
    merged["rho_g"]  = merged["rho_g"].fillna(RHO_PRIOR)
    merged["beta_g"] = merged["beta_g"].fillna(_beta_fallback)
    merged["prot_resid"] = merged["prot_z"] - merged["beta_g"] * merged["rna_z"]

    prot_resid_wide = merged.pivot(index="model_id", columns="gene_id", values="prot_resid")

    # ── Standardise protein residual: shrink SD toward panel-wide, clip ────
    # n_shared here is EVERY gene's actual shared-line count (not gated by
    # N_MIN_RHO=30 the way rho_g's own estimation is) -- a gene with, say, 10
    # shared lines still needs a shrinkage weight even though it fell back to
    # the population beta_g above; w=n/(n+n0) naturally shrinks it hard
    # (mostly toward the panel-wide SD) at that n, which is the correct
    # behaviour, not a special case to code around.
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

    # ── Regime 3 ("genuine disagreement") detection + substitution ─────────
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

    # ── Combined genes (RNA + Protein) — lines with protein data ───────────
    core_long = core_pct.stack().reset_index()
    core_long.columns = ["model_id", "ensg_id", "core_score"]
    core_long["lineage"] = lineage_map.reindex(core_long["model_id"]).values
    core_long["n_layers"] = 2
    variant_long = variant_grid.stack().reset_index()
    variant_long.columns = ["model_id", "ensg_id", "combination_variant"]
    core_long = core_long.merge(variant_long, on=["model_id","ensg_id"], how="left")

    # ── n_layers=1 FALLBACK for shared genes at protein-free lines ──────────
    # Without this, 48.6% of RNA pairs (lines not in shared_models) get no
    # score for the 10,871 shared-gene set. Per-pair fallback to RNA-only
    # for those lines, so every RNA-profiled line gets a score for every gene.
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
        fallback_long["n_layers"] = 1   # RNA-only for these lines (gene has protein, line doesn't)
        fallback_long["combination_variant"] = "A"
        core_long = pd.concat([core_long, fallback_long], ignore_index=True)
        print(f"n_layers=1 fallback (protein-free lines): {len(rna_models_no_prot):,} lines, "
              f"{len(fallback_long):,} pairs")

    # ── RNA-only genes (not in any proteomics platform) ─────────────────────
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

    # ── Stratum rank (within lineage) ──────────────────────────────────────
    # core_long's row order is not guaranteed stable across process runs:
    # shared_genes/shared_models above are built via set intersection, whose
    # iteration order depends on Python's per-process string hash
    # randomization. rank(method="first") breaks ties by row order, so an
    # unstable input order made stratum_rank non-deterministic on tied
    # core_score values (measured: 37.1% of rows differed between two runs on
    # identical data). (ensg_id, model_id) is unique in core_long, so sorting
    # on model_id as an explicit final key removes all remaining ties
    # regardless of incoming row order. Rule: ties in core_score are broken
    # by model_id, ascending.
    core_long = core_long.sort_values(
        ["ensg_id", "lineage", "core_score", "model_id"],
        ascending=[True, True, False, True],
    ).reset_index(drop=True)
    core_long["stratum_rank"] = (
        core_long.groupby(["ensg_id", "lineage"])["core_score"]
        .rank(ascending=False, method="first")
    )

    # ── Gene dispersion ─────────────────────────────────────────────────────
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
