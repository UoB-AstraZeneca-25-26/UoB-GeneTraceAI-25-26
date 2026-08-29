"""
Scoring/core_score.py
----------------------
beta_g SD-shrinkage fix ("C3", promoted 2026-08-29): shrink sd_rna/sd_prot
before forming beta_g.

Prior to this promotion, beta_g = rho_g * (sd_prot_raw/sd_rna_raw), where
rho_g is carefully empirical-Bayes shrunk (lam=K_RHO/(n+K_RHO)) but sd_rna_raw
and sd_prot_raw were RAW, unshrunk per-gene estimates -- even for genes barely
above N_MIN_RHO=30, where a single-gene SD estimate is itself noisy. This was
the same "shrink one factor, not the other" pattern already fixed once in
this file for the protein-residual arm (_standardize_prot_resid_by_gene) --
the fourth confirmed instance of this pattern in this project (see
docs/WHY_REFERENCE.md for the other three).

Fixed by applying the IDENTICAL shrinkage shape already used elsewhere in
this file (w = n/(n+n0), blend toward the panel-wide SD) to sd_rna and
sd_prot before they enter the beta_g ratio:
    sd_final = sqrt(w * sd_raw^2 + (1-w) * sd_global^2)
sd_global is each SD's panel-wide value across the full genes_both population
(the same _sd_rna_global/_sd_prot_global already computed in this file
for the small-n beta fallback -- reused here as the shrinkage target, not
newly invented).

n0 choice: reuses SD_SHRINK_N0=25, the same value as RNA_STRATUM_N0/
PROT_RESID_N0/PROT_STD_STRATUM_N0 elsewhere in this file. This is a
DELIBERATE reuse for consistency, not an independent derivation for this
specific application (shrinking an SD that feeds a ratio inside a slope
calculation, not a direct location/scale standardization) -- flagged as
first-pass and UNVALIDATED for this context, same discipline as every other
n0 in this file. No evidence was found suggesting a different n0 would serve
this specific use better; reusing the existing convention was chosen over
inventing an unrelated new number with equally no justification.

Verified before promotion (core_score_v8_beta_shrinkage.py): isolated finding
reproduced against the real implementation -- near N_MIN_RHO=30 (30-60 shared
lines), beta_g shifted by median 19.9%, mean 72.1%, p99 1074%; far from the
threshold (n>200), shift was <1%. Downstream full-pipeline check (NOT assumed
to help just because the isolated beta_g shift was well-motivated, per this
project's established pattern of upstream improvements not always composing
well downstream -- this is the pattern's first CLEAN composition, not a fifth
failure): full-panel spike 1.4790% -> 1.4769%, n_layers=2 spike 3.0313% ->
3.0257% (both flat-to-slightly-improved, not worse); full-panel core_score
mean/quantiles unchanged to 3 decimal places (aggregate effect is small
because only ~339/11,460 qualifying genes sit near the N_MIN_RHO threshold),
while individual near-threshold genes shifted meaningfully at the per-line
level (e.g. ENSG00000137463/ach-000939: 0.8609 -> 0.8244). See
docs/WHY_REFERENCE.md for the full verification record, and
core_score_v5_pre_beta_shrinkage_fix_*.py.bak for the pre-promotion snapshot.

Full history below (unchanged from live core_score.py):
Combines pre-scored RNA (3 sources) and protein (2 platforms) into a single
core_score per (gene, cell-line), then ranks within lineage strata.

Architecture:
  RNA z-scores   <- 01_Transcriptomics/rna_scorer.py  -> bulk_rna_z.parquet
  Protein z-scores <- 02_Proteinomics/protein_scorer.py -> bulk_prot_z.parquet
  Tier filter    <- 02_Proteinomics/platform_tier.py  -> protein_platform_tier.parquet
  Metadata       <- Stage-0 warehouse (lineage, model labels)
  -> core_score.parquet, gene_dispersion.parquet

Combination (default, Variant A):
  RNA std z       = (rna_z - median_by_gene_and_nsources) / SD_by_gene_and_nsources,
                    winsorized at +/-RNA_WINSOR_BOUND. Replaces the prior
                    within-lineage percentile + clipped norm.ppf treatment,
                    which double-conditioned RNA (already lineage-conditioned
                    upstream) and left it structurally suppressed relative to
                    protein's uncapped per-gene treatment. See
                    _standardize_rna_by_stratum() for the thin-stratum
                    shrinkage detail.
  Prot residual   = prot_z_t - beta_g * rna_z_t        (OLS residual; beta_g = rho_g * SD_prot/SD_rna)
  prot_resid_z    = (prot_resid - med_gene) / SD_final_gene, winsorized at
                    +/-PROT_RESID_WINSOR_BOUND. SD_final_gene shrinks each
                    gene's own residual SD toward the panel-wide residual SD
                    (w = n_gene/(n_gene+PROT_RESID_N0)) -- the protein
                    analogue of RNA's per-stratum shrinkage above. Replaces a
                    prior version with neither shrinkage nor a bound (git-
                    diffed against core_score_v1_pre_rna_fix_*.py.bak: the RNA
                    fix added both to the RNA arm and left this block
                    untouched -- an inherited gap, not a considered decision).
                    See _standardize_prot_resid_by_gene() for detail.
  W_RNA = sqrt(1.431)  [Kish n_eff: N=3 sources, rho_avg=0.548]
  W_PROT = sqrt(1.45)  [2 sources, ProCAN-CCLE rho=0.373 -> n_eff=1.45; DISTINCT from RHO_PRIOR=0.353]
  core_z = (W_RNA * rna_std_z + W_PROT * prot_resid_z) / sqrt(W_RNA^2 + W_PROT^2)
    -- W_RNA/W_PROT left UNCHANGED by both the RNA and protein standardization
       fixes; whether they still make sense given the corrected arms is an
       explicitly separate, out-of-scope question for a later task.
  core_score = Phi(core_z)  [monotone relabelling to [0,1]]

Regime-3 substitution (promoted on top of the above, 2026-08-28):
  Regime-discovery analysis (unsupervised-free, interpretable binning of
  |rna_std_z|/|prot_std_z| against a "signaled" threshold) found 5 candidate
  regimes in the n_layers=2 population. One of them -- "genuine disagreement"
  (both arms confidently signaled, |z|>REGIME3_SIG_THRESH, pointing OPPOSITE
  directions; ~1.1% of n_layers=2 rows) -- showed A's residualization
  uniquely MORE extreme than a plain no-residualization average there
  (reversed from every other regime), and the lowest per-gene rho_g of any
  regime (weakest RNA-protein coupling, so the regression A relies on is
  least trustworthy exactly where disagreement occurs). For rows matching
  this regime ONLY, core_z is computed as the direct average of
  rna_std_z/prot_std_z (Variant C's no-residualization formula) instead of
  A's default. This is mechanistically safe here specifically because the
  regime requires protein to ALSO be confidently signaled -- the "C lets RNA
  run unchecked" failure mode (measured when C was tried pipeline-wide, and
  again when a whole-bucket RNA-only substitution was tried on a DIFFERENT
  regime, "protein-neutral") cannot occur in this regime by construction.
  Verified before promotion: bit-identical to the pre-regime3 scorer on every
  non-matching row (max diff 0.0); population 75,117/6,762,902 rows (1.111%,
  reproduced exactly); n_layers=2 spike 3.0389% -> 3.0313%; full-panel spike
  2.8473% -> 2.8444%. See docs/SCORING_METHODS_FULL.md for the regime-
  discovery task and this promotion's full verification record.
  Rows carry a `combination_variant` audit column ("A" default, "C_regime3"
  for the substituted rows; n_layers=1 and RNA-only-gene rows are always "A")
  so the substitution is auditable, not silent.

  A second regime from the same discovery task -- "protein-neutral" (RNA
  signaled, protein reads near-neutral; ~21% of n_layers=2 rows, the KRT86
  pattern) -- was tested three separate ways (a precision-targeted 8,500-row
  C-substitution; a whole-bucket RNA-only substitution; a precision-targeted
  RNA-only substitution) and NONE improved the aggregate spike -- the
  whole-bucket RNA-only attempt made it substantially worse (+1.58pp),
  because ~92% of that regime was already being correctly damped by A's
  weighted combination and had no suppression problem to fix. This regime is
  NOT addressed here; it remains a documented open limitation. See
  docs/SCORING_METHODS_FULL.md.

  Known limitation (measured, not eliminated by this fix): core_score's
  distribution has a small pile-up at exactly 0.0/1.0. Part of this is
  irreducible from the protein side: even with a hypothetically perfect
  protein arm, RNA's own +/-5 clip alone pushes core_z into Phi's saturated
  region on ~0.48% of n_layers=2 rows (rna_std_z only needs to reach ~3.30,
  well under its own clip, given W_RNA's ~0.70 fractional share of the
  combination) -- against a measured combined-arms rate of ~3.10% pre-fix,
  i.e. RNA's clip alone accounts for roughly 15% of the pile-up on its own.
  The protein shrinkage+clip fix measurably reduced the remainder
  (n_layers=2: 3.098% -> 3.039%; full panel: 2.870% -> 2.847%); the regime-3
  fix above reduces it further (3.039% -> 3.031%) but does not eliminate it.
  A follow-up investigation found that shrinking the earlier, per-(gene,
  lineage) z-scoring stage instead (02_Proteinomics/06_protein_score_M2.py's
  robust_z_matrix()) substantially improves the raw protein z-scores (|z|>=5
  fraction: -68% on the full panel) but, tested end to end, made this
  stage's OWN combined core_score spike slightly WORSE (it interacts with
  the per-gene rho_g/beta_g regression this file fits fresh from whatever
  protein input it receives) -- that fix was evaluated and explicitly NOT
  promoted; see docs/SCORING_METHODS_FULL.md.

n_layers=1 fallback scale fix ("C7", promoted 2026-08-29):
  The full mathematical audit found that the n_layers=1 fallback path
  (protein-free lines; RNA-only genes) computed core_score = Phi(rna_std_z)
  directly, with NO W_RNA/W_NORM scaling -- while n_layers=2 rows apply that
  ~0.7048 downweight to rna_std_z before Phi. The same RNA evidence therefore
  scored systematically more extreme via the fallback path than it would have
  via the combined path with a neutral (zero) protein residual, and both row
  types are ranked together within core_long.groupby(["ensg_id","lineage"]).
  Fixed by treating "no protein arm" as "protein residual = 0" rather than
  "skip arm weighting": core_score = Phi((W_RNA/W_NORM) * rna_std_z) for both
  the protein-free-line and RNA-only-gene fallback paths.
  Verified before promotion (core_score_v7_fallback_fix.py, re-confirmed
  against live inputs unchanged since): n_layers=2 rows bit-identical (max
  diff 0.0); full-panel spike 2.8444% -> 1.4790%; n_layers=1-only spike
  2.7252% -> 0.4894%. Worked example: ENSG00000214309/ach-001396 (n_layers=1)
  0.8858 -> 0.8021, crossing below its (gene, lineage) stratum's best
  n_layers=2 score (0.8025), correcting a ranking inversion that existed
  under the unscaled fallback. See docs/WHY_REFERENCE.md §3.6 for the full
  verification record, and core_score_v4_pre_fallback_fix_*.py.bak for the
  pre-promotion snapshot.

Reads from:  final_pipeline/outputs/{bulk_rna_z, bulk_prot_z, protein_platform_tier}.parquet
             celllineselector.db  (lineage, harmonised_enriched)
Writes to:   final_pipeline/outputs/{core_score, gene_dispersion}.parquet
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

W_RNA  = np.sqrt(1.431)   # Kish n_eff: N=3 sources, ρ_avg=0.548 → n_eff=3/(1+2×0.548)=1.431; replaces provisional sqrt(2.00)
W_PROT = np.sqrt(1.45)    # 2 sources, ProCAN-CCLE rho=0.373 -> n_eff = 2/(1+0.373)
W_NORM = np.sqrt(W_RNA**2 + W_PROT**2)

RHO_PRIOR = 0.353   # median per-gene r(prot_z, rna_z) — RNA-protein EB prior
N_MIN_RHO = 30      # minimum shared lines to estimate per-gene rho
K_RHO     = 30      # EB pseudo-count: prior/data get equal weight at n = K_RHO

# C3 candidate: shrink sd_rna/sd_prot toward their panel-wide values before
# forming beta_g = rho_g*(sd_prot/sd_rna), same w=n/(n+n0) shape used
# throughout this file. n0=25 REUSES the existing convention (RNA_STRATUM_N0/
# PROT_RESID_N0/PROT_STD_STRATUM_N0) deliberately, not independently derived
# for this context -- flagged first-pass/UNVALIDATED, see module docstring.
SD_SHRINK_N0 = 25

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

    # Panel-wide SDs -- computed once, used both as the small-n beta fallback
    # (as before) AND as the C3 shrinkage target inside the loop below.
    _sd_rna_global  = merged["rna_z"].std(ddof=1)
    _sd_prot_global = merged["prot_z"].std(ddof=1)

    rho_rows = []
    for g, grp in merged.groupby("gene_id"):
        if len(grp) < N_MIN_RHO:
            continue
        # rho_raw: per-gene Pearson r(prot_z, rna_z) across shared lines.
        # NOTE: this is NOT the platform-agreement rho (ProCAN vs CCLE, W_PROT=√1.45).
        # Those are distinct quantities; 0.373 applies to W_PROT, 0.353 to this prior.
        rho_raw, _ = stats.pearsonr(grp["prot_z"].fillna(0), grp["rna_z"].fillna(0))
        n   = len(grp)
        lam = K_RHO / (n + K_RHO)   # decreases with n: more data → less shrinkage
        rho_g = lam * RHO_PRIOR + (1 - lam) * rho_raw
        # β = ρ × SD(Y)/SD(X): the OLS regression coefficient, not the correlation.
        # β = ρ only when SD(prot_z) = SD(rna_z); they differ substantially here
        # (measured: SD_rna=1.525, SD_prot=0.951), so using ρ as-is overcorrects by ~1.88×.
        # C3 fix: sd_rna/sd_prot are themselves per-gene sample estimates and can
        # be noisy for genes just above N_MIN_RHO -- shrink each toward its
        # panel-wide value with the same w=n/(n+n0) shape used elsewhere in this
        # file, same as rho_g already is, before forming the ratio.
        sd_rna_raw  = grp["rna_z"].std(ddof=1)
        sd_prot_raw = grp["prot_z"].std(ddof=1)
        w_sd = n / (n + SD_SHRINK_N0)
        sd_rna  = np.sqrt(w_sd * sd_rna_raw**2  + (1 - w_sd) * _sd_rna_global**2)
        sd_prot = np.sqrt(w_sd * sd_prot_raw**2 + (1 - w_sd) * _sd_prot_global**2)
        beta_g  = rho_g * (sd_prot / sd_rna) if sd_rna > 0 else rho_g
        rho_rows.append({"gene_id": g, "rho_g": rho_g, "beta_g": beta_g, "n_shared": n})
    rho_df = pd.DataFrame(rho_rows)

    # Fallback β for genes with too few shared lines (uses global SD ratio)
    _beta_fallback  = (RHO_PRIOR * (_sd_prot_global / _sd_rna_global)
                       if _sd_rna_global > 0 else RHO_PRIOR)

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
        # C7 fix: scale by W_RNA/W_NORM, same as the n_layers=2 combination --
        # "no protein arm" is treated as "protein residual = 0", not as
        # "skip arm weighting". See module docstring.
        fallback_core_z = (W_RNA / W_NORM) * rna_std_z_wide.loc[rna_models_no_prot, shared_genes].values
        fallback_pct = pd.DataFrame(
            stats.norm.cdf(fallback_core_z),
            index=rna_models_no_prot, columns=shared_genes,
        )
        fallback_long = fallback_pct.stack().reset_index()
        fallback_long.columns = ["model_id", "ensg_id", "core_score"]
        fallback_long["lineage"] = lineage_map.reindex(fallback_long["model_id"]).values
        fallback_long["n_layers"] = 1   # RNA-only for these lines (gene has protein, line doesn't)
        fallback_long["combination_variant"] = "A_fallback_scaled"
        core_long = pd.concat([core_long, fallback_long], ignore_index=True)
        print(f"n_layers=1 fallback (protein-free lines, C7-scaled): {len(rna_models_no_prot):,} lines, "
              f"{len(fallback_long):,} pairs")

    # ── RNA-only genes (not in any proteomics platform) ─────────────────────
    rna_only_genes = [g for g in rna_std_z_wide.columns if g not in prot_resid_z_wide.columns]
    if rna_only_genes:
        # C7 fix: same W_RNA/W_NORM scaling as above.
        rna_only_core_z = (W_RNA / W_NORM) * rna_std_z_wide[rna_only_genes].values
        rna_only_pct = pd.DataFrame(
            stats.norm.cdf(rna_only_core_z),
            index=rna_std_z_wide.index, columns=rna_only_genes,
        )
        rna_only_long = rna_only_pct.stack().reset_index()
        rna_only_long.columns = ["model_id", "ensg_id", "core_score"]
        rna_only_long["lineage"] = lineage_map.reindex(rna_only_long["model_id"]).values
        rna_only_long["n_layers"] = 1
        rna_only_long["combination_variant"] = "A_fallback_scaled"
        core_long = pd.concat([core_long, rna_only_long], ignore_index=True)
        print(f"RNA-only genes added (C7-scaled): {len(rna_only_genes):,}")

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

    core_long.to_parquet(CORE_SCORE, index=False)
    disp.to_parquet(GENE_DISP, index=False)

    elapsed = time.time() - t0
    print(f"\nDone in {elapsed:.0f}s")
    print(f"core_score: {len(core_long):,} rows  ->  {CORE_SCORE}")
    print(f"  genes: {core_long.ensg_id.nunique():,}  |  lines: {core_long.model_id.nunique():,}")
    print(f"  core_score mean: {core_long.core_score.mean():.3f}")
    print(f"  combination_variant counts: {core_long.combination_variant.value_counts().to_dict()}")


if __name__ == "__main__":
    run()
