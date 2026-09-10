"""
Scoring/core_score.py
----------------------
Combines pre-scored RNA (3 sources) and protein (2 platforms) into a single
core_score per (gene, cell-line), then ranks within lineage strata.

Architecture:
  RNA z-scores   <- 01_Transcriptomics/rna_scorer.py  -> bulk_rna_z.parquet
  Protein z-scores <- 02_Proteinomics/protein_scorer.py -> bulk_prot_z.parquet
  Tier filter    <- 02_Proteinomics/platform_tier.py  -> protein_platform_tier.parquet
  Metadata       <- Stage-0 warehouse (lineage, model labels)
  -> core_score.parquet, gene_dispersion.parquet

Combination (Cell 7i logic):
  RNA percentile  = rank(rna_z_t, ascending=True) / n  within lineage
  Prot residual   = prot_z_t - beta_g * rna_z_t        (OLS residual; beta_g = rho_g * SD_prot/SD_rna)
  prot_resid_z    = (prot_resid − median) / SD(prot_resid)   per gene  (unit-variance; MAD traded away)
  W_RNA = sqrt(2.0)  [provisional — replace after T3 covariance test]
  W_PROT = sqrt(1.45)  [2 sources, ProCAN-CCLE rho=0.373 -> n_eff=1.45; DISTINCT from RHO_PRIOR=0.353]
  core_z = (W_RNA * rna_pct_z + W_PROT * prot_resid_z) / sqrt(W_RNA^2 + W_PROT^2)
  core_score = Phi(core_z)  [monotone relabelling to [0,1]]

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


def _percentile_within_lineage(z_wide: pd.DataFrame, lineage_map: pd.Series) -> pd.DataFrame:
    z_wide = z_wide.copy()
    z_wide["lineage"] = lineage_map.reindex(z_wide.index)
    pct = z_wide.groupby("lineage", group_keys=False).apply(
        lambda g: g.rank(pct=True, ascending=True, na_option="keep"),
        include_groups=False,
    )
    return pct


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

    # ── RNA: percentile within lineage ─────────────────────────────────────
    rna_wide = rna.pivot(index="model_id", columns="gene_id", values="rna_z")
    rna_pct  = _percentile_within_lineage(rna_wide, lineage_map)
    rna_pct_z_wide = pd.DataFrame(
        stats.norm.ppf(np.clip(rna_pct.values, 0.001, 0.999)),
        index=rna_pct.index, columns=rna_pct.columns,
    )

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
        sd_rna  = grp["rna_z"].std(ddof=1)
        sd_prot = grp["prot_z"].std(ddof=1)
        beta_g  = rho_g * (sd_prot / sd_rna) if sd_rna > 0 else rho_g
        rho_rows.append({"gene_id": g, "rho_g": rho_g, "beta_g": beta_g})
    rho_df = pd.DataFrame(rho_rows)

    # Fallback β for genes with too few shared lines (uses global SD ratio)
    _sd_rna_global  = merged["rna_z"].std(ddof=1)
    _sd_prot_global = merged["prot_z"].std(ddof=1)
    _beta_fallback  = (RHO_PRIOR * (_sd_prot_global / _sd_rna_global)
                       if _sd_rna_global > 0 else RHO_PRIOR)

    merged = merged.merge(rho_df, on="gene_id", how="left")
    merged["rho_g"]  = merged["rho_g"].fillna(RHO_PRIOR)
    merged["beta_g"] = merged["beta_g"].fillna(_beta_fallback)
    merged["prot_resid"] = merged["prot_z"] - merged["beta_g"] * merged["rna_z"]

    prot_resid_wide = merged.pivot(index="model_id", columns="gene_id", values="prot_resid")

    # ── Standardise protein residual: (x − median) / SD per gene ──────────
    # Unit variance is required by the Stouffer weight derivation; both arms
    # must have SD = 1 at the stratum where tier is decided (within-gene).
    # A two-step MAD-z → SD-rescale would cancel algebraically:
    #   (x−med)/(MAD×1.4826) × (MAD×1.4826)/SD(x) = (x−med)/SD(x)
    # We compute the one-step form directly.  Deliberate trade-off: robustness
    # (MAD scale) is yielded in favour of unit variance.  MAD_FLOOR no longer
    # affects the protein path.
    _prot_med = prot_resid_wide.median(axis=0)
    _prot_sd  = prot_resid_wide.std(ddof=1, axis=0).clip(lower=1e-6)
    prot_resid_z_wide = prot_resid_wide.subtract(_prot_med, axis=1).div(_prot_sd, axis=1)
    print(f"  prot_resid_z SD (post-standardise): "
          f"mean={prot_resid_z_wide.std(ddof=1).mean():.4f}  (target 1.00)")

    # ── Orthogonal combination ──────────────────────────────────────────────
    shared_genes   = list(set(rna_pct_z_wide.columns) & set(prot_resid_z_wide.columns))
    shared_models  = list(set(rna_pct_z_wide.index)   & set(prot_resid_z_wide.index))

    rna_m  = rna_pct_z_wide.loc[shared_models, shared_genes].values
    prot_m = prot_resid_z_wide.loc[shared_models, shared_genes].values

    core_z = (W_RNA * rna_m + W_PROT * prot_m) / W_NORM
    core_pct = pd.DataFrame(
        stats.norm.cdf(core_z),
        index=shared_models, columns=shared_genes,
    )

    # ── Combined genes (RNA + Protein) — lines with protein data ───────────
    core_long = core_pct.stack().reset_index()
    core_long.columns = ["model_id", "ensg_id", "core_score"]
    core_long["lineage"] = lineage_map.reindex(core_long["model_id"]).values
    core_long["n_layers"] = 2

    # ── n_layers=1 FALLBACK for shared genes at protein-free lines ──────────
    # Without this, 48.6% of RNA pairs (lines not in shared_models) get no
    # score for the 10,871 shared-gene set. Per-pair fallback to RNA-only
    # for those lines, so every RNA-profiled line gets a score for every gene.
    shared_models_set = set(shared_models)
    rna_models_no_prot = [m for m in rna_pct_z_wide.index if m not in shared_models_set]
    if rna_models_no_prot and shared_genes:
        fallback_pct = pd.DataFrame(
            stats.norm.cdf(rna_pct_z_wide.loc[rna_models_no_prot, shared_genes].values),
            index=rna_models_no_prot, columns=shared_genes,
        )
        fallback_long = fallback_pct.stack().reset_index()
        fallback_long.columns = ["model_id", "ensg_id", "core_score"]
        fallback_long["lineage"] = lineage_map.reindex(fallback_long["model_id"]).values
        fallback_long["n_layers"] = 1   # RNA-only for these lines (gene has protein, line doesn't)
        core_long = pd.concat([core_long, fallback_long], ignore_index=True)
        print(f"n_layers=1 fallback (protein-free lines): {len(rna_models_no_prot):,} lines, "
              f"{len(fallback_long):,} pairs")

    # ── RNA-only genes (not in any proteomics platform) ─────────────────────
    rna_only_genes = [g for g in rna_pct_z_wide.columns if g not in prot_resid_z_wide.columns]
    if rna_only_genes:
        rna_only_pct = pd.DataFrame(
            stats.norm.cdf(rna_pct_z_wide[rna_only_genes].values),
            index=rna_pct_z_wide.index, columns=rna_only_genes,
        )
        rna_only_long = rna_only_pct.stack().reset_index()
        rna_only_long.columns = ["model_id", "ensg_id", "core_score"]
        rna_only_long["lineage"] = lineage_map.reindex(rna_only_long["model_id"]).values
        rna_only_long["n_layers"] = 1
        core_long = pd.concat([core_long, rna_only_long], ignore_index=True)
        print(f"RNA-only genes added: {len(rna_only_genes):,}")

    # ── Stratum rank (within lineage) ──────────────────────────────────────
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


if __name__ == "__main__":
    run()
