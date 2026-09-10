"""
expression_fusion/test_run_t2_protein_residual.py
--------------------------------------------------
T2: Protein residual reproducibility  — PRE-REGISTRATION GATE TEST

Pre-registration: docs/EXPR_FUSION_PREREGISTRATION.md §4, T2
RNG seeds: BASE_SEED + replicate_index  (BASE_SEED = 20260813, per §3)
N1 replicates: B = 200
Bootstrap resamples: 10 000  (over genes, per §3)

QUESTION
--------
Is prot_resid = prot_z - rho_g * rna_z post-transcriptional biology,
or platform noise?

If ProCAN and CCLE capture the same post-transcriptional biology, their
residuals (after removing RNA-explained variance) should correlate across
the ~288 shared cell lines, gene by gene.

PROTOCOL (verbatim from pre-registration)
-----------------------------------------
For each gene g with both platforms and RNA:
  rho_g^P  = shrunk Pearson(rna_z, procan_z)   on shared lines
  resid_P  = procan_z - rho_g^P * rna_z
  rho_g^C  = shrunk Pearson(rna_z, ccle_z)     on shared lines
  resid_C  = ccle_z  - rho_g^C  * rna_z
  rho_resid[g] = Spearman(resid_P, resid_C)    across shared lines

Direction: ACROSS LINES within gene (not across genes within line).

PRIMARY STATISTIC
-----------------
rho_tilde = median rho_resid[g] for genes with n_overlap >= 30.
Bootstrap CI: 10 000 resamples over genes.
N1 null: 200 replicates, both platforms permuted within gene independently.

PASS: rho_tilde >= 0.20 AND lower CI > 0.10 AND rho_tilde > N1 97.5th pctile
FAIL: rho_tilde < 0.10 OR CI includes null
INDETERMINATE: 0.10 <= rho_tilde < 0.20 with CI excluding null -> treated as FAIL

SHRINKAGE PARAMETERS (from cell 7i, unchanged)
----------------------------------------------
RHO_PRIOR = 0.373  (median ProCAN-CCLE Spearman from overlap test)
N_MIN_RHO = 30     (shrink toward prior below this n_overlap)
"""
from __future__ import annotations

import sys, time, json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import common as C

# import protein loaders from scorer (avoids code duplication)
sys.path.insert(0, str(HERE))
import importlib.util as _ilu
_spec = _ilu.spec_from_file_location("prot_scorer", HERE / "06_bulk_protein_scorer.py")
_prot_mod = _ilu.module_from_spec(_spec)
_spec.loader.exec_module(_prot_mod)
_build_uniprot_ensg = _prot_mod._build_uniprot_ensg
_load_procan        = _prot_mod._load_procan
_load_ccle          = _prot_mod._load_ccle

OUT    = HERE / "outputs"
RESULT = OUT / "t2_protein_residual_result.json"

BASE_SEED  = 20260813
B          = 200
N_BOOT     = 10_000
N_MIN_OVL  = 30       # pre-registration: genes with n_overlap >= 30
RHO_PRIOR  = 0.373
N_MIN_RHO  = 30
PASS_MED   = 0.20
PASS_CI_LO = 0.10


# ── per-gene Pearson (vectorised, NaN-safe) ──────────────────────────────────
def _pearson_per_gene(X: np.ndarray, Y: np.ndarray):
    """X, Y: (n_lines, n_genes). Returns (rho, n_valid) per gene."""
    valid = np.isfinite(X) & np.isfinite(Y)
    n     = valid.sum(axis=0)
    Xm = np.where(valid, X, np.nan)
    Ym = np.where(valid, Y, np.nan)
    mx = np.nanmean(Xm, axis=0)
    my = np.nanmean(Ym, axis=0)
    Xc = Xm - mx;  Yc = Ym - my
    num = np.nansum(Xc * Yc, axis=0)
    dX  = np.sqrt(np.nansum(Xc ** 2, axis=0))
    dY  = np.sqrt(np.nansum(Yc ** 2, axis=0))
    den = dX * dY
    rho = np.where(den > 0, num / np.where(den > 0, den, 1.0), np.nan)
    rho = np.where(n >= 3, rho, np.nan)
    return rho, n


# ── per-gene Spearman (rank then Pearson) ────────────────────────────────────
def _rank_matrix(X: np.ndarray) -> np.ndarray:
    """Rank each column independently, preserving NaN positions."""
    from scipy.stats import rankdata
    R = np.full_like(X, np.nan)
    for g in range(X.shape[1]):
        col   = X[:, g]
        valid = np.isfinite(col)
        if valid.sum() >= 2:
            R[valid, g] = rankdata(col[valid])
    return R


def _spearman_per_gene(resid_P: np.ndarray, resid_C: np.ndarray) -> np.ndarray:
    """Spearman rho[g] = Spearman(resid_P[:, g], resid_C[:, g])."""
    rP = _rank_matrix(resid_P)
    rC = _rank_matrix(resid_C)
    rho, _ = _pearson_per_gene(rP, rC)
    return rho


# ── T2 core computation ───────────────────────────────────────────────────────
def _compute_t2(rna_M: np.ndarray,
                prot_P: np.ndarray,
                prot_C: np.ndarray,
                n_valid_P: np.ndarray,
                n_valid_C: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """
    All arrays: (n_shared_lines, n_genes).
    Returns rho_resid (n_genes,) and n_overlap (n_genes,).
    """
    # Per-platform rho with James-Stein shrinkage toward RHO_PRIOR
    def _shrunk_rho(prot_M, n_valid):
        rho_raw, _ = _pearson_per_gene(rna_M, prot_M)
        lam     = np.where(n_valid >= N_MIN_RHO,
                           N_MIN_RHO / (n_valid + N_MIN_RHO), 1.0)
        rho_raw = np.where(np.isfinite(rho_raw), rho_raw, RHO_PRIOR)
        return np.clip(lam * RHO_PRIOR + (1 - lam) * rho_raw, -0.99, 0.99)

    rho_P = _shrunk_rho(prot_P, n_valid_P)
    rho_C = _shrunk_rho(prot_C, n_valid_C)

    # Residuals
    resid_P = prot_P - rho_P[np.newaxis, :] * rna_M
    resid_C = prot_C - rho_C[np.newaxis, :] * rna_M

    # Propagate NaN where either platform or RNA is missing
    both = np.isfinite(prot_P) & np.isfinite(prot_C) & np.isfinite(rna_M)
    resid_P = np.where(both, resid_P, np.nan)
    resid_C = np.where(both, resid_C, np.nan)

    n_overlap = both.sum(axis=0)
    rho_resid = _spearman_per_gene(resid_P, resid_C)
    return rho_resid, n_overlap


def main():
    t0 = time.time()
    print("=" * 72)
    print("T2 -- Protein residual reproducibility  (B=200 N1, 10k bootstrap)")
    print("Pre-registration: docs/EXPR_FUSION_PREREGISTRATION.md section 4")
    print("=" * 72)

    # ── 1. Load protein source-level data from DB ─────────────────────────
    print("\n[1] Loading protein data from DB ...", flush=True)
    con         = C.connect()
    lineage_map = C.load_lineage(con)
    u2e         = _build_uniprot_ensg(con)

    procan_wide = _load_procan(con, u2e)   # model_id (index) x ensg_id (cols)
    ccle_wide   = _load_ccle(con, u2e)
    con.close()

    # Identify shared model_ids (ProCAN ∩ CCLE)
    shared_models = sorted(set(procan_wide.index) & set(ccle_wide.index))
    print(f"  ProCAN: {procan_wide.shape[0]} lines  "
          f"CCLE: {ccle_wide.shape[0]} lines  "
          f"Shared: {len(shared_models)}")

    # Shared genes (ProCAN ∩ CCLE)
    shared_genes_prot = sorted(set(procan_wide.columns) & set(ccle_wide.columns))
    print(f"  Shared genes (ProCAN+CCLE): {len(shared_genes_prot)}")

    # ── 2. Lineage-conditioned z-scores per platform ─────────────────────
    print("\n[2] Computing lineage z-scores ...", flush=True)
    procan_z = C.score_source_lineage(
        procan_wide, shared_genes_prot, lineage_map, "procan")
    ccle_z = C.score_source_lineage(
        ccle_wide, shared_genes_prot, lineage_map, "ccle_prot")
    print(f"  ProCAN z: {len(procan_z):,} rows  "
          f"CCLE z: {len(ccle_z):,} rows  ({time.time()-t0:.0f}s)")

    # ── 3. Load RNA z-scores ──────────────────────────────────────────────
    print("\n[3] Loading RNA z-scores from bulk_rna_z.parquet ...", flush=True)
    rna_long = pd.read_parquet(OUT / "bulk_rna_z.parquet")
    # Keep only RNA genes that are also in prot shared set
    rna_long = rna_long[rna_long.gene_id.isin(shared_genes_prot)]
    # Keep only shared model_ids
    rna_long = rna_long[rna_long.model_id.isin(shared_models)]
    print(f"  RNA rows (shared lines+genes): {len(rna_long):,}")

    # ── 4. Build aligned matrices ─────────────────────────────────────────
    print("\n[4] Pivoting to matrices (shared_lines x genes) ...", flush=True)

    def _to_wide_matrix(long_df, val_col, model_ids, genes):
        w = (long_df.pivot(index="model_id", columns="gene_id", values=val_col)
                    .reindex(index=model_ids, columns=genes))
        return w.values.astype(float)  # (n_lines, n_genes)

    procan_z_wide = (procan_z[procan_z.model_id.isin(shared_models)]
                     .pivot(index="model_id", columns="gene_id", values="z")
                     .reindex(index=shared_models, columns=shared_genes_prot)
                     .values.astype(float))
    ccle_z_wide   = (ccle_z[ccle_z.model_id.isin(shared_models)]
                     .pivot(index="model_id", columns="gene_id", values="z")
                     .reindex(index=shared_models, columns=shared_genes_prot)
                     .values.astype(float))
    rna_z_wide    = _to_wide_matrix(rna_long, "rna_z_t",
                                    shared_models, shared_genes_prot)

    genes_arr  = np.array(shared_genes_prot)
    n_lines    = len(shared_models)
    n_genes    = len(shared_genes_prot)

    # NaN coverage
    n_valid_P = (np.isfinite(procan_z_wide) & np.isfinite(rna_z_wide)).sum(axis=0)
    n_valid_C = (np.isfinite(ccle_z_wide)   & np.isfinite(rna_z_wide)).sum(axis=0)
    all_valid  = np.isfinite(procan_z_wide) & np.isfinite(ccle_z_wide) & np.isfinite(rna_z_wide)
    n_overlap  = all_valid.sum(axis=0)

    print(f"  Matrix shape: {n_lines} lines x {n_genes} genes")
    print(f"  Genes with n_overlap >= {N_MIN_OVL}: "
          f"{(n_overlap >= N_MIN_OVL).sum():,}")
    print(f"  Median n_overlap: {np.median(n_overlap[n_overlap>0]):.0f}  "
          f"({time.time()-t0:.0f}s)")

    # ── 5. Real-data T2 statistic ─────────────────────────────────────────
    print("\n[5] Computing real-data rho_resid per gene ...", flush=True)
    rho_resid_real, n_ovl_real = _compute_t2(
        rna_z_wide, procan_z_wide, ccle_z_wide, n_valid_P, n_valid_C)

    eligible    = (n_ovl_real >= N_MIN_OVL) & np.isfinite(rho_resid_real)
    rho_el_real = rho_resid_real[eligible]
    rho_tilde   = float(np.nanmedian(rho_el_real))
    print(f"  Eligible genes (n_overlap>={N_MIN_OVL}, finite rho): {eligible.sum():,}")
    print(f"  rho_tilde (real) = {rho_tilde:.4f}")
    print(f"  raw protein-protein rho (ceiling) = 0.373")
    print(f"  Distribution: min={rho_el_real.min():.3f}  "
          f"p25={np.percentile(rho_el_real,25):.3f}  "
          f"median={np.median(rho_el_real):.3f}  "
          f"p75={np.percentile(rho_el_real,75):.3f}  "
          f"max={rho_el_real.max():.3f}")

    # ── 6. Bootstrap CI over genes (10 000 resamples) ────────────────────
    print("\n[6] Bootstrap CI (10 000 resamples over genes) ...", flush=True)
    rng_boot = np.random.default_rng(BASE_SEED)
    boot_med = np.empty(N_BOOT)
    idx_pool = np.where(eligible)[0]
    for b in range(N_BOOT):
        idx = rng_boot.choice(idx_pool, size=len(idx_pool), replace=True)
        boot_med[b] = np.nanmedian(rho_resid_real[idx])
    ci_lo = float(np.percentile(boot_med,  2.5))
    ci_hi = float(np.percentile(boot_med, 97.5))
    print(f"  95% bootstrap CI = [{ci_lo:.4f}, {ci_hi:.4f}]")

    # ── 7. N1 null (B = 200 replicates) ──────────────────────────────────
    print(f"\n[7] N1 null ({B} replicates) ...", flush=True)
    # Permute ProCAN and CCLE independently within each gene across shared lines
    n1_medians = np.empty(B)
    for rep in range(B):
        rng_rep = np.random.default_rng(BASE_SEED + rep + 1)  # +1 to separate from bootstrap seed

        P_perm = procan_z_wide.copy()
        C_perm = ccle_z_wide.copy()
        for g in range(n_genes):
            # permute ProCAN across lines where valid
            vP = np.where(np.isfinite(P_perm[:, g]))[0]
            if len(vP) > 1:
                P_perm[vP, g] = rng_rep.permutation(P_perm[vP, g])
            # permute CCLE across lines where valid (independent)
            vC = np.where(np.isfinite(C_perm[:, g]))[0]
            if len(vC) > 1:
                C_perm[vC, g] = rng_rep.permutation(C_perm[vC, g])

        rho_null, _ = _compute_t2(
            rna_z_wide, P_perm, C_perm, n_valid_P, n_valid_C)
        eligible_null = (n_ovl_real >= N_MIN_OVL) & np.isfinite(rho_null)
        n1_medians[rep] = float(np.nanmedian(rho_null[eligible_null]))

        if (rep + 1) % 20 == 0:
            print(f"  rep {rep+1:3d}/{B}  N1 median={np.nanmedian(n1_medians[:rep+1]):.4f}"
                  f"  ({time.time()-t0:.0f}s)", flush=True)

    n1_p975 = float(np.nanpercentile(n1_medians, 97.5))
    n1_med  = float(np.nanmedian(n1_medians))
    print(f"  N1 median = {n1_med:.4f}  97.5th pctile = {n1_p975:.4f}")

    # ── 8. Verdict ────────────────────────────────────────────────────────
    above_n1      = rho_tilde > n1_p975
    ci_above_null = ci_lo > n1_p975

    if rho_tilde >= PASS_MED and ci_lo > PASS_CI_LO and above_n1:
        verdict = "PASS"
    elif rho_tilde < 0.10 or (ci_lo <= n1_p975):
        verdict = "FAIL"
    else:
        verdict = "INDETERMINATE"

    k2 = verdict in ("FAIL", "INDETERMINATE")

    print(f"\n{'='*72}")
    print(f"rho_tilde  = {rho_tilde:.4f}  (PASS threshold: >= {PASS_MED})")
    print(f"CI 95%     = [{ci_lo:.4f}, {ci_hi:.4f}]  (PASS: lower > {PASS_CI_LO})")
    print(f"N1 97.5pct = {n1_p975:.4f}  rho_tilde > N1: {above_n1}")
    print(f"\nT2 VERDICT: {verdict}")
    if k2:
        print("KILL SWITCH K2 TRIGGERED.")
        print("Protein arm removed permanently for this thesis.")
        print("Grounds: reliability ceiling (sqrt(0.373) = 0.611) and MNAR missingness,")
        print("with T2 as direct measurement. Orthogonalisation step withdrawn.")
    print(f"Total elapsed: {time.time()-t0:.0f}s")
    print("=" * 72)

    output = {
        "test": "T2",
        "preregistration": "docs/EXPR_FUSION_PREREGISTRATION.md section 4",
        "B": B, "n_boot": N_BOOT, "base_seed": BASE_SEED,
        "n_shared_lines": len(shared_models),
        "n_shared_genes_prot": n_genes,
        "n_eligible_genes": int(eligible.sum()),
        "overall_verdict": verdict,
        "kill_switch_K2": k2,
        "rho_tilde_real": rho_tilde,
        "bootstrap_ci_95": [ci_lo, ci_hi],
        "n1_median": n1_med,
        "n1_p975": n1_p975,
        "rho_tilde_above_n1": bool(above_n1),
        "ci_lo_above_n1_p975": bool(ci_above_null),
        "pass_thresholds": {
            "rho_tilde_min": PASS_MED,
            "ci_lo_min": PASS_CI_LO,
        },
        "raw_platform_rho_ceiling": 0.373,
    }
    RESULT.write_text(json.dumps(output, indent=2), encoding="utf-8")
    print(f"\nResult saved -> {RESULT}")


if __name__ == "__main__":
    main()
