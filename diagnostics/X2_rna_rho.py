"""
X2 — Measure rho_rna between RNA sources; assess W_RNA.
COWORK v4.0 spec: compute-only, no changes to W_RNA.
Run from repo root: python diagnostics/X2_rna_rho.py
"""
import sys, time
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import stats

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "final_pipeline"))
from config import DB, RNA_Z
sys.path.insert(0, str(REPO / "final_pipeline" / "utils"))
import common as C

N_GENE_SAMPLE = 2000
SEED          = 42
MIN_SHARED    = 30

print("=" * 72)
print("X2 -- rho_rna measurement (COWORK v4.0)")
print("=" * 72)

t0 = time.time()
con = C.connect()
lineage_map = C.load_lineage(con)

# Gene universe: intersection of all 3 sources
gene_univ = C.gene_universe(con)
all_genes = gene_univ["gene_id"].tolist()
print(f"[COUNT] Gene universe (all 3 sources): {len(all_genes):,}")

rng = np.random.default_rng(SEED)
if len(all_genes) > N_GENE_SAMPLE:
    sample_genes = list(rng.choice(all_genes, N_GENE_SAMPLE, replace=False))
else:
    sample_genes = all_genes
print(f"[COUNT] Sampled genes for rho computation: {len(sample_genes):,}  (seed={SEED})")

# Fetch per-source data for sampled genes
print("[STATUS] Fetching per-source data from DuckDB...")
dm = C.fetch_depmap(con, sample_genes)
hp = C.fetch_hpa(con, sample_genes)
ge = C.fetch_geo(con, sample_genes)
con.close()

print(f"[COUNT] DepMap rows: {len(dm):,}  genes: {dm.gene_id.nunique():,}  models: {dm.model_id.nunique():,}")
print(f"[COUNT] HPA rows:    {len(hp):,}  genes: {hp.gene_id.nunique():,}  models: {hp.model_id.nunique():,}")
print(f"[COUNT] GEO rows:    {len(ge):,}  genes: {ge.gene_id.nunique():,}  models: {ge.model_id.nunique():,}")

# Compute lineage-conditioned z-scores per source
print("[STATUS] Computing z-scores per source...")
z_dm = C.score_source_lineage(
    dm.pivot_table(index="model_id", columns="gene_id", values="value", aggfunc="first"),
    sample_genes, lineage_map, "depmap"
)
z_hp = C.score_source_lineage(
    hp.pivot_table(index="model_id", columns="gene_id", values="value", aggfunc="first"),
    sample_genes, lineage_map, "hpa"
)
z_ge = C.score_source_lineage(
    ge.pivot_table(index="model_id", columns="gene_id", values="value", aggfunc="first"),
    sample_genes, lineage_map, "geo"
)
print(f"[COUNT] z_depmap: {len(z_dm):,}  z_hpa: {len(z_hp):,}  z_geo: {len(z_ge):,}")

# Per-gene Spearman correlation for each pair
def per_gene_spearman(za, zb, name_a, name_b):
    joined = za.rename(columns={"z": f"z_{name_a}"})[["gene_id","model_id",f"z_{name_a}"]].merge(
        zb.rename(columns={"z": f"z_{name_b}"})[["gene_id","model_id",f"z_{name_b}"]],
        on=["gene_id","model_id"], how="inner"
    ).dropna()
    print(f"[COUNT] {name_a} x {name_b}: {len(joined):,} shared (gene,model) pairs  "
          f"genes: {joined.gene_id.nunique():,}")
    rhos = []
    for gid, grp in joined.groupby("gene_id"):
        if len(grp) < MIN_SHARED:
            continue
        r, _ = stats.spearmanr(grp[f"z_{name_a}"], grp[f"z_{name_b}"])
        if np.isfinite(r):
            rhos.append(r)
    rhos = np.array(rhos)
    if len(rhos) == 0:
        print(f"[RESULT] {name_a} x {name_b}: NO eligible genes")
        return np.nan
    print(f"[RESULT] {name_a} x {name_b}: n_genes={len(rhos):,}")
    print(f"  median={np.median(rhos):.4f}  mean={np.mean(rhos):.4f}  SD={np.std(rhos):.4f}")
    print(f"  IQR=[{np.percentile(rhos,25):.4f}, {np.percentile(rhos,75):.4f}]")
    print(f"  p5={np.percentile(rhos,5):.4f}  p95={np.percentile(rhos,95):.4f}")
    print(f"  >0.3: {(rhos>0.3).mean()*100:.1f}%  >0.5: {(rhos>0.5).mean()*100:.1f}%")
    return float(np.median(rhos))

print("\n[STATUS] Computing pairwise Spearman rho_rna...")
rho_dm_hp = per_gene_spearman(z_dm, z_hp, "depmap", "hpa")
rho_dm_ge = per_gene_spearman(z_dm, z_ge, "depmap", "geo")
rho_hp_ge = per_gene_spearman(z_hp, z_ge, "hpa", "geo")

valid_rhos = [r for r in [rho_dm_hp, rho_dm_ge, rho_hp_ge] if np.isfinite(r)]
if valid_rhos:
    rho_avg = np.mean(valid_rhos)
    print(f"\n[RESULT] Average pairwise rho_rna: {rho_avg:.4f}  "
          f"(mean of {len(valid_rhos)} pairwise medians)")
else:
    rho_avg = np.nan
    print("[RESULT] Could not compute rho_rna")

print(f"\n--- W_RNA ASSESSMENT (N_SOURCES=3) ---")
print(f"[CURRENT] W_RNA = sqrt(2.00) = {np.sqrt(2.00):.6f}")
print(f"[NOTE] Code comment says 'n_eff=2 RNA sources' but N_SOURCES=3 in rna_scorer.py line 29")
print(f"[NOTE] Under independence (rho=0), W_RNA should be sqrt(3) = {np.sqrt(3):.6f}")
print(f"[NOTE] Current W_RNA = sqrt(2) implies n_eff=2, which requires rho_avg such that")
print(f"       n_eff = N / (1 + (N-1)*rho_avg) = 3 / (1+2*rho_avg) = 2")
print(f"       => rho_implied = 0.25 (not measured, just assumed)")

if np.isfinite(rho_avg):
    n_eff_corrected = 3.0 / (1.0 + 2.0 * rho_avg)
    W_RNA_corrected = np.sqrt(n_eff_corrected)
    pct_change = 100 * (W_RNA_corrected - np.sqrt(2.00)) / np.sqrt(2.00)
    print(f"\n[MEASURED] rho_avg = {rho_avg:.4f}")
    print(f"[MEASURED] n_eff_corrected = 3/(1+2*{rho_avg:.4f}) = {n_eff_corrected:.4f}")
    print(f"[MEASURED] W_RNA_corrected = sqrt({n_eff_corrected:.4f}) = {W_RNA_corrected:.6f}")
    print(f"[MEASURED] % change from current sqrt(2.00) = {pct_change:+.2f}%")

    # Variance shares
    W_RNA_curr  = np.sqrt(2.00)
    W_PROT      = np.sqrt(1.45)
    W_NORM_curr = np.sqrt(W_RNA_curr**2 + W_PROT**2)
    W_NORM_corr = np.sqrt(W_RNA_corrected**2 + W_PROT**2)

    share_rna_curr  = 100 * W_RNA_curr**2 / W_NORM_curr**2
    share_prot_curr = 100 * W_PROT**2     / W_NORM_curr**2
    share_rna_corr  = 100 * W_RNA_corrected**2 / W_NORM_corr**2
    share_prot_corr = 100 * W_PROT**2          / W_NORM_corr**2

    print(f"\n[MEASURED] Current variance shares: RNA {share_rna_curr:.1f}%  Protein {share_prot_curr:.1f}%")
    print(f"[MEASURED] Corrected variance shares: RNA {share_rna_corr:.1f}%  Protein {share_prot_corr:.1f}%")

    invert = W_RNA_corrected < W_PROT
    print(f"\n[MEASURED] Arms inverted (RNA < PROT)? {invert}  "
          f"(RNA {W_RNA_corrected:.4f} vs PROT {W_PROT:.4f})")

print(f"\n[ELAPSED] {time.time()-t0:.0f}s")
print("\n" + "=" * 72)
print("ESCALATION REQUIRED — see COWORK v4.0 X2 escalation block")
print("Options: (A) Correct and re-run  (B) Document and defer")
print("Decision: Fiona")
print("=" * 72)
