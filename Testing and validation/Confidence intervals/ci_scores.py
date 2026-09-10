"""
Scoring/ci_scores.py
--------------------
Batch confidence intervals + rank-stability for core_score, written as a parquet
that fits alongside the other pipeline outputs.

For every (gene, cell-line) it emits:
  core_score            the stored point score (read verbatim from core_score.parquet)
  ci_lo, ci_hi, ci_width   95% interval from a Monte-Carlo of the input noise
  topk_global_pct       % of draws the line is in the global top-k
  topk_lineage_pct      % of draws the line is in its tissue's top-k
  top1_winrate          % of draws the line is the single global best
  n_layers, n_sources   coverage flags

Design notes
  * The POINT score is read from core_score.parquet, never recomputed — so this
    column matches the pipeline exactly by construction, with zero drift.
  * Only the DRAWS are recomputed (vectorised, all N at once per gene).
  * Noise sizes come from the data: RNA between-source disagreement + coverage,
    protein platform-agreement tier. Same knobs as the per-query notebook.

Inputs (all already produced by earlier stages):
  bulk_rna_z.parquet, bulk_prot_z.parquet, protein_platform_tier.parquet,
  core_score.parquet, prot_resid_stats.parquet, and (built here once)
  rna_dispersion.parquet.
Writes: architecture/outputs/ci_scores.parquet
"""
import sys, time
from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
from scipy import stats

# for running in path Scoring/ use the commented code below so file can read config.py and utils/common.py
# sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
# from config import RNA_Z, PROT_Z, PROT_TIER, CORE_SCORE
# sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "utils"))
# import common as C

PIPELINE_ROOT = Path(__file__).resolve().parent.parent.parent / "architecture"  # <-- the folder with config.py
for p in (PIPELINE_ROOT, PIPELINE_ROOT / "utils", PIPELINE_ROOT / "02_Proteinomics"):
    sp = str(p)
    if sp not in sys.path:
        sys.path.insert(0, sp)

from config import RNA_Z, PROT_Z, PROT_TIER, CORE_SCORE
import common as C

OUT_DIR   = CORE_SCORE.parent
CI_SCORES = OUT_DIR / "ci_scores.parquet"
RNA_DISP  = OUT_DIR / "rna_dispersion.parquet"
RESID_ST  = OUT_DIR / "prot_resid_stats.parquet"

# ── constants (must match core_score.py) ─────────────────────────────────────
W_RNA  = np.sqrt(1.431)
W_PROT = np.sqrt(1.45)
W_NORM = np.sqrt(W_RNA**2 + W_PROT**2)
RNA_SOURCES_TOTAL = 3

# ── noise knobs (same as the notebook) ───────────────────────────────────────
RNA_BASE_FLOOR = 0.30
RNA_TAU_PRIOR  = 0.25
COVERAGE_CAP   = 1.5
PROT_BASE      = 0.30
PROT_TIER_MULT = {"consistent": 1.0, "cautious": 1.5, "conflicting": 2.2}
PROT_TAU_PRIOR = 0.45

# ── run parameters ───────────────────────────────────────────────────────────
N_DRAWS = 500     # draws per gene; 500 is plenty for stable percentiles, halves runtime vs 1000
TOPK    = 10
SEED    = 42


# ═════════════════════════════════════════════════════════════════════════════
# One-time: RNA between-source dispersion (the only input not already persisted)
# ═════════════════════════════════════════════════════════════════════════════
def build_rna_dispersion(con):
    """Per (gene, line): SD of the per-source RNA z's — the data-driven noise.
    Heavy but run ONCE; cached to rna_dispersion.parquet. Mirrors rna_scorer's
    per-source step, kept before the Stouffer combine."""
    print("Building rna_dispersion.parquet (one-time) …")
    lineage_map = C.load_lineage(con)
    gene_ids = list(pd.read_parquet(RNA_Z, columns=["gene_id"])["gene_id"].str.lower().unique())

    def _wide(long):
        return long.pivot_table(index="model_id", columns="gene_id", values="value", aggfunc="median")

    frames = []
    for fetch, name in [(C.fetch_depmap, "depmap_expr"),
                        (C.fetch_geo,    "geo_expr"),
                        (C.fetch_hpa,    "hpa_rna")]:
        long = fetch(con, gene_ids)
        if len(long):
            frames.append(C.score_source_lineage(_wide(long), gene_ids, lineage_map, name))
    per = pd.concat(frames, ignore_index=True)[["gene_id", "model_id", "z"]].dropna()

    disp = (per.groupby(["gene_id", "model_id"])["z"]
               .agg(s_between="std", k_obs="count").reset_index())
    disp["gene_id"] = disp["gene_id"].str.upper()
    disp["model_id"] = disp["model_id"].str.lower()
    disp.to_parquet(RNA_DISP, index=False)
    print(f"  wrote {RNA_DISP}  ({len(disp):,} rows)")
    return disp


# ═════════════════════════════════════════════════════════════════════════════
# Load shared inputs once, indexed for fast per-gene slicing
# ═════════════════════════════════════════════════════════════════════════════
def load_inputs(con):
    lineage_map = C.load_lineage(con); lineage_map.index = lineage_map.index.str.lower()

    rna = pd.read_parquet(RNA_Z).rename(columns={"z_t": "rna_z"})
    rna["gene_id"] = rna["gene_id"].str.upper(); rna["model_id"] = rna["model_id"].str.lower()

    prot = pd.read_parquet(PROT_Z).rename(columns={"z_t": "prot_z"})
    prot["gene_id"] = prot["gene_id"].str.upper(); prot["model_id"] = prot["model_id"].str.lower()

    core = pd.read_parquet(CORE_SCORE)[["ensg_id", "model_id", "core_score"]]
    core["ensg_id"] = core["ensg_id"].str.upper(); core["model_id"] = core["model_id"].str.lower()

    stats_df = pd.read_parquet(RESID_ST); stats_df["ensg_id"] = stats_df["ensg_id"].str.upper()
    stats_df = stats_df.set_index("ensg_id")

    disp = pd.read_parquet(RNA_DISP) if RNA_DISP.exists() else build_rna_dispersion(con)

    # protein tier -> per-gene multiplier (via uniprot map)
    import protein_scorer as PS
    u2e = PS._build_uniprot_ensg(con)
    e2u = {}
    for u, e in u2e.items():
        e2u.setdefault(e, set()).add(u)
    tier_df = pd.read_parquet(PROT_TIER); tier_df["uniprot"] = tier_df["uniprot"].str.lower()
    tier_lookup = dict(zip(tier_df["uniprot"], tier_df["platform_tier"]))
    rank = {"conflicting": 0, "cautious": 1, "consistent": 2}
    gene_tier = {}
    for e, us in e2u.items():
        tiers = [tier_lookup[u] for u in us if u in tier_lookup]
        gene_tier[e] = min(tiers, key=lambda t: rank[t]) if tiers else "cautious"

    # merge RNA noise once
    rna = rna.merge(disp[["gene_id", "model_id", "s_between", "k_obs"]],
                    left_on=["gene_id", "model_id"], right_on=["gene_id", "model_id"], how="left")
    rna["lineage"] = rna["model_id"].map(lineage_map)
    return rna, prot, core, stats_df, gene_tier


def _rna_se(k, s_between):
    se_between = np.where(k >= 2, s_between / np.sqrt(np.maximum(k, 1)), RNA_TAU_PRIOR)
    se_between = np.where(np.isfinite(se_between), se_between, RNA_TAU_PRIOR)
    se_within  = RNA_BASE_FLOOR * np.minimum(np.sqrt(RNA_SOURCES_TOTAL / np.maximum(k, 1)), COVERAGE_CAP)
    return np.sqrt(se_between**2 + se_within**2)


# ═════════════════════════════════════════════════════════════════════════════
# Vectorised Monte-Carlo for ONE gene (validated against the per-gene reference)
# ═════════════════════════════════════════════════════════════════════════════
def gene_ci(g_rna, g_prot, point, beta, med, sd, N=N_DRAWS, k=TOPK, seed=SEED, ci=0.95):
    base = g_rna.merge(g_prot, on="model_id", how="left").reset_index(drop=True)
    n = len(base)
    if n == 0:
        return None
    ids = base["model_id"].to_numpy(); lin = base["lineage"].to_numpy()
    rna_mu = base["rna_z"].to_numpy()
    rna_sd = _rna_se(base["n_sources"].to_numpy(), base["s_between"].to_numpy())
    has_p = base["prot_z"].notna().to_numpy()
    prot_mu = base["prot_z"].to_numpy(float)
    prot_sd = np.where(np.isfinite(base["se_prot"].to_numpy(float)), base["se_prot"].to_numpy(float), 0.0)

    rng = np.random.default_rng(seed)
    R = rna_mu[None, :] + rng.normal(size=(N, n)) * rna_sd[None, :]
    P = np.full((N, n), np.nan)
    if has_p.any():
        P[:, has_p] = prot_mu[has_p][None, :] + rng.normal(size=(N, int(has_p.sum()))) * prot_sd[has_p][None, :]

    pct = np.empty((N, n))
    codes, uniq = pd.factorize(lin)
    groups = [np.where(codes == c)[0] for c in range(len(uniq))]
    for idx in groups:
        sub = R[:, idx]
        pct[:, idx] = (sub.argsort(axis=1).argsort(axis=1) + 1) / len(idx)
    rna_pct_z = stats.norm.ppf(np.clip(pct, 0.001, 0.999))

    resid_z = np.full((N, n), np.nan)
    if has_p.any():
        resid_z[:, has_p] = (P[:, has_p] - beta * R[:, has_p] - med) / sd
    core_z2 = (W_RNA * rna_pct_z + W_PROT * resid_z) / W_NORM
    core = np.where(has_p[None, :], stats.norm.cdf(core_z2), stats.norm.cdf(rna_pct_z))

    top_global = np.zeros(n); top_lineage = np.zeros(n); win1 = np.zeros(n)
    kk = min(k, n)
    topk_idx = np.argpartition(-core, kth=kk - 1, axis=1)[:, :kk]
    np.add.at(top_global, topk_idx.ravel(), 1)
    np.add.at(win1, core.argmax(axis=1), 1)
    for idx in groups:
        sub = core[:, idx]; m = len(idx); kg = min(k, m)
        member = ((-sub).argsort(axis=1).argsort(axis=1) < kg)
        np.add.at(top_lineage, np.tile(idx, N)[member.ravel()], 1)

    lo, hi = (1 - ci) / 2, 1 - (1 - ci) / 2
    ci_lo = np.quantile(core, lo, axis=0); ci_hi = np.quantile(core, hi, axis=0)
    pts = np.array([point.get(m, np.nan) for m in ids])
    pts = np.where(np.isfinite(pts), pts, np.median(core, axis=0))   # fallback (shouldn't trigger)
    ci_lo = np.minimum(ci_lo, pts); ci_hi = np.maximum(ci_hi, pts)

    return pd.DataFrame({
        "model_id": ids, "lineage": lin,
        "core_score": pts.astype("float32"),
        "ci_lo": ci_lo.astype("float32"), "ci_hi": ci_hi.astype("float32"),
        "ci_width": (ci_hi - ci_lo).astype("float32"),
        "topk_global_pct": (100 * top_global / N).astype("float32"),
        "topk_lineage_pct": (100 * top_lineage / N).astype("float32"),
        "top1_winrate": (100 * win1 / N).astype("float32"),
        "n_layers": np.where(has_p, 2, 1).astype("int8"),
        "n_sources": base["n_sources"].fillna(0).astype("int8").to_numpy(),
    })


# ═════════════════════════════════════════════════════════════════════════════
def run(genes=None, limit=None):
    t0 = time.time()
    print("=" * 70); print("Scoring — confidence intervals (batch)"); print("=" * 70)
    con = C.connect()
    rna, prot, core, stats_df, gene_tier = load_inputs(con)
    con.close()

    # protein se depends only on the gene tier + n_platforms
    prot = prot.rename(columns={"n_sources": "n_platforms"})

    all_genes = genes or sorted(core["ensg_id"].unique())
    if limit:
        all_genes = all_genes[:limit]
    print(f"Genes to score: {len(all_genes):,}  |  draws/gene: {N_DRAWS}")

    rna_g   = {g: d for g, d in rna.groupby("gene_id")}
    prot_g  = {g: d for g, d in prot.groupby("gene_id")}
    point_g = {g: dict(zip(d.model_id, d.core_score)) for g, d in core.groupby("ensg_id")}

    writer = None
    written = 0
    for i, g in enumerate(all_genes, 1):
        g_rna = rna_g.get(g)
        if g_rna is None:
            continue
        g_prot = prot_g.get(g, prot.iloc[0:0]).copy()
        if len(g_prot):
            mult = PROT_TIER_MULT.get(gene_tier.get(g, "cautious"), 1.5)
            prior = np.where(g_prot["n_platforms"] < 2, PROT_TAU_PRIOR, 0.0)
            g_prot = g_prot.assign(se_prot=np.sqrt((mult * PROT_BASE) ** 2 + prior ** 2))
            g_prot = g_prot[["model_id", "prot_z", "n_platforms", "se_prot"]]
        else:
            g_prot = pd.DataFrame(columns=["model_id", "prot_z", "n_platforms", "se_prot"])

        st = stats_df.loc[g] if g in stats_df.index else None
        beta = float(st["beta_g"]) if st is not None else 0.0
        med  = float(st["resid_med"]) if st is not None else 0.0
        sd   = float(st["resid_sd"]) if st is not None else 1.0

        out = gene_ci(g_rna[["model_id", "rna_z", "n_sources", "s_between", "lineage"]],
                      g_prot, point_g.get(g, {}), beta, med, sd)
        if out is None or out.empty:
            continue
        out.insert(0, "ensg_id", g)

        table = pa.Table.from_pandas(out, preserve_index=False)
        if writer is None:
            writer = pq.ParquetWriter(CI_SCORES, table.schema, compression="zstd")
        writer.write_table(table)
        written += len(out)

        if i % 500 == 0:
            el = time.time() - t0
            print(f"  {i:,}/{len(all_genes):,} genes  |  {written:,} rows  |  "
                  f"{el:.0f}s  |  {el/i*1000:.0f} ms/gene")
    if writer:
        writer.close()

    print(f"\nDone in {time.time()-t0:.0f}s")
    print(f"ci_scores: {written:,} rows  ->  {CI_SCORES}")


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--limit", type=int, default=None, help="score only the first N genes (smoke test)")
    ap.add_argument("--genes", nargs="*", default=None, help="specific ENSG/symbols (ENSG form)")
    args = ap.parse_args()
    run(genes=[g.upper() for g in args.genes] if args.genes else None, limit=args.limit)
