"""
ci_checks.py
--------------------
Verification and robustness checks for the confidence-interval / rank-stability
layer. Each function is self-contained and prints the headline result quoted in
the write-up; most run against the real warehouse over a gene sample.

Every check reuses the CI machinery defined in the confidence-intervals notebook:
  run_gene_ci, _gene_leaves, build_anchor, core_score_chain, monte_carlo_ci,
  _lookup, _resolve_gid, and the module-level constants (W_RNA, RNA_WINSOR_BOUND,
  REGIME3_SIG_THRESH, PROT_STD_WINSOR_BOUND, RNA_BASE_FLOOR, RNA_TAU_PRIOR,
  COVERAGE_CAP, PROT_BASE, ...).

Intended use: import these names into the notebook (or run the notebook's
definition cells first), then call the checks below. They are collected here so
the method and its validation live in one place.

Results obtained on the real pipeline (30-gene sample unless noted):
  centre_draw_identity        max|Δ| ≈ 1.8e-08   (EGFR, TP53, MYC)
  order_consistency           100% of 570 disjoint-CI pairs order-stable at ≥95%
  self_consistency            mean |median − centre| = 0.006
  convergence                 win-rate change ≤ 0.9 pp from N=1e3 to 1e4
  knob_sweep                  mean Spearman 0.91 under ±50% knob moves
  variance_decomposition      two-layer lines: protein ≈ 91% / RNA ≈ 9%
  regime3_flip                ≈ 67% of two-layer lines flip regime ≥ once
  correlated_error_sensitivity mean Spearman 0.93 (ρ=.25), 0.88 (ρ=.5)
  winsor_test                 saturation correlation unchanged ±5 → ±10 (Φ, not clip)

These functions do NOT define the CI method; they validate it. They expect the
notebook's names to be importable (e.g. `from confidence_intervals import *`) or
to be run inside the notebook kernel after the definition cells.
"""
import numpy as np
import pandas as pd
from scipy import stats


# ═════════════════════════════════════════════════════════════════════════════
# helpers shared by several checks
# ═════════════════════════════════════════════════════════════════════════════
def _sample_genes(CORE_SCORE, n=30, seed=0):
    """Draw a reproducible random sample of gene IDs from the core-score table.

    Parameters
    ----------
    CORE_SCORE : pathlib.Path
        Path to the pipeline's core-score parquet (must contain an "ensg_id" column).
    n : int, optional
        Number of genes to sample (default 30).
    seed : int, optional
        RNG seed for a reproducible draw (default 0).

    Returns
    -------
    list[str]
        Upper-cased Ensembl gene IDs, sampled without replacement.
    """
    core = pd.read_parquet(CORE_SCORE)
    core["ensg_id"] = core["ensg_id"].str.upper()
    return list(np.random.default_rng(seed).choice(core["ensg_id"].unique(), n, replace=False))


def _base_anchor(con, gene, run_gene_ci_ctx):
    """Build the (base, anchor) pair for a gene, exactly as run_gene_ci does internally.

    Reads the gene's leaves and noise, loads the panel-wide globals from
    ci_globals.json, builds the anchor, attaches lineage, and assembles the per-line
    `base` frame — the shared setup that the checks below drive their own simulations
    from, without repeating run_gene_ci's Monte-Carlo step.

    Parameters
    ----------
    con : database connection
        Open warehouse connection.
    gene : str
        Gene symbol or Ensembl ID.
    run_gene_ci_ctx : dict
        Notebook names needed to reproduce the wiring, with keys "_gene_leaves",
        "build_anchor" and "CORE_SCORE".

    Returns
    -------
    tuple[pandas.DataFrame, dict]
        (base, anchor), ready for the Monte-Carlo routines.
    """
    import json
    import common as C
    _gene_leaves = run_gene_ci_ctx["_gene_leaves"]
    build_anchor = run_gene_ci_ctx["build_anchor"]
    CORE_SCORE   = run_gene_ci_ctx["CORE_SCORE"]
    rz, pz, gid = _gene_leaves(con, gene)
    g_globals = json.loads((CORE_SCORE.parent / "ci_globals.json").read_text())
    anchor = build_anchor(rz, pz, g_globals)
    lm = C.load_lineage(con); lm.index = lm.index.str.lower()
    base = (rz.rename(columns={"n_sources": "n_sources_rna"})
              .merge(pz.rename(columns={"n_sources": "n_sources_prot"}), on="model_id", how="left"))
    base["lineage"] = base.model_id.map(lm)
    return base, anchor


# ═════════════════════════════════════════════════════════════════════════════
# 1. Centre-draw identity — recompute equals the stored core_score
# ═════════════════════════════════════════════════════════════════════════════
def centre_draw_identity(con, run_gene_ci, _resolve_gid, CORE_SCORE,
                         genes=("EGFR", "TP53", "MYC"), N=300):
    """Check 1 — the noise-free recompute reproduces the stored core_score.

    For each gene, runs the CI machinery and compares the centre (unperturbed) score
    against the pipeline's persisted core_score line-by-line, printing the maximum
    absolute deviation. A clean result is max|Δ| of order 1e-8, confirming the
    notebook's recompute path matches core_score.py exactly. This is the single most
    important check: everything downstream assumes the reconstruction is faithful.

    Parameters
    ----------
    con : database connection
        Open warehouse connection.
    run_gene_ci : callable
        The notebook's run_gene_ci(con, gene, N=...) entry point.
    _resolve_gid : callable
        The notebook's symbol→gene_id resolver.
    CORE_SCORE : pathlib.Path
        Path to the core-score parquet holding the reference scores.
    genes : tuple[str, ...], optional
        Genes to check (default EGFR, TP53, MYC).
    N : int, optional
        Monte-Carlo draws per gene; only the centre is used here (default 300).

    Returns
    -------
    None
        Results are printed, one line per gene.
    """
    core = pd.read_parquet(CORE_SCORE); core["ensg_id"] = core.ensg_id.str.upper()
    for gsym in genes:
        res, _, _ = run_gene_ci(con, gsym, N=N)
        gid = _resolve_gid(con, gsym).upper()
        ref = core[core.ensg_id == gid].assign(model_id=lambda d: d.model_id.str.lower())
        got = res.assign(model_id=lambda d: d.model_id.str.lower())
        m = got.merge(ref[["model_id", "core_score"]].rename(columns={"core_score": "r"}), on="model_id")
        dmax = (m.core_score - m.r).abs().max()
        n1 = int((res.n_layers == 1).sum())
        print(f"{gsym:6s} compared={len(m):5d}  n_layers1={n1:5d}  max|Δ|={dmax:.2e}")
    print("  (target: max|Δ| ~1e-8 => recompute reproduces the pipeline)")


# ═════════════════════════════════════════════════════════════════════════════
# 2. Order-consistency — disjoint-CI pairs keep their order across draws
# ═════════════════════════════════════════════════════════════════════════════
def order_consistency(con, genes, run_gene_ci, N=1000, k=10, level=0.95, topn=30):
    """Check 2 — disjoint-CI pairs keep their order across Monte-Carlo draws.

    Among the top-`topn` lines per gene, finds every pair whose 95% intervals do not
    overlap and measures how often the higher-scoring line also outscores the other
    within the raw per-draw scores. A pair "holds" if its order is preserved in at
    least `level` of draws; the pooled fraction across genes is printed. If separated
    intervals reliably imply a stable ranking, the intervals mean what they claim.

    Parameters
    ----------
    con : database connection
        Open warehouse connection.
    genes : iterable of str
        Genes to evaluate.
    run_gene_ci : callable
        The notebook's run_gene_ci entry point (returns table, score matrix, ids).
    N : int, optional
        Monte-Carlo draws per gene (default 1000).
    k : int, optional
        Top-k passed through to run_gene_ci (default 10).
    level : float, optional
        Order-preservation threshold a pair must clear to count as stable (default 0.95).
    topn : int, optional
        Number of top lines per gene to form pairs from (default 30).

    Returns
    -------
    pandas.DataFrame
        One row per gene with "n_pairs" (disjoint-CI pairs) and "frac" (share stable
        at >= level).
    """
    rows = []
    for gid in genes:
        res, scores, ids = run_gene_ci(con, gid, N=N, k=k)
        cs = res.set_index("model_id").loc[ids, ["core_score", "ci_lo", "ci_hi"]].to_numpy()
        sc, lo, hi = cs[:, 0], cs[:, 1], cs[:, 2]
        order = np.argsort(-sc)[:topn]
        preserved = []
        for i in range(len(order)):
            for j in range(i + 1, len(order)):
                a, b = order[i], order[j]
                if hi[a] < lo[b] or hi[b] < lo[a]:
                    hi_, lo_ = (a, b) if sc[a] >= sc[b] else (b, a)
                    preserved.append(np.mean(scores[:, hi_] > scores[:, lo_]) >= level)
        if preserved:
            rows.append({"gene": gid, "n_pairs": len(preserved), "frac": np.mean(preserved)})
    df = pd.DataFrame(rows)
    if len(df):
        tot = (df.n_pairs * df.frac).sum() / df.n_pairs.sum()
        print(f"disjoint-CI pairs: {int(df.n_pairs.sum())} across {len(df)} genes")
        print(f"pairs whose order holds in >={level:.0%} of draws: {tot:.1%}")
    return df


# ═════════════════════════════════════════════════════════════════════════════
# 3. Self-consistency — draw median vs stored centre, by score bin
# ═════════════════════════════════════════════════════════════════════════════
def self_consistency(con, genes, run_gene_ci, N=1000):
    """Check 3 — the Monte-Carlo median tracks the stored centre, by score bin.

    Pools every line across `genes` and compares the per-line draw median to the
    stored centre score, bucketed into score bins. Small, roughly unbiased gaps show
    the perturbation does not systematically shift scores away from the point
    estimate (e.g. through the non-linearity of the normal-CDF transform near 0/1).

    Parameters
    ----------
    con : database connection
        Open warehouse connection.
    genes : iterable of str
        Genes to pool over.
    run_gene_ci : callable
        The notebook's run_gene_ci entry point.
    N : int, optional
        Monte-Carlo draws per gene (default 1000).

    Returns
    -------
    pandas.DataFrame
        Per-line "core_score", "gap" (median − centre) and the assigned score "bin".
        A by-bin summary and the overall mean |gap| are printed.
    """
    pool = []
    for gid in genes:
        res, scores, ids = run_gene_ci(con, gid, N=N)
        med = np.median(scores, axis=0)
        centre = res.set_index("model_id").loc[ids, "core_score"].to_numpy()
        pool.append(pd.DataFrame({"core_score": centre, "gap": med - centre}))
    P = pd.concat(pool, ignore_index=True)
    P["bin"] = pd.cut(P.core_score, [0, .2, .4, .6, .8, 1.0])
    print("draw-median minus stored centre, by score bin:")
    print(P.groupby("bin", observed=True)["gap"].agg(["mean", "std", "count"]).round(4).to_string())
    print(f"\noverall mean |gap|: {P.gap.abs().mean():.4f}")
    return P


# ═════════════════════════════════════════════════════════════════════════════
# 4. Convergence — estimates stable as N grows
# ═════════════════════════════════════════════════════════════════════════════
def convergence(con, gene, run_gene_ci, N_grid=(100, 300, 1000, 3000, 10000), k=10, seed=42, plot=True):
    """Check 4 — win-rate and CI width stabilise as the draw count N grows.

    Re-runs a single gene over an increasing grid of draw counts, tracking the top-1
    win-rate and CI width of its top-5 lines. Optionally plots both against N (log
    scale) and prints the mean absolute win-rate change from N=1e3 to 1e4, which
    should be small (order ~1 pp) once the Monte-Carlo estimates have converged.

    Parameters
    ----------
    con : database connection
        Open warehouse connection.
    gene : str
        Gene symbol or Ensembl ID.
    run_gene_ci : callable
        The notebook's run_gene_ci entry point.
    N_grid : tuple[int, ...], optional
        Draw counts to sweep (default 100, 300, 1000, 3000, 10000).
    k : int, optional
        Top-k passed through to run_gene_ci (default 10).
    seed : int, optional
        RNG seed, held fixed across the grid (default 42).
    plot : bool, optional
        Whether to draw the win-rate/width-vs-N figure (default True).

    Returns
    -------
    pandas.DataFrame
        Long-format "N", "model_id", "top1_winrate", "ci_width" for the tracked lines.
    """
    rows = []
    for N in N_grid:
        res, _, _ = run_gene_ci(con, gene, N=N, k=k, seed=seed)
        for _, r in res.head(5).iterrows():
            rows.append({"N": N, "model_id": r.model_id, "top1_winrate": r.top1_winrate,
                         "ci_width": r.ci_width})
    df = pd.DataFrame(rows)
    if plot:
        import matplotlib.pyplot as plt
        fig, (a1, a2) = plt.subplots(1, 2, figsize=(11, 4.2))
        for mid, gg in df.groupby("model_id"):
            a1.plot(gg.N, gg.top1_winrate, "o-", label=mid, alpha=.8)
            a2.plot(gg.N, gg.ci_width, "o-", alpha=.8)
        a1.set_xscale("log"); a2.set_xscale("log")
        a1.set_xlabel("N draws"); a1.set_ylabel("top-1 win-rate (%)"); a1.set_title("Win-rate vs N")
        a2.set_xlabel("N draws"); a2.set_ylabel("95% CI width"); a2.set_title("CI width vs N")
        a1.legend(fontsize=7, frameon=False); plt.tight_layout(); plt.show()
    piv = df.pivot_table(index="model_id", columns="N", values="top1_winrate")
    if 1000 in piv and 10000 in piv:
        print(f"mean |win-rate change| 1000->10000: {(piv[10000]-piv[1000]).abs().mean():.2f} pp")
    return df


# ═════════════════════════════════════════════════════════════════════════════
# 5. Knob sweep — win-rate ranking under ±50% moves of each noise knob
# ═════════════════════════════════════════════════════════════════════════════
def knob_sweep(con, genes, run_gene_ci, knob_names=("RNA_BASE_FLOOR", "PROT_BASE", "RNA_TAU_PRIOR"),
               N=400, k=10, ns=None):
    """Check 5 — win-rate ranking is stable under ±50% moves of each noise knob.

    For each gene and each named knob, perturbs the knob to 0.5× and 1.5× its value,
    re-runs the CI, and correlates (Spearman) the resulting top-1 win-rate ranking
    against the unperturbed one; the knob is restored after each move. High mean
    Spearman means the deliberately conservative knob settings are not what drives
    the ranking, so the exact values are not load-bearing.

    Parameters
    ----------
    con : database connection
        Open warehouse connection.
    genes : iterable of str
        Genes to evaluate.
    run_gene_ci : callable
        The notebook's run_gene_ci entry point.
    knob_names : tuple[str, ...], optional
        Noise knobs to sweep (default RNA_BASE_FLOOR, PROT_BASE, RNA_TAU_PRIOR).
    N : int, optional
        Monte-Carlo draws per run (default 400).
    k : int, optional
        Top-k passed through to run_gene_ci (default 10).
    ns : dict
        The namespace holding the knobs, usually the notebook's globals(); required
        so the knobs can be perturbed in place and restored.

    Returns
    -------
    pandas.DataFrame
        One row per (gene, knob, factor) with the Spearman correlation to baseline;
        by-knob means and the overall mean are printed.

    Raises
    ------
    ValueError
        If `ns` is not provided.
    """
    if ns is None:
        raise ValueError("pass ns=globals() so the noise knobs can be perturbed/restored")
    base_vals = {kk: ns[kk] for kk in knob_names}
    def _wr(gid):
        """Top-1 win-rate per line for one gene, under the current knob values."""
        r, _, _ = run_gene_ci(con, gid, N=N, k=k)
        return r.set_index("model_id")["top1_winrate"]
    out = []
    for gid in genes:
        base = _wr(gid)
        for kk in knob_names:
            for factor in (0.5, 1.5):
                ns[kk] = base_vals[kk] * factor
                try:
                    alt = _wr(gid).reindex(base.index)
                    out.append({"gene": gid, "knob": kk, "factor": factor,
                                "spearman": base.corr(alt, method="spearman")})
                finally:
                    ns[kk] = base_vals[kk]
    df = pd.DataFrame(out)
    print("win-rate ranking stability (Spearman) under ±50% knob moves:")
    print(df.groupby(["knob", "factor"])["spearman"].mean().round(3).to_string())
    print(f"\noverall mean Spearman: {df.spearman.mean():.3f}")
    return df


# ═════════════════════════════════════════════════════════════════════════════
# 6. Variance decomposition — RNA vs protein share of score variance
# ═════════════════════════════════════════════════════════════════════════════
def variance_decomposition(con, gene, ctx, N=2000, seed=42):
    """Check 6 — split score variance into RNA and protein contributions.

    For one gene, runs the Monte-Carlo perturbation three ways — RNA-only,
    protein-only and both — and, on two-layer lines, reports the first-order variance
    shares normalised to sum to 100%. Quantifies how much of a line's score
    uncertainty each evidence layer contributes (typically protein-dominated).

    Parameters
    ----------
    con : database connection
        Open warehouse connection.
    gene : str
        Gene symbol or Ensembl ID.
    ctx : dict
        Notebook names, with keys "_gene_leaves", "build_anchor", "core_score_chain"
        and "CORE_SCORE".
    N : int, optional
        Monte-Carlo draws per variance estimate (default 2000).
    seed : int, optional
        RNG seed (default 42).

    Returns
    -------
    dict
        "v_rna", "v_prot", "v_both" (per-line score variance under each perturbation
        scheme) and "has_p" (two-layer mask). The normalised shares are printed.
    """
    import json
    _gene_leaves = ctx["_gene_leaves"]; build_anchor = ctx["build_anchor"]
    core_score_chain = ctx["core_score_chain"]; CORE_SCORE = ctx["CORE_SCORE"]
    g_globals = json.loads((CORE_SCORE.parent / "ci_globals.json").read_text())
    rz, pz, gid = _gene_leaves(con, gene)
    anchor = build_anchor(rz, pz, g_globals)
    base = (rz.rename(columns={"n_sources": "n_sources_rna"})
              .merge(pz.rename(columns={"n_sources": "n_sources_prot"}), on="model_id", how="left"))
    n = len(base); has_p = base.prot_z.notna().to_numpy()
    rna_mu, rna_sd = base.rna_z.to_numpy(), base.se_rna.to_numpy()
    prot_mu, prot_sd = base.prot_z.to_numpy(float), np.nan_to_num(base.se_prot.to_numpy(float))
    nsr, nsp = base.n_sources_rna.to_numpy(float), base.n_sources_prot.to_numpy(float)
    rng = np.random.default_rng(seed)

    def var_under(pr_rna, pr_prot):
        """Per-line score variance with RNA and/or protein perturbation switched on."""
        S = np.empty((N, n))
        for i in range(N):
            rz_ = rna_mu + (rng.normal(0, 1, n) * rna_sd if pr_rna else 0)
            pz_ = prot_mu.copy()
            if pr_prot:
                pz_[has_p] = prot_mu[has_p] + rng.normal(0, 1, int(has_p.sum())) * prot_sd[has_p]
            S[i] = core_score_chain(rz_, pz_, nsr, nsp, has_p, anchor)
        return S.var(axis=0)

    v_rna, v_prot, v_both = var_under(True, False), var_under(False, True), var_under(True, True)
    two = has_p; tot = np.where(v_both[two] > 0, v_both[two], np.nan)
    # normalise the two first-order terms to sum to 100%
    s_rna = np.nanmean(v_rna[two] / (v_rna[two] + v_prot[two]))
    s_prot = 1 - s_rna
    print(f"two-layer lines: {int(two.sum())}")
    print(f"  normalised share — RNA: {s_rna:.1%}  |  Prot: {s_prot:.1%}")
    return {"v_rna": v_rna, "v_prot": v_prot, "v_both": v_both, "has_p": has_p}


# ═════════════════════════════════════════════════════════════════════════════
# 7. Regime-3 flip — fraction of two-layer lines that flip regime >= once
# ═════════════════════════════════════════════════════════════════════════════
def _regime(rna_z, prot_z, nsrc_rna, nsrc_prot, has_p, a, ctx):
    """Flag each two-layer line as conflicting (regime 3) or not, for given z-scores.

    Reproduces the regime-3 test inside core_score_chain: a line is flagged when its
    standardised, winsorised RNA and protein signals both exceed REGIME3_SIG_THRESH
    in magnitude and disagree in sign. Used by `regime3_flip` to detect regime changes
    under perturbation.

    Parameters
    ----------
    rna_z, prot_z : numpy.ndarray
        Per-line RNA and protein z-scores.
    nsrc_rna, nsrc_prot : numpy.ndarray
        Per-line source counts (select the standardisation stratum).
    has_p : numpy.ndarray of bool
        Two-layer mask.
    a : dict
        Per-gene anchor.
    ctx : dict
        Notebook names/constants, with keys "_lookup", "RNA_WINSOR_BOUND",
        "PROT_STD_WINSOR_BOUND" and "REGIME3_SIG_THRESH".

    Returns
    -------
    numpy.ndarray of bool
        True where the line falls in the conflicting (regime-3) branch.
    """
    _lookup = ctx["_lookup"]
    RNA_W = ctx["RNA_WINSOR_BOUND"]; PROT_W = ctx["PROT_STD_WINSOR_BOUND"]; THR = ctx["REGIME3_SIG_THRESH"]
    rmed, rsd = _lookup(nsrc_rna, a["rna_str"])
    rna_std = np.clip((rna_z - rmed) / rsd, -RNA_W, RNA_W)
    reg = np.zeros(len(rna_z), bool)
    if has_p.any():
        pmed, psd = _lookup(np.where(has_p, nsrc_prot, np.nan), a["protstd_str"])
        with np.errstate(invalid="ignore"):
            pstd = np.clip((prot_z - pmed) / psd, -PROT_W, PROT_W)
            reg = ((np.abs(rna_std) > THR) & (np.abs(pstd) > THR)
                   & (np.sign(rna_std) != np.sign(pstd)) & np.isfinite(pstd) & has_p)
    return reg


def regime3_flip(con, genes, ctx, N=800, seed=42):
    """Check 7 — how often two-layer lines switch scoring regime under perturbation.

    For each gene, records each two-layer line's regime at the centre, then perturbs
    RNA and protein across N draws and flags any line whose regime ever differs from
    its centre value. Reports the mean fraction of two-layer lines that flip at least
    once — a high fraction shows the conflicting/consistent branch boundary sits near
    the data for many lines, so the regime logic is genuinely exercised.

    Parameters
    ----------
    con : database connection
        Open warehouse connection.
    genes : iterable of str
        Genes to evaluate.
    ctx : dict
        Notebook names/constants: must include "_gene_leaves", "build_anchor",
        "_lookup", "CORE_SCORE" and the winsor/threshold constants used by `_regime`.
    N : int, optional
        Monte-Carlo draws per gene (default 800).
    seed : int, optional
        RNG seed (default 42).

    Returns
    -------
    pandas.DataFrame
        One row per gene with "ever_flip_frac"; the mean across genes is printed.
    """
    import json
    _gene_leaves = ctx["_gene_leaves"]; build_anchor = ctx["build_anchor"]; CORE_SCORE = ctx["CORE_SCORE"]
    g_globals = json.loads((CORE_SCORE.parent / "ci_globals.json").read_text())
    rng = np.random.default_rng(seed); rows = []
    for gid in genes:
        rz, pz, _ = _gene_leaves(con, gid)
        anchor = build_anchor(rz, pz, g_globals)
        base = (rz.rename(columns={"n_sources": "n_sources_rna"})
                  .merge(pz.rename(columns={"n_sources": "n_sources_prot"}), on="model_id", how="left"))
        has_p = base.prot_z.notna().to_numpy()
        if has_p.sum() == 0:
            continue
        rna_mu, rna_sd = base.rna_z.to_numpy(), base.se_rna.to_numpy()
        prot_mu, prot_sd = base.prot_z.to_numpy(float), np.nan_to_num(base.se_prot.to_numpy(float))
        nsr, nsp = base.n_sources_rna.to_numpy(float), base.n_sources_prot.to_numpy(float)
        n = len(base)
        centre_reg = _regime(rna_mu, prot_mu, nsr, nsp, has_p, anchor, ctx)
        ever = np.zeros(n, bool)
        for _ in range(N):
            rz_ = rna_mu + rng.normal(0, 1, n) * rna_sd
            pz_ = prot_mu.copy()
            pz_[has_p] = prot_mu[has_p] + rng.normal(0, 1, int(has_p.sum())) * prot_sd[has_p]
            ever |= (_regime(rz_, pz_, nsr, nsp, has_p, anchor, ctx) != centre_reg)
        rows.append({"gene": gid, "ever_flip_frac": float(ever[has_p].mean())})
    df = pd.DataFrame(rows)
    print(f"mean fraction of two-layer lines that flip regime >= once: {df.ever_flip_frac.mean():.1%}")
    return df


# ═════════════════════════════════════════════════════════════════════════════
# 8. Correlated-error sensitivity — win-rate ranking under correlated draws
# ═════════════════════════════════════════════════════════════════════════════
def _mc_corr(base, anchor, ctx, N=500, k=10, seed=1, rho_err=0.0, ci=0.95):
    """Monte-Carlo top-1 win-rate per line with correlated RNA/protein errors.

    A trimmed variant of the CI simulation that returns only the win-rate, used by
    `correlated_error_sensitivity`. The protein perturbation is drawn with correlation
    `rho_err` to the RNA perturbation via e_prot = rho·e_rna + sqrt(1 − rho²)·e_ind, so
    rho_err=0 reproduces the independent baseline.

    Parameters
    ----------
    base : pandas.DataFrame
        Per-line frame from `_base_anchor`.
    anchor : dict
        Per-gene anchor.
    ctx : dict
        Notebook names; requires "core_score_chain".
    N : int, optional
        Monte-Carlo draws (default 500).
    k : int, optional
        Unused; kept for signature parity with the full simulation (default 10).
    seed : int, optional
        RNG seed (default 1).
    rho_err : float, optional
        Target error correlation, clipped to (−0.999, 0.999). Default 0.0.
    ci : float, optional
        Unused; kept for signature parity with the full simulation (default 0.95).

    Returns
    -------
    pandas.Series
        Top-1 win-rate (%) indexed by model_id.
    """
    core_score_chain = ctx["core_score_chain"]
    rng = np.random.default_rng(seed)
    ids = base.model_id.to_numpy(); n = len(base); lineage = base.lineage.to_numpy()
    rna_mu, rna_sd = base.rna_z.to_numpy(), base.se_rna.to_numpy()
    nsr = base.n_sources_rna.to_numpy(float)
    has_p = base.prot_z.notna().to_numpy()
    prot_mu = base.prot_z.to_numpy(float); prot_sd = np.nan_to_num(base.se_prot.to_numpy(float))
    nsp = base.n_sources_prot.to_numpy(float)
    win1 = np.zeros(n); rho = float(np.clip(rho_err, -0.999, 0.999))
    for i in range(N):
        e_rna = rng.normal(0, 1, n); e_ind = rng.normal(0, 1, n)
        e_prot = rho * e_rna + np.sqrt(1 - rho**2) * e_ind
        rz = rna_mu + e_rna * rna_sd
        pz = prot_mu.copy(); pz[has_p] = prot_mu[has_p] + e_prot[has_p] * prot_sd[has_p]
        s = core_score_chain(rz, pz, nsr, nsp, has_p, anchor)
        win1[np.argmax(s)] += 1
    return pd.Series(100 * win1 / N, index=ids)


def correlated_error_sensitivity(con, genes, ctx, rhos=(0.0, 0.25, 0.5), N=400):
    """Check 8 — win-rate ranking is stable under correlated measurement errors.

    For each gene, compares the top-1 win-rate ranking at each `rho_err` against the
    independent (rho=0) baseline by Spearman correlation, also tracking the largest
    single-line win-rate shift in percentage points. Spearman near 1 with small shifts
    means the independence assumption in the main method is not load-bearing. Genes
    that fail to load are skipped with a message.

    Parameters
    ----------
    con : database connection
        Open warehouse connection.
    genes : iterable of str
        Genes to sweep.
    ctx : dict
        Notebook names, with keys "_gene_leaves", "build_anchor", "core_score_chain"
        and "CORE_SCORE".
    rhos : tuple[float, ...], optional
        Error correlations to test; must include 0.0 as the baseline
        (default 0.0, 0.25, 0.5).
    N : int, optional
        Monte-Carlo draws per run (default 400).

    Returns
    -------
    pandas.DataFrame
        One row per (gene, non-zero rho) with "spearman" and "max_shift_pp"; a by-rho
        summary is printed.
    """
    rows = []
    for gid in genes:
        try:
            base, anchor = _base_anchor(con, gid, ctx)
        except Exception as e:
            print("skip", gid, str(e)[:40]); continue
        base_wr = None
        for rho in rhos:
            wr = _mc_corr(base, anchor, ctx, N=N, seed=1, rho_err=rho)
            if rho == 0.0:
                base_wr = wr
            else:
                rho_s = base_wr.corr(wr.reindex(base_wr.index), method="spearman")
                shift = (base_wr - wr.reindex(base_wr.index)).abs().max()
                rows.append({"gene": gid, "rho_err": rho, "spearman": rho_s, "max_shift_pp": shift})
    df = pd.DataFrame(rows)
    print("win-rate ranking stability vs the independent baseline:")
    print(df.groupby("rho_err").agg(mean_spearman=("spearman", "mean"),
                                    min_spearman=("spearman", "min"),
                                    mean_max_shift_pp=("max_shift_pp", "mean")).round(3).to_string())
    return df


# ═════════════════════════════════════════════════════════════════════════════
# 9. Winsor test — is the saturation dominance the Φ transform or the clip?
# ═════════════════════════════════════════════════════════════════════════════
def winsor_test(con, genes, run_gene_ci, RNA_Z, _resolve_gid, ns, bounds=(5.0, 10.0), N=300):
    """Check 9 — is the score saturation driven by the Φ transform or the winsor clip?

    Temporarily widens all three winsor bounds (e.g. ±5 → ±10) and checks whether the
    relationship between CI width and score saturation, and between CI width and
    measured source disagreement, survives. If those correlations are essentially
    unchanged, saturation near 0/1 comes from the normal-CDF transform rather than the
    clipping. The original bounds are always restored in a finally block.

    Parameters
    ----------
    con : database connection
        Open warehouse connection.
    genes : iterable of str
        Genes to pool over.
    run_gene_ci : callable
        The notebook's run_gene_ci entry point.
    RNA_Z : pathlib.Path
        Path to the RNA z-score parquet; its sibling bulk_rna_z_by_source.parquet
        supplies the per-source disagreement (src_sd).
    _resolve_gid : callable
        The notebook's symbol→gene_id resolver.
    ns : dict
        Namespace holding the winsor bounds (usually the notebook globals()); the
        three bounds are rebound in place and restored on exit.
    bounds : tuple[float, ...], optional
        Winsor bounds to test (default 5.0, 10.0).
    N : int, optional
        Monte-Carlo draws per run (default 300).

    Returns
    -------
    dict
        Maps each bound to {"corr_saturation", "corr_srcsd", "mean_width"}; a summary
        table is printed.
    """
    orig = (ns["RNA_WINSOR_BOUND"], ns["PROT_STD_WINSOR_BOUND"], ns["PROT_RESID_WINSOR_BOUND"])
    rbs = pd.read_parquet(RNA_Z.parent / "bulk_rna_z_by_source.parquet")
    rbs["gene_id"] = rbs.gene_id.str.upper(); rbs["model_id"] = rbs.model_id.str.lower()
    rbs["src_sd"] = rbs[["z_depmap", "z_hpa_rna", "z_geo"]].std(axis=1, ddof=1)
    results = {}
    try:
        for B in bounds:
            ns["RNA_WINSOR_BOUND"] = ns["PROT_STD_WINSOR_BOUND"] = ns["PROT_RESID_WINSOR_BOUND"] = B
            pool = []
            for gid in genes:
                r, _, _ = run_gene_ci(con, gid, N=N)
                r = r.assign(model_id=r.model_id.str.lower())
                g = rbs[rbs.gene_id == _resolve_gid(con, gid).upper()][["model_id", "src_sd"]]
                pool.append(r.merge(g, on="model_id", how="left"))
            P = pd.concat(pool, ignore_index=True); P["dist"] = (P.core_score - 0.5).abs()
            results[B] = {"corr_saturation": P.dist.corr(-P.ci_width),
                          "corr_srcsd": P[P.src_sd.notna()].ci_width.corr(P[P.src_sd.notna()].src_sd),
                          "mean_width": P.ci_width.mean()}
    finally:
        ns["RNA_WINSOR_BOUND"], ns["PROT_STD_WINSOR_BOUND"], ns["PROT_RESID_WINSOR_BOUND"] = orig
    print("winsor | corr(width,saturation) | corr(width,src_sd) | mean width")
    for B, v in results.items():
        print(f"  ±{B:<5} | {v['corr_saturation']:+.3f} | {v['corr_srcsd']:+.3f} | {v['mean_width']:.3f}")
    return results