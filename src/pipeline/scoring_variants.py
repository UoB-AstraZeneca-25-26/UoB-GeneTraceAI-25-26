"""
Parameterised scoring variants for GeneTraceAI Layer 2.

Produces 12 (normalisation x combination) core_score outputs from the same
(expr_wide, prot_wide) inputs. Used by the GDSC validation harness to pick
the best-performing pair empirically.

Invariants across all variants:
  - Missing layer contributes no term (not zero) — reindex to union of
    models/genes so structural absences remain NaN
  - 1-layer cells return that layer's normalised score directly
  - stratum_rank computed per (ensg_id, n_layers)
  - Output schema: model_id, ensg_id, core_score, n_layers, stratum_rank
"""
import numpy as np
import pandas as pd

NORMALISATIONS = ["minmax", "percentile", "zscore", "hill"]
COMBINATIONS   = ["noisy_or", "weighted_sum", "product_floor"]


def normalise(df, method):
    """Per-column normalisation to [0,1], NaN preserved.
    df: (model_id x ensg_id) wide matrix."""
    if method == "minmax":
        mins = df.min(axis=0)
        maxs = df.max(axis=0)
        span = (maxs - mins).replace(0, np.nan)
        return (df - mins) / span
    if method == "percentile":
        return df.rank(axis=0, pct=True, na_option="keep")
    if method == "zscore":
        mu = df.mean(axis=0)
        sd = df.std(axis=0).replace(0, np.nan)
        z = (df - mu) / sd
        return 1 / (1 + np.exp(-z))            # logistic -> [0,1]
    if method == "hill":
        # Percentile first, then Hill with placeholder k=2, p0=0.5
        # Uncalibrated -- post-validation tuning would fit k, p0
        pct = df.rank(axis=0, pct=True, na_option="keep")
        k, p0 = 2.0, 0.5
        return (pct ** k) / (pct ** k + p0 ** k)
    raise ValueError(f"unknown normalisation: {method}")


def combine(E, P, method):
    """Combine two aligned (model x gene) matrices, both in [0,1] with NaN
    marking absent layers. Returns (core_score, n_layers).
    Missing-layer cells fall through to the present layer's value."""
    n_layers = (~E.isna()).astype(int) + (~P.isna()).astype(int)
    both   = ~E.isna() & ~P.isna()
    e_only = ~E.isna() &  P.isna()
    p_only =  E.isna() & ~P.isna()

    core = pd.DataFrame(np.nan, index=E.index, columns=E.columns, dtype=float)
    core[e_only] = E[e_only]
    core[p_only] = P[p_only]

    if method == "noisy_or":
        core[both] = 1 - (1 - E[both]) * (1 - P[both])
    elif method == "weighted_sum":
        core[both] = 0.5 * E[both] + 0.5 * P[both]
    elif method == "product_floor":
        # Product but never below half the weaker layer -- prevents one
        # weak layer from annihilating a strong one
        prod = E[both] * P[both]
        floor = np.minimum(E[both], P[both]) * 0.5
        core[both] = np.maximum(prod, floor)
    else:
        raise ValueError(f"unknown combination: {method}")

    return core, n_layers


def combine_three(E, P, C, rho_ep=0.46, rho_ec=0.005, rho_pc=0.005):
    """Combine three aligned (model x gene) matrices, all in [0,1] with NaN
    marking absent layers. Uses correlation-penalised weighted sum for cells
    with all three layers; falls back to two-layer weighted_sum or single
    layer as appropriate.

    Weights are derived from mean absolute pairwise correlation:
        w_i = 1 / (1 + mean|rho with other layers|), then normalised.

    Default rho values from empirical measurement:
        rho_ep = 0.46  (expr-prot, from prior validation)
        rho_ec = 0.005 (expr-chronos, measured in Improvement 2)
        rho_pc = 0.005 (prot-chronos, estimated equal to rho_ec)
    """
    # Per-layer mean abs correlation with others
    mean_rho_e = (abs(rho_ep) + abs(rho_ec)) / 2
    mean_rho_p = (abs(rho_ep) + abs(rho_pc)) / 2
    mean_rho_c = (abs(rho_ec) + abs(rho_pc)) / 2

    w_e_raw = 1 / (1 + mean_rho_e)
    w_p_raw = 1 / (1 + mean_rho_p)
    w_c_raw = 1 / (1 + mean_rho_c)
    total   = w_e_raw + w_p_raw + w_c_raw
    w_e, w_p, w_c = w_e_raw / total, w_p_raw / total, w_c_raw / total

    all_models = sorted(set(E.index) | set(P.index) | set(C.index))
    all_genes  = sorted(set(E.columns) | set(P.columns) | set(C.columns))
    E = E.reindex(index=all_models, columns=all_genes)
    P = P.reindex(index=all_models, columns=all_genes)
    C = C.reindex(index=all_models, columns=all_genes)

    has_e = ~E.isna()
    has_p = ~P.isna()
    has_c = ~C.isna()
    n_layers = has_e.astype(int) + has_p.astype(int) + has_c.astype(int)

    core = pd.DataFrame(np.nan, index=E.index, columns=E.columns, dtype=float)

    # 1-layer fallbacks
    core[has_e & ~has_p & ~has_c] = E[has_e & ~has_p & ~has_c]
    core[~has_e & has_p & ~has_c] = P[~has_e & has_p & ~has_c]
    core[~has_e & ~has_p & has_c] = C[~has_e & ~has_p & has_c]

    # 2-layer: equal-weight among present layers (re-normalise pair weights)
    ep = has_e & has_p & ~has_c
    ec = has_e & ~has_p & has_c
    pc = ~has_e & has_p & has_c
    w_ep_e = w_e_raw / (w_e_raw + w_p_raw)
    w_ep_p = w_p_raw / (w_e_raw + w_p_raw)
    w_ec_e = w_e_raw / (w_e_raw + w_c_raw)
    w_ec_c = w_c_raw / (w_e_raw + w_c_raw)
    w_pc_p = w_p_raw / (w_p_raw + w_c_raw)
    w_pc_c = w_c_raw / (w_p_raw + w_c_raw)
    core[ep] = w_ep_e * E[ep] + w_ep_p * P[ep]
    core[ec] = w_ec_e * E[ec] + w_ec_c * C[ec]
    core[pc] = w_pc_p * P[pc] + w_pc_c * C[pc]

    # 3-layer: correlation-penalised weighted sum
    all3 = has_e & has_p & has_c
    core[all3] = w_e * E[all3] + w_p * P[all3] + w_c * C[all3]

    return core, n_layers, (w_e, w_p, w_c)


def score(expr_wide, prot_wide, normalisation, combination):
    """Full scoring pipeline for one (normalisation, combination) pair.

    expr_wide, prot_wide: (model_id x ensg_id) with model_id UPPERCASE,
                          ensg_id bare UPPER. Duplicate model_ids should
                          already be collapsed upstream (dedup_mean).

    Returns long-form DataFrame with columns:
      model_id, ensg_id, core_score, n_layers, stratum_rank
    """
    E = normalise(expr_wide, normalisation)
    P = normalise(prot_wide, normalisation)

    all_models = sorted(set(E.index)   | set(P.index))
    all_genes  = sorted(set(E.columns) | set(P.columns))
    E = E.reindex(index=all_models, columns=all_genes)
    P = P.reindex(index=all_models, columns=all_genes)

    core, n_layers = combine(E, P, combination)

    core.index.name   = "model_id"
    core.columns.name = "ensg_id"
    n_layers.index.name   = "model_id"
    n_layers.columns.name = "ensg_id"

    core_long = core.stack().reset_index().rename(columns={0: "core_score"})
    n_long    = n_layers.stack().reset_index().rename(columns={0: "n_layers"})
    out = core_long.merge(n_long, on=["model_id", "ensg_id"], how="left")
    assert len(out) == len(core_long), "fan-out on n_layers merge"

    out = out.dropna(subset=["core_score"])
    out["n_layers"] = out["n_layers"].astype(int)

    out["stratum_rank"] = (
        out.groupby(["ensg_id", "n_layers"])["core_score"]
           .rank(pct=True, method="average")
    )
    return out
