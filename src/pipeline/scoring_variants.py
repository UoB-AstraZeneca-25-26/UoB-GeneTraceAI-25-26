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
