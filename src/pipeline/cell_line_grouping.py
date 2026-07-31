"""
cell_line_grouping.py
------------------------
Cell-line context and within-list grouping for a gene's returned ranking.
Two separate things, on purpose (see module docstring in the review that
motivated this file):

  (a) Descriptive context -- a display join on cell_line_lookup.parquet,
      already done in rank_cell_lines.py (diseases/category/sex/etc, keyed
      by model_id). Zero architectural risk, does not touch core_score.

  (b) Grouping the top-N returned lines -- NOT a re-ranking by similarity
      (that failed before, dominated by lineage). Instead: use lineage
      domination as the *feature* -- group and explain, don't reorder by a
      similarity score. Metadata grouping first (free, explainable); only
      escalate to expression-profile clustering if it adds something
      metadata grouping doesn't already give, and test that against the
      metadata grouping rather than assuming it's better.
"""
import re
import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import squareform


def disease_short_label(diseases_raw):
    """cell_line_lookup.diseases is 'ncit; c4004; gastric adenocarcinoma ||
    ordo; orphanet_...; ...' -- take the first NCIT term as a short label."""
    if pd.isna(diseases_raw) or not diseases_raw:
        return "unknown"
    first = diseases_raw.split("||")[0].strip()
    parts = [p.strip() for p in first.split(";")]
    return parts[-1] if parts else "unknown"


def group_top_n_by_metadata(ranking_rows, n=50):
    """ranking_rows: list of dicts from rank_cell_lines()['ranking'], already
    carrying 'diseases' from the cell_line_lookup join. Returns counts, most
    common first -- e.g. "of your top 20, 11 are melanoma, 4 lung, 3 colorectal"."""
    top = ranking_rows[:n]
    labels = [disease_short_label(r.get("diseases")) for r in top]
    counts = pd.Series(labels).value_counts()
    return {
        "n": len(top),
        "groups": counts.to_dict(),
        "model_id_to_group": dict(zip([r["model_id"] for r in top], labels)),
        "summary": ", ".join(f"{c} {g}" for g, c in counts.items()),
    }


def cluster_by_expression_profile(model_ids, expr_wide, n_clusters=None):
    """Hierarchical clustering on Pearson correlation of expression profiles,
    restricted to the returned set only (never used to reorder outside it).
    expr_wide: (model_id x ensg_id) DataFrame, values = log2 TPM.
    n_clusters: if None, cut the dendrogram at the same group count metadata
    grouping found, for a fair comparison."""
    sub = expr_wide.reindex(index=model_ids).dropna(how="all")
    if len(sub) < 3:
        return {mid: 0 for mid in model_ids}
    corr = sub.T.corr(method="pearson")
    dist = 1 - corr.values
    dist = (dist + dist.T) / 2
    np.fill_diagonal(dist, 0)
    condensed = squareform(dist, checks=False)
    Z = linkage(condensed, method="average")
    k = n_clusters or max(2, len(sub) // 5)
    labels = fcluster(Z, t=k, criterion="maxclust")
    return dict(zip(sub.index, labels))


def compare_clusterings(meta_labels: dict, profile_labels: dict):
    """Adjusted Rand Index between metadata grouping and profile clustering,
    over the model_ids both cover. High ARI = profile clustering just
    recovers the disease labels -- use the labels, they're cheaper."""
    from sklearn.metrics import adjusted_rand_score
    common = sorted(set(meta_labels) & set(profile_labels))
    if len(common) < 3:
        return {"ari": None, "n": len(common), "note": "too few overlapping lines"}
    ari = adjusted_rand_score([meta_labels[m] for m in common], [profile_labels[m] for m in common])
    return {"ari": ari, "n": len(common)}


def diversify_top_n(ranking_rows, group_key="diseases", pool_size=50, k=20):
    """Best-line-per-cluster-then-fill. Never touches core_score -- purely a
    presentation re-order of the existing ranked pool. Re-rank algorithm:
      1. Take the top `pool_size` candidates (already core_score-ranked).
      2. First pass: walk clusters in order of their best member's rank,
         taking one representative per cluster.
      3. Second pass: fill remaining slots with the next-best remaining
         candidates in original rank order, until k slots are filled.
    """
    pool = ranking_rows[:pool_size]
    for r in pool:
        r["_group"] = disease_short_label(r.get(group_key)) if group_key == "diseases" else r.get(group_key)

    seen_groups = set()
    diversified = []
    remaining = []
    for r in pool:  # already rank-ordered
        if r["_group"] not in seen_groups:
            diversified.append(r)
            seen_groups.add(r["_group"])
        else:
            remaining.append(r)
        if len(diversified) >= k:
            break

    for r in remaining:
        if len(diversified) >= k:
            break
        diversified.append(r)

    return diversified[:k]
