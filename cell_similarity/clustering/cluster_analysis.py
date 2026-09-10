"""
cell_similarity/clustering/cluster_analysis.py
-----------------------------------------------
CoWork v8.0 — Exploratory clustering on the validated RNA similarity graph.

Question: do RNA-similarity clusters match tissue of origin, or do they cut
across tissue in a biologically informative way?

READ-ONLY with respect to the parent module. Writes only to
  cell_similarity/clustering/results/

Stages:
  C1 — confirm inputs (matrix shape, lineage source, tissue distribution)
  C2 — hierarchical + Louvain clustering
  C3 — ARI/NMI vs tissue labels, null baseline, per-cluster composition
  C4 — sanity checks for any cross-tissue clusters flagged in C3
  C5 — write results note and architecture doc paragraph
"""
from __future__ import annotations

import json
import sys
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd

# ── path setup ──────────────────────────────────────────────────────────────
ROOT = Path(__file__).resolve().parent.parent.parent
HERE = Path(__file__).resolve().parent
RESULTS = HERE / "results"
RESULTS.mkdir(exist_ok=True)

SIM_PATH = ROOT / "cell_similarity" / "outputs" / "rna_similarity.parquet"

# Add cell_similarity to path so we can import common.py
sys.path.insert(0, str(ROOT / "cell_similarity"))
import common as C                     # noqa: E402

# ── constants ────────────────────────────────────────────────────────────────
SEED = 42
N_NULLS = 1000
SMALL_TISSUE_THRESHOLD = 10
CROSS_TISSUE_MAJORITY = 0.50          # cluster flagged if no tissue > this
PURE_TISSUE_THRESHOLD = 0.90          # cluster flagged as nearly-pure above this
LOUVAIN_K_VALUES = [10, 15, 20, 30]   # k for the kNN graph in community detection
HC_LINKAGE = "average"                # Ward requires Euclidean; average works on correlation-dist


def banner(s: str):
    print("=" * 70)
    print(s)
    print("=" * 70)


# ═══════════════════════════════════════════════════════════════════════ C1
banner("C1 — CONFIRM INPUTS")

print("\nLoading similarity matrix …")
S = pd.read_parquet(SIM_PATH)
assert S.shape[0] == S.shape[1], "Matrix is not square"
n_lines = S.shape[0]
print(f"  Matrix: {n_lines} × {n_lines}  (diagonal = NaN)")
assert n_lines == 1673, f"Expected 1673 lines, got {n_lines}"
print("  ✓ 1,673 lines confirmed")

print("\nLoading lineage labels …")
lin = C.load_lineage()
lin.index = lin.index.str.lower()

# Restrict to lines in the matrix
shared_idx = S.index.intersection(lin.index)
n_labelled = len(shared_idx)
n_unlabelled = n_lines - n_labelled
print(f"  Lines with lineage label : {n_labelled:,}")
print(f"  Lines without label      : {n_unlabelled:,}  (excluded from ARI/NMI)")
print(f"  Lineage source           : sample_info.parquet (20Q2) + "
      f"Model.csv (24Q4 supplement)")

# Build aligned arrays for labelled lines
S_sub = S.loc[shared_idx, shared_idx]
tissue_series = lin.loc[shared_idx, "lineage"].fillna("unknown")
unique_tissues = sorted(tissue_series.unique())
n_tissues = len(unique_tissues)
tissue_counts = tissue_series.value_counts()

print(f"\n  Distinct tissue labels   : {n_tissues}")
print(f"\n  {'Tissue':<35} {'Lines':>6}")
print(f"  {'-'*35} {'-'*6}")
small_tissues = []
for tissue, count in tissue_counts.items():
    flag = "  ← <10 lines" if count < SMALL_TISSUE_THRESHOLD else ""
    print(f"  {tissue:<35} {count:>6}{flag}")
    if count < SMALL_TISSUE_THRESHOLD:
        small_tissues.append((tissue, count))

if small_tissues:
    print(f"\n  Tissues with < {SMALL_TISSUE_THRESHOLD} lines "
          f"(noted, not dropped): {[t for t, _ in small_tissues]}")

# Integer tissue labels for metric computation
tissue_list = tissue_series.tolist()
tissue2int = {t: i for i, t in enumerate(unique_tissues)}
y_true = np.array([tissue2int[t] for t in tissue_list])

c1_summary = {
    "n_lines_matrix": n_lines,
    "n_labelled": n_labelled,
    "n_unlabelled": n_unlabelled,
    "n_tissues": n_tissues,
    "tissue_counts": tissue_counts.to_dict(),
    "small_tissues": small_tissues,
}
print("\nC1 complete — inputs confirmed.\n")


# ═══════════════════════════════════════════════════════════════════════ C2
banner("C2 — CLUSTERING")

from scipy.cluster.hierarchy import linkage, fcluster, dendrogram  # noqa: E402
from scipy.spatial.distance import squareform                       # noqa: E402
import networkx as nx                                               # noqa: E402
from networkx.algorithms.community import louvain_communities       # noqa: E402

# ── C2a: Hierarchical clustering ─────────────────────────────────────────────
print("\n[C2a] Hierarchical clustering …")
print(f"  Linkage method : {HC_LINKAGE}")
print("  Distance       : 1 − similarity  (sets diagonal 0, so squareform works)")

# Replace NaN diagonal with 0 for distance
D_vals = 1.0 - S_sub.values
np.fill_diagonal(D_vals, 0.0)

# Clip negative distances (similarity can be slightly below -1 due to float precision)
D_vals = np.clip(D_vals, 0.0, 2.0)

condensed = squareform(D_vals, checks=False)
Z = linkage(condensed, method=HC_LINKAGE)
print("  Linkage computed.")

# Try cluster counts: n_tissues ±30%, plus some round numbers
base_k = n_tissues
hc_k_values = sorted(set([
    max(5, round(base_k * 0.70)),
    max(5, round(base_k * 0.85)),
    base_k,
    round(base_k * 1.15),
    round(base_k * 1.30),
    round(base_k * 1.50),
]))
print(f"  Cluster counts to try    : {hc_k_values}")

hc_assignments = {}
for k in hc_k_values:
    labels = fcluster(Z, k, criterion="maxclust")
    hc_assignments[k] = labels
    print(f"  k={k:3d}  actual clusters: "
          f"{len(np.unique(labels)):3d}")

# ── C2b: Louvain community detection ─────────────────────────────────────────
print(f"\n[C2b] Louvain community detection …")
print(f"  k values: {LOUVAIN_K_VALUES}")
print("  Graph    : k-NN built from similarity matrix (directed -> undirected)")

louvain_assignments = {}
idx_arr = np.array(shared_idx.tolist())
sim_vals = S_sub.values.copy()

for k in LOUVAIN_K_VALUES:
    print(f"  Building k={k} graph …", end="", flush=True)
    G = nx.Graph()
    G.add_nodes_from(range(n_labelled))
    for i in range(n_labelled):
        row = sim_vals[i].copy()
        row[i] = -np.inf
        top_k = np.argsort(-row)[:k]
        for j in top_k:
            if not np.isnan(sim_vals[i, j]):
                G.add_edge(i, j, weight=float(sim_vals[i, j]))

    comms = louvain_communities(G, seed=SEED, weight="weight")
    labels = np.zeros(n_labelled, dtype=int)
    for c_idx, comm in enumerate(comms):
        for node in comm:
            labels[node] = c_idx
    n_comms = len(comms)
    louvain_assignments[k] = labels
    print(f"  {n_comms} communities")

print("C2 complete.\n")


# ═══════════════════════════════════════════════════════════════════════ C3
banner("C3 — TISSUE CONCORDANCE")

from sklearn.metrics import adjusted_rand_score as ARI        # noqa: E402
from sklearn.metrics import normalized_mutual_info_score as NMI  # noqa: E402

rng = np.random.default_rng(SEED)


def null_distribution(y_pred: np.ndarray, y_true: np.ndarray,
                       n: int = N_NULLS) -> tuple[np.ndarray, np.ndarray]:
    """Shuffle y_true labels N times, compute ARI and NMI each time."""
    aris, nmis = [], []
    for _ in range(n):
        y_shuf = rng.permutation(y_true)
        aris.append(ARI(y_shuf, y_pred))
        nmis.append(NMI(y_shuf, y_pred, average_method="arithmetic"))
    return np.array(aris), np.array(nmis)


def score_clustering(name: str, y_pred: np.ndarray) -> dict:
    ari = ARI(y_true, y_pred)
    nmi = NMI(y_true, y_pred, average_method="arithmetic")
    null_aris, null_nmis = null_distribution(y_pred, y_true, N_NULLS)
    ari_pct = float(np.mean(null_aris < ari))
    nmi_pct = float(np.mean(null_nmis < nmi))
    result = {
        "name": name,
        "n_clusters": int(len(np.unique(y_pred))),
        "ari": float(ari),
        "nmi": float(nmi),
        "null_ari_mean": float(null_aris.mean()),
        "null_ari_sd": float(null_aris.std()),
        "null_nmi_mean": float(null_nmis.mean()),
        "null_nmi_sd": float(null_nmis.std()),
        "ari_pct_above_null": ari_pct,
        "nmi_pct_above_null": nmi_pct,
    }
    print(f"\n  {name}")
    print(f"    ARI  = {ari:.4f}  (null mean {null_aris.mean():.4f} ± "
          f"{null_aris.std():.4f}, above null {ari_pct:.1%})")
    print(f"    NMI  = {nmi:.4f}  (null mean {null_nmis.mean():.4f} ± "
          f"{null_nmis.std():.4f}, above null {nmi_pct:.1%})")
    return result


def cluster_composition(y_pred: np.ndarray) -> list[dict]:
    """Per-cluster tissue breakdown, flagging cross-tissue and pure clusters."""
    comp = []
    for c in sorted(np.unique(y_pred)):
        mask = y_pred == c
        n = int(mask.sum())
        counts = Counter(tissue_list[i] for i in range(len(tissue_list)) if mask[i])
        top_tissue, top_n = counts.most_common(1)[0]
        top_frac = top_n / n
        breakdown = {t: cnt / n for t, cnt in counts.most_common()}
        cross_tissue = top_frac < CROSS_TISSUE_MAJORITY
        nearly_pure = top_frac >= PURE_TISSUE_THRESHOLD
        comp.append({
            "cluster": int(c),
            "n_lines": n,
            "top_tissue": top_tissue,
            "top_tissue_frac": float(top_frac),
            "cross_tissue": cross_tissue,
            "nearly_pure": nearly_pure,
            "breakdown": {t: round(f, 4) for t, f in breakdown.items()},
        })
    return comp


print("\nScoring all clusterings against tissue labels …")
print(f"({N_NULLS} null permutations per clustering — this takes a few minutes)\n")

all_scores = []
all_compositions = {}

# Hierarchical
print("--- Hierarchical clustering ---")
for k in hc_k_values:
    name = f"HC_avg_k{k}"
    sc = score_clustering(name, hc_assignments[k])
    all_scores.append(sc)
    all_compositions[name] = cluster_composition(hc_assignments[k])

# Louvain
print("\n--- Louvain community detection ---")
for k in LOUVAIN_K_VALUES:
    name = f"Louvain_k{k}"
    sc = score_clustering(name, louvain_assignments[k])
    all_scores.append(sc)
    all_compositions[name] = cluster_composition(louvain_assignments[k])

print("\nC3 — per-cluster composition (flagged cross-tissue clusters)")
cross_tissue_flagged = {}
for method_name, comp in all_compositions.items():
    ct = [cl for cl in comp if cl["cross_tissue"] and cl["n_lines"] >= 5]
    pure = [cl for cl in comp if cl["nearly_pure"]]
    cross_tissue_flagged[method_name] = ct
    print(f"\n  {method_name}:")
    print(f"    Total clusters      : {len(comp)}")
    print(f"    Nearly-pure (≥90%)  : {len(pure)}")
    print(f"    Cross-tissue (<50%) : {len(ct)}")
    if ct:
        print(f"    Cross-tissue detail:")
        for cl in ct[:10]:
            top_items = list(cl["breakdown"].items())[:4]
            breakdown_str = ", ".join(f"{t}: {f:.0%}" for t, f in top_items)
            print(f"      cluster {cl['cluster']:3d}  n={cl['n_lines']:4d}  "
                  f"top={cl['top_tissue']} ({cl['top_tissue_frac']:.0%})  "
                  f"[{breakdown_str}]")

print("\nC3 complete.\n")


# ═══════════════════════════════════════════════════════════════════════ C4
banner("C4 — SANITY CHECKS ON FLAGGED CROSS-TISSUE CLUSTERS")

print("\nChecking flagged cross-tissue clusters for artefacts …")

# Load expression data to check top-variance genes in flagged clusters
expr_path = ROOT / "cell_similarity" / "outputs" / "expression_24q4.parquet"

print(f"\n  Expression file: {expr_path}")
if expr_path.exists():
    expr = pd.read_parquet(expr_path)
    expr.index = expr.index.str.lower()
    print(f"  Loaded: {expr.shape}")
    expr_available = True
else:
    print("  expression_24q4.parquet not found — will skip gene-driver check")
    expr_available = False

# Y-linked genes (from the sex-guard work) — lowercase to match expression columns
Y_LINKED_GENES = {
    "ensg00000129824",  # RPS4Y1
    "ensg00000067048",  # DDX3Y
    "ensg00000012817",  # KDM5D
    "ensg00000114374",  # USP9Y
    "ensg00000183878",  # UTY
    "ensg00000154620",  # TMSB4Y
    "ensg00000229807",  # XIST (X-inactivation sentinel)
}

c4_results = {}
idx_list = shared_idx.tolist()

for method_name, flagged in cross_tissue_flagged.items():
    if not flagged:
        c4_results[method_name] = []
        continue

    method_c4 = []
    if method_name.startswith("HC"):
        pred_labels = hc_assignments[int(method_name.split("_k")[1])]
    else:
        pred_labels = louvain_assignments[int(method_name.split("_k")[1])]

    for cl in flagged:
        cid = cl["cluster"]
        mask = pred_labels == cid
        cluster_lines = [idx_list[i] for i in range(len(idx_list)) if mask[i]]
        result = {
            "cluster": cid,
            "n_lines": cl["n_lines"],
            "breakdown": cl["breakdown"],
            "checks": {},
        }

        # Check 1: small cluster (< 5 lines — already filtered, but note if <10)
        if cl["n_lines"] < 10:
            result["checks"]["small_cluster"] = True
            result["checks"]["small_cluster_note"] = (
                f"Only {cl['n_lines']} lines — low confidence cluster")
        else:
            result["checks"]["small_cluster"] = False

        # Check 2: gene-driver check (top variance genes in cluster vs full panel)
        if expr_available:
            cl_expr = expr.reindex(cluster_lines).dropna(how="all")
            if len(cl_expr) >= 3:
                cl_var = cl_expr.var(axis=0)
                panel_var = expr.var(axis=0)
                # Top 20 genes by variance excess in cluster vs panel
                excess = (cl_var / (panel_var + 1e-6)).sort_values(ascending=False)
                top20 = excess.head(20).index.tolist()
                y_hits = [g for g in top20 if g in Y_LINKED_GENES]
                result["checks"]["top20_excess_var_genes"] = top20[:10]
                result["checks"]["y_linked_in_top20"] = y_hits
                result["checks"]["y_linked_driven"] = len(y_hits) > 0
            else:
                result["checks"]["expr_coverage"] = "< 3 lines in expr data"
        else:
            result["checks"]["expr_check"] = "skipped — no expression_24q4.parquet"

        method_c4.append(result)
        status = "⚠ Y-LINKED DRIVEN" if result["checks"].get("y_linked_driven") else "OK"
        print(f"\n  [{method_name}] cluster {cid}  n={cl['n_lines']}  {status}")
        print(f"    tissues: {', '.join(f'{t}:{f:.0%}' for t, f in list(cl['breakdown'].items())[:5])}")
        if result["checks"].get("y_linked_driven"):
            print(f"    Y-linked in top variance: {result['checks']['y_linked_in_top20']}")

    c4_results[method_name] = method_c4

print("\nC4 complete.\n")


# ═══════════════════════════════════════════════════════════════════════
# Save raw artefacts
banner("SAVING ARTEFACTS")

# Save ARI/NMI table
scores_df = pd.DataFrame(all_scores)
scores_df.to_csv(RESULTS / "concordance_scores.csv", index=False)
print(f"  concordance_scores.csv  ({len(scores_df)} rows)")

# Save cluster assignments
assign_rows = []
for method_name, comp in all_compositions.items():
    if method_name.startswith("HC"):
        k = int(method_name.split("_k")[1])
        pred = hc_assignments[k]
    else:
        k = int(method_name.split("_k")[1])
        pred = louvain_assignments[k]
    for i, model_id in enumerate(idx_list):
        assign_rows.append({
            "model_id": model_id,
            "method": method_name,
            "cluster": int(pred[i]),
            "tissue": tissue_list[i],
        })
assign_df = pd.DataFrame(assign_rows)
assign_df.to_csv(RESULTS / "cluster_assignments.csv", index=False)
print(f"  cluster_assignments.csv ({len(assign_df):,} rows)")

# Save compositions
comp_rows = []
for method_name, comp in all_compositions.items():
    for cl in comp:
        comp_rows.append({
            "method": method_name,
            "cluster": cl["cluster"],
            "n_lines": cl["n_lines"],
            "top_tissue": cl["top_tissue"],
            "top_tissue_frac": cl["top_tissue_frac"],
            "cross_tissue": cl["cross_tissue"],
            "nearly_pure": cl["nearly_pure"],
            "breakdown_json": json.dumps(cl["breakdown"]),
        })
comp_df = pd.DataFrame(comp_rows)
comp_df.to_csv(RESULTS / "cluster_compositions.csv", index=False)
print(f"  cluster_compositions.csv ({len(comp_df)} rows)")

# Save full results JSON
full_results = {
    "c1": c1_summary,
    "c3_scores": all_scores,
    "c3_compositions": {k: v for k, v in all_compositions.items()},
    "c4_sanity": c4_results,
}
(RESULTS / "full_results.json").write_text(
    json.dumps(full_results, indent=2, default=str), encoding="utf-8")
print("  full_results.json")

print("\nDone. Results in cell_similarity/clustering/results/")
