"""
cell_similarity/clustering/followup_analysis.py
------------------------------------------------
CoWork v8.1 follow-up: F1 (bootstrap), F2 (duplicate check), F3 (cluster #2 vet).

Bootstrap criterion (pre-stated, per spec §2.3):
  For each reference cluster R, a bootstrap replicate is counted as "reproducing"
  R if there exists a replicate cluster C' such that:
    |R_obs ∩ C'| / |R_obs| >= 0.50
  where R_obs = R ∩ {lines that appear in this bootstrap replicate}.
  Threshold = 0.50 (chosen before results were seen).

Methods:
  - Louvain_k20  (primary reference, same parameters as v8.0)
  - HC_avg_k34   (secondary, HC at k=n_tissues)
  200 bootstrap replicates each.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent.parent
HERE = Path(__file__).resolve().parent
RESULTS = HERE / "results"

sys.path.insert(0, str(ROOT / "cell_similarity"))
import common as C  # noqa: E402

SEED = 42
N_BOOT = 200
REPRODUCE_THRESHOLD = 0.50
LOUVAIN_K = 20
HC_K = 34

rng = np.random.default_rng(SEED)


def banner(s: str):
    print("=" * 70)
    print(s)
    print("=" * 70)


# ─────────────────────────────────────────────────────────── load reference
banner("LOADING REFERENCE DATA")

S = pd.read_parquet(ROOT / "cell_similarity" / "outputs" / "rna_similarity.parquet")
lin = C.load_lineage()
lin.index = lin.index.str.lower()
shared_idx = S.index.intersection(lin.index).tolist()
tissue_map = lin.loc[shared_idx, "lineage"].fillna("unknown").to_dict()
n = len(shared_idx)
line2i = {m: i for i, m in enumerate(shared_idx)}
i2line = {i: m for i, m in enumerate(shared_idx)}
sim_vals = S.loc[shared_idx, shared_idx].values.copy()
print(f"  Lines: {n}, Tissues: {len(set(tissue_map.values()))}")

# Reference clusters (from v8.0 Louvain_k20 run)
assign_df = pd.read_csv(RESULTS / "cluster_assignments.csv")
lk20 = assign_df[assign_df["method"] == "Louvain_k20"].set_index("model_id")

# Define the four reference cluster line sets (by Louvain_k20 cluster ID)
# Hematological (cluster 13 in Louvain_k20): blood+lymphocyte+plasma_cell
# Gynecological epithelial (cluster 9): ovary+kidney+uterus
# Squamous carcinoma (cluster 4): upper_aerodigestive+esophagus+urinary_tract
# Mesenchymal/neuroectodermal (cluster 13): CNS+fibroblast+soft_tissue
# First: identify which Louvain_k20 clusters correspond to the four groups
print("\nLouvain_k20 cluster compositions (top 3 tissues per cluster):")
for cid in sorted(lk20["cluster"].unique()):
    sub = lk20[lk20["cluster"] == cid]
    vc = sub["tissue"].value_counts()
    top3 = ", ".join(f"{t}:{c}" for t, c in vc.head(3).items())
    print(f"  cluster {cid:2d}  n={len(sub):4d}  {top3}")

# Hard-code the reference cluster IDs based on v8.0 Louvain_k20 output
REF_CLUSTERS = {
    "hematological":         {"louvain_k20_id": None},
    "gynecological":         {"louvain_k20_id": 9},
    "squamous":              {"louvain_k20_id": None},
    "mesenchymal":           {"louvain_k20_id": None},
}


def find_cluster_by_tissue(top_tissues: list[str]) -> int:
    """Find Louvain_k20 cluster whose top tissue is in top_tissues."""
    for cid in sorted(lk20["cluster"].unique()):
        sub = lk20[lk20["cluster"] == cid]
        top = sub["tissue"].value_counts().idxmax()
        if top in top_tissues:
            return int(cid)
    return -1


heme_id = find_cluster_by_tissue(["blood", "lymphocyte"])
gyno_id = 9  # confirmed above
squa_id = find_cluster_by_tissue(["upper_aerodigestive"])
meso_id = find_cluster_by_tissue(["central_nervous_system"])

ref_sets = {
    "hematological":  set(lk20[lk20["cluster"] == heme_id].index.tolist()),
    "gynecological":  set(lk20[lk20["cluster"] == gyno_id].index.tolist()),
    "squamous":       set(lk20[lk20["cluster"] == squa_id].index.tolist()),
    "mesenchymal":    set(lk20[lk20["cluster"] == meso_id].index.tolist()),
}

print("\nReference cluster sizes:")
for name, lines in ref_sets.items():
    top_tissues = pd.Series([tissue_map.get(l, "?") for l in lines]).value_counts()
    top3 = ", ".join(f"{t}:{c}" for t, c in top_tissues.head(3).items())
    print(f"  {name:<20} n={len(lines):4d}  {top3}")


# ═══════════════════════════════════════════════════════════════════════ F1
banner("F1 — BOOTSTRAP STABILITY (200 replicates × 2 methods)")

from scipy.cluster.hierarchy import linkage, fcluster  # noqa: E402
from scipy.spatial.distance import squareform          # noqa: E402
import networkx as nx                                   # noqa: E402
from networkx.algorithms.community import louvain_communities  # noqa: E402


def check_reproduction(replicate_labels: dict[str, int],
                        ref_set: set[str],
                        threshold: float = REPRODUCE_THRESHOLD) -> bool:
    """
    Returns True if any replicate cluster captures >= threshold fraction
    of the observable reference cluster members.
    replicate_labels: {model_id -> cluster_id} for the resampled unique lines.
    ref_set: set of model_ids in the reference cluster.
    """
    observable = ref_set & set(replicate_labels.keys())
    if len(observable) < 5:    # too few observable members to test
        return None            # mark as untestable, excluded from rate
    cluster_counts: dict[int, int] = {}
    for m in observable:
        c = replicate_labels[m]
        cluster_counts[c] = cluster_counts.get(c, 0) + 1
    best_overlap = max(cluster_counts.values())
    return (best_overlap / len(observable)) >= threshold


def run_louvain_on_submatrix(unique_lines: list[str], k: int,
                              seed: int) -> dict[str, int]:
    """Build kNN graph on a subset of lines, run Louvain, return {line: cluster}."""
    idx = [line2i[m] for m in unique_lines if m in line2i]
    sub = sim_vals[np.ix_(idx, idx)]
    m = len(unique_lines)
    G = nx.Graph()
    G.add_nodes_from(range(m))
    for i in range(m):
        row = sub[i].copy()
        row[i] = -np.inf
        top_k = np.argsort(-row)[:k]
        for j in top_k:
            if not np.isnan(sub[i, j]):
                G.add_edge(i, j, weight=float(sub[i, j]))
    comms = louvain_communities(G, seed=seed, weight="weight")
    labels_arr = np.zeros(m, dtype=int)
    for c_idx, comm in enumerate(comms):
        for node in comm:
            labels_arr[node] = c_idx
    return {unique_lines[i]: int(labels_arr[i]) for i in range(m)}


def run_hc_on_submatrix(unique_lines: list[str], k: int) -> dict[str, int]:
    """Build HC clustering on a subset of lines."""
    idx = [line2i[m] for m in unique_lines if m in line2i]
    sub = sim_vals[np.ix_(idx, idx)]
    D = 1.0 - sub
    np.fill_diagonal(D, 0.0)
    D = np.clip(D, 0.0, 2.0)
    condensed = squareform(D, checks=False)
    Z = linkage(condensed, method="average")
    labels_arr = fcluster(Z, k, criterion="maxclust")
    return {unique_lines[i]: int(labels_arr[i]) for i in range(len(unique_lines))}


print(f"\nBootstrap criterion: reproduced if max_C'(|R_obs ∩ C'| / |R_obs|) >= {REPRODUCE_THRESHOLD}")
print(f"Replicates: {N_BOOT} per method")
print(f"Methods: Louvain k={LOUVAIN_K}, HC average k={HC_K}\n")

# Track reproduction counts
louvain_reps = {name: 0 for name in ref_sets}
hc_reps = {name: 0 for name in ref_sets}
louvain_valid = {name: 0 for name in ref_sets}  # replicates where cluster is testable
hc_valid = {name: 0 for name in ref_sets}

for boot_idx in range(N_BOOT):
    if (boot_idx + 1) % 20 == 0:
        print(f"  Bootstrap {boot_idx + 1}/{N_BOOT} …")

    # Sample with replacement
    sampled = rng.choice(shared_idx, size=n, replace=True)
    unique_lines = list(dict.fromkeys(sampled))  # ordered unique, preserving first occurrence

    # Louvain
    louvain_labels = run_louvain_on_submatrix(
        unique_lines, k=LOUVAIN_K, seed=int(rng.integers(0, 2**31)))
    for name, ref_set in ref_sets.items():
        result = check_reproduction(louvain_labels, ref_set)
        if result is not None:
            louvain_valid[name] += 1
            if result:
                louvain_reps[name] += 1

    # HC
    hc_labels = run_hc_on_submatrix(unique_lines, k=HC_K)
    for name, ref_set in ref_sets.items():
        result = check_reproduction(hc_labels, ref_set)
        if result is not None:
            hc_valid[name] += 1
            if result:
                hc_reps[name] += 1

print("\nF1 RESULTS — Bootstrap reproduction rates")
print(f"  (criterion: >= {REPRODUCE_THRESHOLD:.0%} of observable members in one replicate cluster)\n")
print(f"  {'Cluster':<22} {'Louvain k=20':>14} {'HC k=34':>14}")
print(f"  {'-'*22} {'-'*14} {'-'*14}")

f1_results = {}
escalate = False
for name in ref_sets:
    lv_rate = louvain_reps[name] / louvain_valid[name] if louvain_valid[name] else float("nan")
    hc_rate = hc_reps[name] / hc_valid[name] if hc_valid[name] else float("nan")
    lv_str = f"{lv_rate:.1%} ({louvain_reps[name]}/{louvain_valid[name]})"
    hc_str = f"{hc_rate:.1%} ({hc_reps[name]}/{hc_valid[name]})"
    print(f"  {name:<22} {lv_str:>14} {hc_str:>14}")
    f1_results[name] = {"louvain_rate": lv_rate, "louvain_n_valid": louvain_valid[name],
                        "hc_rate": hc_rate, "hc_n_valid": hc_valid[name]}
    if name == "gynecological" and (lv_rate < 0.50 or hc_rate < 0.50):
        escalate = True

if escalate:
    print("\n⚠ ESCALATE: gynecological cluster bootstrap rate < 50% in at least one method.")
    print("  Review F1 results before finalising F4 characterisation.")
else:
    print("\n  No ESCALATE condition triggered.")

print()


# ═══════════════════════════════════════════════════════════════════════ F2
banner("F2 — DUPLICATE / TECHNICAL REPLICATE CHECK")

print("\nChecking DepMap Model.csv for contamination/problematic flags …")

import importlib.util as _ilu
_ga_spec = _ilu.spec_from_file_location("gate_audit_common",
                                        ROOT / "gate_audit" / "common.py")
GA = _ilu.module_from_spec(_ga_spec)
_ga_spec.loader.exec_module(GA)

model_csv = GA.DEPMAP_24Q4 / "Model.csv"

if model_csv.exists():
    # Read all columns to find any duplicate/contamination flags
    model_df = pd.read_csv(model_csv, low_memory=False)
    model_df["model_id"] = model_df["ModelID"].astype(str).str.lower()
    print(f"  Model.csv columns: {[c for c in model_df.columns if any(k in c.lower() for k in ['contam', 'duplicate', 'problem', 'flag', 'exclude', 'curated', 'blacklist', 'quality'])]}")
    print(f"  All columns: {list(model_df.columns)}")
else:
    print("  Model.csv not found — checking cleaned_track_data")
    model_df = None

# Check for duplicate ModelID within the file (identical DepMap IDs)
f2_results = {}
for cluster_name, ref_lines in ref_sets.items():
    ref_list = sorted(ref_lines)
    issues = []

    if model_df is not None:
        # Any lines that appear more than once in Model.csv (technical replicates)
        sub_model = model_df[model_df["model_id"].isin(ref_lines)]
        dup_models = sub_model[sub_model.duplicated("ModelID", keep=False)]["ModelID"].unique()
        if len(dup_models) > 0:
            issues.append(f"duplicate ModelID in Model.csv: {list(dup_models)}")

        # Check for any "Problematic" or similar column
        prob_cols = [c for c in model_df.columns
                     if any(k in c.lower() for k in ["problem", "exclude", "contam", "flag"])]
        for col in prob_cols:
            flagged = sub_model[sub_model[col].notna() & (sub_model[col] != "")]["model_id"].tolist()
            if flagged:
                issues.append(f"{col} flagged: {flagged}")

    # Check for model_id collisions within cluster (same line appearing twice)
    if len(ref_list) != len(set(ref_list)):
        issues.append("model_id collision within cluster — data error")

    f2_results[cluster_name] = {
        "n_lines": len(ref_lines),
        "issues": issues,
        "clean": len(issues) == 0,
    }
    print(f"\n  {cluster_name} (n={len(ref_lines)}): {'CLEAN' if not issues else 'ISSUES: ' + str(issues)}")

print()


# ═══════════════════════════════════════════════════════════════════════ F3
banner("F3 — VET CLUSTER #2 (GYNECOLOGICAL EPITHELIAL)")

# F3.1: Exact composition
gyno_lines = sorted(ref_sets["gynecological"])
gyno_tissues = [tissue_map.get(m, "unknown") for m in gyno_lines]
tissue_counts = pd.Series(gyno_tissues).value_counts()

print(f"\nF3.1 — Exact composition")
print(f"  Total lines: {len(gyno_lines)}")
print(f"  {'Tissue':<25} {'Count':>6} {'%':>7}")
print(f"  {'-'*25} {'-'*6} {'-'*7}")
for tissue, count in tissue_counts.items():
    print(f"  {tissue:<25} {count:>6} {count/len(gyno_lines):>7.1%}")

core3 = tissue_counts[tissue_counts.index.isin(["ovary", "kidney", "uterus"])].sum()
print(f"\n  Core 3 (ovary+kidney+uterus): {core3} lines ({core3/len(gyno_lines):.1%})")
print(f"  Is it roughly even three-way? ovary/kidney ratio = "
      f"{tissue_counts.get('ovary', 0)}/{tissue_counts.get('kidney', 0)} = "
      f"{tissue_counts.get('ovary', 0)/max(tissue_counts.get('kidney', 0), 1):.2f}")

print(f"\n  Exact line IDs: {gyno_lines}")

# F3.2: Check three-way balance
print(f"\nF3.2 — Three-way balance assessment")
ovary_n = tissue_counts.get("ovary", 0)
kidney_n = tissue_counts.get("kidney", 0)
uterus_n = tissue_counts.get("uterus", 0)
total_core = ovary_n + kidney_n + uterus_n
if total_core > 0:
    print(f"  Within core 3: ovary {ovary_n/total_core:.1%}, "
          f"kidney {kidney_n/total_core:.1%}, "
          f"uterus {uterus_n/total_core:.1%}")
    dom = tissue_counts.head(1).index[0]
    dom_n = tissue_counts.head(1).iloc[0]
    print(f"  Largest single tissue: {dom} ({dom_n} lines, {dom_n/len(gyno_lines):.1%} of cluster)")
    if dom_n / len(gyno_lines) > 0.50:
        print("  → DOMINATED by one tissue (>50%). Not a three-way mix.")
    else:
        print("  → No single tissue dominates. Three-way mix confirmed.")

# F3.3: Defining genes — mean expression cluster vs rest
print(f"\nF3.3 — Cluster-defining genes")
expr_path = ROOT / "cell_similarity" / "outputs" / "expression_24q4.parquet"
expr = pd.read_parquet(expr_path)
expr.index = expr.index.str.lower()

gyno_set = set(gyno_lines)
rest_set = set(shared_idx) - gyno_set

gyno_expr = expr.reindex(sorted(gyno_set & set(expr.index))).dropna(how="all")
rest_expr = expr.reindex(sorted(rest_set & set(expr.index))).dropna(how="all")

print(f"  Cluster lines with expression: {len(gyno_expr)}")
print(f"  Rest-of-panel lines: {len(rest_expr)}")

# Mean difference (cluster - rest), normalised by pooled std
gyno_mean = gyno_expr.mean(axis=0)
rest_mean = rest_expr.mean(axis=0)
pooled_std = np.sqrt(
    (gyno_expr.var(axis=0) * (len(gyno_expr) - 1) +
     rest_expr.var(axis=0) * (len(rest_expr) - 1)) /
    (len(gyno_expr) + len(rest_expr) - 2)
)
pooled_std = pooled_std.replace(0, np.nan)
cohen_d = (gyno_mean - rest_mean) / pooled_std

# Top genes by |Cohen's d|
top_up = cohen_d.dropna().sort_values(ascending=False).head(30)
top_down = cohen_d.dropna().sort_values(ascending=True).head(30)

# Try to get gene symbols from the pipeline lookup
gene_lookup_path = ROOT / "final_pipeline" / "reference" / "gene_lookup.parquet"
ensg2sym = {}
if gene_lookup_path.exists():
    gl = pd.read_parquet(gene_lookup_path)
    if "ensg_id" in gl.columns and "gene_name" in gl.columns:
        gl["ensg_id_lower"] = gl["ensg_id"].str.lower()
        ensg2sym = gl.set_index("ensg_id_lower")["gene_name"].to_dict()
    elif "ensembl_id" in gl.columns and "symbol" in gl.columns:
        gl["ensembl_id_lower"] = gl["ensembl_id"].str.lower()
        ensg2sym = gl.set_index("ensembl_id_lower")["symbol"].to_dict()
    print(f"  Gene lookup loaded: {len(ensg2sym)} mappings")
else:
    print("  Gene lookup not found — reporting ENSG IDs only")

def ensg_to_sym(ensg_id: str) -> str:
    return ensg2sym.get(ensg_id.lower(), ensg_id)

print(f"\n  Top 20 UPREGULATED in cluster #2 (vs rest of panel, by Cohen's d):")
print(f"  {'Gene':<20} {'ENSG':<20} {'Cohen d':>10}")
for ensg, d in top_up.head(20).items():
    sym = ensg_to_sym(ensg)
    print(f"  {sym:<20} {ensg:<20} {d:>10.3f}")

print(f"\n  Top 20 DOWNREGULATED in cluster #2 (vs rest of panel):")
print(f"  {'Gene':<20} {'ENSG':<20} {'Cohen d':>10}")
for ensg, d in top_down.head(20).items():
    sym = ensg_to_sym(ensg)
    print(f"  {sym:<20} {ensg:<20} {d:>10.3f}")

f3_results = {
    "n_lines": len(gyno_lines),
    "tissue_counts": tissue_counts.to_dict(),
    "core3_n": int(core3),
    "core3_frac": float(core3 / len(gyno_lines)),
    "top_up_genes": {ensg_to_sym(e): float(d) for e, d in top_up.head(20).items()},
    "top_down_genes": {ensg_to_sym(e): float(d) for e, d in top_down.head(20).items()},
}

print(f"\nF3.4 — Literature context note")
print("  NOTE: Literature cross-check cannot be automated (no web access).")
print("  Key known facts for the examiner:")
print("  - Ovarian clear cell carcinoma (OCCC) and renal clear cell carcinoma (RCC)")
print("    share clear cell histology and VHL/HIF pathway activity in a subset.")
print("  - Endometrioid histology is shared by ovary and uterus (Mullerian origin).")
print("  - TCGA PanCancer Atlas (2018) identified a 'clear cell' cluster spanning")
print("    kidney and ovary in multi-omics analyses.")
print("  - Kidney does NOT share Mullerian/paramesonephric developmental origin")
print("    with ovary and uterus — this is the puzzle this clustering raises.")
print("  - Whether this grouping is driven by clear-cell lines specifically (a")
print("    subset of the ovary and kidney lineages) or the full lineage is")
print("    unresolvable from lineage labels alone and would require histology metadata.")


# ═══════════════════════════════════════════════════════════════════════
banner("SAVING FOLLOW-UP ARTEFACTS")

followup = {
    "criterion": {
        "description": "reproduced if max_C'(|R_obs ∩ C'| / |R_obs|) >= threshold",
        "threshold": REPRODUCE_THRESHOLD,
        "min_observable": 5,
        "n_bootstrap": N_BOOT,
    },
    "f1_bootstrap": f1_results,
    "f2_duplicates": f2_results,
    "f3_cluster2": f3_results,
    "escalate": escalate,
}
out_path = RESULTS / "followup_results.json"
out_path.write_text(json.dumps(followup, indent=2, default=str), encoding="utf-8")
print(f"  followup_results.json")

# Save defining genes as CSV
gene_rows = []
for ensg, d in top_up.head(20).items():
    gene_rows.append({"ensg": ensg, "symbol": ensg_to_sym(ensg), "cohen_d": float(d), "direction": "up"})
for ensg, d in top_down.head(20).items():
    gene_rows.append({"ensg": ensg, "symbol": ensg_to_sym(ensg), "cohen_d": float(d), "direction": "down"})
pd.DataFrame(gene_rows).to_csv(RESULTS / "cluster2_defining_genes.csv", index=False)
print(f"  cluster2_defining_genes.csv")

print("\nDone. Run complete.")
print(f"  Escalate flag: {escalate}")
