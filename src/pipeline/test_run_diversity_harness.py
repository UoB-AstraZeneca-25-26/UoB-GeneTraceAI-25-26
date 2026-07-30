"""
test_run_diversity_harness.py
--------------------------------
Tests the diversity-aware top-20 (best-line-per-cluster-then-fill) against
raw top-20 on hit@20, for the 9 genes with real GDSC ground truth. Also
checks whether expression-profile clustering just recovers metadata
grouping (if so, use the metadata labels -- cheaper, fully explainable).

Read-only. Does not touch core_score.parquet or any production file.
"""
import sys
import numpy as np
import pandas as pd
from pathlib import Path

sys.path.insert(0, "src/pipeline")
from explain_pair import load_data
from rank_cell_lines import rank_cell_lines
from cell_line_grouping import (
    group_top_n_by_metadata, cluster_by_expression_profile,
    compare_clusterings, diversify_top_n, disease_short_label,
)

REF = Path("reference")
DATA_CLEAN = Path("data/parquet/data_clean")
VAL = Path("validation/prepared")

GENES = {
    "EGFR": "ENSG00000146648", "MAP2K1": "ENSG00000169032", "MAP2K2": "ENSG00000126934",
    "CDK4": "ENSG00000135446", "CDK6": "ENSG00000105810", "BRAF": "ENSG00000157764",
    "MET": "ENSG00000105976", "ERBB2": "ENSG00000141736", "BCL2": "ENSG00000171791",
}

print("Loading data...")
pred, flags, variants, harm, gene_lookup, chronos_val = load_data()
curated = pd.read_parquet(VAL / "curated_validation_pairs.parquet")
curated["model_id"] = curated["model_id"].str.upper()
curated["ensg_id"] = curated["ensg_id"].str.upper()

# Full expression matrix for profile clustering (same construction as 02_core_score.ipynb cell 4)
print("Building expression matrix for profile clustering...")
prof = pd.read_parquet(REF / "depmap_profiles.parquet")
rna_prof = prof[prof["datatype"] == "rna"][["profileid", "modelid"]].copy()
rna_prof["model_id"] = rna_prof["modelid"].str.upper()
rna_prof = rna_prof.drop(columns=["modelid"])
expr_raw = pd.read_parquet(DATA_CLEAN / "depmap_expr_clean.parquet")
expr_raw.index.name = "profileid"
expr = expr_raw.reset_index().merge(rna_prof, on="profileid", how="inner")
expr = expr.drop(columns=["profileid"]).set_index("model_id")
expr.columns = expr.columns.str.upper()
if expr.index.duplicated().any():
    expr = expr.groupby(level=0).mean()
print(f"expr matrix: {expr.shape}")


def hit_at_k(model_ids, sensitive_model_ids, k=20):
    if len(sensitive_model_ids) == 0:
        return np.nan
    top_k = set(model_ids[:k])
    return len(top_k & sensitive_model_ids) / len(sensitive_model_ids)


print(f"\n{'Gene':<8}{'N pool':<8}{'ARI(meta,profile)':<20}{'Raw hit@20':<14}{'Diverse hit@20':<16}{'Delta'}")
results = {}
for gene, ensg in GENES.items():
    r = rank_cell_lines(ensg, pred, chronos_val)
    if "error" in r:
        print(f"{gene}: {r['error']}")
        continue
    ranking = r["ranking"]
    pool = ranking[:50]
    pool_ids = [row["model_id"] for row in pool]

    meta = group_top_n_by_metadata(ranking, n=50)
    profile_labels = cluster_by_expression_profile(
        pool_ids, expr, n_clusters=max(2, len(set(meta["model_id_to_group"].values())))
    )
    cmp = compare_clusterings(meta["model_id_to_group"], profile_labels)

    diversified = diversify_top_n(ranking, group_key="diseases", pool_size=50, k=20)
    div_ids = [row["model_id"] for row in diversified]
    raw_ids = [row["model_id"] for row in ranking]

    cur_gene = curated[curated.ensg_id == ensg]
    sens = set(cur_gene[cur_gene.ground_truth < 0]["model_id"])

    h_raw = hit_at_k(raw_ids, sens, k=20)
    h_div = hit_at_k(div_ids, sens, k=20)

    results[gene] = {
        "n_sensitive": len(sens), "ari": cmp["ari"], "hit_raw": h_raw, "hit_diverse": h_div,
        "delta": (h_div - h_raw) if not (np.isnan(h_raw) or np.isnan(h_div)) else np.nan,
        "top20_group_summary": group_top_n_by_metadata(ranking, n=20)["summary"],
        "diversified_group_summary": ", ".join(
            f"{v}" for v in pd.Series([disease_short_label(row.get("diseases")) for row in diversified])
            .value_counts().index.tolist()
        ),
    }
    ari_str = f"{cmp['ari']:.3f}" if cmp["ari"] is not None else "n/a"
    print(f"{gene:<8}{len(pool):<8}{ari_str:<20}{h_raw:<14.4f}{h_div:<16.4f}{results[gene]['delta']:+.4f}")

print("\n" + "=" * 70)
print("Top-20 composition (metadata grouping) -- raw vs diversified")
print("=" * 70)
for gene, r in results.items():
    print(f"\n{gene}:")
    print(f"  raw top-20:         {r['top20_group_summary']}")
    print(f"  diversified top-20: {r['diversified_group_summary']}")

import json
with open("src/pipeline/outputs/test_run_diversity_harness_results.json", "w") as f:
    json.dump(results, f, indent=2, default=float)
print(f"\nWritten: src/pipeline/outputs/test_run_diversity_harness_results.json")
