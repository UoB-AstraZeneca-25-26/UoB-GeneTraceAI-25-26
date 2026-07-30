"""
test_run_metabolomics_lineage_confound.py
-------------------------------------------
Cheap diagnostic (per DECISIONS.md 2026-07-30 entry): does a per-cell-line
global metabolomics scalar just proxy lineage? Read-only, does not touch
core_score.parquet or any production file.

Step 1: melt metabolomics -> per-line mean z-score across 225 metabolites.
Step 2: ANOVA against lineage (cell_line_lookup / sample_info_clean).
Step 3: add the global scalar as a constant term to the same GDSC harness
        used in test_run_metabolomics_mirna.py (BCL2, MET -- only genes with
        real GDSC ground truth), report delta-rho with a PAIRED bootstrap CI.

NOT a candidate for core_score -- see DECISIONS.md for why a positive
result here would be a lineage-confound finding, not a validation.
"""
import numpy as np
import pandas as pd
from pathlib import Path
from scipy.stats import spearmanr, f_oneway
import sys
sys.path.insert(0, "src/pipeline")
from scoring_variants import combine

REF = Path("reference")
CLEANED = Path("cleaned_track_data")
DATA_CLEAN = Path("data/parquet/data_clean")
OUT = Path("src/pipeline/outputs")
VAL = Path("validation/prepared")

GENE_SYMBOL_TO_ENSG = {"BCL2": "ENSG00000171791", "MET": "ENSG00000105976"}

# ── Step 1: per-line mean z-score across all 225 metabolites ────────────
metab = pd.read_parquet(CLEANED / "metabolomics.parquet")
metab["model_id"] = metab["model_id"].str.upper()
value_cols = [c for c in metab.columns if c not in ("ccle_id", "model_id")]
print(f"Step 1: {len(value_cols)} metabolite columns, {len(metab)} cell lines")

z = metab[value_cols].apply(lambda col: (col - col.mean()) / col.std())
metab["metab_z"] = z.mean(axis=1)
per_line = metab[["model_id", "ccle_id", "metab_z"]].copy()
print(f"metab_z: mean={per_line.metab_z.mean():.4f} std={per_line.metab_z.std():.4f} "
      f"range=[{per_line.metab_z.min():.4f}, {per_line.metab_z.max():.4f}]")

# ── Step 2: ANOVA against lineage ────────────────────────────────────────
si = pd.read_parquet(DATA_CLEAN / "sample_info_clean.parquet")
si["model_id"] = si["depmap_id"].str.upper()
per_line = per_line.merge(si[["model_id", "lineage"]], on="model_id", how="left")
missing_lineage = per_line["lineage"].isna().sum()
print(f"\nStep 2: {missing_lineage}/{len(per_line)} lines missing lineage (dropped)")
per_line_l = per_line.dropna(subset=["lineage"])

groups = [g["metab_z"].values for _, g in per_line_l.groupby("lineage") if len(g) >= 2]
n_groups = len(groups)
fstat, pval = f_oneway(*groups)

# R^2 = SS_between / SS_total (eta-squared, standard ANOVA effect size)
grand_mean = per_line_l["metab_z"].mean()
ss_total = ((per_line_l["metab_z"] - grand_mean) ** 2).sum()
ss_between = per_line_l.groupby("lineage")["metab_z"].apply(
    lambda g: len(g) * (g.mean() - grand_mean) ** 2
).sum()
r2 = ss_between / ss_total

print(f"Lineage groups (n>=2): {n_groups}")
print(f"ANOVA: F={fstat:.3f}  p={pval:.4g}")
print(f"R^2 (eta-squared) = {r2:.4f}  -> lineage explains {r2*100:.1f}% of variance in metab_z")

lineage_means = per_line_l.groupby("lineage")["metab_z"].agg(["mean", "count"]).sort_values("mean")
print("\nPer-lineage mean metab_z (sorted):")
print(lineage_means.to_string())

# ── Step 3: add as constant term to the harness, paired bootstrap CI on delta-rho ──
print("\n" + "=" * 70)
print("STEP 3 -- constant-term blend into 2-layer core_score (BCL2, MET)")
print("=" * 70)

curated = pd.read_parquet(VAL / "curated_validation_pairs.parquet")
curated["model_id"] = curated["model_id"].str.upper()
curated["ensg_id"] = curated["ensg_id"].str.upper()

global_pct = per_line.set_index("model_id")["metab_z"].rank(pct=True)  # [0,1] per line, gene-agnostic

results_b = {}
for gene, ensg in GENE_SYMBOL_TO_ENSG.items():
    base = pd.read_parquet(OUT / "core_score.parquet")
    base["model_id"] = base["model_id"].str.upper()
    base = base[base.ensg_id == ensg][["model_id", "core_score"]].rename(columns={"core_score": "baseline"})
    base["global_pct"] = base["model_id"].map(global_pct)
    base = base.dropna(subset=["global_pct"])
    # constant-term blend: equal-weight with the existing 2-layer score
    base["candidate"] = 0.5 * base["baseline"] + 0.5 * base["global_pct"]

    cur_gene = curated[curated.ensg_id == ensg][["model_id", "ground_truth"]]
    m = base.merge(cur_gene, on="model_id", how="inner").dropna()

    rho_base, _ = spearmanr(m["baseline"], m["ground_truth"])
    rho_cand, _ = spearmanr(m["candidate"], m["ground_truth"])

    # paired bootstrap: resample cell lines (rows of m) together for both scores
    rng = np.random.default_rng(42)
    n = len(m)
    deltas = []
    for _ in range(2000):
        idx = rng.integers(0, n, n)
        s = m.iloc[idx]
        rb, _ = spearmanr(s["baseline"], s["ground_truth"])
        rc, _ = spearmanr(s["candidate"], s["ground_truth"])
        deltas.append(rc - rb)
    deltas = np.array(deltas)
    ci_lo, ci_hi = np.percentile(deltas, [2.5, 97.5])

    results_b[gene] = dict(n=n, rho_baseline=rho_base, rho_candidate=rho_cand,
                            delta=rho_cand - rho_base, ci_lo=ci_lo, ci_hi=ci_hi)
    print(f"\n{gene} ({ensg}): n={n}")
    print(f"  baseline rho  : {rho_base:.4f}")
    print(f"  candidate rho : {rho_cand:.4f}")
    print(f"  delta rho     : {rho_cand - rho_base:+.4f}  95% paired-bootstrap CI [{ci_lo:+.4f}, {ci_hi:+.4f}]")
    sig = "CI excludes 0 -> significant shift" if (ci_lo > 0 or ci_hi < 0) else "CI includes 0 -> not significant"
    print(f"  {sig}")

import json
with open(OUT / "test_run_metabolomics_lineage_confound_results.json", "w") as f:
    json.dump({
        "anova_r2": r2, "anova_f": fstat, "anova_p": pval, "n_lineage_groups": n_groups,
        "lineage_means": lineage_means.reset_index().to_dict("records"),
        "harness_constant_term": results_b,
    }, f, indent=2, default=float)
print(f"\nWritten: {OUT / 'test_run_metabolomics_lineage_confound_results.json'}")
