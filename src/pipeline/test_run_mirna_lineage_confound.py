"""
test_run_mirna_lineage_confound.py
-------------------------------------
Same diagnostic as test_run_metabolomics_lineage_confound.py, applied to
miRNA (see DECISIONS.md 2026-07-30 entry). Read-only, does not touch
core_score.parquet or any production file.

Step 1: per-line mean z-score across all 734 miRNAs (z-scored per miRNA
        across cell lines first, then averaged per line).
Step 2: ANOVA against lineage.
Step 3: add the global scalar as a constant term to the same GDSC harness
        (BCL2, MET), report delta-rho with a paired bootstrap CI.
"""
import numpy as np
import pandas as pd
from pathlib import Path
from scipy.stats import spearmanr, f_oneway

CLEANED = Path("cleaned_track_data")
DATA_CLEAN = Path("data/parquet/data_clean")
OUT = Path("src/pipeline/outputs")
VAL = Path("validation/prepared")

GENE_SYMBOL_TO_ENSG = {"BCL2": "ensg00000171791", "MET": "ensg00000105976"}

# ── Step 1: per-line mean z-score across all miRNAs ──────────────────────
mirna = pd.read_parquet(CLEANED / "mirna_model_level.parquet")
mirna["model_id"] = mirna["model_id"].str.lower()
n_mirnas = mirna["mirna_symbol"].nunique()
n_lines = mirna["model_id"].nunique()
print(f"Step 1: {n_mirnas} miRNAs, {n_lines} cell lines, {len(mirna):,} long rows")

# z-score each miRNA across cell lines (per mirna_symbol group), then
# average per line -- mirrors the metabolomics per-column z-score + row-mean
mirna["z"] = mirna.groupby("mirna_symbol")["expression_value"].transform(
    lambda s: (s - s.mean()) / s.std()
)
per_line = mirna.groupby("model_id")["z"].mean().rename("mirna_z").reset_index()
print(f"mirna_z: mean={per_line.mirna_z.mean():.4f} std={per_line.mirna_z.std():.4f} "
      f"range=[{per_line.mirna_z.min():.4f}, {per_line.mirna_z.max():.4f}]")

# ── Step 2: ANOVA against lineage ────────────────────────────────────────
si = pd.read_parquet(DATA_CLEAN / "sample_info_clean.parquet")
si["model_id"] = si["depmap_id"].str.lower()
per_line = per_line.merge(si[["model_id", "lineage"]], on="model_id", how="left")
missing_lineage = per_line["lineage"].isna().sum()
print(f"\nStep 2: {missing_lineage}/{len(per_line)} lines missing lineage (dropped)")
per_line_l = per_line.dropna(subset=["lineage"])

groups = [g["mirna_z"].values for _, g in per_line_l.groupby("lineage") if len(g) >= 2]
n_groups = len(groups)
fstat, pval = f_oneway(*groups)

grand_mean = per_line_l["mirna_z"].mean()
ss_total = ((per_line_l["mirna_z"] - grand_mean) ** 2).sum()
ss_between = per_line_l.groupby("lineage")["mirna_z"].apply(
    lambda g: len(g) * (g.mean() - grand_mean) ** 2
).sum()
r2 = ss_between / ss_total

print(f"Lineage groups (n>=2): {n_groups}")
print(f"ANOVA: F={fstat:.3f}  p={pval:.4g}")
print(f"R^2 (eta-squared) = {r2:.4f}  -> lineage explains {r2*100:.1f}% of variance in mirna_z")

lineage_means = per_line_l.groupby("lineage")["mirna_z"].agg(["mean", "count"]).sort_values("mean")
print("\nPer-lineage mean mirna_z (sorted):")
print(lineage_means.to_string())

# ── Step 3: constant-term blend, paired bootstrap CI on delta-rho ────────
print("\n" + "=" * 70)
print("STEP 3 -- constant-term blend into 2-layer core_score (BCL2, MET)")
print("=" * 70)

curated = pd.read_parquet(VAL / "curated_validation_pairs.parquet")
curated["model_id"] = curated["model_id"].str.lower()
curated["ensg_id"] = curated["ensg_id"].str.lower()

global_pct = per_line.set_index("model_id")["mirna_z"].rank(pct=True)

results_b = {}
for gene, ensg in GENE_SYMBOL_TO_ENSG.items():
    base = pd.read_parquet(OUT / "core_score.parquet")
    base["model_id"] = base["model_id"].str.lower()
    base = base[base.ensg_id == ensg][["model_id", "core_score"]].rename(columns={"core_score": "baseline"})
    base["global_pct"] = base["model_id"].map(global_pct)
    base = base.dropna(subset=["global_pct"])
    base["candidate"] = 0.5 * base["baseline"] + 0.5 * base["global_pct"]

    cur_gene = curated[curated.ensg_id == ensg][["model_id", "ground_truth"]]
    m = base.merge(cur_gene, on="model_id", how="inner").dropna()

    rho_base, _ = spearmanr(m["baseline"], m["ground_truth"])
    rho_cand, _ = spearmanr(m["candidate"], m["ground_truth"])

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
with open(OUT / "test_run_mirna_lineage_confound_results.json", "w") as f:
    json.dump({
        "anova_r2": r2, "anova_f": fstat, "anova_p": pval, "n_lineage_groups": n_groups,
        "lineage_means": lineage_means.reset_index().to_dict("records"),
        "harness_constant_term": results_b,
    }, f, indent=2, default=float)
print(f"\nWritten: {OUT / 'test_run_mirna_lineage_confound_results.json'}")
