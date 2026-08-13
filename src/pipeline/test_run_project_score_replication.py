"""
test_run_project_score_replication.py
----------------------------------------
Read-only diagnostic. Writes no production file.

QUESTION
--------
The project's headline finding is that gene ABUNDANCE barely tracks functional
DEPENDENCY: rho(expression, Chronos essentiality) ~ 0.005. That rests on a single
essentiality dataset (Broad DepMap Chronos). Does it replicate on an INDEPENDENT
CRISPR screen?

Project Score (Sanger v2 + Broad 21Q2 combined) has been sitting unused in
data/GDSC/ -- 17,649 genes x 1,109 cell lines of genome-wide CRISPR fitness,
referenced only in 03_gdsc_eda.ipynb. It is independent of Chronos in screen,
library, and analysis pipeline.

If the near-zero correlation replicates, the central claim is on much firmer
ground than one dataset. If it does not, that must be known before it is relied
on -- which is the point of running this.

METHOD
------
1. Load scaled Bayesian factors. File layout:
     row 1 model_name | row 2 model_id (SIDM) | row 3 source | row 4 qc_pass
     row 5 header: gene_id, symbol, ensembl_gene_id, <cell line columns>
     row 6+ data, indexed by SIDG with ENSG in column 3
2. Columns are SIDM ids; bridge SIDM -> ACH via gdsc_scored_ready.parquet, which
   carries both. Lines screened by BOTH Broad and Sanger appear twice -- averaged.
3. Per gene, Spearman rho(core_score, essentiality) across shared cell lines,
   using the SAME thresholds as build_chronos_validation.py
   (rho>0.1 & p<0.05 = validated, rho<-0.1 & p<0.05 = inverted, MIN_N = 30).
4. Compare the resulting rho distribution, and the per-gene classifications,
   against chronos_validation.parquet.

SIGN CONVENTION
---------------
Scaled Bayesian factor: HIGHER = stronger evidence the gene is essential in that
line. chronos_validation.parquet already inverted Chronos so that HIGHER =
more essential. Both are therefore oriented the same way and rho values are
directly comparable. A sanity check on known core-essential genes (ribosomal
proteins) confirms the orientation rather than assuming it.

Output: src/pipeline/outputs/test_run_project_score_replication_results.json

Run:
    python src/pipeline/test_run_project_score_replication.py
"""
import json
import numpy as np
import pandas as pd
from pathlib import Path
from scipy.stats import spearmanr

OUTPUTS = Path("src/pipeline/outputs")
PS_FILE = ("data/GDSC/Project_score_combined_Sanger_v2_Broad_21Q2_"
           "fitness_scores_scaled_bayesian_factors_20250624.tsv")

RHO_THRESHOLD = 0.1      # matches build_chronos_validation.py
P_THRESHOLD   = 0.05
MIN_N         = 30

# ---------------------------------------------------------------- Step 1
print("=" * 72)
print("STEP 1 -- load Project Score fitness matrix")
print("=" * 72)

meta = pd.read_csv(PS_FILE, sep="\t", nrows=4, header=None, low_memory=False)
sidm_ids = meta.iloc[1, 3:].tolist()      # row 2 = model_id
sources  = meta.iloc[2, 3:].tolist()      # row 3 = source (Broad / Sanger)
print(f"  cell-line columns: {len(sidm_ids):,}  (distinct SIDM: {len(set(sidm_ids)):,})")

ps = pd.read_csv(PS_FILE, sep="\t", skiprows=5, header=None,
                 names=["gene_id", "symbol", "ensembl_gene_id"] + list(range(len(sidm_ids))),
                 low_memory=False)
print(f"  genes: {len(ps):,}")

ps = ps.dropna(subset=["ensembl_gene_id"])
ps["ensg_id"] = ps["ensembl_gene_id"].astype(str).str.split(".").str[0].str.lower()  # canonical case

vals = ps[list(range(len(sidm_ids)))].astype("float32")
vals.columns = sidm_ids
# average duplicate screens of the same line (Broad + Sanger)
vals = vals.T.groupby(level=0).mean().T
vals.index = ps["ensg_id"].values
print(f"  after averaging duplicate screens: {vals.shape[1]:,} distinct cell lines")

# ---------------------------------------------------------------- Step 2
print()
print("=" * 72)
print("STEP 2 -- bridge SIDM -> ACH model_id")
print("=" * 72)

gd = pd.read_parquet("validation/prepared/gdsc_scored_ready.parquet",
                     columns=["sidm_id", "model_id"]).drop_duplicates()
bridge = dict(zip(gd.sidm_id, gd.model_id.str.lower()))
mapped = [bridge.get(c) for c in vals.columns]
keep = [i for i, m in enumerate(mapped) if m is not None]
vals = vals.iloc[:, keep]
vals.columns = [mapped[i] for i in keep]
vals = vals.T.groupby(level=0).mean().T          # collapse any ACH duplicates
print(f"  bridged to ACH: {vals.shape[1]:,} cell lines")

# ---------------------------------------------------------------- Step 3
print()
print("=" * 72)
print("STEP 3 -- sanity check: are core-essential genes flagged essential?")
print("=" * 72)

gl = pd.read_parquet("reference/gene_lookup.parquet", columns=["ensg_id", "hgnc_symbol"])
gl["ensg_id"] = gl["ensg_id"].str.split(".").str[0].str.lower()  # canonical case
gl["hgnc_symbol"] = gl["hgnc_symbol"].astype("string").str.upper()  # symbol keys below are uppercase literals
sym2ensg = dict(zip(gl.hgnc_symbol, gl.ensg_id))

core_ess = ["RPL13A", "RPL5", "RPS6", "RPS3", "POLR2A", "EIF4A3"]
controls = ["GFP", "OR2T4", "KRTAP1-1", "DEFB126"]
gene_means = vals.mean(axis=1)
print("  known core-essential genes (expect HIGH if higher = more essential):")
ess_vals = []
for s in core_ess:
    e = sym2ensg.get(s)
    if e in gene_means.index:
        v = float(np.atleast_1d(gene_means.loc[e]).mean())
        ess_vals.append(v)
        print(f"    {s:9s} mean scaled BF = {v:8.3f}")
overall_mean = float(gene_means.mean())
print(f"  genome-wide mean scaled BF = {overall_mean:.3f}")
orientation_ok = bool(np.mean(ess_vals) > overall_mean) if ess_vals else None
print(f"  orientation confirmed (essential > genome mean): {orientation_ok}")

# ---------------------------------------------------------------- Step 4
print()
print("=" * 72)
print("STEP 4 -- per-gene rho(core_score, Project Score essentiality)")
print("=" * 72)

core = pd.read_parquet(OUTPUTS / "core_score.parquet",
                       columns=["model_id", "ensg_id", "core_score"])
core["model_id"] = core["model_id"].str.lower()

shared_lines = sorted(set(core.model_id) & set(vals.columns))
shared_genes = sorted(set(core.ensg_id) & set(vals.index))
print(f"  shared cell lines: {len(shared_lines):,}")
print(f"  shared genes     : {len(shared_genes):,}")

vals = vals.loc[~vals.index.duplicated(keep="first")]
ess = vals.loc[[g for g in shared_genes if g in vals.index], shared_lines]

core_w = (core[core.model_id.isin(shared_lines) & core.ensg_id.isin(ess.index)]
          .pivot_table(index="ensg_id", columns="model_id", values="core_score", aggfunc="mean")
          .reindex(index=ess.index, columns=shared_lines))

print("  computing per-gene Spearman (vectorised by ranking within gene)...")
A = core_w.to_numpy(dtype="float64")
B = ess.to_numpy(dtype="float64")
valid = ~np.isnan(A) & ~np.isnan(B)
n_per_gene = valid.sum(axis=1)

rows = []
for i in range(A.shape[0]):
    m = valid[i]
    n = int(m.sum())
    if n < MIN_N:
        continue
    r, p = spearmanr(A[i, m], B[i, m])
    if np.isnan(r):
        continue
    rows.append((ess.index[i], n, float(r), float(p)))

ps_res = pd.DataFrame(rows, columns=["ensg_id", "n", "rho", "pval"])
ps_res["check"] = np.where((ps_res.rho > RHO_THRESHOLD) & (ps_res.pval < P_THRESHOLD), "validated",
                    np.where((ps_res.rho < -RHO_THRESHOLD) & (ps_res.pval < P_THRESHOLD), "inverted", "none"))
print(f"  genes with n>={MIN_N}: {len(ps_res):,}")
print()
print(f"  PROJECT SCORE  mean rho = {ps_res.rho.mean():+.4f}   median = {ps_res.rho.median():+.4f}")
print(ps_res.check.value_counts().to_string())

# ---------------------------------------------------------------- Step 5
print()
print("=" * 72)
print("STEP 5 -- compare against Chronos")
print("=" * 72)

ch = pd.read_parquet(OUTPUTS / "chronos_validation.parquet")
print(f"  CHRONOS        mean rho = {ch.rho.mean():+.4f}   median = {ch.rho.median():+.4f}   n_genes={len(ch):,}")
print(f"  PROJECT SCORE  mean rho = {ps_res.rho.mean():+.4f}   median = {ps_res.rho.median():+.4f}   n_genes={len(ps_res):,}")

m = ch[["ensg_id", "rho", "chronos_check"]].merge(
    ps_res[["ensg_id", "rho", "check"]], on="ensg_id", suffixes=("_chronos", "_ps"))
print(f"  genes in both  : {len(m):,}")
agree_r, agree_p = spearmanr(m.rho_chronos, m.rho_ps)
print(f"  rho(chronos_rho, projectscore_rho) = {agree_r:+.4f}  p={agree_p:.2e}")

xtab = pd.crosstab(m.chronos_check, m.check)
print()
print("  classification agreement (rows=Chronos, cols=Project Score):")
print(xtab.to_string())
concord = float((m.chronos_check == m.check).mean())
print(f"  exact label agreement: {concord*100:.1f}%")

replicates = bool(abs(ps_res.rho.mean()) < 0.05)
print()
if replicates:
    verdict = (f"REPLICATES. Project Score mean rho = {ps_res.rho.mean():+.4f} on "
               f"{len(ps_res):,} genes, independently confirming that abundance does not "
               "track functional dependency. The headline finding no longer rests on a "
               "single essentiality dataset.")
else:
    verdict = (f"DOES NOT REPLICATE. Project Score mean rho = {ps_res.rho.mean():+.4f}, "
               "materially different from the near-zero Chronos result. The abundance-vs-"
               "dependency claim must be re-examined before it is relied on.")
print(f"  VERDICT: {verdict}")

results = {
    "question": "Does the near-zero abundance-vs-dependency correlation replicate on an independent CRISPR screen?",
    "read_only": True,
    "dataset": "Project Score combined Sanger v2 + Broad 21Q2, scaled Bayesian factors",
    "orientation_sanity_check": {
        "core_essential_mean": float(np.mean(ess_vals)) if ess_vals else None,
        "genome_wide_mean": overall_mean,
        "orientation_confirmed": orientation_ok,
    },
    "coverage": {
        "project_score_cell_lines_bridged": int(vals.shape[1]),
        "shared_cell_lines": len(shared_lines),
        "shared_genes": len(shared_genes),
        "genes_with_min_n": int(len(ps_res)),
        "min_n": MIN_N,
    },
    "project_score": {
        "mean_rho": round(float(ps_res.rho.mean()), 6),
        "median_rho": round(float(ps_res.rho.median()), 6),
        "classification": {k: int(v) for k, v in ps_res.check.value_counts().items()},
    },
    "chronos": {
        "mean_rho": round(float(ch.rho.mean()), 6),
        "median_rho": round(float(ch.rho.median()), 6),
        "n_genes": int(len(ch)),
    },
    "cross_dataset_agreement": {
        "genes_in_both": int(len(m)),
        "rho_of_rhos": round(float(agree_r), 6),
        "p": float(agree_p),
        "exact_label_agreement_pct": round(concord * 100, 2),
    },
    "replicates": replicates,
    "verdict": verdict,
}

out = OUTPUTS / "test_run_project_score_replication_results.json"
with open(out, "w") as f:
    json.dump(results, f, indent=2)
print()
print(f"Saved: {out}")
print("No production file was modified.")
