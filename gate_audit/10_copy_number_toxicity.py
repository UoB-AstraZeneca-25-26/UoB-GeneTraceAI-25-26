"""
gate_audit/10_copy_number_toxicity.py
--------------------------------------
AUDIT 10.  COPY-NUMBER-DRIVEN CUTTING TOXICITY.

ORIGIN
------
Audit 09 found: the gate appears under CRISPR knockout, does NOT appear under
RNAi knockdown (d1 ratio inverted, pAUC 0.477, p=4.8e-29), and the verdict was
CRISPR_only__gate_does_not_generalise_beyond_knockout.

Audit 09 identified the leading mechanistic hypothesis: the surviving CRISPR-
specific effect in the <5% prevalence / 20-29% depletion band may not be
abundance-dependency biology at all, but a copy-number-driven CRISPR cutting-
toxicity artifact (Aguirre 2016, Hart 2019, Meyers 2017). When a guide targets
a genomic region present in many copies (amplified locus), the excess DSBs cause
apparent dependency unrelated to the target gene's function. Critically, this
artifact scales with expression (expressed, amplified genes are the easiest CRISPR
targets) -- exactly the pattern the gate detects.

FALSIFIABLE PREDICTION
----------------------
If cutting toxicity drives the effect:
  For each gene, lines where that gene is amplified should show MORE negative
  CRISPR gene effects (appear more essential) than lines where that gene is not
  amplified -- even after controlling for actual dependency.

If the effect is genuine biology:
  Amplification status should NOT predict CRISPR gene effect after expression
  is accounted for (the gate is about expression, not copy number).

TEST DESIGN
-----------
CNA source: final_pipeline/outputs/cna_flags.parquet
  -- COSMIC CNA harmonised to ModelID / ENSG, produced by cna_layer.py.
  -- amplified = cna_call recorded as "amplification" in COSMIC (ploidy-relative
     by construction; COSMIC applies an internal threshold ~2.5x ploidy).
     Operationalised here as total_cn_mean >= CNA_AMP_THRESHOLD.

Per-gene analysis (requires MIN_PER_GROUP lines in each arm):
  dep_rate_amp   = fraction of amplified lines with CRISPR <= DEP_THRESHOLD
  dep_rate_noamp = fraction of non-amplified lines with CRISPR <= DEP_THRESHOLD
  amp_ratio      = dep_rate_amp / dep_rate_noamp      (>1 = cutting toxicity)
  effect_delta   = mean_crispr_amp - mean_crispr_noamp  (<0 = cutting toxicity)

Aggregate: Wilcoxon one-sided test of whether amp_ratio > 1 across genes.

Placebo shuffle: randomly permute amplification labels within each gene
(preserving per-gene marginal counts), repeat the analysis N_SHUFFLE times.
If the real aggregate median amp_ratio is not distinguishable from the shuffle
distribution, the observed signal is a labelling artifact, not a real CNA effect.

STEP 5 DECISION RULE
---------------------
If fewer than MIN_GENES_REQUIRED genes pass the MIN_PER_GROUP filter in both
arms, the test is underpowered. This script halts and reports the overlap size
rather than returning a spurious result.

Read-only. Writes only into gate_audit/outputs/.

Run:
    python gate_audit/10_copy_number_toxicity.py
    python gate_audit/10_copy_number_toxicity.py --n-genes 6000 --min-per-group 3
"""
import argparse
import json

import numpy as np
import pandas as pd
from scipy.stats import wilcoxon

import common as C

ap = argparse.ArgumentParser()
ap.add_argument("--n-genes", type=int, default=6000,
                help="gene sample cap (default 6000, 0 = all)")
ap.add_argument("--min-per-group", type=int, default=5,
                help="min lines per CNA arm per gene (default 5)")
ap.add_argument("--n-shuffle", type=int, default=200,
                help="placebo shuffle replicates (default 200)")
ap.add_argument("--amp-threshold", type=float, default=6.0,
                help="total_cn_mean >= this => amplified (default 6.0)")
args = ap.parse_args()

MIN_PER_GROUP = args.min_per_group
N_SHUFFLE = args.n_shuffle
AMP_THRESHOLD = args.amp_threshold
MIN_GENES_REQUIRED = 10     # halt threshold for Step 5

res = {}

C.banner("AUDIT 10 -- COPY-NUMBER-DRIVEN CRISPR CUTTING TOXICITY")
print(f"  min_per_group={MIN_PER_GROUP}  amp_threshold={AMP_THRESHOLD}  "
      f"n_shuffle={N_SHUFFLE}")

uni, valid_ensg, symbol_of, sym2ensg, u2e = C.load_gene_lookup()

# ------------------------------------------------------------------ load data
C.banner("10A -- DATA LOADING AND COVERAGE AUDIT (Step 5)")

W_cr, chk_cr = C.load_crispr_24q4(sym2ensg)
E = C.load_expression_24q4(sym2ensg)
res["scale_check"] = chk_cr

CNA_PATH = C.ROOT / "final_pipeline" / "outputs" / "cna_flags.parquet"
if not CNA_PATH.exists():
    raise SystemExit(f"{CNA_PATH} not found. "
                     "Run final_pipeline/03_Altercations/cna_layer.py first.")

cna = pd.read_parquet(CNA_PATH, columns=["model_id", "ensg_id", "total_cn_mean"])
cna["amplified"] = cna.total_cn_mean >= AMP_THRESHOLD

print(f"\n  CRISPR 24Q4: {W_cr.shape[0]} lines x {W_cr.shape[1]:,} genes")
print(f"  Expression 24Q4: {E.shape[0]} lines x {E.shape[1]:,} genes")
print(f"  CNA flags: {cna.model_id.nunique()} lines x "
      f"{cna.ensg_id.nunique():,} genes  (rows: {len(cna):,})")

cna_lines = set(cna.model_id)
cna_genes = set(cna.ensg_id)
crispr_lines = set(W_cr.index)
expr_lines = set(E.index)

line_overlap = crispr_lines & cna_lines
gene_overlap = set(W_cr.columns) & cna_genes & set(E.columns) & valid_ensg

print(f"\n  line overlap (CRISPR ∩ CNA): {len(line_overlap)}")
print(f"  gene overlap (CRISPR ∩ CNA ∩ expression ∩ valid): {len(gene_overlap):,}")

res["coverage"] = {
    "crispr_lines": int(W_cr.shape[0]),
    "cna_lines": int(cna.model_id.nunique()),
    "line_overlap": int(len(line_overlap)),
    "gene_overlap": int(len(gene_overlap)),
    "amp_threshold": float(AMP_THRESHOLD),
}

if len(line_overlap) < 50:
    msg = (f"INSUFFICIENT COVERAGE: only {len(line_overlap)} lines overlap "
           f"between CRISPR 24Q4 and CNA data. Minimum 50 required. HALTING.")
    print(f"\n  {msg}")
    res["verdict"] = "HALTED_insufficient_line_overlap"
    res["halt_reason"] = msg
    (C.OUT / "10_copy_number_toxicity_results.json").write_text(
        json.dumps(res, indent=2, default=float), encoding="utf-8")
    raise SystemExit(msg)

# ------------------------------------------------------------------ build CNA matrix
# Pivot cna_flags into model x gene matrix of total_cn_mean, restricted to
# the overlapping lines/genes so memory stays bounded.
C.banner("10B -- BUILD ANALYSIS SET")

genes = sorted(gene_overlap)
rng = np.random.default_rng(C.SEED)
if args.n_genes and args.n_genes < len(genes):
    genes = sorted(rng.choice(genes, args.n_genes, replace=False))
print(f"  analysis genes: {len(genes):,}")

# Filter cna to the overlap region, then pivot to model x gene
cna_sub = cna[cna.model_id.isin(line_overlap) & cna.ensg_id.isin(genes)]
cna_wide = cna_sub.pivot_table(
    index="model_id", columns="ensg_id",
    values="total_cn_mean", aggfunc="mean"
)
print(f"  CNA matrix (lines x genes): {cna_wide.shape[0]} x {cna_wide.shape[1]:,}")
print(f"  fill rate: {cna_wide.notna().mean().mean():.3f}  "
      f"(fraction of gene-line pairs with a CNA record)")

# Amplification fraction among lines with a CNA record
amp_frac = (cna_wide >= AMP_THRESHOLD).sum().sum() / cna_wide.notna().sum().sum()
print(f"  amplified fraction (total_cn >= {AMP_THRESHOLD}): {amp_frac:.3f}")
res["cna_matrix"] = {
    "n_lines": int(cna_wide.shape[0]),
    "n_genes": int(cna_wide.shape[1]),
    "fill_rate": float(cna_wide.notna().mean().mean()),
    "amp_frac": float(amp_frac),
}

# ------------------------------------------------------------------ per-gene analysis
C.banner("10C -- PER-GENE DEPLETION IN AMPLIFIED vs NON-AMPLIFIED LINES")

rows = []
for i, g in enumerate(genes):
    if i and i % 500 == 0:
        print(f"    {i:,}/{len(genes):,}")
    if g not in cna_wide.columns or g not in W_cr.columns or g not in E.columns:
        continue

    # Lines that have a CNA record for this gene AND CRISPR data
    cn_series = cna_wide[g].dropna()
    shared = sorted(set(cn_series.index) & set(W_cr.index))
    if len(shared) < 2 * MIN_PER_GROUP:
        continue

    cn_vals = cn_series.reindex(shared).values
    cr_vals = W_cr[g].reindex(shared).values
    amp_mask = cn_vals >= AMP_THRESHOLD
    noamp_mask = ~amp_mask

    n_amp = int(amp_mask.sum())
    n_noamp = int(noamp_mask.sum())
    if n_amp < MIN_PER_GROUP or n_noamp < MIN_PER_GROUP:
        continue

    dep_amp = float((cr_vals[amp_mask] <= C.DEP_THRESHOLD).mean())
    dep_noamp = float((cr_vals[noamp_mask] <= C.DEP_THRESHOLD).mean())
    eff_amp = float(np.nanmean(cr_vals[amp_mask]))
    eff_noamp = float(np.nanmean(cr_vals[noamp_mask]))

    amp_ratio = dep_amp / dep_noamp if dep_noamp > 0 else np.nan
    effect_delta = eff_amp - eff_noamp     # < 0 => amplified lines more essential

    rows.append({
        "ensg_id": g,
        "n_amp": n_amp, "n_noamp": n_noamp,
        "dep_amp": dep_amp, "dep_noamp": dep_noamp,
        "amp_ratio": amp_ratio,
        "eff_amp": eff_amp, "eff_noamp": eff_noamp,
        "effect_delta": effect_delta,
    })

R = pd.DataFrame(rows)
R_valid = R.dropna(subset=["amp_ratio"])
print(f"\n  genes with ≥{MIN_PER_GROUP} lines in each arm: {len(R):,}")
print(f"  genes with computable amp_ratio: {len(R_valid):,}")
res["per_gene_counts"] = {"eligible": int(len(R)), "with_ratio": int(len(R_valid))}

# ------------------------------------------------------------------ STEP 5 halt check
if len(R_valid) < MIN_GENES_REQUIRED:
    msg = (f"INSUFFICIENT COVERAGE: only {len(R_valid)} genes have ≥{MIN_PER_GROUP} "
           f"lines in each CNA arm. Minimum {MIN_GENES_REQUIRED} required "
           f"for any meaningful aggregate. HALTING. "
           f"Line overlap was {len(line_overlap)} lines.")
    print(f"\n  {msg}")
    res["verdict"] = "HALTED_insufficient_gene_coverage"
    res["halt_reason"] = msg
    if len(R):
        R.to_parquet(C.OUT / "10_copy_number_toxicity_per_gene.parquet", index=False)
    (C.OUT / "10_copy_number_toxicity_results.json").write_text(
        json.dumps(res, indent=2, default=float), encoding="utf-8")
    raise SystemExit(msg)

# ------------------------------------------------------------------ aggregate stats
C.banner("10D -- AGGREGATE: IS amp_ratio SYSTEMATICALLY > 1?")

ar = R_valid.amp_ratio
ed = R.effect_delta.dropna()

med_ratio = float(ar.median())
frac_above_1 = float((ar > 1.0).mean())
try:
    p_ratio = float(wilcoxon(ar - 1.0).pvalue)
except Exception:
    p_ratio = np.nan

med_delta = float(ed.median())
frac_negative = float((ed < 0).mean())
try:
    p_delta = float(wilcoxon(ed).pvalue)
except Exception:
    p_delta = np.nan

print(f"  {'metric':<30} {'value':>10}  {'interpretation'}")
print(f"  {'amp_ratio median':<30} {med_ratio:>10.4f}  (>1 = amplified lines more dependent)")
print(f"  {'frac genes amp_ratio > 1':<30} {frac_above_1:>10.3f}  (>0.5 = majority show artifact)")
print(f"  {'Wilcoxon p (ratio ≠ 1)':<30} {p_ratio:>10.4g}")
print(f"  {'effect_delta median':<30} {med_delta:>10.4f}  (<0 = amplified lines appear more essential)")
print(f"  {'frac genes delta < 0':<30} {frac_negative:>10.3f}")
print(f"  {'Wilcoxon p (delta ≠ 0)':<30} {p_delta:>10.4g}")

print(f"\n  dep_rate amplified median:   {R_valid.dep_amp.median():.4f}")
print(f"  dep_rate non-amplified med:  {R_valid.dep_noamp.median():.4f}")
print(f"  CRISPR eff amplified med:    {R.eff_amp.median():.4f}")
print(f"  CRISPR eff non-amplified:    {R.eff_noamp.median():.4f}")

res["aggregate"] = {
    "n_genes": int(len(R_valid)),
    "amp_ratio_median": med_ratio,
    "amp_ratio_frac_above_1": frac_above_1,
    "amp_ratio_wilcoxon_p": p_ratio,
    "effect_delta_median": med_delta,
    "effect_delta_frac_negative": frac_negative,
    "effect_delta_wilcoxon_p": p_delta,
    "dep_rate_amp_median": float(R_valid.dep_amp.median()),
    "dep_rate_noamp_median": float(R_valid.dep_noamp.median()),
    "crispr_eff_amp_median": float(R.eff_amp.median()),
    "crispr_eff_noamp_median": float(R.eff_noamp.median()),
}

# ------------------------------------------------------------------ placebo shuffle
C.banner("10E -- PLACEBO SHUFFLE (amplification labels randomised within each gene)")
print(f"  {N_SHUFFLE} replicates. If real result is indistinguishable from shuffle,")
print(f"  the signal is not a real CNA effect.")

# Build compact per-gene arrays for fast shuffling
gene_data = []
for g, d in R[["ensg_id","n_amp","n_noamp","dep_amp","dep_noamp","amp_ratio","effect_delta"]].iterrows():
    gene_data.append(d)

shuf_medians_ratio = []
shuf_medians_delta = []
rng_s = np.random.default_rng(C.SEED + 1)

# We shuffle at the combined-lines level per gene
# Build per-gene raw arrays
gene_arrays = {}
for i, g in enumerate(genes):
    if g not in cna_wide.columns or g not in W_cr.columns:
        continue
    cn_series = cna_wide[g].dropna()
    shared = sorted(set(cn_series.index) & set(W_cr.index))
    if len(shared) < 2 * MIN_PER_GROUP:
        continue
    cn_vals = cn_series.reindex(shared).values
    cr_vals = W_cr[g].reindex(shared).values
    amp_mask_orig = cn_vals >= AMP_THRESHOLD
    if amp_mask_orig.sum() < MIN_PER_GROUP or (~amp_mask_orig).sum() < MIN_PER_GROUP:
        continue
    gene_arrays[g] = (cr_vals, amp_mask_orig)

print(f"  genes available for shuffle: {len(gene_arrays):,}")

for rep in range(N_SHUFFLE):
    ratios, deltas = [], []
    for g, (cr_vals, amp_mask_orig) in gene_arrays.items():
        n = len(cr_vals)
        perm = rng_s.permutation(n)
        amp_shuf = amp_mask_orig[perm]
        n_amp_s = amp_shuf.sum()
        n_noamp_s = (~amp_shuf).sum()
        if n_amp_s < MIN_PER_GROUP or n_noamp_s < MIN_PER_GROUP:
            continue
        dep_amp_s = float((cr_vals[amp_shuf] <= C.DEP_THRESHOLD).mean())
        dep_noamp_s = float((cr_vals[~amp_shuf] <= C.DEP_THRESHOLD).mean())
        if dep_noamp_s > 0:
            ratios.append(dep_amp_s / dep_noamp_s)
        eff_amp_s = float(np.nanmean(cr_vals[amp_shuf]))
        eff_noamp_s = float(np.nanmean(cr_vals[~amp_shuf]))
        deltas.append(eff_amp_s - eff_noamp_s)
    shuf_medians_ratio.append(float(np.nanmedian(ratios)) if ratios else np.nan)
    shuf_medians_delta.append(float(np.nanmedian(deltas)) if deltas else np.nan)

shuf_ratio_arr = np.array([x for x in shuf_medians_ratio if not np.isnan(x)])
shuf_delta_arr = np.array([x for x in shuf_medians_delta if not np.isnan(x)])

shuf_ratio_med = float(np.median(shuf_ratio_arr)) if len(shuf_ratio_arr) else np.nan
shuf_delta_med = float(np.median(shuf_delta_arr)) if len(shuf_delta_arr) else np.nan

# One-sided percentile: fraction of shuffle replicates at or above the real value
p_shuf_ratio = float(np.mean(shuf_ratio_arr >= med_ratio)) if len(shuf_ratio_arr) else np.nan
p_shuf_delta = float(np.mean(shuf_delta_arr <= med_delta)) if len(shuf_delta_arr) else np.nan

print(f"\n  {'metric':<30} {'REAL':>10}  {'SHUFFLE med':>12}  {'p (shuffle ≥ real)':>18}")
print(f"  {'amp_ratio median':<30} {med_ratio:>10.4f}  {shuf_ratio_med:>12.4f}  {p_shuf_ratio:>18.4f}")
print(f"  {'effect_delta median':<30} {med_delta:>10.4f}  {shuf_delta_med:>12.4f}  {p_shuf_delta:>18.4f}")

res["placebo"] = {
    "n_shuffle": N_SHUFFLE,
    "n_genes_shuffled": int(len(gene_arrays)),
    "real_amp_ratio_median": med_ratio,
    "shuffle_amp_ratio_median": shuf_ratio_med,
    "shuffle_amp_ratio_p": p_shuf_ratio,
    "real_effect_delta_median": med_delta,
    "shuffle_effect_delta_median": shuf_delta_med,
    "shuffle_effect_delta_p": p_shuf_delta,
}

# ------------------------------------------------------------------ verdict
C.banner("VERDICT")

# Criteria for confirming the artifact:
# 1. amp_ratio > 1 (amplified lines show higher dep rate)
# 2. effect_delta < 0 (amplified lines appear more essential)
# 3. Both statistically significant (p < 0.05)
# 4. Effect is clearly larger than the placebo shuffle

confirms_ratio = (med_ratio > 1.05 and not np.isnan(p_ratio) and p_ratio < 0.05
                  and not np.isnan(p_shuf_ratio) and p_shuf_ratio < 0.1)
confirms_delta = (med_delta < -0.02 and not np.isnan(p_delta) and p_delta < 0.05
                  and not np.isnan(p_shuf_delta) and p_shuf_delta < 0.1)

print(f"  amp_ratio > 1 (amplified more dependent):  {med_ratio:.4f} "
      f"p={p_ratio:.4g}  shuffle_p={p_shuf_ratio:.4f}")
print(f"  effect_delta < 0 (amplified more essential): {med_delta:.4f} "
      f"p={p_delta:.4g}  shuffle_p={p_shuf_delta:.4f}")
print()

if confirms_ratio and confirms_delta:
    verdict = "CONFIRMS_cutting_toxicity_artifact__amplified_lines_systematically_more_CRISPR_dependent"
    print("  Both metrics consistent with cutting-toxicity artifact:")
    print("  amplified (gene, line) pairs show significantly higher CRISPR")
    print("  dependency rate AND more negative gene effect, above the placebo.")
    print("  The CRISPR-specific gate signal is at least partly explained by")
    print("  copy-number-driven excess DSBs in amplified loci.")
elif confirms_ratio or confirms_delta:
    verdict = "PARTIAL_evidence__one_metric_consistent_with_cutting_toxicity"
    print("  One of two metrics consistent with cutting-toxicity, but not both.")
    print("  Evidence is incomplete; alternative explanations remain viable.")
elif med_ratio < 1.0 and med_delta > 0:
    verdict = "DOES_NOT_CONFIRM__amplified_lines_show_LESS_apparent_dependency"
    print("  Amplified lines show LESS apparent CRISPR dependency than non-amplified.")
    print("  This is opposite to the cutting-toxicity prediction.")
    print("  Cutting-toxicity artifact does not explain the CRISPR-only gate.")
else:
    verdict = "INCONCLUSIVE__no_systematic_amplification_effect_detected"
    print("  No systematic relationship between amplification and CRISPR dependency.")
    print("  The test cannot confirm or rule out the cutting-toxicity hypothesis.")
    print("  Possible reasons: insufficient overlap, confounders, or true null.")

print(f"\n  VERDICT: {verdict}")
res["verdict"] = verdict

# ------------------------------------------------------------------ write outputs
R.to_parquet(C.OUT / "10_copy_number_toxicity_per_gene.parquet", index=False)
out_json = C.OUT / "10_copy_number_toxicity_results.json"
out_json.write_text(json.dumps(res, indent=2, default=float), encoding="utf-8")
print(f"\nwrote {out_json}")
print(f"wrote {C.OUT / '10_copy_number_toxicity_per_gene.parquet'}")
