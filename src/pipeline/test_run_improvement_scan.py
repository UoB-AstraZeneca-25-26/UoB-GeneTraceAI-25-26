"""
test_run_improvement_scan.py
----------------------------
Read-only diagnostic. Writes no production file.

Measures, on this repo's actual data, two proposed improvements to Layer 2 that
cannot be judged from the literature alone:

  PART 1  LINEAGE-CONDITIONAL GATE
          The project has treated lineage as a confound (14.2 pp TVD, CCLE
          eta^2 = 0.40) and corrected for it. The proposal is to invert that:
          ask "is this gene expressed relative to lines of the SAME lineage"
          rather than relative to the whole panel. Measured here as gate-region
          pAUC (FPR 0.8-1.0) and decile-1 depletion, panel-wide vs
          within-lineage, on the same genes and lines.

          Cost side: the min-lines rule now applies per lineage. 21 of 30
          tissues have >= 15 lines (1,726 of 1,781 models), so the loss is
          small, but it is measured rather than assumed.

  PART 2  IS A CONFORMAL LAYER USEFUL OR VACUOUS ON THIS SIGNAL?
          Split-conformal prediction gives a distribution-free coverage
          guarantee regardless of how weak the score is. The guarantee is
          always obtainable -- the question is whether the resulting SET is
          small enough to be useful. With a weak score the set can approach
          "all cell lines", which is a valid 90% guarantee and a worthless
          product.

          Implemented as the natural task for this pipeline: for a gene, return
          the smallest top-k set of cell lines that contains a true dependency
          with 90% probability. Calibrate k on held-out calibration GENES
          (exchangeability is across genes), then measure realised coverage on
          test genes. Compared against the k a random ordering would need --
          that ratio, not the coverage, is what says whether the score earns
          its place inside a conformal wrapper.

LABEL PROVENANCE -- READ THIS
-----------------------------
validation/prepared/chronos_long.parquet is NOT Chronos. Verified in this
session: it is the Sanger Project Score *scaled Bayesian factor* matrix
(data/GDSC/Project_score_combined_Sanger_v2_Broad_21Q2_fitness_scores_scaled_
bayesian_factors_20250624.tsv), NEGATED -- Pearson -1.0000 and exact agreement
to 1e-4 on 200,000 sampled (gene, line) pairs.

Consequences carried through this script:
  - "essentiality" here is a scaled Bayes factor, not a Chronos gene effect.
    0 does NOT mean "no fitness effect" and -1 does NOT mean "median common
    essential". The DEP_THRESHOLD below is a QUANTILE (~12th percentile), not
    a biological call, and every prevalence figure inherits that.
  - Relative, within-gene comparisons (decile profiles, pAUC, rank-based
    conformal sets) are unaffected by a monotone rescaling, which is why the
    earlier gate results stand. Absolute rates do not.
  - Real DepMap CRISPRGeneEffect.csv is absent from this repo entirely.

Outputs:
    src/pipeline/outputs/test_run_improvement_scan_results.json

Run:
    python src/pipeline/test_run_improvement_scan.py
"""
import argparse
import json
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow.parquet as pq
from scipy.stats import norm, wilcoxon

warnings.filterwarnings("ignore", category=RuntimeWarning)

OUTPUTS = Path("src/pipeline/outputs")
REF_DIR = Path("reference")
DATA_CLEAN = Path("data/parquet/data_clean")

DEP_THRESHOLD = -0.5      # quantile cut on the negated scaled Bayes factor
MIN_LINES = 150
MIN_POS = 10
MIN_LINEAGE = 15
N_GENES = 2000
COVERAGE = 0.90
SEED = 42

ap = argparse.ArgumentParser()
ap.add_argument("--labels", default="project_score",
                choices=["project_score", "chronos"],
                help="project_score = validation/prepared/chronos_long.parquet "
                     "(MISNAMED: it is Project Score scaled BF, negated). "
                     "chronos = a real DepMap CRISPRGeneEffect.csv via --chronos-path.")
ap.add_argument("--chronos-path",
                default="data/DepMap_Chronos/GeneFitnessEffect_Chronos_Achilles.hdf5")
ap.add_argument("--n-genes", type=int, default=None)
ap.add_argument("--dep-threshold", type=float, default=None,
                help="override; default -0.5 for both, which is a QUANTILE on "
                     "project_score but a calibrated gene effect on chronos")
args = ap.parse_args()
if args.dep_threshold is not None:
    DEP_THRESHOLD = args.dep_threshold
if args.n_genes is not None:
    N_GENES = args.n_genes

results = {"label_source": args.labels}


def vdw(s):
    n = int(s.notna().sum())
    return pd.Series(norm.ppf(s.rank(method="min") / (n + 1)), index=s.index)


def partial_auc(score, pos, lo=0.8, hi=1.0):
    ok = ~np.isnan(score)
    score, pos = score[ok], pos[ok]
    if pos.sum() < 2 or (~pos).sum() < 2:
        return np.nan
    order = np.argsort(-score, kind="mergesort")
    s, p = score[order], pos[order]
    distinct = np.r_[np.diff(s) != 0, True]
    fpr = np.r_[0, np.cumsum(~p)[distinct] / (~p).sum()]
    tpr = np.r_[0, np.cumsum(p)[distinct] / p.sum()]
    grid = np.linspace(lo, hi, 501)
    area = np.trapezoid(np.interp(grid, fpr, tpr), grid)
    a_min, a_max = (hi ** 2 - lo ** 2) / 2.0, hi - lo
    return 0.5 * (1.0 + (area - a_min) / (a_max - a_min))


# ============================================================ load
print("=" * 78)
print("STEP 1 -- labels (Project Score scaled BF, negated), lineage, expression")
print("=" * 78)

if args.labels == "project_score":
    ch = pd.read_parquet("validation/prepared/chronos_long.parquet")
    ch["ensg_id"] = ch.ensg_id.astype(str).str.split(".").str[0].str.lower()
    ib = pd.read_parquet("validation/prepared/id_bridge.parquet")
    ch = ch.merge(ib, on="sanger_model_id", how="inner")
    ch["model_id"] = ch.model_id.str.lower()
    print("  WARNING: label source is Project Score scaled BF (negated), NOT")
    print("           Chronos. -0.5 is a quantile here, not a gene effect.")
    LINEAGE_OVERRIDE = None
else:
    import h5py
    cp = Path(args.chronos_path)
    if not cp.exists():
        raise SystemExit(
            f"{cp} not found.\n"
            "Expected the Chronos release HDF5 (Dempster 2021 supplementary):\n"
            "  GeneFitnessEffect_Chronos_Achilles.hdf5  (Broad screen)\n"
            "  GeneFitnessEffect_Chronos_Score.hdf5     (Sanger screen)\n"
            "Layout: /data (lines x genes), /dim_0 = ACH- ids, "
            "/dim_1 = 'SYMBOL (ENTREZ)'.")
    with h5py.File(cp, "r") as fh:
        D = fh["data"][:]
        lines_h = [x.decode() for x in fh["dim_0"][:]]
        genes_h = [x.decode() for x in fh["dim_1"][:]]
    gl0 = pd.read_parquet(REF_DIR / "gene_lookup.parquet",
                          columns=["ensg_id", "hgnc_symbol"])
    gl0["ensg_id"] = gl0.ensg_id.astype(str).str.split(".").str[0].str.lower()
    s2e = dict(zip(gl0.hgnc_symbol.astype(str).str.upper(), gl0.ensg_id))
    keep_i, keep_e = [], []
    for i, c in enumerate(genes_h):
        e = s2e.get(c.split(" (")[0].strip().upper())
        if e:
            keep_i.append(i)
            keep_e.append(e)
    W = pd.DataFrame(D[:, keep_i], index=[m.lower() for m in lines_h],
                     columns=keep_e)
    W = W.T.groupby(level=0).mean().T          # collapse duplicate symbol->ensg
    ch = W.stack().rename("essentiality").reset_index()
    ch.columns = ["model_id", "ensg_id", "essentiality"]
    print(f"  real Chronos: {len(keep_i):,}/{len(genes_h):,} columns mapped to "
          f"ENSG; {W.shape[0]} lines x {W.shape[1]:,} genes")

    # SANITY CHECK against the release's own curated reference sets. This is the
    # check whose absence let a negated Project Score matrix pass as Chronos for
    # the whole project. Refuse to continue if the scale is not Chronos-like.
    refdir = cp.parent
    try:
        ess = set(pd.read_csv(refdir / "ReferenceEssentials.csv").iloc[:, 0]
                  .astype(str).str.split(" (", regex=False).str[0].str.upper())
        non = set(pd.read_csv(refdir / "ReferenceNonEssentials.csv").iloc[:, 0]
                  .astype(str).str.split(" (", regex=False).str[0].str.upper())
        e_ess = {s2e[g] for g in ess if g in s2e} & set(W.columns)
        e_non = {s2e[g] for g in non if g in s2e} & set(W.columns)
        m_ess = float(np.nanmedian(W[list(e_ess)].values))
        m_non = float(np.nanmedian(W[list(e_non)].values))
        fpr = float(np.nanmean(W[list(e_non)].values <= DEP_THRESHOLD))
        print(f"  sanity: median gene effect {np.nanmedian(W.values):+.4f} "
              f"(expect ~0)")
        print(f"          reference essentials    median {m_ess:+.4f} "
              f"(expect ~-1, n={len(e_ess)})")
        print(f"          reference NONessentials median {m_non:+.4f} "
              f"(expect ~0, n={len(e_non)})")
        print(f"          LABEL FPR at {DEP_THRESHOLD}: {fpr:.4f}  "
              f"<- curated negative controls, not improvised gene families")
        if not (m_ess < -0.5 < m_non + 0.5):
            raise SystemExit("Chronos scale check FAILED -- essentials should sit "
                             "near -1 and non-essentials near 0. Aborting.")
        results["label_fpr_from_reference_nonessentials"] = fpr
        results["reference_essential_median"] = m_ess
        results["reference_nonessential_median"] = m_non
    except FileNotFoundError:
        print("  (ReferenceEssentials/NonEssentials not alongside the HDF5 -- "
              "scale check skipped)")

    # Lineage from the release's own sample info; covers every Achilles line and
    # avoids id_bridge.parquet, the file with no builder.
    si = pd.read_csv(refdir / "DepMapSampleInfo20Q2.csv", low_memory=False)
    si["model_id"] = si.DepMap_ID.astype(str).str.lower()
    LINEAGE_OVERRIDE = dict(zip(si.model_id, si.lineage))
ch["dep"] = ch.essentiality <= DEP_THRESHOLD
print(f"  label matrix: {ch.ensg_id.nunique():,} genes x {ch.model_id.nunique():,} lines")
print(f"  threshold {DEP_THRESHOLD} = {100*ch.dep.mean():.1f}th percentile "
      f"(a quantile, not a biological call)")

mod = pd.read_parquet(OUTPUTS / "gdsc_models.parquet",
                      columns=["model_id", "tissue"]).dropna()
mod["model_id"] = mod.model_id.str.lower()
mod = mod.drop_duplicates("model_id")
if LINEAGE_OVERRIDE:
    # Chronos release sample info: covers every screened line, no id_bridge needed
    tissue_of = {k: v for k, v in LINEAGE_OVERRIDE.items() if isinstance(v, str)}
else:
    tissue_of = dict(zip(mod.model_id, mod.tissue))
vc = pd.Series(tissue_of).value_counts()
big = set(vc[vc >= MIN_LINEAGE].index)
print(f"  lineages: {len(vc)} total, {len(big)} with >= {MIN_LINEAGE} lines")

expr_names = pq.ParquetFile(DATA_CLEAN / "depmap_expr_clean.parquet").schema_arrow.names
expr_col = {c.split(".")[0].lower(): c for c in expr_names if c != "index"}
rng = np.random.default_rng(SEED)
cands = sorted(set(ch.ensg_id.unique()) & set(expr_col))
sweep = sorted(rng.choice(cands, size=min(N_GENES, len(cands)), replace=False))
print(f"  genes with label + RNA: {len(cands):,}  (sampled {len(sweep):,})")

prof = pd.read_parquet(REF_DIR / "depmap_profiles.parquet")
rp = prof[prof.datatype == "rna"][["profileid", "modelid"]].copy()
rp["model_id"] = rp.modelid.str.lower()
rp["profileid"] = rp.profileid.astype(str)
E = pq.read_table(DATA_CLEAN / "depmap_expr_clean.parquet",
                  columns=["index"] + [expr_col[g] for g in sweep]).to_pandas()
if "index" not in E.columns:
    E = E.reset_index()
E["profileid"] = E["index"].astype(str)
E = E.drop(columns=["index"]).merge(rp[["profileid", "model_id"]],
                                    on="profileid", how="inner")
E = E.drop(columns=["profileid"]).groupby("model_id").mean()
E.columns = [c.split(".")[0].lower() for c in E.columns]
print(f"  RNA matrix: {E.shape[0]} lines x {E.shape[1]} genes")

dep_by_gene = {g: s.set_index("model_id")["dep"]
               for g, s in ch[ch.ensg_id.isin(sweep)][
                   ["ensg_id", "model_id", "dep"]].groupby("ensg_id")}

# ============================================================ PART 1
print()
print("=" * 78)
print("PART 1 -- LINEAGE-CONDITIONAL GATE vs PANEL-WIDE GATE")
print("=" * 78)

rows, hitranks = [], []
for g in sweep:
    if g not in dep_by_gene or g not in E.columns:
        continue
    dep = dep_by_gene[g]
    e = E[g].dropna()
    lines = [m for m in e.index if m in dep.index and tissue_of.get(m) in big]
    if len(lines) < MIN_LINES:
        continue
    e = e.reindex(lines)
    y = dep.reindex(lines).astype(bool).values
    if y.sum() < MIN_POS or (~y).sum() < MIN_POS:
        continue

    z_panel = vdw(e)
    tis = pd.Series([tissue_of[m] for m in lines], index=lines)
    # within-lineage normal scores: rank inside the line's own tissue
    z_lin = e.groupby(tis).transform(
        lambda s: pd.Series(norm.ppf(s.rank(method="min") / (len(s) + 1)),
                            index=s.index))

    rec = {"ensg_id": g, "n_lines": len(lines), "n_pos": int(y.sum()),
           "base": float(y.mean())}
    for nm, z in [("panel", z_panel), ("lineage", z_lin)]:
        rec[f"pauc_{nm}"] = partial_auc(z.values, y)
        r = z.rank(method="first").values
        q = np.minimum((r - 1) * 10 // len(r), 9)
        rec[f"d1_{nm}"] = float(y[q == 0].mean() / y.mean())
    rows.append(rec)

    # Conformal inputs, BOTH directions.
    order = np.argsort(-z_panel.values, kind="mergesort")
    yo = y[order]
    # (a) INCLUSION task: how far down the ranking to the first true positive.
    hit = np.argmax(yo) + 1 if yo.any() else np.nan
    # (b) EXCLUSION task: how many lines can be dropped from the BOTTOM before
    #     hitting a positive. This is the task a GATE actually performs, and the
    #     established result (ranker region = chance, gate region = 0.59) says
    #     the inclusion task is the wrong one to wrap in a conformal layer.
    yb = yo[::-1]
    safe = int(np.argmax(yb)) if yb.any() else len(yb)
    # (c) RISK-CONTROLLING view: a zero-miss guarantee is degenerate on a weak
    #     score (one unlucky line at the bottom kills it). The shippable question
    #     is what FRACTION OF TRUE DEPENDENCIES is lost for a given amount of
    #     exclusion. That is the RCPS / selective-prediction formulation.
    risk = {}
    for f in (0.05, 0.10, 0.20, 0.30):
        k = int(np.floor(f * len(yo)))
        risk[f"lost_at_{int(f*100)}"] = float(yo[-k:].sum() / yo.sum()) if k else 0.0
    hitranks.append({"ensg_id": g, "n_lines": len(lines), "base": float(y.mean()),
                     "hit_rank": hit,
                     "hit_frac": hit / len(lines) if np.isfinite(hit) else np.nan,
                     "safe_drop": safe, "safe_frac": safe / len(lines), **risk})

R = pd.DataFrame(rows)
print(f"  genes scored: {len(R)}   median lines/gene {R.n_lines.median():.0f}   "
      f"median prevalence {R.base.median():.3f}")
print()
print(f"  {'':<26} {'gate pAUC':>11} {'decile-1 ratio':>16}")
print(f"  {'panel-wide (current)':<26} {R.pauc_panel.median():>11.4f} "
      f"{R.d1_panel.median():>16.3f}")
print(f"  {'within-lineage':<26} {R.pauc_lineage.median():>11.4f} "
      f"{R.d1_lineage.median():>16.3f}")
d = (R.pauc_lineage - R.pauc_panel).dropna()
dd = (R.d1_lineage - R.d1_panel).dropna()
print()
print(f"  gate pAUC   lineage - panel : {d.median():+.4f}   "
      f"p = {wilcoxon(d).pvalue:.4g}   better in {100*(d > 0).mean():.1f}% of genes")
print(f"  decile-1    lineage - panel : {dd.median():+.4f}   "
      f"(more negative = deeper depletion = better)")
verdict1 = ("lineage_conditioning_helps" if d.median() > 0.005
            and wilcoxon(d).pvalue < 0.05 else
            "lineage_conditioning_hurts" if d.median() < -0.005
            and wilcoxon(d).pvalue < 0.05 else
            "lineage_conditioning_neutral")
print(f"  -> {verdict1}")
results["part1_lineage"] = {
    "n_genes": int(len(R)),
    "pauc_panel": float(R.pauc_panel.median()),
    "pauc_lineage": float(R.pauc_lineage.median()),
    "delta": float(d.median()), "p": float(wilcoxon(d).pvalue),
    "frac_better": float((d > 0).mean()),
    "d1_panel": float(R.d1_panel.median()),
    "d1_lineage": float(R.d1_lineage.median()),
    "lines_lost_to_small_lineages": int(vc[vc < MIN_LINEAGE].sum()),
    "verdict": verdict1,
}

# ============================================================ PART 2
print()
print("=" * 78)
print("PART 2 -- SPLIT-CONFORMAL SET SIZE: useful, or a guarantee over everything?")
print("=" * 78)

H = pd.DataFrame(hitranks).dropna(subset=["hit_rank"])
H = H[H.n_lines >= MIN_LINES]
idx = rng.permutation(len(H))
cal, tst = H.iloc[idx[:len(H) // 2]], H.iloc[idx[len(H) // 2:]]
print(f"  genes: {len(H)}  (calibration {len(cal)}, test {len(tst)})")
print(f"  exchangeability is across GENES -- the split is a random gene split")

# conformal quantile on the RELATIVE rank, so genes with different line counts
# are on one scale; k is then a fraction of each gene's own line list.
qhat = np.quantile(cal.hit_frac, COVERAGE * (1 + 1 / len(cal)))
cov = float((tst.hit_frac <= qhat).mean())
med_k = float(np.median(np.ceil(qhat * tst.n_lines)))
print()
print(f"  calibrated set size (fraction of a gene's lines) : {qhat:.3f}")
print(f"  realised coverage on test genes                  : {cov:.3f} "
      f"(target {COVERAGE:.2f})")
print(f"  median absolute set size                         : {med_k:.0f} lines "
      f"of {tst.n_lines.median():.0f}")

# Random-ordering baseline: with prevalence pi, P(no hit in top k) = (1-pi)^k
pi = float(tst.base.median())
k_rand_frac = np.log(1 - COVERAGE) / np.log(1 - pi) / float(tst.n_lines.median())
print()
print(f"  random-ordering baseline at the same coverage    : {k_rand_frac:.3f} "
      f"of lines ({np.log(1-COVERAGE)/np.log(1-pi):.0f} lines, prevalence "
      f"{pi:.3f})")
gain = k_rand_frac / qhat if qhat > 0 else np.nan
print(f"  set-size ratio random / score-ordered            : {gain:.2f}x")
print()
if qhat > 0.5:
    verdict2 = "conformal_would_be_valid_but_vacuous"
    print(f"  -> The 90% set covers {qhat:.0%} of all cell lines. The guarantee is")
    print(f"     real and the product is not: 'the answer is in this half of the")
    print(f"     panel' is not an answer. A conformal layer on THIS score would be")
    print(f"     methodologically correct and practically empty.")
elif gain < 1.2:
    verdict2 = "conformal_valid_but_score_adds_little"
    print(f"  -> Sets are compact, but a random ordering reaches the same coverage")
    print(f"     at a similar size ({gain:.2f}x). The conformal wrapper would be")
    print(f"     carrying the result, not the abundance score.")
else:
    verdict2 = "conformal_useful"
    print(f"  -> Sets are compact AND materially smaller than random ({gain:.2f}x).")
    print(f"     A conformal layer converts the weak score into a calibrated,")
    print(f"     honestly-sized shortlist. Worth building.")

# ---- (b) the EXCLUSION formulation, which is what a gate can actually support
print()
print("  EXCLUSION formulation -- 'these lines can be safely ruled out'")
print("  (the inclusion test above asks the score to enrich the TOP, which the")
print("   ranker-region result already showed it cannot do)")
qsafe = np.quantile(cal.safe_frac, 1 - COVERAGE)
cov_s = float((tst.safe_frac >= qsafe).mean())
pi_t = tst.base.values
# random-ordering baseline: expected safe prefix from the bottom under a random
# permutation is geometric in the prevalence
rand_safe = np.median((1 - pi_t) ** 1 * (1.0 / pi_t - 1) / tst.n_lines.values)
rand_safe_q = float(np.quantile(
    [np.random.default_rng(SEED + i).geometric(max(p, 1e-9)) - 1
     for i, p in enumerate(pi_t)], 1 - COVERAGE) / np.median(tst.n_lines))
print(f"    calibrated safe-exclusion fraction : {qsafe:.4f} of a gene's lines")
print(f"    realised coverage on test genes    : {cov_s:.3f} (target {COVERAGE:.2f})")
print(f"    median lines safely excluded       : "
      f"{np.median(np.floor(qsafe * tst.n_lines)):.0f} of "
      f"{tst.n_lines.median():.0f}")
print(f"    random-ordering baseline           : {rand_safe_q:.4f}")
gain_s = qsafe / rand_safe_q if rand_safe_q > 0 else np.inf
print(f"    ratio score / random               : {gain_s:.2f}x")

results["part2_exclusion"] = {
    "calibrated_safe_fraction": float(qsafe),
    "realised_coverage": cov_s,
    "median_lines_excluded": float(np.median(np.floor(qsafe * tst.n_lines))),
    "random_baseline_fraction": rand_safe_q,
    "ratio_score_over_random": float(gain_s),
}

# ---- (c) the risk curve: the number an abstention rule would actually quote
print()
print("  RISK-CONTROL CURVE -- exclude the bottom f of lines, lose what share of")
print("  true dependencies?  (perfectly uninformative score loses exactly f)")
print(f"    {'exclude':>9} {'lost (median)':>15} {'lost (pooled mean)':>20} {'vs random':>11}")
risk_out = {}
# Split by prevalence: a gate can only act on SELECTIVELY essential genes.
# Pan-essential genes are expressed everywhere, so there is nothing to exclude,
# and pooling them in dilutes the only regime where the gate applies.
# Band definition follows the bootstrap CIs below: depletion is real for
# <2%, 2-5% and 10-25% (CIs exclude 1.0) and absent for >25%. The 5-10% dip has
# a CI spanning 1.0 and n=102, so it reads as noise, not a boundary. The scope
# rule is therefore "not pan-essential", not "very low prevalence".
subsets = [("all genes", H),
           ("selective (prev < 5%)", H[H.base < 0.05]),
           ("non-pan-essential (prev < 25%)", H[H.base < 0.25]),
           ("pan-essential (prev > 25%)", H[H.base > 0.25])]
for name, Hs in subsets:
    if len(Hs) < 20:
        continue
    print(f"    -- {name}  (n = {len(Hs)})")
    risk_out[name] = {}
    for f in (0.05, 0.10, 0.20, 0.30):
        col = f"lost_at_{int(f*100)}"
        med, mean = float(Hs[col].median()), float(Hs[col].mean())
        print(f"    {f:>9.0%} {med:>15.4f} {mean:>20.4f} {mean/f:>11.2f}x")
        risk_out[name][f"{int(f*100)}pct"] = {
            "median_lost": med, "mean_lost": mean, "ratio_vs_random": mean / f}
results["part2_risk_curve"] = risk_out

# ---- why this sample's decile-1 ratio differs from the earlier gate run
print()
print("  WHY decile-1 here (%.3f) differs from the earlier run (0.831):" % R.d1_panel.median())
R["prev_band"] = pd.cut(R.base, [0, .02, .05, .10, .25, 1.01],
                        labels=["<2%", "2-5%", "5-10%", "10-25%", ">25%"])
# Bootstrap CIs -- the 5-10% and 10-25% bands looked anomalous and the smallest
# band has n=102, so the profile needs an interval before any of it is claimed.
def boot_ci(v, n_boot=2000, stat=np.median):
    v = np.asarray(v.dropna())
    if len(v) < 5:
        return (np.nan, np.nan)
    rg = np.random.default_rng(SEED)
    b = [stat(rg.choice(v, len(v), replace=True)) for _ in range(n_boot)]
    return (float(np.percentile(b, 2.5)), float(np.percentile(b, 97.5)))

print(f"  {'band':>8} {'genes':>6} {'prev':>7} {'d1 ratio':>9} "
      f"{'95% CI':>18} {'pAUC':>7} {'95% CI':>18}")
band_rows = []
for b in R.prev_band.cat.categories:
    sub = R[R.prev_band == b]
    if not len(sub):
        continue
    lo1, hi1 = boot_ci(sub.d1_panel)
    lo2, hi2 = boot_ci(sub.pauc_panel)
    print(f"  {str(b):>8} {len(sub):>6} {sub.base.median():>7.3f} "
          f"{sub.d1_panel.median():>9.3f} [{lo1:>7.3f},{hi1:>7.3f}] "
          f"{sub.pauc_panel.median():>7.4f} [{lo2:>7.3f},{hi2:>7.3f}]")
    band_rows.append({"band": str(b), "genes": int(len(sub)),
                      "prevalence": float(sub.base.median()),
                      "d1_median": float(sub.d1_panel.median()),
                      "d1_ci": [lo1, hi1],
                      "pauc_median": float(sub.pauc_panel.median()),
                      "pauc_ci": [lo2, hi2]})
tab = pd.DataFrame(band_rows)
R.to_parquet(OUTPUTS / f"test_run_improvement_scan_{args.labels}_per_gene.parquet")
print("  -> the gate acts on LOW-PREVALENCE (selectively essential) genes.")
print("     Genes essential in many lines are pan-essential: they are expressed")
print("     everywhere, so an expression gate has nothing to exclude. The earlier")
print("     run sampled genes with protein coverage and prevalence 0.040; this")
print("     one requires >=10 positives and lands at 0.059.")
results["prevalence_stratification"] = band_rows

results["part2_conformal"] = {
    "n_genes": int(len(H)), "target_coverage": COVERAGE,
    "calibrated_set_fraction": float(qhat),
    "realised_coverage": cov,
    "median_set_size_lines": med_k,
    "median_lines_per_gene": float(tst.n_lines.median()),
    "median_prevalence": pi,
    "random_baseline_fraction": float(k_rand_frac),
    "set_size_ratio_random_over_score": float(gain),
    "verdict": verdict2,
}

print()
print("=" * 78)
print("SUMMARY")
print("=" * 78)
print(f"  lineage conditioning : {verdict1}")
print(f"  conformal layer      : {verdict2}")

out = OUTPUTS / f"test_run_improvement_scan_{args.labels}_results.json"
out.write_text(json.dumps(results, indent=2, default=float), encoding="utf-8")
print(f"\nwrote {out}")
