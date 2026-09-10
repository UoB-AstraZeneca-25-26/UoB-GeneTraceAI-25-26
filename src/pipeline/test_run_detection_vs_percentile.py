"""
test_run_detection_vs_percentile.py
-----------------------------------
Read-only diagnostic. Writes no production file.

QUESTION
--------
test_run_abundance_gate.py established that core_score works as a GATE and not as
a ranker: gate-region pAUC 0.596, ranker-region pAUC 0.500, decile-1 dependency
ratio 0.84, decile-10 ratio 0.99. All of the measured signal sits at the bottom
of the score.

If that is the whole story, the continuous within-gene percentile is doing no
work beyond a binary call -- "is this gene on or off in this line" -- and could
be replaced by an ABSOLUTE detection threshold with no loss. That would be a
large architectural simplification:

  * an absolute threshold is panel-independent. The percentile's denominator is
    whichever cell lines happen to be in the build, so every score shifts when
    lines are added. A detection call does not.
  * GEO and HPA could be added to the expression layer without renormalising
    anything, because each source would be thresholded on its own scale.
  * `n_layers`, `stratum_rank` and the percentile machinery in 02_core_score
    would collapse to a boolean.

So: does the percentile carry graded information, or is it an expensive encoding
of a boolean?

A NOTE THAT DECIDES HOW THIS TEST HAD TO BE BUILT
-------------------------------------------------
Within one gene, the percentile is a strictly monotone transform of the raw
log2(TPM+1) value. Every rank-based metric -- AUROC, partial AUROC, decile
profile -- is invariant under monotone transforms, so

    pAUC(percentile)  ==  pAUC(raw log2TPM)      EXACTLY, by construction.

The percentile therefore cannot win or lose against raw expression on any of
those metrics. Its only claim over the raw value is cross-gene comparability.
The thing that genuinely differs is a PER-SAMPLE normalised score, because that
reorders cell lines: each line gets its own reference distribution. That is what
an absolute detection call is, and it is what this script compares against.

WHY THE FIRST VERSION'S THRESHOLD WAS WRONG  (v1, superseded)
-------------------------------------------------------------
v1 implemented Hart et al. 2013 zFPKM: fit a kernel density to per-line
log-expression across genes, take the mode `mu` as the centre of the unexpressed
population, estimate sigma from its right-hand half-normal, call expressed at
z >= -3.

On this data the mode landed at mu = +4.13 log2TPM (TPM 17.5) -- the EXPRESSED
peak, not the unexpressed one. depmap_expr_clean.parquet is log2(TPM+1), so the
+1 pseudocount collapses the unexpressed population onto a spike at exactly 0
(17.4% of protein-coding cells). Inverting the pseudocount recovers the spread,
but the low-TPM tail is then smeared across a wide negative range in log space
and its density never exceeds the expressed peak, so argmax picks the wrong one.
The rule v1 actually tested was TPM >= 0.232, arrived at by accident.

THE FIX -- A CURATED UNEXPRESSED REFERENCE, NOT A MODE SEARCH
-------------------------------------------------------------
Rather than infer the unexpressed population from the shape of the mixture, seed
it from genes known not to be expressed in cultured cell lines: the olfactory
(OR*), taste (TAS2R*) and vomeronasal (VN1R*) receptor families. These are the
same families test_run_gate_mechanism.py used as its negative control, so the
calibration set is the one this repo has already characterised rather than a new
one. Their expression is verified here per run, not assumed.

Per cell line, those genes give an empirical unexpressed distribution. A gene is
called expressed in that line when its log2TPM exceeds a quantile of it:

    expressed(g, line)  <=>  log2TPM(g, line)  >  Q_q( reference | line )

`q` is a per-line FALSE POSITIVE RATE on a real negative control set, which is
the same logic the harness family already uses for LABEL_FPR. That makes the
threshold interpretable and the sweep meaningful: q = 0.90 admits 10% of known-
unexpressed genes as "expressed", q = 0.999 admits 0.1%.

Reported alongside:
  * a MODE-CONSTRAINED Hart variant, argmax restricted left of the per-line
    median, to show what the corrected version of v1's estimator gives
  * TPM >= 1, the dead-simple rule, since "much simpler" is part of the claim
  * a continuous per-line z, z = (log2TPM - mu_ref) / sigma_ref, for pAUC --
    a binary call has a two-point ROC whose pAUC over FPR in [0.8, 1.0] is
    mostly interpolation, so binary rules are NOT evaluated by pAUC here

FIVE READOUTS
-------------
  1. GATE-REGION pAUC, percentile vs continuous per-line z.

  2. MATCHED-EXCLUSION, per rule. The rule excludes n lines for a gene; compare
     the dependency rate in those n against the bottom n by percentile. Same
     cost, directly comparable, no ROC geometry.

  3. PER-GENE 2x2 ACROSS THE WHOLE SWEEP. Cross-tabulate {rule says not
     expressed} x {percentile decile 1} and report per-gene median dependency
     ratios in each cell, at every threshold. v1 anchored this to one accidental
     cut; the conclusion must not depend on the cut.

     The cell that answers the question is "in decile 1 but the rule says
     EXPRESSED": if its dependency rate sits at the gene's base rate, the
     percentile adds nothing detection does not already have.

  4. STRUCTURAL OVERLAP. v1 found only 21 genes where BOTH incremental cells had
     >= 5 lines, and reported an underpowered paired test without explaining
     why. The reason is geometric, not statistical: decile 1 is exactly 10% of a
     gene's lines by construction, so when a rule excludes far fewer than 10% its
     exclusions are swallowed by decile 1 (abs_only empty), and when it excludes
     far more, decile 1 is swallowed by the exclusions (pct_only empty). Both
     cells are populated only where the rule's exclusion fraction sits near 10%
     without equalling it. This section measures that containment directly and
     reports how many genes can support a paired comparison at each threshold,
     so the n is explained rather than discovered.

  5. FLOOR DEPLETION, judged against the Chronos LABEL FPR from the release's own
     curated non-essentials. A floor rate below the label's own false positive
     rate cannot be read as residual dependency.

EFFECT-SIZE FLOORS
------------------
v1's verdict logic tested significance only. Pooled over 1.15M (gene, line)
pairs, an odds ratio of 0.87 returned p = 3e-60 -- entirely a large-n artefact.
Every claim below must now clear BOTH a p-value and an effect-size floor:

    MIN_RATIO_EFFECT  0.10   |1 - dependency ratio| must reach this
    MIN_PAUC_DELTA    0.01   paired pAUC difference must reach this
    MAX_OR_EFFECT     0.80   Fisher odds ratio must be at or below this

Output:
    src/pipeline/outputs/test_run_detection_vs_percentile_<label>_results.json
    src/pipeline/outputs/test_run_detection_vs_percentile_<label>_per_gene.parquet
        (long format: one row per gene x threshold rule)

Run:
    python src/pipeline/test_run_detection_vs_percentile.py --label chronos
    python src/pipeline/test_run_detection_vs_percentile.py --label gdsc
"""
import argparse
import json
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow.parquet as pq
from scipy.ndimage import gaussian_filter1d
from scipy.stats import hypergeom, rankdata, wilcoxon, fisher_exact

warnings.filterwarnings("ignore", category=RuntimeWarning)

OUTPUTS = Path("src/pipeline/outputs")
REF_DIR = Path("reference")
DATA_CLEAN = Path("data/parquet/data_clean")

MIN_POS = 10
MIN_NEG = 10
MIN_LINES = 100
MIN_CELL = 5                 # lines needed in a 2x2 cell to quote a rate
DEP_THRESHOLD = -0.5         # chronos: essentiality at or below this = dependent

# --- effect-size floors (see docstring) -----------------------------------
MIN_RATIO_EFFECT = 0.10
MIN_PAUC_DELTA = 0.01
MAX_OR_EFFECT = 0.80

# --- the sweep ------------------------------------------------------------
REF_QUANTILES = (0.90, 0.95, 0.99, 0.999)
HART_CUT = -3.0
TPM_CUT = 1.0
KDE_BINS = 1024
SEED = 42

# --- half-sample held-out centre selection --------------------------------
# The gate-region pAUC of the per-line score depends on which quantile centre it
# is built around: 0.5940 at Q0.90 up to 0.6486 at Q0.99, with the paired delta
# against the percentile changing SIGN across that range. Picking the best centre
# and then quoting its pAUC on the same labels is selection on the test set. This
# splits the cell lines, selects the centre on one half by gate pAUC, and
# evaluates the fixed choice on the other half. The percentile is RECOMPUTED
# within the evaluation half, because it is panel-dependent by construction and
# the honest comparison gives both scores the same lines.
HALF_SAMPLE_REPS = 20
MIN_HALF_POS = 5
MIN_HALF_NEG = 5

# --- paralog stratification -----------------------------------------------
# Paralog buffering attenuates single-knockout phenotypes: a gene with a close
# functional paralog can be knocked out with little fitness cost because the
# paralog covers it. That depresses the Chronos label for exactly those genes,
# independently of expression, so it can bias the incremental cell this test
# turns on.
#
# PROXY, and stated as one: HGNC gene-group co-membership, not sequence identity.
# The repo carries no paralog annotation (gene_lookup.parquet has no such column;
# the only prior handling is test_run_gate_mechanism.py excluding KRTAP and VN1R
# as paralogy-contaminated control families). HGNC gene_group is the closest
# thing on disk. Its large groups are DOMAIN superfamilies rather than paralog
# families -- "Zinc fingers C2H2-type" has 763 members, "CD molecules" 385 --
# and co-membership there implies no paralogy at all, so groups above
# PARALOG_MAX_GROUP are dropped. The annotation is then validated against the
# prediction that paralogous genes show LOWER dependency prevalence; if that
# fails, it is not capturing paralogy and the stratification is not interpreted.
HGNC_FILE = Path("data/gene_with_protein_product.txt")
PARALOG_MAX_GROUP = 30
PARALOG_MAX_GROUP_STRICT = 10

# Curated unexpressed reference families. Same families test_run_gate_mechanism.py
# used as its negative control set.
REF_PREFIXES = ("OR", "TAS2R", "VN1R")
# OR* must not catch ORAI1, ORC1, ORM1 ... olfactory receptor symbols are
# OR<digit><Letter><digits>, e.g. OR4F5, OR10A2. Anchor on a digit.
import re
OR_RE = re.compile(r"^OR\d+[A-Z]+\d*$")
TAS_RE = re.compile(r"^TAS2R\d+[A-Z]*$")
VN_RE = re.compile(r"^VN1R\d+[A-Z]*$")


def is_ref_family(sym):
    s = str(sym).upper()
    return bool(OR_RE.match(s) or TAS_RE.match(s) or VN_RE.match(s))


ap = argparse.ArgumentParser()
ap.add_argument("--label", choices=["gdsc", "chronos"], default="chronos")
ap.add_argument("--chronos-path",
                default="data/DepMap_Chronos/GeneFitnessEffect_Chronos_Achilles.hdf5")
ap.add_argument("--n-genes", type=int, default=3000)
args = ap.parse_args()
LABEL = args.label

results = {
    "label_source": LABEL,
    "method_note": ("percentile is a monotone transform of raw log2TPM within a "
                    "gene, so pAUC(percentile) == pAUC(raw) exactly; the contrast "
                    "is against a per-line absolute detection call"),
    "effect_size_floors": {"min_ratio_effect": MIN_RATIO_EFFECT,
                           "min_pauc_delta": MIN_PAUC_DELTA,
                           "max_or_effect": MAX_OR_EFFECT},
}


# ============================================================ helpers
def roc_curve(score, pos):
    order = np.argsort(-score, kind="mergesort")
    s, p = score[order], pos[order]
    P, N = p.sum(), (~p).sum()
    distinct = np.r_[np.diff(s) != 0, True]
    tp = np.cumsum(p)[distinct]
    fp = np.cumsum(~p)[distinct]
    return np.r_[0, fp / N], np.r_[0, tp / P]


def partial_auc(score, pos, lo, hi):
    """McClish-standardised partial AUROC over FPR in [lo, hi]. 0.5 = chance."""
    ok = np.isfinite(score)
    score, pos = score[ok], pos[ok]
    if pos.sum() < 2 or (~pos).sum() < 2:
        return np.nan
    fpr, tpr = roc_curve(score, pos)
    grid = np.linspace(lo, hi, 501)
    t = np.interp(grid, fpr, tpr)
    area = np.trapezoid(t, grid)
    a_min = (hi ** 2 - lo ** 2) / 2.0
    a_max = hi - lo
    return 0.5 * (1.0 + (area - a_min) / (a_max - a_min))


def auroc(score, pos):
    ok = np.isfinite(score)
    score, pos = score[ok], pos[ok]
    n1, n0 = int(pos.sum()), int((~pos).sum())
    if n1 < 1 or n0 < 1:
        return np.nan
    r = rankdata(score)
    return (r[pos].sum() - n1 * (n1 + 1) / 2.0) / (n1 * n0)


def boot_ci(v, n_boot=2000, stat=np.median):
    v = np.asarray(pd.Series(v).dropna())
    if len(v) < 5:
        return (np.nan, np.nan)
    rg = np.random.default_rng(SEED)
    b = [stat(rg.choice(v, len(v), replace=True)) for _ in range(n_boot)]
    return (float(np.percentile(b, 2.5)), float(np.percentile(b, 97.5)))


def wilx(v, mu=0.0):
    v = pd.Series(v).dropna()
    if len(v) < 6 or np.allclose(v.values, mu):
        return np.nan
    try:
        return float(wilcoxon(v - mu).pvalue)
    except ValueError:
        return np.nan


# ============================================================ Step 1
print("=" * 78)
print("STEP 1 -- gene universe and expression matrix")
print("=" * 78)

gl = pd.read_parquet(REF_DIR / "gene_lookup.parquet",
                     columns=["ensg_id", "hgnc_symbol", "biotype", "hgnc_status"])
gl["ensg_id"] = gl.ensg_id.astype("string").str.split(".").str[0].str.lower()
pc = gl[(gl.biotype == "protein_coding") & (gl.hgnc_status == "Approved")]
valid_ensg = set(pc.ensg_id.dropna())
symbol_of = dict(zip(pc.ensg_id, pc.hgnc_symbol))
s2e = dict(zip(pc.hgnc_symbol.astype(str).str.upper(), pc.ensg_id))
print(f"  protein-coding, HGNC-approved genes : {len(valid_ensg):,}")

expr_names = pq.ParquetFile(DATA_CLEAN / "depmap_expr_clean.parquet").schema_arrow.names
expr_col = {}
for c in expr_names:
    if c == "index":
        continue
    e = c.split(".")[0].lower()
    if e in valid_ensg:
        expr_col.setdefault(e, c)
print(f"  ... with an expression column       : {len(expr_col):,}")

prof = pd.read_parquet(REF_DIR / "depmap_profiles.parquet")
rna_prof = prof[prof.datatype == "rna"][["profileid", "modelid"]].copy()
rna_prof["model_id"] = rna_prof.modelid.str.lower()
rna_prof["profileid"] = rna_prof.profileid.astype(str)

ordered_ensg = sorted(expr_col)
E = pq.read_table(DATA_CLEAN / "depmap_expr_clean.parquet",
                  columns=["index"] + [expr_col[g] for g in ordered_ensg]).to_pandas()
# pyarrow restores "index" as the pandas index, not a column.
if "index" not in E.columns:
    E = E.reset_index()
E["profileid"] = E["index"].astype(str)
E = E.drop(columns=["index"]).merge(rna_prof[["profileid", "model_id"]],
                                   on="profileid", how="inner")
E = E.drop(columns=["profileid"]).groupby("model_id").mean()
E.columns = ordered_ensg
E = E.astype("float32")
n_lines, n_genes_all = E.shape
print(f"  expression matrix                   : {n_lines} lines x "
      f"{n_genes_all:,} genes")

# ------------------------------------------------------- paralog annotation
print()
print("  PARALOG ANNOTATION (HGNC gene-group co-membership -- a PROXY)")
if not HGNC_FILE.exists():
    raise SystemExit(f"{HGNC_FILE} not found; needed for paralog stratification.")
hg = pd.read_csv(HGNC_FILE, sep="\t", low_memory=False,
                 usecols=["symbol", "ensembl_gene_id", "gene_group"])
hg = hg.dropna(subset=["ensembl_gene_id"])
hg["ensg_id"] = hg.ensembl_gene_id.astype(str).str.split(".").str[0].str.lower()
hg = hg[hg.ensg_id.isin(set(ordered_ensg))]
mem = (hg.assign(grp=hg.gene_group.fillna("").str.split("|"))
         .explode("grp"))
mem["grp"] = mem.grp.str.strip()
mem = mem[mem.grp != ""][["ensg_id", "grp"]].drop_duplicates()
gsize = mem.groupby("grp").ensg_id.nunique()


def paralog_counts(cap):
    """Number of OTHER expression-matrix genes sharing a small HGNC group."""
    keep = set(gsize[gsize <= cap].index)
    m = mem[mem.grp.isin(keep)]
    partners = {}
    for grp, sub in m.groupby("grp").ensg_id:
        members = list(sub)
        for g in members:
            partners.setdefault(g, set()).update(members)
    return {g: len(v - {g}) for g, v in partners.items()}, len(keep)


par_n, n_grp = paralog_counts(PARALOG_MAX_GROUP)
par_n_strict, n_grp_s = paralog_counts(PARALOG_MAX_GROUP_STRICT)
print(f"    groups retained at size <= {PARALOG_MAX_GROUP:<3} : {n_grp} "
      f"of {len(gsize)} groups total")
print(f"    groups retained at size <= {PARALOG_MAX_GROUP_STRICT:<3} : {n_grp_s}")
print(f"    genes with >=1 paralog (cap {PARALOG_MAX_GROUP}) : "
      f"{sum(1 for v in par_n.values() if v):,} of {len(ordered_ensg):,}")
print(f"    genes with >=1 paralog (cap {PARALOG_MAX_GROUP_STRICT}) : "
      f"{sum(1 for v in par_n_strict.values() if v):,}")
results["paralog_annotation"] = {
    "source": "HGNC gene_group co-membership (PROXY, not sequence identity)",
    "max_group_size": PARALOG_MAX_GROUP,
    "max_group_size_strict": PARALOG_MAX_GROUP_STRICT,
    "groups_retained": int(n_grp), "groups_retained_strict": int(n_grp_s),
    "genes_with_paralog": int(sum(1 for v in par_n.values() if v)),
    "genes_with_paralog_strict": int(sum(1 for v in par_n_strict.values() if v)),
}

# --- invert the pseudocount ------------------------------------------------
X = E.values
TPM = np.expm1(X * np.log(2.0))
TPM[TPM < 0] = 0.0
zero_mask = TPM <= 0
frac_zero = float(zero_mask.mean())
with np.errstate(divide="ignore"):
    L = np.log2(TPM, where=~zero_mask, out=np.full_like(TPM, -np.inf))
print(f"  exact-zero cells (log2TPM = -inf)   : {frac_zero:.4f}")


# ============================================================ Step 2
print()
print("=" * 78)
print("STEP 2 -- CURATED UNEXPRESSED REFERENCE  (the v1 fix)")
print("=" * 78)

ref_ensg = [g for g in ordered_ensg if is_ref_family(symbol_of.get(g, ""))]
ref_idx = np.array([ordered_ensg.index(g) for g in ref_ensg])
print(f"  reference families {REF_PREFIXES} : {len(ref_ensg)} genes with an "
      f"expression column")
ex_syms = [symbol_of.get(g) for g in ref_ensg[:6]]
print(f"  e.g. {ex_syms}")

# Verify -- do not assume -- that these are actually unexpressed here.
REF = L[:, ref_idx]
ref_finite = REF[np.isfinite(REF)]
all_finite = L[np.isfinite(L)]
print(f"  reference genes : median log2TPM {np.median(ref_finite):+.3f}, "
      f"{100*np.mean(~np.isfinite(REF)):.1f}% exact zero")
print(f"  all genes       : median log2TPM {np.median(all_finite):+.3f}, "
      f"{100*frac_zero:.1f}% exact zero")
ref_med = float(np.median(ref_finite))
all_med = float(np.median(all_finite))
if not (ref_med < all_med - 1.0):
    raise SystemExit("reference families are NOT clearly lower than the "
                     "transcriptome. Calibration set invalid. Aborting.")
print(f"  -> reference sits {all_med - ref_med:.2f} log2 units below the "
      f"transcriptome median. Usable.")

# Per-line reference quantiles. Exact zeros count as the lowest possible value,
# so they are included via a finite sentinel below the observed minimum.
FLOOR_SENTINEL = float(np.min(all_finite)) - 1.0
REF_f = np.where(np.isfinite(REF), REF, FLOOR_SENTINEL)
Lf = np.where(np.isfinite(L), L, FLOOR_SENTINEL)

ref_mu = np.median(REF_f, axis=1)
ref_sd = 1.4826 * np.median(np.abs(REF_f - ref_mu[:, None]), axis=1)
ref_sd[ref_sd <= 0] = np.nan
print(f"  per-line reference mu    : median {np.nanmedian(ref_mu):+.3f} "
      f"[{np.nanpercentile(ref_mu, 5):+.3f}, {np.nanpercentile(ref_mu, 95):+.3f}]")
print(f"  per-line reference sigma : median {np.nanmedian(ref_sd):.3f}")

ref_q = {q: np.quantile(REF_f, q, axis=1) for q in REF_QUANTILES}
for q in REF_QUANTILES:
    print(f"    Q{q:<6} of reference : median cut log2TPM "
          f"{np.median(ref_q[q]):+.3f}  (TPM {2**np.median(ref_q[q]):.3f})")

# --- mode-constrained Hart variant, for comparison with the corrected v1 ----
lo_edge = float(np.percentile(all_finite, 0.05))
hi_edge = float(np.percentile(all_finite, 99.95))
grid = np.linspace(lo_edge, hi_edge, KDE_BINS)
binw = grid[1] - grid[0]
h_mu = np.full(n_lines, np.nan)
h_sd = np.full(n_lines, np.nan)
for i in range(n_lines):
    v = L[i][np.isfinite(L[i])]
    if len(v) < 1000:
        continue
    med_i = np.median(v)
    hist, _ = np.histogram(v, bins=KDE_BINS, range=(lo_edge, hi_edge))
    sd = v.std()
    iqr = np.subtract(*np.percentile(v, [75, 25]))
    bw = 0.9 * min(sd, iqr / 1.34) * len(v) ** (-0.2)
    sm = gaussian_filter1d(hist.astype(float), max(1.0, bw / binw))
    # CONSTRAIN the mode search left of the per-line median -- this is the v1 fix
    left = grid < med_i
    if not left.any():
        continue
    mu = grid[np.argmax(np.where(left, sm, -np.inf))]
    # sigma from the right-hand half-normal, as Hart specifies. An earlier
    # version also bounded this band above by the per-line median, which left
    # fewer than 100 values for almost every line, so h_mu stayed NaN and
    # `Lf >= NaN` evaluated to False -- the rule then read as "nothing is
    # expressed" rather than "this line never fitted".
    right = v[v > mu]
    if len(right) < 100:
        continue
    sigma = (right.mean() - mu) * np.sqrt(np.pi / 2.0)
    if np.isfinite(sigma) and sigma > 0:
        h_mu[i], h_sd[i] = mu, sigma
hart_fitted = float(np.isfinite(h_mu).mean())
print(f"  mode-constrained Hart mu : median {np.nanmedian(h_mu):+.3f} "
      f"(v1, unconstrained, gave +4.134 -- the expressed peak)")
print(f"  mode-constrained Hart sd : median {np.nanmedian(h_sd):.3f}")
print(f"  lines successfully fitted: {hart_fitted:.3f}")

# --- continuous per-line score, for pAUC ----------------------------------
# NOT a z-score. sigma of the reference is MAD-based and collapses to 0 on every
# line where most reference genes are exact zeros, which made an earlier version
# all-NaN. Distance above the line's OWN detection threshold needs no spread
# estimate, is always defined, and is directly interpretable as "log2 units
# above the point where this line stops calling a gene off".
#
# Computed at EVERY quantile centre, not just one. The centre enters as a
# per-line offset, so it reorders cell lines within a gene and the pAUC genuinely
# depends on it -- anchoring the comparison to a single centre would leave the
# headline number resting on one arbitrary choice, which is the same mistake the
# 2x2 made in v1. `Lf` is kept whole and the offset subtracted per gene inside
# the loop rather than materialising four full matrices.

# --- assemble the rule set -------------------------------------------------
RULES = {}
for q in REF_QUANTILES:
    RULES[f"ref_q{q}"] = Lf > ref_q[q][:, None]
if hart_fitted >= 0.90:
    RULES["hart_constrained"] = Lf >= (h_mu + HART_CUT * h_sd)[:, None]
else:
    print(f"  DROPPING hart_constrained: fitted only {hart_fitted:.1%} of lines. "
          f"An unfitted line yields NaN, and `Lf >= NaN` is False, which would "
          f"masquerade as 'nothing expressed'.")
RULES[f"tpm_ge_{TPM_CUT:g}"] = TPM >= TPM_CUT

# A rule that excludes nothing, or excludes nearly everything, cannot support any
# incremental comparison against decile 1. Flag rather than silently average in.
DEGENERATE = {}
print()
print(f"  {'rule':<22} {'frac EXPRESSED':>15} {'status':>12}")
rule_frac = {}
for k, M_ in RULES.items():
    f = float(np.nanmean(M_.astype(float)))
    rule_frac[k] = f
    deg = (f >= 0.995) or (f <= 0.05)
    DEGENERATE[k] = deg
    print(f"  {k:<22} {f:>15.4f} {'DEGENERATE' if deg else 'ok':>12}")

results["reference_calibration"] = {
    "families": list(REF_PREFIXES), "n_ref_genes": len(ref_ensg),
    "ref_median_log2tpm": ref_med, "all_median_log2tpm": all_med,
    "separation_log2": all_med - ref_med,
    "frac_exact_zero_cells": frac_zero,
    "ref_mu_median": float(np.nanmedian(ref_mu)),
    "ref_sigma_median": float(np.nanmedian(ref_sd)),
    "hart_constrained_mu_median": float(np.nanmedian(h_mu)),
    "v1_unconstrained_mu": 4.1336,
    "cut_log2tpm_by_quantile": {str(q): float(np.median(ref_q[q]))
                                for q in REF_QUANTILES},
    "frac_expressed_by_rule": rule_frac,
    "hart_lines_fitted": hart_fitted,
    "degenerate_rules": {k: bool(v) for k, v in DEGENERATE.items()},
}

RULE_DF = {k: pd.DataFrame(v, index=E.index, columns=E.columns)
           for k, v in RULES.items()}
LFDF = pd.DataFrame(Lf, index=E.index, columns=E.columns)
del X, TPM, L, Lf, REF, REF_f, RULES


# ============================================================ Step 3
print()
print("=" * 78)
print(f"STEP 3 -- labels ({LABEL})")
print("=" * 78)

LABEL_FPR = None
if LABEL == "gdsc":
    gdsc = pd.read_parquet("validation/prepared/gdsc_scored_ready.parquet",
                           columns=["model_id", "target_ensg", "sensitive"])
    gdsc["model_id"] = gdsc.model_id.str.lower()
    gdsc["target_ensg"] = (gdsc.target_ensg.astype("string")
                           .str.split(".").str[0].str.lower())
    sens = gdsc[gdsc.sensitive][["target_ensg", "model_id"]].drop_duplicates()
    pos_by_gene = sens.groupby("target_ensg")["model_id"].apply(set).to_dict()
    target_genes = sorted(pos_by_gene)
    print(f"  GDSC target genes with >=1 sensitive line : {len(target_genes)}")
else:
    import h5py
    cp = Path(args.chronos_path)
    if not cp.exists():
        raise SystemExit(f"{cp} not found. Expected the Chronos release HDF5.")
    with h5py.File(cp, "r") as fh:
        D = fh["data"][:]
        lines_h = [x.decode() for x in fh["dim_0"][:]]
        genes_h = [x.decode() for x in fh["dim_1"][:]]
    keep_i, keep_e = [], []
    for i, c in enumerate(genes_h):
        e = s2e.get(c.split(" (")[0].strip().upper())
        if e:
            keep_i.append(i)
            keep_e.append(e)
    W = pd.DataFrame(D[:, keep_i], index=[m.lower() for m in lines_h],
                     columns=keep_e)
    W = W.T.groupby(level=0).mean().T
    print(f"  real Chronos: {W.shape[0]} lines x {W.shape[1]:,} genes")

    ess = set(pd.read_csv(cp.parent / "ReferenceEssentials.csv").iloc[:, 0]
              .astype(str).str.split(" (", regex=False).str[0].str.upper())
    non = set(pd.read_csv(cp.parent / "ReferenceNonEssentials.csv").iloc[:, 0]
              .astype(str).str.split(" (", regex=False).str[0].str.upper())
    e_ess = {s2e[g] for g in ess if g in s2e} & set(W.columns)
    e_non = {s2e[g] for g in non if g in s2e} & set(W.columns)
    m_ess = float(np.nanmedian(W[list(e_ess)].values))
    m_non = float(np.nanmedian(W[list(e_non)].values))
    LABEL_FPR = float(np.nanmean(W[list(e_non)].values <= DEP_THRESHOLD))
    print(f"  scale check: essentials {m_ess:+.4f} (~-1); "
          f"non-essentials {m_non:+.4f} (~0)")
    if not (m_ess < -0.5 < m_non + 0.5):
        raise SystemExit("Chronos scale check FAILED. Aborting.")
    print(f"  LABEL FPR at {DEP_THRESHOLD}: {LABEL_FPR:.4f} "
          f"({len(e_non)} curated non-essential genes)")
    results["label_fpr"] = LABEL_FPR

    ch = W.stack().rename("essentiality").reset_index()
    ch.columns = ["model_id", "ensg_id", "essentiality"]
    rng0 = np.random.default_rng(SEED)
    all_ch = sorted(ch.ensg_id.unique())
    if args.n_genes < len(all_ch):
        keep = set(rng0.choice(all_ch, size=args.n_genes, replace=False))
        ch = ch[ch.ensg_id.isin(keep)]
    ch["dep"] = ch.essentiality <= DEP_THRESHOLD
    print(f"  sampled to {ch.ensg_id.nunique():,} genes; "
          f"prevalence {ch.dep.mean():.3f}")
    pos_by_gene = ch[ch.dep].groupby("ensg_id")["model_id"].apply(set).to_dict()
    target_genes = sorted(pos_by_gene)
    print(f"  genes with >=1 dependent line             : {len(target_genes)}")
    del W, D

genes = [g for g in target_genes if g in E.columns]
print(f"  ... with an expression column             : {len(genes)}")


# ============================================================ Step 4
print()
print("=" * 78)
print("STEP 4 -- per gene x rule: pAUC, matched exclusion, 2x2, containment")
print("=" * 78)

rows = []
pooled = {k: {c: [0, 0] for c in ("both", "pct_only", "abs_only", "neither")}
          for k in RULE_DF}
base_pool = [0, 0]

for g in genes:
    raw = E[g].values
    idx = E.index
    pos = idx.isin(pos_by_gene[g])
    n = len(idx)
    if n < MIN_LINES or pos.sum() < MIN_POS or (~pos).sum() < MIN_NEG:
        continue
    base = float(pos.mean())
    base_pool[0] += n
    base_pool[1] += int(pos.sum())

    pct = pd.Series(raw).rank(pct=True, method="min").values
    r = rankdata(raw, method="ordinal")
    q10 = np.minimum((r - 1) * 10 // n, 9)
    dec1 = q10 == 0
    n_d1 = int(dec1.sum())

    # Continuous per-line score at every quantile centre. The centre is a
    # per-line offset, so each one induces a different within-gene ordering.
    lf_g = LFDF[g].values
    pauc_pct = partial_auc(pct, pos, 0.8, 1.0)
    rank_pct = partial_auc(pct, pos, 0.0, 0.2)
    auc_pct_g = auroc(pct, pos)
    cont_metrics = {}
    for qq in REF_QUANTILES:
        cont = lf_g - ref_q[qq]
        cont_metrics[f"pauc_gate_cont_q{qq}"] = partial_auc(cont, pos, 0.8, 1.0)
        cont_metrics[f"pauc_rank_cont_q{qq}"] = partial_auc(cont, pos, 0.0, 0.2)
        cont_metrics[f"auc_cont_q{qq}"] = auroc(cont, pos)

    for rname, RM in RULE_DF.items():
        not_expr = ~RM[g].values
        n_ne = int(not_expr.sum())
        overlap = int((dec1 & not_expr).sum())
        rec = {"ensg_id": g, "symbol": symbol_of.get(g, "?"), "rule": rname,
               "n_lines": n, "n_pos": int(pos.sum()), "base_rate": base,
               "n_d1": n_d1, "n_not_expressed": n_ne,
               "frac_not_expressed": n_ne / n,
               "overlap": overlap,
               # containment: how much of each set the other swallows
               "contain_d1_in_ne": overlap / n_d1 if n_d1 else np.nan,
               "contain_ne_in_d1": overlap / n_ne if n_ne else np.nan,
               "cell_pct_only": n_d1 - overlap,
               "cell_abs_only": n_ne - overlap,
               "paired_possible": (n_d1 - overlap >= MIN_CELL
                                   and n_ne - overlap >= MIN_CELL),
               "n_paralogs": par_n.get(g, 0),
               "has_paralog": par_n.get(g, 0) >= 1,
               "n_paralogs_strict": par_n_strict.get(g, 0),
               "has_paralog_strict": par_n_strict.get(g, 0) >= 1,
               "pauc_gate_pct": pauc_pct, "pauc_rank_pct": rank_pct,
               "auc_pct": auc_pct_g, **cont_metrics}

        # matched exclusion
        if 0 < n_ne < n:
            order = np.argsort(pct, kind="mergesort")[:n_ne]
            ra = float(pos[not_expr].mean())
            rp = float(pos[order].mean())
            rec.update({"excl_rate_abs": ra, "excl_rate_pct": rp,
                        "excl_ratio_abs": ra / base, "excl_ratio_pct": rp / base,
                        "excl_delta": rp - ra,
                        "npv_base": 1 - base, "npv_abs": 1 - ra, "npv_pct": 1 - rp,
                        "excl_p_abs": float(hypergeom.cdf(
                            int(pos[not_expr].sum()), n, int(pos.sum()), n_ne)),
                        "excl_p_pct": float(hypergeom.cdf(
                            int(pos[order].sum()), n, int(pos.sum()), n_ne))})

        # 2x2
        for key, mask in (("both", dec1 & not_expr),
                          ("pct_only", dec1 & ~not_expr),
                          ("abs_only", ~dec1 & not_expr),
                          ("neither", ~dec1 & ~not_expr)):
            pooled[rname][key][0] += int(mask.sum())
            pooled[rname][key][1] += int(pos[mask].sum())
            if mask.sum() >= MIN_CELL:
                rt = float(pos[mask].mean())
                rec[f"rate_{key}"] = rt
                rec[f"ratio_{key}"] = rt / base if base else np.nan

        # floor
        if n_ne >= 20:
            rec["floor_abs_rate"] = float(pos[not_expr].mean())
            rec["floor_abs_ratio"] = rec["floor_abs_rate"] / base
        fl_min = raw == raw.min()
        if fl_min.sum() >= 20:
            rec["floor_min_rate"] = float(pos[fl_min].mean())
            rec["floor_min_ratio"] = rec["floor_min_rate"] / base
        rows.append(rec)

G = pd.DataFrame(rows)
G.to_parquet(OUTPUTS / f"test_run_detection_vs_percentile_{LABEL}_per_gene.parquet")
n_scored = G.ensg_id.nunique()
base_all = base_pool[1] / base_pool[0]
print(f"  genes scored : {n_scored}   rules : {len(RULE_DF)}   "
      f"rows : {len(G):,}")
print(f"  median lines/gene {G.n_lines.median():.0f}   "
      f"median per-gene prevalence {G.groupby('ensg_id').base_rate.first().median():.4f}   "
      f"pooled {base_all:.4f}")


# ============================================================ Step 5
print()
print("=" * 78)
print("STEP 5 -- GATE-REGION pAUC.  percentile vs per-line detection distance")
print("=" * 78)
one = G.drop_duplicates("ensg_id")
print("  (percentile == raw log2TPM here, exactly, by construction, so the")
print("   percentile row is fixed and only the per-line centre varies)")
print()
print(f"  {'score':<34} {'gate pAUC':>10} {'rank pAUC':>10} {'full AUC':>10}")
print(f"  {'within-gene percentile':<34} {one.pauc_gate_pct.median():>10.4f} "
      f"{one.pauc_rank_pct.median():>10.4f} {one.auc_pct.median():>10.4f}")
for qq in REF_QUANTILES:
    print(f"  {'per-line dist, centre Q' + str(qq):<34} "
          f"{one[f'pauc_gate_cont_q{qq}'].median():>10.4f} "
          f"{one[f'pauc_rank_cont_q{qq}'].median():>10.4f} "
          f"{one[f'auc_cont_q{qq}'].median():>10.4f}")

print()
print(f"  PAIRED DELTA (percentile - per-line), gate region, per centre")
print(f"  negative = the per-line absolute score WINS")
print(f"  {'centre':<12} {'delta':>9} {'95% CI':>20} {'p':>11} "
      f"{'pct wins':>9} {'floor':>6}")
pauc_by_q = {}
for qq in REF_QUANTILES:
    d = (one.pauc_gate_pct - one[f"pauc_gate_cont_q{qq}"]).dropna()
    if not len(d):
        continue
    lo_q, hi_q = boot_ci(d)
    p_q = wilx(d)
    clears = bool(np.isfinite(d.median()) and abs(d.median()) >= MIN_PAUC_DELTA
                  and np.isfinite(p_q) and p_q < 0.05)
    print(f"  Q{qq:<11} {d.median():>+9.4f} "
          f"[{lo_q:>+8.4f},{hi_q:>+8.4f}] {p_q:>11.3g} "
          f"{100*(d > 0).mean():>8.1f}% {'YES' if clears else 'no':>6}")
    pauc_by_q[str(qq)] = {
        "gate_pct": float(one.pauc_gate_pct.median()),
        "gate_cont": float(one[f"pauc_gate_cont_q{qq}"].median()),
        "rank_cont": float(one[f"pauc_rank_cont_q{qq}"].median()),
        "paired_delta": float(d.median()), "paired_delta_ci": [lo_q, hi_q],
        "paired_p": p_q, "pct_better_frac": float((d > 0).mean()),
        "clears_effect_floor": clears,
    }

deltas = [v["paired_delta"] for v in pauc_by_q.values()]
n_clear = sum(v["clears_effect_floor"] for v in pauc_by_q.values())
n_neg = sum(d < 0 for d in deltas)
spread = max(deltas) - min(deltas) if deltas else np.nan
print()
print(f"  delta across the {len(deltas)} centres : {min(deltas):+.4f} to "
      f"{max(deltas):+.4f}   spread {spread:.4f}")
print(f"  clears the effect floor at      : {n_clear}/{len(deltas)} centres")
print(f"  per-line score wins at          : {n_neg}/{len(deltas)} centres")
stable = bool(spread < MIN_PAUC_DELTA and n_neg == len(deltas))
print(f"  direction STABLE across centres : {'YES' if stable else 'NO'}")
results["pauc_by_centre"] = pauc_by_q
results["pauc_centre_stability"] = {
    "delta_min": min(deltas) if deltas else None,
    "delta_max": max(deltas) if deltas else None,
    "delta_spread": float(spread) if deltas else None,
    "n_centres": len(deltas), "n_clearing_floor": n_clear,
    "n_favouring_per_line": n_neg, "direction_stable": stable,
}
# Conservative pair kept for the verdict: percentile against the centre where the
# per-line score does WORST, so a per-line win here holds at every centre.
worst_cont = one[[f"pauc_gate_cont_q{q}" for q in REF_QUANTILES]].min(axis=1)
dg = (one.pauc_gate_pct - worst_cont).dropna()
p_gate = wilx(dg)
print(f"  conservative delta (vs the per-line score's worst centre) : "
      f"{dg.median():+.4f}  p = {p_gate:.3g}")


# ============================================================ Step 6
print()
print("=" * 78)
print("STEP 5b -- HALF-SAMPLE HELD-OUT CENTRE SELECTION")
print("=" * 78)
print(f"  {HALF_SAMPLE_REPS} random splits of the {n_lines} cell lines.")
print("  select the centre on half A by median gate pAUC, evaluate on half B.")
print("  the percentile is RECOMPUTED within half B -- it is panel-dependent, so")
print("  giving it the full-panel ranking on a half-panel evaluation would flatter it.")
print()

# Pre-extract per-gene arrays once; the rep loop is O(reps x genes).
gene_arrays = []
for g in G.ensg_id.unique():
    gene_arrays.append((g, E[g].values, LFDF[g].values,
                        E.index.isin(pos_by_gene[g])))

rng = np.random.default_rng(SEED)
half = n_lines // 2
sel_counter = {q: 0 for q in REF_QUANTILES}
rep_rows = []
for rep in range(HALF_SAMPLE_REPS):
    perm = rng.permutation(n_lines)
    A, B = perm[:half], perm[half:]

    # --- select the centre on half A -------------------------------------
    a_med = {}
    for qq in REF_QUANTILES:
        vals = []
        for g, raw, lf, pos in gene_arrays:
            pa, la = pos[A], lf[A]
            if pa.sum() < MIN_HALF_POS or (~pa).sum() < MIN_HALF_NEG:
                continue
            vals.append(partial_auc(la - ref_q[qq][A], pa, 0.8, 1.0))
        a_med[qq] = float(np.nanmedian(vals)) if vals else np.nan
    q_star = max(REF_QUANTILES, key=lambda q: (a_med[q] if np.isfinite(a_med[q])
                                               else -np.inf))
    sel_counter[q_star] += 1

    # --- evaluate the FIXED choice on half B ------------------------------
    b_pct, b_cont = [], []
    for g, raw, lf, pos in gene_arrays:
        pb = pos[B]
        if pb.sum() < MIN_HALF_POS or (~pb).sum() < MIN_HALF_NEG:
            continue
        # percentile recomputed on half B only
        pct_b = pd.Series(raw[B]).rank(pct=True, method="min").values
        b_pct.append(partial_auc(pct_b, pb, 0.8, 1.0))
        b_cont.append(partial_auc(lf[B] - ref_q[q_star][B], pb, 0.8, 1.0))
    bp = np.array(b_pct, dtype=float)
    bc = np.array(b_cont, dtype=float)
    ok = np.isfinite(bp) & np.isfinite(bc)
    rep_rows.append({
        "rep": rep, "selected_centre": q_star,
        "selection_pauc_A": a_med[q_star],
        "pauc_B_pct": float(np.nanmedian(bp[ok])),
        "pauc_B_cont": float(np.nanmedian(bc[ok])),
        "delta_B": float(np.nanmedian(bp[ok] - bc[ok])),
        "n_genes_B": int(ok.sum()),
    })

RS = pd.DataFrame(rep_rows)
print(f"  centre selected on half A : "
      + ", ".join(f"Q{q}: {c}/{HALF_SAMPLE_REPS}"
                  for q, c in sel_counter.items() if c))
print()
print(f"  {'quantity':<44} {'median':>9} {'min':>9} {'max':>9}")
print(f"  {'selection pAUC on half A (the winner)':<44} "
      f"{RS.selection_pauc_A.median():>9.4f} {RS.selection_pauc_A.min():>9.4f} "
      f"{RS.selection_pauc_A.max():>9.4f}")
print(f"  {'HELD-OUT pAUC on half B, per-line score':<44} "
      f"{RS.pauc_B_cont.median():>9.4f} {RS.pauc_B_cont.min():>9.4f} "
      f"{RS.pauc_B_cont.max():>9.4f}")
print(f"  {'HELD-OUT pAUC on half B, percentile':<44} "
      f"{RS.pauc_B_pct.median():>9.4f} {RS.pauc_B_pct.min():>9.4f} "
      f"{RS.pauc_B_pct.max():>9.4f}")
print(f"  {'held-out delta (percentile - per-line)':<44} "
      f"{RS.delta_B.median():>+9.4f} {RS.delta_B.min():>+9.4f} "
      f"{RS.delta_B.max():>+9.4f}")
shrink = RS.selection_pauc_A.median() - RS.pauc_B_cont.median()
per_line_wins_ho = int((RS.delta_B < 0).sum())
ho_clears = bool(abs(RS.delta_B.median()) >= MIN_PAUC_DELTA)
print()
print(f"  OPTIMISM (selection pAUC on A - held-out pAUC on B) : {shrink:+.4f}")
print(f"  per-line score wins on held-out data in             : "
      f"{per_line_wins_ho}/{HALF_SAMPLE_REPS} reps")
print(f"  held-out delta clears the effect floor              : "
      f"{'YES' if ho_clears else 'NO'}")
results["half_sample"] = {
    "reps": HALF_SAMPLE_REPS, "lines_per_half": int(half),
    "centre_selection_counts": {str(q): c for q, c in sel_counter.items()},
    "selection_pauc_A_median": float(RS.selection_pauc_A.median()),
    "heldout_pauc_B_cont_median": float(RS.pauc_B_cont.median()),
    "heldout_pauc_B_pct_median": float(RS.pauc_B_pct.median()),
    "heldout_delta_median": float(RS.delta_B.median()),
    "heldout_delta_min": float(RS.delta_B.min()),
    "heldout_delta_max": float(RS.delta_B.max()),
    "optimism": float(shrink),
    "per_line_wins_reps": per_line_wins_ho,
    "heldout_clears_effect_floor": ho_clears,
}


print()
print("=" * 78)
print("STEP 6 -- PER-GENE 2x2 ACROSS THE SWEEP  (not anchored to one cut)")
print("=" * 78)
print("  ratios are PER-GENE MEDIANS of (cell dependency rate / gene base rate).")
print("  the deciding cell is 'dec1 & EXPRESSED' -- the percentile's increment")
print("  over detection. 1.00 = no depletion = the percentile adds nothing.")
print()
hdr = (f"  {'rule':<18} {'%excl':>6} {'both':>18} {'dec1&EXPR':>18} "
       f"{'NE&not-dec1':>18} {'nPair':>6}")
print(hdr)
sweep_out = {}
for rname in RULE_DF:
    s = G[G.rule == rname]
    def cell(col):
        v = s[col].dropna() if col in s.columns else pd.Series(dtype=float)
        if len(v) < 6:
            return np.nan, np.nan, 0
        return float(v.median()), wilx(v, 1.0), len(v)
    m_b, p_b, n_b = cell("ratio_both")
    m_p, p_p, n_p = cell("ratio_pct_only")
    m_a, p_a, n_a = cell("ratio_abs_only")
    npair = int(s.paired_possible.sum())
    print(f"  {rname:<18} {s.frac_not_expressed.median()*100:>5.1f}% "
          f"{m_b:>8.3f} (n={n_b:<4}) {m_p:>8.3f} (n={n_p:<4}) "
          f"{m_a:>8.3f} (n={n_a:<4}) {npair:>6}")
    sweep_out[rname] = {
        "median_frac_excluded": float(s.frac_not_expressed.median()),
        "ratio_both": m_b, "p_both": p_b, "n_both": n_b,
        "ratio_pct_only": m_p, "p_pct_only": p_p, "n_pct_only": n_p,
        "ratio_abs_only": m_a, "p_abs_only": p_a, "n_abs_only": n_a,
        "n_paired_possible": npair,
        "pct_only_clears_floor": bool(np.isfinite(m_p)
                                      and abs(1 - m_p) >= MIN_RATIO_EFFECT
                                      and np.isfinite(p_p) and p_p < 0.05),
        "abs_only_clears_floor": bool(np.isfinite(m_a)
                                      and abs(1 - m_a) >= MIN_RATIO_EFFECT
                                      and np.isfinite(p_a) and p_a < 0.05),
    }
print()
print(f"  effect floor: |1 - ratio| >= {MIN_RATIO_EFFECT} AND p < 0.05")
for rname, v in sweep_out.items():
    print(f"    {rname:<18} percentile increment real: "
          f"{'YES' if v['pct_only_clears_floor'] else 'no ':<4}   "
          f"detection increment real: "
          f"{'YES' if v['abs_only_clears_floor'] else 'no'}")
results["sweep_per_gene_2x2"] = sweep_out

# pooled, reported only to show it disagrees and why
print()
print("  POOLED OVER PAIRS (shown to document the Simpson reversal, not to read)")
pool_out = {}
for rname in RULE_DF:
    row = {}
    for k in ("both", "pct_only", "abs_only", "neither"):
        npair, npos = pooled[rname][k]
        row[k] = {"pairs": npair,
                  "dep_rate": npos / npair if npair else np.nan,
                  "ratio": (npos / npair) / base_all if npair else np.nan}
    pool_out[rname] = row
    print(f"    {rname:<18} dec1&EXPR ratio {row['pct_only']['ratio']:>6.3f}  "
          f"NE&not-dec1 ratio {row['abs_only']['ratio']:>6.3f}")
results["sweep_pooled_2x2"] = pool_out
results["sweep_pooled_base_rate"] = base_all


# ============================================================ Step 7
print()
print("=" * 78)
print("STEP 7 -- WHY n IS SMALL, AND WHETHER THE RULES ARE STRUCTURALLY DISJOINT")
print("=" * 78)
print("  decile 1 is exactly 10% of a gene's lines BY CONSTRUCTION. So:")
print("    rule excludes << 10%  ->  its exclusions sit inside decile 1")
print("                          ->  'NE & not-dec1' is empty")
print("    rule excludes >> 10%  ->  decile 1 sits inside its exclusions")
print("                          ->  'dec1 & EXPRESSED' is empty")
print("  Both cells populate only where the exclusion fraction is NEAR 10%.")
print()
print(f"  {'rule':<18} {'excl=0':>7} {'0-5%':>7} {'5-15%':>7} {'>15%':>7} "
      f"{'d1 in NE':>9} {'NE in d1':>9} {'paired':>7}")
struct = {}
for rname in RULE_DF:
    s = G[G.rule == rname]
    f = s.frac_not_expressed
    b0 = int((f == 0).sum())
    b1 = int(((f > 0) & (f < 0.05)).sum())
    b2 = int(((f >= 0.05) & (f <= 0.15)).sum())
    b3 = int((f > 0.15).sum())
    print(f"  {rname:<18} {b0:>7} {b1:>7} {b2:>7} {b3:>7} "
          f"{s.contain_d1_in_ne.median():>9.3f} {s.contain_ne_in_d1.median():>9.3f} "
          f"{int(s.paired_possible.sum()):>7}")
    struct[rname] = {"n_genes": int(len(s)), "excl_zero": b0, "excl_0_5pct": b1,
                     "excl_5_15pct": b2, "excl_over_15pct": b3,
                     "median_contain_d1_in_ne": float(s.contain_d1_in_ne.median()),
                     "median_contain_ne_in_d1": float(s.contain_ne_in_d1.median()),
                     "n_paired_possible": int(s.paired_possible.sum()),
                     "frac_paired_possible": float(s.paired_possible.mean())}
results["structural_overlap"] = struct
print()
print("  d1 in NE = fraction of decile-1 lines the rule also calls not-expressed")
print("  NE in d1 = fraction of the rule's exclusions that are inside decile 1")
print("  paired   = genes where BOTH incremental cells have >= "
      f"{MIN_CELL} lines")


# ============================================================ Step 8
print()
print("=" * 78)
print("STEP 7b -- PARALOG STRATIFICATION of the deciding cell")
print("=" * 78)
one_g = G.drop_duplicates("ensg_id")
n_par = int(one_g.has_paralog.sum())
print(f"  genes with >=1 HGNC-group paralog : {n_par} of {len(one_g)} "
      f"({100*n_par/len(one_g):.1f}%)   strict cap: "
      f"{int(one_g.has_paralog_strict.sum())}")

# --- validate the annotation before interpreting it ------------------------
# Paralog buffering predicts LOWER dependency prevalence among paralogous genes.
pv_par = one_g[one_g.has_paralog].base_rate
pv_sing = one_g[~one_g.has_paralog].base_rate
print()
print(f"  VALIDATION -- paralog buffering predicts lower dependency prevalence")
print(f"    prevalence, genes WITH paralogs    : {pv_par.median():.4f} "
      f"(n={len(pv_par)})")
print(f"    prevalence, genes WITHOUT paralogs : {pv_sing.median():.4f} "
      f"(n={len(pv_sing)})")
try:
    from scipy.stats import mannwhitneyu
    mw = mannwhitneyu(pv_par.dropna(), pv_sing.dropna(), alternative="less")
    mw_p = float(mw.pvalue)
except Exception:
    mw_p = np.nan
print(f"    Mann-Whitney (paralog < singleton) : p = {mw_p:.4g}")
annot_valid = bool(np.isfinite(mw_p) and mw_p < 0.05
                   and pv_par.median() < pv_sing.median())
print(f"    -> annotation behaves like paralogy : "
      f"{'YES' if annot_valid else 'NO -- stratification NOT interpreted'}")
results["paralog_validation"] = {
    "prevalence_with_paralog": float(pv_par.median()),
    "prevalence_without_paralog": float(pv_sing.median()),
    "n_with": int(len(pv_par)), "n_without": int(len(pv_sing)),
    "mannwhitney_p": mw_p, "annotation_behaves_like_paralogy": annot_valid,
}

print()
print(f"  DECIDING CELL ('dec1 & EXPRESSED') BY PARALOG STRATUM")
print(f"  {'rule':<18} {'paralog':>9} {'n':>5} {'singleton':>10} {'n':>5} "
       f"{'delta':>8} {'p(strat)':>9}")
par_out = {}
for rname in RULE_DF:
    s = G[G.rule == rname]
    a = s[s.has_paralog].ratio_pct_only.dropna()
    b = s[~s.has_paralog].ratio_pct_only.dropna()
    if len(a) < 6 or len(b) < 6:
        print(f"  {rname:<18}   (too few genes in a stratum)")
        continue
    try:
        sp = float(mannwhitneyu(a, b, alternative="two-sided").pvalue)
    except Exception:
        sp = np.nan
    print(f"  {rname:<18} {a.median():>9.3f} {len(a):>5} {b.median():>10.3f} "
          f"{len(b):>5} {a.median()-b.median():>+8.3f} {sp:>9.3g}")
    par_out[rname] = {
        "ratio_pct_only_paralog": float(a.median()), "n_paralog": int(len(a)),
        "ratio_pct_only_singleton": float(b.median()), "n_singleton": int(len(b)),
        "delta": float(a.median() - b.median()), "strat_p": sp,
        # the headline finding recomputed on SINGLETONS ONLY, i.e. with paralog
        # buffering removed rather than merely adjusted for
        "singleton_clears_floor": bool(abs(1 - b.median()) >= MIN_RATIO_EFFECT
                                       and np.isfinite(wilx(b, 1.0))
                                       and wilx(b, 1.0) < 0.05),
    }
results["paralog_stratified_deciding_cell"] = par_out
sing_clear = sum(v["singleton_clears_floor"] for v in par_out.values())
print()
print(f"  the 0.81 finding recomputed on SINGLETONS ONLY clears the effect")
print(f"  floor in {sing_clear}/{len(par_out)} rules "
      f"(paralog buffering excluded, not adjusted)")
results["paralog_singleton_rules_clearing"] = sing_clear
results["paralog_n_rules"] = len(par_out)


print()
print("=" * 78)
print("STEP 8 -- MATCHED EXCLUSION per rule  (equal number of lines removed)")
print("=" * 78)
print(f"  {'rule':<18} {'genes':>6} {'absolute':>10} {'percentile':>11} "
      f"{'base':>8} {'delta':>9} {'p':>9}")
match_out = {}
for rname in RULE_DF:
    s = G[(G.rule == rname)].dropna(subset=["excl_rate_abs", "excl_rate_pct"])
    if len(s) < 6:
        print(f"  {rname:<18} {len(s):>6}   (too few genes)")
        match_out[rname] = {"n_genes": int(len(s))}
        continue
    d = (s.excl_rate_pct - s.excl_rate_abs).dropna()
    pv = wilx(d)
    print(f"  {rname:<18} {len(s):>6} {s.excl_rate_abs.median():>10.4f} "
          f"{s.excl_rate_pct.median():>11.4f} {s.base_rate.median():>8.4f} "
          f"{d.median():>+9.5f} {pv:>9.3g}")
    match_out[rname] = {
        "n_genes": int(len(s)),
        "dep_rate_absolute": float(s.excl_rate_abs.median()),
        "dep_rate_percentile": float(s.excl_rate_pct.median()),
        "ratio_absolute": float(s.excl_ratio_abs.median()),
        "ratio_percentile": float(s.excl_ratio_pct.median()),
        "base_rate": float(s.base_rate.median()),
        "paired_delta": float(d.median()), "paired_p": pv,
        "pct_cleaner_frac": float((d < 0).mean()),
        "npv_absolute": float(s.npv_abs.median()),
        "npv_percentile": float(s.npv_pct.median()),
        "npv_base": float(s.npv_base.median()),
    }
print("  delta = percentile - absolute; negative = percentile excludes cleaner")
results["matched_exclusion_by_rule"] = match_out


# ============================================================ Step 9
print()
print("=" * 78)
print("STEP 9 -- FLOOR DEPLETION, judged against the label's own FPR")
print("=" * 78)
fm = G.drop_duplicates("ensg_id").dropna(subset=["floor_min_ratio"])
print(f"  {'definition':<30} {'genes':>6} {'rate':>9} {'ratio':>7} {'resid':>8}")
if len(fm):
    fr = float(fm.floor_min_rate.median())
    print(f"  {'exact minimum (current)':<30} {len(fm):>6} {fr:>9.4f} "
          f"{fm.floor_min_ratio.median():>7.3f} "
          f"{(max(0, fr-LABEL_FPR) if LABEL_FPR else np.nan):>8.4f}")
floor_out = {"exact_min": ({"n_genes": int(len(fm)),
                            "rate": float(fm.floor_min_rate.median()),
                            "ratio": float(fm.floor_min_ratio.median())}
                           if len(fm) else None)}
for rname in RULE_DF:
    s = G[G.rule == rname].dropna(subset=["floor_abs_ratio"])
    if not len(s):
        continue
    fr = float(s.floor_abs_rate.median())
    print(f"  {rname:<30} {len(s):>6} {fr:>9.4f} "
          f"{s.floor_abs_ratio.median():>7.3f} "
          f"{(max(0, fr-LABEL_FPR) if LABEL_FPR else np.nan):>8.4f}")
    floor_out[rname] = {"n_genes": int(len(s)), "rate": fr,
                        "ratio": float(s.floor_abs_ratio.median()),
                        "residual_above_label_fpr":
                            (max(0.0, fr - LABEL_FPR) if LABEL_FPR else None)}
floor_out["label_fpr"] = LABEL_FPR
results["floor"] = floor_out
if LABEL_FPR:
    print(f"  label FPR (curated non-essentials) = {LABEL_FPR:.4f}; a rate at or "
          f"below it is\n  complete depletion to within what this label can "
          f"resolve.")


# ============================================================ Step 10
print()
print("=" * 78)
print("STEP 10 -- SWEEP INSENSITIVITY")
print("=" * 78)
# Degenerate rules are excluded: a rule that calls everything expressed (or
# nothing) has no incremental cell to measure, and averaging its NaN-driven
# ratio into the spread would manufacture sensitivity that is not there.
valid = [k for k in sweep_out if not DEGENERATE.get(k)]
print(f"  rules in the sweep : {len(sweep_out)}   "
      f"non-degenerate : {len(valid)}  ({', '.join(valid)})")
if set(sweep_out) - set(valid):
    print(f"  excluded as degenerate : "
          f"{', '.join(sorted(set(sweep_out) - set(valid)))}")
print()
ins = {}
for cellname, lbl in (("ratio_both", "'both'"),
                      ("ratio_pct_only", "'dec1 & EXPRESSED'"),
                      ("ratio_abs_only", "'NE & not-dec1'")):
    vals = [sweep_out[k][cellname] for k in valid
            if np.isfinite(sweep_out[k][cellname])]
    if not vals:
        continue
    spread = max(vals) - min(vals)
    ins[cellname] = {"min": min(vals), "max": max(vals), "spread": spread,
                     "insensitive": bool(spread < 0.15)}
    print(f"  {lbl:<22} across {len(vals)} rules : {min(vals):.3f} - "
          f"{max(vals):.3f}   spread {spread:.3f}   "
          f"{'INSENSITIVE' if spread < 0.15 else 'sensitive'}")
key = ins.get("ratio_pct_only", {})
insensitive = bool(key.get("insensitive"))
print()
print(f"  the deciding cell ('dec1 & EXPRESSED') is "
      f"{'INSENSITIVE' if insensitive else 'SENSITIVE'} to the threshold")
results["sweep_insensitivity"] = {"by_cell": ins, "valid_rules": valid,
                                  "degenerate_excluded":
                                      sorted(set(sweep_out) - set(valid)),
                                  "deciding_cell_insensitive": insensitive}


# ============================================================ Step 11
print()
print("=" * 78)
print("STEP 11 -- VERDICT  (significance AND effect size required)")
print("=" * 78)

pct_real = sum(sweep_out[k]["pct_only_clears_floor"] for k in valid)
abs_real = sum(sweep_out[k]["abs_only_clears_floor"] for k in valid)
n_rules = len(valid)
# The pAUC difference counts as real only if it clears the effect floor at EVERY
# quantile centre -- a difference that appears at one centre and vanishes at
# another is a property of the centre, not of the score.
pauc_real = bool(n_clear == len(deltas) and len(deltas) > 0)
per_line_wins = bool(n_neg == len(deltas) and pauc_real)
paired_ok = max(sweep_out[k]["n_paired_possible"] for k in valid) >= 30

print(f"  non-degenerate rules assessed            : {n_rules}")
print(f"  percentile increment clears the floor in {pct_real}/{n_rules} rules")
print(f"  detection increment clears the floor in  {abs_real}/{n_rules} rules")
print(f"  pAUC gap clears the floor at             : {n_clear}/{len(deltas)} centres")
print(f"  per-line score wins at                   : {n_neg}/{len(deltas)} centres")
print(f"  gate-region pAUC difference is real      : {pauc_real}")
print(f"  any rule supports a paired test (n>=30)  : {paired_ok}")
print(f"  HELD-OUT: per-line wins in               : "
      f"{per_line_wins_ho}/{HALF_SAMPLE_REPS} reps, "
      f"delta {RS.delta_B.median():+.4f}, clears floor: {ho_clears}")
print(f"  PARALOG: singleton-only finding holds in : "
      f"{sing_clear}/{len(par_out)} rules "
      f"(annotation valid: {annot_valid})")
results["per_line_beats_percentile_at_all_centres"] = per_line_wins

if pct_real >= n_rules / 2 and pauc_real:
    v = "percentile_carries_graded_information_beyond_detection"
elif pct_real >= n_rules / 2:
    v = "percentile_increment_survives_effect_floor_but_pauc_is_tied"
elif abs_real >= n_rules / 2 and pct_real == 0:
    v = "detection_only__percentile_is_an_encoding_of_a_boolean"
elif not paired_ok:
    v = "rules_structurally_disjoint__not_directly_comparable"
else:
    v = "both_increments_weak__neither_rule_is_load_bearing"
results["verdict"] = v
print(f"\n  VERDICT: {v}")

out = OUTPUTS / f"test_run_detection_vs_percentile_{LABEL}_results.json"
out.write_text(json.dumps(results, indent=2, default=float), encoding="utf-8")
print(f"\nwrote {out}")
print(f"wrote {OUTPUTS / f'test_run_detection_vs_percentile_{LABEL}_per_gene.parquet'}")
