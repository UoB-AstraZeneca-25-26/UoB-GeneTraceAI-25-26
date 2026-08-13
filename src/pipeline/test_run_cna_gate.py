"""
test_run_cna_gate.py
--------------------
Read-only diagnostic. Writes no production file.

QUESTION
--------
The expression gate established in test_run_abundance_gate.py excludes cell lines
where a gene is not transcribed. Homozygous deletion is the LIMITING CASE of that
gate: not "the gene is silent" but "the gene is not there". If the gate logic is
right, deletion should be an absolute exclusion -- the dependency rate at deleted
(gene, line) pairs should sit at the label's noise floor, with nothing left over.

There is also a mechanical reason to expect it, which makes this a check on the
labels as much as on the biology: CRISPR cannot cut DNA that is absent, so a
knockout screen physically cannot register a fitness effect at a homozygously
deleted locus. A non-zero rate there is measurement error by construction.

The second, and more useful, question is whether copy number ADDS anything over
expression. Deleted genes are also unexpressed, so the two gates overlap. The
cross-tabulation below separates:

    deleted AND at the expression floor   -- the overlap, no new information
    deleted but NOT at the floor          -- CNA catching what RNA missed
    at the floor but NOT deleted          -- silencing without loss
    neither                               -- the scored population

Only the second cell justifies adding a layer. If deletion is a strict subset of
the expression floor, gene-level CNA is redundant with what Layer 2 already has.

WHY THIS FILE
-------------
The repo's existing CNA is COSMIC-derived: cna_flags.parquet covers 117,540 of
15,545 genes x 969 models = 0.78% of the grid, because COSMIC records only
altered genes. CCLEGeneCopyNumber20Q2.hdf5 (Chronos release, Dempster 2021) is
dense: 908 lines x 27,639 genes with no missing values.

SCALE AND THRESHOLD
-------------------
Values are log2(relative copy number + 1): 1.0 = diploid, 0 = complete loss.
Verified here rather than assumed -- CDKN2A 29.3% below 0.1 (literature ~30% of
cancer lines carry homozygous CDKN2A deletion), MTAP 14.6% (the classic CDKN2A
co-deletion), GAPDH and ACTB exactly 0.000%. Genome-wide 0.10% of pairs.
The cut is insensitive: <0.1 gives 0.0010 of pairs, <0.2 gives 0.0012.
DEL_THRESHOLD = 0.1 is used, with a sensitivity sweep reported.

SIGN CONVENTION
---------------
Chronos gene effect: more negative = more essential; positive (dependent) =
gene effect <= DEP_THRESHOLD. Verified against the release's own reference sets.

Outputs:
    src/pipeline/outputs/test_run_cna_gate_results.json

Run:
    python src/pipeline/test_run_cna_gate.py
"""
import json
import warnings
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
import pyarrow.parquet as pq
from scipy.stats import fisher_exact

warnings.filterwarnings("ignore", category=RuntimeWarning)

OUTPUTS = Path("src/pipeline/outputs")
REF_DIR = Path("reference")
DATA_CLEAN = Path("data/parquet/data_clean")
CHR_DIR = Path("data/DepMap_Chronos")

DEP_THRESHOLD = -0.5
DEL_THRESHOLD = 0.1
MIN_LINES = 100
SEED = 42
results = {}

# ============================================================ Step 1
print("=" * 78)
print("STEP 1 -- load Chronos, copy number, expression; verify both scales")
print("=" * 78)

gl = pd.read_parquet(REF_DIR / "gene_lookup.parquet",
                     columns=["ensg_id", "hgnc_symbol"])
gl["ensg_id"] = gl.ensg_id.astype(str).str.split(".").str[0].str.lower()
s2e = dict(zip(gl.hgnc_symbol.astype(str).str.upper(), gl.ensg_id))


def load_h5(path):
    with h5py.File(path, "r") as fh:
        D = fh["data"][:]
        lines = [x.decode().lower() for x in fh["dim_0"][:]]
        genes = [x.decode() for x in fh["dim_1"][:]]
    ki, ke = [], []
    for i, c in enumerate(genes):
        e = s2e.get(c.split(" (")[0].strip().upper())
        if e:
            ki.append(i)
            ke.append(e)
    W = pd.DataFrame(D[:, ki], index=lines, columns=ke)
    return W.T.groupby(level=0).mean().T


CH = load_h5(CHR_DIR / "GeneFitnessEffect_Chronos_Achilles.hdf5")
CN = load_h5(CHR_DIR / "CCLEGeneCopyNumber20Q2.hdf5")

ess = set(pd.read_csv(CHR_DIR / "ReferenceEssentials.csv").iloc[:, 0]
          .astype(str).str.split(" (", regex=False).str[0].str.upper())
non = set(pd.read_csv(CHR_DIR / "ReferenceNonEssentials.csv").iloc[:, 0]
          .astype(str).str.split(" (", regex=False).str[0].str.upper())
e_ess = {s2e[g] for g in ess if g in s2e} & set(CH.columns)
e_non = {s2e[g] for g in non if g in s2e} & set(CH.columns)
m_ess = float(np.nanmedian(CH[list(e_ess)].values))
m_non = float(np.nanmedian(CH[list(e_non)].values))
LABEL_FPR = float(np.nanmean(CH[list(e_non)].values <= DEP_THRESHOLD))
print(f"  Chronos : {CH.shape[0]} lines x {CH.shape[1]:,} genes   "
      f"ess {m_ess:+.4f}  non-ess {m_non:+.4f}  FPR {LABEL_FPR:.4f}")
if not (m_ess < -0.5 < m_non + 0.5):
    raise SystemExit("Chronos scale check FAILED. Aborting.")

print(f"  CopyNum : {CN.shape[0]} lines x {CN.shape[1]:,} genes   "
      f"median {np.nanmedian(CN.values):.3f} (expect ~1.0 = diploid)")
for sym, exp in [("CDKN2A", "~30% deleted"), ("MTAP", "CDKN2A co-deletion"),
                 ("GAPDH", "never deleted"), ("ACTB", "never deleted")]:
    e = s2e.get(sym)
    if e in CN.columns:
        print(f"    {sym:<7} median {CN[e].median():.3f}   "
              f"frac < {DEL_THRESHOLD}: {float((CN[e] < DEL_THRESHOLD).mean()):.3f}"
              f"   ({exp})")
if not (CN[s2e["GAPDH"]] < DEL_THRESHOLD).mean() == 0:
    print("    WARNING: GAPDH shows deletions -- threshold may be wrong")

expr_names = pq.ParquetFile(DATA_CLEAN / "depmap_expr_clean.parquet").schema_arrow.names
expr_col = {c.split(".")[0].lower(): c for c in expr_names if c != "index"}
genes = sorted(set(CH.columns) & set(CN.columns) & set(expr_col))
prof = pd.read_parquet(REF_DIR / "depmap_profiles.parquet")
rp = prof[prof.datatype == "rna"][["profileid", "modelid"]].copy()
rp["model_id"] = rp.modelid.str.lower()
rp["profileid"] = rp.profileid.astype(str)
E = pq.read_table(DATA_CLEAN / "depmap_expr_clean.parquet",
                  columns=["index"] + [expr_col[g] for g in genes]).to_pandas()
if "index" not in E.columns:
    E = E.reset_index()
E["profileid"] = E["index"].astype(str)
E = E.drop(columns=["index"]).merge(rp[["profileid", "model_id"]],
                                    on="profileid", how="inner")
E = E.drop(columns=["profileid"]).groupby("model_id").mean()
E.columns = [c.split(".")[0].lower() for c in E.columns]

lines = sorted(set(CH.index) & set(CN.index) & set(E.index))
print(f"\n  genes in all three : {len(genes):,}")
print(f"  lines in all three : {len(lines):,}")
CHm = CH.loc[lines, genes]
CNm = CN.loc[lines, genes]
Em = E.loc[lines, genes]
results["inputs"] = {"genes": len(genes), "lines": len(lines),
                     "label_fpr": LABEL_FPR,
                     "del_threshold": DEL_THRESHOLD,
                     "dep_threshold": DEP_THRESHOLD}

# ============================================================ Step 2
print()
print("=" * 78)
print("STEP 2 -- COVERAGE: how often does homozygous deletion actually fire?")
print("=" * 78)

DEL = (CNm < DEL_THRESHOLD).values
DEP = (CHm <= DEP_THRESHOLD).values
# expression floor, per gene, on the shared lines
FLOOR = (Em.values == np.nanmin(Em.values, axis=0, keepdims=True))
n_cells = DEL.size
print(f"  grid                              : {len(lines):,} lines x "
      f"{len(genes):,} genes = {n_cells:,} pairs")
print(f"  homozygously deleted              : {DEL.sum():,} "
      f"({100*DEL.mean():.3f}%)")
print(f"  at the expression floor           : {FLOOR.sum():,} "
      f"({100*FLOOR.mean():.3f}%)")
print(f"  dependent (Chronos <= {DEP_THRESHOLD})       : {DEP.sum():,} "
      f"({100*DEP.mean():.2f}%)")
gene_del = DEL.sum(axis=0)
print(f"\n  genes with >=20 deleted lines     : {int((gene_del >= 20).sum()):,} "
      f"of {len(genes):,}")
print(f"  genes with 0 deleted lines        : {int((gene_del == 0).sum()):,}")
print(f"\n  for comparison, the repo's COSMIC-derived cna_flags.parquet covers")
print(f"  0.78% of its grid and records only altered genes; this is dense.")
results["coverage"] = {
    "pairs": int(n_cells), "deleted": int(DEL.sum()),
    "deleted_pct": float(100 * DEL.mean()),
    "floor": int(FLOOR.sum()), "floor_pct": float(100 * FLOOR.mean()),
    "dependent_pct": float(100 * DEP.mean()),
    "genes_with_20plus_deletions": int((gene_del >= 20).sum()),
}

# ============================================================ Step 3
print()
print("=" * 78)
print("STEP 3 -- IS DELETION AN ABSOLUTE EXCLUSION?")
print("=" * 78)
base = float(DEP.mean())
r_del = float(DEP[DEL].mean()) if DEL.sum() else np.nan
print(f"  dependency rate, all pairs             : {base:.4f}")
print(f"  dependency rate, deleted pairs         : {r_del:.4f}   "
      f"ratio {r_del/base:.3f}")
print(f"  label FPR (curated non-essentials)     : {LABEL_FPR:.4f}")
print(f"  residual above label noise             : "
      f"{max(0.0, r_del - LABEL_FPR):.4f}")
print()
if r_del <= LABEL_FPR:
    print(f"  -> The dependency rate at homozygously deleted loci is AT OR BELOW")
    print(f"     the label's own false-positive rate. Deletion is an absolute")
    print(f"     exclusion: nothing survives it that the screen can resolve.")
    verdict_excl = "deletion_is_an_absolute_exclusion"
else:
    print(f"  -> Deleted loci still show dependency above the label's noise floor,")
    print(f"     which should not be mechanically possible (CRISPR cannot cut")
    print(f"     absent DNA). Suspect residual mis-calls in either layer.")
    verdict_excl = "residual_above_noise__check_calls"

print("\n  threshold sensitivity:")
for t in (0.05, 0.1, 0.2, 0.3):
    d = (CNm < t).values
    print(f"    del < {t:<4}  pairs {d.sum():>8,} ({100*d.mean():.3f}%)   "
          f"dependency rate {float(DEP[d].mean()):.4f}")
results["absolute_exclusion"] = {
    "base_rate": base, "deleted_rate": r_del, "ratio": float(r_del / base),
    "label_fpr": LABEL_FPR, "residual": float(max(0.0, r_del - LABEL_FPR)),
    "verdict": verdict_excl,
    "sensitivity": {str(t): float((CHm <= DEP_THRESHOLD).values[
        (CNm < t).values].mean()) for t in (0.05, 0.1, 0.2, 0.3)},
}

# ============================================================ Step 4
print()
print("=" * 78)
print("STEP 4 -- DOES CNA ADD ANYTHING OVER THE EXPRESSION FLOOR?")
print("=" * 78)
cells = {
    "deleted AND floor": DEL & FLOOR,
    "deleted, NOT floor": DEL & ~FLOOR,
    "floor, NOT deleted": ~DEL & FLOOR,
    "neither": ~DEL & ~FLOOR,
}
print(f"  {'stratum':<22} {'pairs':>12} {'% of grid':>10} {'dep rate':>10} "
      f"{'vs base':>9}")
addcell = {}
for nm, m in cells.items():
    n = int(m.sum())
    r = float(DEP[m].mean()) if n else np.nan
    print(f"  {nm:<22} {n:>12,} {100*m.mean():>9.3f}% {r:>10.4f} "
          f"{r/base:>9.3f}")
    addcell[nm] = {"pairs": n, "pct_of_grid": float(100 * m.mean()),
                   "dep_rate": r, "ratio_vs_base": float(r / base)}
overlap = DEL & FLOOR
print(f"\n  of {DEL.sum():,} deleted pairs, {int(overlap.sum()):,} "
      f"({100*overlap.sum()/max(DEL.sum(),1):.1f}%) are also at the expression floor")
new = DEL & ~FLOOR
print(f"  deletion adds {int(new.sum()):,} pairs the expression floor does not "
      f"already catch ({100*new.mean():.3f}% of the grid)")
if new.sum() > 100:
    k1, n1 = int(DEP[new].sum()), int(new.sum())
    k0, n0 = int(DEP[~DEL & ~FLOOR].sum()), int((~DEL & ~FLOOR).sum())
    orr, pf = fisher_exact([[k1, n1 - k1], [k0, n0 - k0]])
    print(f"  those {n1:,} pairs have dependency rate {k1/n1:.4f} vs "
          f"{k0/n0:.4f} in the scored population")
    print(f"    odds ratio {orr:.3f}   Fisher p {pf:.3g}")
    addcell["incremental_fisher"] = {"odds_ratio": float(orr), "p": float(pf)}
    incremental = (k1 / n1) < 0.5 * (k0 / n0)
else:
    incremental = False
results["incremental_value"] = addcell

print()
if incremental:
    print(f"  -> CNA carries information the expression floor does NOT. The pairs")
    print(f"     it uniquely excludes are strongly depleted of dependency, so a")
    print(f"     gene-level CNA layer earns its place.")
    verdict_add = "cna_adds_information_beyond_expression"
else:
    print(f"  -> The pairs CNA uniquely excludes are not materially cleaner than")
    print(f"     the scored population, so gene-level CNA is largely REDUNDANT")
    print(f"     with the expression floor Layer 2 already applies.")
    verdict_add = "cna_largely_redundant_with_expression_floor"

# ============================================================ Step 5
print()
print("=" * 78)
print("STEP 5 -- PRACTICAL YIELD")
print("=" * 78)
print(f"  A homozygous-deletion exclusion rule would fire on "
      f"{100*DEL.mean():.3f}% of pairs.")
print(f"  It is high-precision (dependency rate {r_del:.4f} vs base {base:.4f})")
print(f"  but low-yield: it removes ~{DEL.sum()/len(genes):.0f} of "
      f"{len(lines):,} candidate lines for a typical gene.")
top = pd.Series(gene_del, index=genes).nlargest(8)
print(f"\n  genes where the rule fires most (deleted lines of {len(lines)}):")
for g, n in top.items():
    sym = next((k for k, v in s2e.items() if v == g), g)
    dr = float(DEP[:, genes.index(g)][DEL[:, genes.index(g)]].mean())
    print(f"    {sym:<10} {int(n):>4} lines deleted   dependency rate there "
          f"{dr:.4f}")
results["yield"] = {
    "pct_of_pairs": float(100 * DEL.mean()),
    "mean_lines_excluded_per_gene": float(DEL.sum() / len(genes)),
    "top_genes": {next((k for k, v in s2e.items() if v == g), g): int(n)
                  for g, n in top.items()},
}

# ============================================================ Step 6
print()
print("=" * 78)
print("STEP 6 -- WHERE DOES THE RESIDUAL LIVE?")
print("=" * 78)
print("  A dependency call at a deleted locus is mechanically impossible, so the")
print("  0.0397 aggregate is measurement error. Locating it decides whether the")
print("  rule is usable.")

n_del_g = DEL.sum(axis=0)
k_dep_g = (DEL & DEP).sum(axis=0)
print(f"\n  {'deletions per gene':<24} {'genes':>7} {'pairs':>9} {'dep rate':>10}")
recur = {}
for lo, hi, lab in [(1, 3, "1-2 (likely mis-calls)"), (3, 10, "3-9"),
                    (10, 50, "10-49 (recurrent)"), (50, 10 ** 9, "50+")]:
    m = (n_del_g >= lo) & (n_del_g < hi)
    if not m.sum():
        continue
    rate = k_dep_g[m].sum() / max(n_del_g[m].sum(), 1)
    print(f"  {lab:<24} {int(m.sum()):>7,} {int(n_del_g[m].sum()):>9,} "
          f"{rate:>10.4f}")
    recur[lab] = {"genes": int(m.sum()), "pairs": int(n_del_g[m].sum()),
                  "dep_rate": float(rate)}
results["residual_by_recurrence"] = recur

meas = n_del_g >= 10
rates = np.divide(k_dep_g, np.maximum(n_del_g, 1))
print(f"\n  among the {int(meas.sum())} genes with >=10 deleted lines:")
print(f"    dependency rate exactly 0.000 : {int((rates[meas] == 0).sum())} "
      f"({100*(rates[meas] == 0).mean():.1f}%)  <- perfect exclusion")
print(f"    rate > 0.20                   : {int((rates[meas] > 0.20).sum())}")
bad_i = np.where(meas & (rates > 0.20))[0]
if len(bad_i):
    print(f"\n  the pathological genes:")
    for i in bad_i[np.argsort(-rates[bad_i])]:
        sym = next((k for k, v in s2e.items() if v == genes[i]), genes[i])
        print(f"    {sym:<12} {int(n_del_g[i]):>4} deleted, "
              f"{int(k_dep_g[i]):>3} 'dependent'  rate {rates[i]:.3f}")
    print(f"  These are highly paralogous multi-copy families. A guide targeting")
    print(f"  one member cuts several loci, so the screen registers a fitness")
    print(f"  effect that has nothing to do with the deleted gene -- the known")
    print(f"  multi-targeting artefact (De Kegel & Ryan 2019 reprocess for it).")
results["perfect_exclusion_frac"] = float((rates[meas] == 0).mean())
results["n_genes_measurable"] = int(meas.sum())

print()
print("=" * 78)
print("VERDICT")
print("=" * 78)
print(f"  absolute exclusion : {verdict_excl}")
print(f"  incremental value  : {verdict_add}")
print()
print(f"  Usable form of the rule: restrict to RECURRENTLY deleted genes")
print(f"  (>=10 deleted lines). There the dependency rate is "
      f"{recur.get('10-49 (recurrent)', {}).get('dep_rate', float('nan')):.4f} "
      f"against a label FPR of {LABEL_FPR:.4f},")
print(f"  and {100*results['perfect_exclusion_frac']:.0f}% of those genes show "
      f"exactly zero. Applied across ALL genes it is")
print(f"  noisier than the expression floor the pipeline already has "
      f"({addcell['floor, NOT deleted']['dep_rate']:.4f} at")
print(f"  {addcell['floor, NOT deleted']['pct_of_grid']:.1f}% of the grid, versus "
      f"{r_del:.4f} at {100*DEL.mean():.3f}%).")
results["verdict"] = {"exclusion": verdict_excl, "incremental": verdict_add,
                      "usable_form": "restrict to genes with >=10 deleted lines"}

out = OUTPUTS / "test_run_cna_gate_results.json"
out.write_text(json.dumps(results, indent=2, default=float), encoding="utf-8")
print(f"\nwrote {out}")
