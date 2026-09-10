"""
test_run_abundance_gate.py
--------------------------
Read-only diagnostic. Writes no production file.

QUESTION
--------
test_run_stouffer_regime_split.py measured Layer 2 against GDSC drug sensitivity
with AUROC and got 0.52-0.56 for every scoring rule tried. That was read as a
weak result. It may instead be the wrong readout.

core_score is an abundance / expressibility PRIOR -- a necessary-but-not-
sufficient condition -- not a dependency estimate. AUROC integrates over all
thresholds and rewards enrichment at the TOP of the ranking. A gate's job is at
the BOTTOM: the DepMap principle the pipeline leans on is that a gene which is
not expressed in a cell line will behave as non-essential there regardless of
its function. The claim being made is therefore NOT

    "high-abundance lines are enriched for drug sensitivity"       (AUROC's question)

but

    "low-abundance lines are DEPLETED of drug sensitivity"         (a gate's question)

A filter measured with a ranking metric will score near chance even when it is
working perfectly, because the metric spends most of its mass on a region the
filter makes no claim about.

METHOD
------
Same construction as test_run_stouffer_regime_split.py -- same GDSC labels, same
van der Waerden layers, same Stouffer Z -- with three gate-shaped readouts:

  1. DECILE PROFILE. Per gene, sensitivity rate in each decile of the score,
     expressed as a ratio to that gene's base rate. A working gate shows
     ratio << 1 in decile 1 and ratio ~ 1 (flat) across the upper deciles.
     A working RANKER shows a monotone gradient across all ten.
     These are different signatures and the decile profile separates them.

  2. NEGATIVE PREDICTIVE VALUE at candidate abstention thresholds (bottom 5%,
     10%, 20%). NPV = P(not sensitive | score below threshold). Reported against
     the base NPV (= 1 - prevalence), because NPV is high whenever prevalence is
     low and the uplift over base is the only honest quantity.

  3. PARTIAL AUROC over the low-specificity region, FPR in [0.8, 1.0], with
     McClish (1989) standardisation onto [0.5, 1]. That is the region a gate
     actually operates in: a rule that excludes the bottom decile is sitting at
     roughly FPR 0.9. The symmetric high-specificity region FPR in [0.0, 0.2] is
     reported alongside, since the contrast between them is the whole point.

  Plus a FLOOR test, which is the DepMap principle in its most literal form:
  positive rate among lines whose RNA sits exactly at the gene's expression
  floor, versus that gene's base rate. No thresholding, no deciles -- these are
  the lines where the gene is simply not expressed.

TWO LABEL SOURCES  (--label)
----------------------------
  gdsc     GDSC drug sensitivity. 141 genes. A line is positive when it is
           `sensitive` to a drug whose ANNOTATED target is this gene.
  chronos  CRISPR knockout fitness (validation/prepared/chronos_long.parquet).
           ~17.6k genes. A line is positive when knocking the gene out impairs
           fitness in that line.

The GDSC run was done first and returned floor ratio 0.898 -- barely any
depletion. That is not a fair test of the DepMap principle, and the reason is
causal distance. The principle is a claim about KNOCKOUT: remove a gene the cell
does not express and nothing happens. GDSC interposes a whole chain between the
score and the label -- drug reaches the cell, drug hits the annotated target
rather than one of its many off-targets, inhibiting that target kills the cell.
Clinical kinase inhibitors are extensively polypharmacological (Klaeger et al.
2017, Science -- kinobead profiling of 243 clinical drugs), so a line can be
genuinely drug-sensitive through a protein other than the annotated one. That
puts a ceiling on any depletion GDSC can show, regardless of whether the
underlying principle holds.

Chronos removes every link in that chain: the perturbation IS the gene, and the
readout IS fitness. It is also ~125x more genes and a dense continuous label.
It is therefore the test that decides the question; GDSC decides only whether
drug annotation is informative enough to validate against.

SIGN CONVENTION (chronos)
-------------------------
Raw `essentiality` in chronos_long.parquet is more NEGATIVE = more essential,
verified here rather than assumed: EIF4A3 -16.5, RPS6 -15.7, RPL13A -12.5,
PSMB2 -12.4 (core-essential) against ALB +4.2, KLK4 +3.1, CD3E +3.1
(tissue-restricted). Step 1 re-runs that check every time and refuses to
continue if it fails. Binary positive = essentiality <= DEP_THRESHOLD, which
yields 12.2% prevalence, close to GDSC's 13.5% so the two runs are comparable.
A threshold-free continuous readout (within-gene essentiality percentile per
decile) is reported alongside, since the threshold is a convention.

Significance: the null for the decile-1 count is exactly hypergeometric --
"draw n1 of this gene's N lines without replacement, how often do we see at most
the observed number of positives". That is the same null a within-gene label
permutation would approximate, holding prevalence and line count fixed, but it
is exact and closed-form, so it scales to thousands of genes.

READING THE RESULT
------------------
  bottom decile depleted, top deciles flat -> the gate works; AUROC 0.56 was a
      measurement error and Layer 2 should be evaluated and DOCUMENTED as a
      filter, not a ranker
  bottom decile depleted AND a top gradient -> it is both, and AUROC understates
      it but is not meaningless
  bottom decile not depleted -> the real negative result, and a far more serious
      one than any combination-rule difference measured so far

SIGN CONVENTION
---------------
Higher score = higher abundance = stronger candidate. GDSC `sensitive` = True is
the positive class. A ratio below 1.0 in decile 1 is DEPLETION and is the
outcome the gate hypothesis predicts.

Outputs (suffixed by label):
    src/pipeline/outputs/test_run_abundance_gate_<label>_results.json
    src/pipeline/outputs/test_run_abundance_gate_<label>_per_gene.parquet

Run:
    python src/pipeline/test_run_abundance_gate.py --label gdsc
    python src/pipeline/test_run_abundance_gate.py --label chronos --n-genes 3000
"""
import argparse
import json
import warnings
from pathlib import Path

import duckdb
import numpy as np
import pandas as pd
import pyarrow.parquet as pq
from scipy.stats import hypergeom, norm, rankdata, wilcoxon

warnings.filterwarnings("ignore", category=RuntimeWarning)

OUTPUTS = Path("src/pipeline/outputs")
REF_DIR = Path("reference")
DATA_CLEAN = Path("data/parquet/data_clean")
CLEANED_TRACK = Path("cleaned_track_data")

RHO_GLOBAL = 0.46
W_E, W_P = 0.987, 0.381
MIN_POS = 10
MIN_NEG = 10
MIN_LINES = 100          # deciles need enough lines to be meaningful
DEP_THRESHOLD = -0.5     # chronos: essentiality at or below this = dependent
SEED = 42

ap = argparse.ArgumentParser()
ap.add_argument("--label", choices=["gdsc", "project_score", "chronos"],
                default="gdsc",
                help="gdsc = drug sensitivity; project_score = the MISNAMED "
                     "chronos_long.parquet (Project Score scaled BF, negated); "
                     "chronos = the real Chronos gene effect HDF5")
ap.add_argument("--chronos-path",
                default="data/DepMap_Chronos/GeneFitnessEffect_Chronos_Achilles.hdf5")
ap.add_argument("--n-genes", type=int, default=3000,
                help="essentiality labels only: gene sample size "
                     "(GDSC targets always included)")
args = ap.parse_args()
LABEL = args.label

results = {"label_source": LABEL}


# ============================================================ helpers
def vdw(s: pd.Series) -> pd.Series:
    n = int(s.notna().sum())
    return pd.Series(norm.ppf(s.rank(method="min") / (n + 1)), index=s.index)


def stouffer(zE, zP, rho=RHO_GLOBAL, wE=W_E, wP=W_P):
    out = np.full(len(zE), np.nan)
    e_ok, p_ok = zE.notna().values, zP.notna().values
    both = e_ok & p_ok
    ze, zp = zE.values, zP.values
    out[e_ok & ~p_ok] = ze[e_ok & ~p_ok]
    out[p_ok & ~e_ok] = zp[p_ok & ~e_ok]
    out[both] = (wE * ze[both] + wP * zp[both]) / np.sqrt(
        wE ** 2 + wP ** 2 + 2 * wE * wP * rho)
    return pd.Series(out, index=zE.index)


def roc_curve(score, pos):
    """FPR, TPR arrays. Ties handled by grouping equal scores into one step."""
    order = np.argsort(-score, kind="mergesort")
    s, p = score[order], pos[order]
    P, N = p.sum(), (~p).sum()
    # step at each distinct score value
    distinct = np.r_[np.diff(s) != 0, True]
    tp = np.cumsum(p)[distinct]
    fp = np.cumsum(~p)[distinct]
    return np.r_[0, fp / N], np.r_[0, tp / P]


def partial_auc(score, pos, lo, hi):
    """McClish-standardised partial AUROC over FPR in [lo, hi]. 0.5 = chance."""
    ok = ~np.isnan(score)
    score, pos = score[ok], pos[ok]
    if pos.sum() < 2 or (~pos).sum() < 2:
        return np.nan
    fpr, tpr = roc_curve(score, pos)
    grid = np.linspace(lo, hi, 501)
    t = np.interp(grid, fpr, tpr)
    area = np.trapezoid(t, grid)
    a_min = (hi ** 2 - lo ** 2) / 2.0        # chance line over the same interval
    a_max = hi - lo
    return 0.5 * (1.0 + (area - a_min) / (a_max - a_min))


def auroc(score, pos):
    ok = ~np.isnan(score)
    score, pos = score[ok], pos[ok]
    n1, n0 = int(pos.sum()), int((~pos).sum())
    if n1 < 1 or n0 < 1:
        return np.nan
    r = rankdata(score)
    return (r[pos].sum() - n1 * (n1 + 1) / 2.0) / (n1 * n0)


# ============================================================ Step 1
print("=" * 78)
print(f"STEP 1 -- labels ({LABEL}), gene set, layers")
print("=" * 78)

gdsc = pd.read_parquet("validation/prepared/gdsc_scored_ready.parquet",
                       columns=["model_id", "target_ensg", "sensitive"])
gdsc["model_id"] = gdsc.model_id.str.lower()
gdsc["target_ensg"] = gdsc.target_ensg.astype("string").str.split(".").str[0].str.lower()
sens = gdsc[gdsc.sensitive][["target_ensg", "model_id"]].drop_duplicates()
gdsc_genes = sorted(sens.target_ensg.unique())

chr_pct = {}          # gene -> Series(model_id -> within-gene essentiality percentile)
LABEL_FPR = None
if LABEL == "gdsc":
    pos_by_gene = sens.groupby("target_ensg")["model_id"].apply(set).to_dict()
    target_genes = gdsc_genes
    print(f"  GDSC target genes with >=1 sensitive line : {len(target_genes)}")
else:
    gl_chk = pd.read_parquet(REF_DIR / "gene_lookup.parquet",
                             columns=["ensg_id", "hgnc_symbol"])
    gl_chk["ensg_id"] = gl_chk.ensg_id.astype(str).str.split(".").str[0].str.lower()
    s2e = dict(zip(gl_chk.hgnc_symbol.astype(str).str.upper(), gl_chk.ensg_id))

    if LABEL == "project_score":
        ch = pd.read_parquet("validation/prepared/chronos_long.parquet")
        ch["ensg_id"] = ch.ensg_id.astype(str).str.split(".").str[0].str.lower()
        core = [s2e[s] for s in ["RPL13A", "RPS6", "EIF4A3", "PSMB2"] if s in s2e]
        rest = [s2e[s] for s in ["ALB", "KLK4", "CD3E", "MYOD1"] if s in s2e]
        med = ch.groupby("ensg_id").essentiality.median()
        m_core, m_rest = med.reindex(core).median(), med.reindex(rest).median()
        print(f"  sign check: core-essential median {m_core:+.2f} vs "
              f"tissue-restricted median {m_rest:+.2f}")
        if not (m_core < m_rest and m_core < 0):
            raise SystemExit("sign convention check FAILED. Aborting.")
        ib = pd.read_parquet("validation/prepared/id_bridge.parquet")
        ch = ch.merge(ib, on="sanger_model_id", how="inner")
        ch["model_id"] = ch.model_id.str.lower()
        print("  WARNING: this is Project Score scaled BF (negated), NOT Chronos.")
        print(f"           {DEP_THRESHOLD} is a quantile here, not a gene effect.")
        LABEL_FPR = None
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
        ch = W.stack().rename("essentiality").reset_index()
        ch.columns = ["model_id", "ensg_id", "essentiality"]
        print(f"  real Chronos ({cp.name}): {W.shape[0]} lines x "
              f"{W.shape[1]:,} genes")

        # Scale check against the release's own curated reference sets, and the
        # LABEL FPR the floor test has to be judged against. Its absence is what
        # let a negated Project Score matrix pass as Chronos project-wide.
        ess = set(pd.read_csv(cp.parent / "ReferenceEssentials.csv").iloc[:, 0]
                  .astype(str).str.split(" (", regex=False).str[0].str.upper())
        non = set(pd.read_csv(cp.parent / "ReferenceNonEssentials.csv").iloc[:, 0]
                  .astype(str).str.split(" (", regex=False).str[0].str.upper())
        e_ess = {s2e[g] for g in ess if g in s2e} & set(W.columns)
        e_non = {s2e[g] for g in non if g in s2e} & set(W.columns)
        m_ess = float(np.nanmedian(W[list(e_ess)].values))
        m_non = float(np.nanmedian(W[list(e_non)].values))
        LABEL_FPR = float(np.nanmean(W[list(e_non)].values <= DEP_THRESHOLD))
        print(f"  scale check: median {np.nanmedian(W.values):+.4f} (expect ~0); "
              f"essentials {m_ess:+.4f} (~-1); non-essentials {m_non:+.4f} (~0)")
        if not (m_ess < -0.5 < m_non + 0.5):
            raise SystemExit("Chronos scale check FAILED. Aborting.")
        print(f"  LABEL FPR at {DEP_THRESHOLD}: {LABEL_FPR:.4f}  "
              f"({len(e_non)} curated non-essential genes)")
        results["label_fpr"] = LABEL_FPR
        results["reference_medians"] = {"essential": m_ess, "nonessential": m_non}

    # Sample genes for tractability; always keep the GDSC targets for comparability.
    rng0 = np.random.default_rng(SEED)
    all_ch = sorted(ch.ensg_id.unique())
    keep = set(gdsc_genes) & set(all_ch)
    pool = [g for g in all_ch if g not in keep]
    if args.n_genes < len(pool):
        keep |= set(rng0.choice(pool, size=args.n_genes - len(keep), replace=False))
    else:
        keep |= set(pool)
    ch = ch[ch.ensg_id.isin(keep)]
    print(f"  sampled to {ch.ensg_id.nunique():,} genes "
          f"({len(set(gdsc_genes) & keep)} GDSC targets forced in)")

    ch["dep"] = ch.essentiality <= DEP_THRESHOLD
    print(f"  binary positive = essentiality <= {DEP_THRESHOLD}  -> "
          f"prevalence {ch.dep.mean():.3f}")
    pos_by_gene = ch[ch.dep].groupby("ensg_id")["model_id"].apply(set).to_dict()
    # threshold-free companion: within-gene percentile, high = more essential
    ch["pct"] = 1.0 - ch.groupby("ensg_id")["essentiality"].rank(pct=True,
                                                                method="average")
    chr_pct = {g: s.set_index("model_id")["pct"]
               for g, s in ch[["ensg_id", "model_id", "pct"]].groupby("ensg_id")}
    target_genes = sorted(pos_by_gene)
    print(f"  genes with >=1 dependent line             : {len(target_genes)}")

regime = pd.read_parquet(OUTPUTS / "gene_regime.parquet", columns=["ensg_id", "class"])
regime["ensg_id"] = regime.ensg_id.astype(str).str.lower()
class_of = dict(zip(regime.ensg_id, regime["class"]))

gl = pd.read_parquet(REF_DIR / "gene_lookup.parquet",
                     columns=["ensg_id", "hgnc_symbol", "biotype",
                              "hgnc_status", "uniprot_ids"])
uni = gl[(gl.biotype == "protein_coding") & (gl.hgnc_status == "Approved")].copy()
uni["ensg_id"] = uni.ensg_id.astype("string").str.split(".").str[0].str.lower()
valid_ensg = set(uni.ensg_id)
symbol_of = dict(zip(uni.ensg_id, uni.hgnc_symbol))

u2e = {}
for e, u in zip(uni.ensg_id, uni.uniprot_ids):
    if isinstance(u, str):
        for acc in u.split("|"):
            acc = acc.strip().lower()
            if acc:
                u2e.setdefault(acc, e)

expr_names = pq.ParquetFile(DATA_CLEAN / "depmap_expr_clean.parquet").schema_arrow.names
expr_col = {c.split(".")[0].lower(): c for c in expr_names if c != "index"}

prot_names = pq.ParquetFile(CLEANED_TRACK / "proteomics.parquet").schema_arrow.names
key_col = "model_id" if "model_id" in prot_names else "depmap_id"
ccle_by_ensg = {}
for c in prot_names:
    if c == key_col:
        continue
    a = c.strip().lower()
    e = u2e.get(a) or u2e.get(a.split("-")[0])
    if e in valid_ensg:
        ccle_by_ensg.setdefault(e, []).append(c)

con = duckdb.connect(str(OUTPUTS / "celllineselector.db"), read_only=True)
pc_names = [r[0] for r in con.execute("DESCRIBE procan_proteomics").fetchall()]
PC_META = {"gdsc_model_name", "sanger_model_id", "model_id",
           "matched_via", "n_model_id", "is_ambiguous"}
procan_by_ensg = {}
for c in pc_names:
    if c in PC_META:
        continue
    e = u2e.get(c.strip().lower())
    if e in valid_ensg:
        procan_by_ensg.setdefault(e, []).append(c)

genes = [g for g in target_genes if g in expr_col]
print(f"  GDSC target genes                 : {len(target_genes)}")
print(f"  ... with an RNA column            : {len(genes)}")

prof = pd.read_parquet(REF_DIR / "depmap_profiles.parquet")
rna_prof = prof[prof.datatype == "rna"][["profileid", "modelid"]].copy()
rna_prof["model_id"] = rna_prof.modelid.str.lower()
rna_prof["profileid"] = rna_prof.profileid.astype(str)

E = pq.read_table(DATA_CLEAN / "depmap_expr_clean.parquet",
                  columns=["index"] + [expr_col[g] for g in genes]).to_pandas()
if "index" not in E.columns:
    E = E.reset_index()
E["profileid"] = E["index"].astype(str)
E = E.drop(columns=["index"]).merge(rna_prof[["profileid", "model_id"]],
                                    on="profileid", how="inner")
E = E.drop(columns=["profileid"]).groupby("model_id").mean()
E.columns = [c.split(".")[0].lower() for c in E.columns]

need_ccle = sorted({a for g in genes for a in ccle_by_ensg.get(g, [])})
C = pq.read_table(CLEANED_TRACK / "proteomics.parquet",
                  columns=[key_col] + need_ccle).to_pandas()
C.index = C[key_col].str.lower().values
C = C.drop(columns=[key_col]).astype(float).groupby(level=0).mean()

need_pc = sorted({a for g in genes for a in procan_by_ensg.get(g, [])})
q = ", ".join(f'"{a}"' for a in need_pc)
P = con.execute(f"SELECT model_id, {q} FROM procan_proteomics "
                "WHERE model_id IS NOT NULL").df()
P.index = P.model_id.str.lower().values
P = P.drop(columns=["model_id"]).astype(float).groupby(level=0).mean()
con.close()
print(f"  RNA {E.shape[0]} lines | CCLE {C.shape[0]} | ProCan {P.shape[0]}")


def layers_for(g):
    e = E[g].dropna()
    zE = vdw(e) if len(e) else pd.Series(dtype=float)
    parts = {}
    accs = [a for a in ccle_by_ensg.get(g, []) if a in C.columns]
    if accs:
        s = C[accs].mean(axis=1).dropna()
        if len(s):
            parts["ccle"] = vdw(s)
    accs = [a for a in procan_by_ensg.get(g, []) if a in P.columns]
    if accs:
        s = P[accs].mean(axis=1).dropna()
        if len(s):
            parts["procan"] = vdw(s)
    if not parts:
        return e, zE, pd.Series(dtype=float)
    idx = sorted(set().union(*[set(v.index) for v in parts.values()]))
    zP = pd.DataFrame({k: v.reindex(idx) for k, v in parts.items()}).mean(axis=1)
    return e, zE, zP


# ============================================================ Step 2
print()
print("=" * 78)
print("STEP 2 -- decile profile, NPV, partial AUROC, floor test")
print("=" * 78)

rows, dec_Z, dec_R, dec_cont = [], [], [], []
for g in genes:
    e_raw, zE, zP = layers_for(g)
    if not len(zE):
        continue
    idx = sorted(set(zE.index) | set(zP.index)) if len(zP) else list(zE.index)
    d = pd.DataFrame({"zE": zE.reindex(idx),
                      "zP": zP.reindex(idx) if len(zP) else np.nan})
    d = d[d.notna().any(axis=1)]
    d["Z"] = stouffer(d.zE, d.zP) if len(zP) else d.zE
    d["pos"] = d.index.isin(pos_by_gene[g])
    if len(d) < MIN_LINES or d.pos.sum() < MIN_POS or (~d.pos).sum() < MIN_NEG:
        continue

    base = float(d.pos.mean())
    rec = {"ensg_id": g, "symbol": symbol_of.get(g, "?"), "class": class_of.get(g),
           "n_lines": int(len(d)), "n_pos": int(d.pos.sum()), "base_rate": base}

    for tag, col in [("Z", "Z"), ("R", "zE")]:
        sc = d[col].values
        ok = ~np.isnan(sc)
        s, p = sc[ok], d.pos.values[ok]
        # deciles by rank so tie blocks distribute deterministically
        r = rankdata(s, method="ordinal")
        q = np.minimum((r - 1) * 10 // len(s), 9)
        rates = np.array([p[q == i].mean() if (q == i).sum() else np.nan
                          for i in range(10)])
        ratios = rates / base
        (dec_Z if tag == "Z" else dec_R).append(ratios)
        rec[f"d1_rate_{tag}"] = float(rates[0])
        rec[f"d1_ratio_{tag}"] = float(ratios[0])
        rec[f"d10_ratio_{tag}"] = float(ratios[9])
        rec[f"upper_flat_{tag}"] = float(np.nanstd(ratios[4:]))
        # NPV at candidate abstention cuts
        for pctl in (5, 10, 20):
            k = max(1, int(round(len(s) * pctl / 100)))
            below = np.argsort(s, kind="mergesort")[:k]
            rec[f"npv{pctl}_{tag}"] = float(1.0 - p[below].mean())
        rec[f"npv_base_{tag}"] = float(1.0 - base)
        # partial AUROC, gate region and ranker region
        rec[f"pauc_low_{tag}"] = partial_auc(s, p, 0.8, 1.0)
        rec[f"pauc_high_{tag}"] = partial_auc(s, p, 0.0, 0.2)
        rec[f"auc_{tag}"] = auroc(s, p)
        # permutation null for the decile-1 ratio
        # One-sided (depletion is the directional hypothesis). Drawing n1 of the
        # gene's N lines without replacement holds prevalence and size fixed, so
        # the null count is exactly hypergeometric -- no permutation needed.
        n1 = int((q == 0).sum())
        obs = int(p[q == 0].sum())
        rec[f"d1_perm_p_{tag}"] = float(
            hypergeom.cdf(obs, len(p), int(p.sum()), n1))

    # Threshold-free companion (chronos only): median within-gene essentiality
    # percentile per decile of Z. Says the same thing without a cut-off.
    if LABEL != "gdsc" and g in chr_pct:
        pc = chr_pct[g].reindex(d.index)
        sc = d["Z"].values
        ok = ~np.isnan(sc) & pc.notna().values
        if ok.sum() >= MIN_LINES:
            r = rankdata(sc[ok], method="ordinal")
            qq = np.minimum((r - 1) * 10 // ok.sum(), 9)
            v = pc.values[ok]
            dec_cont.append([np.nanmedian(v[qq == i]) if (qq == i).sum() else np.nan
                             for i in range(10)])
            rec["cont_d1"] = float(np.nanmedian(v[qq == 0]))
            rec["cont_d10"] = float(np.nanmedian(v[qq == 9]))

    # FLOOR test -- the DepMap principle in its most literal form
    fl = e_raw[e_raw == e_raw.min()].index
    fl = [m for m in fl if m in d.index]
    rec["n_floor"] = len(fl)
    rec["floor_share"] = len(fl) / len(d)
    if len(fl) >= 20:
        rec["floor_rate"] = float(d.loc[fl, "pos"].mean())
        rec["floor_ratio"] = rec["floor_rate"] / base
    rows.append(rec)

G = pd.DataFrame(rows)
G.to_parquet(OUTPUTS / f"test_run_abundance_gate_{LABEL}_per_gene.parquet")
DZ = np.vstack(dec_Z)
DR = np.vstack(dec_R)
print(f"  genes scored: {len(G)}   median lines/gene {G.n_lines.median():.0f}   "
      f"median prevalence {G.base_rate.median():.3f}")

# ============================================================ Step 3
print()
print("=" * 78)
print("STEP 3 -- DECILE PROFILE.  sensitivity rate / base rate, median across genes")
print("=" * 78)
print("  decile 1 = lowest abundance.  1.00 = exactly the base rate.")
print()
print(f"  {'decile':>7} {'Stouffer Z':>12} {'RNA alone':>12}   {'interpretation':<30}")
for i in range(10):
    z, r = np.nanmedian(DZ[:, i]), np.nanmedian(DR[:, i])
    bar = "#" * max(0, int(round(z * 20)))
    print(f"  {i+1:>7} {z:>12.3f} {r:>12.3f}   {bar}")
med_prof = np.array([np.nanmedian(DZ[:, i]) for i in range(10)])
d1z, d10z = med_prof[0], med_prof[9]
# Flatness must be measured on the MEDIAN PROFILE, not per gene. The per-gene
# spread is dominated by sampling noise at low prevalence and says nothing about
# whether the profile slopes -- an earlier version used it and mislabelled a flat
# chronos profile as "sloped".
upper = float(np.nanstd(med_prof[1:]))              # deciles 2-10 of the profile
slope = float(np.polyfit(np.arange(2, 11), med_prof[1:], 1)[0])
per_gene_noise = float(np.nanmedian(np.nanstd(DZ[:, 4:], axis=1)))
print()
print(f"  decile 1  ratio (median across genes) : {d1z:.3f}   "
      f"({'depleted' if d1z < 1 else 'enriched'} by {abs(1-d1z)*100:.1f}% vs base)")
print(f"  decile 10 ratio                       : {d10z:.3f}")
print(f"  deciles 2-10 of the median profile    : sd {upper:.3f}, "
      f"slope {slope:+.4f}/decile")
print(f"     (flat = gate signature; sloping = ranker signature)")
print(f"  per-gene noise in deciles 5-10        : sd {per_gene_noise:.3f}  "
      f"(sampling noise, NOT a slope)")
print(f"  genes with decile-1 ratio < 1         : "
      f"{100*(G.d1_ratio_Z < 1).mean():.1f}%")
print(f"  genes with decile-1 permutation p<.05 : "
      f"{100*(G.d1_perm_p_Z < 0.05).mean():.1f}%")
if dec_cont:
    DC = np.vstack(dec_cont)
    print()
    print("  THRESHOLD-FREE COMPANION -- median within-gene essentiality percentile")
    print("  (0 = least essential line for this gene, 1 = most essential)")
    print(f"  {'decile':>7} {'ess. percentile':>17}")
    for i in range(10):
        print(f"  {i+1:>7} {np.nanmedian(DC[:, i]):>17.3f}")
    results["continuous_decile_profile"] = [float(np.nanmedian(DC[:, i]))
                                            for i in range(10)]
    print()

wsr = wilcoxon(G.d1_ratio_Z.dropna() - 1.0)
print(f"  Wilcoxon signed-rank, decile-1 ratio vs 1.0 : p = {wsr.pvalue:.4g}")

results["decile_profile"] = {
    "median_ratio_by_decile_Z": [float(np.nanmedian(DZ[:, i])) for i in range(10)],
    "median_ratio_by_decile_RNA": [float(np.nanmedian(DR[:, i])) for i in range(10)],
    "decile1_ratio_Z": float(d1z), "decile10_ratio_Z": float(d10z),
    "upper_profile_sd": float(upper),
    "upper_profile_slope_per_decile": float(slope),
    "per_gene_noise_sd": per_gene_noise,
    "frac_genes_depleted": float((G.d1_ratio_Z < 1).mean()),
    "frac_genes_perm_sig": float((G.d1_perm_p_Z < 0.05).mean()),
    "wilcoxon_p": float(wsr.pvalue),
}

# ============================================================ Step 4
print()
print("=" * 78)
print("STEP 4 -- NEGATIVE PREDICTIVE VALUE at candidate abstention thresholds")
print("=" * 78)
print(f"  base NPV (= 1 - prevalence), median : {G.npv_base_Z.median():.4f}")
for pctl in (5, 10, 20):
    npv = G[f"npv{pctl}_Z"]
    up = (npv - G.npv_base_Z)
    print(f"  cut at bottom {pctl:>2}% of Z : NPV {npv.median():.4f}   "
          f"uplift over base {up.median():+.4f}   "
          f"positive in {100*(up > 0).mean():.1f}% of genes")
results["npv"] = {
    "base_median": float(G.npv_base_Z.median()),
    **{f"npv{p}_median": float(G[f"npv{p}_Z"].median()) for p in (5, 10, 20)},
    **{f"npv{p}_uplift_median": float((G[f"npv{p}_Z"] - G.npv_base_Z).median())
       for p in (5, 10, 20)},
}

# ============================================================ Step 5
print()
print("=" * 78)
print("STEP 5 -- PARTIAL AUROC (McClish standardised, 0.5 = chance)")
print("=" * 78)
print(f"  {'region':<38} {'Stouffer Z':>12} {'RNA alone':>12}")
print(f"  {'FPR 0.8-1.0  (the gate region)':<38} "
      f"{G.pauc_low_Z.median():>12.4f} {G.pauc_low_R.median():>12.4f}")
print(f"  {'FPR 0.0-0.2  (the ranker region)':<38} "
      f"{G.pauc_high_Z.median():>12.4f} {G.pauc_high_R.median():>12.4f}")
print(f"  {'full AUROC (for reference)':<38} "
      f"{G.auc_Z.median():>12.4f} {G.auc_R.median():>12.4f}")
pl = wilcoxon(G.pauc_low_Z.dropna() - 0.5)
ph = wilcoxon(G.pauc_high_Z.dropna() - 0.5)
print(f"\n  gate region vs chance   : p = {pl.pvalue:.4g}   "
      f"above 0.5 in {100*(G.pauc_low_Z > 0.5).mean():.1f}% of genes")
print(f"  ranker region vs chance : p = {ph.pvalue:.4g}   "
      f"above 0.5 in {100*(G.pauc_high_Z > 0.5).mean():.1f}% of genes")
results["partial_auc"] = {
    "gate_region_Z": float(G.pauc_low_Z.median()),
    "gate_region_RNA": float(G.pauc_low_R.median()),
    "ranker_region_Z": float(G.pauc_high_Z.median()),
    "ranker_region_RNA": float(G.pauc_high_R.median()),
    "full_auc_Z": float(G.auc_Z.median()),
    "gate_p": float(pl.pvalue), "ranker_p": float(ph.pvalue),
    "gate_frac_above_chance": float((G.pauc_low_Z > 0.5).mean()),
    "ranker_frac_above_chance": float((G.pauc_high_Z > 0.5).mean()),
}

# ============================================================ Step 6
print()
print("=" * 78)
print("STEP 6 -- FLOOR TEST.  Lines where the gene is simply not expressed.")
print("=" * 78)
F = G.dropna(subset=["floor_ratio"])
print(f"  genes with a floor block of >=20 lines : {len(F)} of {len(G)}")
if len(F):
    print(f"  median floor block size                : {F.n_floor.median():.0f} lines "
          f"({F.floor_share.median():.1%} of the gene's lines)")
    print(f"  median positive rate at the floor      : {F.floor_rate.median():.4f}")
    print(f"  median base rate for the same genes    : {F.base_rate.median():.4f}")
    print(f"  median ratio (floor / base)            : {F.floor_ratio.median():.3f}")
    print(f"  genes with ratio < 1                   : "
          f"{100*(F.floor_ratio < 1).mean():.1f}%")
    if len(F) >= 6:
        wf = wilcoxon(F.floor_ratio - 1.0)
        print(f"  Wilcoxon vs 1.0                        : p = {wf.pvalue:.4g}")
        results["floor_test"] = {
            "n_genes": int(len(F)), "median_ratio": float(F.floor_ratio.median()),
            "median_floor_rate": float(F.floor_rate.median()),
            "median_base_rate": float(F.base_rate.median()),
            "frac_below_1": float((F.floor_ratio < 1).mean()),
            "wilcoxon_p": float(wf.pvalue)}

    # --- stratified by prevalence. The gate decays with prevalence (measured in
    #     test_run_improvement_scan.py on real Chronos), so a pooled floor ratio
    #     mixes bands where the gate works with bands where it provably cannot.
    def boot_ci(v, n_boot=2000, stat=np.median):
        v = np.asarray(pd.Series(v).dropna())
        if len(v) < 5:
            return (np.nan, np.nan)
        rg = np.random.default_rng(SEED)
        b = [stat(rg.choice(v, len(v), replace=True)) for _ in range(n_boot)]
        return (float(np.percentile(b, 2.5)), float(np.percentile(b, 97.5)))

    F = F.copy()
    F["band"] = pd.cut(F.base_rate, [0, .02, .05, .10, .25, 1.01],
                       labels=["<2%", "2-5%", "5-10%", "10-25%", ">25%"])
    print()
    print(f"  STRATIFIED BY PREVALENCE")
    print(f"  {'band':>8} {'genes':>6} {'floor rate':>11} {'base':>8} "
          f"{'ratio':>7} {'95% CI':>18}")
    band_rows = []
    for b in F.band.cat.categories:
        sub = F[F.band == b]
        if len(sub) < 5:
            continue
        lo, hi = boot_ci(sub.floor_ratio)
        print(f"  {str(b):>8} {len(sub):>6} {sub.floor_rate.median():>11.4f} "
              f"{sub.base_rate.median():>8.4f} {sub.floor_ratio.median():>7.3f} "
              f"[{lo:>7.3f},{hi:>7.3f}]")
        band_rows.append({"band": str(b), "genes": int(len(sub)),
                          "floor_rate": float(sub.floor_rate.median()),
                          "base_rate": float(sub.base_rate.median()),
                          "ratio": float(sub.floor_ratio.median()),
                          "ci": [lo, hi]})
    results["floor_test_by_prevalence"] = band_rows

    # --- the decisive comparison: is the floor rate at the label's noise floor?
    if LABEL_FPR is not None:
        print()
        print(f"  IS THE FLOOR AT THE LABEL'S NOISE FLOOR?")
        print(f"    label FPR (curated non-essentials)   : {LABEL_FPR:.4f}")
        sel = F[F.base_rate < 0.10]        # bands where the gate demonstrably acts
        for nm, sub in [("all floor genes", F), ("gate-eligible (prev<10%)", sel)]:
            if not len(sub):
                continue
            fr = float(sub.floor_rate.median())
            resid = max(0.0, fr - LABEL_FPR)
            print(f"    {nm:<34} floor rate {fr:.4f}   "
                  f"residual above noise {resid:.4f}   "
                  f"({100*resid/float(sub.base_rate.median()):.1f}% of base)")
            results.setdefault("floor_vs_label_noise", {})[nm] = {
                "floor_rate": fr, "label_fpr": LABEL_FPR,
                "residual_above_noise": resid,
                "residual_as_pct_of_base":
                    100 * resid / float(sub.base_rate.median())}
        fr_sel = float(sel.floor_rate.median()) if len(sel) else np.nan
        if np.isfinite(fr_sel) and fr_sel <= LABEL_FPR * 1.5:
            print(f"    -> the floor rate is within ~1.5x the label's own false")
            print(f"       positive rate. Depletion at the floor is COMPLETE to")
            print(f"       within what this label can resolve; quoting it as a")
            print(f"       percentage depletion understates the gate.")
            results["floor_noise_verdict"] = "at_label_noise_floor"
        else:
            print(f"    -> the floor rate is well above the label's FPR, so most")
            print(f"       of it is real residual dependency, not screen noise.")
            results["floor_noise_verdict"] = "real_residual_signal"
else:
    print("  no gene has a floor block that large")
    results["floor_test"] = {"n_genes": 0}

# ============================================================ Step 7
print()
print("=" * 78)
print("STEP 7 -- VERDICT")
print("=" * 78)
depleted = d1z < 0.90 and wsr.pvalue < 0.05
flat_top = upper < 0.05 and abs(slope) < 0.01
gradient = results["partial_auc"]["ranker_region_Z"] > 0.52

if depleted and flat_top:
    v = "gate_works_auroc_was_the_wrong_metric"
elif depleted and gradient:
    v = "both_filter_and_ranker"
elif depleted:
    v = "gate_works_weakly"
else:
    v = "NO_DEPLETION__the_real_negative_result"
results["verdict"] = v

print(f"  decile-1 depletion   : ratio {d1z:.3f}, p {wsr.pvalue:.3g}  "
      f"-> {'DEPLETED' if depleted else 'NOT depleted'}")
print(f"  upper-decile profile : sd {upper:.3f}, slope {slope:+.4f}  "
      f"-> {'FLAT (filter-like)' if flat_top else 'sloped (ranker-like)'}")
print(f"  ranker region pAUC   : {results['partial_auc']['ranker_region_Z']:.4f}")
print(f"\n  VERDICT: {v}")

out = OUTPUTS / f"test_run_abundance_gate_{LABEL}_results.json"
out.write_text(json.dumps(results, indent=2, default=float), encoding="utf-8")
print(f"\nwrote {out}")
print(f"wrote {OUTPUTS / f'test_run_abundance_gate_{LABEL}_per_gene.parquet'}")
