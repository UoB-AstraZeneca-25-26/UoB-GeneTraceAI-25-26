"""
test_run_stouffer_regime_split.py
---------------------------------
Read-only diagnostic. Writes no production file.

QUESTION
--------
Correlation-corrected weighted Stouffer Z was proposed to replace the current
two-layer rule in 02_core_score.ipynb Cell 8b (`0.5*E + 0.5*P`). The scale
argument for it is settled: Phi^-1 of a within-gene rank is standard normal, so
every n_layers stratum lands on one scale and `stratum_rank` becomes unnecessary
-- for dense genes. The measured suppression fix is real (two-layer share of the
top 20 moves 0.100 -> 0.350 against a 0.349 base rate, 600-gene sweep).

None of that says combining the layers PREDICTS BETTER. It says the combined
number is on a defensible scale. Those are different claims and only the second
one has been tested.

The open question is whether combining helps at all where the two layers do not
measure the same thing. RNA and protein are treated as two noisy readings of one
latent abundance. Per-gene rho(E,P) ranges from 0.04 to 0.72 across the genome
(600-gene sweep: median 0.418, 2.7% negative). At the low end they are plausibly
measuring two different biological quantities -- transcript level and, after
post-transcriptional regulation, functional protein level. Combining
measurements of two different things is not a statistics problem that a better
divisor fixes.

    If combining beats RNA-alone only for high-rho genes, rho is a fourth regime
    axis and the combination rule should be gated on it.
    If combining beats RNA-alone everywhere, the layers are complementary even
    when they disagree, and the single rule stands.
    If combining beats RNA-alone nowhere, the protein layer is not earning its
    place in the score at all, and that needs to be known.

METHOD
------
Ground truth: GDSC drug sensitivity, the same labels stage4_eval_save.py uses.
A (gene, cell line) pair is positive when the line is `sensitive` to a drug whose
annotated target is that gene. 141 target genes, 54-587 sensitive lines each.

Per gene, AUROC of each scoring rule at separating sensitive from non-sensitive
cell lines. Two populations, because they answer different questions:

  RESTRICTED  only lines carrying BOTH layers. Isolates the combination rule
              from coverage effects -- every line is scored the same way, so
              this is a clean RNA-alone vs protein-alone vs combined contest.
  FULL        every line carrying at least one layer, one-layer and two-layer
              pooled into a single ranking. This is what the output spec
              actually proposes to ship, and the only place the rho correction
              can change anything.

    IMPORTANT, and the reason both populations are needed: within the restricted
    population the Stouffer divisor is a positive constant, so Z is a monotone
    transform of (w_E*z_E + w_P*z_P) and AUROC is INVARIANT to rho. The rho
    correction is a cross-stratum calibration device, not a discrimination
    improvement. Any AUROC difference it produces can only appear in FULL, by
    re-ordering one-layer lines against two-layer lines. This is asserted in the
    maths and verified numerically in Step 5.

Genes are then split at the median per-gene rho and the AUROC lift of combining
over RNA-alone is compared across halves (Mann-Whitney, plus a permutation test
on the difference of medians). The same split is repeated on the existing
`class` regime from gene_regime.parquet, since abundance_tracking is the regime
that already claims RNA tracks function.

WEIGHTS
-------
The current weights (w_E=0.987, w_P=0.381) are Q5 measurement reliabilities --
they encode which layer is measured more precisely, not which layer predicts
sensitivity better. Only the second is what the score is for. Step 6 sweeps
w_P/w_E and reports the AUROC-optimal ratio, fitted on TRAIN genes and reported
on TEST genes using the same rng(42) split as stage4_eval_save.py, so the
reported number is not the fitted one.

NORMALISATION
-------------
van der Waerden normal scores, Phi^-1(rank / (n+1)), with method="min" to match
the tie rule of Cell 7b. This cannot emit +/-inf, unlike percentile-then-clip
(100% of genes in the 600-gene sweep produced pct == 1.0 from a unique maximum).
AUROC is rank-based, so the transform choice cannot affect single-layer results
and affects combined results only through the relative spacing of the layers.

SIGN CONVENTION
---------------
Higher score = stronger candidate. GDSC `sensitive` = True is the positive class.
An AUROC of 0.5 is chance; below 0.5 means the score is anti-correlated with
sensitivity, which is a real finding for a gene, not a bug to be flipped.

Outputs:
    src/pipeline/outputs/test_run_stouffer_regime_split_results.json
    src/pipeline/outputs/test_run_stouffer_regime_split_per_gene.parquet

Run:
    python src/pipeline/test_run_stouffer_regime_split.py
"""
import json
import warnings
from pathlib import Path

import duckdb
import numpy as np
import pandas as pd
import pyarrow.parquet as pq
from scipy.stats import norm, rankdata, spearmanr, mannwhitneyu, wilcoxon

warnings.filterwarnings("ignore", category=RuntimeWarning)

OUTPUTS = Path("src/pipeline/outputs")
REF_DIR = Path("reference")
DATA_CLEAN = Path("data/parquet/data_clean")
CLEANED_TRACK = Path("cleaned_track_data")

RHO_GLOBAL = 0.46          # 02_core_score Cell 8b / scoring_variants.combine_three
W_E, W_P = 0.987, 0.381    # Q5 reliabilities: RNA technical r, protein cross-platform rho
MIN_POS = 10               # minimum sensitive lines in a population to score a gene
MIN_NEG = 10
MIN_BOTH = 30              # minimum two-layer lines before a per-gene rho is trusted
SEED = 42                  # same split seed as stage4_eval_save.py
N_PERM = 10_000

results = {}


# ============================================================ helpers
def vdw(s: pd.Series) -> pd.Series:
    """van der Waerden normal scores. method='min' matches Cell 7b's tie rule."""
    n = int(s.notna().sum())
    return pd.Series(norm.ppf(s.rank(method="min") / (n + 1)), index=s.index)


def auroc(score: np.ndarray, pos: np.ndarray) -> float:
    """Mann-Whitney AUROC with correct mid-rank handling of ties."""
    ok = ~np.isnan(score)
    score, pos = score[ok], pos[ok]
    n1, n0 = int(pos.sum()), int((~pos).sum())
    if n1 < 1 or n0 < 1:
        return np.nan
    r = rankdata(score)
    return (r[pos].sum() - n1 * (n1 + 1) / 2.0) / (n1 * n0)


def stouffer(zE: pd.Series, zP: pd.Series, rho: float,
             wE: float = W_E, wP: float = W_P) -> pd.Series:
    """Correlation-corrected weighted Stouffer Z over the layers actually present.

    One layer  -> w*z / sqrt(w^2) = z, so the stratum is standard normal.
    Two layers -> (wE*zE + wP*zP) / sqrt(wE^2 + wP^2 + 2*wE*wP*rho).
    """
    out = np.full(len(zE), np.nan)
    e_ok, p_ok = zE.notna().values, zP.notna().values
    both = e_ok & p_ok
    ze, zp = zE.values, zP.values
    out[e_ok & ~p_ok] = ze[e_ok & ~p_ok]
    out[p_ok & ~e_ok] = zp[p_ok & ~e_ok]
    out[both] = (wE * ze[both] + wP * zp[both]) / np.sqrt(
        wE ** 2 + wP ** 2 + 2 * wE * wP * rho)
    return pd.Series(out, index=zE.index)


def perm_test_median_diff(a: np.ndarray, b: np.ndarray, n_perm=N_PERM, seed=SEED):
    """Two-sided permutation test on the difference of medians."""
    a, b = a[~np.isnan(a)], b[~np.isnan(b)]
    if len(a) < 3 or len(b) < 3:
        return np.nan, np.nan
    obs = np.median(a) - np.median(b)
    pool = np.concatenate([a, b])
    rng = np.random.default_rng(seed)
    na = len(a)
    cnt = 0
    for _ in range(n_perm):
        rng.shuffle(pool)
        if abs(np.median(pool[:na]) - np.median(pool[na:])) >= abs(obs):
            cnt += 1
    return float(obs), float((cnt + 1) / (n_perm + 1))


# ============================================================ Step 1
print("=" * 78)
print("STEP 1 -- GDSC sensitivity labels and the candidate gene set")
print("=" * 78)

gdsc = pd.read_parquet("validation/prepared/gdsc_scored_ready.parquet",
                       columns=["model_id", "target_ensg", "sensitive"])
gdsc["model_id"] = gdsc["model_id"].str.lower()
gdsc["target_ensg"] = gdsc["target_ensg"].astype("string").str.split(".").str[0].str.lower()

sens = (gdsc[gdsc["sensitive"]][["target_ensg", "model_id"]]
        .drop_duplicates())
sens_by_gene = sens.groupby("target_ensg")["model_id"].apply(set).to_dict()
target_genes = sorted(sens_by_gene)
print(f"  GDSC target genes with >=1 sensitive line : {len(target_genes)}")
print(f"  sensitive (gene, line) pairs              : {len(sens):,}")

regime = pd.read_parquet(OUTPUTS / "gene_regime.parquet", columns=["ensg_id", "class"])
regime["ensg_id"] = regime["ensg_id"].astype(str).str.lower()
class_of = dict(zip(regime.ensg_id, regime["class"]))
print(f"  gene_regime classes available for         : "
      f"{len(set(target_genes) & set(class_of))} of them")

# ============================================================ Step 2
print()
print("=" * 78)
print("STEP 2 -- rebuild the RNA and protein layers for those genes")
print("=" * 78)

gl = pd.read_parquet(REF_DIR / "gene_lookup.parquet",
                     columns=["ensg_id", "hgnc_symbol", "biotype",
                              "hgnc_status", "uniprot_ids"])
uni = gl[(gl.biotype == "protein_coding") & (gl.hgnc_status == "Approved")].copy()
uni["ensg_id"] = uni.ensg_id.astype("string").str.split(".").str[0].str.lower()
valid_ensg = set(uni.ensg_id)
symbol_of = dict(zip(uni.ensg_id, uni.hgnc_symbol))

# UniProt -> ENSG, multi-accession aware (02_core_score Cell 3b)
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
    e = u2e.get(a) or u2e.get(a.split("-")[0])          # isoform suffix, Cell 3b
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

genes = [g for g in target_genes
         if g in expr_col and (g in ccle_by_ensg or g in procan_by_ensg)]
print(f"  target genes with an RNA column           : "
      f"{sum(g in expr_col for g in target_genes)}")
print(f"  ... and at least one protein column       : {len(genes)}")

# RNA
prof = pd.read_parquet(REF_DIR / "depmap_profiles.parquet")
rna_prof = prof[prof.datatype == "rna"][["profileid", "modelid"]].copy()
rna_prof["model_id"] = rna_prof.modelid.str.lower()          # rule 1
rna_prof["profileid"] = rna_prof.profileid.astype(str)

E = pq.read_table(DATA_CLEAN / "depmap_expr_clean.parquet",
                  columns=["index"] + [expr_col[g] for g in genes]).to_pandas()
if "index" not in E.columns:
    E = E.reset_index()
E["profileid"] = E["index"].astype(str)
E = E.drop(columns=["index"]).merge(rna_prof[["profileid", "model_id"]],
                                    on="profileid", how="inner")
E = E.drop(columns=["profileid"]).groupby("model_id").mean()   # Cell 7 dedup_mean
E.columns = [c.split(".")[0].lower() for c in E.columns]
print(f"  RNA matrix    : {E.shape[0]} lines x {E.shape[1]} genes")

# protein, both platforms
need_ccle = sorted({a for g in genes for a in ccle_by_ensg.get(g, [])})
C = pq.read_table(CLEANED_TRACK / "proteomics.parquet",
                  columns=[key_col] + need_ccle).to_pandas()
C.index = C[key_col].str.lower().values                        # rule 1
C = C.drop(columns=[key_col]).astype(float).groupby(level=0).mean()
print(f"  CCLE protein  : {C.shape[0]} lines x {C.shape[1]} accessions")

need_pc = sorted({a for g in genes for a in procan_by_ensg.get(g, [])})
if need_pc:
    q = ", ".join(f'"{a}"' for a in need_pc)
    P = con.execute(f"SELECT model_id, {q} FROM procan_proteomics "
                    "WHERE model_id IS NOT NULL").df()
    P.index = P.model_id.str.lower().values
    P = P.drop(columns=["model_id"]).astype(float).groupby(level=0).mean()
else:
    P = pd.DataFrame()
con.close()
print(f"  ProCan protein: {P.shape[0]} lines x {P.shape[1]} accessions")

results["inputs"] = {
    "gdsc_target_genes": len(target_genes),
    "genes_with_rna_and_protein": len(genes),
    "rna_lines": int(E.shape[0]),
    "ccle_lines": int(C.shape[0]),
    "procan_lines": int(P.shape[0]),
    "rho_global": RHO_GLOBAL,
    "weights": {"rna": W_E, "protein": W_P},
}


# ============================================================ Step 3
print()
print("=" * 78)
print("STEP 3 -- per-gene AUROC, restricted (both layers) and full (>=1 layer)")
print("=" * 78)


def layers_for(g):
    """Return (z_rna, z_protein) as van der Waerden scores, indexed by model_id."""
    e = E[g].dropna()
    zE = vdw(e) if len(e) else pd.Series(dtype=float)

    parts = {}
    accs = [a for a in ccle_by_ensg.get(g, []) if a in C.columns]
    if accs:
        s = C[accs].mean(axis=1).dropna()               # average duplicate accessions
        if len(s):
            parts["ccle"] = vdw(s)
    accs = [a for a in procan_by_ensg.get(g, []) if a in P.columns]
    if accs:
        s = P[accs].mean(axis=1).dropna()
        if len(s):
            parts["procan"] = vdw(s)
    if not parts:
        return zE, pd.Series(dtype=float)
    idx = sorted(set().union(*[set(v.index) for v in parts.values()]))
    # Cell 7e merge rule: mean where both platforms cover, whichever exists otherwise.
    # Applied on normal scores rather than percentiles -- monotone equivalent per
    # platform, and it keeps one protein layer (n_layers is unaffected).
    zP = pd.DataFrame({k: v.reindex(idx) for k, v in parts.items()}).mean(axis=1)
    return zE, zP


rows = []
for g in genes:
    zE, zP = layers_for(g)
    if not len(zE) or not len(zP):
        continue
    idx = sorted(set(zE.index) | set(zP.index))
    d = pd.DataFrame({"zE": zE.reindex(idx), "zP": zP.reindex(idx)})
    d = d[d.notna().any(axis=1)]
    d["pos"] = d.index.isin(sens_by_gene[g])

    both = d.zE.notna() & d.zP.notna()
    if int(both.sum()) < MIN_BOTH:
        continue
    rho_g = spearmanr(d.loc[both, "zE"], d.loc[both, "zP"]).statistic
    rho_use = float(np.clip(rho_g, -0.90, 0.99))

    # percentile form of the current rule, for a like-for-like comparison
    d["pE"] = norm.cdf(d.zE)
    d["pP"] = norm.cdf(d.zP)
    d["current"] = np.where(both, 0.5 * d.pE + 0.5 * d.pP,
                            d.pE.fillna(0) + d.pP.fillna(0))
    d["Z46"] = stouffer(d.zE, d.zP, RHO_GLOBAL)
    d["Zg"] = stouffer(d.zE, d.zP, rho_use)
    d["Zeq"] = stouffer(d.zE, d.zP, RHO_GLOBAL, 1.0, 1.0)

    rec = {"ensg_id": g, "symbol": symbol_of.get(g, "?"),
           "class": class_of.get(g), "rho": float(rho_g),
           "n_lines_full": int(len(d)), "n_both": int(both.sum()),
           "n_rna_only": int((d.zE.notna() & d.zP.isna()).sum()),
           "n_prot_only": int((d.zP.notna() & d.zE.isna()).sum())}

    # -------- RESTRICTED: two-layer lines only
    r = d[both]
    pos = r.pos.values
    rec["n_pos_restricted"] = int(pos.sum())
    rec["n_neg_restricted"] = int((~pos).sum())
    if pos.sum() >= MIN_POS and (~pos).sum() >= MIN_NEG:
        rec["auc_rna_r"] = auroc(r.zE.values, pos)
        rec["auc_prot_r"] = auroc(r.zP.values, pos)
        rec["auc_Z46_r"] = auroc(r.Z46.values, pos)
        rec["auc_Zg_r"] = auroc(r.Zg.values, pos)
        rec["auc_Zeq_r"] = auroc(r.Zeq.values, pos)
        rec["auc_cur_r"] = auroc(r.current.values, pos)
        rec["lift_r"] = rec["auc_Z46_r"] - rec["auc_rna_r"]

    # -------- FULL: every line with >=1 layer, one ranking
    pos = d.pos.values
    rec["n_pos_full"] = int(pos.sum())
    rec["n_neg_full"] = int((~pos).sum())
    if pos.sum() >= MIN_POS and (~pos).sum() >= MIN_NEG:
        rec["auc_rna_f"] = auroc(d.zE.values, pos)      # NaN protein lines drop out
        rec["auc_Z46_f"] = auroc(d.Z46.values, pos)
        rec["auc_Zg_f"] = auroc(d.Zg.values, pos)
        rec["auc_cur_f"] = auroc(d.current.values, pos)
        rec["lift_f"] = rec["auc_Z46_f"] - rec["auc_rna_f"]
        rec["vs_current_f"] = rec["auc_Z46_f"] - rec["auc_cur_f"]
    rows.append(rec)

G = pd.DataFrame(rows)
G.to_parquet(OUTPUTS / "test_run_stouffer_regime_split_per_gene.parquet")
Gr = G.dropna(subset=["auc_Z46_r"]).copy()
Gf = G.dropna(subset=["auc_Z46_f"]).copy()
print(f"  genes scored (restricted population): {len(Gr)}")
print(f"  genes scored (full population)      : {len(Gf)}")
print(f"  median two-layer lines per gene     : {G.n_both.median():.0f}")
print(f"  per-gene rho: median {G.rho.median():.3f}  "
      f"IQR [{G.rho.quantile(.25):.3f}, {G.rho.quantile(.75):.3f}]  "
      f"negative {100*(G.rho < 0).mean():.1f}%")

# ============================================================ Step 4
print()
print("=" * 78)
print("STEP 4 -- does combining beat RNA alone?")
print("=" * 78)


def line(label, s):
    s = pd.Series(s).dropna()
    return (f"  {label:34s} median {s.median():.4f}   mean {s.mean():.4f}   "
            f">0.5 in {100*(s > 0.5).mean():5.1f}% of genes")


print("RESTRICTED (two-layer lines only -- clean contest, no coverage effects)")
for lab, col in [("RNA alone", "auc_rna_r"), ("protein alone", "auc_prot_r"),
                 ("current 0.5E+0.5P", "auc_cur_r"),
                 ("Stouffer Z (rho=0.46)", "auc_Z46_r"),
                 ("Stouffer Z (per-gene rho)", "auc_Zg_r"),
                 ("Stouffer Z (equal weights)", "auc_Zeq_r")]:
    print(line(lab, Gr[col]))
lift_r = (Gr.auc_Z46_r - Gr.auc_rna_r).dropna()
# Paired: both AUROCs come from the same gene on the same lines, so the signed-rank
# test on the per-gene difference is the correct one. An unpaired test here would
# throw away the pairing and badly understate power.
w_r = wilcoxon(Gr.auc_Z46_r.dropna(), Gr.auc_rna_r.dropna())
w_rc = wilcoxon(Gr.auc_Z46_r.dropna(), Gr.auc_cur_r.dropna())
print(f"\n  lift of Stouffer over RNA alone: median {lift_r.median():+.4f}   "
      f"positive in {100*(lift_r > 0).mean():.1f}% of genes")
print(f"  Wilcoxon signed-rank p (paired by gene) = {w_r.pvalue:.4g}")
print(f"  vs the CURRENT rule: median {(Gr.auc_Z46_r - Gr.auc_cur_r).median():+.4f}"
      f"   signed-rank p = {w_rc.pvalue:.4g}")

print()
print("FULL (one-layer and two-layer lines pooled -- what the output spec ships)")
for lab, col in [("RNA alone", "auc_rna_f"), ("current 0.5E+0.5P", "auc_cur_f"),
                 ("Stouffer Z (rho=0.46)", "auc_Z46_f"),
                 ("Stouffer Z (per-gene rho)", "auc_Zg_f")]:
    print(line(lab, Gf[col]))
lift_f = (Gf.auc_Z46_f - Gf.auc_rna_f).dropna()
w_f = wilcoxon(Gf.auc_Z46_f.dropna(), Gf.auc_rna_f.dropna())
w_fc = wilcoxon(Gf.auc_Z46_f.dropna(), Gf.auc_cur_f.dropna())
print(f"\n  lift of Stouffer over RNA alone: median {lift_f.median():+.4f}   "
      f"positive in {100*(lift_f > 0).mean():.1f}% of genes"
      f"   signed-rank p = {w_f.pvalue:.4g}")
print(f"  lift over the CURRENT rule     : median "
      f"{(Gf.auc_Z46_f - Gf.auc_cur_f).median():+.4f}   positive in "
      f"{100*((Gf.auc_Z46_f - Gf.auc_cur_f) > 0).mean():.1f}% of genes"
      f"   signed-rank p = {w_fc.pvalue:.4g}")

results["auroc"] = {
    "n_genes_restricted": int(len(Gr)), "n_genes_full": int(len(Gf)),
    "restricted": {c: float(Gr[c].median()) for c in
                   ["auc_rna_r", "auc_prot_r", "auc_cur_r", "auc_Z46_r",
                    "auc_Zg_r", "auc_Zeq_r"] if c in Gr},
    "full": {c: float(Gf[c].median()) for c in
             ["auc_rna_f", "auc_cur_f", "auc_Z46_f", "auc_Zg_f"] if c in Gf},
    "lift_restricted_median": float(lift_r.median()),
    "lift_restricted_frac_positive": float((lift_r > 0).mean()),
    "lift_full_median": float(lift_f.median()),
    "lift_full_frac_positive": float((lift_f > 0).mean()),
    "vs_current_full_median": float((Gf.auc_Z46_f - Gf.auc_cur_f).median()),
    "wilcoxon_p": {
        "restricted_vs_rna": float(w_r.pvalue),
        "restricted_vs_current": float(w_rc.pvalue),
        "full_vs_rna": float(w_f.pvalue),
        "full_vs_current": float(w_fc.pvalue),
    },
}

# ============================================================ Step 5
print()
print("=" * 78)
print("STEP 5 -- rho invariance check (the claim the method rests on)")
print("=" * 78)
inv = (Gr.auc_Z46_r - Gr.auc_Zg_r).abs()
print(f"  RESTRICTED: max |AUROC(rho=0.46) - AUROC(per-gene rho)| = {inv.max():.2e}")
print(f"    -> confirms rho cannot change within-stratum discrimination.")
print(f"       The divisor is a positive constant on a fixed layer set, so it")
print(f"       rescales Z without reordering. rho is a CALIBRATION parameter.")
invf = (Gf.auc_Z46_f - Gf.auc_Zg_f).abs()
print(f"  FULL      : median |difference| = {invf.median():.4f}   "
      f"max = {invf.max():.4f}")
print(f"    -> here rho DOES matter: it sets where two-layer lines sit relative")
print(f"       to one-layer lines. This is the only place it can pay off.")
results["rho_invariance"] = {
    "restricted_max_abs_diff": float(inv.max()),
    "full_median_abs_diff": float(invf.median()),
    "full_max_abs_diff": float(invf.max()),
}

# ============================================================ Step 6
print()
print("=" * 78)
print("STEP 6 -- THE SPLIT: does the lift depend on rho?")
print("=" * 78)

for pop, Gx, lift_col, auc_c, auc_r in [
        ("RESTRICTED", Gr, "lift_r", "auc_Z46_r", "auc_rna_r"),
        ("FULL", Gf, "lift_f", "auc_Z46_f", "auc_rna_f")]:
    d = Gx.dropna(subset=[lift_col]).copy()
    med = d.rho.median()
    hi = d[d.rho >= med][lift_col].values
    lo = d[d.rho < med][lift_col].values
    obs, pperm = perm_test_median_diff(hi, lo)
    u = mannwhitneyu(hi[~np.isnan(hi)], lo[~np.isnan(lo)])
    print(f"\n{pop}  (split at median rho = {med:.3f})")
    print(f"  high-rho genes n={len(hi):>3}  median lift {np.nanmedian(hi):+.4f}  "
          f"positive in {100*np.nanmean(hi > 0):.1f}%")
    print(f"  low-rho  genes n={len(lo):>3}  median lift {np.nanmedian(lo):+.4f}  "
          f"positive in {100*np.nanmean(lo > 0):.1f}%")
    print(f"  difference of medians {obs:+.4f}   permutation p = {pperm:.4f}   "
          f"Mann-Whitney p = {u.pvalue:.4f}")
    sp = spearmanr(d.rho, d[lift_col])
    print(f"  Spearman(rho, lift) = {sp.statistic:+.3f}  p = {sp.pvalue:.4f}")
    results[f"rho_split_{pop.lower()}"] = {
        "median_rho": float(med), "n_high": int(len(hi)), "n_low": int(len(lo)),
        "median_lift_high": float(np.nanmedian(hi)),
        "median_lift_low": float(np.nanmedian(lo)),
        "diff_of_medians": obs, "perm_p": pperm,
        "mannwhitney_p": float(u.pvalue),
        "spearman_rho_vs_lift": float(sp.statistic),
        "spearman_p": float(sp.pvalue),
    }

# quartile view -- a median split hides a monotone trend
print("\n  rho QUARTILES (restricted population):")
d = Gr.dropna(subset=["lift_r"]).copy()
d["q"] = pd.qcut(d.rho, 4, labels=["Q1 low", "Q2", "Q3", "Q4 high"])
qa = d.groupby("q", observed=True).agg(
    genes=("ensg_id", "size"), rho=("rho", "median"),
    auc_rna=("auc_rna_r", "median"), auc_prot=("auc_prot_r", "median"),
    auc_Z=("auc_Z46_r", "median"), lift=("lift_r", "median"))
print(qa.round(4).to_string())
results["rho_quartiles_restricted"] = json.loads(qa.reset_index().to_json(orient="records"))

# ============================================================ Step 7
print()
print("=" * 78)
print("STEP 7 -- the same split on the EXISTING regime class")
print("=" * 78)
d = Gr.dropna(subset=["lift_r"])
d = d[d["class"].notna()]
if d["class"].nunique() > 1:
    ca = d.groupby("class", observed=True).agg(
        genes=("ensg_id", "size"), rho=("rho", "median"),
        auc_rna=("auc_rna_r", "median"), auc_prot=("auc_prot_r", "median"),
        auc_Z=("auc_Z46_r", "median"), lift=("lift_r", "median"))
    print(ca.round(4).to_string())
    grp = [d[d["class"] == c]["lift_r"].values for c in ca.index]
    if len(grp) == 2:
        obs, pperm = perm_test_median_diff(grp[0], grp[1])
        print(f"\n  {ca.index[0]} vs {ca.index[1]}: difference of median lift "
              f"{obs:+.4f}, permutation p = {pperm:.4f}")
        results["class_split"] = {"diff_of_medians": obs, "perm_p": pperm}
    results["class_table"] = json.loads(ca.reset_index().to_json(orient="records"))
else:
    print("  only one regime class present among scored genes -- split not run")

# ============================================================ Step 8
print()
print("=" * 78)
print("STEP 8 -- weight sweep, fitted on TRAIN genes, reported on TEST genes")
print("=" * 78)
print("  Reliability weights answer 'which layer is measured better'.")
print("  This asks 'which layer predicts better' -- the thing the score is for.")

curated = sorted(Gr.ensg_id)
rng = np.random.default_rng(SEED)                    # same split as stage4_eval_save.py
test_g = set(rng.choice(curated, size=max(1, len(curated) // 5), replace=False))
train_g = [g for g in curated if g not in test_g]
print(f"\n  train genes {len(train_g)}   test genes {len(test_g)}")

RATIOS = [0.0, 0.1, 0.2, 0.386, 0.6, 1.0, 1.5, 2.5, 100.0]   # w_P / w_E
cache = {}
for g in Gr.ensg_id:
    zE, zP = layers_for(g)
    idx = sorted(set(zE.index) | set(zP.index))
    d = pd.DataFrame({"zE": zE.reindex(idx), "zP": zP.reindex(idx)})
    d = d[d.zE.notna() & d.zP.notna()]
    cache[g] = (d.zE.values, d.zP.values, d.index.isin(sens_by_gene[g]))

sweep = []
for r in RATIOS:
    wE, wP = (0.0, 1.0) if r == 100.0 else (1.0, r)
    tr = [auroc(wE * a + wP * b, p) for g in train_g for a, b, p in [cache[g]]]
    te = [auroc(wE * a + wP * b, p) for g in test_g for a, b, p in [cache[g]]]
    sweep.append({"w_ratio": r, "train_median_auc": float(np.nanmedian(tr)),
                  "test_median_auc": float(np.nanmedian(te))})
    tag = "  <- current (0.381/0.987)" if r == 0.386 else ""
    tag += "  <- RNA only" if r == 0.0 else ""
    tag += "  <- protein only" if r == 100.0 else ""
    print(f"    w_P/w_E = {r:>6.3f}   train {sweep[-1]['train_median_auc']:.4f}   "
          f"test {sweep[-1]['test_median_auc']:.4f}{tag}")

best = max(sweep, key=lambda s: s["train_median_auc"])
cur = [s for s in sweep if s["w_ratio"] == 0.386][0]
print(f"\n  best ratio on TRAIN            : {best['w_ratio']:.3f} "
      f"(train {best['train_median_auc']:.4f})")
print(f"  its held-out TEST median AUROC : {best['test_median_auc']:.4f}")
print(f"  current ratio on TEST          : {cur['test_median_auc']:.4f}")
print(f"  honest gain from refitting     : "
      f"{best['test_median_auc'] - cur['test_median_auc']:+.4f}")
results["weight_sweep"] = {
    "ratios": sweep, "n_train": len(train_g), "n_test": len(test_g),
    "best_train_ratio": best["w_ratio"],
    "best_test_auc": best["test_median_auc"],
    "current_test_auc": cur["test_median_auc"],
    "honest_gain": best["test_median_auc"] - cur["test_median_auc"],
}

# ============================================================ Step 9
print()
print("=" * 78)
print("STEP 9 -- VERDICT")
print("=" * 78)

sp_r = results["rho_split_restricted"]
lift_med = results["auroc"]["lift_restricted_median"]
lift_pos = results["auroc"]["lift_restricted_frac_positive"]

if lift_med > 0.005 and lift_pos > 0.55:
    combine = "combining_helps"
elif lift_med < -0.005:
    combine = "combining_hurts"
else:
    combine = "combining_neutral"

if sp_r["perm_p"] < 0.05 and abs(sp_r["diff_of_medians"]) > 0.005:
    gate = "rho_is_a_regime_axis"
elif sp_r["spearman_p"] < 0.05:
    gate = "rho_trend_present_but_weak"
else:
    gate = "no_rho_dependence_detected"

results["verdict"] = {"combination": combine, "rho_gating": gate}
print(f"  combination : {combine}")
print(f"                median lift {lift_med:+.4f}, positive in {lift_pos:.1%} of genes")
print(f"  rho gating  : {gate}")
print(f"                high-rho lift {sp_r['median_lift_high']:+.4f} vs "
      f"low-rho {sp_r['median_lift_low']:+.4f}, perm p {sp_r['perm_p']:.4f}")
print(f"\n  NOTE: rho cannot affect discrimination within a stratum (Step 5).")
print(f"  A verdict of '{gate}' concerns whether the two layers are measuring")
print(f"  the same thing, not whether the divisor is set correctly.")

out = OUTPUTS / "test_run_stouffer_regime_split_results.json"
out.write_text(json.dumps(results, indent=2, default=float), encoding="utf-8")
print(f"\nwrote {out}")
print(f"wrote {OUTPUTS / 'test_run_stouffer_regime_split_per_gene.parquet'}")
