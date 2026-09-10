"""
test_run_stratum_centring.py
----------------------------
Read-only diagnostic. Writes no production file.

QUESTION
--------
Three independent measurements say the protein layer costs discrimination on the
population the pipeline actually ships:

    RESTRICTED (two-layer lines only)   Stouffer vs RNA alone   +0.0011, p=0.23
    FULL       (strata pooled)          Stouffer vs RNA alone   -0.0026, p=0.014
    FULL, gate region (pAUC 0.8-1.0)    RNA 0.6176 vs Z 0.6009

The only mechanism consistent with all three is the coverage-selection lift, not
the protein data. In RESTRICTED every line carries both layers, strata cannot
mix, and protein is harmless. The deficit appears only when one-layer and
two-layer lines are ranked against each other -- and two-layer lines carry a
measured centre shift of about +0.21 Z (600-gene sweep: median mean(z_E | 2L)
= +0.151, positive in 81.5% of genes). Lines sent for proteomics are
systematically higher-expressing ones. That lift has no merit basis: it is a
fact about which cell lines a lab put through a mass spectrometer.

If that is the whole story, removing the shift should recover the deficit.

METHOD
------
Same construction as test_run_abundance_gate.py / test_run_stouffer_regime_split.py.
On the FULL population, per gene, compare four rankings:

    rna        z_E alone (covers every line, so this is the coverage ceiling)
    current    0.5*p_E + 0.5*p_P on two-layer lines, single percentile otherwise
    Z          Stouffer, uncentred                     <- what was proposed
    Z_cent     Stouffer, per-stratum centred           <- what is being tested

Centring subtracts, within each gene, the mean of Z over each n_layers stratum:

    Z_cent = Z - mean(Z | same n_layers, same gene)

The Stouffer divisor is left alone, so the scale argument is untouched -- this
removes a location shift, not a spread correction. Centring is applied within
gene because the shift is a property of which lines were profiled for that gene.

BINARY OUTCOME
--------------
    Z_cent >= rna   the coverage lift was the entire deficit. Centring is
                    justified on measurement rather than aesthetics, and the
                    protein layer is vindicated as neutral-to-helpful.
    Z_cent <  rna   something else is wrong with the two-layer stratum, and the
                    protein layer is a net negative on the shipped population.
                    It should then be reported as displayed evidence rather than
                    folded into the ranked score.

Reported against both GDSC drug sensitivity (where the deficit was measured) and
Chronos knockout fitness (the denser, causally closer label), because a result
that only holds on one of them is not a result.

Metrics: full AUROC and McClish-standardised partial AUROC over FPR in
[0.8, 1.0] -- the operating region established by test_run_abundance_gate.py,
where the abundance signal was shown to concentrate.

Outputs:
    src/pipeline/outputs/test_run_stratum_centring_results.json
    src/pipeline/outputs/test_run_stratum_centring_per_gene.parquet

Run:
    python src/pipeline/test_run_stratum_centring.py
"""
import json
import warnings
from pathlib import Path

import duckdb
import numpy as np
import pandas as pd
import pyarrow.parquet as pq
from scipy.stats import norm, rankdata, wilcoxon

warnings.filterwarnings("ignore", category=RuntimeWarning)

OUTPUTS = Path("src/pipeline/outputs")
REF_DIR = Path("reference")
DATA_CLEAN = Path("data/parquet/data_clean")
CLEANED_TRACK = Path("cleaned_track_data")

RHO_GLOBAL = 0.46
W_E, W_P = 0.987, 0.381
DEP_THRESHOLD = -0.5
MIN_POS = 10
MIN_NEG = 10
MIN_BOTH = 30
N_CHRONOS_GENES = 3000
SEED = 42

results = {}


def vdw(s):
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


def auroc(score, pos):
    ok = ~np.isnan(score)
    score, pos = score[ok], pos[ok]
    n1, n0 = int(pos.sum()), int((~pos).sum())
    if n1 < 1 or n0 < 1:
        return np.nan
    r = rankdata(score)
    return (r[pos].sum() - n1 * (n1 + 1) / 2.0) / (n1 * n0)


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


# ============================================================ Step 1
print("=" * 78)
print("STEP 1 -- layers and both label sets")
print("=" * 78)

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

gdsc = pd.read_parquet("validation/prepared/gdsc_scored_ready.parquet",
                       columns=["model_id", "target_ensg", "sensitive"])
gdsc["model_id"] = gdsc.model_id.str.lower()
gdsc["target_ensg"] = gdsc.target_ensg.astype("string").str.split(".").str[0].str.lower()
gdsc_pos = (gdsc[gdsc.sensitive].groupby("target_ensg")["model_id"]
            .apply(set).to_dict())

ch = pd.read_parquet("validation/prepared/chronos_long.parquet")
ch["ensg_id"] = ch.ensg_id.astype(str).str.split(".").str[0].str.lower()
ib = pd.read_parquet("validation/prepared/id_bridge.parquet")
ch = ch.merge(ib, on="sanger_model_id", how="inner")
ch["model_id"] = ch.model_id.str.lower()
rng0 = np.random.default_rng(SEED)
all_ch = sorted(ch.ensg_id.unique())
keep = set(gdsc_pos) & set(all_ch)
pool = [g for g in all_ch if g not in keep]
keep |= set(rng0.choice(pool, size=min(N_CHRONOS_GENES - len(keep), len(pool)),
                        replace=False))
ch = ch[ch.ensg_id.isin(keep)]
chronos_pos = (ch[ch.essentiality <= DEP_THRESHOLD].groupby("ensg_id")["model_id"]
               .apply(set).to_dict())
print(f"  GDSC genes {len(gdsc_pos)}   chronos genes {len(chronos_pos)}")

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

genes = sorted({g for g in (set(gdsc_pos) | set(chronos_pos))
                if g in expr_col and (g in ccle_by_ensg or g in procan_by_ensg)})
print(f"  genes with RNA + at least one protein column: {len(genes)}")

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
        return zE, pd.Series(dtype=float)
    idx = sorted(set().union(*[set(v.index) for v in parts.values()]))
    return zE, pd.DataFrame({k: v.reindex(idx)
                             for k, v in parts.items()}).mean(axis=1)


# ============================================================ Step 2
print()
print("=" * 78)
print("STEP 2 -- FULL population: uncentred vs per-stratum centred")
print("=" * 78)

rows = []
for g in genes:
    zE, zP = layers_for(g)
    if not len(zE) or not len(zP):
        continue
    idx = sorted(set(zE.index) | set(zP.index))
    d = pd.DataFrame({"zE": zE.reindex(idx), "zP": zP.reindex(idx)})
    d = d[d.notna().any(axis=1)]
    both = d.zE.notna() & d.zP.notna()
    if int(both.sum()) < MIN_BOTH:
        continue

    d["Z"] = stouffer(d.zE, d.zP)
    # Equal-weight Stouffer. The current rule is equal-weight and Stouffer as
    # proposed is 2.6:1 RNA-dominant, so `current` vs `Z` confounds the TRANSFORM
    # (percentile vs probit) with the WEIGHTS. This isolates the transform.
    d["Z_eq"] = stouffer(d.zE, d.zP, wE=1.0, wP=1.0)
    d["pE"], d["pP"] = norm.cdf(d.zE), norm.cdf(d.zP)
    d["current"] = np.where(both, 0.5 * d.pE + 0.5 * d.pP,
                            d.pE.fillna(0) + d.pP.fillna(0))
    # per-stratum centring: remove the location shift, leave the divisor alone
    shift = d.loc[both, "Z"].mean() - d.loc[~both, "Z"].mean()
    d["Z_cent"] = d["Z"] - both.map({True: d.loc[both, "Z"].mean(),
                                     False: d.loc[~both, "Z"].mean()})

    rec = {"ensg_id": g, "symbol": symbol_of.get(g, "?"),
           "n_lines": int(len(d)), "n_both": int(both.sum()),
           "stratum_shift": float(shift)}
    for lab, posmap in [("gdsc", gdsc_pos), ("chr", chronos_pos)]:
        if g not in posmap:
            continue
        pos = d.index.isin(posmap[g])
        if pos.sum() < MIN_POS or (~pos).sum() < MIN_NEG:
            continue
        for nm, col in [("rna", "zE"), ("cur", "current"),
                        ("Z", "Z"), ("Zeq", "Z_eq"), ("Zc", "Z_cent")]:
            rec[f"auc_{lab}_{nm}"] = auroc(d[col].values, pos)
            rec[f"pauc_{lab}_{nm}"] = partial_auc(d[col].values, pos)
    rows.append(rec)

G = pd.DataFrame(rows)
G.to_parquet(OUTPUTS / "test_run_stratum_centring_per_gene.parquet")
print(f"  genes scored: {len(G)}")
print(f"  measured stratum shift (2-layer mean - 1-layer mean), in Z units:")
print(f"    median {G.stratum_shift.median():+.4f}   "
      f"IQR [{G.stratum_shift.quantile(.25):+.3f}, "
      f"{G.stratum_shift.quantile(.75):+.3f}]   "
      f"positive in {100*(G.stratum_shift > 0).mean():.1f}% of genes")
results["stratum_shift"] = {
    "median": float(G.stratum_shift.median()),
    "frac_positive": float((G.stratum_shift > 0).mean()),
    "n_genes": int(len(G)),
}

# ============================================================ Step 3
for lab, name in [("gdsc", "GDSC drug sensitivity"),
                  ("chr", "Chronos knockout fitness")]:
    sub = G.dropna(subset=[f"auc_{lab}_Z"])
    if len(sub) < 10:
        continue
    print()
    print("=" * 78)
    print(f"STEP 3 -- {name}   (n = {len(sub)} genes)")
    print("=" * 78)
    print(f"  {'rule':<34} {'full AUROC':>11} {'pAUC 0.8-1.0':>14}")
    for nm, disp in [("rna", "RNA alone (coverage ceiling)"),
                     ("cur", "current 0.5E+0.5P (equal weight)"),
                     ("Zeq", "Stouffer Z, EQUAL weights"),
                     ("Z", "Stouffer Z, reliability weights"),
                     ("Zc", "Stouffer Z, per-stratum centred")]:
        print(f"  {disp:<34} {sub[f'auc_{lab}_{nm}'].median():>11.4f} "
              f"{sub[f'pauc_{lab}_{nm}'].median():>14.4f}")

    d_unc = (sub[f"auc_{lab}_Z"] - sub[f"auc_{lab}_rna"]).dropna()
    d_cen = (sub[f"auc_{lab}_Zc"] - sub[f"auc_{lab}_rna"]).dropna()
    d_gain = (sub[f"auc_{lab}_Zc"] - sub[f"auc_{lab}_Z"]).dropna()
    print()
    print(f"  uncentred Z  vs RNA alone : {d_unc.median():+.4f}   "
          f"p = {wilcoxon(d_unc).pvalue:.4g}")
    print(f"  centred   Z  vs RNA alone : {d_cen.median():+.4f}   "
          f"p = {wilcoxon(d_cen).pvalue:.4g}")
    print(f"  gain from centring        : {d_gain.median():+.4f}   "
          f"p = {wilcoxon(d_gain).pvalue:.4g}   "
          f"positive in {100*(d_gain > 0).mean():.1f}% of genes")
    # THE DECISIVE COMPARISON: does equal-weight Stouffer close the gate-region
    # gap to the mean rule? If so the reversal was about weights, not Stouffer.
    g_eq_cur = (sub[f"pauc_{lab}_Zeq"] - sub[f"pauc_{lab}_cur"]).dropna()
    g_eq_z = (sub[f"pauc_{lab}_Zeq"] - sub[f"pauc_{lab}_Z"]).dropna()
    print()
    print(f"  GATE REGION (pAUC 0.8-1.0):")
    print(f"    equal-weight Z  vs current : {g_eq_cur.median():+.4f}   "
          f"p = {wilcoxon(g_eq_cur).pvalue:.4g}   "
          f"positive in {100*(g_eq_cur > 0).mean():.1f}%")
    print(f"    equal-weight Z  vs reliability-weight Z : {g_eq_z.median():+.4f}   "
          f"p = {wilcoxon(g_eq_z).pvalue:.4g}")
    results.setdefault(f"{lab}_equal_weight", {}).update({
        "pauc_eq_vs_current": float(g_eq_cur.median()),
        "pauc_eq_vs_current_p": float(wilcoxon(g_eq_cur).pvalue),
        "pauc_eq_vs_reliability": float(g_eq_z.median()),
        "pauc_eq_vs_reliability_p": float(wilcoxon(g_eq_z).pvalue),
    })

    recovered = d_cen.median() >= 0
    print(f"\n  -> centred Z {'MEETS OR BEATS' if recovered else 'STILL TRAILS'} "
          f"RNA alone")
    results[f"{lab}_result"] = {
        "n_genes": int(len(sub)),
        "auc": {nm: float(sub[f"auc_{lab}_{nm}"].median())
                for nm in ["rna", "cur", "Z", "Zeq", "Zc"]},
        "pauc": {nm: float(sub[f"pauc_{lab}_{nm}"].median())
                 for nm in ["rna", "cur", "Z", "Zeq", "Zc"]},
        "uncentred_vs_rna": float(d_unc.median()),
        "uncentred_vs_rna_p": float(wilcoxon(d_unc).pvalue),
        "centred_vs_rna": float(d_cen.median()),
        "centred_vs_rna_p": float(wilcoxon(d_cen).pvalue),
        "centring_gain": float(d_gain.median()),
        "centring_gain_p": float(wilcoxon(d_gain).pvalue),
        "recovered": bool(recovered),
    }

# ============================================================ Step 4
print()
print("=" * 78)
print("STEP 4 -- VERDICT")
print("=" * 78)
rg = results.get("gdsc_result", {})
rc = results.get("chr_result", {})
both_rec = rg.get("recovered") and rc.get("recovered")
either = rg.get("recovered") or rc.get("recovered")
verdict = ("coverage_lift_explains_the_deficit" if both_rec else
           "partial__recovers_on_one_label_only" if either else
           "centring_does_not_recover__protein_is_a_net_negative")
results["verdict"] = verdict
print(f"  GDSC    : centred vs RNA {rg.get('centred_vs_rna', float('nan')):+.4f} "
      f"(was {rg.get('uncentred_vs_rna', float('nan')):+.4f})")
print(f"  Chronos : centred vs RNA {rc.get('centred_vs_rna', float('nan')):+.4f} "
      f"(was {rc.get('uncentred_vs_rna', float('nan')):+.4f})")
print(f"\n  VERDICT: {verdict}")

out = OUTPUTS / "test_run_stratum_centring_results.json"
out.write_text(json.dumps(results, indent=2, default=float), encoding="utf-8")
print(f"\nwrote {out}")
print(f"wrote {OUTPUTS / 'test_run_stratum_centring_per_gene.parquet'}")
