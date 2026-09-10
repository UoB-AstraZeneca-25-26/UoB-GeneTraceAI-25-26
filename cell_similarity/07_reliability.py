"""
cell_similarity/07_reliability.py
---------------------------------
POWER OR CONFOUND?  Turn the 0.25x ratio into a number.

THE ARGUMENT
------------
`graph - tissue` is a difference of correlations between two dependency
profiles. If each profile carries reliability r_xx, the observed correlation is
attenuated by r_xx, and so is the difference. Therefore

    delta_true  ~=  delta_observed / r_xx

and the CRISPR-vs-RNAi comparison has an exact interpretation rather than a
rhetorical one:

    RNAi reliability   disattenuated delta   vs CRISPR 0.080
    0.90               0.022                 0.28x
    0.50               0.040                 0.50x
    0.35               0.057                 0.72x
    0.25               0.080                 1.00x   <- ratio is ALL noise

So: if RNAi's reliability is near 0.25, the 0.25x ratio is entirely measurement
error, there is no confound left to explain, and the component passes outright.

THREE ESTIMATES OF RNAi RELIABILITY, PLUS A POWER CONTROL
----------------------------------------------------------
  A. TEST-RETEST ACROSS INDEPENDENT SCREENS.  DEMETER2 is fitted to separate
     source datasets (Achilles shRNA, Novartis DRIVE). Lines screened in BOTH
     give a direct test-retest reliability for the same perturbation type. This
     is the cleanest estimate because it holds perturbation constant -- unlike
     CRISPR-vs-RNAi, which confounds reliability with mechanism.

     Reliability of a PROFILE is what matters here, not of a single gene, so it
     is measured as the correlation between a line's Achilles profile and its
     DRIVE profile over the same genes. Spearman-Brown is then used to express
     it on the scale of the combined (two-screen) estimate that the main
     analysis actually uses.

  B. THE PUBLISHED POSTERIOR SDs.  DEMETER2 emits per (gene, line) posterior
     standard deviations. Restricting to low-SD observations should RAISE the
     observed delta if attenuation is the explanation, and leave it flat if a
     confound is.

  C. POWER CONTROL.  CRISPR restricted to the same 659 lines as the RNAi arm.
     This separates "fewer lines" from "different perturbation" at zero cost. If
     CRISPR at n=659 still gives ~0.080, the gap is not sample size.

Read-only. Writes only into cell_similarity/outputs/.

Run:
    python cell_similarity/07_reliability.py
"""
import json
from pathlib import Path

import numpy as np
import pandas as pd

import common as C

D2 = C.GA.ROOT / "data" / "DEMETER2"
res = {}
CRISPR_REF = 0.0802          # RNA vs CRISPR, common essentials removed (05)

C.banner("POWER OR CONFOUND -- ESTIMATING RNAi RELIABILITY")

_, _, _, sym2ensg, _ = C.GA.load_gene_lookup()
lin = C.load_lineage()
ce = C.common_essentials(sym2ensg, "inferred")
S_rna = C.load_axis_similarity("rna")


def load_d2_matrix(fname):
    """DEMETER2 per-screen matrix: genes x CCLE-name -> lines x ensg, ModelID."""
    f = D2 / fname
    if not f.exists():
        return None
    df = pd.read_csv(f, index_col=0, low_memory=False).T
    W = C.GA._wide_from_symbol_matrix(df, sym2ensg)
    md = pd.read_csv(C.GA.DEPMAP_24Q4 / "Model.csv", low_memory=False)
    cc = "CCLEName" if "CCLEName" in md.columns else "CCLE_Name"
    bridge = {str(c).strip().lower(): str(m).lower()
              for c, m in zip(md[cc], md["ModelID"]) if isinstance(c, str)}
    mapped = [bridge.get(x) for x in W.index]
    ok = [i for i, v in enumerate(mapped) if v]
    W = W.iloc[ok]
    W.index = [mapped[i] for i in ok]
    return W.groupby(level=0).mean()


# ------------------------------------------------------------ A
C.banner("7A -- TEST-RETEST ACROSS TWO INDEPENDENT RNAi SCREENS")
A = load_d2_matrix("D2_Achilles_gene_dep_scores.csv")
Dr = load_d2_matrix("D2_DRIVE_gene_dep_scores.csv")
if A is None or Dr is None:
    print("  per-screen DEMETER2 files not present -- skipping")
    res["test_retest"] = None
else:
    shared_lines = sorted(set(A.index) & set(Dr.index))
    shared_genes = sorted(set(A.columns) & set(Dr.columns) - ce)
    print(f"  Achilles shRNA : {A.shape[0]} lines x {A.shape[1]:,} genes")
    print(f"  Novartis DRIVE : {Dr.shape[0]} lines x {Dr.shape[1]:,} genes")
    print(f"  screened in BOTH: {len(shared_lines)} lines, "
          f"{len(shared_genes):,} shared genes (common essentials removed)")
    if len(shared_lines) < 20:
        print("  too few doubly-screened lines")
        res["test_retest"] = None
    else:
        Aa = A.loc[shared_lines, shared_genes]
        Dd = Dr.loc[shared_lines, shared_genes]
        # variable genes only -- matching the main analysis's selective subset
        sd = pd.concat([Aa, Dd]).std(axis=0)
        keep = sd[sd > sd.median()].index
        Aa, Dd = Aa[keep], Dd[keep]
        rel = []
        for m in shared_lines:
            a, b = Aa.loc[m].values, Dd.loc[m].values
            ok = ~np.isnan(a) & ~np.isnan(b)
            if ok.sum() > 200:
                rel.append(np.corrcoef(a[ok], b[ok])[0, 1])
        rel = np.array(rel)
        r1 = float(np.median(rel))
        # Spearman-Brown: the analysis uses the COMBINED (2-screen) estimate,
        # whose reliability is higher than one screen's.
        r2 = 2 * r1 / (1 + r1) if r1 > 0 else np.nan
        print(f"\n  per-line profile test-retest (single screen): "
              f"median r = {r1:.3f}")
        print(f"    IQR [{np.percentile(rel, 25):.3f}, "
              f"{np.percentile(rel, 75):.3f}]   n = {len(rel)} lines")
        print(f"  Spearman-Brown to the 2-screen combined estimate: "
              f"r_xx = {r2:.3f}")
        res["test_retest"] = {
            "n_lines": int(len(rel)), "n_genes": int(len(keep)),
            "single_screen_r": r1,
            "iqr": [float(np.percentile(rel, 25)), float(np.percentile(rel, 75))],
            "combined_r_xx": float(r2)}

        print()
        print("  DISATTENUATION")
        obs = 0.0201
        for nm, r in (("single screen", r1), ("combined (used in 05)", r2)):
            if r and r > 0:
                dis = obs / r
                print(f"    r_xx = {r:.3f} ({nm:<22}) -> delta_true "
                      f"{dis:+.4f}   {dis / CRISPR_REF:.2f}x CRISPR")
        res["disattenuation"] = {
            "observed_rnai_delta": obs,
            "crispr_reference": CRISPR_REF,
            "disattenuated_single": float(obs / r1) if r1 > 0 else None,
            "disattenuated_combined": float(obs / r2) if r2 and r2 > 0 else None}

# ------------------------------------------------------------ C  (power control)
C.banner("7C -- POWER CONTROL:  CRISPR RESTRICTED TO THE RNAi LINES")
print("  Separates 'fewer lines' from 'different perturbation'. Same estimator,")
print("  same genes, same 659 lines.")
print()

W_ri = C.load_rnai(verbose=False)
sd = W_ri.std(axis=0)
W_ri = W_ri.loc[:, sd > sd.median()]
W_ri = W_ri[[g for g in W_ri.columns if g not in ce]]
D_ri = C.dependency_correlation(W_ri)

W_cr = C.load_dependency(selective_only=True, verbose=False)
W_cr = W_cr[[g for g in W_cr.columns if g not in ce]]
D_cr = C.dependency_correlation(W_cr)

rnai_lines = sorted(set(S_rna.index) & set(D_ri.index) & set(lin.index))


def delta(S, D, lines, n_boot=600, top_frac=0.01, max_pairs=300_000, seed=C.SEED):
    Sx = S.loc[lines, lines].values
    Dx = D.loc[lines, lines].values
    lv = lin.reindex(lines).lineage.values.astype(object)
    n = len(lines)

    def stat(idx, orig=None, rng=None):
        iu, ju = np.triu_indices(len(idx), k=1)
        if orig is not None:
            g = orig[iu] != orig[ju]
            iu, ju = iu[g], ju[g]
        if len(iu) > max_pairs:
            sel = rng.choice(len(iu), max_pairs, replace=False)
            iu, ju = iu[sel], ju[sel]
        ai, aj = idx[iu], idx[ju]
        s, d = Sx[ai, aj], Dx[ai, aj]
        same = lv[ai] == lv[aj]
        ok = ~np.isnan(s) & ~np.isnan(d)
        s, d, same = s[ok], d[ok], same[ok]
        if len(s) < 200 or same.sum() < 50:
            return np.nan
        k = max(10, int(len(s) * top_frac))
        part = np.argpartition(-s, k - 1)[:k]
        return float(np.median(d[part]) - np.median(d[same]))

    rng = np.random.default_rng(seed)
    pt = stat(np.arange(n), rng=rng)
    b = []
    for _ in range(n_boot):
        dr = rng.integers(0, n, n)
        v = stat(dr, orig=dr, rng=rng)
        if np.isfinite(v):
            b.append(v)
    lo, hi = (np.percentile(b, [2.5, 97.5]) if len(b) > 20 else (np.nan, np.nan))
    return pt, float(lo), float(hi), (np.array(b) <= 0).mean() if b else np.nan


print(f"  {'arm':<38} {'lines':>6} {'delta':>9} {'95% CI':>22} {'p':>8}")
cr_full = sorted(set(S_rna.index) & set(D_cr.index) & set(lin.index))
for nm, D_, lines_ in [("CRISPR, all lines", D_cr, cr_full),
                       ("CRISPR, restricted to RNAi lines", D_cr,
                        sorted(set(cr_full) & set(rnai_lines))),
                       ("RNAi", D_ri, rnai_lines)]:
    pt, lo, hi, p = delta(S_rna, D_, lines_)
    print(f"  {nm:<38} {len(lines_):>6,} {pt:>+9.4f} "
          f"{f'[{lo:+.4f}, {hi:+.4f}]':>22} {p:>8.3f}")
    res.setdefault("power_control", {})[nm] = {
        "n_lines": len(lines_), "delta": pt, "ci": [lo, hi], "boot_p": float(p)}

# ------------------------------------------------------------ B  (low-SD)
C.banner("7B -- LOW-SD SUBSET  (DEMETER2 posterior standard deviations)")
sdf = D2 / "D2_combined_gene_dep_score_SDs.csv"
if not sdf.exists():
    print("  posterior SD file not present -- skipping")
    res["low_sd"] = None
else:
    SD = load_d2_matrix("D2_combined_gene_dep_score_SDs.csv")
    common_g = [g for g in W_ri.columns if g in SD.columns]
    common_l = [m for m in W_ri.index if m in SD.index]
    SDx = SD.loc[common_l, common_g]
    print(f"  posterior SDs: {SDx.shape[0]} lines x {SDx.shape[1]:,} genes")
    print(f"  median posterior SD: {float(np.nanmedian(SDx.values)):.4f}")
    print()
    print("  If attenuation is the explanation, restricting to precisely-estimated")
    print("  genes should RAISE the delta. If it is a confound, it stays flat.")
    print()
    print(f"  {'gene subset':<34} {'genes':>7} {'delta':>9} {'95% CI':>22}")
    res["low_sd"] = {}
    gene_sd = SDx.median(axis=0).sort_values()
    for q, nm in ((1.00, "all genes"), (0.50, "lowest 50% posterior SD"),
                  (0.25, "lowest 25% posterior SD"),
                  (0.10, "lowest 10% posterior SD")):
        keep = list(gene_sd.index[:max(200, int(len(gene_sd) * q))])
        Wq = W_ri[keep]
        Dq = C.dependency_correlation(Wq)
        lines_q = sorted(set(S_rna.index) & set(Dq.index) & set(lin.index))
        pt, lo, hi, p = delta(S_rna, Dq, lines_q, n_boot=400)
        print(f"  {nm:<34} {len(keep):>7,} {pt:>+9.4f} "
              f"{f'[{lo:+.4f}, {hi:+.4f}]':>22}")
        res["low_sd"][nm] = {"n_genes": len(keep), "delta": pt, "ci": [lo, hi]}

# ------------------------------------------------------------ D
C.banner("7D -- THE SYMMETRY CHECK:  CRISPR's RELIABILITY IS NOT 1.0 EITHER")
print("  Disattenuating only the RNAi arm would rig the comparison. CRISPR has")
print("  measurement error too, and the fair contrast is disattenuated against")
print("  disattenuated. Estimated the same way: two independent CRISPR screens,")
print("  Broad (Achilles) and Sanger (Project Score), on lines screened in both.")
print()
W_b, _ = C.GA.load_chronos("achilles", sym2ensg, verbose=False)
W_s, _ = C.GA.load_chronos("score", sym2ensg, verbose=False)
sl = sorted(set(W_b.index) & set(W_s.index))
sg = sorted((set(W_b.columns) & set(W_s.columns)) - ce)
print(f"  screened in BOTH CRISPR screens: {len(sl)} lines, {len(sg):,} genes")
if len(sl) >= 20:
    Bb, Ss = W_b.loc[sl, sg], W_s.loc[sl, sg]
    sdc = pd.concat([Bb, Ss]).std(axis=0)
    kc = sdc[sdc > sdc.median()].index
    Bb, Ss = Bb[kc], Ss[kc]
    relc = []
    for m in sl:
        a, b = Bb.loc[m].values, Ss.loc[m].values
        ok = ~np.isnan(a) & ~np.isnan(b)
        if ok.sum() > 200:
            relc.append(np.corrcoef(a[ok], b[ok])[0, 1])
    relc = np.array(relc)
    rc1 = float(np.median(relc))
    rc2 = 2 * rc1 / (1 + rc1) if rc1 > 0 else np.nan
    print(f"  per-line profile test-retest (single screen): r = {rc1:.3f}   "
          f"n = {len(relc)}")
    print(f"  Spearman-Brown to a 2-screen estimate       : r_xx = {rc2:.3f}")
    res["crispr_test_retest"] = {"n_lines": int(len(relc)),
                                 "single_screen_r": rc1, "combined_r_xx": rc2}
    tr_ = res.get("test_retest")
    if tr_:
        # The main CRISPR arm is 24Q4 Avana, a SINGLE screen, so its reliability
        # is the single-screen figure. The RNAi arm is DEMETER2 combined, which
        # pools two screens, so its reliability is the Spearman-Brown one. Match
        # each arm to what it actually used.
        dis_cr = CRISPR_REF / rc1 if rc1 > 0 else np.nan
        dis_ri = 0.0201 / tr_["combined_r_xx"] if tr_["combined_r_xx"] > 0 else np.nan
        print()
        print(f"  {'arm':<28} {'r_xx':>7} {'observed':>10} {'disattenuated':>15}")
        print(f"  {'CRISPR (24Q4, 1 screen)':<28} {rc1:>7.3f} "
              f"{CRISPR_REF:>+10.4f} {dis_cr:>+15.4f}")
        print(f"  {'RNAi (DEMETER2, 2 screens)':<28} "
              f"{tr_['combined_r_xx']:>7.3f} {0.0201:>+10.4f} {dis_ri:>+15.4f}")
        print()
        print(f"  ratio of DISATTENUATED estimates, RNAi / CRISPR: "
              f"{dis_ri / dis_cr:.2f}x")
        res["symmetric_disattenuation"] = {
            "crispr_r_xx": rc1, "crispr_disattenuated": float(dis_cr),
            "rnai_r_xx": tr_["combined_r_xx"], "rnai_disattenuated": float(dis_ri),
            "ratio": float(dis_ri / dis_cr)}
else:
    print("  too few doubly-screened CRISPR lines")
    res["crispr_test_retest"] = None

# ------------------------------------------------------------ verdict
C.banner("VERDICT")
tr = res.get("test_retest")
if tr:
    r2 = tr["combined_r_xx"]
    dis = 0.0201 / r2 if r2 > 0 else np.nan
    ratio = dis / CRISPR_REF
    print(f"  RNAi reliability (2-screen, Spearman-Brown) : r_xx = {r2:.3f}")
    print(f"  observed RNAi delta                         : +0.0201")
    print(f"  disattenuated                               : {dis:+.4f}")
    print(f"  vs CRISPR {CRISPR_REF:+.4f}                        : {ratio:.2f}x")
    print()
    sym = res.get("symmetric_disattenuation")
    if sym:
        ratio = sym["ratio"]
        print(f"  (using the SYMMETRIC comparison -- both arms disattenuated)")
    if ratio >= 0.75:
        v = "gap_is_measurement_error__component_passes"
        print("  The gap between the CRISPR and RNAi arms is accounted for by")
        print("  RNAi's measurement error. There is no residual confound to")
        print("  explain. This is a FULL pass that was hiding behind attenuation,")
        print("  and the cautious version of the write-up is the wrong one.")
    elif ratio >= 0.4:
        v = "gap_is_part_noise_part_confound"
        print("  Attenuation explains a substantial part of the gap but not all")
        print("  of it. Report the disattenuated estimate WITH its assumption")
        print("  stated; some CRISPR-specific inflation remains.")
    else:
        v = "gap_survives_disattenuation__confound_remains"
        print("  The gap survives correction for RNAi's reliability. Most of the")
        print("  CRISPR advantage is CRISPR-specific and the cautious reading")
        print("  stands.")
    print(f"\n  VERDICT: {v}")
    res["verdict"] = v

out = C.OUT / "07_reliability_results.json"
out.write_text(json.dumps(res, indent=2, default=float), encoding="utf-8")
print(f"\nwrote {out}")
