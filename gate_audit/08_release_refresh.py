"""
gate_audit/08_release_refresh.py
--------------------------------
AUDIT 8.  Does the corrected result hold on a current DepMap release?

WHY
---
Audits 1-7 ran on `data/DepMap_Chronos/`, which is the Zenodo 14067047 benchmark
bundle built on DepMap 20Q2: 767 Achilles lines, 754 after intersecting with the
20Q2-vintage expression matrix. That is a five-year-old snapshot and it caps the
power of every band in audit 6 and 7.

This rebuilds on DepMap 24Q4 -- CRISPRGeneEffect and
OmicsExpressionProteinCodingGenesTPMLogp1 from the SAME release, so no model-ID
drift between vintages is introduced. Both are keyed by ModelID directly, so
there is no profile-to-model collapse step and no chance of the profile join
losing lines.

WHAT IT TESTS
-------------
  A. How many lines the refresh actually buys, after the expression intersection.
  B. Whether the audit-3 denominator conclusion holds at the larger n: the gate
     computed the published way (unscreened lines default to negative) versus on
     screened lines only.
  C. Whether the audit-2 validity conclusion holds.
  D. The banded gate at the larger n, which is the number audits 6 and 7 were
     short of power for.

A refresh is a POWER gain, not a fix. If the denominator artefact reproduces at
1,100 lines, it was never a small-sample effect.

Read-only. Writes only into gate_audit/outputs/.

Run:
    python gate_audit/08_release_refresh.py
    python gate_audit/08_release_refresh.py --n-genes 4000
"""
import argparse
import json

import numpy as np
import pandas as pd
from scipy.stats import wilcoxon

import common as C

ap = argparse.ArgumentParser()
ap.add_argument("--n-genes", type=int, default=4000)
args = ap.parse_args()

BANDS = [0, .02, .05, .10, .25, 1.01]
BAND_LABELS = ["<2%", "2-5%", "5-10%", "10-25%", ">25%"]
res = {}

C.banner("AUDIT 8 -- DEPMAP 24Q4 REFRESH")
uni, valid_ensg, symbol_of, sym2ensg, u2e = C.load_gene_lookup()
W, chk = C.load_crispr_24q4(sym2ensg)
E = C.load_expression_24q4(sym2ensg)
res["scale_check"] = chk

# ------------------------------------------------------------ A
C.banner("8A -- WHAT THE REFRESH BUYS")
old_W, old_chk = C.load_chronos("achilles", sym2ensg, verbose=False)
old_E = C.load_expression(list(old_W.columns)[:50])
old_n = len(set(old_W.index) & set(old_E.index))
new_n = len(set(W.index) & set(E.index))
print(f"  20Q2 bundle : {old_W.shape[0]:>5} screened, {old_E.shape[0]:>5} expression,"
      f" {old_n:>5} usable")
print(f"  24Q4 release: {W.shape[0]:>5} screened, {E.shape[0]:>5} expression,"
      f" {new_n:>5} usable")
print(f"  gain        : {new_n - old_n:+} lines  ({new_n/old_n:.2f}x)")
print()
print(f"  lines still WITHOUT a screen in 24Q4: "
      f"{len(set(E.index) - set(W.index)):,} of {E.shape[0]:,} expression lines")
print(f"  -> the denominator hazard from audit 3 does NOT go away with the")
print(f"     refresh; it shrinks. It still has to be handled explicitly.")
res["line_counts"] = {
    "old_screened": int(old_W.shape[0]), "old_usable": int(old_n),
    "new_screened": int(W.shape[0]), "new_expression": int(E.shape[0]),
    "new_usable": int(new_n),
    "new_expression_without_screen": int(len(set(E.index) - set(W.index))),
    "gain_factor": float(new_n / old_n)}

# ------------------------------------------------------------ build
C.banner("8B -- SCORING ON 24Q4")
genes = sorted(set(W.columns) & set(E.columns) & valid_ensg)
rng = np.random.default_rng(C.SEED)
if args.n_genes < len(genes):
    genes = sorted(rng.choice(genes, args.n_genes, replace=False))
print(f"  genes: {len(genes):,}")

fe = (E[genes] > C.EXPRESSED_MIN).sum(axis=0) / E[genes].notna().sum(axis=0)
is_valid = fe >= C.SILENT_FRAC
print(f"  valid (frac_expressed >= {C.SILENT_FRAC}): {int(is_valid.sum()):,} "
      f"({100*is_valid.mean():.1f}%)   invalid {int((~is_valid).sum()):,} "
      f"({100*(~is_valid).mean():.1f}%)")
res["validity"] = {"n_valid": int(is_valid.sum()),
                   "n_invalid": int((~is_valid).sum()),
                   "invalid_frac": float((~is_valid).mean())}

screened = set(W.index)
rows = []
for i, g in enumerate(genes):
    if i and i % 1000 == 0:
        print(f"    {i:,}/{len(genes):,}")
    e = E[g].dropna()
    if len(e) < C.MIN_LINES:
        continue
    z = C.vdw(e)
    col = W[g]
    pos = set(col.index[(col <= C.DEP_THRESHOLD).fillna(False)])
    meas = set(col.index[col.notna()])
    rec = {"ensg_id": g, "is_valid": bool(is_valid[g])}
    for mode in ("original", "screened_only"):
        idx = list(z.index) if mode == "original" else [m for m in z.index
                                                        if m in meas]
        if len(idx) < C.MIN_LINES:
            continue
        s = z.reindex(idx).values
        p = np.fromiter((m in pos for m in idx), bool, len(idx))
        if p.sum() < C.MIN_POS or (~p).sum() < C.MIN_NEG:
            continue
        base = p.mean()
        q, keep = C.decile_index(s, 10, "ordinal")
        qv, pv = q[keep], p[keep]
        rec[f"{mode}_base"] = float(base)
        rec[f"{mode}_d1"] = float(pv[qv == 0].mean() / base)
        rec[f"{mode}_pauc"] = C.partial_auc(s, p, *C.GATE_REGION)
        rec[f"{mode}_n"] = int(len(idx))
    rows.append(rec)

R = pd.DataFrame(rows)
R.to_parquet(C.OUT / "08_release_refresh_per_gene.parquet", index=False)
print(f"  scored genes: {len(R):,}")

# ------------------------------------------------------------ C
C.banner("8C -- DOES THE DENOMINATOR ARTEFACT REPRODUCE AT 24Q4?")
print(f"  {'quantity':<28} {'ORIGINAL':>13} {'SCREENED ONLY':>16} {'delta':>10}")
res["denominator"] = {}
for col, lab in [("d1", "decile-1 ratio"), ("pauc", "gate-region pAUC"),
                 ("base", "base rate"), ("n", "lines per gene")]:
    a = R[f"original_{col}"].dropna()
    b = R[f"screened_only_{col}"].dropna()
    if not len(a) or not len(b):
        continue
    print(f"  {lab:<28} {a.median():>13.4f} {b.median():>16.4f} "
          f"{b.median()-a.median():>+10.4f}")
    res["denominator"][col] = {"original": float(a.median()),
                               "screened_only": float(b.median()),
                               "delta": float(b.median() - a.median())}
both = R.dropna(subset=["original_d1", "screened_only_d1"])
if len(both) > 5:
    p = wilcoxon(both.screened_only_d1 - both.original_d1).pvalue
    print(f"\n  paired shift: "
          f"{(both.screened_only_d1 - both.original_d1).median():+.4f}   "
          f"p {p:.3g}   ({len(both):,} genes)")
    res["denominator"]["paired_p"] = float(p)

# ------------------------------------------------------------ D
C.banner("8D -- VALIDITY AT 24Q4  (audit 2, re-tested)")
print(f"  {'stratum':<10} {'genes':>7} {'d1 ratio':>10} {'pAUC gate':>11}")
res["validity_split"] = {}
for nm, sub in [("all", R), ("valid", R[R.is_valid]), ("invalid", R[~R.is_valid])]:
    s = sub.dropna(subset=["screened_only_d1"])
    if len(s) < 6:
        continue
    print(f"  {nm:<10} {len(s):>7,} {s.screened_only_d1.median():>10.4f} "
          f"{s.screened_only_pauc.median():>11.4f}")
    res["validity_split"][nm] = {"n": int(len(s)),
                                 "d1_ratio": float(s.screened_only_d1.median()),
                                 "pauc": float(s.screened_only_pauc.median())}

# ------------------------------------------------------------ E
C.banner("8E -- THE BANDED GATE AT 24Q4  (the number audits 6-7 lacked power for)")
V = R[R.is_valid].dropna(subset=["screened_only_d1", "screened_only_base"]).copy()
V["band"] = pd.cut(V.screened_only_base, BANDS, labels=BAND_LABELS)
print(f"  {'band':>8} {'genes':>7} {'d1 ratio':>10} {'p':>11} {'pAUC gate':>11} "
      f"{'p':>11}")
res["banded"] = {}
for b in BAND_LABELS:
    s = V[V.band == b]
    if len(s) < 6:
        print(f"  {b:>8} {len(s):>7}   too few")
        continue
    p1 = wilcoxon(s.screened_only_d1 - 1.0).pvalue
    pa = s.screened_only_pauc.dropna()
    p2 = wilcoxon(pa - 0.5).pvalue if len(pa) > 5 else np.nan
    print(f"  {b:>8} {len(s):>7,} {s.screened_only_d1.median():>10.4f} "
          f"{p1:>11.3g} {pa.median():>11.4f} {p2:>11.3g}")
    res["banded"][b] = {"n": int(len(s)),
                        "d1_ratio": float(s.screened_only_d1.median()),
                        "d1_p": float(p1), "pauc": float(pa.median()),
                        "pauc_p": float(p2),
                        "frac_depleted": float((s.screened_only_d1 < 1).mean())}

# ------------------------------------------------------------ verdict
C.banner("VERDICT")
o = res["denominator"]["d1"]["original"]
f = res["denominator"]["d1"]["screened_only"]
print(f"  denominator artefact at 24Q4: original {o:.4f} -> screened {f:.4f} "
      f"({f-o:+.4f})")
b510 = res["banded"].get("5-10%", {})
if b510:
    print(f"  5-10% band at n={res['line_counts']['new_usable']}: "
          f"d1 {b510['d1_ratio']:.4f} (p {b510['d1_p']:.3g}), "
          f"pAUC {b510['pauc']:.4f} (p {b510['pauc_p']:.3g}), "
          f"{b510['n']} genes")
reproduces = (f - o) > 0.05
verdict = ("denominator_artefact_reproduces_at_larger_n" if reproduces
           else "denominator_artefact_smaller_at_24Q4")
print(f"\n  VERDICT: {verdict}")
res["verdict"] = verdict

out = C.OUT / "08_release_refresh_results.json"
out.write_text(json.dumps(res, indent=2, default=float), encoding="utf-8")
print(f"\nwrote {out}")
