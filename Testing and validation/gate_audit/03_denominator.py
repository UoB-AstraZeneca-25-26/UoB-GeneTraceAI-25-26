"""
gate_audit/03_denominator.py
----------------------------
AUDIT 3 of 6.  THE PROFILING CONFOUND, in its direct form.

This audit was not on the original list. It was found while building 01, and it
supersedes the protein-detection version of the profiling question, because it
acts on the label rather than on a covariate.

THE FINDING
-----------
`src/pipeline/test_run_abundance_gate.py` builds each gene's line set from the
EXPRESSION and PROTEOMICS indices:

    idx = sorted(set(zE.index) | set(zP.index))
    d   = pd.DataFrame({...}, index=idx)
    d["pos"] = d.index.isin(pos_by_gene[g])          # <-- here

`pos_by_gene[g]` contains only lines that were actually CRISPR-screened. Any
line with expression but no screen therefore receives `pos = False` — it is
counted as a confirmed non-dependency rather than excluded as unmeasured.

    expression lines                1,479
    Chronos (Achilles) lines          767
    overlap                           754
    expression lines never screened   725   (49% of the denominator)

So roughly half of every gene's "negatives" are lines where nobody looked.

WHY THAT MANUFACTURES THE GATE
------------------------------
Guaranteed-negative lines are not spread evenly across the score. They
concentrate in the bottom decile — the unscreened share of decile 1 runs about
7 percentage points above the unscreened share overall. Decile 1 is therefore
enriched for lines that CANNOT be positive, and its positive rate falls for a
reason that has nothing to do with abundance.

That is the same confound raised about protein profiling, but one step closer to
the outcome: it is not that profiled lines differ biologically, it is that
unprofiled lines were scored as negatives.

WHAT THIS SCRIPT MEASURES
-------------------------
  A. The two line sets, and the unscreened share by decile.
  B. The decile-1 ratio and the gate-region pAUC computed both ways, on the same
     genes and the same scores. `original` reproduces the published construction;
     `screened_only` is the corrected one.
  C. A placebo: relabel every unscreened line as negative but REMOVE the real
     labels, leaving only the screening indicator. Any depletion that survives
     is manufactured entirely by the denominator, with no dependency information
     in the data at all.

Read-only. Writes only into gate_audit/outputs/.

Run:
    python gate_audit/03_denominator.py
    python gate_audit/03_denominator.py --n-genes 3000
"""
import argparse
import json

import numpy as np
import pandas as pd
from scipy.stats import wilcoxon

import common as C

ap = argparse.ArgumentParser()
ap.add_argument("--n-genes", type=int, default=2000,
                help="genes to evaluate; this audit re-loads expression for the "
                     "UNSCREENED lines, which the panel deliberately excludes")
args = ap.parse_args()

res = {}
C.banner("AUDIT 3 -- IS THE DENOMINATOR HALF UNMEASURED LINES?")

uni, valid_ensg, symbol_of, sym2ensg, u2e = C.load_gene_lookup()
W, chk = C.load_chronos("achilles", sym2ensg)
panel = C.read_panel()
G = pd.read_parquet(C.OUT / "genes.parquet")
valid_genes = set(G.loc[G.is_valid, "ensg_id"])

genes = sorted(panel.ensg_id.unique())[: args.n_genes]
E = C.load_expression(genes)
genes = [g for g in genes if g in E.columns]
screened = set(W.index)

# ------------------------------------------------------------ A
C.banner("3A -- THE TWO LINE SETS")
print(f"  expression lines                : {E.shape[0]:,}")
print(f"  Chronos (Achilles) lines        : {W.shape[0]:,}")
print(f"  overlap                         : {len(set(E.index) & screened):,}")
print(f"  expression lines NEVER screened : {len(set(E.index) - screened):,} "
      f"({100*len(set(E.index)-screened)/E.shape[0]:.1f}% of the denominator)")
res["line_sets"] = {
    "expression_lines": int(E.shape[0]),
    "screened_lines": int(W.shape[0]),
    "overlap": int(len(set(E.index) & screened)),
    "never_screened": int(len(set(E.index) - screened)),
    "never_screened_frac": float(len(set(E.index) - screened) / E.shape[0]),
}

# ------------------------------------------------------------ B
C.banner("3B -- THE GATE, BOTH DENOMINATORS, SAME GENES AND SAME SCORES")
rows, unsc_by_dec = [], []
for g in genes:
    e = E[g].dropna()
    if len(e) < C.MIN_LINES:
        continue
    z = C.vdw(e)
    posset = set(W.index[(W[g] <= C.DEP_THRESHOLD).values])
    rec = {"ensg_id": g, "is_valid": g in valid_genes}

    for mode in ("original", "screened_only"):
        idx = list(z.index) if mode == "original" else [m for m in z.index
                                                        if m in screened]
        if len(idx) < C.MIN_LINES:
            continue
        s = z.reindex(idx).values
        p = np.fromiter((m in posset for m in idx), bool, len(idx))
        if p.sum() < C.MIN_POS or (~p).sum() < C.MIN_NEG:
            continue
        base = p.mean()
        q, keep = C.decile_index(s, 10, "ordinal")
        qv, pv = q[keep], p[keep]
        rec[f"{mode}_base"] = float(base)
        rec[f"{mode}_d1"] = float(pv[qv == 0].mean() / base)
        rec[f"{mode}_pauc"] = C.partial_auc(s, p, *C.GATE_REGION)
        rec[f"{mode}_n"] = int(len(idx))

    # unscreened share by decile, on the ORIGINAL construction
    s = z.values
    q, keep = C.decile_index(s, 10, "ordinal")
    unsc = np.fromiter((m not in screened for m in z.index), bool, len(z))
    unsc_by_dec.append([unsc[keep][q[keep] == i].mean() if (q[keep] == i).sum()
                        else np.nan for i in range(10)])
    rows.append(rec)

D = pd.DataFrame(rows)
D.to_parquet(C.OUT / "03_denominator_per_gene.parquet", index=False)
U = np.vstack(unsc_by_dec)

print("  Share of each decile that was NEVER CRISPR-screened.")
print("  A flat row would mean the denominator is neutral; a falling row means")
print("  decile 1 is packed with lines that cannot be positive.")
print()
print(f"  {'decile':>7} {'unscreened share':>18}")
for i in range(10):
    v = np.nanmean(U[:, i])
    print(f"  {i+1:>7} {v:>18.4f}   {'#' * int(round(v * 60))}")
overall = float(np.nanmean(U))
d1_unsc = float(np.nanmean(U[:, 0]))
print()
print(f"  overall unscreened share : {overall:.4f}")
print(f"  decile-1 unscreened share: {d1_unsc:.4f}   "
      f"({100*(d1_unsc-overall):+.1f} percentage points)")
res["unscreened_share_by_decile"] = [float(np.nanmean(U[:, i])) for i in range(10)]
res["unscreened_share_overall"] = overall
res["unscreened_enrichment_decile1_pp"] = float(100 * (d1_unsc - overall))

print()
print(f"  {'quantity':<32} {'ORIGINAL':>14} {'SCREENED ONLY':>16} {'delta':>10}")
res["comparison"] = {}
for col, lab in [("d1", "decile-1 ratio"), ("pauc", "gate-region pAUC"),
                 ("base", "base rate"), ("n", "lines per gene")]:
    a = D[f"original_{col}"].dropna()
    b = D[f"screened_only_{col}"].dropna()
    if not len(a) or not len(b):
        continue
    print(f"  {lab:<32} {a.median():>14.4f} {b.median():>16.4f} "
          f"{b.median()-a.median():>+10.4f}")
    res["comparison"][col] = {"original": float(a.median()),
                              "screened_only": float(b.median()),
                              "delta": float(b.median() - a.median())}
both = D.dropna(subset=["original_d1", "screened_only_d1"])
if len(both) > 5:
    w = wilcoxon(both.screened_only_d1 - both.original_d1)
    print(f"\n  paired shift in decile-1 ratio: "
          f"{(both.screened_only_d1 - both.original_d1).median():+.4f}   "
          f"p {w.pvalue:.3g}   ({len(both)} genes)")
    res["comparison"]["paired_d1_shift_p"] = float(w.pvalue)
    res["comparison"]["n_paired"] = int(len(both))

# ------------------------------------------------------------ C
C.banner("3C -- PLACEBO:  no dependency information at all")
print("  Keep the construction, discard the biology. Every screened line is")
print("  relabelled NEGATIVE at random with the gene's own prevalence, so the")
print("  only real structure left is 'unscreened lines are always negative'.")
print("  Depletion that survives here is manufactured by the denominator.")
print()
rng = np.random.default_rng(C.SEED)
plac = []
for g in genes[:1000]:
    e = E[g].dropna()
    if len(e) < C.MIN_LINES:
        continue
    z = C.vdw(e)
    n_pos = int((W[g] <= C.DEP_THRESHOLD).sum())
    if n_pos < C.MIN_POS:
        continue
    scr = [m for m in z.index if m in screened]
    if len(scr) < C.MIN_LINES:
        continue
    fake_pos = set(rng.choice(scr, size=min(n_pos, len(scr)), replace=False))
    idx = list(z.index)
    s = z.reindex(idx).values
    p = np.fromiter((m in fake_pos for m in idx), bool, len(idx))
    if p.sum() < C.MIN_POS:
        continue
    q, keep = C.decile_index(s, 10, "ordinal")
    plac.append(p[keep][q[keep] == 0].mean() / p.mean())
plac = pd.Series(plac).dropna()
print(f"  placebo decile-1 ratio (median over {len(plac)} genes): {plac.median():.4f}")
print(f"  depleted in {100*(plac < 1).mean():.1f}% of genes   "
      f"p {wilcoxon(plac - 1.0).pvalue:.3g}")
res["placebo"] = {"n_genes": int(len(plac)),
                  "decile1_ratio": float(plac.median()),
                  "frac_depleted": float((plac < 1).mean()),
                  "wilcoxon_p": float(wilcoxon(plac - 1.0).pvalue)}

# ------------------------------------------------------------ verdict
C.banner("VERDICT")
orig = res["comparison"]["d1"]["original"]
fix = res["comparison"]["d1"]["screened_only"]
pl = res["placebo"]["decile1_ratio"]
recovered = (1 - pl) / (1 - orig) if orig < 1 else np.nan
print(f"  published construction      : decile-1 ratio {orig:.4f}  "
      f"({100*(1-orig):.1f}% depletion)")
print(f"  screened lines only         : decile-1 ratio {fix:.4f}  "
      f"({100*(1-fix):.1f}% depletion)")
print(f"  placebo, no real labels     : decile-1 ratio {pl:.4f}  "
      f"({100*(1-pl):.1f}% depletion)")
print(f"  share of the published depletion reproduced with NO biology: "
      f"{100*recovered:.0f}%")
print()
if fix < 0.90:
    verdict = "denominator_inflated_the_effect_but_a_real_gate_remains"
    print("  Correcting the denominator shrinks the effect but leaves a")
    print("  depletion above the project's 10% effect floor. The gate is real and")
    print("  was overstated.")
elif fix < 0.99 and res["comparison"].get("paired_d1_shift_p", 1) < 0.05:
    verdict = "effect_is_mostly_denominator__residual_gate_is_below_the_effect_floor"
    print("  Correcting the denominator removes most of the effect. What remains")
    print("  is statistically detectable but sits BELOW the project's own 10%")
    print("  effect floor, so it does not support the headline as written.")
    print("  The decile-1 depletion must be requoted on screened lines only.")
else:
    verdict = "effect_does_not_survive_the_denominator_correction"
    print("  The effect does not survive. It was an artefact of scoring")
    print("  never-screened cell lines as confirmed non-dependencies.")
print(f"\n  VERDICT: {verdict}")
res["verdict"] = verdict
res["share_of_depletion_without_biology"] = float(recovered)

out = C.OUT / "03_denominator_results.json"
out.write_text(json.dumps(res, indent=2, default=float), encoding="utf-8")
print(f"\nwrote {out}")
print(f"wrote {C.OUT / '03_denominator_per_gene.parquet'}")
