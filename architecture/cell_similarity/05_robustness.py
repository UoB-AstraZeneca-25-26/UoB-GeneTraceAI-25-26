"""
cell_similarity/05_robustness.py
--------------------------------
Three checks on the +0.056, run together because they are the same estimator
with different inputs and the interactions matter.

    1. CLUSTER BOOTSTRAP OVER LINES.  Every quantity here gets an interval.
    2. RNAi INSTEAD OF CRISPR.        Does the RNA advantage survive a
                                      perturbation that does not share CRISPR's
                                      expression confound?
    3. COMMON ESSENTIALS REMOVED.     What does the graph capture once the
                                      housekeeping floor is taken out?

THE ESTIMATOR
-------------
    graph - tissue  =  median dep-corr of the top 1% most similar pairs
                     - median dep-corr of same-lineage pairs

The subtraction is the point: same-lineage is the free baseline you get from a
tissue lookup, so this is what building a graph buys over not building one.

WHY THE BOOTSTRAP RESAMPLES LINES, NOT PAIRS
---------------------------------------------
Each cell line appears in ~1,600 pairs. Pair-level statistics therefore treat
one line's idiosyncrasies as ~1,600 independent observations, and the nominal
precision is fiction. This is the same fault class as gaps 7 and 17 in the main
pipeline, and the same one that inflated the gate's p-value by a design factor
of 2.4 in `gate_audit/04`.

So the resampling unit is the CELL LINE. Lines are drawn with replacement, the
pair set is rebuilt from the drawn lines, and the statistic is recomputed.
Pairs formed between two draws of the SAME original line are dropped — they
would have similarity 1.0 and dependency correlation 1.0 by construction and
would bias every replicate upward.

WHY RNAi IS THE DECIDING TEST
------------------------------
The graph is built from RNA. It is being validated against CRISPR dependency
correlation. But `gate_audit` established that CRISPR dependency is confounded
with expression — a gene that is not expressed scores as non-essential
regardless of function, and audit 09 found the abundance-dependency effect is
CRISPR-specific (RNAi ran OPPOSITE, Spearman +0.08 between the two).

So two lines with similar expression could show similar CRISPR profiles partly
through that artefact rather than through shared biology. RNAi does not share
the confound. If RNA still beats tissue on RNAi profiles, the similarity result
is real. If the advantage disappears, it is the expression-CRISPR artefact
found for a third time.

Read-only. Writes only into cell_similarity/outputs/.

Run:
    python cell_similarity/05_robustness.py
    python cell_similarity/05_robustness.py --n-boot 1000
"""
import argparse
import json

import numpy as np
import pandas as pd

import common as C

ap = argparse.ArgumentParser()
ap.add_argument("--n-boot", type=int, default=400)
ap.add_argument("--top-frac", type=float, default=0.01)
ap.add_argument("--max-pairs", type=int, default=300_000)
args = ap.parse_args()

res = {"params": {"n_boot": args.n_boot, "top_frac": args.top_frac}}
rng_global = np.random.default_rng(C.SEED)

C.banner("ROBUSTNESS -- CLUSTER BOOTSTRAP, RNAi, AND COMMON ESSENTIALS")

_, _, _, sym2ensg, _ = C.GA.load_gene_lookup()
lin = C.load_lineage()
ce = C.common_essentials(sym2ensg, "inferred")
print(f"\n  DepMap inferred common essentials: {len(ce):,} genes (as ensg)")

# ------------------------------------------------------------ variants
print("\n  building dependency-profile variants")
variants = {}

W_cr = C.load_dependency(selective_only=True, verbose=True)
variants["crispr_selective"] = W_cr
keep = [g for g in W_cr.columns if g not in ce]
variants["crispr_no_common_ess"] = W_cr[keep]
print(f"    crispr_selective        : {W_cr.shape[1]:,} genes")
print(f"    crispr_no_common_ess    : {len(keep):,} genes "
      f"({W_cr.shape[1]-len(keep):,} common essentials removed)")

W_ri_raw = C.load_rnai(verbose=True)
sd = W_ri_raw.std(axis=0)
sel = sd > sd.median()          # RNAi is compressed; use a relative cut
W_ri = W_ri_raw.loc[:, sel]
variants["rnai_selective"] = W_ri
keep_r = [g for g in W_ri.columns if g not in ce]
variants["rnai_no_common_ess"] = W_ri[keep_r]
print(f"    rnai_selective          : {W_ri.shape[1]:,} genes "
      f"(sd above median; RNAi effect sizes are compressed so an absolute "
      f"cut would not transfer)")
print(f"    rnai_no_common_ess      : {len(keep_r):,} genes")

DEPS = {k: C.dependency_correlation(v) for k, v in variants.items()}
res["variants"] = {k: {"n_genes": int(v.shape[1]), "n_lines": int(v.shape[0])}
                   for k, v in variants.items()}


# ------------------------------------------------------------ estimator
def stat_on(S_vals, D_vals, lv, idx, top_frac, max_pairs, rng, orig=None):
    """(random, same-lineage, graph top-frac, graph - tissue) on the given line
    indices. `orig` marks which ORIGINAL line each drawn index came from, so
    pairs between two draws of the same line can be dropped."""
    n = len(idx)
    iu, ju = np.triu_indices(n, k=1)
    if orig is not None:
        good = orig[iu] != orig[ju]
        iu, ju = iu[good], ju[good]
    if max_pairs and len(iu) > max_pairs:
        sel = rng.choice(len(iu), max_pairs, replace=False)
        iu, ju = iu[sel], ju[sel]
    ai, aj = idx[iu], idx[ju]
    s = S_vals[ai, aj]
    d = D_vals[ai, aj]
    same = lv[ai] == lv[aj]
    ok = ~np.isnan(s) & ~np.isnan(d)
    s, d, same = s[ok], d[ok], same[ok]
    if len(s) < 200 or same.sum() < 50:
        return (np.nan,) * 4
    k = max(10, int(len(s) * top_frac))
    part = np.argpartition(-s, k - 1)[:k]
    rnd = float(np.median(d))
    tis = float(np.median(d[same]))
    gph = float(np.median(d[part]))
    return rnd, tis, gph, gph - tis


# ------------------------------------------------------------ run
C.banner("5A -- POINT ESTIMATES AND CLUSTER-BOOTSTRAP INTERVALS")
print(f"  resampling LINES with replacement, {args.n_boot} replicates")
print(f"  top {args.top_frac:.0%} of pairs by similarity = 'graph'")
print()
print(f"  {'axis':<7} {'dependency':<22} {'lines':>6} {'random':>8} "
      f"{'tissue':>8} {'graph':>8} {'graph-tissue':>13} {'95% CI':>20}")

res["results"] = {}
for axis in C.available_axes():
    S = C.load_axis_similarity(axis)
    for vname, D in DEPS.items():
        lines = sorted(set(S.index) & set(D.index) & set(lin.index))
        if len(lines) < 120:
            print(f"  {axis:<7} {vname:<22} {len(lines):>6}   too few shared lines")
            continue
        Sx = S.loc[lines, lines].values
        Dx = D.loc[lines, lines].values
        lv = lin.reindex(lines).lineage.values.astype(object)
        n = len(lines)
        base_idx = np.arange(n)

        rnd, tis, gph, delta = stat_on(Sx, Dx, lv, base_idx, args.top_frac,
                                       args.max_pairs, rng_global)
        boots = []
        rng = np.random.default_rng(C.SEED)
        for _ in range(args.n_boot):
            draw = rng.integers(0, n, n)
            v = stat_on(Sx, Dx, lv, draw, args.top_frac, args.max_pairs, rng,
                        orig=draw)
            if np.isfinite(v[3]):
                boots.append(v[3])
        if len(boots) < 20:
            lo = hi = np.nan
        else:
            lo, hi = np.percentile(boots, [2.5, 97.5])
        sig = "" if (np.isnan(lo) or (lo <= 0 <= hi)) else "  *"
        print(f"  {axis:<7} {vname:<22} {n:>6,} {rnd:>8.4f} {tis:>8.4f} "
              f"{gph:>8.4f} {delta:>+13.4f} "
              f"{f'[{lo:+.4f}, {hi:+.4f}]':>20}{sig}")
        res["results"].setdefault(axis, {})[vname] = {
            "n_lines": int(n), "random": rnd, "tissue": tis, "graph": gph,
            "graph_minus_tissue": delta, "ci": [float(lo), float(hi)],
            "n_boot_ok": len(boots),
            "excludes_zero": bool(not np.isnan(lo) and not (lo <= 0 <= hi))}
print()
print("  * = 95% CI excludes zero")
print()
print("  Compare against the uncorrected numbers this replaces: RNA +0.056,")
print("  metabolomics -0.020, miRNA -0.015, all with no intervals at all.")

# ------------------------------------------------------------ 5B
C.banner("5B -- WHAT THE HOUSEKEEPING FLOOR WAS HIDING")
print("  The random-pair baseline is high because every line depends on")
print("  ribosome, proteasome and polymerase. Two unrelated lines look similar")
print("  because they are both alive.")
print()
print(f"  {'dependency':<22} {'random-pair floor':>19} {'graph top 1%':>14} "
      f"{'headroom used':>15}")
res["floor"] = {}
for vname in DEPS:
    r = res["results"].get("rna", {}).get(vname)
    if not r:
        continue
    head = (r["graph"] - r["random"]) / (1.0 - r["random"])
    print(f"  {vname:<22} {r['random']:>19.4f} {r['graph']:>14.4f} "
          f"{head:>14.1%}")
    res["floor"][vname] = {"random_floor": r["random"], "graph": r["graph"],
                           "headroom_used": float(head)}
print()
print("  'headroom used' = how far the graph moves from the floor toward 1.0.")
print("  It is the effect size that the shared-housekeeping baseline was")
print("  compressing.")

# ------------------------------------------------------------ 5C
C.banner("5C -- THE DECIDING TEST:  DOES RNA BEAT TISSUE ON RNAi PROFILES?")
r_cr = res["results"].get("rna", {}).get("crispr_no_common_ess") \
    or res["results"].get("rna", {}).get("crispr_selective")
r_ri = res["results"].get("rna", {}).get("rnai_no_common_ess") \
    or res["results"].get("rna", {}).get("rnai_selective")
if r_cr and r_ri:
    print(f"  CRISPR  graph - tissue {r_cr['graph_minus_tissue']:+.4f}   "
          f"95% CI [{r_cr['ci'][0]:+.4f}, {r_cr['ci'][1]:+.4f}]   "
          f"n = {r_cr['n_lines']:,}")
    print(f"  RNAi    graph - tissue {r_ri['graph_minus_tissue']:+.4f}   "
          f"95% CI [{r_ri['ci'][0]:+.4f}, {r_ri['ci'][1]:+.4f}]   "
          f"n = {r_ri['n_lines']:,}")
    print()
    ratio = (r_ri["graph_minus_tissue"] / r_cr["graph_minus_tissue"]
             if r_cr["graph_minus_tissue"] else np.nan)
    print(f"  RNAi advantage is {ratio:.2f}x the CRISPR one")
    print()
    # Excluding zero is NOT sufficient. An advantage that survives in sign but
    # collapses in magnitude means most of what CRISPR measured is
    # CRISPR-specific -- which is the confound, just not all of it.
    if r_ri["excludes_zero"] and r_ri["graph_minus_tissue"] > 0 and ratio >= 0.5:
        verdict = "rna_advantage_survives_rnai_at_similar_magnitude"
        print("  The RNA advantage survives on a perturbation that does not share")
        print("  CRISPR's expression confound, at comparable magnitude. The")
        print("  similarity result is real. This component can ship.")
    elif r_ri["excludes_zero"] and r_ri["graph_minus_tissue"] > 0:
        verdict = "rna_advantage_survives_in_sign_but_collapses_in_magnitude"
        print("  The advantage survives in SIGN but collapses in MAGNITUDE.")
        print("  Both things are true and both must be reported:")
        print("    - it is not purely the expression-CRISPR artefact, because it")
        print("      is still positive on a perturbation that lacks that confound;")
        print(f"    - but only {ratio:.0%} of the CRISPR-measured advantage")
        print("      replicates, so the CRISPR figure is substantially inflated")
        print("      by something CRISPR-specific.")
        print("  Quote the RNAi number as the defensible estimate of what the")
        print("  graph adds, and the CRISPR one only alongside it.")
    elif r_ri["graph_minus_tissue"] > 0:
        verdict = "rna_advantage_positive_on_rnai_but_interval_includes_zero"
        print("  Directionally positive on RNAi but the interval includes zero.")
        print("  Not the artefact, but not established either -- RNAi has fewer")
        print("  lines and a noisier readout, so this is underpowered rather")
        print("  than negative. Report both intervals; do not claim the RNAi arm.")
    else:
        verdict = "rna_advantage_vanishes_on_rnai__expression_confound_again"
        print("  The advantage DISAPPEARS on RNAi. The RNA graph was beating the")
        print("  tissue baseline by exploiting the same expression-CRISPR")
        print("  confound that killed the gate. Found for the third time, and")
        print("  that is now the thesis rather than a caveat.")
    res["deciding_test"] = {"crispr": r_cr, "rnai": r_ri, "verdict": verdict}
    print(f"\n  VERDICT: {verdict}")

# ------------------------------------------------------------ 5D
C.banner("5D -- ARE THE TWO DEMOTIONS ACTUALLY NEGATIVE?")
print("  -0.020 and -0.015 were used to demote metabolomics and miRNA from")
print("  rankable. With intervals, is either distinguishable from zero?")
print()
print(f"  {'axis':<8} {'dependency':<22} {'delta':>9} {'95% CI':>22} {'call':>12}")
res["demotions"] = {}
for axis in ("metab", "mirna"):
    for vname in DEPS:
        r = res["results"].get(axis, {}).get(vname)
        if not r:
            continue
        lo, hi = r["ci"]
        call = ("negative" if hi < 0 else
                "positive" if lo > 0 else "indistinguishable")
        print(f"  {axis:<8} {vname:<22} {r['graph_minus_tissue']:>+9.4f} "
              f"{f'[{lo:+.4f}, {hi:+.4f}]':>22} {call:>12}")
        res["demotions"].setdefault(axis, {})[vname] = call

out = C.OUT / "05_robustness_results.json"
out.write_text(json.dumps(res, indent=2, default=float), encoding="utf-8")
print(f"\nwrote {out}")
