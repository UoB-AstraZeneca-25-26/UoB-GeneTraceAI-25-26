"""
gate_audit/01_tie_mass.py
-------------------------
AUDIT 1 of 5.  "Report the tie mass in the expression percentile and the
tie-handling method used."

THE METHOD, ESTABLISHED BY INSPECTION NOT ASSUMPTION
----------------------------------------------------
There are TWO tie-handling decisions in the gate, in series, and they are not
the same decision. Both are in `src/pipeline/test_run_abundance_gate.py`:

  (1) THE PERCENTILE.   `vdw()` -> `s.rank(method="min")`
      Ties share a rank, so every line in a tie block receives the SAME z.
      This is the "average rank" branch of the concern in spirit -- the block
      does not spread, it collapses to one point -- except the shared rank is
      the MINIMUM of the block, not its mean, so a floor block sits at the
      BOTTOM of the score, not in the middle. (min-rank puts an 800-line floor
      block at z = ppf(1/(n+1)); average-rank would put it at ppf(400/(n+1)),
      i.e. near the median.)

  (2) THE DECILE.       `rankdata(s, method="ordinal")`
      This is where the block spreads. Ordinal ranking assigns 1..k to k tied
      values in ARRAY ORDER, so a tie block larger than a decile is smeared
      across as many deciles as it spans, and which line lands in decile 1
      versus decile 5 is determined by cell-line index order and nothing else.

So the answer to "which of the two is it" is: BOTH, at different stages. The
percentile is average-like (collapsing), the decile assignment is ordinal
(spreading). The score is not contaminated; the STRATIFICATION is.

That matters because the headline decile-1 ratio is computed on the ordinal
strata. This script measures how much.

WHAT IT REPORTS
---------------
  A. Tie structure of the expression vector, per gene and pooled: tie mass,
     largest block, floor block, and the same split by gene validity.
  B. Decile contamination: what fraction of decile 1 is a tie block, and how
     many deciles the largest block spans under the production `ordinal` rule.
  C. The decile-1 ratio recomputed under all three tie rules. If `ordinal` and
     `average` agree, the smearing is cosmetic. If they diverge, the headline
     depends on cell-line array order.

Read-only. Writes only into gate_audit/outputs/.

Run:
    python gate_audit/01_tie_mass.py
"""
import json

import numpy as np
import pandas as pd
from scipy.stats import wilcoxon

import common as C

panel = C.read_panel()
G = pd.read_parquet(C.OUT / "genes.parquet")
res = {"tie_handling_method": {
    "percentile": "pandas.Series.rank(method='min') inside vdw() -- ties COLLAPSE "
                  "to one shared z at the block minimum",
    "decile": "scipy.stats.rankdata(method='ordinal') -- ties SPREAD across "
              "deciles in array order",
    "roc_and_pauc": "equal scores collapse into a single ROC step, so partial "
                    "AUROC is tie-safe and is NOT affected by (2)",
}}

# ============================================================ A
C.banner("AUDIT 1A -- TIE MASS IN THE EXPRESSION PERCENTILE")
print("  tie mass = share of a gene's lines that share their value with >=1 other")
print()
for col, lab in [("tie_tie_mass", "tie mass"),
                 ("tie_largest_block_frac", "largest tie block (share of lines)"),
                 ("tie_floor_block_frac", "floor block (share of lines)")]:
    v = G[col]
    print(f"  {lab:<40} median {v.median():.3f}   mean {v.mean():.3f}   "
          f"p90 {v.quantile(.90):.3f}   max {v.max():.3f}")
print()
print(f"  genes with >50% of lines tied      : {100*(G.tie_tie_mass > .50).mean():5.1f}%")
print(f"  genes with >20% of lines tied      : {100*(G.tie_tie_mass > .20).mean():5.1f}%")
print(f"  genes whose largest block > 1 decile: "
      f"{100*(G.tie_largest_block_frac > .10).mean():5.1f}%")
print(f"  genes with NO ties at all          : {100*(G.tie_tie_mass == 0).mean():5.1f}%")

print()
print("  BY VALIDITY  (valid = frac_expressed >= %.2f)" % C.SILENT_FRAC)
print(f"  {'stratum':<12} {'genes':>6} {'tie mass':>10} {'largest blk':>12} "
      f"{'floor blk':>10} {'base rate':>10}")
by_valid = {}
for flag, name in [(True, "valid"), (False, "invalid")]:
    s = G[G.is_valid == flag]
    if not len(s):
        continue
    print(f"  {name:<12} {len(s):>6} {s.tie_tie_mass.median():>10.3f} "
          f"{s.tie_largest_block_frac.median():>12.3f} "
          f"{s.tie_floor_block_frac.median():>10.3f} "
          f"{s.base_rate.median():>10.4f}")
    by_valid[name] = {"n_genes": int(len(s)),
                      "tie_mass_median": float(s.tie_tie_mass.median()),
                      "largest_block_frac_median": float(s.tie_largest_block_frac.median()),
                      "floor_block_frac_median": float(s.tie_floor_block_frac.median())}
res["tie_mass"] = {
    "median": float(G.tie_tie_mass.median()),
    "mean": float(G.tie_tie_mass.mean()),
    "p90": float(G.tie_tie_mass.quantile(.90)),
    "frac_genes_over_50pct": float((G.tie_tie_mass > .50).mean()),
    "frac_genes_over_20pct": float((G.tie_tie_mass > .20).mean()),
    "frac_genes_block_over_one_decile": float((G.tie_largest_block_frac > .10).mean()),
    "frac_genes_no_ties": float((G.tie_tie_mass == 0).mean()),
    "by_validity": by_valid,
}

# ============================================================ B
C.banner("AUDIT 1B -- HOW MUCH OF DECILE 1 IS AN ARBITRARILY-ORDERED TIE BLOCK?")
print("  Under the production rule (ordinal), a tie block larger than n/10 lines")
print("  is split across deciles by cell-line array order.")
print()

span_rows = []
for g, d in panel.groupby("ensg_id", sort=False):
    s = d.Z.values
    q, keep = C.decile_index(s, 10, "ordinal")
    if keep.sum() < C.MIN_LINES:
        continue
    sv = s[keep]
    qv = q[keep]
    vals, counts = np.unique(sv, return_counts=True)
    big = vals[counts.argmax()]
    blk = sv == big
    d1 = qv == 0
    span_rows.append({
        "ensg_id": g,
        "largest_block_n": int(blk.sum()),
        "largest_block_spans_deciles": int(len(np.unique(qv[blk]))),
        "d1_share_that_is_largest_block": float((d1 & blk).sum() / max(1, d1.sum())),
        "d1_share_tied": float(np.isin(sv[d1], vals[counts > 1]).mean()),
        "largest_block_is_floor": bool(big == sv.min()),
    })
S = pd.DataFrame(span_rows)
S.to_parquet(C.OUT / "01_decile_contamination.parquet", index=False)

print(f"  {'quantity':<52} {'median':>8} {'mean':>8} {'p90':>8}")
for col, lab in [("largest_block_spans_deciles", "deciles spanned by the largest tie block"),
                 ("d1_share_that_is_largest_block", "share of decile 1 from that one block"),
                 ("d1_share_tied", "share of decile 1 sitting in ANY tie")]:
    v = S[col]
    print(f"  {lab:<52} {v.median():>8.3f} {v.mean():>8.3f} {v.quantile(.9):>8.3f}")
print()
print(f"  genes whose largest block spans >1 decile : "
      f"{100*(S.largest_block_spans_deciles > 1).mean():5.1f}%")
print(f"  genes whose largest block spans >3 deciles: "
      f"{100*(S.largest_block_spans_deciles > 3).mean():5.1f}%")
print(f"  genes where decile 1 is >50% one tie block: "
      f"{100*(S.d1_share_that_is_largest_block > .5).mean():5.1f}%")
print(f"  ... and that block is the expression floor : "
      f"{100*S.largest_block_is_floor.mean():5.1f}%")
res["decile_contamination"] = {
    "median_deciles_spanned": float(S.largest_block_spans_deciles.median()),
    "mean_deciles_spanned": float(S.largest_block_spans_deciles.mean()),
    "frac_spanning_gt1": float((S.largest_block_spans_deciles > 1).mean()),
    "frac_spanning_gt3": float((S.largest_block_spans_deciles > 3).mean()),
    "median_d1_share_from_one_block": float(S.d1_share_that_is_largest_block.median()),
    "median_d1_share_tied": float(S.d1_share_tied.median()),
    "frac_largest_block_is_floor": float(S.largest_block_is_floor.mean()),
}

# ============================================================ C
C.banner("AUDIT 1C -- DOES THE DECILE-1 RATIO DEPEND ON THE TIE RULE?")
print("  Same score, same labels, same genes. Only the stratification changes.")
print("  If these agree, the smearing is cosmetic. If they diverge, the headline")
print("  depends on cell-line array order.")
print()

prof = {m: [] for m in ("ordinal", "average", "drop")}
d1 = {m: [] for m in ("ordinal", "average", "drop")}
kept = {m: [] for m in ("ordinal", "average", "drop")}
for g, d in panel.groupby("ensg_id", sort=False):
    s = d.Z.values
    p = d.dep.values.astype(bool)
    base = p.mean()
    if base <= 0:
        continue
    for m in prof:
        q, keep = C.decile_index(s, 10, m)
        if keep.sum() < C.MIN_LINES or p[keep].sum() < C.MIN_POS:
            prof[m].append([np.nan] * 10)
            d1[m].append(np.nan)
            kept[m].append(np.nan)
            continue
        qv, pv = q[keep], p[keep]
        rates = np.array([pv[qv == i].mean() if (qv == i).sum() else np.nan
                          for i in range(10)])
        prof[m].append(rates / base)
        d1[m].append(rates[0] / base)
        kept[m].append(keep.mean())

print(f"  {'decile':>7} {'ordinal':>10} {'average':>10} {'drop-ties':>11}")
med = {}
for m in prof:
    med[m] = np.nanmedian(np.vstack(prof[m]), axis=0)
for i in range(10):
    print(f"  {i+1:>7} {med['ordinal'][i]:>10.3f} {med['average'][i]:>10.3f} "
          f"{med['drop'][i]:>11.3f}")
print()
res["tie_rule_sensitivity"] = {}
for m in prof:
    v = pd.Series(d1[m]).dropna()
    w = wilcoxon(v - 1.0) if len(v) > 5 else None
    print(f"  {m:<10} decile-1 ratio {v.median():.4f}   "
          f"depleted in {100*(v < 1).mean():.1f}% of genes   "
          f"lines retained {100*np.nanmedian(kept[m]):.1f}%   "
          f"p {w.pvalue if w else np.nan:.3g}")
    res["tie_rule_sensitivity"][m] = {
        "decile1_ratio_median": float(v.median()),
        "profile": [float(x) for x in med[m]],
        "frac_genes_depleted": float((v < 1).mean()),
        "median_lines_retained": float(np.nanmedian(kept[m])),
        "wilcoxon_p": float(w.pvalue) if w else None,
        "n_genes": int(len(v)),
    }

o = pd.Series(d1["ordinal"])
a = pd.Series(d1["average"])
both = o.notna() & a.notna()
delta = (a[both] - o[both])
print()
print(f"  paired shift (average - ordinal) : median {delta.median():+.4f}   "
      f"p {wilcoxon(delta).pvalue:.3g}" if len(delta) > 5 else "")
res["tie_rule_sensitivity"]["paired_average_minus_ordinal"] = {
    "median": float(delta.median()),
    "p": float(wilcoxon(delta).pvalue) if len(delta) > 5 else None,
    "n": int(both.sum()),
}

# ============================================================ verdict
C.banner("VERDICT")
lo, hi = min(v["decile1_ratio_median"] for v in
             (res["tie_rule_sensitivity"][m] for m in prof)), \
         max(v["decile1_ratio_median"] for v in
             (res["tie_rule_sensitivity"][m] for m in prof))
spread = hi - lo
if spread < 0.02:
    verdict = "tie_rule_immaterial__headline_is_not_an_ordering_artefact"
elif spread < 0.05:
    verdict = "tie_rule_minor__report_the_rule_and_the_spread"
else:
    verdict = "tie_rule_MATERIAL__the_headline_depends_on_array_order"
print(f"  decile-1 ratio across the three tie rules: {lo:.4f} .. {hi:.4f}  "
      f"(spread {spread:.4f})")
print(f"  VERDICT: {verdict}")
res["verdict"] = verdict
res["decile1_ratio_spread_across_tie_rules"] = float(spread)

out = C.OUT / "01_tie_mass_results.json"
out.write_text(json.dumps(res, indent=2, default=float), encoding="utf-8")
print(f"\nwrote {out}")
print(f"wrote {C.OUT / '01_decile_contamination.parquet'}")
