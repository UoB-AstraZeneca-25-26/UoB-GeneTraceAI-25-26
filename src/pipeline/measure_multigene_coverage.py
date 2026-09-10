"""
measure_multigene_coverage.py
-----------------------------
Measurements behind MULTIGENE_SPEC.md (A2, A4, A9) and MULTIGENE_COVERAGE.md
(A4, A5). Read-only; writes one JSON of results.

    python src/pipeline/measure_multigene_coverage.py
"""
import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import multi_gene_kleene as K  # noqa: E402

OUT = Path(__file__).resolve().parents[2] / "src" / "pipeline" / "outputs"
OUT.mkdir(parents=True, exist_ok=True)
RNG = np.random.default_rng(42)
res = {}

print("loading...")
E = K.expression()
fe = K.gene_validity()
sequenced, alt, _mod = K.alterations()
lin = K.lineage()
lines = list(E.index)
valid_genes = [g for g in E.columns if fe.get(g, 0) >= K.SILENT_FRAC]
print(f"  {E.shape[0]:,} lines x {E.shape[1]:,} genes; "
      f"{len(valid_genes):,} pass the validity guard "
      f"({len(valid_genes)/E.shape[1]:.1%})")
res["panel"] = {"n_lines": int(E.shape[0]), "n_genes": int(E.shape[1]),
                "n_valid_genes": len(valid_genes),
                "valid_frac": len(valid_genes) / E.shape[1],
                "n_sequenced_lines": len(sequenced & set(lines))}

# ---------------------------------------------------------------- A2 band
print("\nA2 -- how many calls does the uncertainty band move?")
sub = RNG.choice(valid_genes, size=min(2000, len(valid_genes)), replace=False)
V = E[list(sub)]
meas = V.notna()
near = meas & ((V - K.EXPRESSED_MIN).abs() <= K.BAND)
tot = int(meas.values.sum())
res["uncertainty_band"] = {
    "band_log2": K.BAND, "measurement_sd": K.MEAS_SD,
    "n_genes_sampled": len(sub),
    "gene_line_pairs": tot,
    "pairs_in_band": int(near.values.sum()),
    "frac_in_band": float(near.values.sum() / tot),
    "per_gene_frac_median": float(near.sum(axis=0).div(meas.sum(axis=0)).median()),
}
b = res["uncertainty_band"]
print(f"  band = +/-{K.BAND:.3f} log2 (1.96 x measured SD {K.MEAS_SD:.4f})")
print(f"  {b['pairs_in_band']:,} of {b['gene_line_pairs']:,} measured gene x line "
      f"calls fall inside it ({b['frac_in_band']:.2%})")
print(f"  median per gene: {b['per_gene_frac_median']:.2%}")

# ---------------------------------------------------------------- A4 joint
print("\nA4 -- real joint coverage vs the independence prediction")


def term_resolved_mask(kind, gene, negate=False):
    t = K.Term(kind, gene, negate)
    st, _, _ = t.evaluate(lines)
    return (st != K.U).values


single = {}
for g in RNG.choice(valid_genes, size=300, replace=False):
    single[g] = term_resolved_mask("expressed", g)
S = np.vstack(list(single.values()))
marg = S.mean(axis=1)
res["single_term_coverage"] = {"mean": float(marg.mean()),
                               "median": float(np.median(marg)),
                               "n_genes": int(S.shape[0])}
print(f"  single expressed() term: mean {marg.mean():.1%}, "
       f"median {np.median(marg):.1%} of lines resolved")

joint = {}
gl = list(single)
for n_terms in (1, 2, 3, 4):
    obs, ind = [], []
    for _ in range(200):
        pick = RNG.choice(len(gl), size=n_terms, replace=False)
        M = S[pick]
        obs.append(M.all(axis=0).mean())
        ind.append(np.prod(M.mean(axis=1)))
    joint[n_terms] = {"observed_mean": float(np.mean(obs)),
                      "independence_pred": float(np.mean(ind)),
                      "ratio": float(np.mean(obs) / np.mean(ind))}
    print(f"  {n_terms} term(s): observed {np.mean(obs):.1%}   "
          f"independence predicts {np.mean(ind):.1%}   "
          f"ratio {np.mean(obs)/np.mean(ind):.2f}x")
res["joint_coverage_expressed"] = joint

# mixed: expressed(A) AND NOT altered(B)
print("\n  mixed expressed(A) AND NOT altered(B):")
alt_genes = [g for g, s in alt.items() if len(s) >= 10 and g in E.columns]
mix = []
for _ in range(150):
    a = gl[RNG.integers(len(gl))]
    b = alt_genes[RNG.integers(len(alt_genes))]
    ma = single[a]
    mb = term_resolved_mask("altered", b, negate=True)
    mix.append((ma.mean(), mb.mean(), (ma & mb).mean()))
mix = np.array(mix)
res["joint_coverage_mixed"] = {
    "expressed_mean": float(mix[:, 0].mean()),
    "not_altered_mean": float(mix[:, 1].mean()),
    "joint_observed": float(mix[:, 2].mean()),
    "independence_pred": float((mix[:, 0] * mix[:, 1]).mean()),
    "n_pairs": len(mix)}
m = res["joint_coverage_mixed"]
print(f"    expressed {m['expressed_mean']:.1%}  x  NOT altered "
      f"{m['not_altered_mean']:.1%}")
print(f"    joint observed {m['joint_observed']:.1%}   independence predicts "
      f"{m['independence_pred']:.1%}")

# ---------------------------------------------------------------- A5 abstention
print("\nA5 -- abstention trigger rates across the gene universe")
all_genes = list(E.columns)
samp = RNG.choice(all_genes, size=3000, replace=False)
silent = sum(1 for g in samp if fe.get(g, 0) < K.SILENT_FRAC)
near_univ = sum(1 for g in samp if fe.get(g, 0) > 0.95)
res["abstention_single"] = {
    "n_sampled": len(samp),
    "SILENT_GENE_frac": silent / len(samp),
    "NEAR_UNIVERSAL_frac": near_univ / len(samp),
}
print(f"  expressed(G) on a random gene:")
print(f"    SILENT_GENE (fails validity guard) : {silent/len(samp):.1%}")
print(f"    NEAR_UNIVERSAL (>95% of lines)     : {near_univ/len(samp):.1%}")

alt_counts = {g: len(s) for g, s in alt.items()}
samp_alt = RNG.choice(all_genes, size=3000, replace=False)
trivial = sum(1 for g in samp_alt if alt_counts.get(g, 0) < 10)
res["abstention_not_altered"] = {
    "n_sampled": len(samp_alt),
    "TRIVIALLY_SATISFIED_frac": trivial / len(samp_alt),
    "n_genes_with_10plus_altered_lines": int(sum(1 for v in alt_counts.values()
                                                 if v >= 10)),
}
print(f"  NOT altered(G) on a random gene:")
print(f"    TRIVIALLY_SATISFIED (<10 altered lines): {trivial/len(samp_alt):.1%}")
print(f"    genes with >=10 altered lines          : "
      f"{res['abstention_not_altered']['n_genes_with_10plus_altered_lines']:,}")

two = 0
N2 = 1000
for _ in range(N2):
    a, b = RNG.choice(all_genes, size=2, replace=False)
    if fe.get(a, 0) < K.SILENT_FRAC or fe.get(b, 0) < K.SILENT_FRAC:
        two += 1
res["abstention_two_gene"] = {"n_sampled": N2, "any_trigger_frac": two / N2}
print(f"  two-gene expressed(A) AND expressed(B): at least one abstention "
      f"in {two/N2:.1%} of random pairs")

# ---------------------------------------------------------------- per lineage
print("\nper-lineage resolution (expressed(EGFR) AND NOT altered(KRAS))")
terms = K.parse_terms(expressed=["EGFR"], not_altered=["KRAS"])
q = K.query(terms, k_alternatives=0)
cols = q["columns"]
rows = []
for lg, g in cols.groupby("lineage"):
    if len(g) < 10:
        continue
    rows.append({"lineage": lg, "n": int(len(g)),
                 "matches": int((g.state == "TRUE").sum()),
                 "possible": int((g.state == "UNKNOWN").sum()),
                 "excluded": int((g.state == "FALSE").sum()),
                 "resolved_frac": float((g.state != "UNKNOWN").mean())})
rows.sort(key=lambda r: r["resolved_frac"])
res["per_lineage_example"] = rows
print(f"  {'lineage':<26} {'n':>5} {'match':>6} {'poss':>6} {'excl':>6} "
      f"{'resolved':>9}")
for r in rows[:6] + [{"lineage": "...", "n": 0, "matches": 0, "possible": 0,
                      "excluded": 0, "resolved_frac": float('nan')}] + rows[-4:]:
    if r["lineage"] == "...":
        print("    ...")
        continue
    print(f"  {r['lineage']:<26} {r['n']:>5} {r['matches']:>6} "
          f"{r['possible']:>6} {r['excluded']:>6} {r['resolved_frac']:>8.1%}")
res["example_query"] = {"counts": q["counts"], "coverage": q["coverage"]}

# ---------------------------------------------------------------- A9
print("\nA9 -- timing")
res["timing"] = K.benchmark()
print(f"  cold load {res['timing']['cold_load_ms']:.0f} ms   "
      f"warm query {res['timing']['warm_query_ms_median']:.0f} ms "
      f"(median of {res['timing']['n_rep']})")

p = OUT / "multigene_coverage_results.json"
p.write_text(json.dumps(res, indent=2, default=float), encoding="utf-8")
print(f"\nwrote {p}")
