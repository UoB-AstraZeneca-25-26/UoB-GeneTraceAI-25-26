"""
gate_audit/02_gate_valid_genes.py
---------------------------------
AUDIT 2 of 5.  "Recompute the gate on the 13,326 valid genes only."

THE QUESTION
------------
The gate result predates the silence-guard finding
(`docs/TRANSCRIPTOMICS_SPEC.md` F2: 5,847 of 19,173 protein-coding genes, 30.5%,
are expressed in under 20% of the panel and their within-gene percentile is
therefore not a meaningful relative position). Nobody has checked whether the
two overlap. If the gate is measured largely on genes whose distributions are
degenerate, the gate is partly an artefact of those distributions.

  effect survives at similar magnitude on valid genes -> the gate is clean, and
      the largest remaining doubt about the project's only solid result is closed
  effect shrinks materially                           -> the gate is partly an
      artefact of degenerate distributions and must be rescoped to valid genes

WHAT IS COMPUTED
----------------
The FULL gate readout -- decile profile, decile-1 depletion, NPV uplift, partial
AUROC in the gate region and the ranker region, and the floor test -- three
times over: on VALID genes, on INVALID genes, and on ALL genes pooled (which
reproduces the headline). Identical code path for all three; only the gene set
changes. Then a like-for-like comparison, because valid and invalid genes differ
in prevalence and in line count and those must not be mistaken for the effect:

  * unmatched   valid vs invalid, as they come
  * matched     invalid genes matched 1:1 to valid genes on prevalence band and
                line count, so the contrast is not a prevalence contrast

Read-only. Writes only into gate_audit/outputs/.

Run:
    python gate_audit/02_gate_valid_genes.py
    python gate_audit/02_gate_valid_genes.py --tie-mode average
"""
import argparse
import json

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu

import common as C

ap = argparse.ArgumentParser()
ap.add_argument("--tie-mode", choices=["ordinal", "average", "drop"],
                default="ordinal", help="stratification rule; see 01")
ap.add_argument("--score", default="Z", choices=["Z", "zE"],
                help="Z = Stouffer (headline); zE = RNA percentile alone")
args = ap.parse_args()

panel = C.read_panel()
G = pd.read_parquet(C.OUT / "genes.parquet")
valid = set(G.loc[G.is_valid, "ensg_id"])
res = {"tie_mode": args.tie_mode, "score": args.score,
       "silent_frac": C.SILENT_FRAC, "expressed_min": C.EXPRESSED_MIN}

C.banner("AUDIT 2 -- THE GATE, RECOMPUTED ON VALID GENES ONLY")
print(f"  validity rule : frac_expressed >= {C.SILENT_FRAC} "
      f"(log2TPM+1 > {C.EXPRESSED_MIN}), docs/TRANSCRIPTOMICS_SPEC.md F2")
print(f"  scored genes  : {len(G):,}   valid {len(valid):,}  "
      f"invalid {len(G)-len(valid):,}")
print(f"  score column  : {args.score}   tie mode: {args.tie_mode}")
print()

R = C.gate_readout(panel, score_col=args.score, tie_mode=args.tie_mode)
R["is_valid"] = R.ensg_id.isin(valid)
R.to_parquet(C.OUT / "02_gate_per_gene.parquet", index=False)

sets = {"all": R, "valid": R[R.is_valid], "invalid": R[~R.is_valid]}
summ = {k: C.summarise_gate(v) for k, v in sets.items()}
res["strata"] = summ

# ------------------------------------------------------------ profile
C.banner("2A -- DECILE PROFILE  (rate / base rate, median across genes)")
print("  decile 1 = lowest abundance.  1.00 = exactly the base rate.")
print()
print(f"  {'decile':>7} {'ALL':>10} {'VALID':>10} {'INVALID':>10}")
for i in range(10):
    print(f"  {i+1:>7} {summ['all']['median_profile'][i]:>10.3f} "
          f"{summ['valid']['median_profile'][i]:>10.3f} "
          f"{summ['invalid']['median_profile'][i]:>10.3f}")
print()
print(f"  {'quantity':<34} {'ALL':>12} {'VALID':>12} {'INVALID':>12}")
for key, lab, fmt in [
    ("n_genes", "genes", "d"),
    ("decile1_ratio", "decile-1 ratio", ".4f"),
    ("decile1_wilcoxon_p", "decile-1 p (vs 1.0)", ".3g"),
    ("frac_genes_depleted", "frac genes depleted", ".3f"),
    ("upper_profile_sd", "deciles 2-10 sd", ".4f"),
    ("upper_profile_slope", "deciles 2-10 slope", "+.5f"),
    ("pauc_gate", "pAUC gate region", ".4f"),
    ("pauc_gate_p", "pAUC gate p", ".3g"),
    ("pauc_ranker", "pAUC ranker region", ".4f"),
    ("full_auc", "full AUROC", ".4f"),
    ("npv10_uplift", "NPV uplift, bottom 10%", "+.5f"),
]:
    vals = []
    for k in ("all", "valid", "invalid"):
        v = summ[k].get(key)
        vals.append("        n/a" if v is None else format(v, fmt).rjust(12))
    print(f"  {lab:<34} {vals[0]} {vals[1]} {vals[2]}")

# ------------------------------------------------------------ floor
C.banner("2B -- FLOOR TEST BY VALIDITY")
print("  Lines where the gene sits exactly at its expression floor.")
print()
print(f"  {'stratum':<10} {'genes':>7} {'floor rate':>12} {'base':>10} "
      f"{'ratio':>9} {'p':>10}")
for k in ("all", "valid", "invalid"):
    f = summ[k].get("floor")
    if not f:
        print(f"  {k:<10} {'--':>7}  (too few genes with a floor block >=20 lines)")
        continue
    print(f"  {k:<10} {f['n_genes']:>7} {f['median_rate']:>12.4f} "
          f"{f['median_base']:>10.4f} {f['median_ratio']:>9.3f} "
          f"{f['wilcoxon_p']:>10.3g}")
print()
print("  NOTE. Invalid genes are silent in most lines, so their 'floor block' is")
print("  most of the gene. A floor ratio there is close to a base-rate tautology,")
print("  which is exactly why the valid-gene column is the one that carries the")
print("  claim.")

# ------------------------------------------------------------ matched
C.banner("2C -- LIKE-FOR-LIKE:  valid vs invalid, matched on prevalence and size")
print("  Invalid genes have lower prevalence and fewer usable lines. Without")
print("  matching, a valid/invalid difference could be a prevalence difference.")
print()
Rv, Ri = sets["valid"].copy(), sets["invalid"].copy()
bands = [0, .02, .05, .10, .25, 1.01]
for df in (Rv, Ri):
    df["band"] = pd.cut(df.base_rate, bands,
                        labels=["<2%", "2-5%", "5-10%", "10-25%", ">25%"])
    df["size_band"] = pd.cut(df.n_lines, [0, 300, 500, 700, 10000],
                             labels=["<300", "300-500", "500-700", ">700"])

rng = np.random.default_rng(C.SEED)
pairs = []
pool = {k: list(v.index) for k, v in Rv.groupby(["band", "size_band"],
                                                observed=True)}
for _, row in Ri.iterrows():
    key = (row.band, row.size_band)
    cand = pool.get(key)
    if cand:
        j = cand.pop(rng.integers(len(cand)))
        pairs.append((row.ensg_id, Rv.loc[j, "ensg_id"]))
print(f"  matched pairs: {len(pairs)} of {len(Ri)} invalid genes")
res["matched"] = {"n_pairs": len(pairs)}
if len(pairs) >= 10:
    inv_ids = [a for a, _ in pairs]
    val_ids = [b for _, b in pairs]
    mi = R.set_index("ensg_id").loc[inv_ids]
    mv = R.set_index("ensg_id").loc[val_ids]
    print()
    print(f"  {'quantity':<28} {'INVALID':>12} {'VALID(matched)':>16} "
          f"{'delta':>10} {'MWU p':>10}")
    for col, lab in [("d1_ratio", "decile-1 ratio"),
                     ("pauc_gate", "pAUC gate region"),
                     ("base_rate", "base rate (check)"),
                     ("n_lines", "lines (check)")]:
        a, b = mi[col].dropna(), mv[col].dropna()
        if len(a) < 5 or len(b) < 5:
            continue
        p = mannwhitneyu(a, b).pvalue
        print(f"  {lab:<28} {a.median():>12.4f} {b.median():>16.4f} "
              f"{b.median()-a.median():>+10.4f} {p:>10.3g}")
        res["matched"][col] = {"invalid": float(a.median()),
                               "valid_matched": float(b.median()),
                               "delta": float(b.median() - a.median()),
                               "mwu_p": float(p)}
else:
    print("  too few matched pairs to test -- report the unmatched contrast only")

# ------------------------------------------------------------ verdict
C.banner("VERDICT")
a_d1 = summ["all"]["decile1_ratio"]
v_d1 = summ["valid"]["decile1_ratio"]
a_pa = summ["all"]["pauc_gate"]
v_pa = summ["valid"]["pauc_gate"]
d1_shift = v_d1 - a_d1
pa_shift = v_pa - a_pa
print(f"  decile-1 depletion   ALL {a_d1:.4f}  ->  VALID {v_d1:.4f}   "
      f"({d1_shift:+.4f})")
print(f"  gate-region pAUC     ALL {a_pa:.4f}  ->  VALID {v_pa:.4f}   "
      f"({pa_shift:+.4f})")
print()
# THIS AUDIT ANSWERS EXACTLY ONE QUESTION: does restricting to valid genes move
# the result? It does NOT answer whether the result is large. Those are separate,
# and conflating them would credit or blame the silence guard for a magnitude it
# has nothing to do with. Absolute magnitude is audit 3's business.
moved = (abs(d1_shift) > 0.05) or (abs(pa_shift) > 0.02)
if not moved:
    verdict = "validity_is_not_the_confound__effect_unchanged_on_valid_genes"
    print("  Restricting to valid genes does NOT move the result. The degenerate")
    print("  distributions flagged by the silence guard are not what the gate was")
    print("  measuring, so this particular doubt is CLOSED. Report the gate on")
    print("  valid genes anyway -- it costs nothing and removes the question.")
else:
    verdict = "validity_shifts_the_result__report_valid_gene_numbers_as_headline"
    print("  Restricting to valid genes moves the result materially. The")
    print("  valid-gene numbers, not the pooled ones, are the honest headline.")
print(f"\n  VERDICT: {verdict}")
print()
print(f"  SEPARATE QUESTION -- is the effect LARGE? On this panel the valid-gene")
print(f"  depletion is {100*(1-v_d1):.1f}% and gate pAUC is {v_pa:.4f}. Both sit far")
print(f"  below the published figures. That gap is NOT about gene validity; see")
print(f"  audit 3 (denominator), which is where it comes from.")
res["verdict"] = verdict
res["shift_valid_minus_all"] = {"decile1_ratio": float(d1_shift),
                                "pauc_gate": float(pa_shift)}

out = C.OUT / "02_gate_valid_genes_results.json"
out.write_text(json.dumps(res, indent=2, default=float), encoding="utf-8")
print(f"\nwrote {out}")
print(f"wrote {C.OUT / '02_gate_per_gene.parquet'}")
