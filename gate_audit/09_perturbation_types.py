"""
gate_audit/09_perturbation_types.py
-----------------------------------
AUDIT 9.  THREE PERTURBATION TYPES.  The tiebreak.

THE LOGIC
---------
  audit 6   Broad vs Sanger CRISPR agree            -> not screen-specific noise
  audit 7   GDSC2 drug label is FLAT                -> either the gate is specific
                                                       to CRISPR, or GDSC2's
                                                       causal distance hides it
  audit 9   RNAi (DEMETER2) decides between those two

RNAi is the discriminating third case because of where it sits:

                       genetic?   complete?   acts on the gene product directly?
    CRISPR knockout      yes         yes                    yes
    RNAi knockdown       yes         NO (partial)           yes
    Drug (GDSC2)         no          no          NO (annotated target may be wrong)

If the gate reproduces on RNAi, the drug null was causal distance and the gate is
a statement about gene abundance and dependency. If the gate is flat on RNAi too,
then it is specific to CRISPR knockout, and the most likely reason is mechanical
rather than biological -- CRISPR cutting, copy-number-driven cutting toxicity,
and guide-level effects all scale with things that correlate with expression.

THE SCALE PROBLEM, AND WHY PREVALENCE MATCHING IS THE ONLY HONEST COMPARISON
----------------------------------------------------------------------------
DEP_THRESHOLD = -0.5 is a Chronos quantity. It does NOT transfer:

    CRISPR 24Q4  curated essentials median  -0.98   (anchor is ~-1)
    DEMETER2     curated essentials median  -0.35   (knockdown is PARTIAL)
    GDSC2        no gene-effect scale at all

Applying one absolute cut across three assays would compare three different
prevalences and call the difference biology. So every label here is binarised the
SAME way: per gene, the most-dependent q fraction of that gene's own measured
lines is positive, for q in {0.02, 0.05, 0.10, 0.25}. That makes prevalence
identical by construction, which is exactly what the banded result needs.

Each label's native threshold is reported alongside for reference, but the
cross-label claim rests on the matched binarisation.

COVERAGE IS PER LABEL, PER GENE -- NEVER POOLED
------------------------------------------------
A line absent from a label is not a negative for that label; it is not in that
label's denominator. Enforced per (label, gene), and the coverage table is
written out.

Read-only. Writes only into gate_audit/outputs/.

Run:
    python gate_audit/09_perturbation_types.py
    python gate_audit/09_perturbation_types.py --n-genes 6000
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

QS = (0.02, 0.05, 0.10, 0.25)
GDSC_READY = C.ROOT / "validation" / "prepared" / "gdsc_scored_ready.parquet"
res = {}

C.banner("AUDIT 9 -- CRISPR vs RNAi vs DRUG")
uni, valid_ensg, symbol_of, sym2ensg, u2e = C.load_gene_lookup()

W_cr, chk_cr = C.load_crispr_24q4(sym2ensg)
W_ri, chk_ri = C.load_demeter2(sym2ensg, target_fpr=chk_cr["label_fpr"])
E = C.load_expression_24q4(sym2ensg)

gdc = pd.read_parquet(GDSC_READY, columns=["model_id", "target_ensg", "AUC"])
gdc["model_id"] = gdc.model_id.astype(str).str.lower()
gdc["target_ensg"] = (gdc.target_ensg.astype("string").str.split(".").str[0]
                      .str.lower())
gd = (gdc.groupby(["target_ensg", "model_id"], observed=True)["AUC"]
      .min().rename("v").reset_index())
gd_by_gene = {g: s.set_index("model_id")["v"] for g, s in gd.groupby("target_ensg")}
res["scale_checks"] = {"crispr_24q4": chk_cr, "rnai_demeter2": chk_ri}

print(f"\n  NOTE: DEMETER2 curated essentials sit at "
      f"{chk_ri['median_essential']:+.3f}, not ~-1. Knockdown is partial.")
print(f"  DEP_THRESHOLD {C.DEP_THRESHOLD} is NOT applied to it. All three labels")
print(f"  are binarised by matched per-gene prevalence instead.")

# ------------------------------------------------------------ gene sets
cr_g, ri_g, gd_g = set(W_cr.columns), set(W_ri.columns), set(gd_by_gene)
genetic = sorted(cr_g & ri_g & set(E.columns) & valid_ensg)
rng = np.random.default_rng(C.SEED)
if args.n_genes < len(genetic):
    genetic = sorted(rng.choice(genetic, args.n_genes, replace=False))
all3 = sorted(set(genetic) & gd_g)
print(f"\n  genes in CRISPR AND RNAi AND expression : {len(genetic):,}")
print(f"  ... also in GDSC2 (drug-target genes)   : {len(all3)}")
res["gene_sets"] = {"crispr_rnai": len(genetic), "all_three": len(all3)}


def values_for(g, label):
    """(measured Series indexed by model_id, direction) -- LOWER = more dependent
    for all three, so one comparator works throughout."""
    if label == "crispr":
        s = W_cr[g] if g in W_cr.columns else None
    elif label == "rnai":
        s = W_ri[g] if g in W_ri.columns else None
    else:
        s = gd_by_gene.get(g)
    return None if s is None else s.dropna()


# ------------------------------------------------------------ coverage + gate
C.banner("9A -- COVERAGE PER LABEL, PER GENE")
rows = []
for i, g in enumerate(genetic):
    if i and i % 1000 == 0:
        print(f"    {i:,}/{len(genetic):,}")
    e = E[g].dropna()
    if len(e) < C.MIN_LINES:
        continue
    z = C.vdw(e)
    for label in ("crispr", "rnai", "gdsc2"):
        v = values_for(g, label)
        if v is None:
            continue
        used = sorted(set(z.index) & set(v.index))
        rec = {"ensg_id": g, "label": label, "n_measured": len(v),
               "n_expr": len(z), "n_used": len(used)}
        if len(used) >= C.MIN_LINES:
            s = z.reindex(used).values
            vv = v.reindex(used).values
            for q in QS:
                k = int(round(len(used) * q))
                if k < C.MIN_POS or len(used) - k < C.MIN_NEG:
                    continue
                cut = np.partition(vv, k - 1)[k - 1]
                p = vv <= cut
                if p.sum() < C.MIN_POS or (~p).sum() < C.MIN_NEG:
                    continue
                base = float(p.mean())
                qq, keep = C.decile_index(s, 10, "ordinal")
                qv, pv = qq[keep], p[keep]
                rec[f"d1_q{q}"] = float(pv[qv == 0].mean() / base)
                rec[f"pauc_q{q}"] = C.partial_auc(s, p, *C.GATE_REGION)
                rec[f"base_q{q}"] = base
        rows.append(rec)

L = pd.DataFrame(rows)
L.to_parquet(C.OUT / "09_perturbation_per_gene.parquet", index=False)

print(f"\n  {'label':<10} {'genes':>7} {'n_measured med [p10,p90]':>28} "
      f"{'n_used med [p10,p90]':>26}")
res["coverage"] = {}
for label in ("crispr", "rnai", "gdsc2"):
    s = L[L.label == label]
    if not len(s):
        continue
    m = s.n_measured.quantile([.1, .5, .9])
    u = s.n_used.quantile([.1, .5, .9])
    print(f"  {label:<10} {len(s):>7,} "
          f"{f'{m[.5]:.0f} [{m[.1]:.0f}, {m[.9]:.0f}]':>28} "
          f"{f'{u[.5]:.0f} [{u[.1]:.0f}, {u[.9]:.0f}]':>26}")
    res["coverage"][label] = {"n_genes": int(len(s)),
                              "n_measured_median": float(m[.5]),
                              "n_used_median": float(u[.5]),
                              "n_used_p10": float(u[.1]),
                              "n_used_p90": float(u[.9])}

# ------------------------------------------------------------ matched
C.banner("9B -- MATCHED PREVALENCE.  Same q, same genes where possible.")
print("  q is the fraction of each gene's own measured lines called dependent,")
print("  so base rate is identical across labels by construction.")
print()
print(f"  {'q':>6} {'label':<10} {'genes':>7} {'d1 ratio':>10} {'p':>11} "
      f"{'pAUC gate':>11} {'p':>11}")
res["matched"] = {}
for q in QS:
    for label in ("crispr", "rnai", "gdsc2"):
        s = L[(L.label == label)].dropna(subset=[f"d1_q{q}"])
        if len(s) < 6:
            print(f"  {q:>6.2f} {label:<10} {len(s):>7}   too few")
            continue
        d1 = s[f"d1_q{q}"]
        pa = s[f"pauc_q{q}"].dropna()
        p1 = wilcoxon(d1 - 1.0).pvalue
        p2 = wilcoxon(pa - 0.5).pvalue if len(pa) > 5 else np.nan
        print(f"  {q:>6.2f} {label:<10} {len(s):>7,} {d1.median():>10.4f} "
              f"{p1:>11.3g} {pa.median():>11.4f} {p2:>11.3g}")
        res["matched"].setdefault(f"{q:.2f}", {})[label] = {
            "n": int(len(s)), "d1_ratio": float(d1.median()), "d1_p": float(p1),
            "pauc": float(pa.median()), "pauc_p": float(p2),
            "frac_depleted": float((d1 < 1).mean())}
    print()

# ------------------------------------------------------------ paired
C.banner("9C -- PAIRED WITHIN GENE:  CRISPR vs RNAi")
print("  Same genes, same prevalence, same score. Only the perturbation differs.")
print()
res["paired_crispr_rnai"] = {}
for q in QS:
    piv = L.pivot_table(index="ensg_id", columns="label", values=f"d1_q{q}")
    if "crispr" not in piv or "rnai" not in piv:
        continue
    b = piv.dropna(subset=["crispr", "rnai"])
    if len(b) < 6:
        continue
    p = wilcoxon(b.crispr - b.rnai).pvalue
    rho = b.crispr.corr(b.rnai, method="spearman")
    print(f"  q={q:.2f}  n={len(b):>5,}  CRISPR {b.crispr.median():.4f}   "
          f"RNAi {b.rnai.median():.4f}   paired p {p:.3g}   Spearman {rho:+.4f}")
    res["paired_crispr_rnai"][f"{q:.2f}"] = {
        "n": int(len(b)), "crispr": float(b.crispr.median()),
        "rnai": float(b.rnai.median()), "paired_p": float(p),
        "spearman": float(rho)}

# ------------------------------------------------------------ verdict
C.banner("VERDICT")
key = "0.05"
m = res["matched"].get(key, {})
cr = m.get("crispr", {})
ri = m.get("rnai", {})
gd_ = m.get("gdsc2", {})
print(f"  at matched prevalence q=0.05:")
for nm, d in (("CRISPR", cr), ("RNAi", ri), ("drug (GDSC2)", gd_)):
    if d:
        print(f"    {nm:<14} d1 {d['d1_ratio']:.4f} (p {d['d1_p']:.3g})   "
              f"pAUC {d['pauc']:.4f} (p {d['pauc_p']:.3g})   n {d['n']:,}")
print()
rnai_works = bool(ri) and ri["d1_ratio"] < 0.95 and ri["d1_p"] < 0.05
crispr_works = bool(cr) and cr["d1_ratio"] < 0.95 and cr["d1_p"] < 0.05
if crispr_works and rnai_works:
    verdict = "reproduces_on_a_second_genetic_perturbation__gate_is_biology"
    print("  The gate reproduces under RNAi knockdown as well as CRISPR knockout.")
    print("  Two independent genetic perturbations agree, so the drug null in")
    print("  audit 7 is causal distance, not a refutation. The gate is a")
    print("  statement about abundance and dependency.")
elif crispr_works and not rnai_works:
    verdict = "CRISPR_only__gate_does_not_generalise_beyond_knockout"
    print("  The gate appears under CRISPR knockout and NOT under RNAi knockdown,")
    print("  having already failed under drug. Two of three perturbations say no.")
    print("  The most likely explanation is mechanical rather than biological --")
    print("  cutting toxicity and guide-level effects scale with quantities that")
    print("  correlate with expression. The claim must be restricted to")
    print("  CRISPR-derived dependency and that restriction stated prominently.")
elif rnai_works and not crispr_works:
    verdict = "rnai_only__unexpected__investigate_before_reporting"
    print("  RNAi shows the gate and CRISPR does not, at matched prevalence.")
    print("  That inverts the expected ordering and should be investigated")
    print("  before any of it is reported.")
else:
    verdict = "no_perturbation_shows_the_gate_at_matched_prevalence"
    print("  At matched prevalence no perturbation shows the gate. The banded")
    print("  result does not survive putting all three labels on equal footing.")
print(f"\n  VERDICT: {verdict}")
res["verdict"] = verdict

out = C.OUT / "09_perturbation_types_results.json"
out.write_text(json.dumps(res, indent=2, default=float), encoding="utf-8")
print(f"\nwrote {out}")
print(f"wrote {C.OUT / '09_perturbation_per_gene.parquet'}")
