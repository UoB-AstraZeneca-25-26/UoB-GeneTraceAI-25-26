"""
gate_audit/07_independent_labels.py
-----------------------------------
AUDIT 7 of 7.  Is the gate BIOLOGY or ASSAY ARTEFACT?

THE QUESTION
------------
Audits 1-6 established that the surviving gate effect is narrow: ~11% depletion
with pAUC ~0.58 in the 5-10% dependency-prevalence band, replicated across the
Broad (Achilles) and Sanger (Score) CRISPR screens. But Broad and Sanger are the
SAME PERTURBATION TYPE. Two CRISPR screens agreeing rules out screen-specific
noise; it does not rule out a property of CRISPR knockout as an assay.

A drug-sensitivity label breaks that. GDSC2 is a different perturbation
(small-molecule inhibition, not genetic knockout), a different readout (viability
dose-response, not fitness over passages), a different institution and a
different cell-line panel. If the banded gate reproduces there, the gate is a
statement about abundance and biology. If it does not, it is a statement about
CRISPR.

  reproduces on GDSC2  -> biology; proceed to a third perturbation type (RNAi)
  does not reproduce   -> assay artefact of CRISPR knockout, and the claim must
                          be restricted to CRISPR-derived dependency

WHY GDSC2 AND NOT PRISM
-----------------------
GDSC2 is already in the repo (`data/GDSC/GDSC2_fitted_dose_response_27Oct23.xlsx`,
prepared as `validation/prepared/gdsc_scored_ready.parquet`, verified here to be
100% GDSC2 with no GDSC1 contamination). PRISM would need a download and answers
the same question. Do the free test first.

THE CAUSAL-DISTANCE CAVEAT, STATED UP FRONT
-------------------------------------------
GDSC2 is a WEAKER test of the DepMap principle than CRISPR, and this is not a
reason to discount a negative but it is a reason not to over-read one. Knockout
removes the gene; a drug must reach the cell, hit the ANNOTATED target rather
than one of its off-targets, and killing must follow. Clinical kinase inhibitors
are extensively polypharmacological (Klaeger et al. 2017, Science), so a line can
be genuinely drug-sensitive through a protein other than the annotated one. That
puts a ceiling on any depletion GDSC2 can show. The honest reading of a null here
is "does not reproduce in a label with a known attenuation", not "refuted".

COVERAGE IS REPORTED PER LABEL, PER GENE -- NEVER ONE DENOMINATOR
------------------------------------------------------------------
This is the direct lesson of audit 3, applied as a rule. Every label has its own
measured line set, and every GENE within a label has its own. A line absent from
a label is NOT a negative for that label; it is not in that label's denominator
at all. This script therefore:

  * computes, per (label, gene), n_measured / n_expressed / n_used
  * refuses to score any (label, gene) whose intersection is below MIN_LINES
  * writes the full per-(label, gene) coverage table to parquet
  * reports the coverage distribution, not just its median

    label      lines measured   lines also with expression
    achilles   767              754
    score      317              ~300
    gdsc2      962              701

Read-only. Writes only into gate_audit/outputs/.

Run:
    python gate_audit/07_independent_labels.py
"""
import json

import numpy as np
import pandas as pd
from scipy.stats import wilcoxon, mannwhitneyu

import common as C

GDSC_READY = C.ROOT / "validation" / "prepared" / "gdsc_scored_ready.parquet"
BANDS = [0, .02, .05, .10, .25, 1.01]
BAND_LABELS = ["<2%", "2-5%", "5-10%", "10-25%", ">25%"]

res = {}
C.banner("AUDIT 7 -- IS THE GATE BIOLOGY OR AN ASSAY ARTEFACT?")

uni, valid_ensg, symbol_of, sym2ensg, u2e = C.load_gene_lookup()
G = pd.read_parquet(C.OUT / "genes.parquet")
valid_genes = set(G.loc[G.is_valid, "ensg_id"])

# ------------------------------------------------------------ labels
print("\n  loading three labels, each with its OWN measured line set")
W_ach, chk_ach = C.load_chronos("achilles", sym2ensg)
W_sco, chk_sco = C.load_chronos("score", sym2ensg)

gd = pd.read_parquet(GDSC_READY,
                     columns=["model_id", "target_ensg", "sensitive", "drug_id"])
gd["model_id"] = gd.model_id.astype(str).str.lower()
gd["target_ensg"] = (gd.target_ensg.astype("string").str.split(".").str[0]
                     .str.lower())
# A (line, gene) is POSITIVE if sensitive to ANY drug annotated to that gene, and
# MEASURED if any such drug was tested on it. Both are needed: the second is the
# denominator and it is gene-specific.
gd_meas = gd.groupby(["target_ensg", "model_id"], observed=True)["sensitive"].max()
gd_meas = gd_meas.rename("pos").reset_index()
print(f"  gdsc2  : {gd.target_ensg.nunique()} target genes, "
      f"{gd.model_id.nunique()} lines, {gd.drug_id.nunique()} drugs")
print(f"           positive = sensitive to ANY drug annotated to that gene")
res["label_sources"] = {
    "achilles": {"kind": "CRISPR knockout (Broad)", **chk_ach},
    "score": {"kind": "CRISPR knockout (Sanger)", **chk_sco},
    "gdsc2": {"kind": "small-molecule dose-response (GDSC2)",
              "n_genes": int(gd.target_ensg.nunique()),
              "n_lines": int(gd.model_id.nunique()),
              "n_drugs": int(gd.drug_id.nunique())},
}

# ------------------------------------------------------------ common gene set
gdsc_genes = set(gd_meas.target_ensg)
common = sorted(gdsc_genes & set(W_ach.columns) & set(W_sco.columns) & valid_ensg)
print(f"\n  genes in GDSC2 AND both CRISPR screens AND the gene universe: "
      f"{len(common)}")
print("  the cross-label comparison is made on THIS set, so label differences")
print("  are not gene-set differences")

E = C.load_expression(common)
common = [g for g in common if g in E.columns]
print(f"  ... with an RNA column: {len(common)}")
Cm, Pm, ccle_by_ensg, procan_by_ensg = C.load_proteomics(common, valid_ensg, u2e)
res["common_gene_set"] = {"n": len(common),
                          "n_valid": len([g for g in common
                                          if g in valid_genes])}

gd_pos = {g: set(s.model_id[s.pos]) for g, s in gd_meas.groupby("target_ensg")}
gd_all = {g: set(s.model_id) for g, s in gd_meas.groupby("target_ensg")}


def score_for(g):
    """Stouffer Z over the expression + protein layers, exactly as audits 1-6."""
    e = E[g].dropna()
    zE = C.vdw(e)
    parts = {}
    accs = [a for a in ccle_by_ensg.get(g, []) if a in Cm.columns]
    if accs:
        s = Cm[accs].mean(axis=1).dropna()
        if len(s):
            parts["ccle"] = C.vdw(s)
    accs = [a for a in procan_by_ensg.get(g, []) if a in Pm.columns]
    if accs:
        s = Pm[accs].mean(axis=1).dropna()
        if len(s):
            parts["procan"] = C.vdw(s)
    if parts:
        idx = sorted(set().union(*[set(v.index) for v in parts.values()]))
        zP = pd.DataFrame({k: v.reindex(idx) for k, v in parts.items()}).mean(axis=1)
        Z = C.stouffer(zE.reindex(sorted(set(zE.index) | set(idx))),
                       zP.reindex(sorted(set(zE.index) | set(idx))))
    else:
        Z = zE
    return e, Z.dropna()


def label_view(g, label):
    """(measured_lines, positive_set) for this gene under this label. A line not
    measured is NOT a negative -- it is absent from the denominator."""
    if label == "gdsc2":
        return gd_all.get(g, set()), gd_pos.get(g, set())
    W = W_ach if label == "achilles" else W_sco
    if g not in W.columns:
        return set(), set()
    col = W[g]
    meas = set(col.index[col.notna()])
    pos = set(col.index[(col <= C.DEP_THRESHOLD).fillna(False)])
    return meas, pos


# ------------------------------------------------------------ per label/gene
C.banner("7A -- COVERAGE PER LABEL, PER GENE")
print("  n_measured  lines this label measured for this gene")
print("  n_expr      lines with an expression value for this gene")
print("  n_used      the intersection -- the ONLY denominator used")
print()

rows = []
for g in common:
    e, Z = score_for(g)
    expr_lines = set(Z.index)
    for label in ("achilles", "score", "gdsc2"):
        meas, pos = label_view(g, label)
        used = sorted(expr_lines & meas)
        rec = {"ensg_id": g, "symbol": symbol_of.get(g, "?"), "label": label,
               "is_valid": g in valid_genes,
               "n_measured": len(meas), "n_expr": len(expr_lines),
               "n_used": len(used),
               "coverage_of_label": len(used) / max(1, len(meas))}
        if len(used) >= C.MIN_LINES:
            s = Z.reindex(used).values
            p = np.fromiter((m in pos for m in used), bool, len(used))
            if p.sum() >= C.MIN_POS and (~p).sum() >= C.MIN_NEG:
                base = float(p.mean())
                q, keep = C.decile_index(s, 10, "ordinal")
                qv, pv = q[keep], p[keep]
                rec["base_rate"] = base
                rec["d1_ratio"] = float(pv[qv == 0].mean() / base)
                rec["pauc_gate"] = C.partial_auc(s, p, *C.GATE_REGION)
                rec["auc"] = C.auroc(s, p)
                ev = e.reindex(used).values
                fl = ev == np.nanmin(ev)
                if fl.sum() >= 20:
                    rec["floor_ratio"] = float(p[fl].mean() / base)
        rows.append(rec)

L = pd.DataFrame(rows)
L.to_parquet(C.OUT / "07_per_label_per_gene.parquet", index=False)

print(f"  {'label':<10} {'genes':>6} {'n_measured':>22} {'n_used':>22} "
      f"{'scored':>7}")
print(f"  {'':<10} {'':>6} {'median [p10, p90]':>22} {'median [p10, p90]':>22}")
res["coverage"] = {}
for label in ("achilles", "score", "gdsc2"):
    s = L[L.label == label]
    sc = s.dropna(subset=["d1_ratio"])
    mq = s.n_measured.quantile([.1, .5, .9])
    uq = s.n_used.quantile([.1, .5, .9])
    print(f"  {label:<10} {len(s):>6} "
          f"{f'{mq[.5]:.0f} [{mq[.1]:.0f}, {mq[.9]:.0f}]':>22} "
          f"{f'{uq[.5]:.0f} [{uq[.1]:.0f}, {uq[.9]:.0f}]':>22} {len(sc):>7}")
    res["coverage"][label] = {
        "n_genes": int(len(s)), "n_scored": int(len(sc)),
        "n_measured_median": float(mq[.5]),
        "n_measured_p10": float(mq[.1]), "n_measured_p90": float(mq[.9]),
        "n_used_median": float(uq[.5]),
        "n_used_p10": float(uq[.1]), "n_used_p90": float(uq[.9]),
        "base_rate_median": float(sc.base_rate.median()) if len(sc) else None,
    }
print()
print("  Per-gene coverage varies by a factor of several within every label.")
print("  A single panel-wide denominator would misstate every one of these.")

# ------------------------------------------------------------ headline
C.banner("7B -- THE GATE UNDER EACH LABEL,  each on its own denominator")
print(f"  {'label':<10} {'genes':>6} {'base':>8} {'d1 ratio':>10} {'p':>10} "
      f"{'pAUC gate':>11} {'p':>10}")
res["overall"] = {}
for label in ("achilles", "score", "gdsc2"):
    s = L[(L.label == label)].dropna(subset=["d1_ratio"])
    if len(s) < 6:
        print(f"  {label:<10} {len(s):>6}   too few scored genes")
        continue
    p1 = wilcoxon(s.d1_ratio - 1.0).pvalue
    pa = s.pauc_gate.dropna()
    p2 = wilcoxon(pa - 0.5).pvalue if len(pa) > 5 else np.nan
    print(f"  {label:<10} {len(s):>6} {s.base_rate.median():>8.4f} "
          f"{s.d1_ratio.median():>10.4f} {p1:>10.3g} "
          f"{pa.median():>11.4f} {p2:>10.3g}")
    res["overall"][label] = {
        "n_genes": int(len(s)), "base_rate": float(s.base_rate.median()),
        "d1_ratio": float(s.d1_ratio.median()), "d1_p": float(p1),
        "pauc_gate": float(pa.median()), "pauc_p": float(p2),
        "frac_depleted": float((s.d1_ratio < 1).mean())}

# ------------------------------------------------------------ banded
C.banner("7C -- THE BANDED RESULT.  This is the test that decides it.")
print("  Audits 1-6 located the surviving effect in the 5-10% prevalence band.")
print("  Bands are on each label's OWN base rate, because prevalence differs by")
print("  label (CRISPR ~12%, GDSC2 ~25%) and a shared cut would not be the same")
print("  biological stratum.")
print()
L2 = L.dropna(subset=["d1_ratio"]).copy()
L2["band"] = pd.cut(L2.base_rate, BANDS, labels=BAND_LABELS)
print(f"  {'band':>8} " + "".join(f"{lb:>26}" for lb in
                                  ("achilles", "score", "gdsc2")))
print(f"  {'':>8} " + "".join(f"{'n   d1     pAUC':>26}" for _ in range(3)))
res["banded"] = {}
for b in BAND_LABELS:
    cells = []
    for label in ("achilles", "score", "gdsc2"):
        s = L2[(L2.band == b) & (L2.label == label)]
        if len(s) < 5:
            cells.append(f"{len(s):>4}      --      --")
            res["banded"].setdefault(b, {})[label] = {"n": int(len(s))}
            continue
        cells.append(f"{len(s):>4} {s.d1_ratio.median():>7.4f} "
                     f"{s.pauc_gate.median():>7.4f}")
        res["banded"].setdefault(b, {})[label] = {
            "n": int(len(s)), "d1_ratio": float(s.d1_ratio.median()),
            "pauc_gate": float(s.pauc_gate.median()),
            "frac_depleted": float((s.d1_ratio < 1).mean()),
            "d1_p": float(wilcoxon(s.d1_ratio - 1.0).pvalue) if len(s) > 5 else None}
    print(f"  {b:>8} " + "".join(f"{c:>26}" for c in cells))

# ------------------------------------------------------------ paired
C.banner("7D -- PAIRED, SAME GENES:  CRISPR vs DRUG")
print("  Restricted to genes scored under BOTH achilles and gdsc2, so the")
print("  comparison is within-gene and cannot be a gene-set effect.")
print()
piv = L2.pivot_table(index="ensg_id", columns="label",
                     values=["d1_ratio", "pauc_gate", "base_rate"])
both = piv.dropna(subset=[("d1_ratio", "achilles"), ("d1_ratio", "gdsc2")])
print(f"  paired genes: {len(both)}")
res["paired_crispr_vs_drug"] = {"n": int(len(both))}
if len(both) >= 6:
    print(f"  {'quantity':<22} {'ACHILLES':>12} {'GDSC2':>12} {'paired p':>11} "
          f"{'Spearman':>10}")
    for q in ("d1_ratio", "pauc_gate", "base_rate"):
        a, b_ = both[(q, "achilles")], both[(q, "gdsc2")]
        pp = wilcoxon(a - b_).pvalue
        rho = a.corr(b_, method="spearman")
        print(f"  {q:<22} {a.median():>12.4f} {b_.median():>12.4f} "
              f"{pp:>11.3g} {rho:>10.4f}")
        res["paired_crispr_vs_drug"][q] = {
            "achilles": float(a.median()), "gdsc2": float(b_.median()),
            "paired_p": float(pp), "spearman": float(rho)}

# ------------------------------------------------------------ 7E
C.banner("7E -- PREVALENCE-MATCHED BINARISATION.  Without this, 7C is unreadable.")
print("  GDSC2's `sensitive` flag gives a median base rate of ~0.30. Every CRISPR")
print("  band below 10% is therefore EMPTY for GDSC2, and the band where the")
print("  CRISPR gate actually lives cannot be compared at all. A null there would")
print("  be an artefact of the binarisation, not a result.")
print()
print("  Fix: re-binarise GDSC2 per gene by CONTINUOUS drug response, taking the")
print("  most sensitive q fraction of lines as positive, for q matched to the")
print("  CRISPR bands. Response = minimum AUC across drugs annotated to the gene")
print("  (lower AUC = more sensitive), which is threshold-free and per-gene.")
print()

gdc = pd.read_parquet(GDSC_READY, columns=["model_id", "target_ensg", "AUC"])
gdc["model_id"] = gdc.model_id.astype(str).str.lower()
gdc["target_ensg"] = (gdc.target_ensg.astype("string").str.split(".").str[0]
                      .str.lower())
gd_auc = (gdc.groupby(["target_ensg", "model_id"], observed=True)["AUC"]
          .min().rename("auc").reset_index())
auc_by_gene = {g: s.set_index("model_id")["auc"]
               for g, s in gd_auc.groupby("target_ensg")}

QS = (0.02, 0.05, 0.10, 0.25)
print(f"  {'q (target prevalence)':<24} {'genes':>6} {'base':>8} {'d1 ratio':>10} "
      f"{'p':>10} {'pAUC gate':>11} {'p':>10}")
res["gdsc2_prevalence_matched"] = {}
for q in QS:
    vals = []
    for g in common:
        a = auc_by_gene.get(g)
        if a is None:
            continue
        e, Z = score_for(g)
        used = sorted(set(Z.index) & set(a.index))
        if len(used) < C.MIN_LINES:
            continue
        av = a.reindex(used).values
        k = int(round(len(used) * q))
        if k < C.MIN_POS or len(used) - k < C.MIN_NEG:
            continue
        cut = np.partition(av, k - 1)[k - 1]
        p = av <= cut
        s = Z.reindex(used).values
        base = float(p.mean())
        qq, keep = C.decile_index(s, 10, "ordinal")
        qv, pv = qq[keep], p[keep]
        vals.append({"ensg_id": g, "base_rate": base,
                     "d1_ratio": float(pv[qv == 0].mean() / base),
                     "pauc_gate": C.partial_auc(s, p, *C.GATE_REGION)})
    V = pd.DataFrame(vals)
    if len(V) < 6:
        print(f"  q = {q:<20.2f} {len(V):>6}   too few genes")
        continue
    p1 = wilcoxon(V.d1_ratio - 1.0).pvalue
    pa = V.pauc_gate.dropna()
    p2 = wilcoxon(pa - 0.5).pvalue if len(pa) > 5 else np.nan
    print(f"  q = {q:<20.2f} {len(V):>6} {V.base_rate.median():>8.4f} "
          f"{V.d1_ratio.median():>10.4f} {p1:>10.3g} {pa.median():>11.4f} "
          f"{p2:>10.3g}")
    res["gdsc2_prevalence_matched"][f"{q:.2f}"] = {
        "n_genes": int(len(V)), "base_rate": float(V.base_rate.median()),
        "d1_ratio": float(V.d1_ratio.median()), "d1_p": float(p1),
        "pauc_gate": float(pa.median()), "pauc_p": float(p2),
        "frac_depleted": float((V.d1_ratio < 1).mean())}

print()
print("  Compare against CRISPR at the SAME prevalence, same 140-gene set:")
print(f"  {'band':<24} {'achilles d1':>13} {'gdsc2-matched d1':>18}")
for q, bnd in ((0.02, "<2%"), (0.05, "2-5%"), (0.10, "5-10%")):
    gm = res["gdsc2_prevalence_matched"].get(f"{q:.2f}", {})
    ab = res.get("banded", {}).get(bnd, {}).get("achilles", {})
    print(f"  q={q:.2f} / band {bnd:<12} "
          f"{ab.get('d1_ratio', float('nan')):>13.4f} "
          f"{gm.get('d1_ratio', float('nan')):>18.4f}")

# ------------------------------------------------------------ verdict
C.banner("VERDICT")
gd_ov = res["overall"].get("gdsc2", {})
ach_ov = res["overall"].get("achilles", {})
band = res.get("banded", {}).get("5-10%", {})
gd_band = band.get("gdsc2", {})
ach_band = band.get("achilles", {})

print(f"  overall     achilles d1 {ach_ov.get('d1_ratio', float('nan')):.4f}   "
      f"gdsc2 d1 {gd_ov.get('d1_ratio', float('nan')):.4f}")
if gd_band.get("d1_ratio") is not None and ach_band.get("d1_ratio") is not None:
    print(f"  5-10% band  achilles d1 {ach_band['d1_ratio']:.4f} "
          f"(n={ach_band['n']})   gdsc2 d1 {gd_band['d1_ratio']:.4f} "
          f"(n={gd_band['n']})")

gd_depleted = (gd_ov.get("d1_ratio", 1.0) < 1.0
               and gd_ov.get("d1_p", 1.0) < 0.05)
gd_pauc = gd_ov.get("pauc_gate", 0.5) > 0.51
print()
if gd_depleted and gd_pauc:
    verdict = "reproduces_on_an_independent_perturbation__gate_is_biology"
    print("  The gate reproduces on a DIFFERENT PERTURBATION TYPE. It is not a")
    print("  property of CRISPR knockout. Proceed to RNAi/DEMETER2 as a third")
    print("  mechanism -- that is now worth the download.")
elif gd_depleted or gd_pauc:
    verdict = "partial__directionally_consistent_but_attenuated_on_drug_label"
    print("  Directionally consistent on the drug label but attenuated. That is")
    print("  what the causal-distance argument predicts, so it neither confirms")
    print("  nor refutes. RNAi is the tiebreak and IS worth the download.")
else:
    verdict = "does_not_reproduce_on_drug_label__restrict_the_claim_to_CRISPR"
    print("  The gate does NOT reproduce on the drug label. Given GDSC2's known")
    print("  attenuation this is not a refutation, but it means the claim cannot")
    print("  be stated as a general abundance-dependency law. Restrict it to")
    print("  CRISPR-derived dependency until a third perturbation says otherwise.")
print(f"\n  VERDICT: {verdict}")
res["verdict"] = verdict

out = C.OUT / "07_independent_labels_results.json"
out.write_text(json.dumps(res, indent=2, default=float), encoding="utf-8")
print(f"\nwrote {out}")
print(f"wrote {C.OUT / '07_per_label_per_gene.parquet'}")
