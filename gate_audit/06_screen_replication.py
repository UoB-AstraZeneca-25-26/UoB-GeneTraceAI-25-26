"""
gate_audit/06_screen_replication.py
-----------------------------------
AUDIT 6 of 6.  "Re-verify on the real labels" — and use the two screens as
independent replications.

WHAT WAS ACTUALLY ALREADY DONE, AND WHAT WAS NOT
------------------------------------------------
The list this audit came from says the gate has not been re-run on real Chronos.
That is not quite right, and the record should be corrected rather than
repeated: `src/pipeline/outputs/test_run_abundance_gate_chronos_results.json`
carries `label_fpr = 0.00928` and reference-set medians, which only the real
Chronos branch of that script produces. The gate HAS been run on real Chronos.

What has NOT been done is the part that matters here:

  * the gate on real Chronos with the DENOMINATOR CORRECTED (audit 3)
  * the gate on real Chronos restricted to VALID genes (audit 2)
  * Broad (Achilles) versus Sanger (Score) as two independent screens, under
    those corrections, which is what answers the post-hoc region challenge in
    audit 5

This script does those three things.

THE TWO SCREENS
---------------
    data/DepMap_Chronos/GeneFitnessEffect_Chronos_Achilles.hdf5   Broad
    data/DepMap_Chronos/GeneFitnessEffect_Chronos_Score.hdf5      Sanger

Different institutions, different libraries, different cell-line panels, same
Chronos processing. Genes and lines are intersected so the comparison is paired:
the same gene measured twice. A paired test is the right one — an unpaired
comparison of two medians would confound the screens with their gene sets.

THE QUESTION A PANEL WILL ASK, AND THE ANSWER
----------------------------------------------
  "Isn't this just saying you can't knock out a gene that isn't expressed?"

Partly yes, and that is why it replicates. Near-tautological findings are
reliable and unexciting in equal measure. The answer is that the percentile is
NOT merely an encoded detection boolean: bottom decile INTERSECT
rule-says-EXPRESSED still shows depletion, stable across calibration rules
(`test_run_detection_vs_percentile.py`, ratio 0.810-0.827 across six rules).
This script re-runs that intersection here, on both screens and on the corrected
denominator, so the answer travels with the result instead of sitting in a
robustness appendix.

Read-only. Writes only into gate_audit/outputs/.

Run:
    python gate_audit/06_screen_replication.py
"""
import json

import numpy as np
import pandas as pd
from scipy.stats import wilcoxon

import common as C

res = {}
C.banner("AUDIT 6 -- BROAD vs SANGER ON THE REAL CHRONOS RELEASE")

uni, valid_ensg, symbol_of, sym2ensg, u2e = C.load_gene_lookup()
G = pd.read_parquet(C.OUT / "genes.parquet")
valid_genes = set(G.loc[G.is_valid, "ensg_id"])
panel = C.read_panel()
genes = sorted(panel.ensg_id.unique())

Ws, chks = {}, {}
for which in ("achilles", "score"):
    Ws[which], chks[which] = C.load_chronos(which, sym2ensg)
res["scale_checks"] = chks

E = C.load_expression(genes)
genes = [g for g in genes if g in E.columns
         and g in Ws["achilles"].columns and g in Ws["score"].columns]
print(f"\n  genes present in BOTH screens and RNA : {len(genes):,}")
print(f"    of which valid (silence guard)      : "
      f"{len([g for g in genes if g in valid_genes]):,}")
print(f"  Broad lines {Ws['achilles'].shape[0]}  |  "
      f"Sanger lines {Ws['score'].shape[0]}  |  "
      f"shared {len(set(Ws['achilles'].index) & set(Ws['score'].index))}")
res["overlap"] = {
    "n_genes": len(genes),
    "n_valid_genes": len([g for g in genes if g in valid_genes]),
    "broad_lines": int(Ws["achilles"].shape[0]),
    "sanger_lines": int(Ws["score"].shape[0]),
    "shared_lines": int(len(set(Ws["achilles"].index) & set(Ws["score"].index))),
}

# ------------------------------------------------------------ per-screen gate
C.banner("6A -- THE GATE IN EACH SCREEN,  denominator corrected, valid genes only")
print("  Every line in the denominator was actually screened in THAT screen.")
print()

rows = []
for g in genes:
    e = E[g].dropna()
    z = C.vdw(e)
    rec = {"ensg_id": g, "is_valid": g in valid_genes,
           "symbol": symbol_of.get(g, "?")}
    for which in ("achilles", "score"):
        W = Ws[which]
        idx = [m for m in z.index if m in W.index]
        if len(idx) < C.MIN_LINES:
            continue
        s = z.reindex(idx).values
        ess = W.loc[idx, g].values
        ok = ~np.isnan(ess)
        s, ess, idxk = s[ok], ess[ok], np.array(idx)[ok]
        p = ess <= C.DEP_THRESHOLD
        if len(s) < C.MIN_LINES or p.sum() < C.MIN_POS or (~p).sum() < C.MIN_NEG:
            continue
        base = p.mean()
        q, keep = C.decile_index(s, 10, "ordinal")
        qv, pv = q[keep], p[keep]
        tag = "broad" if which == "achilles" else "sanger"
        rec[f"{tag}_base"] = float(base)
        rec[f"{tag}_d1"] = float(pv[qv == 0].mean() / base)
        rec[f"{tag}_pauc"] = C.partial_auc(s, p, *C.GATE_REGION)
        rec[f"{tag}_n"] = int(len(s))
        # detection-intersection: bottom decile AND the rule says EXPRESSED
        expr_k = e.reindex(idxk).values
        sel = (qv == 0) & (expr_k[keep] > C.EXPRESSED_MIN)
        if sel.sum() >= 10:
            rec[f"{tag}_d1_expressed_ratio"] = float(pv[sel].mean() / base)
            rec[f"{tag}_d1_expressed_n"] = int(sel.sum())
    rows.append(rec)

R = pd.DataFrame(rows)
R.to_parquet(C.OUT / "06_screen_replication_per_gene.parquet", index=False)

print(f"  {'quantity':<30} {'BROAD':>12} {'SANGER':>12} {'paired p':>12} "
      f"{'n pairs':>9}")
res["per_screen"] = {}
for col, lab in [("d1", "decile-1 ratio"),
                 ("pauc", "gate-region pAUC"),
                 ("base", "base rate"),
                 ("d1_expressed_ratio", "decile-1 & EXPRESSED ratio")]:
    for sub, nm in [(R, "all"), (R[R.is_valid], "valid")]:
        a = f"broad_{col}"
        b = f"sanger_{col}"
        if a not in sub or b not in sub:
            continue
        both = sub.dropna(subset=[a, b])
        if len(both) < 6:
            continue
        p = wilcoxon(both[a] - both[b]).pvalue
        if nm == "valid":
            print(f"  {lab:<30} {both[a].median():>12.4f} "
                  f"{both[b].median():>12.4f} {p:>12.3g} {len(both):>9,}")
        res["per_screen"].setdefault(nm, {})[col] = {
            "broad": float(both[a].median()), "sanger": float(both[b].median()),
            "paired_p": float(p), "n_pairs": int(len(both))}
print("  (rows above are VALID genes; the 'all genes' figures are in the JSON)")

# ------------------------------------------------------------ prevalence bands
C.banner("6A2 -- A SELECTION EFFECT THAT MUST BE STATED, NOT BURIED")
print("  Sanger's panel is 317 lines against Broad's 767. Requiring >=10 dependent")
print("  AND >=10 non-dependent lines in BOTH screens therefore selects genes that")
print("  are essential in a large share of lines. The paired set's median base")
print("  rate is far above the panel's.")
print()
pair = R.dropna(subset=["broad_base", "sanger_base"])
print(f"  paired genes                    : {len(pair):,}")
print(f"  their median base rate (Broad)  : {pair.broad_base.median():.4f}")
print(f"  whole-panel median base rate    : "
      f"{pd.read_parquet(C.OUT / 'genes.parquet').base_rate.median():.4f}")
print()
print("  A gate cannot act on a gene that is essential in 44% of lines. The")
print("  replication must therefore be read BY PREVALENCE BAND, and the band")
print("  where the gate is claimed to operate is the low one.")
print()
pair = pair.copy()
pair["band"] = pd.cut(pair.broad_base, [0, .02, .05, .10, .25, 1.01],
                      labels=["<2%", "2-5%", "5-10%", "10-25%", ">25%"])
print(f"  {'band':>8} {'genes':>7} {'BROAD d1':>10} {'SANGER d1':>11} "
      f"{'BROAD pAUC':>12} {'SANGER pAUC':>13}")
res["by_prevalence"] = {}
for b in pair.band.cat.categories:
    s = pair[pair.band == b]
    if len(s) < 10:
        continue
    print(f"  {str(b):>8} {len(s):>7,} {s.broad_d1.median():>10.4f} "
          f"{s.sanger_d1.median():>11.4f} {s.broad_pauc.median():>12.4f} "
          f"{s.sanger_pauc.median():>13.4f}")
    res["by_prevalence"][str(b)] = {
        "n": int(len(s)),
        "broad_d1": float(s.broad_d1.median()),
        "sanger_d1": float(s.sanger_d1.median()),
        "broad_pauc": float(s.broad_pauc.median()),
        "sanger_pauc": float(s.sanger_pauc.median())}
res["paired_set_selection"] = {
    "n_paired": int(len(pair)),
    "median_base_rate_broad": float(pair.broad_base.median()),
    "note": ("Sanger's 317-line panel plus MIN_POS/MIN_NEG selects high-prevalence "
             "genes; the paired comparison is therefore NOT representative of the "
             "panel and must be read by prevalence band.")}

# ------------------------------------------------------------ agreement
C.banner("6B -- DO THE TWO SCREENS AGREE GENE BY GENE?")
print("  A difference test that fails to reject is weak evidence of agreement.")
print("  The per-gene correlation is the direct question.")
print()
res["agreement"] = {}
for col, lab in [("d1", "decile-1 ratio"), ("pauc", "gate-region pAUC")]:
    for sub, nm in [(R, "all"), (R[R.is_valid], "valid")]:
        both = sub.dropna(subset=[f"broad_{col}", f"sanger_{col}"])
        if len(both) < 10:
            continue
        rho = both[f"broad_{col}"].corr(both[f"sanger_{col}"], method="spearman")
        agree = float(((both[f"broad_{col}"] < 1) ==
                       (both[f"sanger_{col}"] < 1)).mean()) if col == "d1" else \
                float(((both[f"broad_{col}"] > .5) ==
                       (both[f"sanger_{col}"] > .5)).mean())
        if nm == "valid":
            print(f"  {lab:<26} Spearman rho {rho:>7.4f}   "
                  f"same side of the null in {100*agree:.1f}% of genes   "
                  f"n {len(both):,}")
        res["agreement"].setdefault(nm, {})[col] = {
            "spearman": float(rho), "sign_agreement": agree,
            "n": int(len(both))}

# ------------------------------------------------------------ tautology
C.banner("6C -- 'ISN'T THIS JUST SAYING YOU CAN'T KNOCK OUT A SILENT GENE?'")
print("  The bottom decile INTERSECTED with lines the detection rule calls")
print(f"  EXPRESSED (log2TPM+1 > {C.EXPRESSED_MIN:.1f}). If depletion survives")
print("  there, the percentile is not merely an encoded detection boolean.")
print()
print(f"  {'screen':<10} {'genes':>7} {'d1 ratio':>10} {'d1 & EXPRESSED':>16} "
      f"{'p':>11}")
res["tautology_test"] = {}
for tag in ("broad", "sanger"):
    sub = R[R.is_valid].dropna(subset=[f"{tag}_d1", f"{tag}_d1_expressed_ratio"])
    if len(sub) < 6:
        print(f"  {tag:<10} too few genes")
        continue
    v = sub[f"{tag}_d1_expressed_ratio"]
    p = wilcoxon(v - 1.0).pvalue
    print(f"  {tag:<10} {len(sub):>7,} {sub[f'{tag}_d1'].median():>10.4f} "
          f"{v.median():>16.4f} {p:>11.3g}")
    res["tautology_test"][tag] = {
        "n_genes": int(len(sub)), "d1_ratio": float(sub[f"{tag}_d1"].median()),
        "d1_expressed_ratio": float(v.median()), "p": float(p),
        "frac_depleted": float((v < 1).mean())}

# ------------------------------------------------------------ verdict
C.banner("VERDICT")
pv = res.get("per_screen", {}).get("valid", {})
d1 = pv.get("d1", {})
pa = pv.get("pauc", {})
if d1:
    print(f"  decile-1 ratio   Broad {d1['broad']:.4f}   Sanger {d1['sanger']:.4f}   "
          f"difference p {d1['paired_p']:.3g}")
if pa:
    print(f"  gate pAUC        Broad {pa['broad']:.4f}   Sanger {pa['sanger']:.4f}   "
          f"difference p {pa['paired_p']:.3g}")
taut = res["tautology_test"]
surv = all(v["d1_expressed_ratio"] < 1.0 for v in taut.values()) and taut
consistent = bool(pa) and abs(pa["broad"] - pa["sanger"]) < 0.05

if consistent and surv:
    verdict = "replicates_across_screens_and_is_not_a_detection_tautology"
elif consistent:
    verdict = "replicates_across_screens_but_reduces_to_detection"
else:
    verdict = "does_not_replicate_across_screens"
print(f"\n  VERDICT: {verdict}")
print()
print("  FOR THE WRITE-UP. The gate region was selected post-hoc on the")
print("  discovery data (audit 5). This comparison is the answer to that: the")
print("  region was fixed before these two screens were compared, and screen")
print("  identity was not used to choose it. State that explicitly.")
res["verdict"] = verdict

out = C.OUT / "06_screen_replication_results.json"
out.write_text(json.dumps(res, indent=2, default=float), encoding="utf-8")
print(f"\nwrote {out}")
print(f"wrote {C.OUT / '06_screen_replication_per_gene.parquet'}")
