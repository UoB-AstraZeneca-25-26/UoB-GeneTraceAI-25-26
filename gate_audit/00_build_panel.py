"""
gate_audit/00_build_panel.py
----------------------------
STEP 0 of 5. Build the analysis panel every other script in this folder reads.
Read-only with respect to the pipeline; writes only into gate_audit/outputs/.

WHY A PANEL STEP EXISTS
-----------------------
The original diagnostics each re-loaded the expression matrix, the two protein
layers and the Chronos HDF5 from scratch, and each re-derived the layers. That
is slow, and worse, it is five chances for the five audits to disagree about
what "the gate" is. This builds the scored panel ONCE, and 01-05 answer their
questions on that fixed object. If two audits disagree, they disagree about
statistics, not about inputs.

OUTPUTS
-------
  outputs/panel.parquet       long, one row per (gene, line):
                              ensg_id, model_id, expr, zE, zP, Z, essentiality,
                              dep, lineage
  outputs/genes.parquet       one row per gene: symbol, n_lines, base_rate,
                              frac_expressed, is_valid, and the full tie
                              structure of the expression vector
  outputs/panel_meta.json     scale checks, label FPR, validity counts, the
                              exact constants used

LABEL
-----
Real DepMap Chronos (`GeneFitnessEffect_Chronos_Achilles.hdf5`, Broad), NOT the
misnamed `validation/prepared/chronos_long.parquet`, which is Project Score
scaled BF negated. The scale check in common.load_chronos() is what distinguishes
them and it aborts rather than warns.

Run:
    python gate_audit/00_build_panel.py                 # 3000 genes, default
    python gate_audit/00_build_panel.py --n-genes 6000  # slower, tighter CIs
    python gate_audit/00_build_panel.py --rna-only      # score = zE, no protein
"""
import argparse
import json

import numpy as np
import pandas as pd

import common as C

ap = argparse.ArgumentParser()
ap.add_argument("--n-genes", type=int, default=3000,
                help="gene sample size; the full panel is ~17k and is tractable "
                     "but slow. 3000 matches the prior runs, so results are "
                     "comparable to the numbers in docs/.")
ap.add_argument("--rna-only", action="store_true",
                help="score = zE only. The tie question in 01 is about the "
                     "EXPRESSION percentile, so this isolates it.")
ap.add_argument("--screen", choices=["achilles", "score"], default="achilles")
args = ap.parse_args()

meta = {"n_genes_requested": args.n_genes, "rna_only": args.rna_only,
        "screen": args.screen,
        "constants": {"DEP_THRESHOLD": C.DEP_THRESHOLD,
                      "EXPRESSED_MIN": C.EXPRESSED_MIN,
                      "SILENT_FRAC": C.SILENT_FRAC,
                      "RHO_GLOBAL": C.RHO_GLOBAL,
                      "W_E": C.W_E, "W_P": C.W_P,
                      "MIN_LINES": C.MIN_LINES, "SEED": C.SEED}}

# ------------------------------------------------------------ 1. universe
C.banner("STEP 0.1 -- gene universe and labels")
uni, valid_ensg, symbol_of, sym2ensg, u2e = C.load_gene_lookup()
print(f"  protein-coding + HGNC Approved: {len(valid_ensg):,} genes")

W, chk = C.load_chronos(args.screen, sym2ensg)
meta["chronos_scale_check"] = chk
meta["label_fpr"] = chk["label_fpr"]

ctrl = C.control_fpr(W, symbol_of)
meta["control_family_fpr"] = ctrl
if ctrl:
    for k, v in ctrl.items():
        print(f"  control FPR [{k}]: {v['dep_rate']:.4f}  ({v['n_genes']} genes)")

# ------------------------------------------------------------ 2. gene sample
expr_col = C.expression_columns()
cands = sorted(set(W.columns) & set(expr_col) & valid_ensg)
print(f"  genes with Chronos + an RNA column: {len(cands):,}")

rng = np.random.default_rng(C.SEED)
if args.n_genes < len(cands):
    genes = sorted(rng.choice(cands, size=args.n_genes, replace=False))
else:
    genes = cands
print(f"  sampled: {len(genes):,}")

# ------------------------------------------------------------ 3. layers
C.banner("STEP 0.2 -- expression, protein layers, validity")
E = C.load_expression(genes, expr_col)
genes = [g for g in genes if g in E.columns]
print(f"  RNA matrix: {E.shape[0]} lines x {E.shape[1]:,} genes")

# Validity is defined on the FULL expression panel, before any label join, so it
# reconciles with docs/TRANSCRIPTOMICS_SPEC.md F2 (5,847 of 19,173 flagged).
fe = C.frac_expressed(E)
is_valid = fe >= C.SILENT_FRAC
print(f"  frac_expressed: median {fe.median():.3f}  IQR "
      f"{fe.quantile(.25):.3f}-{fe.quantile(.75):.3f}")
print(f"  VALID (frac_expressed >= {C.SILENT_FRAC}) : {int(is_valid.sum()):,} "
      f"({100*is_valid.mean():.1f}%)")
print(f"  INVALID / UNINFORMATIVE               : {int((~is_valid).sum()):,} "
      f"({100*(~is_valid).mean():.1f}%)")
meta["validity"] = {
    "n_genes": int(len(fe)),
    "n_valid": int(is_valid.sum()),
    "n_invalid": int((~is_valid).sum()),
    "invalid_frac": float((~is_valid).mean()),
    "frac_expressed_median": float(fe.median()),
}

if args.rna_only:
    Cm = Pm = pd.DataFrame()
    ccle_by_ensg = procan_by_ensg = {}
else:
    Cm, Pm, ccle_by_ensg, procan_by_ensg = C.load_proteomics(genes, valid_ensg, u2e)
    print(f"  CCLE {Cm.shape[0]} lines x {Cm.shape[1]:,} cols | "
          f"ProCan {Pm.shape[0]} lines x {Pm.shape[1]:,} cols")

lineage = C.load_lineage()
print(f"  lineage labels: {lineage.nunique()} distinct over {len(lineage):,} lines")

# ------------------------------------------------------------ 4. score + label
C.banner("STEP 0.3 -- scoring")
rows, gene_rows = [], []
label_lines = set(W.index)

for i, g in enumerate(genes):
    if i and i % 500 == 0:
        print(f"    {i:,}/{len(genes):,}")
    e_raw = E[g].dropna()
    if not len(e_raw):
        continue
    zE = C.vdw(e_raw)

    parts = {}
    if not args.rna_only:
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
    else:
        zP = pd.Series(dtype=float)

    idx = sorted((set(zE.index) | set(zP.index)) & label_lines)
    if len(idx) < C.MIN_LINES:
        continue
    d = pd.DataFrame({"zE": zE.reindex(idx),
                      "zP": zP.reindex(idx) if len(zP) else np.nan})
    d = d[d.notna().any(axis=1)]
    if len(d) < C.MIN_LINES:
        continue
    d["Z"] = C.stouffer(d.zE, d.zP) if len(zP) else d.zE
    d["expr"] = e_raw.reindex(d.index)
    ess = W.loc[d.index, g]
    d["essentiality"] = ess.values
    d = d[d.essentiality.notna()]
    if len(d) < C.MIN_LINES:
        continue
    d["dep"] = d.essentiality <= C.DEP_THRESHOLD
    if int(d.dep.sum()) < C.MIN_POS or int((~d.dep).sum()) < C.MIN_NEG:
        continue

    d = d.reset_index().rename(columns={"index": "model_id"})
    d["ensg_id"] = g
    d["lineage"] = d.model_id.map(lineage).fillna("unknown")
    rows.append(d[["ensg_id", "model_id", "lineage", "expr",
                   "zE", "zP", "Z", "essentiality", "dep"]])

    ts = C.tie_stats(d.expr.values)
    gene_rows.append({
        "ensg_id": g, "symbol": symbol_of.get(g, "?"),
        "n_lines": int(len(d)), "n_pos": int(d.dep.sum()),
        "base_rate": float(d.dep.mean()),
        "frac_expressed": float(fe[g]), "is_valid": bool(is_valid[g]),
        "has_protein": bool(len(zP) > 0),
        **{f"tie_{k}": v for k, v in ts.items()},
    })

panel = pd.concat(rows, ignore_index=True)
G = pd.DataFrame(gene_rows)

panel.to_parquet(C.PANEL, index=False)
G.to_parquet(C.OUT / "genes.parquet", index=False)

C.banner("STEP 0.4 -- panel built")
print(f"  genes scored          : {len(G):,}")
print(f"    valid               : {int(G.is_valid.sum()):,}")
print(f"    invalid             : {int((~G.is_valid).sum()):,}")
print(f"  (gene, line) rows     : {len(panel):,}")
print(f"  median lines per gene : {G.n_lines.median():.0f}")
print(f"  median prevalence     : {G.base_rate.median():.4f}")
print(f"  lineages represented  : {panel.lineage.nunique()}")

meta["panel"] = {
    "n_genes_scored": int(len(G)),
    "n_valid_scored": int(G.is_valid.sum()),
    "n_invalid_scored": int((~G.is_valid).sum()),
    "n_rows": int(len(panel)),
    "median_lines_per_gene": float(G.n_lines.median()),
    "median_base_rate": float(G.base_rate.median()),
    "n_lineages": int(panel.lineage.nunique()),
}
C.PANEL_META.write_text(json.dumps(meta, indent=2, default=float), encoding="utf-8")
print(f"\nwrote {C.PANEL}")
print(f"wrote {C.OUT / 'genes.parquet'}")
print(f"wrote {C.PANEL_META}")
