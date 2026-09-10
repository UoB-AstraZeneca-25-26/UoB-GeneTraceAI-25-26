"""
cell_similarity/04_query.py
---------------------------
THE QUERY PATH.  One route: gene -> what is known -> candidate lines ->
alternatives.

    python cell_similarity/04_query.py --gene EGFR
    python cell_similarity/04_query.py --gene EGFR --line ACH-000552
    python cell_similarity/04_query.py --line ACH-000552 --k 10

WHAT THIS TOOL CLAIMS, AND WHAT IT DOES NOT
--------------------------------------------
CLAIMS: for a cell line, it returns other lines that behave similarly under
drug perturbation. That claim is validated (`08_drug_arm.py`): RNA-similarity
neighbours share drug-response profiles at **+0.186 [+0.149, +0.227]** above a
same-tissue baseline, cluster-bootstrapped over cell lines, on 649 lines and 150
compounds. It survives dropping haematological lines (+0.190) and holds within
lineage (+0.252), so it is not a haem/solid shortcut.

DOES NOT CLAIM: anything about whether knocking out a gene will kill a line.
The abundance-gate work (`gate_audit/`) established that the dependency claim
this project was originally built around does not survive; it is absent here on
purpose and must not reappear.

The gene view is therefore DESCRIPTIVE, not predictive. It reports what was
measured. "Candidate lines" means "lines where this gene is expressed", which is
an observation. Ranking them by predicted dependency would be reintroducing the
dead claim under a new name.

THE THREE AXES
--------------
    rna      VALIDATED RANKER. +0.186 against drug response over tissue.
             Split-half reliability 83.7%.
    metab    within-tissue context only. Does not beat a tissue lookup
             (+0.008, n.s. against drug response). Split-half reliability 42.9%
             -- 225 metabolites cannot pin down a stable 928-line graph.
    mirna    within-tissue context only (+0.019, n.s.). Split-half 58.8% after
             the log1p rebuild.

Only RNA is offered as a ranking. The other two are shown as context and
labelled with what they are worth.

DESIGN RULES, EACH INHERITED FROM SOMETHING MEASURED
-----------------------------------------------------
1. ABSENCE IS REPORTED AS ABSENCE. A line with no metabolomics gets
   "NOT COVERED - this line has no metabolomics data (928 lines do have it)",
   never an empty list. This is the friendly form of the mistake `gate_audit`
   found, where unmeasured lines were silently counted as confirmed negatives.
2. AXES STAY SEPARATE. Top-25 neighbour overlap between axes is 7-15%, so a
   single fused list would imply an agreement that does not exist.
3. LOW-CONFIDENCE NEIGHBOURS ARE FLAGGED. Lines whose nearest neighbour is below
   0.40 similarity have no real neighbour; `03` named those cases.
4. NO DEPENDENCY PREDICTION ANYWHERE.

Read-only. Reads `cell_similarity/outputs/` and the pipeline's reference tables.
"""
import argparse
import json

import numpy as np
import pandas as pd

import common as C

RANKABLE_AXES = ("rna",)
CONTEXT_AXES = ("metab", "mirna")
LOW_SIM_WARN = 0.40

# Headline validation, quoted in the output so the claim travels with the answer.
VALIDATION = {
    "readout": "GDSC2 drug response (150 compounds)",
    "delta": 0.1856,
    "ci": (0.1492, 0.2270),
    "n_lines": 649,
    "baseline": "same-tissue lookup",
    "reliability": 0.707,
}
AXIS_ROLE = {
    "rna": ("validated ranker",
            "+0.186 drug-response concordance over a same-tissue baseline"),
    "metab": ("within-tissue context",
              "does not beat a tissue lookup (+0.008, n.s.); reliability 42.9%"),
    "mirna": ("within-tissue context",
              "does not beat a tissue lookup (+0.019, n.s.); reliability 58.8%"),
}

_cache = {}


def _lineage():
    if "lin" not in _cache:
        _cache["lin"] = C.load_lineage()
    return _cache["lin"]


def _sim(axis):
    k = f"sim_{axis}"
    if k not in _cache:
        _cache[k] = C.load_axis_similarity(axis)
    return _cache[k]


def _rna_matrix():
    if "rna_mat" not in _cache:
        _cache["rna_mat"] = C.load_rna()
    return _cache["rna_mat"]


def _gene_lookup():
    if "gl" not in _cache:
        gl = pd.read_parquet(C.GA.ROOT / "reference" / "gene_lookup.parquet",
                             columns=["ensg_id", "hgnc_symbol", "biotype"])
        gl["ensg_id"] = gl.ensg_id.astype("string").str.split(".").str[0].str.lower()
        gl["hgnc_symbol"] = gl.hgnc_symbol.astype(str).str.upper()
        _cache["gl"] = gl
    return _cache["gl"]


def resolve_gene(q: str):
    gl = _gene_lookup()
    q = str(q).strip()
    hit = (gl[gl.ensg_id == q.split(".")[0].lower()] if q.lower().startswith("ensg")
           else gl[gl.hgnc_symbol == q.upper()])
    return (None, None) if not len(hit) else (hit.iloc[0].ensg_id,
                                              hit.iloc[0].hgnc_symbol)


# =============================================================== gene view
def describe_gene(q: str) -> dict:
    """What is measured for this gene, and what is not. Descriptive only."""
    ensg, symbol = resolve_gene(q)
    out = {"query": q, "ensg_id": ensg, "symbol": symbol}
    if ensg is None:
        return {**out, "state": "UNKNOWN_GENE",
                "message": f"'{q}' does not resolve to a gene in the lookup."}

    E = _rna_matrix()
    axes = {}
    if ensg in E.columns:
        v = E[ensg]
        axes["rna"] = {
            "gene_resolved": True,
            "lines_measured": int(v.notna().sum()),
            "lines_total": int(len(v)),
            "lines_expressed": int((v > C.GA.EXPRESSED_MIN).sum()),
            "frac_expressed": float((v > C.GA.EXPRESSED_MIN).mean()),
        }
    else:
        axes["rna"] = {"gene_resolved": False, "lines_measured": 0,
                       "note": "no RNA column for this gene in DepMap 24Q4"}
    for a in CONTEXT_AXES:
        try:
            axes[a] = {"gene_resolved": False,
                       "lines_with_axis": int(_sim(a).shape[0]),
                       "note": ("this assay has no gene axis; it contributes "
                                "cell-line similarity only")}
        except SystemExit:
            axes[a] = {"gene_resolved": False, "lines_with_axis": 0,
                       "note": "axis not built"}
    out["axes"] = axes

    r = axes["rna"]
    if not r.get("gene_resolved"):
        out["state"], out["message"] = "NO_EVIDENCE", \
            "No RNA measurement exists for this gene."
    elif r.get("frac_expressed", 0) < C.GA.SILENT_FRAC:
        out["state"] = "UNINFORMATIVE"
        out["message"] = (
            f"Expressed in only {r['frac_expressed']:.1%} of measured lines "
            f"(< {C.GA.SILENT_FRAC:.0%}). A within-gene percentile is not a "
            f"meaningful relative position here.")
    else:
        out["state"] = "OK"
        out["message"] = (f"Measured in {r['lines_measured']:,} lines, expressed "
                          f"in {r['frac_expressed']:.1%} of them.")
    return out


def candidate_lines(q: str, k: int = 10, lineage: str = None) -> dict:
    """Lines where this gene is EXPRESSED, highest first.

    This is an observation, not a prediction. It does not say the line depends
    on the gene, and it is not ordered by any dependency estimate -- that claim
    did not survive `gate_audit/` and is deliberately absent.
    """
    ensg, symbol = resolve_gene(q)
    out = {"query": q, "ensg_id": ensg, "symbol": symbol,
           "basis": "measured RNA expression, descending",
           "not_a_prediction": ("these lines express the gene; this is not a "
                                "claim that they depend on it")}
    if ensg is None:
        return {**out, "state": "UNKNOWN_GENE", "candidates": []}
    E = _rna_matrix()
    if ensg not in E.columns:
        return {**out, "state": "NO_EVIDENCE", "candidates": []}
    lin = _lineage()
    v = E[ensg].dropna()
    if lineage:
        keep = [m for m in v.index
                if str(lin.lineage.get(m, "")).lower() == lineage.lower()]
        v = v.loc[keep]
        out["filtered_to_lineage"] = lineage
    expressed = v[v > C.GA.EXPRESSED_MIN].sort_values(ascending=False)
    out["state"] = "OK" if len(expressed) else "UNINFORMATIVE"
    out["n_expressed"] = int(len(expressed))
    out["n_measured"] = int(len(v))
    out["candidates"] = [
        {"model_id": m, "log2_tpm1": float(x),
         "lineage": (str(lin.lineage.get(m)) if m in lin.index else None)}
        for m, x in expressed.head(k).items()]
    return out


# =============================================================== line view
def similar_lines(model_id: str, k: int = 10, axes=None) -> dict:
    """Lines that behave similarly to this one under drug perturbation.

    One ranked list per axis; never fused. Only the RNA list is a validated
    ranking -- the others are labelled as context.
    """
    mid = str(model_id).strip().lower()
    lin = _lineage()
    axes = axes or C.available_axes()
    out = {"model_id": mid,
           "lineage": (str(lin.lineage.get(mid)) if mid in lin.index else None),
           "validated_claim": VALIDATION, "axes": {}, "not_covered": []}

    per_axis = {}
    for axis in axes:
        S = _sim(axis)
        role, note = AXIS_ROLE.get(axis, ("unlabelled", ""))
        if mid not in S.index:
            out["not_covered"].append(axis)
            out["axes"][axis] = {
                "covered": False, "role": role,
                "reason": f"this line has no {C.AXIS_LABEL[axis]} data",
                "lines_on_axis": int(S.shape[0])}
            continue
        v = S.loc[mid].dropna().sort_values(ascending=False).head(k)
        rows = [{"rank": r, "model_id": nb, "similarity": float(s),
                 "lineage": (str(lin.lineage.get(nb)) if nb in lin.index else None)}
                for r, (nb, s) in enumerate(v.items(), start=1)]
        per_axis[axis] = set(v.index)
        out["axes"][axis] = {
            "covered": True, "role": role, "role_note": note,
            "lines_on_axis": int(S.shape[0]), "neighbours": rows,
            "same_lineage_in_topk": sum(1 for x in rows
                                        if x["lineage"] == out["lineage"]),
            "top1_similarity": float(v.iloc[0]) if len(v) else None,
            "low_confidence": bool(len(v) and v.iloc[0] < LOW_SIM_WARN)}

    if len(per_axis) > 1:
        from collections import Counter
        cnt = Counter()
        for s in per_axis.values():
            cnt.update(s)
        out["multi_axis_agreement"] = [
            {"model_id": m, "n_axes": n,
             "axes": sorted(a for a, s in per_axis.items() if m in s),
             "lineage": (str(lin.lineage.get(m)) if m in lin.index else None)}
            for m, n in cnt.most_common() if n > 1]
    return out


def explore(gene: str = None, model_id: str = None, k: int = 10) -> dict:
    """The whole path in one call: gene -> what is known -> candidates ->
    alternatives."""
    out = {}
    if gene:
        out["gene"] = describe_gene(gene)
        out["candidates"] = candidate_lines(gene, k=k)
        if model_id is None and out["candidates"].get("candidates"):
            model_id = out["candidates"]["candidates"][0]["model_id"]
            out["auto_selected_line"] = model_id
    if model_id:
        out["line"] = similar_lines(model_id, k=k)
    return out


# =============================================================== CLI
def _p_gene(g):
    C.banner(f"1. GENE  {g.get('symbol') or g['query']}  "
             f"({g.get('ensg_id') or '?'})")
    print(f"  STATE: {g['state']}   {g['message']}")
    if "axes" not in g:
        return
    r = g["axes"]["rna"]
    print()
    print("  WHAT IS MEASURED")
    if r.get("gene_resolved"):
        print(f"    RNA            {r['lines_measured']:,} of "
              f"{r['lines_total']:,} lines; expressed in {r['lines_expressed']:,}"
              f" ({r['frac_expressed']:.1%})")
    else:
        print("    RNA            NOT MEASURED for this gene")
    for a in CONTEXT_AXES:
        d = g["axes"].get(a, {})
        print(f"    {C.AXIS_LABEL[a]:<14} no gene axis; "
              f"{d.get('lines_with_axis', 0):,} lines carry it")
    print()
    print("  WHAT IS NOT")
    if r.get("gene_resolved"):
        print(f"    {r['lines_total'] - r['lines_measured']:,} lines have no RNA "
              f"value for this gene.")
        print(f"    {r['lines_measured'] - r['lines_expressed']:,} measured lines "
              f"sit below the detection floor -- measured, but silent.")
    print("    Nothing here predicts whether knocking this gene out would kill a")
    print("    line. That claim is not supported and is not offered.")


def _p_cand(c, k):
    C.banner(f"2. CANDIDATE LINES  ({c.get('symbol') or c['query']})")
    print(f"  Basis: {c['basis']}")
    print(f"  NOT a prediction: {c['not_a_prediction']}")
    if not c.get("candidates"):
        print(f"\n  STATE: {c.get('state')} -- no lines to offer.")
        return
    print(f"\n  {c['n_expressed']:,} of {c['n_measured']:,} measured lines "
          f"express this gene. Top {min(k, len(c['candidates']))}:")
    print(f"    {'line':<16} {'log2(TPM+1)':>12}  lineage")
    for r in c["candidates"][:k]:
        print(f"    {r['model_id']:<16} {r['log2_tpm1']:>12.2f}  "
              f"{r['lineage'] or 'unknown'}")


def _p_line(L, k):
    v = L["validated_claim"]
    C.banner(f"3. ALTERNATIVE LINES  {L['model_id']}   "
             f"lineage: {L['lineage'] or 'unknown'}")
    print(f"  Lines that behave similarly under drug perturbation.")
    print(f"  Validated on {v['readout']}: +{v['delta']:.3f} "
          f"[{v['ci'][0]:.3f}, {v['ci'][1]:.3f}] over a {v['baseline']},")
    print(f"  {v['n_lines']} lines, cluster-bootstrapped over cell lines.")
    for axis, d in L["axes"].items():
        print()
        if not d["covered"]:
            print(f"  [{axis}]  NOT COVERED -- {d['reason']} "
                  f"({d['lines_on_axis']:,} lines do have it)")
            continue
        flag = "   << low confidence: nearest neighbour is not near" \
            if d["low_confidence"] else ""
        print(f"  [{axis}]  {C.AXIS_LABEL[axis]} -- {d['role'].upper()}{flag}")
        if d.get("role_note"):
            print(f"        {d['role_note']}")
        print(f"    {'rank':>4} {'line':<16} {'sim':>7}  lineage")
        for r in d["neighbours"][:k]:
            mark = "*" if r["lineage"] == L["lineage"] else " "
            print(f"    {r['rank']:>4} {r['model_id']:<16} "
                  f"{r['similarity']:>7.3f} {mark} {r['lineage'] or 'unknown'}")
        print(f"    {d['same_lineage_in_topk']} of {len(d['neighbours'])} share "
              f"this line's lineage (*)")
    if L.get("multi_axis_agreement"):
        print()
        print("  FOUND BY MORE THAN ONE AXIS")
        for a in L["multi_axis_agreement"][:8]:
            print(f"    {a['model_id']:<16} {a['n_axes']} axes  "
                  f"({', '.join(a['axes'])})  {a['lineage'] or 'unknown'}")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--gene", help="symbol or ENSG")
    ap.add_argument("--line", help="model_id, e.g. ACH-000552")
    ap.add_argument("--k", type=int, default=10)
    ap.add_argument("--json", action="store_true")
    a = ap.parse_args()
    if not a.gene and not a.line:
        ap.error("give --gene and/or --line")

    payload = explore(gene=a.gene, model_id=a.line, k=a.k)
    if a.json:
        print(json.dumps(payload, indent=2, default=float))
    else:
        if "gene" in payload:
            _p_gene(payload["gene"])
            _p_cand(payload["candidates"], a.k)
            if payload.get("auto_selected_line"):
                print(f"\n  (no --line given; using the top candidate "
                      f"{payload['auto_selected_line']} below)")
        if "line" in payload:
            _p_line(payload["line"], a.k)
