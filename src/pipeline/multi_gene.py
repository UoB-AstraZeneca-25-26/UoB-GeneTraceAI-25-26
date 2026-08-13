"""
multi_gene.py
-------------
Cell lines that satisfy SEVERAL genes at once.

    Gene A + Gene B  ->  a short list of cell lines that suit both,
                         not every line that suits either.

THE COMBINATION RULE
    From CellLineSelector - Confidence Engine: Design Record (v1), §5:

        "Conjunctive claims (needs several facts jointly) are bounded by the
         weakest necessary link (a min-like rule / Frechet-Hoeffding lower bound)."

    So the joint score for a cell line is the MINIMUM of its per-gene scores, not
    the mean:

        joint(line) = min_g  core_score(g, line)

    A line that is 0.99 for gene A and 0.10 for gene B scores 0.10, not 0.55. The
    mean would let a strong gene carry a weak one and produce lines that suit
    neither claim jointly -- exactly the failure a conjunctive query must avoid.
    The minimum is the Frechet-Hoeffding lower bound on P(A and B): the highest
    value the conjunction can be guaranteed regardless of how the two depend on
    each other, which is the honest bound when that dependence is unknown.

    core_score is used (not stratum_rank) because it is already a WITHIN-GENE
    percentile in [0,1] and is therefore comparable across genes. stratum_rank is
    computed within (gene, n_layers) strata, so two genes' values are not on the
    same footing.

CONFIDENCE, ALSO WEAKEST-LINK
    The joint verdict is the worst per-gene verdict (§5 again: conjunctive claims
    are bounded by the weakest link). If ANY requested gene is NO_EVIDENCE or
    UNINFORMATIVE, the conjunction inherits it and is reported outright rather
    than silently dropping that gene from the query.

EFFECTIVE INDEPENDENCE (L6)
    The Design Record warns that correlated sources must not be double-counted,
    and leaves n_eff as future work. The same caution applies to genes: asking for
    two co-expressed genes is close to asking for one. This module does not
    reweight -- it MEASURES and reports the pairwise expression correlation, so a
    redundant conjunction is visible instead of looking like independent support.

WHAT IT DOES NOT DO
    It does not claim the returned lines are experimentally validated for the
    combination. It reports which lines the existing evidence ranks highly for
    every gene asked for, with the weakest link named. Role separation (C7) holds:
    nothing here writes to or alters core_score.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from evidence_state import (NO_EVIDENCE, UNINFORMATIVE, CONFLICTING, RANKED,
                            gene_verdict)

DEFAULT_TOP_N = 10
# A line must clear this per-gene percentile for EVERY gene before it is offered
# as a joint hit. Without a floor, min() still returns a ranking when every line
# is mediocre for one of the genes, and the caller cannot tell the difference.
JOINT_FLOOR = 0.50


def _gene_frame(pred, ensg_id):
    cols = [c for c in ("model_id", "core_score", "n_layers", "confidence",
                        "class", "stratum_rank") if c in pred.columns]
    g = pred.loc[pred["ensg_id"] == ensg_id, cols].copy()
    if g.empty:
        return g
    # one row per line: a line can appear once per stratum
    return (g.sort_values("core_score", ascending=False)
             .drop_duplicates("model_id"))


def gene_pair_correlation(ensg_a, ensg_b, expr=None, base="src/pipeline/outputs"):
    """
    Spearman correlation of the two genes' core_score across shared lines.

    Reported, not corrected for -- see L6. A high value means the conjunction is
    close to redundant: the second gene is adding little independent constraint.
    """
    if expr is None or ensg_a not in expr or ensg_b not in expr:
        return None
    a, b = expr[ensg_a], expr[ensg_b]
    ok = a.notna() & b.notna()
    if ok.sum() < 30:
        return None
    from scipy.stats import spearmanr
    return float(spearmanr(a[ok], b[ok]).statistic)


def co_select(gene_inputs, pred, gene_lookup, resolve_gene, dispersion=None,
              chronos_val=None, top_n=DEFAULT_TOP_N, floor=JOINT_FLOOR):
    """
    gene_inputs : list of symbols or ENSG ids
    returns     : dict with per-gene verdicts, the joint verdict, and the shortlist
    """
    if len(gene_inputs) < 2:
        return {"error": "Give at least two genes -- use 'gene <X>' for a single gene."}

    resolved, per_gene, frames = [], [], {}
    for gi in gene_inputs:
        try:
            ensg, sym = resolve_gene(gi, gene_lookup)
        except ValueError as e:
            return {"error": str(e)}
        f = _gene_frame(pred, ensg)
        frames[ensg] = f

        disp = None
        if dispersion is not None and ensg in dispersion.index:
            disp = dispersion.loc[[ensg]]
        contradicted = False
        if chronos_val is not None and "chronos_check" in chronos_val.columns:
            cv = chronos_val[chronos_val["ensg_id"] == ensg]
            if len(cv):
                contradicted = bool((cv["chronos_check"] == "inverted").iloc[0])

        gclass = f["class"].iloc[0] if len(f) and "class" in f.columns else None
        v = gene_verdict(ensg, sym, len(f), disp, contradicted, gene_class=gclass)
        resolved.append((ensg, sym))
        per_gene.append({"ensg_id": ensg, "symbol": sym, "verdict": v,
                         "n_lines": len(f)})

    # ---- weakest link governs the conjunction ---------------------------
    order = {NO_EVIDENCE: 0, UNINFORMATIVE: 1, CONFLICTING: 2, RANKED: 3}
    weakest = min(per_gene, key=lambda d: order[d["verdict"].state])
    joint_state = weakest["verdict"].state

    usable = [d for d in per_gene if d["n_lines"] > 0]
    if len(usable) < len(per_gene):
        missing = [d["symbol"] for d in per_gene if d["n_lines"] == 0]
        return {"genes": per_gene, "joint_state": NO_EVIDENCE, "weakest": weakest,
                "shortlist": pd.DataFrame(),
                "note": (f"No cell line can satisfy this combination: "
                         f"{', '.join(missing)} has no scored data at all.")}

    # ---- intersect on lines scored for EVERY gene -----------------------
    shared = None
    for ensg, _ in resolved:
        s = set(frames[ensg]["model_id"])
        shared = s if shared is None else (shared & s)
    shared = sorted(shared or [])

    if not shared:
        return {"genes": per_gene, "joint_state": NO_EVIDENCE, "weakest": weakest,
                "shortlist": pd.DataFrame(),
                "note": "No cell line is scored for all of these genes at once."}

    wide = pd.DataFrame(index=shared)
    for ensg, sym in resolved:
        f = frames[ensg].set_index("model_id")
        wide[sym] = f.loc[shared, "core_score"]
        wide[f"{sym}_layers"] = f.loc[shared, "n_layers"]

    score_cols = [sym for _, sym in resolved]
    # Frechet-Hoeffding lower bound: the conjunction is only as strong as its
    # weakest component
    wide["joint_score"] = wide[score_cols].min(axis=1)
    wide["limiting_gene"] = wide[score_cols].idxmin(axis=1)
    wide["spread"] = wide[score_cols].max(axis=1) - wide[score_cols].min(axis=1)

    passing = wide[wide["joint_score"] >= floor].copy()
    shortlist = (passing if len(passing) else wide).sort_values(
        "joint_score", ascending=False).head(top_n)

    # ---- redundancy check (L6) ------------------------------------------
    corr = None
    if len(resolved) == 2:
        a, b = resolved[0][0], resolved[1][0]
        expr = pd.DataFrame({resolved[0][0]: frames[a].set_index("model_id")["core_score"],
                             resolved[1][0]: frames[b].set_index("model_id")["core_score"]})
        corr = gene_pair_correlation(a, b, expr)

    return {
        "genes": per_gene,
        "joint_state": joint_state,
        "weakest": weakest,
        "n_shared_lines": len(shared),
        "n_passing_floor": int(len(passing)),
        "floor": floor,
        "shortlist": shortlist.reset_index().rename(columns={"index": "model_id"}),
        "score_cols": score_cols,
        "pair_correlation": corr,
        "note": None,
    }


def format_result(res, indent="  "):
    """Human-readable rendering for the CLI."""
    from evidence_state import render
    out = []
    if "error" in res:
        return res["error"]

    out.append("\nPER-GENE EVIDENCE (the conjunction cannot be stronger than the weakest)")
    for d in res["genes"]:
        v = d["verdict"]
        mark = {NO_EVIDENCE: "NO DATA", UNINFORMATIVE: "NOT DISCRIMINATING",
                CONFLICTING: "CONFLICTING", RANKED: "ok"}[v.state]
        out.append(f"{indent}{d['symbol']:<10} {mark:<20} {d['n_lines']:>6,} lines scored")

    if res["joint_state"] != RANKED:
        w = res["weakest"]
        out.append("")
        out.append(render(w["verdict"]))
        out.append(f"{indent}This is the limiting gene, so the combined answer inherits it.")

    if res.get("note"):
        out.append(f"\n{res['note']}")
        return "\n".join(out)

    corr = res.get("pair_correlation")
    if corr is not None:
        redundant = abs(corr) > 0.7
        out.append(f"\n{indent}Gene-gene score correlation: rho = {corr:+.3f}"
                   + ("   *** the two genes track each other -- asking for both adds "
                      "little beyond asking for one ***" if redundant else
                      "   (largely independent constraints)"))

    out.append(f"\n{indent}{res['n_shared_lines']:,} cell lines are scored for every gene; "
               f"{res['n_passing_floor']:,} clear {res['floor']:.2f} on all of them.")
    if res["n_passing_floor"] == 0:
        out.append(f"{indent}*** NO cell line clears the floor for every gene. The list below "
                   f"is the closest available, NOT a set of good joint matches. ***")

    sl = res["shortlist"]
    cols = res["score_cols"]
    out.append("")
    header = f"{indent}{'model_id':<14}{'joint':<9}" + "".join(f"{c[:11]:<12}" for c in cols) \
             + f"{'limiting':<11}"
    out.append(header)
    out.append(indent + "-" * (len(header) - len(indent)))
    for _, r in sl.iterrows():
        line = f"{indent}{r['model_id']:<14}{r['joint_score']:<9.4f}"
        line += "".join(f"{r[c]:<12.4f}" for c in cols)
        line += f"{r['limiting_gene']:<11}"
        out.append(line)
    return "\n".join(out)
