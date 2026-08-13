"""
cell_similarity_bridge.py
-------------------------
Makes the validated cell-line similarity graph reachable from the pipeline's own
query layer, alongside `rank_cell_lines()` / `explain_pair()`.

Additive and read-only: it does not modify `core_score`, does not re-sort any
ranking, and writes nothing. It attaches two things to an existing result --
alternatives per line, and honest provenance for both claims in the answer.

WHAT IS BEING WIRED IN
----------------------
`cell_similarity/` builds a line x line RNA-expression similarity graph. Its one
validated claim:

    lines that are close on the RNA axis share DRUG-RESPONSE profiles
    +0.186 [+0.149, +0.227] above a same-tissue baseline
    (GDSC2, 649 lines x 150 compounds, cluster bootstrap over cell lines)

It survives dropping haematological lines (+0.190) and holds within lineage
(+0.252), so it is not a tissue shortcut. Split-half feature reliability of the
RNA axis is 83.7%. Metabolomics and miRNA are NOT wired in as rankings: neither
beats a plain tissue lookup (+0.008 and +0.019, both n.s.).

THE PROVENANCE PROBLEM THIS MODULE HAS TO HANDLE
-------------------------------------------------
`rank_cell_lines()` sorts by `core_score`, which is the abundance prior. The
`gate_audit/` work established that the *dependency* claim attached to that
prior does not survive:

  * 98% of the published decile-1 depletion is reproduced by a placebo with the
    dependency labels randomised -- it was an artefact of counting
    never-screened cell lines as confirmed non-dependencies (`gate_audit/03`,
    confirmed at DepMap 24Q4 in `gate_audit/08`)
  * what remains is CRISPR-specific: flat on drug sensitivity, and OPPOSITE on
    RNAi at matched prevalence (`gate_audit/09`)

So this module must not quietly attach a validated claim to an invalidated one
and let the pair read as equally supported. Every result it returns carries a
`claims` block naming what each part rests on. `core_score` remains an
*abundance* ordering, which is what it always measured; it is the dependency
interpretation that is withdrawn, and the caveat says exactly that rather than
implying the column is meaningless.

USE
---
    from explain_pair import load_data
    from rank_cell_lines import rank_cell_lines
    from cell_similarity_bridge import attach_alternatives, similar_lines_for

    pred, flags, variants, harm, gl, cv = load_data()
    res = rank_cell_lines("ensg00000146648", pred, cv)
    res = attach_alternatives(res, k=5, top_n=20)

    similar_lines_for("ach-000552", k=10)

CLI
    python src/pipeline/cell_similarity_bridge.py --gene EGFR --top-n 10
    python src/pipeline/cell_similarity_bridge.py --line ACH-000552
"""
from __future__ import annotations

import argparse
import importlib.util as _ilu
import json
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
SIM_DIR = ROOT / "cell_similarity"
SIM_OUT = SIM_DIR / "outputs"

# The validated claim, quoted with every answer so it travels with the result
# instead of living in a README. Observed value first -- the correction is
# secondary and no claim here depends on it.
SIMILARITY_CLAIM = {
    "what": "lines whose RNA profile is closest, which share drug-response behaviour",
    "readout": "GDSC2 drug response, 150 compounds",
    "delta_observed": 0.1856,
    "ci": [0.1492, 0.2270],
    "baseline": "same-tissue lookup",
    "n_lines": 649,
    "resampling": "cluster bootstrap over cell lines",
    "robustness": {
        "haematological lines dropped": 0.1901,
        "within-lineage pairs only": 0.2519,
    },
    "axis_reliability_split_half": 0.837,
    "delta_disattenuated": 0.2623,
    "disattenuation_note": (
        "secondary; readout reliability is 0.707, high enough that the observed "
        "value carries the claim on its own"),
    "source": "cell_similarity/08_drug_arm.py",
}

CORE_SCORE_CAVEAT = {
    "what_it_is": "an abundance / expressibility prior (RNA + protein percentile)",
    "what_it_is_not": (
        "a dependency estimate. The claim that low abundance predicts "
        "non-dependency does not survive: a placebo with randomised labels "
        "reproduces 98% of the published effect, which came from scoring "
        "never-screened cell lines as confirmed non-dependencies"),
    "also": (
        "the residual effect is CRISPR-specific -- flat on drug sensitivity and "
        "opposite on RNAi at matched prevalence"),
    "source": "gate_audit/03_denominator.py, gate_audit/09_perturbation_types.py",
}

# Only RNA is wired in as a ranking. See module docstring.
RANKABLE_AXIS = "rna"
CONTEXT_AXES = ("metab", "mirna")
LOW_SIM_WARN = 0.40

_cache: dict = {}


def _sim_common():
    """cell_similarity/common.py, loaded by explicit path -- src/pipeline also
    has modules that would shadow a plain import, and cell_similarity/common.py
    itself loads gate_audit/common.py the same way."""
    if "mod" not in _cache:
        spec = _ilu.spec_from_file_location("cell_similarity_common",
                                            SIM_DIR / "common.py")
        mod = _ilu.module_from_spec(spec)
        sys.modules["cell_similarity_common"] = mod
        spec.loader.exec_module(mod)
        _cache["mod"] = mod
    return _cache["mod"]


def graph_available() -> bool:
    return (SIM_OUT / f"{RANKABLE_AXIS}_similarity.parquet").exists()


def _sim(axis: str = RANKABLE_AXIS):
    key = f"sim_{axis}"
    if key not in _cache:
        p = SIM_OUT / f"{axis}_similarity.parquet"
        _cache[key] = pd.read_parquet(p) if p.exists() else None
    return _cache[key]


def _lineage():
    if "lin" not in _cache:
        _cache["lin"] = _sim_common().load_lineage()
    return _cache["lin"]


def similar_lines_for(model_id: str, k: int = 5, axis: str = RANKABLE_AXIS) -> dict:
    """Alternatives for one line. Absence is reported as absence -- a line with
    no data on this axis returns `covered: False` and the size of the axis, never
    an empty list that would read as 'no similar lines exist'."""
    mid = str(model_id).strip().lower()
    S = _sim(axis)
    if S is None:
        return {"model_id": mid, "covered": False,
                "reason": "similarity graph not built; run cell_similarity/01_build_axes.py"}
    lin = _lineage()
    if mid not in S.index:
        return {"model_id": mid, "covered": False,
                "reason": f"this line has no {axis} data",
                "lines_on_axis": int(S.shape[0])}
    v = S.loc[mid].dropna().sort_values(ascending=False).head(k)
    self_lin = str(lin.lineage.get(mid)) if mid in lin.index else None
    rows = [{"rank": r, "model_id": nb, "similarity": float(s),
             "lineage": (str(lin.lineage.get(nb)) if nb in lin.index else None),
             "same_lineage": (str(lin.lineage.get(nb)) == self_lin
                              if nb in lin.index else None)}
            for r, (nb, s) in enumerate(v.items(), start=1)]
    return {"model_id": mid, "covered": True, "axis": axis,
            "lineage": self_lin, "lines_on_axis": int(S.shape[0]),
            "alternatives": rows,
            "top1_similarity": float(v.iloc[0]) if len(v) else None,
            "low_confidence": bool(len(v) and v.iloc[0] < LOW_SIM_WARN)}


def attach_alternatives(result: dict, k: int = 5, top_n: int = 20) -> dict:
    """Attach drug-response alternatives to a `rank_cell_lines()` result.

    Does NOT re-sort. `sorted_position` and `core_score` are untouched; each of
    the first `top_n` rows simply gains an `alternatives` list, and the result
    gains a `claims` block so the two differently-supported statements in the
    answer cannot be mistaken for one another.
    """
    if "ranking" not in result:
        return result
    out = dict(result)
    out["claims"] = {
        "ranking_basis": {
            "column": "core_score", "status": "ABUNDANCE ONLY",
            **CORE_SCORE_CAVEAT},
        "alternatives": {"status": "VALIDATED", **SIMILARITY_CLAIM},
    }
    if not graph_available():
        out["alternatives_available"] = False
        out["alternatives_note"] = (
            "similarity graph not built; run cell_similarity/01_build_axes.py")
        return out

    ranking = [dict(r) for r in out["ranking"]]
    covered = 0
    for row in ranking[:top_n]:
        alt = similar_lines_for(row["model_id"], k=k)
        row["alternatives"] = alt.get("alternatives", [])
        row["alternatives_covered"] = alt["covered"]
        if not alt["covered"]:
            row["alternatives_reason"] = alt.get("reason")
        else:
            covered += 1
            row["alternatives_low_confidence"] = alt["low_confidence"]
    out["ranking"] = ranking
    out["alternatives_available"] = True
    out["alternatives_coverage"] = {
        "rows_annotated": min(top_n, len(ranking)),
        "rows_covered": covered,
        "rows_not_covered": min(top_n, len(ranking)) - covered,
        "context_axes_not_wired": list(CONTEXT_AXES),
        "why": ("metabolomics and miRNA do not beat a same-tissue lookup "
                "(+0.008 and +0.019, both n.s.), so they are not offered as "
                "rankings here"),
    }
    return out


# =============================================================== CLI
def _print(res, top_n, k):
    print("=" * 78)
    print(f"GENE {res.get('gene')}   class={res.get('class')}   "
          f"candidates={res.get('total_candidates')}")
    print("=" * 78)
    c = res.get("claims", {})
    rb, al = c.get("ranking_basis", {}), c.get("alternatives", {})
    print(f"  RANKING BASIS : core_score -- {rb.get('status')}")
    print(f"    {rb.get('what_it_is')}")
    print(f"    NOT a dependency estimate ({rb.get('source')})")
    print(f"  ALTERNATIVES  : {al.get('status')}")
    print(f"    +{al.get('delta_observed'):.3f} "
          f"[{al['ci'][0]:.3f}, {al['ci'][1]:.3f}] over a {al.get('baseline')} "
          f"on {al.get('readout')}")
    cov = res.get("alternatives_coverage")
    if cov:
        print(f"    coverage: {cov['rows_covered']}/{cov['rows_annotated']} "
              f"of the top rows have similarity data")
    print()
    for row in res["ranking"][:top_n]:
        print(f"  #{row['sorted_position']:<4} {row['model_id']:<14} "
              f"abundance={row['core_score']:.3f}  conf={row['confidence']:<8} "
              f"{str(row.get('diseases'))[:34]}")
        if not row.get("alternatives_covered", False):
            print(f"        alternatives: NOT COVERED -- "
                  f"{row.get('alternatives_reason')}")
            continue
        flag = "  << low confidence" if row.get("alternatives_low_confidence") else ""
        alts = ", ".join(f"{a['model_id']}({a['similarity']:.2f})"
                         for a in row["alternatives"][:k])
        print(f"        alternatives: {alts}{flag}")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--gene", help="symbol or ENSG")
    ap.add_argument("--line", help="model_id, e.g. ACH-000552")
    ap.add_argument("--top-n", type=int, default=10)
    ap.add_argument("--k", type=int, default=5)
    ap.add_argument("--json", action="store_true")
    a = ap.parse_args()
    if not a.gene and not a.line:
        ap.error("give --gene and/or --line")

    if a.line and not a.gene:
        out = similar_lines_for(a.line, k=a.k)
        out["claim"] = SIMILARITY_CLAIM
        print(json.dumps(out, indent=2, default=float))
        raise SystemExit

    sys.path.insert(0, str(Path(__file__).resolve().parent))
    from explain_pair import load_data, resolve_gene          # noqa: E402
    from rank_cell_lines import rank_cell_lines               # noqa: E402

    pred, flags, variants, harm, gene_lookup, chronos_val = load_data()
    ensg = a.gene if a.gene.lower().startswith("ensg") else None
    if ensg is None:
        r = resolve_gene(a.gene, gene_lookup)
        ensg = r[0] if isinstance(r, (tuple, list)) else r
        if isinstance(ensg, dict):
            ensg = ensg.get("ensg_id")
    res = rank_cell_lines(str(ensg).lower(), pred, chronos_val)
    if "error" in res:
        raise SystemExit(res["error"])
    res = attach_alternatives(res, k=a.k, top_n=a.top_n)
    if a.json:
        print(json.dumps(res, indent=2, default=float))
    else:
        _print(res, a.top_n, a.k)
