"""
measure_fix1_abstention_cost.py
-------------------------------
What Fix 1 costs, measured before adoption.

Fix 1: in `altered(G)`, when the line is sequenced, no driver call exists, at
least one MODERATE-impact row exists, and G is on the pre-committed uncertainty
list -> return UNKNOWN instead of FALSE.

Three questions, as specified:
  1. How many current Matches become Possible, across ~100 sampled two-term
     queries? Reported twice -- for a RANDOM modifier B (the everyday cost) and
     for B drawn from the 22-gene uncertainty list (the worst case, since the
     rule cannot fire anywhere else).
  2. Do PC9 and HCC827 move to UNKNOWN? They should.
  3. Do the mutation rows that currently pass stay passing? They should -- those
     are HIGH-impact calls and the rule does not touch them.

    python src/pipeline/measure_fix1_abstention_cost.py
"""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(Path(__file__).resolve().parent))
import multi_gene_kleene as K  # noqa: E402

OUT = ROOT / "src" / "pipeline" / "outputs"
RNG = np.random.default_rng(42)
res = {}

print("=" * 78)
print("FIX 1 -- ABSTENTION COST, MEASURED BEFORE ADOPTION")
print("=" * 78)

E = K.expression()
fe = K.gene_validity()
sequenced, alt, moderate = K.alterations()
unc = K.uncertainty_genes()
print(f"  uncertainty list: {len(unc)} genes (pre-committed)")

valid = [g for g in E.columns if fe.get(g, 0) >= K.SILENT_FRAC]
alt_counts = {g: len(s) for g, s in alt.items()}
varying = [g for g, n in alt_counts.items() if n >= 10 and g in E.columns]
unc_varying = [g for g in unc if g in varying]
print(f"  of those, varying (>=10 altered lines): {len(unc_varying)}")


def n_unknown_added(b_ensg):
    """Lines that flip FALSE -> UNKNOWN for NOT altered(B)."""
    if b_ensg not in unc:
        return 0
    mod = moderate.get(b_ensg, set())
    altered = alt.get(b_ensg, set())
    return len((mod & sequenced) - altered)


def run_pair(a, b, use_fix):
    saved = K._cache.get("unc")
    K._cache["unc"] = set(unc) if use_fix else set()
    try:
        terms = K.parse_terms(expressed=[a], not_altered=[b])
        r = K.query(terms, k_alternatives=0)
        return r["counts"]
    finally:
        K._cache["unc"] = saved


# ---------------------------------------------------------------- 1
print("\n1. MATCHES LOST ACROSS SAMPLED TWO-TERM QUERIES")
print("   expressed(A) AND NOT altered(B), 100 samples per arm")
print()
res["queries"] = {}
for arm, pool in (("random B", varying), ("B on the uncertainty list", unc_varying)):
    if not pool:
        continue
    rows = []
    n = 100 if arm == "random B" else min(100, len(pool) * 5)
    for _ in range(n):
        a = valid[RNG.integers(len(valid))]
        b = pool[RNG.integers(len(pool))]
        before = run_pair(a, b, use_fix=False)
        after = run_pair(a, b, use_fix=True)
        rows.append({"A": a, "B": b,
                     "matches_before": before["matches"],
                     "matches_after": after["matches"],
                     "possible_before": before["possible"],
                     "possible_after": after["possible"],
                     "lost": before["matches"] - after["matches"]})
    D = pd.DataFrame(rows)
    D["lost_frac"] = D.lost / D.matches_before.replace(0, np.nan)
    changed = int((D.lost > 0).sum())
    print(f"  [{arm}]  n = {len(D)}")
    print(f"    queries where any Match moved  : {changed} ({changed/len(D):.0%})")
    print(f"    Matches lost, median           : {D.lost.median():.0f}")
    print(f"    Matches lost, mean             : {D.lost.mean():.1f}")
    print(f"    as a share of Matches, median  : {D.lost_frac.median():.2%}")
    print(f"    as a share of Matches, max     : {D.lost_frac.max():.2%}")
    res["queries"][arm] = {
        "n": len(D), "queries_changed": changed,
        "frac_queries_changed": changed / len(D),
        "lost_median": float(D.lost.median()),
        "lost_mean": float(D.lost.mean()),
        "lost_frac_median": float(D.lost_frac.median(skipna=True)),
        "lost_frac_max": float(D.lost_frac.max(skipna=True))}
    print()

# how many lines the rule can touch at all, per uncertainty gene
touch = {g: n_unknown_added(g) for g in unc}
tv = np.array([v for v in touch.values()])
print(f"  lines the rule can move, per uncertainty gene: "
      f"median {np.median(tv):.0f}, max {tv.max()}, total {tv.sum()}")
print(f"  as a share of the 1,424 sequenced lines: median "
      f"{np.median(tv)/1424:.2%}, max {tv.max()/1424:.2%}")
res["per_gene_reach"] = {"median": float(np.median(tv)), "max": int(tv.max()),
                         "total": int(tv.sum()),
                         "detail": {g: int(v) for g, v in touch.items()}}

# ---------------------------------------------------------------- 2 + 3
print("\n2 + 3. GROUND TRUTH -- do the right rows move, and only those?")
sys.path.insert(0, str(ROOT / "src" / "pipeline"))
import importlib.util as ilu
spec = ilu.spec_from_file_location("gt", ROOT / "src" / "pipeline" /
                                   "test_kleene_ground_truth.py")
print("   (re-running test_kleene_ground_truth.py)")
print()

gt_json = OUT / "kleene_ground_truth_results.json"
before_rows = json.loads(gt_json.read_text(encoding="utf-8"))["rows"] \
    if gt_json.exists() else []
before = {(r["gene"], r["line"]): r.get("observed") for r in before_rows
          if r.get("tested_term") == "altered"}

mod_gt = ilu.module_from_spec(spec)
import io, contextlib
buf = io.StringIO()
with contextlib.redirect_stdout(buf):
    spec.loader.exec_module(mod_gt)
after_rows = json.loads(gt_json.read_text(encoding="utf-8"))["rows"]
after = {(r["gene"], r["line"]): r for r in after_rows
         if r.get("tested_term") == "altered"}

moved = [(k, before.get(k), v["observed"]) for k, v in after.items()
         if before.get(k) is not None and before[k] != v["observed"]]
print(f"  altered() ground-truth rows: {len(after)}")
print(f"  rows that changed state    : {len(moved)}")
for (g, l), b_, a_ in moved:
    exp = after[(g, l)]["expected"]
    print(f"    {g:<8} {l:<11} {b_} -> {a_}   (expected {exp})")
res["ground_truth_moves"] = [{"gene": g, "line": l, "before": b_, "after": a_}
                             for (g, l), b_, a_ in moved]

still_pass = sum(1 for v in after.values() if v["status"] == "PASS")
still_fail = sum(1 for v in after.values() if v["status"] == "FAIL")
now_unknown = sum(1 for v in after.values() if v["status"] == "UNKNOWN")
print(f"\n  altered() rows now: PASS {still_pass}  FAIL {still_fail}  "
      f"UNKNOWN {now_unknown}")
res["ground_truth_after"] = {"pass": still_pass, "fail": still_fail,
                             "unknown": now_unknown}

# ---------------------------------------------------------------- verdict
print("\n" + "=" * 78)
print("VERDICT")
print("=" * 78)
worst = res["queries"].get("B on the uncertainty list", {})
rand = res["queries"].get("random B", {})
print(f"  everyday cost (random B)      : "
      f"{rand.get('frac_queries_changed', 0):.0%} of queries affected, "
      f"median {rand.get('lost_frac_median', 0):.2%} of Matches")
print(f"  worst case (B on the list)    : "
      f"{worst.get('frac_queries_changed', 0):.0%} of queries affected, "
      f"median {worst.get('lost_frac_median', 0):.2%} of Matches")
small = (worst.get("lost_frac_median", 1) or 0) < 0.05
if small:
    print("\n  ADOPT. The abstention cost is small, and it buys the guarantee the")
    print("  layer sells: that a Match is real. Trading a few Matches for")
    print("  Possibles to keep that guarantee is the trade the layer exists for.")
else:
    print("\n  The cost is LARGE. That is itself a finding about how much of")
    print("  altered() rests on VEP severity ranking, and it belongs in the spec")
    print("  whichever way the decision goes.")
res["verdict"] = "adopt_small_cost" if small else "large_cost_report_either_way"

p = OUT / "fix1_abstention_cost.json"
p.write_text(json.dumps(res, indent=2, default=float), encoding="utf-8")
print(f"\nwrote {p}")
