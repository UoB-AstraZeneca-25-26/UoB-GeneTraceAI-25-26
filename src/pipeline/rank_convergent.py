"""
rank_convergent.py
------------------
Convergent-evidence ranking, implementing docs/RANKING_DECLARATION.md exactly.

    "how much independent evidence is there that this line is in the state the
     researcher asked for?"

No gene-role classification anywhere. Direction comes from the query.

THE RULE (declaration section 4)
---------------------------------
Per line, build the ordered evidence vector for the declared direction:

    v = (s1, ..., s5)      si in { +1 confirms, 0 unknown, -1 contradicts }

Sort by, in order:
    1. number of confirming layers, descending
    2. the vector v itself, lexicographically  (layer 1 outranks layers 2+3)
    3. number of contradicting layers, ascending
    4. still tied -> DECLARE TIED, do not impose an order

A weighted sum is explicitly rejected: equal weights with nothing to justify
them, and it hides which layers fired. Every layer is its own column.

UNKNOWN IS NOT CONTRADICTION (declaration section 5)
-----------------------------------------------------
si = 0 when the layer did not MEASURE that (gene, line); -1 only when it measured
and disagreed. Requires a `measured` flag per layer:

    CN          membership in cn_gene_complete.parquet  (908 x 18,276, 0% missing)
    mutation    line present in mutations_collapsed
    expression  non-NaN in the 24Q4 matrix
    protein     NOT AVAILABLE -- NaN conflates not-detected with not-run

So the protein layer contributes 0 always. It can raise a rank, never lower one.
Declared, not a bug.

    python src/pipeline/rank_convergent.py --gene TP53 --direction loss
    python src/pipeline/rank_convergent.py --gene BRAF --direction gain
    python src/pipeline/rank_convergent.py --acceptance-test-1
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(Path(__file__).resolve().parent))
import multi_gene_kleene as K  # noqa: E402

OUT = ROOT / "src" / "pipeline" / "outputs"

# ---- thresholds, from docs/RANKING_DECLARATION.md section 7 (CN marginal only)
CN_DEEP_DEL = 0.25      # bottom 0.121% of the CN marginal
CN_LOSS = 0.50          # bottom 1.68%
CN_AMP = 1.50           # top 5.12%
CN_NEUTRAL_HI = 1.15    # p75 -- used only for the CONTRADICTS direction
CN_P25 = 0.885

LOSS_LAYERS = ["cn_deep_deletion", "truncating_mutation", "cn_loss",
               "protein_absent", "expression_low"]
GAIN_LAYERS = ["activating_mutation", "cn_amplification", "protein_high",
               "expression_high"]

_c: dict = {}


def cn_matrix():
    if "cn" not in _c:
        d = pd.read_parquet(OUT / "cn_gene_complete.parquet")
        d.index = d.index.astype(str).str.lower()
        _c["cn"] = d
    return _c["cn"]


def high_impact():
    """(sequenced lines, ensg -> lines with a HIGH-impact (rank 3) row)."""
    if "hi" not in _c:
        m = pd.read_parquet(ROOT / "cleaned_track_data" / "mutations_collapsed.parquet",
                            columns=["ensg_id", "model_id", "max_vep_rank"])
        m["model_id"] = m.model_id.astype(str).str.lower()
        m["ensg_id"] = m.ensg_id.astype(str).str.split(".").str[0].str.lower()
        seq = set(m.model_id.unique())
        d = (m[m.max_vep_rank == 3].groupby("ensg_id")["model_id"]
             .apply(set).to_dict())
        _c["hi"] = (seq, d)
    return _c["hi"]


def evidence(gene: str, lines, direction: str) -> pd.DataFrame:
    """One row per line, one column per layer with value in {+1, 0, -1}, PLUS a
    parallel `<layer>__measured` boolean column.

    BUG FOUND AND FIXED (2026-08-10, prompted by a question about whether 0
    correctly distinguishes measured-negative from unknown). It did not, for the
    mutation layers specifically:

        r["truncating_mutation"] = (1 if m in hi_set else 0)

    never checked whether the line was sequenced at all. A never-sequenced line
    and a sequenced-and-clean line both produced 0, with no way to tell them
    apart afterwards. Confirmed on EGFR: line ach-001979 (never sequenced) and
    ach-002376 (sequenced, no HIGH-impact row) were IDENTICAL in the output.

    This is NOT the Kleene core (`multi_gene_kleene.py`), which already
    distinguishes them correctly -- ground truth verification: all 9 FAIL cases
    in the 53-fact test resolve to a definite TRUE/FALSE, never UNKNOWN: the
    core's UNKNOWN handling is sound. The bug was specific to this ranking
    module's OWN evidence encoder, built independently of the core.

    THE FIX. Every layer now carries `<layer>__measured`, computed the same way
    completeness_for() computes it. The SORT VALUE is unchanged by design: a
    measured-clean call still contributes 0, not -1, matching the CN layers
    (which already do this: "measured, neutral" is 0, not -1 -- only the
    opposite extreme is -1). That design choice stands: absence of one
    mutation class is not evidence AGAINST loss via another mechanism (RB1 in
    SAOS-2 is CN-neutral and mutation-positive; the reverse case must not be
    penalised as if it contradicted). What changes is that measured-vs-unknown
    is now TRACKED and REPORTABLE, which it was not before.
    """
    ensg, _ = K.resolve_gene(gene)
    E = K.expression()
    CN = cn_matrix()
    seq, hi = high_impact()
    seq_all, alt, _mod = K.alterations()

    cn_ok = ensg in CN.columns
    cn_col = CN[ensg] if cn_ok else None
    hi_set = hi.get(ensg, set())
    drv_set = alt.get(ensg, set())
    expr = E[ensg] if ensg in E.columns else None

    rows = []
    for m in lines:
        r = {"model_id": m}
        cn = float(cn_col[m]) if (cn_ok and m in CN.index) else np.nan
        cn_measured = np.isfinite(cn)
        e = float(expr[m]) if (expr is not None and m in expr.index) else np.nan
        e_measured = np.isfinite(e)
        sequenced = m in seq

        if direction == "loss":
            # 1 CN deep deletion
            r["cn_deep_deletion"] = (0 if not cn_measured else
                                     1 if cn < CN_DEEP_DEL else
                                     -1 if cn > CN_AMP else 0)
            r["cn_deep_deletion__measured"] = cn_measured
            # 2 truncating mutation -- absence is not evidence AGAINST loss;
            #   sort value unaffected by the fix, measured status now tracked
            r["truncating_mutation"] = (1 if m in hi_set else 0)
            r["truncating_mutation__measured"] = sequenced
            # 3 CN loss
            r["cn_loss"] = (0 if not cn_measured else
                            1 if cn < CN_LOSS else
                            -1 if cn > CN_NEUTRAL_HI else 0)
            r["cn_loss__measured"] = cn_measured
            # 4 protein -- no measured flag exists yet (declaration s.5);
            #   always unknown, always UNmeasured -- consistent by construction
            r["protein_absent"] = 0
            r["protein_absent__measured"] = False
            # 5 expression low
            if not e_measured:
                r["expression_low"] = 0
            elif e < K.EXPRESSED_MIN - K.BAND:
                r["expression_low"] = 1
            elif e > K.EXPRESSED_MIN + K.BAND:
                r["expression_low"] = -1
            else:
                r["expression_low"] = 0
            r["expression_low__measured"] = e_measured
        else:
            r["activating_mutation"] = (1 if m in drv_set else 0)
            r["activating_mutation__measured"] = sequenced
            r["cn_amplification"] = (0 if not cn_measured else
                                     1 if cn > CN_AMP else
                                     -1 if cn < CN_LOSS else 0)
            r["cn_amplification__measured"] = cn_measured
            r["protein_high"] = 0
            r["protein_high__measured"] = False
            if not e_measured:
                r["expression_high"] = 0
            elif e > K.EXPRESSED_MIN + K.BAND:
                r["expression_high"] = 1
            elif e < K.EXPRESSED_MIN - K.BAND:
                r["expression_high"] = -1
            else:
                r["expression_high"] = 0
            r["expression_high__measured"] = e_measured
        r["_cn"] = cn
        r["_log2tpm"] = e
        r["_sequenced"] = sequenced
        rows.append(r)
    return pd.DataFrame(rows).set_index("model_id")


def rank(gene: str, lines, direction: str) -> pd.DataFrame:
    """Lexicographic convergent-evidence ranking with declared tied bands."""
    layers = LOSS_LAYERS if direction == "loss" else GAIN_LAYERS
    V = evidence(gene, lines, direction)
    L = V[layers]
    V["n_confirm"] = (L == 1).sum(axis=1)
    V["n_contradict"] = (L == -1).sum(axis=1)
    # sort key: confirms desc, then the vector lexicographically desc,
    # then contradictions asc
    V = V.sort_values(["n_confirm"] + layers + ["n_contradict"],
                      ascending=[False] + [False] * len(layers) + [True])
    # tied band = identical (n_confirm, vector, n_contradict)
    key = V[["n_confirm"] + layers + ["n_contradict"]].astype(int).astype(str) \
        .agg("|".join, axis=1)
    V["_key"] = key
    band = (key != key.shift()).cumsum()
    V["band"] = band
    sizes = band.value_counts()
    V["band_size"] = band.map(sizes)
    starts = V.groupby("band").cumcount().eq(0)
    first_pos = {}
    pos = 1
    for b, n in zip(sizes.index.sort_values(), [sizes[b] for b in sizes.index.sort_values()]):
        pass
    p = 1
    for b in V.band.drop_duplicates():
        first_pos[b] = p
        p += int(sizes[b])
    V["rank_from"] = V.band.map(first_pos)
    V["rank_to"] = V.rank_from + V.band_size - 1
    return V.drop(columns=["_key"])


# =============================================================== acceptance 1
def acceptance_test_1(n_queries: int = 12, seed: int = 42) -> dict:
    """Declaration section 6, test 1: DISCRIMINATION.

    On the survivors of a representative Kleene query, fewer than 60% of lines
    may fall in the single largest tied band. Above that, the ranking is mostly
    ties and the filter is doing all the work -- which is a finding, reported as
    one, not patched away.
    """
    E = K.expression()
    fe = K.gene_validity()
    CN = cn_matrix()
    rng = np.random.default_rng(seed)
    # representative genes: valid, CN-measured, and actually varying
    pool = [g for g in E.columns
            if fe.get(g, 0) >= K.SILENT_FRAC and g in CN.columns]
    genes = list(rng.choice(pool, n_queries, replace=False))
    # plus the ground-truth genes, which are the realistic use case
    for s in ["TP53", "PTEN", "RB1", "CDKN2A", "BRAF", "EGFR", "KRAS", "ERBB2"]:
        e, _ = K.resolve_gene(s)
        if e in pool and e not in genes:
            genes.append(e)

    rows = []
    for direction in ("loss", "gain"):
        for g in genes:
            terms = K.parse_terms(**({"not_expressed": [g]} if direction == "loss"
                                     else {"expressed": [g]}))
            q = K.query(terms, k_alternatives=0)
            survivors = q["buckets"]["matches"]
            if len(survivors) < 30:
                continue
            R = rank(g, survivors, direction)
            biggest = int(R.band_size.max())
            rows.append({
                "gene": g, "direction": direction,
                "survivors": len(survivors),
                "n_bands": int(R.band.nunique()),
                "largest_band": biggest,
                "largest_band_frac": biggest / len(survivors),
                "median_band_size": float(R.band_size.median()),
            })
    D = pd.DataFrame(rows)
    res = {
        "threshold": 0.60,
        "n_queries_evaluated": int(len(D)),
        "largest_band_frac_median": float(D.largest_band_frac.median()),
        "largest_band_frac_mean": float(D.largest_band_frac.mean()),
        "frac_queries_failing": float((D.largest_band_frac >= 0.60).mean()),
        "median_bands_per_query": float(D.n_bands.median()),
        "by_direction": {
            d: {"n": int(len(s)),
                "largest_band_frac_median": float(s.largest_band_frac.median()),
                "median_bands": float(s.n_bands.median())}
            for d, s in D.groupby("direction")},
    }
    res["passes"] = bool(res["largest_band_frac_median"] < 0.60)
    return res, D


# =============================================================== CLI
if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--gene")
    ap.add_argument("--direction", choices=["loss", "gain"], default="loss")
    ap.add_argument("--top", type=int, default=12)
    ap.add_argument("--acceptance-test-1", action="store_true")
    ap.add_argument("--sort-by",
                    help="display sort column instead of the declared "
                         "evidence-band order (e.g. a layer name, _cn, "
                         "_log2tpm). Re-sorts the printed rows only; it does "
                         "not change n_confirm/band, and is labelled as a "
                         "custom view, not an evidence claim.")
    a = ap.parse_args()

    if a.acceptance_test_1:
        print("=" * 78)
        print("ACCEPTANCE TEST 1 -- DISCRIMINATION")
        print("=" * 78)
        print("  Declared before building: fewer than 60% of a query's survivors")
        print("  may fall in the single largest tied band.")
        print()
        res, D = acceptance_test_1()
        print(f"  {'gene':<18}{'dir':<6}{'survivors':>10}{'bands':>7}"
              f"{'largest band':>14}{'frac':>8}")
        for _, r in D.sort_values("largest_band_frac", ascending=False).head(14).iterrows():
            print(f"  {r['gene'][:16]:<18}{r['direction']:<6}{r['survivors']:>10,}"
                  f"{r['n_bands']:>7}{r['largest_band']:>14,}"
                  f"{r['largest_band_frac']:>8.1%}")
        print()
        print(f"  queries evaluated            : {res['n_queries_evaluated']}")
        print(f"  largest-band fraction, median: {res['largest_band_frac_median']:.1%}")
        print(f"  queries at or above 60%      : {res['frac_queries_failing']:.0%}")
        print(f"  median distinct bands/query  : {res['median_bands_per_query']:.0f}")
        for d, s in res["by_direction"].items():
            print(f"    {d:<6} median largest band {s['largest_band_frac_median']:.1%}"
                  f"   median bands {s['median_bands']:.0f}   n={s['n']}")
        print()
        print(f"  VERDICT: {'PASS' if res['passes'] else 'FAIL'}")
        if not res["passes"]:
            print("  The ranking is mostly ties. The Kleene filter is doing the")
            print("  work and the ranking room adds little. Report this as the")
            print("  finding; do not tune thresholds to force a pass.")
        p = OUT / "ranking_acceptance_test_1.json"
        p.write_text(json.dumps(res, indent=2, default=float), encoding="utf-8")
        D.to_parquet(OUT / "ranking_acceptance_test_1.parquet", index=False)
        print(f"\nwrote {p}")
        raise SystemExit

    if not a.gene:
        ap.error("--gene required unless --acceptance-test-1")
    terms = K.parse_terms(**({"not_expressed": [a.gene]} if a.direction == "loss"
                             else {"expressed": [a.gene]}))
    q = K.query(terms, k_alternatives=0)
    surv = q["buckets"]["matches"]
    print("=" * 78)
    print("  WARNING: ranked by MEASUREMENT COMPLETENESS, not biological")
    print("  evidence strength. Position 1 is NOT a claim of biological")
    print("  superiority. Completeness correlates with how well-studied a")
    print("  line is: rho=0.60 pooled (p=1.4e-134), rho=0.38 within-lineage")
    print("  (median 0.36 across 18 lineages) -- fame bias survives holding")
    print("  lineage fixed. MUTATION is the worst layer for this: rho=0.72,")
    print("  more than any other -- a mutation-confirmed line is more likely")
    print("  to simply be a heavily-sequenced one (measure_completeness.py).")
    print("  Even restricted to lines with >=1 confirming layer, the largest")
    print("  tied band is 94.3%, close to the 92.9% for the full survivor set")
    print("  -- ordering does not help even there (measure_completeness_2 Q7).")
    print()
    print("  'TIED' HERE MEANS EVIDENCE-EQUIVALENT, NOT MISSING DATA. Most")
    print("  non-confirming cells are measured negatives, not unknowns -- e.g.")
    print("  83.1% of truncating-mutation zeros are lines that WERE sequenced")
    print("  and found clean, only 16.9% were never sequenced. The tie is a")
    print("  property of the biology matching across lines, not of gaps in")
    print("  what was measured.")
    print("=" * 78)
    if not surv:
        print(f"query: {a.direction} of {a.gene}   0 lines qualify")
        raise SystemExit
    R = rank(a.gene, surv, a.direction)
    layers = LOSS_LAYERS if a.direction == "loss" else GAIN_LAYERS
    n_confirmed = int((R.n_confirm >= 1).sum())
    print(f"query: {a.direction} of {a.gene}")
    print(f"  {len(surv):,} lines qualify; {n_confirmed:,} have >=1 additional "
          f"confirming layer beyond the query itself")

    if a.sort_by:
        if a.sort_by not in R.columns:
            ap.error(f"--sort-by {a.sort_by!r} not a column; choices: "
                     f"{layers + ['_cn', '_log2tpm']}")
        print()
        print(f"  CUSTOM VIEW: sorted by {a.sort_by}, descending. This is a")
        print(f"  display convenience only -- it is NOT the evidence-based")
        print(f"  ranking above, and position here is not an evidence claim.")
        R = R.sort_values(a.sort_by, ascending=False)

    print()
    print(f"  {'rank':<10}{'line':<14}" + "".join(f"{l[:11]:>13}" for l in layers)
          + f"{'CN':>8}")
    print(f"  {'':<10}{'':<14}"
          + "".join(f"{'(+1/0m/0u/-1)':>13}" for _ in layers) + f"{'':>8}")
    for mid, r in R.head(a.top).iterrows():
        if a.sort_by:
            rk = f"#{list(R.index).index(mid) + 1}"
        else:
            rk = (f"#{int(r.rank_from)}" if r.band_size == 1
                  else f"#{int(r.rank_from)}-{int(r.rank_to)}")
        cn = f"{r._cn:.3f}" if np.isfinite(r._cn) else "--"

        def sym(layer):
            v = int(r[layer])
            if v == 1:
                return "+1"
            if v == -1:
                return "-1"
            measured = bool(r.get(f"{layer}__measured", False))
            return "0" if measured else "U"   # U = true unknown, never measured

        print(f"  {rk:<10}{mid:<14}"
              + "".join(f"{sym(l):>13}" for l in layers) + f"{cn:>8}")
    print()
    print("  legend: +1 confirms, 0 = measured, did not confirm, "
          "U = never measured (true unknown), -1 contradicts")
