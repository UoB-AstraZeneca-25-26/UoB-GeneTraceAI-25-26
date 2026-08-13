"""
multi_gene_kleene.py
--------------------
Multi-gene set membership under Kleene strong three-valued logic.

    expressed(EGFR) AND NOT altered(KRAS)

Answers with THREE buckets, never two: Matches, Possible, Excluded.

WHY THREE-VALUED, NOT BOOLEAN
------------------------------
`NOT altered(B)` carries two entirely different meanings that boolean logic
collapses into one:

    B was sequenced and carries no driver alteration   -> a real exclusion
    B was never sequenced                              -> nothing is known

Collapsing them is the absence-as-negative fault, which this project has already
found three times: `gate_audit/03` (never-screened lines scored as confirmed
non-dependencies, which manufactured 98% of a headline), `gate_audit/07A`
(unmeasured lines in the label denominator), and the coverage handling in the
similarity query path. Under `AND NOT` the failure is worse than in scoring: a
line with no data for B silently PASSES the filter and appears among confident
matches.

So UNKNOWN is a first-class truth value and it propagates:

    AND     | T  F  U        NOT
    --------|---------       ---------
    T       | T  F  U        T -> F
    F       | F  F  F        F -> T
    U       | U  F  U        U -> U

The one line that matters: `TRUE AND UNKNOWN = UNKNOWN`, so an unresolved line
can never reach Matches.

WHAT RESOLVES EACH TERM
-----------------------
expressed(G) on line L:
    TRUE     L measured for G, value > EXPRESSED_MIN + BAND
    FALSE    L measured for G, value < EXPRESSED_MIN - BAND
    UNKNOWN  L absent from the expression matrix
             OR G fails the validity guard (silent across the panel)
             OR |value - EXPRESSED_MIN| <= BAND  (inside measurement uncertainty)

    BAND is measured, not asserted: 16 DepMap models carry two independent RNA
    profiles; over 64,000 paired gene observations the SD of the paired
    difference is 0.5042 log2 units, so a single measurement has
    SD = 0.5042/sqrt(2) = 0.3565. BAND = 1.96 * 0.3565 = 0.699 log2 units, a 95%
    interval. Provenance: reference/depmap_profiles.parquet +
    data/parquet/data_clean/depmap_expr_clean.parquet. CAVEAT: n = 16 models is
    thin, and both profiles of a model share a library-prep batch, so this is a
    lower bound on true measurement error.

altered(G) on line L:
    TRUE     L sequenced AND a driver alteration is recorded for G
    FALSE    L sequenced AND no driver alteration recorded for G
    UNKNOWN  L never sequenced

    BOOLEAN ONLY. `mutations_collapsed.parquet` carries `variant_count` and
    `max_pathogenicity`, which aggregate by max and drift upward with variant
    count. A continuous alteration score would import that bias; this uses
    `any_driver` (plus in-frame fusions) as a flag and nothing else.

WHAT THIS LAYER DOES NOT DO
---------------------------
It produces no combined score. Each gene's value is a separate column, so the
ordering is transparently "sorted by a number you can see". Nothing is combined,
so there is nothing to validate and nothing to retract later.

Read-only. Writes nothing.

CLI
    python src/pipeline/multi_gene_kleene.py --expressed EGFR --not-altered KRAS
    python src/pipeline/multi_gene_kleene.py --expressed EGFR MET --not-expressed PTEN
"""
from __future__ import annotations

import argparse
import importlib.util as _ilu
import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]

# Kleene truth values
T, F, U = "TRUE", "FALSE", "UNKNOWN"

EXPRESSED_MIN = 1.0        # log2(TPM+1); the TPM >= 1 detection convention
SILENT_FRAC = 0.20         # validity guard, docs/TRANSCRIPTOMICS_SPEC.md F2
MEAS_SD = 0.3565           # measured, see module docstring
BAND = 1.96 * MEAS_SD      # 0.699 log2 units

_cache: dict = {}


# =============================================================== Kleene ops
def k_not(v: str) -> str:
    return U if v == U else (F if v == T else T)


def k_and(a: str, b: str) -> str:
    if a == F or b == F:
        return F                      # FALSE dominates -- one hard exclusion is enough
    if a == U or b == U:
        return U
    return T


def k_and_all(vals) -> str:
    out = T
    for v in vals:
        out = k_and(out, v)
    return out


# =============================================================== data
def _cs_common():
    """cell_similarity/common.py by explicit path -- src/pipeline would shadow a
    plain `import common`."""
    if "cs" not in _cache:
        spec = _ilu.spec_from_file_location("cs_common",
                                            ROOT / "cell_similarity" / "common.py")
        m = _ilu.module_from_spec(spec)
        sys.modules["cs_common"] = m
        spec.loader.exec_module(m)
        _cache["cs"] = m
    return _cache["cs"]


EXPR_CACHE = ROOT / "cell_similarity" / "outputs" / "expression_24q4.parquet"


def expression():
    """lines x ensg, log2(TPM+1), DepMap 24Q4.

    Cached to parquet on first use. The source is a 507 MB CSV whose parse
    dominates everything else here (~14 s); from parquet it is well under a
    second, which is what makes the query interactive. The boolean work itself
    was never the cost -- see `benchmark()`.
    """
    if "E" not in _cache:
        if EXPR_CACHE.exists():
            _cache["E"] = pd.read_parquet(EXPR_CACHE)
        else:
            E = _cs_common().load_rna()
            EXPR_CACHE.parent.mkdir(parents=True, exist_ok=True)
            try:
                E.to_parquet(EXPR_CACHE)
            except Exception:
                pass
            _cache["E"] = E
    return _cache["E"]


def benchmark(n_rep: int = 5):
    """A9: confirm the query is interactive, separating load from evaluation."""
    t0 = time.perf_counter()
    expression(); gene_validity(); alterations(); lineage()
    load_ms = (time.perf_counter() - t0) * 1000.0
    terms = parse_terms(expressed=["EGFR"], not_altered=["KRAS"])
    times = []
    for _ in range(n_rep):
        t1 = time.perf_counter()
        query(terms, k_alternatives=0)
        times.append((time.perf_counter() - t1) * 1000.0)
    E = expression()
    return {"cold_load_ms": load_ms,
            "warm_query_ms_median": float(np.median(times)),
            "warm_query_ms_min": float(np.min(times)),
            "matrix_shape": [int(E.shape[0]), int(E.shape[1])],
            "n_rep": n_rep}


def gene_validity():
    """frac_expressed per gene, and the validity flag."""
    if "fe" not in _cache:
        E = expression()
        fe = (E > EXPRESSED_MIN).sum(axis=0) / E.notna().sum(axis=0)
        _cache["fe"] = fe
    return _cache["fe"]


UNCERTAINTY_LIST = ROOT / "src" / "pipeline" / "outputs" / "altered_uncertainty_list.json"

# Which alteration classes altered() can see. Stated here and printed in the
# output legend, because `altered() == FALSE` reads as "not altered" when what it
# actually means is "no SNV, indel or in-frame fusion found".
ALTERED_SCOPE = {
    "covered": ["point mutations (missense, nonsense) flagged as drivers",
                "small indels flagged as drivers",
                "in-frame gene fusions"],
    "not_covered": ["copy-number amplification",
                    "homozygous deletion",
                    "translocations that juxtapose an enhancer without creating "
                    "a protein fusion (e.g. IGH-MYC)",
                    "exon-level structural deletions absent from the SNV table "
                    "(e.g. CTNNB1 exon 3)"],
    "fix_for_not_covered": ("wiring the copy-number layer would make "
                            "amplification and exon-level deletion visible; not "
                            "done here"),
}


def uncertainty_genes():
    """Pre-committed list, generated by build_altered_uncertainty_list.py from
    COSMIC CGC Tier 1. Committed BEFORE any call was inspected."""
    if "unc" not in _cache:
        if UNCERTAINTY_LIST.exists():
            _cache["unc"] = set(json.loads(
                UNCERTAINTY_LIST.read_text(encoding="utf-8"))["ensg"])
        else:
            _cache["unc"] = set()
    return _cache["unc"]


def tumour_suppressors():
    """gene_role from reference/gene_lookup.parquet, built from COSMIC CGC
    GRCh37 by build_gene_roles.py."""
    if "tsg" not in _cache:
        gl = pd.read_parquet(ROOT / "reference" / "gene_lookup.parquet",
                             columns=["ensg_id", "gene_role"])
        gl["ensg_id"] = gl.ensg_id.astype("string").str.split(".").str[0].str.lower()
        _cache["tsg"] = set(gl.loc[gl.gene_role.isin(["tsg", "both"]), "ensg_id"])
    return _cache["tsg"]


def alterations():
    """(sequenced model_ids, ensg -> altered lines, ensg -> MODERATE-only lines).

    Altered = a recorded driver mutation OR an in-frame fusion. Boolean by
    construction; see the docstring on aggregation bias.

    The third return value supports Fix 1. `any_driver` is effectively
    `max_vep_rank == 3` (100% of HIGH rows carry it, 0.5% of MODERATE), so an
    activating IN-FRAME INDEL -- scored MODERATE by VEP -- returns a confident
    FALSE. `moderate_only[g]` is the set of lines that have a MODERATE row for g
    and no driver call, i.e. evidence of something the annotation could not rank.
    """
    if "alt" not in _cache:
        mut = pd.read_parquet(ROOT / "cleaned_track_data" / "mutations_collapsed.parquet",
                              columns=["ensg_id", "model_id", "any_driver",
                                       "max_vep_rank"])
        mut["model_id"] = mut.model_id.astype(str).str.lower()
        mut["ensg_id"] = mut.ensg_id.astype(str).str.split(".").str[0].str.lower()
        sequenced = set(mut.model_id.unique())
        drv = mut.any_driver.fillna(False).astype(bool)
        d = (mut[drv].groupby("ensg_id")["model_id"].apply(set).to_dict())
        mod = (mut[(~drv) & (mut.max_vep_rank == 2)]
               .groupby("ensg_id")["model_id"].apply(set).to_dict())
        _cache["alt_moderate"] = mod
        try:
            fus = pd.read_parquet(ROOT / "cleaned_track_data" / "fusions_gene_level.parquet",
                                  columns=["ensg_id", "model_id", "any_in_frame"])
            fus["model_id"] = fus.model_id.astype(str).str.lower()
            fus["ensg_id"] = fus.ensg_id.astype(str).str.split(".").str[0].str.lower()
            for g, s in (fus[fus.any_in_frame.fillna(False).astype(bool)]
                         .groupby("ensg_id")["model_id"].apply(set).items()):
                d[g] = d.get(g, set()) | s
        except Exception:
            pass
        _cache["alt"] = (sequenced, d, _cache.get("alt_moderate", {}))
    return _cache["alt"]


def lineage():
    if "lin" not in _cache:
        _cache["lin"] = _cs_common().load_lineage()
    return _cache["lin"]


def resolve_gene(q: str):
    if "gl" not in _cache:
        gl = pd.read_parquet(ROOT / "reference" / "gene_lookup.parquet",
                             columns=["ensg_id", "hgnc_symbol"])
        gl["ensg_id"] = gl.ensg_id.astype("string").str.split(".").str[0].str.lower()
        gl["hgnc_symbol"] = gl.hgnc_symbol.astype(str).str.upper()
        _cache["gl"] = gl
    gl = _cache["gl"]
    q = str(q).strip()
    hit = (gl[gl.ensg_id == q.split(".")[0].lower()] if q.lower().startswith("ensg")
           else gl[gl.hgnc_symbol == q.upper()])
    return (None, None) if not len(hit) else (hit.iloc[0].ensg_id,
                                              hit.iloc[0].hgnc_symbol)


# =============================================================== terms
class Term:
    """One clause. `kind` in {expressed, altered}; `negate` applies Kleene NOT."""

    def __init__(self, kind: str, gene: str, negate: bool = False):
        self.kind = kind
        self.raw = gene
        self.ensg, self.symbol = resolve_gene(gene)
        self.negate = negate

    def label(self):
        s = f"{self.kind}({self.symbol or self.raw})"
        return f"NOT {s}" if self.negate else s

    # ---- evaluation
    def evaluate(self, lines):
        """Returns (states Series over `lines`, reason Series, value Series)."""
        n = len(lines)
        if self.ensg is None:
            return (pd.Series([U] * n, index=lines),
                    pd.Series(["gene not in lookup"] * n, index=lines),
                    pd.Series([np.nan] * n, index=lines))
        if self.kind == "expressed":
            return self._eval_expressed(lines)
        return self._eval_altered(lines)

    def _eval_expressed(self, lines):
        E = expression()
        fe = gene_validity()
        if self.ensg not in E.columns:
            return (pd.Series([U] * len(lines), index=lines),
                    pd.Series(["no RNA column for this gene"] * len(lines), index=lines),
                    pd.Series([np.nan] * len(lines), index=lines))
        if fe.get(self.ensg, 0.0) < SILENT_FRAC:
            return (pd.Series([U] * len(lines), index=lines),
                    pd.Series([f"gene fails the validity guard "
                               f"(expressed in {fe[self.ensg]:.1%} of lines)"]
                              * len(lines), index=lines),
                    E[self.ensg].reindex(lines))
        v = E[self.ensg].reindex(lines)
        st = pd.Series([U] * len(lines), index=lines)
        rs = pd.Series([""] * len(lines), index=lines)
        missing = v.isna()
        st[missing] = U
        rs[missing] = "line absent from the expression matrix"
        near = (~missing) & ((v - EXPRESSED_MIN).abs() <= BAND)
        st[near] = U
        rs[near] = (f"within measurement uncertainty of the threshold "
                    f"(+/-{BAND:.3f} log2)")
        hi = (~missing) & (v > EXPRESSED_MIN + BAND)
        st[hi] = T
        lo = (~missing) & (v < EXPRESSED_MIN - BAND)
        st[lo] = F
        if self.negate:
            st = st.map(k_not)
        return st, rs, v

    def _eval_altered(self, lines):
        sequenced, d, moderate = alterations()
        altered_set = d.get(self.ensg, set())
        mod_set = moderate.get(self.ensg, set())
        uncertain_gene = self.ensg in uncertainty_genes()
        st, rs, val = [], [], []
        for m in lines:
            if m not in sequenced:
                st.append(U)
                rs.append("line was never sequenced")
                val.append(np.nan)
            elif m in altered_set:
                st.append(T)
                rs.append("")
                val.append(1.0)
            elif uncertain_gene and m in mod_set:
                # FIX 1. Sequenced, no driver call, but a MODERATE-impact row
                # exists on a gene where activating in-frame indels are a known
                # mechanism. VEP cannot rank that event, so neither can this
                # layer -- and a confident FALSE here would admit the line to
                # Matches. Abstain instead.
                st.append(U)
                rs.append("only MODERATE-impact variants; this gene's activating "
                          "in-frame indels are not rankable by VEP severity")
                val.append(np.nan)
            else:
                st.append(F)
                rs.append("")
                val.append(0.0)
        st = pd.Series(st, index=lines)
        if self.negate:
            st = st.map(k_not)
        return st, pd.Series(rs, index=lines), pd.Series(val, index=lines)


def parse_terms(expressed=(), not_expressed=(), altered=(), not_altered=()):
    ts = []
    ts += [Term("expressed", g) for g in expressed or ()]
    ts += [Term("expressed", g, negate=True) for g in not_expressed or ()]
    ts += [Term("altered", g) for g in altered or ()]
    ts += [Term("altered", g, negate=True) for g in not_altered or ()]
    return ts


# =============================================================== query
def query(terms, k_alternatives: int = 3, cluster_order: str = "count") -> dict:
    """Evaluate a conjunction of terms across the panel."""
    t0 = time.perf_counter()
    E = expression()
    lines = list(E.index)
    lin = lineage()

    per_term, reasons, values = {}, {}, {}
    for t in terms:
        st, rs, v = t.evaluate(lines)
        per_term[t.label()] = st
        reasons[t.label()] = rs
        values[t.label()] = v

    ST = pd.DataFrame(per_term, index=lines)
    combined = ST.apply(lambda row: k_and_all(row.values), axis=1)

    # ---- buckets. UNKNOWN can never reach Matches; that is the whole point.
    matches = list(combined.index[combined == T])
    possible = list(combined.index[combined == U])
    excluded = list(combined.index[combined == F])

    # ---- coverage, the most important line in the output
    resolved = len(matches) + len(excluded)
    per_term_cov = {}
    for lbl, st in per_term.items():
        per_term_cov[lbl] = {
            "resolved": int((st != U).sum()),
            "unknown": int((st == U).sum()),
            "resolved_frac": float((st != U).mean()),
            "top_unknown_reason": (reasons[lbl][st == U].value_counts().index[0]
                                   if (st == U).any() else None),
        }
    unknown_blame = {}
    for m in possible:
        culprits = [lbl for lbl in per_term if per_term[lbl][m] == U]
        unknown_blame[m] = [{"term": c, "why": reasons[c][m]} for c in culprits]

    # ---- abstention
    fe = gene_validity()
    warnings = []
    for t in terms:
        if t.kind == "expressed" and t.ensg in fe.index:
            f = float(fe[t.ensg])
            if f < SILENT_FRAC:
                warnings.append({
                    "term": t.label(), "condition": "SILENT_GENE",
                    "detail": (f"{t.symbol} is expressed in only {f:.1%} of lines "
                               f"(< {SILENT_FRAC:.0%}); this term resolves to "
                               f"UNKNOWN everywhere and constrains nothing")})
            elif t.negate and f > 0.95:
                warnings.append({
                    "term": t.label(), "condition": "NEAR_UNIVERSAL_GENE",
                    "detail": (f"{t.symbol} is expressed in {f:.1%} of lines, so "
                               f"NOT expressed({t.symbol}) excludes almost "
                               f"everything")})
        if t.kind == "altered" and t.negate:
            _, d, _mod = alterations()
            n_alt = len(d.get(t.ensg, set()))
            if n_alt < 10:
                warnings.append({
                    "term": t.label(), "condition": "TRIVIALLY_SATISFIED",
                    "detail": (f"only {n_alt} lines carry a driver alteration in "
                               f"{t.symbol}, so NOT altered({t.symbol}) is "
                               f"satisfied by almost every sequenced line and "
                               f"barely narrows the result")})

        # FIX 2. `NOT expressed(G)` on a tumour suppressor is a confidently wrong
        # question: truncating mutations leave abundant transcript because
        # nonsense-mediated decay is incomplete, so a protein-null line still
        # reads as expressed. Measured: the deletion_null category scores 25%
        # agreement (MULTIGENE_GROUND_TRUTH.md). Suggest the right term rather
        # than only flagging the problem.
        if t.kind == "expressed" and t.negate and t.ensg in tumour_suppressors():
            warnings.append({
                "term": t.label(), "condition": "LOSS_IS_PROTEIN_LEVEL",
                "detail": (f"{t.symbol} loss is typically protein-level. "
                           f"Truncating mutations leave transcript intact, so "
                           f"this term will not find {t.symbol}-null lines. "
                           f"Consider altered({t.symbol}) instead.")})

    # ---- per-line columns: each gene's value in its own column, never combined
    cols = pd.DataFrame({lbl: values[lbl] for lbl in per_term}, index=lines)
    cols["lineage"] = [str(lin.lineage.get(m, "unknown")) for m in lines]
    cols["state"] = combined

    # ---- clustering, descriptive
    clusters = []
    for lg, grp in cols.groupby("lineage"):
        clusters.append({
            "lineage": lg,
            "matches": int((grp.state == T).sum()),
            "possible": int((grp.state == U).sum()),
            "excluded": int((grp.state == F).sum()),
            "n_lines": int(len(grp)),
            "match_ids": list(grp.index[grp.state == T]),
        })
    if cluster_order == "count":
        clusters.sort(key=lambda c: (-c["matches"], -c["possible"]))
    else:
        clusters.sort(key=lambda c: c["lineage"])

    out = {
        "query": [t.label() for t in terms],
        "altered_scope": ALTERED_SCOPE,
        "coverage": {
            "lines_considered": len(lines),
            "lines_resolved": resolved,
            "resolved_frac": resolved / len(lines) if lines else 0.0,
            "per_term": per_term_cov,
        },
        "buckets": {"matches": matches, "possible": possible, "excluded": excluded},
        "counts": {"matches": len(matches), "possible": len(possible),
                   "excluded": len(excluded)},
        "unknown_blame": unknown_blame,
        "abstention_warnings": warnings,
        "clusters": clusters,
        "columns": cols,
        "elapsed_ms": (time.perf_counter() - t0) * 1000.0,
        "combined_score": None,
        "combined_score_note": (
            "no combined score is produced. Each term is its own column. Task B "
            "Stage 1 did not run: the pre-registered eligibility kill switch "
            "fired at 2 of 10 established pairs. See "
            "docs/COMBINED_SCORE_RESULT.md."),
    }
    if k_alternatives:
        out["alternatives"] = _alternatives(matches, k_alternatives)
    return out


def _alternatives(model_ids, k):
    """The validated similarity claim, per matched line. Strongest WITHIN
    lineage (+0.2519) against +0.1856 overall, so clustered output is the
    setting where this component performs best."""
    try:
        sys.path.insert(0, str(ROOT / "src" / "pipeline"))
        from cell_similarity_bridge import similar_lines_for, SIMILARITY_CLAIM
    except Exception:
        return {"available": False, "reason": "cell_similarity_bridge unavailable"}
    return {
        "available": True,
        "claim": SIMILARITY_CLAIM,
        "within_lineage_delta": 0.2519,
        "within_lineage_note": (
            "the similarity graph is strongest within lineage (+0.2519 vs "
            "+0.1856 overall), so the clustered view is where it performs best"),
        "per_line": {m: similar_lines_for(m, k=k) for m in model_ids[:50]},
    }


# =============================================================== CLI
def _print(res, top_clusters=8, show=6):
    q = " AND ".join(res["query"])
    print("=" * 78)
    print(f"QUERY  {q}")
    print("=" * 78)
    c = res["coverage"]
    print(f"  COVERAGE  {c['lines_resolved']:,} of {c['lines_considered']:,} lines "
          f"resolved ({c['resolved_frac']:.1%})")
    for lbl, d in c["per_term"].items():
        print(f"    {lbl:<34} {d['resolved']:>5,} resolved "
              f"({d['resolved_frac']:>5.1%})   {d['unknown']:>5,} unknown"
              + (f" -- {d['top_unknown_reason']}" if d["top_unknown_reason"] else ""))
    k = res["counts"]
    print(f"\n  MATCHES {k['matches']:,}    POSSIBLE {k['possible']:,}    "
          f"EXCLUDED {k['excluded']:,}")
    print(f"  (POSSIBLE = at least one clause UNKNOWN; these are NOT matches)")
    for w in res["abstention_warnings"]:
        print(f"\n  !! ABSTENTION [{w['condition']}] on {w['term']}")
        print(f"     {w['detail']}")
    print(f"\n  CLUSTERS (by match count; alternative is alphabetical)")
    print(f"    {'lineage':<26} {'match':>6} {'poss':>6} {'excl':>6} {'lines':>6}")
    for cl in res["clusters"][:top_clusters]:
        if cl["matches"] == 0 and cl["possible"] == 0:
            continue
        print(f"    {cl['lineage']:<26} {cl['matches']:>6} {cl['possible']:>6} "
              f"{cl['excluded']:>6} {cl['n_lines']:>6}")
    cols = res["columns"]
    m = cols[cols.state == "TRUE"]
    if len(m):
        valcols = [c for c in cols.columns if c not in ("lineage", "state")]
        # A7: "sorted by a number you can see". Sort by the FIRST term's value,
        # descending, and name it in the header. Previously this printed
        # `m.head(show)` in expression-matrix order under the label "TOP
        # MATCHES" -- an arbitrary slice implying a ranking that did not exist.
        sort_col = valcols[0] if valcols else None
        if sort_col is not None and m[sort_col].notna().any():
            m = m.sort_values(sort_col, ascending=False, na_position="last")
            print(f"\n  MATCHES, sorted by {sort_col} (descending)")
            print(f"  Each term is its own column; nothing is combined. Sort by "
                  f"whichever you care about --")
            print(f"  the order below is one choice, not a ranking the tool "
                  f"endorses.")
        else:
            print(f"\n  MATCHES (unordered -- no numeric column to sort on)")
        print("    " + "".join(f"{c[:22]:>24}" for c in valcols)
              + f"{'lineage':>22}")
        for mid, row in m.head(show).iterrows():
            print(f"    {mid:<12}"
                  + "".join(f"{row[c]:>24.3f}" if pd.notna(row[c]) else f"{'--':>24}"
                            for c in valcols)
                  + f"{row['lineage']:>22}")
    if res["counts"]["possible"]:
        print(f"\n  WHY LINES ARE UNRESOLVED (first 5)")
        for mid, why in list(res["unknown_blame"].items())[:5]:
            for w in why:
                print(f"    {mid:<14} {w['term']:<30} {w['why']}")
    if any(t.startswith("altered") or " altered(" in t for t in res["query"]):
        sc = res["altered_scope"]
        print(f"\n  WHAT altered() CAN SEE")
        print(f"    covered    : " + "; ".join(sc["covered"]))
        print(f"    NOT covered: " + "; ".join(sc["not_covered"][:2]))
        for x in sc["not_covered"][2:]:
            print(f"                 {x}")
        print(f"    altered()==FALSE means 'no SNV, indel or in-frame fusion "
              f"found', not 'not altered'.")
    print(f"\n  elapsed {res['elapsed_ms']:.1f} ms")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    # action="extend" so REPEATED flags accumulate. With the argparse default
    # each occurrence REPLACES the previous one, so
    #     --altered KRAS --altered BRAF
    # silently ran altered(BRAF) alone -- a different query from the one typed,
    # with no warning. Both forms now work and mean the same thing.
    ap.add_argument("--expressed", nargs="*", action="extend", default=[])
    ap.add_argument("--not-expressed", nargs="*", action="extend", default=[])
    ap.add_argument("--altered", nargs="*", action="extend", default=[])
    ap.add_argument("--not-altered", nargs="*", action="extend", default=[])
    ap.add_argument("--json", action="store_true")
    ap.add_argument("--no-alternatives", action="store_true")
    a = ap.parse_args()
    terms = parse_terms(a.expressed, a.not_expressed, a.altered, a.not_altered)
    if not terms:
        ap.error("give at least one term")
    res = query(terms, k_alternatives=0 if a.no_alternatives else 3)
    if a.json:
        res = {k: v for k, v in res.items() if k != "columns"}
        print(json.dumps(res, indent=2, default=str))
    else:
        _print(res)
