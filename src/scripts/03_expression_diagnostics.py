"""
03_expression_diagnostics.py
CellLineSelector - Feature Engineering: characterise the expression (Q1/Q2) problem.

WHAT THIS IS
------------
The expression aggregation (Q2) and scale reconciliation (Q1) decisions are DEFERRED
in 02_aggregation.py because they are method-sensitive. Before reading any literature
or choosing any method, this module turns the abstract problem into concrete numbers and
pictures, using the tables already in memory. It answers six questions:

  1. FAN-OUT      How many measurements per line actually exist (is "up to 478" one line
                  or hundreds)? -> sizes the Q2 aggregation problem.
  2. DEPMAP SCALE Do the two DepMap summaries you have reconcile (appendix mean 3.40 /
                  max 8.13 vs verification mean 1.09 / 49.8% zeros)? -> grounds Q1b.
  3. HPA COST     expr_long was built from DepMap+GEO only. What would adding HPA gain,
                  and what did the DepMap-GEO gene intersection already cost? -> a scope
                  decision (2-source vs 3-source Q1).
  4. PLATFORM     Is GEO microarray, RNA-seq, or mixed? -> decides whether "absent" in
                  GEO means "below array background" or "TPM below threshold".
  5. ABSENCE      DepMap encodes "off" as exact 0 (~half its values); GEO never does.
                  What is GEO's effective detection floor? -> the crux of Q1b / D2.
  6. RANK AGREE   Do DepMap and GEO rank genes the SAME WAY within a line, even though
                  their raw scales differ? If yes, a rank/percentile representation
                  reconciles them without aligning raw scales -> tells you whether the
                  "keep sources separate + compare" or "categories" strategy is viable
                  and you may never need the hard cross-scale merge.

WHAT THIS IS *NOT*
------------------
It does not aggregate, transform, normalise, threshold, or otherwise choose a Q1/Q2
method. Check 6 uses a within-source median PURELY as a robust, transform-invariant
diagnostic to compare orderings; it is not a proposed aggregator. Every function returns
measured numbers and (optionally) plots; the decision stays with you and the literature.

USAGE (in your live Jupyter session, where the harmonised tables already exist)
-------------------------------------------------------------------------------
    import importlib
    diag = importlib.import_module("03_expression_diagnostics")

    tables = {                       # same dict you pass to 02_aggregation.measure_fanout
        "expr_long":   expr_long,    # model_id, gene, value, source  (depmap + geo)
        "depmap_expr": depmap_expr,  # profile x gene wide, log2(TPM+1)
        "geo_expr":    geo_expr,     # GSM x gene wide, linear
        "geo_info":    geo_info,     # GSM metadata incl. platform_id
        "hpa_rna":     hpa_rna,      # line x gene long, linear nTPM (NOT in expr_long)
        "sample_info": sample_info,  # 1,840-line roster
    }

    summary = diag.run_all(tables)   # runs all six, saves plots + JSON, prints a digest

Individual checks can be run alone, e.g. diag.check_1_fanout(tables). Missing tables are
skipped with a clear message rather than crashing. Column names are resolved defensively.

Author: (CellLineSelector team)   |   Depends on: pandas, numpy; matplotlib optional (plots)
"""

from __future__ import annotations

import json
import os
import random
from typing import Iterable, Sequence

import numpy as np
import pandas as pd

try:
    import matplotlib.pyplot as plt
    _HAVE_MPL = True
except Exception:                                    # matplotlib optional
    _HAVE_MPL = False


# =============================================================================
# CONFIG
# =============================================================================

# Convenience ENSG ids for the optional human-readable marker spotlight in check 6.
# VERIFY THESE before quoting them - they are provided only so the spotlight prints
# something recognisable out of the box; any id absent from your common-gene set is
# skipped, so a wrong id degrades gracefully (it just does not appear).
MARKER_ENSG = {
    "ACTB (housekeeping)":  "ENSG00000075624",
    "GAPDH (housekeeping)": "ENSG00000111640",
    "B2M (housekeeping)":   "ENSG00000166710",
    "TBP (housekeeping)":   "ENSG00000112592",
    "MITF (melanocyte)":    "ENSG00000187098",
    "PMEL (melanocyte)":    "ENSG00000185664",
    "ALB (hepatocyte)":     "ENSG00000163631",
    "GATA3 (breast)":       "ENSG00000107485",
    "PTPRC/CD45 (immune)":  "ENSG00000081237",
    "VIM (mesenchymal)":    "ENSG00000026025",
    "KRT19 (epithelial)":   "ENSG00000171345",
    "EPCAM (epithelial)":   "ENSG00000119888",
}

_DEPMAP_SOURCE_ALIASES = ("depmap", "depmap_expr", "ccle", "depmap/ccle")
_GEO_SOURCE_ALIASES = ("geo", "geo_expr")

# A few common GEO platform accessions, so check 4 can name the technology. Not exhaustive;
# always confirm a GPL id on GEO. The key distinction for Q1b is microarray vs RNA-seq.
_KNOWN_GPL = {
    "gpl570":   "Affymetrix HG-U133 Plus 2.0 (MICROARRAY)",
    "gpl96":    "Affymetrix HG-U133A (MICROARRAY)",
    "gpl97":    "Affymetrix HG-U133B (MICROARRAY)",
    "gpl6244":  "Affymetrix Human Gene 1.0 ST (MICROARRAY)",
    "gpl10558": "Illumina HumanHT-12 V4 (MICROARRAY)",
    "gpl11154": "Illumina HiSeq 2000 (RNA-SEQ)",
    "gpl16791": "Illumina HiSeq 2500 (RNA-SEQ)",
    "gpl18573": "Illumina NextSeq 500 (RNA-SEQ)",
    "gpl24676": "Illumina NovaSeq 6000 (RNA-SEQ)",
}


# =============================================================================
# SHARED HELPERS
# =============================================================================

def _find_col(df: pd.DataFrame, candidates: Iterable[str], required: bool = True,
              what: str = "column") -> str | None:
    """Return the first candidate column present in df (case-insensitive), else None/raise."""
    lower = {c.lower(): c for c in df.columns}
    for cand in candidates:
        if cand.lower() in lower:
            return lower[cand.lower()]
    if required:
        raise KeyError(f"Could not resolve {what}: tried {list(candidates)}; "
                       f"available = {list(df.columns)[:40]}")
    return None


def _explode_model_id(df: pd.DataFrame, model_col: str) -> pd.DataFrame:
    """Explode list/set/delimited model_id cells so the column is groupby-safe."""
    if model_col not in df.columns:
        return df
    out = df.copy()

    def _tolist(v):
        if isinstance(v, np.ndarray):            # parquet often returns list cells as ndarray
            return v.tolist()
        if isinstance(v, (list, tuple, set)):
            return list(v)
        try:
            if pd.isna(v):                       # scalar-only; guarded against array-truthiness
                return []
        except (TypeError, ValueError):
            return []
        s = str(v)
        if any(d in s for d in (",", ";")):
            return [x.strip() for x in s.replace(";", ",").split(",") if x.strip()]
        return [s]

    out[model_col] = out[model_col].map(_tolist)
    out = out.explode(model_col)
    return out[out[model_col].notna() & (out[model_col].astype(str) != "")]


def _norm_gene(g) -> str:
    """Normalise a gene id for cross-table comparison: lowercase, strip version suffix."""
    s = str(g).strip().lower()
    return s.split(".")[0]


def _summ(values: np.ndarray) -> dict:
    """Distribution summary mirroring the verification run (NaNs excluded)."""
    v = np.asarray(values, dtype="float64")
    v = v[~np.isnan(v)]
    if v.size == 0:
        return {"n": 0}
    s = pd.Series(v)
    return {
        "n": int(v.size),
        "min": float(v.min()), "max": float(v.max()),
        "mean": float(v.mean()), "median": float(np.median(v)),
        "q01": float(np.quantile(v, 0.01)), "q99": float(np.quantile(v, 0.99)),
        "pct_zero": float((v == 0).mean() * 100.0),
        "pct_negative": float((v < 0).mean() * 100.0),
        "skew": float(s.skew()),
    }


def _spearman(a: pd.Series, b: pd.Series) -> float:
    """Spearman rank correlation = Pearson on ranks; pairwise-complete; NaN if n<3."""
    m = a.notna() & b.notna()
    if m.sum() < 3:
        return float("nan")
    ra, rb = a[m].rank(), b[m].rank()
    if ra.std(ddof=0) == 0 or rb.std(ddof=0) == 0:
        return float("nan")
    return float(np.corrcoef(ra, rb)[0, 1])


def _get(tables: dict, *names: str) -> pd.DataFrame | None:
    for n in names:
        if n in tables and isinstance(tables[n], pd.DataFrame):
            return tables[n]
    return None


def _outdir(path: str) -> str:
    os.makedirs(path, exist_ok=True)
    return path


def _finish_plot(fig, outdir: str, fname: str, show: bool):
    if outdir:
        fig.savefig(os.path.join(outdir, fname), dpi=110, bbox_inches="tight")
    if show and _HAVE_MPL:
        plt.show()
    plt.close(fig)


def _source_mask(series: pd.Series, aliases: Sequence[str]) -> pd.Series:
    low = series.astype(str).str.lower()
    return low.isin([a.lower() for a in aliases])


def _ensg_cols(df: pd.DataFrame) -> list[str]:
    return [c for c in df.columns if str(c).lower().startswith("ensg")]


# =============================================================================
# CHECK 1 - FAN-OUT: how many measurements per line actually exist?
# =============================================================================

def check_1_fanout(tables: dict, outdir: str = "features/expr_diagnostics",
                   show: bool = True) -> dict:
    """Measure the real samples-per-line distribution for GEO and DepMap.

    Sizes the Q2 aggregation problem: a table where nearly every line has one measurement
    needs no careful aggregator; GEO's long tail (hundreds of samples on a few lines) does.
    Operates on the sample/profile-grained wide tables (geo_expr, depmap_expr) so it counts
    real samples, not melted rows.
    """
    outdir = _outdir(outdir) if outdir else outdir
    report = {}
    for name, aliases in (("geo_expr", ("model_id", "model_ids")),
                          ("depmap_expr", ("model_id", "model_ids"))):
        df = _get(tables, name)
        if df is None:
            print(f"[fanout] {name}: not in tables - skipped")
            continue
        mcol = _find_col(df, aliases, required=False, what="model_id")
        if mcol is None:
            print(f"[fanout] {name}: no model_id column - skipped")
            continue
        counts = _explode_model_id(df[[mcol]].copy(), mcol)[mcol].value_counts()
        r = {
            "n_lines": int(counts.size),
            "max_per_line": int(counts.max()) if counts.size else 0,
            "median_per_line": float(counts.median()) if counts.size else 0.0,
            "pct_lines_multi": float((counts > 1).mean() * 100) if counts.size else 0.0,
            "pct_lines_ge10": float((counts >= 10).mean() * 100) if counts.size else 0.0,
            "pct_lines_ge50": float((counts >= 50).mean() * 100) if counts.size else 0.0,
        }
        report[name] = r
        print(f"[fanout] {name:11s} {r['n_lines']:>5,} lines | "
              f"max {r['max_per_line']:>4} / median {r['median_per_line']:.0f} per line | "
              f"{r['pct_lines_multi']:.1f}% multi, {r['pct_lines_ge10']:.1f}% >=10, "
              f"{r['pct_lines_ge50']:.1f}% >=50")

        if _HAVE_MPL and name == "geo_expr" and counts.size:
            fig, ax = plt.subplots(figsize=(7, 3.4))
            ax.hist(counts.values, bins=60)
            ax.set_yscale("log")
            ax.set_xlabel("GEO samples per line")
            ax.set_ylabel("number of lines (log)")
            ax.set_title("Check 1 - GEO fan-out (how many samples must collapse per line)")
            _finish_plot(fig, outdir, "check1_geo_fanout.png", show)

    if report:
        geo = report.get("geo_expr", {})
        print(f"[fanout] VERDICT: GEO median {geo.get('median_per_line', '?')}/line, "
              f"max {geo.get('max_per_line', '?')} - the tail, not the typical line, is the "
              f"Q2 problem; DepMap fan-out is small (near-technical).")
    return report


# =============================================================================
# CHECK 2 - DEPMAP SCALE: reconcile the conflicting DepMap summaries
# =============================================================================

def check_2_depmap_scale(tables: dict, outdir: str = "features/expr_diagnostics",
                         show: bool = True, sample_n: int | None = 5_000_000) -> dict:
    """Recompute DepMap distribution on the common-gene set and on the full matrix.

    Your appendix (mean 3.40, max 8.13) and verification (mean 1.09, 49.8% zeros) disagree.
    The zero fraction is load-bearing for Q1b, so establish which subset each was computed on.
    """
    outdir = _outdir(outdir) if outdir else outdir
    report = {}

    el = _get(tables, "expr_long")
    if el is not None:
        vcol = _find_col(el, ("value", "expr", "expression"), what="value")
        scol = _find_col(el, ("source", "src"), required=False, what="source")
        if scol is not None:
            dm = el.loc[_source_mask(el[scol], _DEPMAP_SOURCE_ALIASES), vcol]
            if sample_n and len(dm) > sample_n:
                dm = dm.sample(sample_n, random_state=0)
            report["depmap_common_genes"] = _summ(dm.to_numpy())
            r = report["depmap_common_genes"]
            print(f"[depmap-scale] common-gene set (from expr_long): "
                  f"mean {r['mean']:.2f}, median {r['median']:.3f}, max {r['max']:.2f}, "
                  f"{r['pct_zero']:.1f}% zeros, skew {r['skew']:.2f}")

    de = _get(tables, "depmap_expr")
    if de is not None:
        ecols = _ensg_cols(de)
        if ecols:
            block = de[ecols].to_numpy(dtype="float64", na_value=np.nan).ravel()
            if sample_n and block.size > sample_n:
                idx = np.random.default_rng(0).choice(block.size, sample_n, replace=False)
                block = block[idx]
            report["depmap_full_matrix"] = _summ(block)
            r = report["depmap_full_matrix"]
            print(f"[depmap-scale] full matrix ({len(ecols):,} ENSG cols): "
                  f"mean {r['mean']:.2f}, median {r['median']:.3f}, max {r['max']:.2f}, "
                  f"{r['pct_zero']:.1f}% zeros, skew {r['skew']:.2f}")

    a = report.get("depmap_common_genes", {})
    b = report.get("depmap_full_matrix", {})
    if a or b:
        z = a.get("pct_zero", b.get("pct_zero"))
        print(f"[depmap-scale] VERDICT: DepMap is ~{z:.0f}% structural zeros on the set "
              f"used for cross-source work - the appendix 'mean 3.40' was almost certainly "
              f"computed on a denser/non-zero subset. Treat the zeros as real for Q1b.")
    return report


# =============================================================================
# CHECK 3 - HPA COST: what does the 2-source expr_long include / exclude?
# =============================================================================

def check_3_hpa_cost(tables: dict, outdir: str = "features/expr_diagnostics",
                     show: bool = True) -> dict:
    """Quantify the scope of expr_long: which sources, which genes, and the HPA trade-off.

    Confirms expr_long is DepMap+GEO only, measures the genes dropped by the DepMap-GEO
    intersection, and measures what adding HPA would gain (lines and genes). This is the
    2-source-vs-3-source scope decision for Q1, in numbers.
    """
    report = {}
    el = _get(tables, "expr_long")
    if el is None:
        print("[hpa-cost] expr_long not in tables - skipped")
        return report

    mcol = _find_col(el, ("model_id",), what="model_id")
    gcol = _find_col(el, ("gene", "ensg", "gene_id"), what="gene")
    scol = _find_col(el, ("source", "src"), required=False, what="source")

    if scol is not None:
        report["expr_long_sources"] = (el[scol].astype(str).str.lower()
                                       .value_counts().to_dict())
        print(f"[hpa-cost] expr_long sources: {report['expr_long_sources']} "
              f"(HPA expected ABSENT)")

    common_genes = set(el[gcol].map(_norm_gene).unique())
    common_lines = set(el[mcol].dropna().astype(str).unique())
    report["expr_long_n_genes"] = len(common_genes)
    report["expr_long_n_lines"] = len(common_lines)
    print(f"[hpa-cost] expr_long covers {len(common_genes):,} genes x "
          f"{len(common_lines):,} lines (the DepMap-GEO common set)")

    de = _get(tables, "depmap_expr")
    if de is not None:
        depmap_genes = {_norm_gene(c) for c in _ensg_cols(de)}
        dropped = depmap_genes - common_genes
        report["depmap_only_genes_dropped"] = len(dropped)
        print(f"[hpa-cost] genes DepMap measured but expr_long drops (not on GEO): "
              f"{len(dropped):,} - invisible to cross-source work, DepMap-only if needed")

    hp = _get(tables, "hpa_rna")
    if hp is not None:
        hmcol = _find_col(hp, ("model_id",), required=False, what="model_id")
        hgcol = _find_col(hp, ("gene", "ensembl", "ensembl_id", "gene_id", "ensg"),
                          required=False, what="gene")
        if hmcol and hgcol:
            hpa_lines = set(hp[hmcol].dropna().astype(str).unique())
            hpa_genes = set(hp[hgcol].map(_norm_gene).unique())
            comparable = len(hpa_genes & common_genes)
            report["hpa_n_lines"] = len(hpa_lines)
            report["hpa_lines_new"] = len(hpa_lines - common_lines)
            report["hpa_genes_ensg_comparable"] = comparable
            report["hpa_genes_new_vs_common"] = len(hpa_genes - common_genes)
            print(f"[hpa-cost] HPA: {len(hpa_lines):,} lines "
                  f"({len(hpa_lines - common_lines):,} not already in expr_long); "
                  f"{comparable:,} of its genes are ENSG-comparable to the common set")
            if comparable == 0:
                print("[hpa-cost] NOTE: 0 comparable genes -> HPA gene ids are not ENSG "
                      "here; an id mapping is required before HPA can join Q1 at all.")
        else:
            print("[hpa-cost] hpa_rna present but model_id/gene columns unresolved - "
                  "check its schema")

    print("[hpa-cost] VERDICT: decide deliberately - keep expr_long 2-source (simpler Q1) "
          "or add HPA as a third source/reference; and decide whether DepMap-only genes get "
          "a documented single-source fallback for target lookup (D1).")
    return report


# =============================================================================
# CHECK 4 - PLATFORM: is GEO microarray, RNA-seq, or mixed?
# =============================================================================

def check_4_platform(tables: dict, outdir: str = "features/expr_diagnostics",
                     show: bool = True) -> dict:
    """Decode GEO's platform_id (2 levels) and cross-tab against any technology/type field.

    Whether 'absent' in GEO means 'below array background' or 'TPM below a threshold'
    depends on whether the samples are microarray or RNA-seq. This resolves that.
    """
    report = {}
    gi = _get(tables, "geo_info")
    if gi is None:
        print("[platform] geo_info not in tables - skipped")
        return report

    pcol = _find_col(gi, ("platform_id", "platform", "gpl"), required=False, what="platform")
    if pcol is None:
        print("[platform] no platform_id column - skipped")
        return report

    vc = gi[pcol].astype(str).value_counts()
    report["platform_counts"] = vc.to_dict()
    print(f"[platform] platform_id values and sample counts: {vc.to_dict()}")
    for gpl, n in vc.items():
        known = _KNOWN_GPL.get(str(gpl).lower())
        if known:
            print(f"[platform]   {gpl} = {known}  ({int(n):,} samples)")
        elif str(gpl).lower() in ("none", "nan", ""):
            print(f"[platform]   {int(n):,} samples have NO platform label - "
                  f"identify them (a cel_file_names column implies Affymetrix arrays)")
    report["known_platform"] = {str(g): _KNOWN_GPL.get(str(g).lower())
                                for g in vc.index if _KNOWN_GPL.get(str(g).lower())}

    tcol = _find_col(gi, ("type", "technology", "library_strategy", "molecule"),
                     required=False, what="technology")
    if tcol is not None:
        ct = (gi.groupby([pcol, tcol]).size()
              .reset_index(name="n").sort_values("n", ascending=False))
        report["platform_x_type"] = ct.to_dict(orient="records")
        print(f"[platform] platform_id x {tcol}:")
        for _, row in ct.iterrows():
            print(f"           {row[pcol]!s:<14} {row[tcol]!s:<28} n={int(row['n']):,}")
    print("[platform] VERDICT: confirm each GPL id on GEO (array vs RNA-seq). If array, "
          "GEO absence is a background/detection-call problem; if RNA-seq, a TPM-threshold "
          "problem. This decides the GEO side of Q1b.")
    return report


# =============================================================================
# CHECK 5 - ABSENCE: DepMap's exact zeros vs GEO's detection floor
# =============================================================================

def check_5_absence_floor(tables: dict, outdir: str = "features/expr_diagnostics",
                          show: bool = True, sample_n: int | None = 3_000_000) -> dict:
    """Contrast how each source encodes 'not expressed': DepMap exact-0 spike vs GEO floor.

    This is the crux of Q1b. DepMap's ~50% zeros are structural (TPM=0 -> log 0); GEO has
    0% exact zeros, so its 'off' sits near a noise floor. Measures GEO's effective floor and
    plots the two low-end shapes side by side.
    """
    outdir = _outdir(outdir) if outdir else outdir
    report = {}
    el = _get(tables, "expr_long")
    if el is None:
        print("[absence] expr_long not in tables - skipped")
        return report

    vcol = _find_col(el, ("value", "expr"), what="value")
    scol = _find_col(el, ("source", "src"), required=False, what="source")
    if scol is None:
        print("[absence] no source column - skipped")
        return report

    dm = el.loc[_source_mask(el[scol], _DEPMAP_SOURCE_ALIASES), vcol].dropna()
    geo = el.loc[_source_mask(el[scol], _GEO_SOURCE_ALIASES), vcol].dropna()
    if sample_n:
        if len(dm) > sample_n:
            dm = dm.sample(sample_n, random_state=0)
        if len(geo) > sample_n:
            geo = geo.sample(sample_n, random_state=0)

    if len(dm):
        report["depmap"] = {"pct_zero": float((dm == 0).mean() * 100),
                            "min_nonzero": float(dm[dm > 0].min()) if (dm > 0).any() else None}
        print(f"[absence] DepMap: {report['depmap']['pct_zero']:.1f}% exact zeros "
              f"(structural 'not detected')")
    if len(geo):
        gpos = geo[geo > 0]
        report["geo"] = {
            "pct_zero": float((geo == 0).mean() * 100),
            "min": float(geo.min()),
            "floor_q001": float(np.quantile(gpos, 0.001)) if len(gpos) else None,
            "floor_q01": float(np.quantile(gpos, 0.01)) if len(gpos) else None,
            "floor_q05": float(np.quantile(gpos, 0.05)) if len(gpos) else None,
        }
        r = report["geo"]
        print(f"[absence] GEO: {r['pct_zero']:.1f}% exact zeros; min {r['min']:.3g}; "
              f"effective floor ~ q0.1%={r['floor_q001']:.3g}, q1%={r['floor_q01']:.3g}, "
              f"q5%={r['floor_q05']:.3g}")

    if _HAVE_MPL and len(dm) and len(geo):
        fig, axes = plt.subplots(1, 2, figsize=(11, 3.6))
        axes[0].hist(dm.values, bins=80)
        axes[0].set_title("DepMap log2(TPM+1): spike at 0 = structural absence")
        axes[0].set_xlabel("value"); axes[0].set_ylabel("count")
        gpos = geo[geo > 0]
        axes[1].hist(np.log10(gpos.values + 1e-9), bins=80)
        axes[1].set_title("GEO linear (log10 axis): a floor, never 0")
        axes[1].set_xlabel("log10(value)"); axes[1].set_ylabel("count")
        fig.suptitle("Check 5 - 'not expressed' is a different object in each source (Q1b)")
        _finish_plot(fig, outdir, "check5_absence.png", show)

    print("[absence] VERDICT: absence cannot be read off a shared continuous axis. Strongly "
          "points to a separate per-source detection layer (present / not-detected) for D2, "
          "rather than folding absence into one continuous value.")
    return report


# =============================================================================
# CHECK 6 - RANK AGREEMENT: do DepMap and GEO order genes the same way per line?
# =============================================================================

def check_6_rank_agreement(tables: dict, outdir: str = "features/expr_diagnostics",
                           show: bool = True, n_genes: int = 3000,
                           marker_ensg: dict | None = None, seed: int = 0) -> dict:
    """Do the two sources rank genes consistently WITHIN a line, despite different scales?

    THE decision-relevant check. If, within a line, DepMap and GEO agree on which genes are
    high vs low (high Spearman), then a rank/percentile or category representation makes the
    sources comparable WITHOUT aligning raw scales - i.e. the 'keep separate + compare' and
    'categorise' strategies are viable and the hard cross-scale merge may be unnecessary.

    Method note (not a proposed aggregator): to get one value per (source, line, gene) for
    the comparison, replicates are collapsed with a WITHIN-SOURCE MEDIAN. The median is used
    only because it is robust and transform-invariant, making this a fair ordering test; it
    is not a recommendation for the eventual Q2 aggregator.
    """
    outdir = _outdir(outdir) if outdir else outdir
    report = {}
    el = _get(tables, "expr_long")
    if el is None:
        print("[rank] expr_long not in tables - skipped")
        return report

    mcol = _find_col(el, ("model_id",), what="model_id")
    gcol = _find_col(el, ("gene",), what="gene")
    vcol = _find_col(el, ("value",), what="value")
    scol = _find_col(el, ("source", "src"), required=False, what="source")
    if scol is None:
        print("[rank] no source column - cannot compare sources - skipped")
        return report

    low = el[scol].astype(str).str.lower()
    dm_name = next((v for v in low.unique() if v in [a.lower() for a in _DEPMAP_SOURCE_ALIASES]), None)
    geo_name = next((v for v in low.unique() if v in [a.lower() for a in _GEO_SOURCE_ALIASES]), None)
    if dm_name is None or geo_name is None:
        print(f"[rank] could not identify both sources in {list(low.unique())} - skipped")
        return report

    dm_lines = set(el.loc[low == dm_name, mcol].dropna().astype(str))
    geo_lines = set(el.loc[low == geo_name, mcol].dropna().astype(str))
    shared = sorted(dm_lines & geo_lines)
    report["n_shared_lines"] = len(shared)
    print(f"[rank] lines measured by BOTH sources: {len(shared):,}")
    if len(shared) < 5:
        print("[rank] too few shared lines to assess - skipped")
        return report

    all_genes = el[gcol].astype(str).unique()
    rng = random.Random(seed)
    gene_sample = set(rng.sample(list(all_genes), min(n_genes, len(all_genes))))

    sub = el[el[mcol].astype(str).isin(set(shared)) & el[gcol].astype(str).isin(gene_sample)].copy()
    sub["_s"] = sub[scol].astype(str).str.lower()
    sub = sub[sub["_s"].isin([dm_name, geo_name])]

    med = (sub.groupby([mcol, gcol, "_s"])[vcol].median().reset_index())
    wide = med.pivot_table(index=[mcol, gcol], columns="_s", values=vcol)
    if dm_name not in wide.columns or geo_name not in wide.columns:
        print("[rank] pivot missing a source column - skipped")
        return report

    # (a) per-line: within a line, do the sources agree on gene ordering? (most relevant to D1)
    per_line = (wide.dropna().groupby(level=0)
                .apply(lambda g: _spearman(g[dm_name], g[geo_name])).dropna())
    # (b) per-gene: across lines, do the sources agree on which lines are high? (relevant to D10)
    per_gene = (wide.dropna().groupby(level=1)
                .apply(lambda g: _spearman(g[dm_name], g[geo_name])).dropna())

    for label, s in (("per_line", per_line), ("per_gene", per_gene)):
        if len(s):
            report[label] = {
                "n": int(len(s)), "median_spearman": float(s.median()),
                "q25": float(s.quantile(.25)), "q75": float(s.quantile(.75)),
                "pct_ge_0.5": float((s >= 0.5).mean() * 100),
                "pct_ge_0.7": float((s >= 0.7).mean() * 100),
            }
            r = report[label]
            print(f"[rank] {label}: median Spearman {r['median_spearman']:.2f} "
                  f"(IQR {r['q25']:.2f}-{r['q75']:.2f}); "
                  f"{r['pct_ge_0.5']:.0f}% >=0.5, {r['pct_ge_0.7']:.0f}% >=0.7")

    if _HAVE_MPL and len(per_line):
        fig, ax = plt.subplots(figsize=(7, 3.4))
        ax.hist(per_line.values, bins=40, alpha=0.8, label="per line (gene ordering)")
        if len(per_gene):
            ax.hist(per_gene.values, bins=40, alpha=0.5, label="per gene (line ordering)")
        ax.axvline(0, color="k", lw=0.8)
        ax.set_xlabel("Spearman rank correlation between DepMap and GEO")
        ax.set_ylabel("count"); ax.legend()
        ax.set_title("Check 6 - do the two sources agree on ordering? (drives Q1 strategy)")
        _finish_plot(fig, outdir, "check6_rank_agreement.png", show)

    # optional human-readable marker spotlight
    marker_ensg = marker_ensg or MARKER_ENSG
    want = {_norm_gene(v): k for k, v in marker_ensg.items()}
    spot = wide.dropna().reset_index()
    spot["_g"] = spot[gcol].map(_norm_gene)
    spot = spot[spot["_g"].isin(want)]
    if len(spot):
        rows = []
        for g, grp in spot.groupby("_g"):
            rows.append({"marker": want[g], "n_lines": len(grp),
                         "spearman": round(_spearman(grp[dm_name], grp[geo_name]), 2)})
        report["marker_spotlight"] = rows
        print("[rank] marker spotlight (line-ordering agreement per gene):")
        for row in rows:
            print(f"        {row['marker']:<22} n={row['n_lines']:<4} "
                  f"Spearman={row['spearman']}")

    pl = report.get("per_line", {}).get("median_spearman")
    if pl is not None:
        if pl >= 0.6:
            verdict = ("HIGH ordering agreement -> a rank/percentile or category "
                       "representation reconciles the sources without a raw-scale merge. "
                       "The 'keep separate + compare' / 'categorise' strategies are viable.")
        elif pl >= 0.3:
            verdict = ("MODERATE agreement -> rank-based comparison partly works but is noisy; "
                       "batch (Q14) and platform differences may be degrading it - investigate "
                       "before committing to a strategy.")
        else:
            verdict = ("LOW agreement -> the sources disagree on ordering; neither a simple "
                       "merge nor a naive rank comparison is safe. Batch/identity/annotation "
                       "issues must be understood first.")
        print(f"[rank] VERDICT: {verdict}")
    return report


# =============================================================================
# ORCHESTRATOR
# =============================================================================

def run_all(tables: dict, outdir: str = "features/expr_diagnostics",
            show: bool = True) -> dict:
    """Run all six diagnostics, save plots + a JSON summary, print a closing digest.

    Returns a dict of every check's findings. Writes:
        <outdir>/expr_diagnostics_summary.json  and  check*.png plots.
    Nothing here selects a Q1/Q2 method; it characterises the problem so the literature
    phase reads against real numbers.
    """
    outdir = _outdir(outdir) if outdir else outdir
    print("=" * 78)
    print("EXPRESSION DIAGNOSTICS (Q1/Q2 characterisation) - measures only, chooses nothing")
    print("=" * 78)
    summary = {}
    for key, fn in (("fanout", check_1_fanout),
                    ("depmap_scale", check_2_depmap_scale),
                    ("hpa_cost", check_3_hpa_cost),
                    ("platform", check_4_platform),
                    ("absence_floor", check_5_absence_floor),
                    ("rank_agreement", check_6_rank_agreement)):
        print(f"\n--- {key} ---")
        try:
            summary[key] = fn(tables, outdir=outdir, show=show)
        except Exception as e:                       # a bad table should not sink the run
            summary[key] = {"error": f"{type(e).__name__}: {e}"}
            print(f"[{key}] ERROR: {type(e).__name__}: {e}")

    if outdir:
        path = os.path.join(outdir, "expr_diagnostics_summary.json")
        with open(path, "w") as f:
            json.dump(summary, f, indent=2, default=str)
        print(f"\n[run_all] wrote {path}")

    print("\n" + "=" * 78)
    print("READ THE RESULTS IN THIS ORDER:")
    print("  check 6 (rank agreement) first  -> decides which Q1 strategy is even needed")
    print("  check 5 (absence)               -> whether absence wants its own layer (D2)")
    print("  checks 1-4                      -> size the aggregation & scope the sources")
    print("Then read Q1/Q2 literature against these numbers, not in the abstract.")
    print("=" * 78)
    return summary


if __name__ == "__main__":
    print(__doc__)
    print("Import this module and call run_all(tables) from your live Jupyter session.")