"""
export_web.py
-------------
Export real pipeline output as JSON for the React app.

    python src/pipeline/export_web.py --genes BRAF ERBB2 EGFR KRAS TP53 GAPDH
    python src/pipeline/export_web.py --genes BRAF --out ../Personal/genetraceai/public/data

WHY A STATIC EXPORT AND NOT AN API
    The React app (`Personal/genetraceai`) is a Vite + React SPA with no backend
    and no fetch calls -- its `GENES` object is a hand-written constant. The
    smallest honest bridge is therefore a static JSON file per gene, dropped into
    `public/data/`, which `vite build` serves as-is. No server to run, no CORS, no
    Python on the deploy target.

    `core_score.parquet` is 29.8M rows and cannot ship to a browser. One gene's
    slice is ~1,500 rows, and the top-N of it is a few KB. Exporting per gene keeps
    every payload small and lets the app fetch only what the user asked for.

WHAT IS AND IS NOT EXPORTED
    Exported: the gene's identity and dispersion, the abstention verdict, and the
    top-N cell lines with their full harmonised metadata and per-layer evidence.

    NOT exported: a 0-1 "confidence" number. The Confidence Engine design record
    (C1/C5) is explicit that the reported unit is an ordinal band with no
    calibration behind it. The React prototype currently renders a fabricated
    `confidence: 0.96` float; this export deliberately sends `confidence_tier`
    instead, plus `core_score` which IS a real quantity. Anything that wants a
    bar-length can use core_score and label it as such.

    The `verdict` block travels with the data, so the web app can refuse to draw a
    ranking for the same reasons the CLI refuses -- the abstention logic lives in
    the payload, not only in the terminal.
"""

from __future__ import annotations

import argparse
import json
import os
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from evidence_state import gene_verdict, load_dispersion            # noqa: E402
from metadata import (load_cell_lines, load_genes, cell_line_meta,   # noqa: E402
                      gene_meta, MODALITY_COLS)

ROOT = Path(__file__).resolve().parents[2]
PIPE = ROOT / "src" / "pipeline" / "outputs"
DEFAULT_OUT = ROOT.parent / "Personal" / "genetraceai" / "public" / "data"


def _pct(s):
    """Within-gene percentile across cell lines, ties to the floor.

    Same convention as 02_core_score Cell 7b (method='min'): a block of cell
    lines tied at zero sits at the bottom of the range it spans, not its midpoint.
    """
    s = pd.to_numeric(s, errors="coerce").dropna()
    if s.empty:
        return pd.Series(dtype=float)
    return s.rank(pct=True, method="min")


def layer_magnitudes(con, ensg, uniprots):
    """
    Real per-layer values for one gene, per cell line, as within-gene percentiles.

    This is what makes the omics map mean something. The alternative -- a
    present/absent flag per layer -- collapses every well-covered line to an
    identical row, which is true but carries no information. These are the actual
    measured quantities, ranked within the gene so the five layers are on one
    comparable 0-1 axis (the same transform the score itself uses).

    Returns {layer: {model_id: {"pct": float, "raw": float}}}.
    """
    out = {}

    def wide(table, col):
        cols = {c[0] for c in con.execute(f'DESCRIBE {table}').fetchall()}
        if col not in cols:
            return pd.Series(dtype=float)
        d = con.execute(f'SELECT model_id, "{col}" AS v FROM {table} '
                        f'WHERE model_id IS NOT NULL').df()
        return d.groupby("model_id")["v"].mean()

    out["depmap_expr"] = wide("depmap_expr", ensg)
    out["geo_expr"] = wide("geo_expr", ensg)

    hpa = con.execute("SELECT model_id, avg(ntpm) v FROM hpa_rna "
                      "WHERE lower(gene)=? AND model_id IS NOT NULL GROUP BY 1",
                      [ensg]).df()
    out["hpa_rna"] = hpa.set_index("model_id")["v"] if len(hpa) else pd.Series(dtype=float)

    for tbl, key in (("proteomics", "ccle_prot"), ("procan_proteomics", "procan_prot")):
        cols = {c[0] for c in con.execute(f"DESCRIBE {tbl}").fetchall()}
        have = [u for u in uniprots if u in cols]
        if not have:
            out[key] = pd.Series(dtype=float)
            continue
        expr = " + ".join(f'coalesce("{u}", 0)' for u in have)
        cnt = " + ".join(f'CASE WHEN "{u}" IS NULL THEN 0 ELSE 1 END' for u in have)
        d = con.execute(f'SELECT model_id, ({expr}) s, ({cnt}) n FROM {tbl} '
                        f'WHERE model_id IS NOT NULL').df()
        d = d[d["n"] > 0]
        out[key] = (d["s"] / d["n"]).groupby(d["model_id"]).mean() if len(d) else pd.Series(dtype=float)

    cna = con.execute("SELECT model_id, avg(total_cn) v FROM cosmic_cna "
                      "WHERE gene_id=? AND model_id IS NOT NULL GROUP BY 1", [ensg]).df()
    out["cna"] = cna.set_index("model_id")["v"] if len(cna) else pd.Series(dtype=float)

    packed = {}
    for k, s in out.items():
        if s.empty:
            packed[k] = {}
            continue
        p = _pct(s)
        packed[k] = {m: {"pct": round(float(p[m]), 4), "raw": round(float(s[m]), 4)}
                     for m in s.index if m in p.index}
    return packed


def _layer_evidence(row, cmeta, mags):
    """Per-layer evidence for one (gene, line): real magnitudes where they exist,
    honest booleans where the layer is genuinely categorical."""
    mods = set(cmeta.get("_modalities") or [])
    mid = row["model_id"]
    driver = bool(row.get("has_driver_alteration"))

    def m(layer):
        return (mags.get(layer) or {}).get(mid)

    dep, geo, hpa = m("depmap_expr"), m("geo_expr"), m("hpa_rna")
    expr_parts = [x for x in (dep, geo, hpa) if x]
    expr_pct = (sum(x["pct"] for x in expr_parts) / len(expr_parts)) if expr_parts else None

    ccle, procan = m("ccle_prot"), m("procan_prot")
    prot_parts = [x for x in (ccle, procan) if x]
    prot_pct = (sum(x["pct"] for x in prot_parts) / len(prot_parts)) if prot_parts else None

    cna = m("cna")

    return {
        # continuous, real: percentile of the measured value within this gene
        "expression": {
            "present": expr_pct is not None,
            "score": round(expr_pct, 4) if expr_pct is not None else 0,
            "detail": {"depmap_log2tpm": dep["raw"] if dep else None,
                       "geo": geo["raw"] if geo else None,
                       "hpa_ntpm": hpa["raw"] if hpa else None},
            "sources": [s for s, v in (("depmap_expr", dep), ("geo_expr", geo),
                                       ("hpa_rna", hpa)) if v],
        },
        "proteomics": {
            "present": prot_pct is not None,
            "score": round(prot_pct, 4) if prot_pct is not None else 0,
            "detail": {"ccle": ccle["raw"] if ccle else None,
                       "procan": procan["raw"] if procan else None},
            "sources": [s for s, v in (("proteomics", ccle),
                                       ("procan_proteomics", procan)) if v],
        },
        # CNA is NOT percentile-ranked. Copy number is heavily tied and capped --
        # for ERBB2, 18 of 34 lines sit at the ceiling of 14, and rank(method='min')
        # hands that tied top block 0.5, i.e. the strongest amplifications score
        # mid-range. Deviation from diploid is both artefact-free and directly
        # readable: CN 14 is a 7x amplification, CN 0 is a homozygous deletion.
        "cna": {
            "present": cna is not None,
            "score": round(min(1.0, abs(cna["raw"] - 2) / 6.0), 4) if cna else 0,
            "total_cn": cna["raw"] if cna else None,
            "direction": (None if not cna else
                          "amplification" if cna["raw"] > 2.5 else
                          "deletion" if cna["raw"] < 1.5 else "neutral"),
        },
        # genuinely categorical -- a driver call is not a magnitude, and
        # inventing one would be the exact failure this pipeline avoids
        "mutation": {"present": "mutations" in mods,
                     "score": 1.0 if driver else (0.35 if "mutations" in mods else 0),
                     "driver_alteration": driver},
        "fusion": {"present": "fusions" in mods,
                   "score": 1.0 if "fusions" in mods else 0},
        "signature": {"present": "signatures" in mods,
                      "score": 1.0 if "signatures" in mods else 0},
    }


def resolve_symbol(genes, symbol):
    """ENSG for a symbol, or None.

    hugo_symbol in gene_enriched is sourced from the mutations table, so a gene
    never seen mutated has it blank (MUC1). Fall back to the alias list, which is
    unioned across every source, then to a direct ENSG.
    """
    sym = str(symbol).lower()
    hit = genes[genes["hugo_symbol"] == sym]
    if len(hit):
        return hit.index[0]
    if sym.startswith("ensg") and sym in genes.index:
        return sym
    def has_alias(v):
        # gene_names is a list column; missing values arrive as pd.NA, which is
        # neither None nor iterable
        if v is None or not hasattr(v, "__iter__") or isinstance(v, str):
            return False
        try:
            return sym in [str(x).lower() for x in v]
        except TypeError:
            return False

    hit = genes[genes["gene_names"].map(has_alias).fillna(False).astype(bool)]
    return hit.index[0] if len(hit) else None


def export_gene(symbol, pred, cells, genes, dispersion, top_n, out_dir, con=None):
    ensg = resolve_symbol(genes, symbol)
    if ensg is None:
        return None, f"{symbol}: no gene in the warehouse carries that symbol"
    grow = genes.loc[[ensg]]

    g = pred[pred["ensg_id"] == ensg]
    if g.empty:
        g = pd.DataFrame(columns=pred.columns)

    disp = dispersion.loc[[ensg]] if dispersion is not None and ensg in dispersion.index else None
    gclass = g["class"].iloc[0] if len(g) and "class" in g.columns else None
    v = gene_verdict(ensg, symbol.lower(), len(g), disp, contradicted=False,
                     gene_class=gclass)

    ranked = (g.sort_values("core_score", ascending=False)
                .drop_duplicates("model_id").head(top_n))

    grow0 = grow.iloc[0]
    # `or []` is unsafe here: uniprot_ids is a numpy array, and a multi-element
    # array raises on truth-testing. CDKN2A is the first gene with two accessions.
    _u = grow0.get("uniprot_ids")
    uniprots = [] if _u is None else [str(x) for x in list(_u)]
    mags = layer_magnitudes(con, ensg, uniprots) if con is not None else {}

    lines = []
    for _, r in ranked.iterrows():
        cm = cell_line_meta(r["model_id"], cells)
        lines.append({
            "model_id": r["model_id"],
            "name": cm.get("name") or r["model_id"],
            "lineage": cm.get("lineage"),
            "disease": cm.get("disease") or cm.get("cancer type (Sanger)"),
            "subtype": cm.get("subtype") or cm.get("lineage subtype"),
            "tissue": cm.get("tissue (Sanger)"),
            # a real quantity, clearly named -- NOT a calibrated confidence
            "core_score": round(float(r["core_score"]), 4),
            "n_layers": int(r.get("n_layers") or 0),
            "confidence_tier": r.get("confidence"),
            "confidence_reason": r.get("confidence_reason"),
            "rank_basis": r.get("rank_basis"),
            "n_modalities": cm.get("_n_modalities"),
            "metadata": {k: v2 for k, v2 in cm.items() if not k.startswith("_")},
            "tracks": _layer_evidence(r, cm, mags),
        })

    payload = {
        "symbol": symbol.upper(),
        "ensg": ensg,
        "gene": gene_meta(ensg, genes),
        "class": gclass,
        "n_candidates": int(len(g)),
        "verdict": {
            "state": v.state,
            "headline": v.headline,
            "detail": v.detail,
            "show_ranking": v.show_ranking,
            "qualify_ranking": v.qualify_ranking,
        },
        "lines": lines,
        "provenance": {
            "core_score_rows": None,
            "note": ("core_score is a WITHIN-GENE percentile; it is comparable across "
                     "cell lines for one gene, not across genes. confidence_tier is an "
                     "ordinal band with no calibration behind it (design record C5)."),
        },
    }

    out_dir.mkdir(parents=True, exist_ok=True)
    path = out_dir / f"{symbol.upper()}.json"
    path.write_text(json.dumps(payload, indent=1, default=str), encoding="utf-8")
    return path, None


def main():
    ap = argparse.ArgumentParser(description="Export gene payloads as JSON for the web app")
    ap.add_argument("--genes", nargs="+", required=True)
    ap.add_argument("--top", type=int, default=25, help="cell lines per gene")
    ap.add_argument("--out", default=str(DEFAULT_OUT))
    args = ap.parse_args()

    out_dir = Path(args.out)
    print(f"reading predictions ...")
    cols = ["model_id", "ensg_id", "core_score", "n_layers", "class",
            "confidence", "confidence_reason", "rank_basis", "has_driver_alteration"]
    genes_l = [g.lower() for g in args.genes]
    gtab = load_genes()
    if gtab is None:
        raise SystemExit("warehouse not found - run build_warehouse_views.py first")
    ensgs = [e for e in (resolve_symbol(gtab, g) for g in args.genes) if e]
    if not ensgs:
        raise SystemExit("none of the requested symbols resolved to a gene")
    pred = pd.read_parquet(PIPE / "predictions_with_confidence.parquet",
                           columns=cols, filters=[("ensg_id", "in", ensgs)])
    cells = load_cell_lines()
    dispersion = load_dispersion()
    import duckdb
    con = duckdb.connect(str(PIPE / "celllineselector.db"), read_only=True)

    written, index = [], []
    for sym in args.genes:
        path, err = export_gene(sym, pred, cells, gtab, dispersion, args.top, out_dir, con)
        if err:
            print(f"  {err}")
            continue
        size = path.stat().st_size / 1024
        print(f"  {path.name:<16} {size:7.1f} KB")
        written.append(path)
        index.append({"symbol": sym.upper(), "file": path.name})

    con.close()

    # Rebuild the index from what is actually on disk, not just this run --
    # otherwise exporting one extra gene silently drops every gene exported
    # before it.
    on_disk = sorted(f.stem for f in out_dir.glob("*.json") if f.stem != "index")
    (out_dir / "index.json").write_text(
        json.dumps({"genes": [{"symbol": s, "file": f"{s}.json"} for s in on_disk]},
                   indent=1), encoding="utf-8")
    print(f"\nwrote {len(written)} gene file(s); index now lists {len(on_disk)} genes")
    print(f"out: {out_dir}")


if __name__ == "__main__":
    main()
