"""
final_pipeline/api/ranking.py
------------------------------
Data loading and co-selection logic for the /v1/rank endpoint.

Loaded once at startup via load_ranking_data(); everything after that is
pure in-memory pandas — no network, no LLM.
"""

from __future__ import annotations

import logging
import math
import os
from collections import Counter
from pathlib import Path

import pandas as pd

logger = logging.getLogger(__name__)

_PIPELINE = Path(__file__).resolve().parent.parent


def _env_path(env_var: str, default: Path) -> Path:
    val = os.getenv(env_var)
    return Path(val) if val else default


def ensure_loaded() -> None:
    """Idempotent data load. Safe to call from either the HTTP lifespan or the
    Bedrock handler — loads once, no-ops thereafter. Paths come from env vars
    (Lambda/EFS overrides) or the bundled defaults."""
    if is_ready():
        return
    load_ranking_data(
        predictions_path=_env_path("PREDICTIONS_PATH",
                                   _PIPELINE / "outputs" / "predictions_with_confidence.parquet"),
        gene_lookup_path=_env_path("GENE_LOOKUP_PATH",
                                   _PIPELINE / "reference" / "gene_lookup.parquet"),
        cell_lookup_path=_env_path("CELL_LOOKUP_PATH",
                                   _PIPELINE / "reference" / "cell_line_lookup.parquet"),
        db_path=_env_path("RANKING_DB_PATH",
                          _PIPELINE / "outputs" / "celllineselector.db"),
    )

_pred: pd.DataFrame | None = None  # index=ensg_id (categorical, sorted), cols: model_id, core_score
_sym_index: dict[str, tuple[str, str]] = {}  # upper(symbol) -> (ensg_id, hgnc_symbol)
_ensg_index: dict[str, tuple[str, str]] = {}  # upper(ensg_id) -> (ensg_id, hgnc_symbol)
_gene_lkp: pd.DataFrame | None = None
_cell_lkp: pd.DataFrame | None = None
_cell_name: dict[str, str] = {}  # lower(model_id) -> cell_line_name
_lineage_map: dict[str, str] = {}  # lower(model_id) -> lineage
_meta_map: dict[str, dict] = {}  # lower(model_id) -> {lineage, subtype, disease, sex}
_pred_path: Path | None = None  # kept so /gene/detail can pull the full row on demand
_db_path: Path | None = None


def load_ranking_data(predictions_path: Path, gene_lookup_path: Path,
                      cell_lookup_path: Path, db_path: Path) -> None:
    global _pred, _sym_index, _ensg_index, _gene_lkp, _cell_lkp, _cell_name
    global _lineage_map, _meta_map, _pred_path, _db_path

    if not predictions_path.exists():
        logger.warning("predictions not found at %s — /v1/rank will 503", predictions_path)
        return

    pred = pd.read_parquet(predictions_path, columns=["model_id", "ensg_id", "core_score"])
    logger.info("predictions loaded: %d rows", len(pred))
    # Single DataFrame indexed by gene, sorted for fast per-gene .loc slicing.
    # float32 halves the score column; categorical model_id + ensg_id replace
    # 40M repeated strings with small int codes — keeps the frame under ~1 GB so
    # it fits in a modest Lambda memory budget.
    pred["core_score"] = pred["core_score"].astype("float32")
    pred["model_id"] = pred["model_id"].astype("category")
    pred["ensg_id"] = pred["ensg_id"].astype("category")
    _pred = pred.set_index("ensg_id").sort_index()
    _pred_path = predictions_path
    logger.info("prediction index built: %d genes", _pred.index.nunique())

    _gene_lkp = pd.read_parquet(gene_lookup_path, columns=["ensg_id", "hgnc_symbol"])
    # Vectorized zip — far faster than iterrows over the gene table.
    for ensg, sym in zip(_gene_lkp["ensg_id"], _gene_lkp["hgnc_symbol"]):
        pair = (ensg, sym)
        _sym_index[sym.upper()] = pair
        _ensg_index[ensg.upper()] = pair

    if cell_lookup_path.exists():
        _cell_lkp = pd.read_parquet(cell_lookup_path, columns=["model_id", "cell_line_name"])
        _cell_lkp["model_id"] = _cell_lkp["model_id"].str.lower()
        _cell_name = {
            mid: name
            for mid, name in zip(_cell_lkp["model_id"], _cell_lkp["cell_line_name"])
            if isinstance(name, str)
        }

    _db_path = db_path if db_path.exists() else None
    if _db_path:
        try:
            import duckdb
            con = duckdb.connect(str(_db_path), read_only=True)
            rows = con.execute(
                "SELECT lower(model_id), lineage, lineage_subtype, "
                "primary_disease, sex FROM sample_info"
            ).fetchall()
            con.close()
            for mid, lineage, subtype, disease, sex in rows:
                meta = {}
                if lineage:
                    meta["lineage"] = lineage
                    _lineage_map[mid] = lineage
                if subtype:
                    meta["lineage_subtype"] = subtype
                if disease:
                    meta["primary_disease"] = disease
                if sex:
                    meta["sex"] = sex
                if meta:
                    _meta_map[mid] = meta
            logger.info("metadata map loaded: %d lines", len(_meta_map))
        except Exception as exc:
            logger.warning("metadata map load failed: %s", exc)


def is_ready() -> bool:
    return _pred is not None and _gene_lkp is not None


def _resolve(query: str) -> tuple[str, str]:
    q = query.strip().upper()
    if q in _sym_index:
        return _sym_index[q]
    if q in _ensg_index:
        return _ensg_index[q]
    raise ValueError(f"Gene not found: {query!r}")


def _lineage_model_ids(terms: list[str]) -> set[str] | None:
    if not terms or _db_path is None:
        return None
    try:
        import duckdb
        con = duckdb.connect(str(_db_path), read_only=True)
        clauses = " OR ".join(
            f"(lower(lineage) LIKE '%{t.lower()}%' OR lower(lineage_subtype) LIKE '%{t.lower()}%')"
            for t in terms
        )
        rows = con.execute(
            f"SELECT lower(model_id) FROM sample_info WHERE {clauses}"
        ).fetchall()
        con.close()
        return {r[0] for r in rows} or None
    except Exception as exc:
        logger.warning("lineage filter failed: %s", exc)
        return None


def _meta(model_id: str) -> dict:
    """Per-line metadata (lineage, subtype, disease, sex) or {} if unknown."""
    return _meta_map.get(model_id.lower(), {})


def _lineage_dist(model_ids) -> dict[str, int]:
    counts: Counter = Counter()
    for mid in model_ids:
        lin = _lineage_map.get(mid.lower())
        if lin:
            counts[lin] += 1
    return dict(counts.most_common())


def _gene_series(ensg: str, lineage_ids: set[str] | None) -> pd.Series:
    """Return core_score Series (index=model_id) for one gene, optionally filtered.

    Uses .loc on the sorted, gene-indexed frame — O(log n) slice, no full scan.
    """
    try:
        sub = _pred.loc[[ensg]]
    except KeyError:
        return pd.Series(dtype="float32", name="core_score",
                         index=pd.Index([], name="model_id"))
    # model_id is categorical in the big frame; cast the small per-gene slice back
    # to plain strings so the index carries only present ids (not all categories)
    # and downstream merges stay object-typed.
    s = sub.set_index(sub["model_id"].astype(str))["core_score"]
    s.index.name = "model_id"
    if lineage_ids is not None:
        s = s[s.index.isin(lineage_ids)]
    return s


def rank(gene_syms: list[str], lineages: list[str],
         floor: float, top_n: int) -> dict:
    resolved = [_resolve(s) for s in gene_syms]
    syms = [s for _, s in resolved]

    lineage_ids = _lineage_model_ids(lineages)

    # ── 1-gene: straight core_score ranking ────────────────────────────────
    if len(resolved) == 1:
        ensg, sym = resolved[0]
        s = _gene_series(ensg, lineage_ids).dropna().sort_values(ascending=False)
        sub = s.head(top_n).reset_index()
        sub.columns = ["model_id", "core_score"]
        if _cell_lkp is not None:
            sub = sub.merge(_cell_lkp, on="model_id", how="left")
        lines = [
            {
                "model_id": r["model_id"].upper(),
                "name": r.get("cell_line_name") if isinstance(r.get("cell_line_name"), str) else None,
                "joint_score": round(float(r["core_score"]), 6),
                "limiting_gene": sym,
                "scores": {sym: round(float(r["core_score"]), 6)},
                "metadata": _meta(r["model_id"]),
            }
            for _, r in sub.iterrows()
        ]
        return {"genes": syms, "lineage": lineages, "floor": floor,
                "total_passing": len(s), "lines": lines,
                "lineage_distribution": _lineage_dist(s.index.tolist())}

    # ── N-gene: inner join + weakest-link min ──────────────────────────────
    wide = None
    for ensg, sym in resolved:
        s = _gene_series(ensg, lineage_ids).rename(sym).reset_index()
        wide = s if wide is None else wide.merge(s, on="model_id", how="inner")

    if wide is None or wide.empty:
        return {"genes": syms, "lineage": lineages, "floor": floor,
                "total_passing": 0, "lines": []}

    wide["joint_score"] = wide[syms].min(axis=1, skipna=False)
    scoreable = wide["joint_score"].notna()
    wide.loc[scoreable, "limiting_gene"] = wide.loc[scoreable, syms].idxmin(axis=1)

    passing = wide[scoreable & (wide["joint_score"] >= floor)]
    shortlist = passing.sort_values("joint_score", ascending=False).head(top_n)

    if _cell_lkp is not None:
        shortlist = shortlist.merge(_cell_lkp[["model_id", "cell_line_name"]],
                                    on="model_id", how="left")

    lines = []
    for _, row in shortlist.iterrows():
        raw_name = row.get("cell_line_name")
        lines.append({
            "model_id": row["model_id"].upper(),
            "name": raw_name if isinstance(raw_name, str) else None,
            "joint_score": round(float(row["joint_score"]), 6),
            "limiting_gene": row.get("limiting_gene"),
            "scores": {s: round(float(row[s]), 6) for s in syms if pd.notna(row[s])},
            "metadata": _meta(row["model_id"]),
        })

    return {"genes": syms, "lineage": lineages, "floor": floor,
            "total_passing": len(passing), "lines": lines,
            "lineage_distribution": _lineage_dist(passing["model_id"].tolist())}


def exclude_many(gene_a: str, gene_bs: list[str], lineages: list[str], top_n: int) -> dict:
    """Like exclude(), generalized to N excluded genes: selectivity is
    score_a * Π(1 - score_bi) — a joint gene must pass, and every excluded
    gene independently pulls the ranking down."""
    ensg_a, sym_a = _resolve(gene_a)
    resolved_b = [_resolve(g) for g in gene_bs]
    syms_b = [s for _, s in resolved_b]

    lineage_ids = _lineage_model_ids(lineages)

    merged = _gene_series(ensg_a, lineage_ids).rename("score_a").reset_index()
    for ensg_b, sym_b in resolved_b:
        b = _gene_series(ensg_b, lineage_ids).rename(sym_b).reset_index()
        merged = merged.merge(b, on="model_id", how="inner")

    selectivity = merged["score_a"].astype("float64").copy()
    for sym_b in syms_b:
        selectivity = selectivity * (1.0 - merged[sym_b])
    merged["selectivity"] = selectivity

    merged = merged.sort_values("selectivity", ascending=False, na_position="last")
    scored = merged[merged["selectivity"].notna()]
    shortlist = scored.head(top_n)

    if _cell_lkp is not None:
        shortlist = shortlist.merge(_cell_lkp, on="model_id", how="left")

    lines = []
    for _, row in shortlist.iterrows():
        raw_name = row.get("cell_line_name")
        lines.append({
            "model_id": row["model_id"].upper(),
            "name": raw_name if isinstance(raw_name, str) else None,
            "score_a": round(float(row["score_a"]), 6),
            "exclusion_scores": {s: round(float(row[s]), 6) for s in syms_b},
            "selectivity": round(float(row["selectivity"]), 6),
            "metadata": _meta(row["model_id"]),
        })

    return {"gene_a": sym_a, "excluded_genes": syms_b, "lineage": lineages,
            "total_ranked": len(scored), "lines": lines,
            "lineage_distribution": _lineage_dist(scored["model_id"].tolist())}


def exclude(gene_a: str, gene_b: str, lineages: list[str], top_n: int) -> dict:
    ensg_a, sym_a = _resolve(gene_a)
    ensg_b, sym_b = _resolve(gene_b)

    lineage_ids = _lineage_model_ids(lineages)

    a = _gene_series(ensg_a, lineage_ids).rename("score_a").reset_index()
    b = _gene_series(ensg_b, lineage_ids).rename("score_b").reset_index()

    merged = a.merge(b, on="model_id", how="inner")
    merged["selectivity"] = merged["score_a"] * (1.0 - merged["score_b"])
    merged = merged.sort_values("selectivity", ascending=False, na_position="last")
    scored = merged[merged["selectivity"].notna()]
    shortlist = scored.head(top_n)

    if _cell_lkp is not None:
        shortlist = shortlist.merge(_cell_lkp, on="model_id", how="left")

    lines = []
    for _, row in shortlist.iterrows():
        raw_name = row.get("cell_line_name")
        lines.append({
            "model_id": row["model_id"].upper(),
            "name": raw_name if isinstance(raw_name, str) else None,
            "score_a": round(float(row["score_a"]), 6),
            "score_b": round(float(row["score_b"]), 6),
            "selectivity": round(float(row["selectivity"]), 6),
            "metadata": _meta(row["model_id"]),
        })

    return {"gene_a": sym_a, "gene_b": sym_b, "lineage": lineages,
            "total_ranked": len(scored), "lines": lines,
            "lineage_distribution": _lineage_dist(scored["model_id"].tolist())}


def _lineage_meta(model_id: str) -> dict:
    if _db_path is None:
        return {}
    try:
        import duckdb
        con = duckdb.connect(str(_db_path), read_only=True)
        row = con.execute(
            "SELECT lineage, lineage_subtype, primary_disease, sex, age, "
            "default_growth_pattern, primary_or_metastasis "
            "FROM main.sample_info WHERE model_id = ?",
            [model_id.lower()]
        ).fetchone()
        con.close()
        if row:
            keys = ["lineage", "lineage_subtype", "primary_disease", "sex",
                    "age", "growth_pattern", "primary_or_metastasis"]
            return {k: v for k, v in zip(keys, row) if v is not None}
    except Exception as exc:
        logger.warning("lineage meta failed: %s", exc)
    return {}


# sample_info columns surfaced in the detail metadata block, in the order the
# UI labels them. Kept flat (single dict) because that's what the UI expects.
_DETAIL_META_COLS = [
    "lineage", "lineage_subtype", "primary_disease", "subtype",
    "cellosaurus_ncit_disease", "primary_or_metastasis",
    "sample_collection_site", "default_growth_pattern", "sex", "age",
]


def _full_metadata(mid: str) -> dict:
    """Flat metadata dict for one line, keyed as the UI's detail view expects."""
    if _db_path is None:
        return {}
    try:
        import duckdb
        con = duckdb.connect(str(_db_path), read_only=True)
        cols = ", ".join(_DETAIL_META_COLS)
        row = con.execute(
            f"SELECT {cols} FROM main.sample_info WHERE lower(model_id) = ?",
            [mid.lower()],
        ).fetchone()
        con.close()
        if row:
            return {k: v for k, v in zip(_DETAIL_META_COLS, row) if v is not None}
    except Exception as exc:
        logger.warning("full metadata failed: %s", exc)
    return {}


def _prediction_row(ensg: str, mid: str) -> dict:
    """Pull the full prediction row for one (gene, line) straight from the
    parquet via DuckDB predicate pushdown — reads only the matching row, so it
    adds nothing to steady-state memory."""
    if _pred_path is None:
        return {}
    try:
        import duckdb
        con = duckdb.connect()
        row = con.execute(
            "SELECT n_layers, p_mutation, p_fusion, has_cna_alteration, "
            "has_driver_alteration, confidence_tier "
            "FROM read_parquet(?) WHERE ensg_id = ? AND lower(model_id) = ?",
            [str(_pred_path), ensg, mid.lower()],
        ).fetchone()
        con.close()
        if row:
            keys = ["n_layers", "p_mutation", "p_fusion", "has_cna_alteration",
                    "has_driver_alteration", "confidence_tier"]
            return dict(zip(keys, row))
    except Exception as exc:
        logger.warning("prediction row lookup failed: %s", exc)
    return {}


def _z_to_level(z) -> float | None:
    """Map an expression/proteomics z-score to a 0-1 level (normal-CDF percentile).
    z=0 -> 0.50, z=+2 -> ~0.98, z=-2 -> ~0.02. Gives the omics map a real gradient."""
    if z is None or pd.isna(z):
        return None
    return round(0.5 * (1.0 + math.erf(float(z) / math.sqrt(2.0))), 6)


def _omics_levels(ensg: str, mid: str) -> dict:
    """Per-line RNA and protein expression level (0-1), joined from the z-score
    parquets that sit beside the predictions file. Pulled on demand via DuckDB."""
    out = {"expression_level": None, "proteomics_level": None}
    if _pred_path is None:
        return out
    base = _pred_path.parent
    sources = {"expression_level": base / "bulk_rna_z.parquet",
               "proteomics_level": base / "bulk_prot_z.parquet"}
    try:
        import duckdb
        con = duckdb.connect()
        for key, path in sources.items():
            if not path.exists():
                continue
            row = con.execute(
                "SELECT z_t FROM read_parquet(?) WHERE gene_id = ? AND lower(model_id) = ?",
                [str(path), ensg, mid.lower()],
            ).fetchone()
            if row:
                out[key] = _z_to_level(row[0])
        con.close()
    except Exception as exc:
        logger.warning("omics levels lookup failed: %s", exc)
    return out


def detail(gene: str, model_id: str) -> dict:
    """Rich per-line detail matching the UI's CellLineDetailApiResponse shape."""
    ensg, sym = _resolve(gene)
    mid = model_id.strip().lower()

    s = _gene_series(ensg, None)
    if s.empty:
        raise ValueError(f"No predictions for gene {gene!r}")

    score = s.get(mid)
    if score is None or pd.isna(score):
        raise ValueError(f"No prediction for {gene!r} + {model_id!r}")

    rank_pos = int((s.dropna() > score).sum()) + 1
    total = int(s.notna().sum())

    name = _cell_name.get(mid) or mid.upper()
    pr = _prediction_row(ensg, mid)
    omics = _omics_levels(ensg, mid)

    def _num(v):
        return float(v) if v is not None and pd.notna(v) else None

    return {
        "gene": sym,
        "ensg": ensg,
        "cell_line": {"model_id": mid.upper(), "name": name},
        "rank": rank_pos,
        "total": total,
        "score": round(float(score), 6),
        "tier": pr.get("confidence_tier") or "unknown",
        "n_layers": int(pr.get("n_layers") or 0),
        "driver_alteration": bool(pr.get("has_driver_alteration") or False),
        "p_mutation": _num(pr.get("p_mutation")),
        "p_fusion": _num(pr.get("p_fusion")),
        "has_cna_alteration": bool(pr.get("has_cna_alteration") or False),
        "expression_level": omics["expression_level"],
        "proteomics_level": omics["proteomics_level"],
        "metadata": _full_metadata(mid),
        "rna_alternatives": [],
    }
