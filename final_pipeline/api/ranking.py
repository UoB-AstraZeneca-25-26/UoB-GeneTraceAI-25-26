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


# Where cold-start S3 downloads land -- Lambda's image filesystem is read-only,
# /tmp is the only writable directory (backed by DATA_BUCKET's EphemeralStorage).
_S3_CACHE_DIR = Path("/tmp") / "genetraceai_data"


def _resolve_data_path(local_path: Path, s3_key: str) -> Path:
    """Returns a local, readable path to a data file. If `local_path` already
    exists (local dev, or any future image that still bakes data in), use it
    directly -- no network. Otherwise download it from S3 (DATA_BUCKET env var,
    key `processed/<s3_key>`) into /tmp once per cold start; every warm
    invocation after that reuses the same container's /tmp, so the download
    cost is paid at most once per execution environment, not per request."""
    if local_path.exists():
        return local_path
    bucket = os.getenv("DATA_BUCKET")
    if not bucket:
        return local_path  # no S3 configured -- let the caller's own "missing file" handling apply
    cached = _S3_CACHE_DIR / s3_key
    if not cached.exists():
        cached.parent.mkdir(parents=True, exist_ok=True)
        import boto3
        key = f"processed/{s3_key}"
        logger.info("downloading s3://%s/%s -> %s", bucket, key, cached)
        boto3.client("s3").download_file(bucket, key, str(cached))
    return cached


def _resolve_rna_neighbours_path() -> Path | None:
    """Same as _resolve_data_path, but for the one genuinely OPTIONAL file:
    RNA-similarity neighbours power a nice-to-have UI panel, not core scoring,
    and load_ranking_data() already degrades gracefully (empty
    rna_alternatives) when this path doesn't exist. _resolve_data_path()
    itself has no such fallback -- a failed S3 download (object not yet
    uploaded, permissions, transient failure, ...) raises straight out of it,
    and since this whole chain runs during the FastAPI lifespan, that
    exception takes down EVERY endpoint's cold start, not just this feature.
    Learned the hard way: rna_neighbours.parquet was deployed in code before
    it existed in the production data bucket, and every request 500'd until
    the object was uploaded. Catching here contains a missing/failed fetch of
    this one optional file to "the panel is empty," instead of "the API is
    down," without changing behaviour for the four required files above --
    those SHOULD crash loudly if actually missing, since the app cannot serve
    correct data at all without them."""
    try:
        return _resolve_data_path(
            _env_path("RNA_NEIGHBOURS_PATH", _PIPELINE.parent / "cell_similarity" / "outputs" / "rna_neighbours.parquet"),
            "rna_neighbours.parquet")
    except Exception as exc:
        logger.warning("RNA neighbours path resolution failed (non-fatal, feature degrades): %s", exc)
        return None


def ensure_loaded() -> None:
    """Idempotent data load. Safe to call from either the HTTP lifespan or the
    Bedrock handler — loads once, no-ops thereafter. Every data file is
    fetched from S3 on a cold start (see _resolve_data_path); local dev reuses
    whatever's already on disk under final_pipeline/outputs|reference without
    touching the network."""
    if is_ready():
        return
    load_ranking_data(
        predictions_path=_resolve_data_path(
            _env_path("PREDICTIONS_PATH", _PIPELINE / "outputs" / "predictions_with_confidence.parquet"),
            "predictions_with_confidence.parquet"),
        gene_lookup_path=_resolve_data_path(
            _env_path("GENE_LOOKUP_PATH", _PIPELINE / "reference" / "gene_lookup.parquet"),
            "gene_lookup.parquet"),
        cell_lookup_path=_resolve_data_path(
            _env_path("CELL_LOOKUP_PATH", _PIPELINE / "reference" / "cell_line_lookup.parquet"),
            "cell_line_lookup.parquet"),
        sample_info_path=_resolve_data_path(
            _env_path("SAMPLE_INFO_PATH", _PIPELINE / "reference" / "sample_info.parquet"),
            "sample_info.parquet"),
        rna_neighbours_path=_resolve_rna_neighbours_path(),
    )
    # The four per-line omics parquets are read lazily, on demand, straight off
    # disk by _omics_levels()/_omics_levels_by_source() (base = _pred_path.parent)
    # -- pre-fetch them here too so the first /gene/detail request doesn't pay a
    # surprise per-request download.
    if _pred_path is not None:
        base = _pred_path.parent
        for fname in ["bulk_rna_z.parquet", "bulk_prot_z.parquet",
                      "bulk_rna_z_by_source.parquet", "bulk_prot_z_by_source.parquet"]:
            _resolve_data_path(base / fname, fname)

_pred: pd.DataFrame | None = None  # index=ensg_id (categorical, sorted), cols: model_id, core_score
_sym_index: dict[str, tuple[str, str]] = {}  # upper(symbol) -> (ensg_id, hgnc_symbol)
_ensg_index: dict[str, tuple[str, str]] = {}  # upper(ensg_id) -> (ensg_id, hgnc_symbol)
_gene_lkp: pd.DataFrame | None = None
_cell_lkp: pd.DataFrame | None = None
_cell_name: dict[str, str] = {}  # lower(model_id) -> cell_line_name
_lineage_map: dict[str, str] = {}  # lower(model_id) -> lineage
_meta_map: dict[str, dict] = {}  # lower(model_id) -> full sample_info row (lineage, subtype, disease, sex, age, ...)
_pred_path: Path | None = None  # kept so /gene/detail can pull the full row on demand
_y_set: set[str] = set()  # ENSG IDs on the Y chromosome (chromosomal_location starts "Y")
_rna_neighbours: dict[str, list[tuple[str, float]]] = {}  # lower(model_id) -> [(neighbour, similarity), ...] sorted by rank
RNA_ALTERNATIVES_TOP_N = 5  # how many similar lines /gene/detail surfaces


def load_ranking_data(predictions_path: Path, gene_lookup_path: Path,
                      cell_lookup_path: Path, sample_info_path: Path,
                      rna_neighbours_path: Path | None = None) -> None:
    global _pred, _sym_index, _ensg_index, _gene_lkp, _cell_lkp, _cell_name
    global _lineage_map, _meta_map, _pred_path, _y_set, _rna_neighbours

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

    _gene_lkp = pd.read_parquet(gene_lookup_path,
                                columns=["ensg_id", "hgnc_symbol", "approved_name",
                                         "chromosomal_location"])
    loc = _gene_lkp["chromosomal_location"].fillna("")
    _y_set = set(_gene_lkp.loc[loc.str.startswith("Y"), "ensg_id"])
    logger.info("Y-linked gene set: %d genes", len(_y_set))
    # Vectorized zip — far faster than iterrows over the gene table.
    for ensg, sym in zip(_gene_lkp["ensg_id"], _gene_lkp["hgnc_symbol"]):
        pair = (ensg, sym)
        _sym_index[sym.upper()] = pair
        _ensg_index[ensg.upper()] = pair

    rrid_map: dict[str, str] = {}
    if cell_lookup_path.exists():
        _cell_lkp = pd.read_parquet(cell_lookup_path, columns=["model_id", "cell_line_name", "rrid"])
        _cell_lkp["model_id"] = _cell_lkp["model_id"].str.lower()
        _cell_name = {
            mid: name
            for mid, name in zip(_cell_lkp["model_id"], _cell_lkp["cell_line_name"])
            if isinstance(name, str)
        }
        rrid_map = {
            mid: rrid
            for mid, rrid in zip(_cell_lkp["model_id"], _cell_lkp["rrid"])
            if isinstance(rrid, str)
        }

    if sample_info_path.exists():
        try:
            info = pd.read_parquet(sample_info_path)
            info["model_id"] = info["model_id"].str.lower()
            for row in info.to_dict("records"):
                mid = row.pop("model_id")
                meta = {k: v for k, v in row.items() if v is not None and v == v}  # v==v excludes NaN
                if "lineage" in meta:
                    _lineage_map[mid] = meta["lineage"]
                if meta:
                    _meta_map[mid] = meta
            logger.info("metadata map loaded: %d lines", len(_meta_map))
        except Exception as exc:
            logger.warning("metadata map load failed: %s", exc)

    # rrid comes from cell_line_lookup, a different source than sample_info --
    # merge it into the same per-line metadata dict every ranking response
    # already returns, so a line with no sample_info row still gets an rrid.
    for mid, rrid in rrid_map.items():
        _meta_map.setdefault(mid, {})["rrid"] = rrid

    # RNA-similarity neighbours (transcriptomic k-NN, precomputed by
    # cell_similarity/) -- powers /gene/detail's "RNA-similar alternatives"
    # panel. Optional: the UI section this feeds simply stays empty if the
    # file is absent, same graceful-degradation as sample_info above.
    if rna_neighbours_path is not None and rna_neighbours_path.exists():
        try:
            nbr = pd.read_parquet(rna_neighbours_path, columns=["model_id", "neighbour", "rank", "similarity"])
            nbr["model_id"] = nbr["model_id"].str.lower()
            nbr["neighbour"] = nbr["neighbour"].str.lower()
            nbr = nbr.sort_values(["model_id", "rank"])
            for mid, grp in nbr.groupby("model_id", sort=False):
                _rna_neighbours[mid] = list(zip(grp["neighbour"], grp["similarity"].astype(float)))
            logger.info("RNA-similarity neighbours loaded: %d lines", len(_rna_neighbours))
        except Exception as exc:
            logger.warning("RNA neighbours load failed: %s", exc)


def is_ready() -> bool:
    return _pred is not None and _gene_lkp is not None


def list_genes() -> dict:
    """The full recognised gene universe (symbol + ENSG + approved name), for
    client-side autocomplete/validation and the Reference lookup page -- so
    the UI can catch an unrecognised symbol before submitting, instead of
    round-tripping to find out. Served from the same in-memory table
    _resolve() already uses; no extra I/O."""
    genes = [
        {"symbol": sym, "ensg": ensg, "name": name if isinstance(name, str) else None}
        for ensg, sym, name in zip(_gene_lkp["ensg_id"], _gene_lkp["hgnc_symbol"], _gene_lkp["approved_name"])
    ]
    return {"genes": genes}


def list_cell_lines() -> dict:
    """The full recognised cell-line universe (model_id + display name +
    lineage + primary disease), for the Reference lookup page -- one bulk
    export computed from the same _cell_lkp / _meta_map tables already loaded
    at startup, no extra I/O. Lines with no sample_info row (lineage-less,
    see the lineage-grouped ranking's unassigned-count note) simply omit
    lineage/primary_disease -- still resolvable by id/name."""
    if _cell_lkp is None:
        return {"cell_lines": []}
    lines = []
    for mid, name in zip(_cell_lkp["model_id"], _cell_lkp["cell_line_name"]):
        meta = _meta_map.get(mid, {})
        lines.append({
            "model_id": mid.upper(),
            "name": name if isinstance(name, str) else None,
            "lineage": meta.get("lineage"),
            "primary_disease": meta.get("primary_disease"),
        })
    return {"cell_lines": lines}


def _resolve(query: str) -> tuple[str, str]:
    q = query.strip().upper()
    if q in _sym_index:
        return _sym_index[q]
    if q in _ensg_index:
        return _ensg_index[q]
    raise ValueError(f"Gene not found: {query!r}")


def _resolve_unique(queries: list[str]) -> list[tuple[str, str]]:
    """_resolve() every query, then drop later duplicates by canonical symbol
    (first occurrence wins). Two different query strings can resolve to the
    SAME gene -- typed twice, or one as a symbol ("BRAF") and one as its own
    Ensembl ID ("ENSG00000157764") -- and every multi-gene endpoint below
    builds a wide DataFrame with one column per resolved symbol; two columns
    that end up with the identical name crash the wide-merge/selection step
    with a pandas KeyError surfaced to the client as a raw 500. Deduping here,
    once, before any of that column-building happens, closes that off for
    every caller instead of guarding each merge site separately."""
    seen: set[str] = set()
    out: list[tuple[str, str]] = []
    for q in queries:
        ensg, sym = _resolve(q)
        if sym in seen:
            continue
        seen.add(sym)
        out.append((ensg, sym))
    return out


def _lineage_model_ids(terms: list[str]) -> set[str] | None:
    """model_ids whose lineage or lineage_subtype contains any of `terms`
    (case-insensitive substring) -- served from the in-memory _meta_map, no
    file I/O per call."""
    if not terms:
        return None
    needles = [t.lower() for t in terms]
    matched = {
        mid for mid, meta in _meta_map.items()
        if any(
            n in meta.get("lineage", "").lower() or n in meta.get("lineage_subtype", "").lower()
            for n in needles
        )
    }
    return matched or None


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


def _sex_of(model_id: str) -> str:
    """Return the sex of a line from _meta_map: 'male', 'female', or 'unknown'.
    Values in sample_info.parquet are already lowercase ('male'/'female').
    Missing entries (line absent from sample_info) → 'unknown'.
    Matches cli.py's _line_sex() / fillna('unknown') convention exactly."""
    return _meta_map.get(model_id.lower(), {}).get("sex", "unknown")


def _apply_sex_guard_wide(wide: pd.DataFrame, y_set: set[str],
                          resolved: list[tuple[str, str]],
                          syms: list[str]) -> tuple[pd.DataFrame, int]:
    """For each Y-linked gene in the query, set its score column to NaN for
    female and unknown-sex lines.  Returns (wide, n_nulled).

    Mirrors cli.py's cmd_genes sex-guard block exactly: NaN propagates through
    min(skipna=False) so those lines get joint_score=NaN and are excluded by
    the scoreable filter before the floor gate."""
    ensgs = [e for e, _ in resolved]
    y_in_query = [e for e in ensgs if e in y_set]
    if not y_in_query:
        return wide, 0

    y_syms = [syms[ensgs.index(e)] for e in y_in_query]
    non_male = wide["model_id"].map(_sex_of).isin(["female", "unknown"])
    n_nulled = int(non_male.sum())
    for sym in y_syms:
        wide = wide.copy()
        wide.loc[non_male, sym] = float("nan")
    logger.info("sex guard: %d female/unknown-sex rows nulled for Y-linked genes %s",
                n_nulled, y_syms)
    return wide, n_nulled


def rank(gene_syms: list[str], lineages: list[str],
         floor: float, top_n: int) -> dict:
    resolved = _resolve_unique(gene_syms)
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

    # Sex guard: Y-linked genes are structurally absent in female/unknown-sex
    # lines — null their score columns before joint_score so those lines get
    # joint_score=NaN and are excluded by the scoreable filter, not penalised
    # as failed co-selection.  Mirrors cli.py cmd_genes guard exactly.
    wide, _ = _apply_sex_guard_wide(wide, _y_set, resolved, syms)

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
    resolved_b = _resolve_unique(gene_bs)
    syms_b = [s for _, s in resolved_b]

    lineage_ids = _lineage_model_ids(lineages)

    merged = _gene_series(ensg_a, lineage_ids).rename("score_a").reset_index()
    for ensg_b, sym_b in resolved_b:
        b = _gene_series(ensg_b, lineage_ids).rename(sym_b).reset_index()
        merged = merged.merge(b, on="model_id", how="inner")

    # Sex guard: for each Y-linked exclusion gene, null its score column for
    # female/unknown-sex lines before computing selectivity.  Mirrors cli.py
    # cmd_exclude guard applied per exclusion gene.
    for ensg_b, sym_b in resolved_b:
        if ensg_b in _y_set:
            non_male = merged["model_id"].map(_sex_of).isin(["female", "unknown"])
            n_nulled = int(non_male.sum())
            merged = merged.copy()
            merged.loc[non_male, sym_b] = float("nan")
            logger.info("sex guard (exclude_many): %d female/unknown-sex rows nulled "
                        "for Y-linked %s", n_nulled, sym_b)

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

    # Sex guard: if gene_b is Y-linked, score_b is structurally absent in
    # female/unknown-sex lines — set it to NaN so selectivity=NaN and those
    # lines are excluded, not rewarded for structural absence of the gene.
    # Mirrors cli.py cmd_exclude guard exactly.
    if ensg_b in _y_set:
        non_male = merged["model_id"].map(_sex_of).isin(["female", "unknown"])
        n_nulled = int(non_male.sum())
        merged = merged.copy()
        merged.loc[non_male, "score_b"] = float("nan")
        logger.info("sex guard (exclude): %d female/unknown-sex rows nulled for Y-linked %s",
                    n_nulled, sym_b)

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


def _group_by_lineage(metric: pd.Series, top_lineages: int, per_lineage: int) -> dict:
    """Shared lineage-grouping core, metric-agnostic: given ANY (model_id ->
    score) Series pooled across all lineages, keep the top `per_lineage` rows
    in each lineage, then keep the top `top_lineages` lineages (ranked by
    their own best candidate's score). Sizing is fixed by this rule -- there
    is no top_n / floor here. Used by single-gene, multi-gene, and both
    selectivity modes -- each builds its own metric Series (core_score /
    joint_score / selectivity) and calls this once.

    rank_within_lineage is a clean 1..per_lineage ordinal within the kept
    lineage's scoreable population, ties broken by model_id ascending -- the
    same tie rule core_score.py uses for stratum_rank (not reused directly:
    see the data-layer note below), computed fresh so a "top 5" is always
    exactly 5 distinct rows.

    rank_global reproduces the pooled, cross-lineage rank formula already
    live on /gene/detail (ties share a rank: 1 + count of strictly higher
    scores) so it means the same thing here as it already does there.

    Data-layer note: core_score.parquet's `stratum_rank` column already
    computes (gene, lineage) rank but is unread anywhere live, and is
    verified NULL for ~768k otherwise-scoreable rows (67 cell lines absent
    from celllineselector.db's sample_info, so lineage-less) -- deliberately
    not reused here; both rank fields are computed fresh, in-memory, per
    request, the same way /gene/detail already computes its rank.

    Returns {lineages_returned, total_scoreable, lineage_unassigned_count,
    kept: DataFrame[model_id, metric, lineage, rank_within_lineage,
    rank_global]} -- callers attach whichever mode-specific fields (scores
    dict, exclusion_scores, ...) belong on top of this shared scaffold.
    """
    s = metric.dropna()
    if s.empty:
        return {"lineages_returned": 0, "total_scoreable": 0,
                "lineage_unassigned_count": 0, "kept": None}

    df = s.rename("metric").reset_index()  # columns: model_id, metric
    df["lineage"] = df["model_id"].str.lower().map(_lineage_map)

    unassigned = df["lineage"].isna()
    lineage_unassigned_count = int(unassigned.sum())
    df = df[~unassigned].copy()

    # rank_global: pooled across every scoreable line (lineage-blind) -- same
    # "ties share a rank" formula detail() already uses.
    global_rank = s.rank(ascending=False, method="min")
    df["rank_global"] = df["model_id"].map(global_rank).astype(int)

    # Deterministic order for within-lineage ranking: metric desc, then
    # model_id asc on ties (core_score.py's own stratum_rank tie rule).
    df = df.sort_values(["lineage", "metric", "model_id"],
                        ascending=[True, False, True])
    df["rank_within_lineage"] = (
        df.groupby("lineage")["metric"].rank(ascending=False, method="first")
    ).astype(int)

    # Each lineage's top candidate score decides lineage order.
    lineage_best = df.groupby("lineage")["metric"].max().sort_values(ascending=False)
    kept_lineages = lineage_best.head(top_lineages).index.tolist()
    lineage_order = {lin: i for i, lin in enumerate(kept_lineages)}

    kept = df[df["lineage"].isin(kept_lineages) & (df["rank_within_lineage"] <= per_lineage)].copy()
    kept["_lineage_order"] = kept["lineage"].map(lineage_order)
    kept = kept.sort_values(["_lineage_order", "rank_within_lineage"])

    return {
        "lineages_returned": len(kept_lineages),
        "total_scoreable": int(len(s)),
        "lineage_unassigned_count": lineage_unassigned_count,
        "kept": kept,
    }


def rank_by_lineage(gene: str, top_lineages: int = 10, per_lineage: int = 5) -> dict:
    """Single-gene lineage ranking: metric is core_score for the one gene."""
    ensg, sym = _resolve(gene)
    grouped = _group_by_lineage(_gene_series(ensg, None), top_lineages, per_lineage)

    kept = grouped.pop("kept")
    if kept is None:
        return {"gene": sym, **grouped, "lines": []}

    if _cell_lkp is not None:
        kept = kept.merge(_cell_lkp, on="model_id", how="left")

    lines = []
    for _, row in kept.iterrows():
        raw_name = row.get("cell_line_name")
        lines.append({
            "model_id": row["model_id"].upper(),
            "name": raw_name if isinstance(raw_name, str) else None,
            "core_score": round(float(row["metric"]), 6),
            "lineage": row["lineage"],
            "rank_within_lineage": int(row["rank_within_lineage"]),
            "rank_global": int(row["rank_global"]),
            "metadata": _meta(row["model_id"]),
        })

    return {"gene": sym, **grouped, "lines": lines}


def rank_by_lineage_multi(genes: list[str], top_lineages: int = 10, per_lineage: int = 5) -> dict:
    """Multi-gene lineage ranking: metric is the joint score (min across
    genes, same weakest-link rule as rank()'s N-gene path)."""
    resolved = _resolve_unique(genes)
    syms = [s for _, s in resolved]

    wide = None
    for ensg, sym in resolved:
        s = _gene_series(ensg, None).rename(sym).reset_index()
        wide = s if wide is None else wide.merge(s, on="model_id", how="inner")

    if wide is None or wide.empty:
        return {"genes": syms, "lineages_returned": 0, "total_scoreable": 0,
                "lineage_unassigned_count": 0, "lines": []}

    wide["joint_score"] = wide[syms].min(axis=1, skipna=False)
    metric = wide.set_index("model_id")["joint_score"]

    grouped = _group_by_lineage(metric, top_lineages, per_lineage)
    kept = grouped.pop("kept")
    if kept is None:
        return {"genes": syms, **grouped, "lines": []}

    kept = kept.merge(wide[["model_id", *syms]], on="model_id", how="left")
    if _cell_lkp is not None:
        kept = kept.merge(_cell_lkp, on="model_id", how="left")

    lines = []
    for _, row in kept.iterrows():
        raw_name = row.get("cell_line_name")
        gene_scores = {sym: (round(float(row[sym]), 6) if pd.notna(row[sym]) else None) for sym in syms}
        limiting = min((s for s in syms if row[s] == row["metric"]), default=None)
        lines.append({
            "model_id": row["model_id"].upper(),
            "name": raw_name if isinstance(raw_name, str) else None,
            "joint_score": round(float(row["metric"]), 6),
            "scores": gene_scores,
            "limiting_gene": limiting,
            "lineage": row["lineage"],
            "rank_within_lineage": int(row["rank_within_lineage"]),
            "rank_global": int(row["rank_global"]),
            "metadata": _meta(row["model_id"]),
        })

    return {"genes": syms, **grouped, "lines": lines}


def exclude_by_lineage(gene_a: str, gene_b: str, top_lineages: int = 10, per_lineage: int = 5) -> dict:
    """Single-exclusion selectivity lineage ranking: metric is selectivity
    (score_a * (1 - score_b)), same formula as exclude()."""
    ensg_a, sym_a = _resolve(gene_a)
    ensg_b, sym_b = _resolve(gene_b)

    a = _gene_series(ensg_a, None).rename("score_a").reset_index()
    b = _gene_series(ensg_b, None).rename("score_b").reset_index()
    merged = a.merge(b, on="model_id", how="inner")
    merged["selectivity"] = merged["score_a"] * (1.0 - merged["score_b"])
    metric = merged.set_index("model_id")["selectivity"]

    grouped = _group_by_lineage(metric, top_lineages, per_lineage)
    kept = grouped.pop("kept")
    if kept is None:
        return {"gene_a": sym_a, "gene_b": sym_b, **grouped, "lines": []}

    kept = kept.merge(merged[["model_id", "score_a", "score_b"]], on="model_id", how="left")
    if _cell_lkp is not None:
        kept = kept.merge(_cell_lkp, on="model_id", how="left")

    lines = []
    for _, row in kept.iterrows():
        raw_name = row.get("cell_line_name")
        lines.append({
            "model_id": row["model_id"].upper(),
            "name": raw_name if isinstance(raw_name, str) else None,
            "score_a": round(float(row["score_a"]), 6),
            "score_b": round(float(row["score_b"]), 6),
            "selectivity": round(float(row["metric"]), 6),
            "lineage": row["lineage"],
            "rank_within_lineage": int(row["rank_within_lineage"]),
            "rank_global": int(row["rank_global"]),
            "metadata": _meta(row["model_id"]),
        })

    return {"gene_a": sym_a, "gene_b": sym_b, **grouped, "lines": lines}


def exclude_many_by_lineage(gene_a: str, gene_bs: list[str], top_lineages: int = 10,
                            per_lineage: int = 5) -> dict:
    """Multi-exclusion selectivity lineage ranking: metric is selectivity
    (score_a * Pi(1 - score_bi)), same formula as exclude_many()."""
    ensg_a, sym_a = _resolve(gene_a)
    resolved_b = _resolve_unique(gene_bs)
    syms_b = [s for _, s in resolved_b]

    merged = _gene_series(ensg_a, None).rename("score_a").reset_index()
    for ensg_b, sym_b in resolved_b:
        b = _gene_series(ensg_b, None).rename(sym_b).reset_index()
        merged = merged.merge(b, on="model_id", how="inner")

    selectivity = merged["score_a"].astype("float64").copy()
    for sym_b in syms_b:
        selectivity = selectivity * (1.0 - merged[sym_b])
    merged["selectivity"] = selectivity
    metric = merged.set_index("model_id")["selectivity"]

    grouped = _group_by_lineage(metric, top_lineages, per_lineage)
    kept = grouped.pop("kept")
    if kept is None:
        return {"gene_a": sym_a, "excluded_genes": syms_b, **grouped, "lines": []}

    kept = kept.merge(merged[["model_id", "score_a", *syms_b]], on="model_id", how="left")
    if _cell_lkp is not None:
        kept = kept.merge(_cell_lkp, on="model_id", how="left")

    lines = []
    for _, row in kept.iterrows():
        raw_name = row.get("cell_line_name")
        lines.append({
            "model_id": row["model_id"].upper(),
            "name": raw_name if isinstance(raw_name, str) else None,
            "score_a": round(float(row["score_a"]), 6),
            "exclusion_scores": {s: round(float(row[s]), 6) for s in syms_b},
            "selectivity": round(float(row["metric"]), 6),
            "lineage": row["lineage"],
            "rank_within_lineage": int(row["rank_within_lineage"]),
            "rank_global": int(row["rank_global"]),
            "metadata": _meta(row["model_id"]),
        })

    return {"gene_a": sym_a, "excluded_genes": syms_b, **grouped, "lines": lines}


def rank_exclude_many(gene_as: list[str], gene_bs: list[str], lineages: list[str],
                      floor: float, top_n: int) -> dict:
    """Combines rank()'s N-target weakest-link joint score with exclude_many()'s
    N-exclusion selectivity penalty: joint_score = min(target scores), gated by
    floor exactly as rank()'s N-gene path does; combined_score = joint_score *
    Pi(1 - score_bi) for every excluded gene, same multiplicative penalty as
    exclude_many(). Requires 2+ target genes -- for one target use
    exclude()/exclude_many() instead."""
    resolved_a = _resolve_unique(gene_as)
    syms_a = [s for _, s in resolved_a]
    resolved_b = _resolve_unique(gene_bs)
    syms_b = [s for _, s in resolved_b]

    lineage_ids = _lineage_model_ids(lineages)

    wide = None
    for ensg, sym in resolved_a:
        s = _gene_series(ensg, lineage_ids).rename(sym).reset_index()
        wide = s if wide is None else wide.merge(s, on="model_id", how="inner")

    if wide is None or wide.empty:
        return {"genes": syms_a, "excluded_genes": syms_b, "lineage": lineages,
                "floor": floor, "total_passing": 0, "lines": [],
                "lineage_distribution": {}}

    wide["joint_score"] = wide[syms_a].min(axis=1, skipna=False)
    scoreable = wide["joint_score"].notna()
    wide.loc[scoreable, "limiting_gene"] = wide.loc[scoreable, syms_a].idxmin(axis=1)

    for ensg_b, sym_b in resolved_b:
        b = _gene_series(ensg_b, lineage_ids).rename(sym_b).reset_index()
        wide = wide.merge(b, on="model_id", how="inner")

    passing = wide[scoreable & (wide["joint_score"] >= floor)].copy()
    combined = passing["joint_score"].astype("float64").copy()
    for sym_b in syms_b:
        combined = combined * (1.0 - passing[sym_b])
    passing["combined_score"] = combined

    shortlist = passing.sort_values("combined_score", ascending=False, na_position="last").head(top_n)

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
            "scores": {s: round(float(row[s]), 6) for s in syms_a if pd.notna(row[s])},
            "exclusion_scores": {s: round(float(row[s]), 6) for s in syms_b},
            "combined_score": round(float(row["combined_score"]), 6),
            "metadata": _meta(row["model_id"]),
        })

    return {"genes": syms_a, "excluded_genes": syms_b, "lineage": lineages,
            "floor": floor, "total_passing": len(passing), "lines": lines,
            "lineage_distribution": _lineage_dist(passing["model_id"].tolist())}


def rank_exclude_many_by_lineage(gene_as: list[str], gene_bs: list[str],
                                 top_lineages: int = 10, per_lineage: int = 5) -> dict:
    """Lineage-grouped counterpart of rank_exclude_many(): metric is
    combined_score, same formula, no floor gate (matches every other
    *_by_lineage function -- sizing comes from _group_by_lineage alone)."""
    resolved_a = _resolve_unique(gene_as)
    syms_a = [s for _, s in resolved_a]
    resolved_b = _resolve_unique(gene_bs)
    syms_b = [s for _, s in resolved_b]

    wide = None
    for ensg, sym in resolved_a:
        s = _gene_series(ensg, None).rename(sym).reset_index()
        wide = s if wide is None else wide.merge(s, on="model_id", how="inner")

    if wide is None or wide.empty:
        return {"genes": syms_a, "excluded_genes": syms_b, "lineages_returned": 0,
                "total_scoreable": 0, "lineage_unassigned_count": 0, "lines": []}

    wide["joint_score"] = wide[syms_a].min(axis=1, skipna=False)

    for ensg_b, sym_b in resolved_b:
        b = _gene_series(ensg_b, None).rename(sym_b).reset_index()
        wide = wide.merge(b, on="model_id", how="inner")

    combined = wide["joint_score"].astype("float64").copy()
    for sym_b in syms_b:
        combined = combined * (1.0 - wide[sym_b])
    wide["combined_score"] = combined
    metric = wide.set_index("model_id")["combined_score"]

    grouped = _group_by_lineage(metric, top_lineages, per_lineage)
    kept = grouped.pop("kept")
    if kept is None:
        return {"genes": syms_a, "excluded_genes": syms_b, **grouped, "lines": []}

    kept = kept.merge(wide[["model_id", "joint_score", *syms_a, *syms_b]], on="model_id", how="left")
    if _cell_lkp is not None:
        kept = kept.merge(_cell_lkp, on="model_id", how="left")

    lines = []
    for _, row in kept.iterrows():
        raw_name = row.get("cell_line_name")
        gene_scores = {sym: (round(float(row[sym]), 6) if pd.notna(row[sym]) else None) for sym in syms_a}
        limiting = min((s for s in syms_a if row[s] == row["joint_score"]), default=None)
        lines.append({
            "model_id": row["model_id"].upper(),
            "name": raw_name if isinstance(raw_name, str) else None,
            "joint_score": round(float(row["joint_score"]), 6),
            "scores": gene_scores,
            "limiting_gene": limiting,
            "exclusion_scores": {s: round(float(row[s]), 6) for s in syms_b},
            "combined_score": round(float(row["metric"]), 6),
            "lineage": row["lineage"],
            "rank_within_lineage": int(row["rank_within_lineage"]),
            "rank_global": int(row["rank_global"]),
            "metadata": _meta(row["model_id"]),
        })

    return {"genes": syms_a, "excluded_genes": syms_b, **grouped, "lines": lines}


# sample_info columns surfaced in the detail metadata block, in the order the
# UI labels them. Kept flat (single dict) because that's what the UI expects.
_DETAIL_META_COLS = [
    "lineage", "lineage_subtype", "primary_disease", "subtype",
    "cellosaurus_ncit_disease", "primary_or_metastasis",
    "sample_collection_site", "default_growth_pattern", "sex", "age",
]


def _full_metadata(mid: str) -> dict:
    """Flat metadata dict for one line, keyed as the UI's detail view expects.
    Served from the in-memory _meta_map, no file I/O per call."""
    meta = _meta_map.get(mid.lower(), {})
    return {k: meta[k] for k in _DETAIL_META_COLS if k in meta}


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


def _omics_levels_by_source(ensg: str, mid: str) -> dict:
    """Per-source expression/proteomics level (0-1 each) -- DepMap/HPA/GEO for
    RNA, ProCan/CCLE for protein. Purely additive: reads two new export files
    (bulk_rna_z_by_source.parquet, bulk_prot_z_by_source.parquet) that sit
    beside the existing combined ones and were built without touching any
    scoring script or re-running the pipeline -- the per-source z-scores were
    already computed as an internal step of the existing RNA/protein scorers,
    just never persisted before now. Does not read or affect
    expression_level/proteomics_level (the existing combined fields) at all."""
    out = {"expression_by_source": {}, "proteomics_by_source": {}}
    if _pred_path is None:
        return out
    base = _pred_path.parent
    try:
        import duckdb
        con = duckdb.connect()
        rna_path = base / "bulk_rna_z_by_source.parquet"
        if rna_path.exists():
            row = con.execute(
                "SELECT z_depmap, z_hpa_rna, z_geo FROM read_parquet(?) "
                "WHERE gene_id = ? AND lower(model_id) = ?",
                [str(rna_path), ensg, mid.lower()],
            ).fetchone()
            if row:
                out["expression_by_source"] = {
                    "DepMap": _z_to_level(row[0]),
                    "HPA": _z_to_level(row[1]),
                    "GEO": _z_to_level(row[2]),
                }
        prot_path = base / "bulk_prot_z_by_source.parquet"
        if prot_path.exists():
            row = con.execute(
                "SELECT z_procan, z_ccle FROM read_parquet(?) "
                "WHERE gene_id = ? AND lower(model_id) = ?",
                [str(prot_path), ensg, mid.lower()],
            ).fetchone()
            if row:
                out["proteomics_by_source"] = {
                    "ProCan": _z_to_level(row[0]),
                    "CCLE": _z_to_level(row[1]),
                }
        con.close()
    except Exception as exc:
        logger.warning("per-source omics levels lookup failed: %s", exc)
    return out


def _rna_alternatives(mid: str, scores_by_gene: dict[str, pd.Series],
                      exclude_scores_by_gene: dict[str, pd.Series] | None = None) -> list[dict]:
    """Top RNA-similar lines for `mid` that also have real (non-NaN) evidence
    for EVERY gene the original query actually needed -- not just the one gene
    currently on screen. `scores_by_gene` covers every TARGET gene (the gene
    being viewed plus any other targets from a multi/jointSelectivity query);
    `exclude_scores_by_gene` covers exclusion genes from a
    selectivity/jointSelectivity query, if any.

    A neighbour missing evidence for any target gene would be a dead-end click
    ("no prediction for ...") for the "click to inspect the same gene in that
    line" affordance the UI promises, so such neighbours are filtered out here
    rather than left for the frontend to discover after a failed fetch. The
    same reasoning extends to target genes beyond the one on screen: showing
    an "alternative" that can't actually stand in for the full multi-gene
    query it was suggested from would be misleading.

    DESIGN DECISION NEEDING SIGN-OFF (Fiona/Daniel): exclusion-gene validity
    here means only "has a real score" (non-NaN), NOT "scores low enough to
    still satisfy the selectivity criterion." A neighbour could have a real
    but HIGH exclusion-gene score and still pass this filter. Requiring an
    actual low-score threshold for exclusion genes was NOT implemented --
    doing so would duplicate the selectivity floor logic from rank()/exclude()
    without an agreed threshold to reuse, so this stops at "has data" pending
    that decision."""
    def _missing(s: pd.Series, neighbour: str) -> bool:
        score = s.get(neighbour)
        return score is None or pd.isna(score)

    candidates = _rna_neighbours.get(mid, [])
    out: list[dict] = []
    for neighbour, similarity in candidates:
        if any(_missing(s, neighbour) for s in scores_by_gene.values()):
            continue
        if exclude_scores_by_gene and any(_missing(s, neighbour) for s in exclude_scores_by_gene.values()):
            continue
        meta = _meta_map.get(neighbour, {})
        out.append({
            "model_id": neighbour.upper(),
            "name": _cell_name.get(neighbour) or neighbour.upper(),
            "similarity": round(float(similarity), 6),
            "lineage": meta.get("lineage"),
            "primary_disease": meta.get("primary_disease"),
        })
        if len(out) >= RNA_ALTERNATIVES_TOP_N:
            break
    return out


def detail(gene: str, model_id: str, other_genes: list[str] | None = None,
          exclude_genes: list[str] | None = None) -> dict:
    """Rich per-line detail matching the UI's CellLineDetailApiResponse shape.

    other_genes / exclude_genes: the REST of the original multi-gene query's
    target/exclusion genes (if any) -- passed through only to scope
    rna_alternatives correctly (see _rna_alternatives). They don't affect any
    other field in this response, which stays single-gene (`gene`/`ensg`)
    exactly as before. A gene symbol here that fails to resolve is silently
    skipped rather than raising -- it's supplementary context for narrowing
    alternatives, not something this endpoint is itself being asked to score."""
    ensg, sym = _resolve(gene)
    mid = model_id.strip().lower()

    s = _gene_series(ensg, None)
    if s.empty:
        raise ValueError(f"No predictions for gene {gene!r}")

    score = s.get(mid)
    if score is None or pd.isna(score):
        raise ValueError(f"No prediction for {gene!r} + {model_id!r}")

    scores_by_gene: dict[str, pd.Series] = {sym: s}
    for g in (other_genes or []):
        try:
            g_ensg, g_sym = _resolve(g)
        except ValueError:
            continue
        scores_by_gene[g_sym] = _gene_series(g_ensg, None)

    exclude_scores_by_gene: dict[str, pd.Series] = {}
    for g in (exclude_genes or []):
        try:
            g_ensg, g_sym = _resolve(g)
        except ValueError:
            continue
        exclude_scores_by_gene[g_sym] = _gene_series(g_ensg, None)

    rank_pos = int((s.dropna() > score).sum()) + 1
    total = int(s.notna().sum())

    name = _cell_name.get(mid) or mid.upper()
    pr = _prediction_row(ensg, mid)
    omics = _omics_levels(ensg, mid)
    omics_by_source = _omics_levels_by_source(ensg, mid)

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
        # 0.5 matches Scoring/driver_routing.py's MUT_DRIVER_THRESHOLD/FUS_DRIVER_THRESHOLD --
        # per-layer split of the combined has_driver_alteration flag, so the UI can attribute
        # the "driver" badge to the layer that actually earned it (mutation vs fusion vs CNA).
        "mutation_driver": bool((pr.get("p_mutation") or 0.0) >= 0.5),
        "fusion_driver": bool((pr.get("p_fusion") or 0.0) >= 0.5),
        "has_cna_alteration": bool(pr.get("has_cna_alteration") or False),
        "expression_level": omics["expression_level"],
        "proteomics_level": omics["proteomics_level"],
        "expression_by_source": omics_by_source["expression_by_source"],
        "proteomics_by_source": omics_by_source["proteomics_by_source"],
        "metadata": _full_metadata(mid),
        "rna_alternatives": _rna_alternatives(mid, scores_by_gene, exclude_scores_by_gene or None),
    }
