"""Routes for the architecture ranking API."""

from __future__ import annotations

import asyncio
from typing import Optional

from fastapi import APIRouter, HTTPException, Query

from .ranking import detail as _detail
from .ranking import exclude as _exclude
from .ranking import exclude_by_lineage as _exclude_by_lineage
from .ranking import exclude_many as _exclude_many
from .ranking import exclude_many_by_lineage as _exclude_many_by_lineage
from .ranking import is_ready, list_cell_lines, list_genes, rank
from .ranking import rank_by_lineage as _rank_by_lineage
from .ranking import rank_by_lineage_multi as _rank_by_lineage_multi
from .ranking import rank_exclude_many as _rank_exclude_many
from .ranking import rank_exclude_many_by_lineage as _rank_exclude_many_by_lineage
from .schemas import (
    CellLineListResponse,
    DetailResponse,
    ExcludeLineageRankedResponse,
    ExcludeManyLineageRankedResponse,
    ExcludeManyResponse,
    ExcludeResponse,
    GeneListResponse,
    HealthResponse,
    LineageRankedResponse,
    MultiLineageRankedResponse,
    RankExcludeManyLineageRankedResponse,
    RankExcludeManyResponse,
    RankRequest,
    RankResponse,
)

router = APIRouter()


def _guard():
    """
    Reject requests before the ranking data has finished loading.

    Parameters
    ----------
    None

    Returns
    -------
    None

    Notes
    -----
    Raises ``HTTPException(503)`` if :func:`ranking.is_ready` is False.
    Called at the top of every endpoint below except ``/health``, so a
    client polling readiness during Lambda cold start gets a clean 503
    instead of a confusing crash from code that assumes data is present.
    """
    if not is_ready():
        raise HTTPException(status_code=503, detail="Ranking data not loaded.")


@router.post("/v1/rank", response_model=RankResponse)
async def rank_endpoint(req: RankRequest) -> RankResponse:
    """
    Rank cell lines for one or more genes (POST body variant).

    Parameters
    ----------
    req : RankRequest
        Request body with ``genes`` (list of symbols/ENSG IDs),
        ``lineage`` (optional filter list), ``floor`` (minimum score to
        include), and ``top_n`` (result cap).

    Returns
    -------
    RankResponse
        Ranked cell lines, per :func:`ranking.rank`.

    Notes
    -----
    Runs the (synchronous, CPU-bound) ranking call in the default
    executor via ``run_in_executor`` so it doesn't block the event
    loop. Raises ``HTTPException(422)`` if ``ranking.rank`` raises
    ``ValueError`` (e.g. an unrecognised gene).
    """
    _guard()
    try:
        loop = asyncio.get_event_loop()
        result = await loop.run_in_executor(
            None, lambda: rank(req.genes, req.lineage, req.floor, req.top_n)
        )
    except ValueError as exc:
        raise HTTPException(status_code=422, detail=str(exc))
    return RankResponse(**result)


@router.get("/gene", response_model=RankResponse)
async def gene_endpoint(
    gene: str = Query(..., description="Gene symbol or ENSG ID"),
    lineage: list[str] = Query(default=[], description="Lineage filter (repeatable)"),
    top_n: int = Query(default=30, ge=1, le=200),
) -> RankResponse:
    """
    Rank cell lines for a single gene (GET, query-string variant).

    Parameters
    ----------
    gene : str
        Gene symbol or ENSG ID.
    lineage : list of str, optional
        Repeatable query param (``?lineage=lung&lineage=breast``) to
        restrict results to these lineages. Empty means no filter.
    top_n : int, default 30
        Maximum number of ranked lines to return (1-200).

    Returns
    -------
    RankResponse
        Ranked cell lines for this one gene, floor fixed at 0.0 (every
        scored line is eligible; only ``top_n`` caps the result size).

    Notes
    -----
    Raises ``HTTPException(422)`` on an unrecognised gene.
    """
    _guard()
    try:
        loop = asyncio.get_event_loop()
        result = await loop.run_in_executor(
            None, lambda: rank([gene], lineage, 0.0, top_n)
        )
    except ValueError as exc:
        raise HTTPException(status_code=422, detail=str(exc))
    return RankResponse(**result)


@router.get("/gene/by-lineage", response_model=LineageRankedResponse)
async def gene_by_lineage_endpoint(
    gene: str = Query(..., description="Gene symbol or ENSG ID"),
) -> LineageRankedResponse:
    """
    Rank cell lines for one gene, grouped by lineage.

    Parameters
    ----------
    gene : str
        Gene symbol or ENSG ID.

    Returns
    -------
    LineageRankedResponse
        Per-lineage top lines for this gene, per
        :func:`ranking.rank_by_lineage`. Powers the "best line per
        tissue type" view rather than one flat ranked list.

    Notes
    -----
    Raises ``HTTPException(422)`` on an unrecognised gene.
    """
    _guard()
    try:
        loop = asyncio.get_event_loop()
        result = await loop.run_in_executor(
            None, lambda: _rank_by_lineage(gene)
        )
    except ValueError as exc:
        raise HTTPException(status_code=422, detail=str(exc))
    return LineageRankedResponse(**result)


@router.get("/genes", response_model=RankResponse)
async def genes_endpoint(
    genes: str = Query(..., description="Comma-separated gene symbols, e.g. BRAF,KRAS"),
    lineage: list[str] = Query(default=[], description="Lineage filter (repeatable)"),
    floor: float = Query(default=0.50, ge=0.0, le=1.0),
    top_n: int = Query(default=30, ge=1, le=200),
) -> RankResponse:
    """
    Rank cell lines that satisfy ALL of several genes at once (joint query).

    Parameters
    ----------
    genes : str
        Comma-separated gene symbols/IDs, e.g. ``"BRAF,KRAS"``. Must
        contain at least 2 after stripping whitespace/empties.
    lineage : list of str, optional
        Repeatable lineage filter, as in :func:`gene_endpoint`.
    floor : float, default 0.50
        Minimum per-gene score a line must clear on every gene to be
        included (unlike the single-gene endpoint, which fixes this
        at 0.0).
    top_n : int, default 30
        Maximum number of ranked lines to return (1-200).

    Returns
    -------
    RankResponse
        Lines ranked by their combined score across all requested
        genes, per :func:`ranking.rank`.

    Notes
    -----
    Raises ``HTTPException(422)`` if fewer than 2 genes are given after
    parsing, or if any gene is unrecognised.
    """
    _guard()
    gene_list = [g.strip() for g in genes.split(",") if g.strip()]
    if len(gene_list) < 2:
        raise HTTPException(status_code=422, detail="Provide at least 2 comma-separated genes.")
    try:
        loop = asyncio.get_event_loop()
        result = await loop.run_in_executor(
            None, lambda: rank(gene_list, lineage, floor, top_n)
        )
    except ValueError as exc:
        raise HTTPException(status_code=422, detail=str(exc))
    return RankResponse(**result)


@router.get("/genes/by-lineage", response_model=MultiLineageRankedResponse)
async def genes_by_lineage_endpoint(
    genes: str = Query(..., description="Comma-separated gene symbols, e.g. BRAF,KRAS"),
) -> MultiLineageRankedResponse:
    """
    Rank cell lines for several genes at once, grouped by lineage.

    Parameters
    ----------
    genes : str
        Comma-separated gene symbols/IDs; must resolve to at least 2
        genes.

    Returns
    -------
    MultiLineageRankedResponse
        Per-lineage joint ranking across all requested genes, per
        :func:`ranking.rank_by_lineage_multi`.

    Notes
    -----
    Raises ``HTTPException(422)`` if fewer than 2 genes are given, or
    on an unrecognised gene.
    """
    _guard()
    gene_list = [g.strip() for g in genes.split(",") if g.strip()]
    if len(gene_list) < 2:
        raise HTTPException(status_code=422, detail="Provide at least 2 comma-separated genes.")
    try:
        loop = asyncio.get_event_loop()
        result = await loop.run_in_executor(
            None, lambda: _rank_by_lineage_multi(gene_list)
        )
    except ValueError as exc:
        raise HTTPException(status_code=422, detail=str(exc))
    return MultiLineageRankedResponse(**result)


@router.get("/exclude", response_model=ExcludeResponse)
async def exclude_endpoint(
    gene_a: str = Query(..., description="Gene that must be essential (high score)"),
    gene_b: str = Query(..., description="Gene that must NOT be essential (low score)"),
    lineage: list[str] = Query(default=[], description="Lineage filter (repeatable)"),
    top_n: int = Query(default=30, ge=1, le=200),
) -> ExcludeResponse:
    """
    Rank cell lines where gene_a is essential and gene_b is not (selectivity query).

    Parameters
    ----------
    gene_a : str
        Target gene that must score high (essential).
    gene_b : str
        Gene that must score low (must NOT be essential) — the
        selectivity constraint.
    lineage : list of str, optional
        Repeatable lineage filter.
    top_n : int, default 30
        Maximum number of ranked lines to return (1-200).

    Returns
    -------
    ExcludeResponse
        Lines ranked by gene_a's score among those where gene_b reads
        low, per :func:`ranking.exclude`.

    Notes
    -----
    Raises ``HTTPException(422)`` on an unrecognised gene.
    """
    _guard()
    try:
        loop = asyncio.get_event_loop()
        result = await loop.run_in_executor(
            None, lambda: _exclude(gene_a, gene_b, lineage, top_n)
        )
    except ValueError as exc:
        raise HTTPException(status_code=422, detail=str(exc))
    return ExcludeResponse(**result)


@router.get("/exclude/by-lineage", response_model=ExcludeLineageRankedResponse)
async def exclude_by_lineage_endpoint(
    gene_a: str = Query(..., description="Gene that must be essential (high score)"),
    gene_b: str = Query(..., description="Gene that must NOT be essential (low score)"),
) -> ExcludeLineageRankedResponse:
    """
    Selectivity ranking (gene_a high, gene_b low), grouped by lineage.

    Parameters
    ----------
    gene_a : str
        Target gene that must score high.
    gene_b : str
        Gene that must score low.

    Returns
    -------
    ExcludeLineageRankedResponse
        Per-lineage selectivity ranking, per
        :func:`ranking.exclude_by_lineage`.

    Notes
    -----
    Raises ``HTTPException(422)`` on an unrecognised gene.
    """
    _guard()
    try:
        loop = asyncio.get_event_loop()
        result = await loop.run_in_executor(
            None, lambda: _exclude_by_lineage(gene_a, gene_b)
        )
    except ValueError as exc:
        raise HTTPException(status_code=422, detail=str(exc))
    return ExcludeLineageRankedResponse(**result)


@router.get("/exclude/many", response_model=ExcludeManyResponse)
async def exclude_many_endpoint(
    gene_a: str = Query(..., description="Gene that must be essential (high score)"),
    exclude: list[str] = Query(
        ..., description="Genes that must NOT be essential (repeatable, e.g. ?exclude=KRAS&exclude=TP53)"
    ),
    lineage: list[str] = Query(default=[], description="Lineage filter (repeatable)"),
    top_n: int = Query(default=30, ge=1, le=200),
) -> ExcludeManyResponse:
    """
    Selectivity ranking against several exclusion genes at once.

    Parameters
    ----------
    gene_a : str
        Target gene that must score high.
    exclude : list of str
        Repeatable query param; every gene here must score low.
        1-10 genes after stripping whitespace/empties.
    lineage : list of str, optional
        Repeatable lineage filter.
    top_n : int, default 30
        Maximum number of ranked lines to return (1-200).

    Returns
    -------
    ExcludeManyResponse
        Lines ranked by gene_a's score among those where every gene
        in ``exclude`` reads low, per :func:`ranking.exclude_many`.

    Notes
    -----
    Raises ``HTTPException(422)`` if ``exclude`` is empty, has more
    than 10 genes, or contains an unrecognised gene.
    """
    _guard()
    exclude_genes = [g.strip() for g in exclude if g.strip()]
    if not exclude_genes:
        raise HTTPException(status_code=422, detail="Provide at least 1 exclusion gene.")
    if len(exclude_genes) > 10:
        raise HTTPException(status_code=422, detail="Provide at most 10 exclusion genes.")
    try:
        loop = asyncio.get_event_loop()
        result = await loop.run_in_executor(
            None, lambda: _exclude_many(gene_a, exclude_genes, lineage, top_n)
        )
    except ValueError as exc:
        raise HTTPException(status_code=422, detail=str(exc))
    return ExcludeManyResponse(**result)


@router.get("/exclude/many/by-lineage", response_model=ExcludeManyLineageRankedResponse)
async def exclude_many_by_lineage_endpoint(
    gene_a: str = Query(..., description="Gene that must be essential (high score)"),
    exclude: list[str] = Query(
        ..., description="Genes that must NOT be essential (repeatable, e.g. ?exclude=KRAS&exclude=TP53)"
    ),
) -> ExcludeManyLineageRankedResponse:
    """
    Selectivity ranking against several exclusion genes, grouped by lineage.

    Parameters
    ----------
    gene_a : str
        Target gene that must score high.
    exclude : list of str
        Repeatable query param; 1-10 genes that must score low.

    Returns
    -------
    ExcludeManyLineageRankedResponse
        Per-lineage selectivity ranking, per
        :func:`ranking.exclude_many_by_lineage`.

    Notes
    -----
    Raises ``HTTPException(422)`` if ``exclude`` is empty, has more
    than 10 genes, or contains an unrecognised gene.
    """
    _guard()
    exclude_genes = [g.strip() for g in exclude if g.strip()]
    if not exclude_genes:
        raise HTTPException(status_code=422, detail="Provide at least 1 exclusion gene.")
    if len(exclude_genes) > 10:
        raise HTTPException(status_code=422, detail="Provide at most 10 exclusion genes.")
    try:
        loop = asyncio.get_event_loop()
        result = await loop.run_in_executor(
            None, lambda: _exclude_many_by_lineage(gene_a, exclude_genes)
        )
    except ValueError as exc:
        raise HTTPException(status_code=422, detail=str(exc))
    return ExcludeManyLineageRankedResponse(**result)


@router.get("/genes/exclude/many", response_model=RankExcludeManyResponse)
async def genes_exclude_many_endpoint(
    genes: str = Query(..., description="Comma-separated target genes, e.g. BRAF,TP53"),
    exclude: list[str] = Query(
        ..., description="Genes that must NOT be essential (repeatable, e.g. ?exclude=KRAS&exclude=EGFR)"
    ),
    lineage: list[str] = Query(default=[], description="Lineage filter (repeatable)"),
    floor: float = Query(default=0.50, ge=0.0, le=1.0),
    top_n: int = Query(default=30, ge=1, le=200),
) -> RankExcludeManyResponse:
    """
    Joint ranking across several target genes, with several exclusion genes.

    The most general ranking query: combines the multi-gene joint
    ranking of :func:`genes_endpoint` with the multi-exclusion
    selectivity of :func:`exclude_many_endpoint` in one call.

    Parameters
    ----------
    genes : str
        Comma-separated target genes; must resolve to at least 2.
    exclude : list of str
        Repeatable query param; 1-10 genes that must score low on the
        same lines.
    lineage : list of str, optional
        Repeatable lineage filter.
    floor : float, default 0.50
        Minimum per-target-gene score a line must clear.
    top_n : int, default 30
        Maximum number of ranked lines to return (1-200).

    Returns
    -------
    RankExcludeManyResponse
        Lines ranked by combined target-gene score among those where
        every exclusion gene reads low, per
        :func:`ranking.rank_exclude_many`.

    Notes
    -----
    Raises ``HTTPException(422)`` if fewer than 2 target genes are
    given, ``exclude`` is empty or has more than 10 genes, or any gene
    is unrecognised.
    """
    _guard()
    gene_list = [g.strip() for g in genes.split(",") if g.strip()]
    if len(gene_list) < 2:
        raise HTTPException(status_code=422, detail="Provide at least 2 comma-separated target genes.")
    exclude_genes = [g.strip() for g in exclude if g.strip()]
    if not exclude_genes:
        raise HTTPException(status_code=422, detail="Provide at least 1 exclusion gene.")
    if len(exclude_genes) > 10:
        raise HTTPException(status_code=422, detail="Provide at most 10 exclusion genes.")
    try:
        loop = asyncio.get_event_loop()
        result = await loop.run_in_executor(
            None, lambda: _rank_exclude_many(gene_list, exclude_genes, lineage, floor, top_n)
        )
    except ValueError as exc:
        raise HTTPException(status_code=422, detail=str(exc))
    return RankExcludeManyResponse(**result)


@router.get("/genes/exclude/many/by-lineage", response_model=RankExcludeManyLineageRankedResponse)
async def genes_exclude_many_by_lineage_endpoint(
    genes: str = Query(..., description="Comma-separated target genes, e.g. BRAF,TP53"),
    exclude: list[str] = Query(
        ..., description="Genes that must NOT be essential (repeatable, e.g. ?exclude=KRAS&exclude=EGFR)"
    ),
) -> RankExcludeManyLineageRankedResponse:
    """
    Joint ranking across several target and exclusion genes, grouped by lineage.

    Parameters
    ----------
    genes : str
        Comma-separated target genes; must resolve to at least 2.
    exclude : list of str
        Repeatable query param; 1-10 genes that must score low.

    Returns
    -------
    RankExcludeManyLineageRankedResponse
        Per-lineage version of :func:`genes_exclude_many_endpoint`'s
        result, per :func:`ranking.rank_exclude_many_by_lineage`.

    Notes
    -----
    Raises ``HTTPException(422)`` under the same conditions as
    :func:`genes_exclude_many_endpoint`.
    """
    _guard()
    gene_list = [g.strip() for g in genes.split(",") if g.strip()]
    if len(gene_list) < 2:
        raise HTTPException(status_code=422, detail="Provide at least 2 comma-separated target genes.")
    exclude_genes = [g.strip() for g in exclude if g.strip()]
    if not exclude_genes:
        raise HTTPException(status_code=422, detail="Provide at least 1 exclusion gene.")
    if len(exclude_genes) > 10:
        raise HTTPException(status_code=422, detail="Provide at most 10 exclusion genes.")
    try:
        loop = asyncio.get_event_loop()
        result = await loop.run_in_executor(
            None, lambda: _rank_exclude_many_by_lineage(gene_list, exclude_genes)
        )
    except ValueError as exc:
        raise HTTPException(status_code=422, detail=str(exc))
    return RankExcludeManyLineageRankedResponse(**result)


@router.get("/gene/detail", response_model=DetailResponse)
async def gene_detail_endpoint(
    gene: str = Query(..., description="Gene symbol or ENSG ID"),
    model_id: str = Query(..., description="Cell line ACH ID, e.g. ACH-000088"),
    other_genes: str = Query(
        default="", description="Comma-separated OTHER target genes from the "
        "original multi-gene query (multi/jointSelectivity), if any -- narrows "
        "rna_alternatives to lines with real evidence for every target gene, "
        "not just this one."),
    exclude_genes: str = Query(
        default="", description="Comma-separated exclusion genes from the "
        "original query (selectivity/jointSelectivity), if any -- narrows "
        "rna_alternatives to lines with real evidence for every exclusion "
        "gene too."),
) -> DetailResponse:
    """
    Full evidence detail for one (gene, cell-line) pair.

    Powers the UI's per-line drill-down: every omics layer's
    contribution to this gene/line's score, plus RNA-similar
    alternative lines.

    Parameters
    ----------
    gene : str
        Gene symbol or ENSG ID.
    model_id : str
        Cell line ACH ID, e.g. ``"ACH-000088"``.
    other_genes : str, optional
        Comma-separated OTHER target genes from the original
        multi-gene query (multi/jointSelectivity), if this detail
        view was reached from one. Narrows the returned
        ``rna_alternatives`` to lines with real evidence for every
        target gene, not just this one.
    exclude_genes : str, optional
        Comma-separated exclusion genes from the original query
        (selectivity/jointSelectivity), if any. Narrows
        ``rna_alternatives`` to lines with real evidence for every
        exclusion gene too.

    Returns
    -------
    DetailResponse
        Per-layer evidence and RNA-similar alternatives for this
        exact (gene, model_id) pair, per :func:`ranking.detail`.

    Notes
    -----
    Raises ``HTTPException(422)`` on an unrecognised gene or model_id.
    """
    _guard()
    other_gene_list = [g.strip() for g in other_genes.split(",") if g.strip()]
    exclude_gene_list = [g.strip() for g in exclude_genes.split(",") if g.strip()]
    try:
        loop = asyncio.get_event_loop()
        result = await loop.run_in_executor(
            None, lambda: _detail(gene, model_id, other_gene_list, exclude_gene_list)
        )
    except ValueError as exc:
        raise HTTPException(status_code=422, detail=str(exc))
    return DetailResponse(**result)


@router.get("/genes/list", response_model=GeneListResponse)
async def genes_list_endpoint() -> GeneListResponse:
    """Full recognised gene universe, for client-side autocomplete/validation
    and the Reference lookup page."""
    _guard()
    return GeneListResponse(**list_genes())


@router.get("/cell_lines/list", response_model=CellLineListResponse)
async def cell_lines_list_endpoint() -> CellLineListResponse:
    """Full recognised cell-line universe, for the Reference lookup page."""
    _guard()
    return CellLineListResponse(**list_cell_lines())


@router.get("/health", response_model=HealthResponse)
async def health() -> HealthResponse:
    """
    Liveness check.

    Parameters
    ----------
    None

    Returns
    -------
    HealthResponse
        Always ``status="ok"`` if the process is running to answer at
        all — deliberately does NOT call :func:`_guard`, so this stays
        reachable during cold start while ranking data is still
        loading (use the 503 from any other endpoint to detect that
        state instead).
    """
    return HealthResponse(status="ok")
