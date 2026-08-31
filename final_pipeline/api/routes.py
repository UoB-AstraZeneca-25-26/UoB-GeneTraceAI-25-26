"""Routes for the final_pipeline ranking API."""

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
    if not is_ready():
        raise HTTPException(status_code=503, detail="Ranking data not loaded.")


@router.post("/v1/rank", response_model=RankResponse)
async def rank_endpoint(req: RankRequest) -> RankResponse:
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
) -> DetailResponse:
    _guard()
    try:
        loop = asyncio.get_event_loop()
        result = await loop.run_in_executor(
            None, lambda: _detail(gene, model_id)
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
    return HealthResponse(status="ok")
