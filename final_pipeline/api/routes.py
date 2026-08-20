"""Routes for the final_pipeline ranking API."""

from __future__ import annotations

import asyncio
from typing import Optional

from fastapi import APIRouter, HTTPException, Query

from .ranking import detail as _detail
from .ranking import exclude as _exclude
from .ranking import is_ready, rank
from .schemas import DetailResponse, ExcludeResponse, HealthResponse, RankRequest, RankResponse

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


@router.get("/health", response_model=HealthResponse)
async def health() -> HealthResponse:
    return HealthResponse(status="ok")
