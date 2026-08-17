"""Tool definitions for the CellLineFinder agent.

Each tool is self-contained: no cross-tool dependencies, no shared state.
Every @tool function returns a JSON string (a serialised Pydantic model).
"""

import asyncio
import json
from pathlib import Path

import httpx
from langchain_core.tools import tool
from pydantic import BaseModel

KNOWLEDGE_DIR = Path(__file__).resolve().parent.parent / "knowledge"


# ---------------------------------------------------------------------------
# Tool 1: gene_alias_lookup
# ---------------------------------------------------------------------------

class GeneAliasResult(BaseModel):
    found: bool
    query: str
    symbol: str | None = None
    ensembl_id: str | None = None
    full_name: str | None = None
    previous_symbols: list[str] = []
    synonyms: list[str] = []
    cross_validated: bool = False
    sources: list[str] = []


async def _query_ensembl(query: str) -> dict | None:
    url = f"https://rest.ensembl.org/lookup/symbol/homo_sapiens/{query}"
    headers = {"Content-Type": "application/json"}
    async with httpx.AsyncClient(timeout=10.0) as client:
        response = await client.get(url, headers=headers)
    if response.status_code == 200:
        return response.json()
    return None


async def _query_ensembl_by_id(ensg_id: str) -> dict | None:
    url = f"https://rest.ensembl.org/lookup/id/{ensg_id}"
    headers = {"Content-Type": "application/json"}
    async with httpx.AsyncClient(timeout=10.0) as client:
        response = await client.get(url, headers=headers)
    if response.status_code == 200:
        return response.json()
    return None


async def _query_hgnc(query: str) -> dict | None:
    url = f"https://rest.genenames.org/fetch/symbol/{query}"
    headers = {"Accept": "application/json"}
    async with httpx.AsyncClient(timeout=10.0) as client:
        response = await client.get(url, headers=headers)
    if response.status_code != 200:
        return None
    docs = response.json().get("response", {}).get("docs", [])
    return docs[0] if docs else None


@tool
async def gene_alias_lookup(query: str) -> str:
    """Look up gene aliases and synonyms. Takes a gene name, symbol,
    or Ensembl ID. Returns validated aliases from Ensembl and HGNC."""
    query = query.strip()

    ensembl_coro = (
        _query_ensembl_by_id(query)
        if query.upper().startswith("ENSG")
        else _query_ensembl(query)
    )
    ensembl_data, hgnc_data = await asyncio.gather(ensembl_coro, _query_hgnc(query))

    if ensembl_data is None and hgnc_data is None:
        return GeneAliasResult(found=False, query=query).model_dump_json()

    sources = []
    ensembl_id = None
    symbol = None
    full_name = None
    previous_symbols: list[str] = []
    synonyms: list[str] = []

    if ensembl_data:
        sources.append("ensembl")
        ensembl_id = ensembl_data.get("id")
        symbol = ensembl_data.get("display_name")
        full_name = ensembl_data.get("description")

    hgnc_ensembl_id = None
    if hgnc_data:
        sources.append("hgnc")
        symbol = hgnc_data.get("symbol", symbol)
        full_name = hgnc_data.get("name", full_name)
        hgnc_ensembl_id = hgnc_data.get("ensembl_gene_id")
        previous_symbols = hgnc_data.get("prev_symbol", []) or []
        synonyms = hgnc_data.get("alias_symbol", []) or []
        ensembl_id = ensembl_id or hgnc_ensembl_id

    cross_validated = bool(
        ensembl_data and hgnc_data and hgnc_ensembl_id
        and ensembl_data.get("id") == hgnc_ensembl_id
    )

    result = GeneAliasResult(
        found=True,
        query=query,
        symbol=symbol,
        ensembl_id=ensembl_id,
        full_name=full_name,
        previous_symbols=previous_symbols,
        synonyms=synonyms,
        cross_validated=cross_validated,
        sources=sources,
    )
    return result.model_dump_json()


# ---------------------------------------------------------------------------
# Tool 2: score_explainer
# ---------------------------------------------------------------------------

class ScoreResult(BaseModel):
    found: bool
    ensg_id: str
    model_id: str
    raw_expression: float | None = None
    raw_proteomics: float | None = None
    pit_expression: float | None = None
    pit_proteomics: float | None = None
    rho_ep: float | None = None
    weight_expression: float | None = None
    weight_proteomics: float | None = None
    core_score: float | None = None
    confidence_tier: str | None = None
    is_tsg: bool = False
    driver_gated: bool = False
    message: str = ""


# Hardcoded mock data for initial scaffold.
# TODO: replace with a parquet read from the pipeline scoring output.
_MOCK_SCORES = {
    ("ENSG00000141510", "ACH-000001"): {
        "raw_expression": 4.32,
        "raw_proteomics": 0.85,
        "pit_expression": 0.72,
        "pit_proteomics": 0.68,
        "rho_ep": 0.41,
        "weight_expression": 0.55,
        "weight_proteomics": 0.45,
        "core_score": 0.702,
        "confidence_tier": "high",
        "is_tsg": True,
        "driver_gated": True,
    },
}


def _load_scores(ensg_id: str, model_id: str) -> ScoreResult:
    row = _MOCK_SCORES.get((ensg_id, model_id))
    if row is None:
        return ScoreResult(
            found=False,
            ensg_id=ensg_id,
            model_id=model_id,
            message=f"No scoring data found for {ensg_id} in {model_id}.",
        )
    return ScoreResult(
        found=True,
        ensg_id=ensg_id,
        model_id=model_id,
        message="Scoring intermediates retrieved.",
        **row,
    )


@tool
def score_explainer(ensg_id: str, model_id: str) -> str:
    """Retrieve scoring intermediates for a gene-cell line pair.
    Returns raw values, PIT ranks, weights, core score, and tier."""
    return _load_scores(ensg_id, model_id).model_dump_json()


# ---------------------------------------------------------------------------
# Tool 3: dataset_info
# ---------------------------------------------------------------------------

_DATASETS_PATH = KNOWLEDGE_DIR / "datasets.json"
_KB: dict = json.loads(_DATASETS_PATH.read_text())


@tool
def dataset_info(dataset_name: str) -> str:
    """Get a description of a CellLineFinder data source.
    Valid: depmap, hpa, geo, cellosaurus, ensembl, cosmic, hgnc"""
    key = dataset_name.lower().strip()

    if key in _KB:
        return json.dumps(_KB[key])

    for name, entry in _KB.items():
        if key in name:
            return json.dumps(entry)

    return json.dumps({"error": "not found", "available": sorted(_KB.keys())})
