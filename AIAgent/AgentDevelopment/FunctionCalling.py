"""Tool definitions for the CellLineFinder agent.

Each tool is self-contained: no cross-tool dependencies, no shared state.
Every @tool function returns a JSON string (a serialised Pydantic model).
"""

import asyncio
import json
import os
from pathlib import Path

from langchain_core.tools import tool
from pydantic import BaseModel

from gene_index import GeneIndex, load_gene_index
from http_utils import CircuitBreaker, safe_get

# NOTE: capital K — this matches the actual on-disk directory name
# (AIAgent/Knowledge/). The lowercase "knowledge/" that used to be here
# raised FileNotFoundError at import on case-sensitive filesystems (gap G1).
KNOWLEDGE_DIR = Path(__file__).resolve().parent.parent / "Knowledge"

# reference/gene_lookup.parquet lives at the repo root, one level above
# AIAgent/. Sibling layout is preserved in the Docker image (see Dockerfile).
_DEFAULT_GENE_LOOKUP_PATH = (
    Path(__file__).resolve().parent.parent.parent / "reference" / "gene_lookup.parquet"
)


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
    degraded: list[str] = []


# Live Ensembl/HGNC are the primary source; this local index is now only the
# last-resort fallback when both are unreachable. Kept prewarmed anyway
# (19,213 genes, ~0 ms lookups) so the fallback path never pays disk I/O.
# `set_gene_index` lets the FastAPI lifespan prewarm it once at startup;
# standalone/test use falls back to a lazy load here.
_GENE_INDEX: GeneIndex | None = None


def set_gene_index(index: GeneIndex) -> None:
    global _GENE_INDEX
    _GENE_INDEX = index


def _get_gene_index() -> GeneIndex:
    global _GENE_INDEX
    if _GENE_INDEX is None:
        _GENE_INDEX = load_gene_index(_DEFAULT_GENE_LOOKUP_PATH)
    return _GENE_INDEX


# Ensembl/HGNC are now the primary source for gene_alias_lookup, so a couple
# of transient failures should not permanently route every query to the
# degraded local fallback — threshold raised from 5 to 8 versus the old
# enrichment-only tuning. Still trips and cools down after sustained outages.
_ensembl_breaker = CircuitBreaker(failure_threshold=8, reset_after=60.0)
_hgnc_breaker = CircuitBreaker(failure_threshold=8, reset_after=60.0)


async def _query_ensembl(query: str) -> dict | None:
    if _ensembl_breaker.is_open():
        return None
    url = f"https://rest.ensembl.org/lookup/symbol/homo_sapiens/{query}"
    # GET with no body -> Accept, not Content-Type (fixes R4).
    result = await safe_get(url, headers={"Accept": "application/json"})
    (_ensembl_breaker.record_success if result is not None else _ensembl_breaker.record_failure)()
    return result


async def _query_ensembl_by_id(ensg_id: str) -> dict | None:
    if _ensembl_breaker.is_open():
        return None
    url = f"https://rest.ensembl.org/lookup/id/{ensg_id}"
    result = await safe_get(url, headers={"Accept": "application/json"})
    (_ensembl_breaker.record_success if result is not None else _ensembl_breaker.record_failure)()
    return result


async def _query_hgnc(query: str) -> dict | None:
    if _hgnc_breaker.is_open():
        return None
    url = f"https://rest.genenames.org/fetch/symbol/{query}"
    data = await safe_get(url, headers={"Accept": "application/json"})
    (_hgnc_breaker.record_success if data is not None else _hgnc_breaker.record_failure)()
    if data is None:
        return None
    docs = data.get("response", {}).get("docs", [])
    return docs[0] if docs else None


async def _lookup_from_network(query: str) -> GeneAliasResult:
    ensembl_coro = (
        _query_ensembl_by_id(query)
        if query.upper().startswith("ENSG")
        else _query_ensembl(query)
    )
    ensembl_data, hgnc_data = await asyncio.gather(ensembl_coro, _query_hgnc(query))

    degraded = []
    if ensembl_data is None:
        degraded.append("ensembl")
    if hgnc_data is None:
        degraded.append("hgnc")

    if ensembl_data is None and hgnc_data is None:
        return GeneAliasResult(found=False, query=query, degraded=degraded)

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

    return GeneAliasResult(
        found=True,
        query=query,
        symbol=symbol,
        ensembl_id=ensembl_id,
        full_name=full_name,
        previous_symbols=previous_symbols,
        synonyms=synonyms,
        cross_validated=cross_validated,
        sources=sources,
        degraded=degraded,
    )


@tool
async def gene_alias_lookup(query: str) -> str:
    """Look up gene aliases and synonyms. Takes a gene name, symbol,
    or Ensembl ID. Returns validated aliases from Ensembl and HGNC."""
    query = query.strip()

    # Resolution order: live Ensembl+HGNC (queried concurrently) -> local
    # parquet index only if BOTH network sources are unreachable ->
    # found=False. The local index still covers the full 19,213-gene panel
    # offline, but is now a last-resort fallback rather than the default.
    result = await _lookup_from_network(query)
    if result.found:
        return result.model_dump_json()

    local = _get_gene_index().lookup(query)
    if local is not None:
        fallback = GeneAliasResult(
            found=True,
            query=query,
            symbol=local.symbol,
            ensembl_id=local.ensg_id,
            full_name=local.full_name,
            previous_symbols=local.previous_symbols,
            synonyms=local.synonyms,
            cross_validated=True,
            sources=["local_index"],
            degraded=result.degraded,
        )
        return fallback.model_dump_json()

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
    is_tsg: bool | None = None
    driver_gated: bool = False
    message: str = ""


# The scoring pipeline's intermediates are served by a separate scoring API
# (SCORING_API_URL). Until that endpoint is wired up, the tool runs in
# methodology-only mode: it returns the pair's identity with every numeric
# field null, and the LLM explains the 7-step method from the mathematical
# reference in the system prompt rather than narrating a specific result.
# No scores are ever synthesised here.

_API_NUMERIC_FIELDS = (
    "raw_expression",
    "raw_proteomics",
    "pit_expression",
    "pit_proteomics",
    "rho_ep",
    "weight_expression",
    "weight_proteomics",
    "core_score",
)

_METHODOLOGY_MESSAGE = (
    "Methodology explanation — scoring API not connected. "
    "Explaining the 7-step pipeline methodology for this gene."
)
_METHODOLOGY_MESSAGE_DEGRADED = (
    "Methodology explanation — the scoring API returned no data for this "
    "pair. Explaining the 7-step pipeline methodology for this gene."
)


def _scoring_api_url() -> str:
    """Resolve the scoring endpoint, empty string meaning methodology mode.

    ``os.environ`` wins so a container env var (or a test) can flip modes
    without touching .env; the validated ``Settings`` object is the fallback
    so a value set only in AIAgent/.env is still picked up. Settings is
    imported lazily because this module also runs standalone, outside the
    FastAPI process that owns the `api` package.
    """
    url = os.getenv("SCORING_API_URL", "").strip()
    if url:
        return url
    try:
        from api.settings import get_settings

        return get_settings().scoring_api_url.strip()
    except Exception:  # noqa: BLE001 - no settings available => methodology mode
        return ""


def _as_float(value: object) -> float | None:
    try:
        return float(value)  # type: ignore[arg-type]
    except (TypeError, ValueError):
        return None


def _result_from_api(
    ensg_id: str, model_id: str, payload: dict, *, is_tsg: bool | None
) -> ScoreResult:
    """Map a scoring-API payload onto ScoreResult.

    The API is expected to return the ScoreResult field names; anything it
    omits stays null rather than being filled in with a guess.
    """
    numerics = {name: _as_float(payload.get(name)) for name in _API_NUMERIC_FIELDS}
    tier = payload.get("confidence_tier")
    payload_is_tsg = payload.get("is_tsg", is_tsg)

    return ScoreResult(
        found=bool(payload.get("found", True)),
        ensg_id=ensg_id,
        model_id=model_id,
        confidence_tier=str(tier) if tier is not None else None,
        is_tsg=bool(payload_is_tsg) if payload_is_tsg is not None else None,
        driver_gated=bool(payload.get("driver_gated", False)),
        message=str(payload.get("message") or "Scoring intermediates retrieved."),
        **numerics,
    )


def _methodology_result(
    ensg_id: str, model_id: str, *, is_tsg: bool | None, message: str
) -> ScoreResult:
    """A ScoreResult carrying only the pair's identity and gene class.

    Every numeric field is null on purpose: the LLM has the formulas in its
    system prompt and explains the method, so there is nothing to invent.
    """
    return ScoreResult(
        found=True,
        ensg_id=ensg_id,
        model_id=model_id,
        raw_expression=None,
        raw_proteomics=None,
        pit_expression=None,
        pit_proteomics=None,
        rho_ep=None,
        weight_expression=None,
        weight_proteomics=None,
        core_score=None,
        confidence_tier=None,
        is_tsg=is_tsg,
        driver_gated=False,
        message=message,
    )


@tool
async def score_explainer(ensg_id: str, model_id: str) -> str:
    """Retrieve scoring intermediates for a gene-cell line pair and explain
    the 7-step scoring methodology. If the scoring API is connected and has
    data for this pair, returns real pipeline values. Otherwise explains the
    mathematical methodology generically using math_reference.md, without
    inventing numbers."""
    ensg_id = ensg_id.strip().upper()
    model_id = model_id.strip().upper()

    scoring_url = _scoring_api_url()
    if scoring_url:
        url = f"{scoring_url}?gene={ensg_id}&model={model_id}"
        data = await safe_get(url, headers={"Accept": "application/json"})
        if data:
            return _result_from_api(ensg_id, model_id, data, is_tsg=None).model_dump_json()
        message = _METHODOLOGY_MESSAGE_DEGRADED
    else:
        message = _METHODOLOGY_MESSAGE

    return _methodology_result(
        ensg_id, model_id, is_tsg=None, message=message
    ).model_dump_json()


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
