"""Tests for the gene_alias_lookup tool.

Unit tests below mock the network layer entirely and force a local-index
miss, so they exercise parsing/cross-validation logic offline —
deterministic, always pass, gate CI. The `live` tests hit real
Ensembl/HGNC (or the real local parquet index) and are excluded from CI;
they may fail on third-party downtime (plan.md R7 — this closes gap G4,
where the previous all-live suite silently depended on Ensembl's uptime).
"""

import json

import pytest

import FunctionCalling as fc
from FunctionCalling import _query_ensembl, _query_hgnc, gene_alias_lookup
from gene_index import GeneIndex


@pytest.fixture(autouse=True)
def _reset_circuit_breaker():
    # Shared module-level breaker; clear it so one test's simulated
    # failures can't trip the circuit for the next test.
    fc._ensembl_breaker.record_success()
    yield


@pytest.fixture()
def empty_index(monkeypatch):
    """Forces a local-index miss so tests exercise the network fallback."""
    monkeypatch.setattr(fc, "_GENE_INDEX", GeneIndex())


# ---------------------------------------------------------------------------
# Unit tests — network mocked or local-index-backed, always pass offline
# ---------------------------------------------------------------------------

@pytest.mark.asyncio
async def test_gene_alias_lookup_uses_local_index_first():
    raw = await gene_alias_lookup.ainvoke({"query": "TP53"})
    data = json.loads(raw)

    assert data["found"] is True
    assert data["symbol"] == "TP53"
    assert data["ensembl_id"] == "ENSG00000141510"
    assert data["sources"] == ["local_index"]
    assert data["cross_validated"] is True


@pytest.mark.asyncio
async def test_gene_alias_lookup_by_ensembl_id_uses_local_index():
    raw = await gene_alias_lookup.ainvoke({"query": "ENSG00000141510"})
    data = json.loads(raw)

    assert data["found"] is True
    assert data["symbol"] == "TP53"


@pytest.mark.asyncio
async def test_gene_alias_lookup_cross_validated_offline(empty_index, monkeypatch):
    async def fake_safe_get(url, *, headers=None, client=None, attempts=3):
        if "ensembl.org" in url:
            return {"id": "ENSG00000141510", "display_name": "TP53", "description": "tumor protein p53"}
        if "genenames.org" in url:
            return {"response": {"docs": [{
                "symbol": "TP53", "name": "tumor protein p53",
                "ensembl_gene_id": "ENSG00000141510",
                "prev_symbol": [], "alias_symbol": ["p53"],
            }]}}
        return None

    monkeypatch.setattr(fc, "safe_get", fake_safe_get)

    raw = await gene_alias_lookup.ainvoke({"query": "TP53"})
    data = json.loads(raw)

    assert data["found"] is True
    assert data["symbol"] == "TP53"
    assert data["ensembl_id"] == "ENSG00000141510"
    assert data["cross_validated"] is True
    assert set(data["sources"]) == {"ensembl", "hgnc"}
    assert data["degraded"] == []


@pytest.mark.asyncio
async def test_gene_alias_lookup_degrades_when_ensembl_down(empty_index, monkeypatch):
    async def fake_safe_get(url, *, headers=None, client=None, attempts=3):
        if "genenames.org" in url:
            return {"response": {"docs": [{
                "symbol": "TP53", "name": "tumor protein p53",
                "ensembl_gene_id": "ENSG00000141510",
                "prev_symbol": [], "alias_symbol": [],
            }]}}
        return None  # ensembl unreachable

    monkeypatch.setattr(fc, "safe_get", fake_safe_get)

    raw = await gene_alias_lookup.ainvoke({"query": "TP53"})
    data = json.loads(raw)

    assert data["found"] is True
    assert data["sources"] == ["hgnc"]
    assert data["degraded"] == ["ensembl"]
    assert data["cross_validated"] is False


@pytest.mark.asyncio
async def test_gene_alias_lookup_not_found_offline(empty_index, monkeypatch):
    async def fake_safe_get(url, *, headers=None, client=None, attempts=3):
        return None

    monkeypatch.setattr(fc, "safe_get", fake_safe_get)

    raw = await gene_alias_lookup.ainvoke({"query": "NONEXISTENT_GENE_XYZ"})
    data = json.loads(raw)

    assert data["found"] is False
    assert set(data["degraded"]) == {"ensembl", "hgnc"}


# ---------------------------------------------------------------------------
# Live tests — real Ensembl/HGNC network calls, excluded from CI
# ---------------------------------------------------------------------------

@pytest.mark.live
@pytest.mark.asyncio
async def test_query_ensembl_valid_gene():
    result = await _query_ensembl("TP53")
    assert result is not None
    assert result["id"] == "ENSG00000141510"


@pytest.mark.live
@pytest.mark.asyncio
async def test_query_hgnc_valid_gene():
    result = await _query_hgnc("TP53")
    assert result is not None
    assert result["symbol"] == "TP53"


@pytest.mark.live
@pytest.mark.asyncio
async def test_query_ensembl_unknown_gene():
    result = await _query_ensembl("NONEXISTENT_GENE")
    assert result is None


@pytest.mark.live
@pytest.mark.asyncio
async def test_gene_alias_lookup_returns_her2_synonym():
    raw = await gene_alias_lookup.ainvoke({"query": "ERBB2"})
    data = json.loads(raw)
    assert data["found"] is True
    assert "HER2" in data["synonyms"]
