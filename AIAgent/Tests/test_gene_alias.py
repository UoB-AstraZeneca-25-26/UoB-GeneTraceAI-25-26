"""Tests for the gene_alias_lookup tool. Hits live Ensembl/HGNC REST APIs."""

import json
import pytest
from FunctionCalling import _query_ensembl, _query_hgnc, gene_alias_lookup


@pytest.mark.asyncio
async def test_query_ensembl_valid_gene():
    result = await _query_ensembl("TP53")
    print(f"\n[test_query_ensembl_valid_gene] {json.dumps(result, indent=2)}")

    assert result is not None
    assert result["id"] == "ENSG00000141510"


@pytest.mark.asyncio
async def test_query_hgnc_valid_gene():
    result = await _query_hgnc("TP53")
    print(f"\n[test_query_hgnc_valid_gene] {json.dumps(result, indent=2)}")

    assert result is not None
    assert result["symbol"] == "TP53"


@pytest.mark.asyncio
async def test_query_ensembl_unknown_gene():
    result = await _query_ensembl("NONEXISTENT_GENE")
    print(f"\n[test_query_ensembl_unknown_gene] {result}")

    assert result is None


@pytest.mark.asyncio
async def test_gene_alias_lookup_cross_validated():
    raw = await gene_alias_lookup.ainvoke({"query": "TP53"})
    data = json.loads(raw)
    print(f"\n[test_gene_alias_lookup_cross_validated] {json.dumps(data, indent=2)}")

    assert data["found"] is True
    assert data["symbol"] == "TP53"
    assert data["ensembl_id"] == "ENSG00000141510"
    assert data["cross_validated"] is True
    assert set(data["sources"]) == {"ensembl", "hgnc"}


@pytest.mark.asyncio
async def test_gene_alias_lookup_by_ensembl_id():
    raw = await gene_alias_lookup.ainvoke({"query": "ENSG00000141510"})
    data = json.loads(raw)
    print(f"\n[test_gene_alias_lookup_by_ensembl_id] {json.dumps(data, indent=2)}")

    assert data["found"] is True
    assert data["symbol"] == "TP53"
    assert data["ensembl_id"] == "ENSG00000141510"


@pytest.mark.asyncio
async def test_gene_alias_lookup_returns_her2_synonym():
    raw = await gene_alias_lookup.ainvoke({"query": "ERBB2"})
    data = json.loads(raw)
    print(f"\n[test_gene_alias_lookup_returns_her2_synonym] {json.dumps(data, indent=2)}")

    assert data["found"] is True
    assert "HER2" in data["synonyms"]


@pytest.mark.asyncio
async def test_gene_alias_lookup_not_found():
    raw = await gene_alias_lookup.ainvoke({"query": "NONEXISTENT_GENE_XYZ"})
    data = json.loads(raw)
    print(f"\n[test_gene_alias_lookup_not_found] {json.dumps(data, indent=2)}")

    assert data["found"] is False
