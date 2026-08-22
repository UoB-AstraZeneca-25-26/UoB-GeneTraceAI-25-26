import pandas as pd

from src.scripts import transcriptomics as tx
from src.scripts.harmonisation_functions import (
    filter_protein_coding_expression_tables,
    summarize_missing_gene_biotypes,
    summarize_missing_model_species,
)


def test_filter_valid_ensg_ids_keeps_only_roster_ids():
    genes = [
        "ENSG00000123456",
        "ENSG00000123456.7",
        "ENSG00000999999",
        "NOT_A_GENE",
        "ENSP00000345678",
    ]
    roster = {"ENSG00000123456", "ENSG00000999999"}

    filtered = tx.filter_valid_ensg_ids(genes, roster)

    assert filtered == ["ENSG00000123456", "ENSG00000999999"]


def test_filter_protein_coding_expression_tables_keeps_only_roster_ensg_ids():
    roster = pd.DataFrame({"gene_id": ["ENSG00000123456", "ENSG00000999999"]})
    tables = {
        "gene_roster": roster,
        "hpa_rna": pd.DataFrame({
            "gene_id": ["ENSG00000123456", "ENSG00000999999.7", "ENSP00000345678"]
        }),
        "depmap_expr": pd.DataFrame({
            "profile_id": ["A", "B", "C"],
            "ENSG00000123456": [1.0, 2.0, 3.0],
            "ENSG00000999999.7": [4.0, 5.0, 6.0],
            "ENSP00000345678": [7.0, 8.0, 9.0],
        }),
        "geo_expr": pd.DataFrame({
            "gsm_id": ["G1", "G2"],
            "ENSG00000999999": [1.0, 2.0],
            "ENSG00000123456.7": [3.0, 4.0],
            "ENST00000000001": [5.0, 6.0],
        }),
    }

    out = filter_protein_coding_expression_tables(tables, tables["gene_roster"])

    assert list(out["hpa_rna"]["gene_id"]) == ["ENSG00000123456", "ENSG00000999999"]
    assert list(out["depmap_expr"].columns) == ["profile_id", "ENSG00000123456", "ENSG00000999999"]
    assert list(out["geo_expr"].columns) == ["gsm_id", "ENSG00000999999", "ENSG00000123456"]


def test_summarize_missing_gene_biotypes_reports_protein_coding_status():
    tables = {
        "gene_roster": pd.DataFrame({"gene_id": ["ENSG00000123456"]}),
        "mutations": pd.DataFrame({
            "gene_id": ["ENSG00000123456", "ENSG00000234567", "ENSG00000345678"],
            "vepbiotype": ["protein_coding", "lncRNA", "protein_coding"],
        }),
        "depmap_expr": pd.DataFrame({
            "profile_id": ["A", "B"],
            "ENSG00000123456": [1.0, 2.0],
            "ENSG00000234567": [3.0, 4.0],
            "ENSG00000456789": [5.0, 6.0],
        }),
    }

    summary = summarize_missing_gene_biotypes(tables)
    stats = summary["depmap_expr"]

    assert stats["n_missing"] == 2
    assert stats["n_protein_coding"] == 1
    assert stats["n_non_protein_coding"] == 1
    assert stats["protein_coding_examples"] == ["ENSG00000456789"]


def test_summarize_missing_model_species_reports_human_status():
    tables = {
        "cell_line_roster": pd.DataFrame({"model_id": ["ACH-000001"]}),
        "sample_info": pd.DataFrame({
            "model_id": ["ACH-000002", "ACH-000003", "ACH-000004"],
            "cell_line_name": ["A", "B", "C"],
        }),
        "cellosaurus": pd.DataFrame({
            "model_id": ["ACH-000002", "ACH-000004"],
            "species_of_origin": ["Homo sapiens", "Mus musculus"],
        }),
        "model_list_20260709": pd.DataFrame({
            "sample_id": ["ACH-000003"],
            "species": ["Homo sapiens"],
        }),
    }

    summary = summarize_missing_model_species(tables)
    stats = summary["sample_info"]

    assert stats["n_missing"] == 3
    assert stats["n_human"] == 2
    assert stats["n_non_human"] == 1
    assert stats["n_unknown"] == 0
