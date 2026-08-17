"""Tests for the score_explainer tool, against the scaffold's hardcoded mock scores."""

import json

from FunctionCalling import ScoreResult, score_explainer


def test_known_pair_returns_all_fields():
    raw = score_explainer.invoke({"ensg_id": "ENSG00000141510", "model_id": "ACH-000001"})
    data = json.loads(raw)
    print(f"\n[test_known_pair_returns_all_fields] {json.dumps(data, indent=2)}")

    assert data["found"] is True
    assert data["ensg_id"] == "ENSG00000141510"
    assert data["model_id"] == "ACH-000001"
    for field in (
        "raw_expression",
        "raw_proteomics",
        "pit_expression",
        "pit_proteomics",
        "rho_ep",
        "weight_expression",
        "weight_proteomics",
        "core_score",
        "confidence_tier",
    ):
        assert data[field] is not None


def test_unknown_pair_returns_found_false():
    raw = score_explainer.invoke({"ensg_id": "ENSG_UNKNOWN", "model_id": "ACH-UNKNOWN"})
    data = json.loads(raw)
    print(f"\n[test_unknown_pair_returns_found_false] {json.dumps(data, indent=2)}")

    assert data["found"] is False
    assert data["ensg_id"] == "ENSG_UNKNOWN"
    assert data["model_id"] == "ACH-UNKNOWN"
    assert data["message"]


def test_pydantic_validation_field_types():
    raw = score_explainer.invoke({"ensg_id": "ENSG00000141510", "model_id": "ACH-000001"})
    result = ScoreResult.model_validate_json(raw)
    print(f"\n[test_pydantic_validation_field_types] {result.model_dump_json(indent=2)}")

    assert isinstance(result.found, bool)
    assert isinstance(result.core_score, float)
    assert isinstance(result.confidence_tier, str)
    assert isinstance(result.is_tsg, bool)
    assert isinstance(result.driver_gated, bool)
