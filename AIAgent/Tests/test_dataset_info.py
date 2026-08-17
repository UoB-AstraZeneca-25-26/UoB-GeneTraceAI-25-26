"""Tests for the dataset_info tool against knowledge/datasets.json."""

import json

import pytest

from FunctionCalling import _KB, dataset_info

ALL_DATASETS = ["depmap", "hpa", "geo", "cellosaurus", "ensembl", "cosmic", "hgnc"]


def test_dataset_info_depmap():
    raw = dataset_info.invoke({"dataset_name": "depmap"})
    data = json.loads(raw)
    print(f"\n[test_dataset_info_depmap] {json.dumps(data, indent=2)}")

    assert data["name"] == "depmap"
    assert "full_name" in data
    assert "url" in data


def test_dataset_info_case_insensitive():
    raw = dataset_info.invoke({"dataset_name": "DepMap"})
    data = json.loads(raw)
    print(f"\n[test_dataset_info_case_insensitive] {json.dumps(data, indent=2)}")

    assert data["name"] == "depmap"


def test_dataset_info_unknown_returns_available_keys():
    raw = dataset_info.invoke({"dataset_name": "unknown"})
    data = json.loads(raw)
    print(f"\n[test_dataset_info_unknown_returns_available_keys] {json.dumps(data, indent=2)}")

    assert data["error"] == "not found"
    assert set(data["available"]) == set(ALL_DATASETS)


@pytest.mark.parametrize("name", ALL_DATASETS)
def test_all_datasets_return_valid_json(name):
    raw = dataset_info.invoke({"dataset_name": name})
    data = json.loads(raw)
    print(f"\n[test_all_datasets_return_valid_json:{name}] {json.dumps(data, indent=2)}")

    assert data["name"] == name
    for field in (
        "full_name",
        "url",
        "what_it_measures",
        "data_types",
        "cell_line_count",
        "version_used",
        "caveats",
        "reference",
    ):
        assert field in data


def test_knowledge_base_has_seven_entries():
    print(f"\n[test_knowledge_base_has_seven_entries] keys={sorted(_KB.keys())}")

    assert set(_KB.keys()) == set(ALL_DATASETS)
