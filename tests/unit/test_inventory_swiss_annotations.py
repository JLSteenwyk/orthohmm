import pytest

from benchmark_tools.inventory_swiss_annotations import features, summarize


def test_repeats_and_distinct_types_are_separate():
    result = features({"length": 100, "pfam": {
        "pfam_B": {"instance": [[80, 90, 1e-5]]},
        "pfam_A": {"instance": [[10, 20, 1e-5], [40, 50, 1e-5]]}}})
    assert result["pfam_type_count"] == 2
    assert result["pfam_instance_count"] == 3
    assert result["has_repeated_pfam_type"]
    assert [r["domain"] for r in result["ordered_pfam_instances"]] == ["pfam_A", "pfam_A", "pfam_B"]


def test_absent_annotation_is_not_zero_domains():
    zero = features({"length": 50, "pfam": {}})
    result = summarize(["present", "absent"], {"present": zero})
    assert result["annotated_genes"] == result["annotated_zero_pfam_genes"] == 1
    assert result["missing_annotation_genes"] == ["absent"]
    assert summarize(["absent"], {})["median_pfam_types_among_annotated"] is None


@pytest.mark.parametrize("entry", [[-1, 10], [10, 9], [1, 101], [1.5, 10]])
def test_bad_coordinates_rejected(entry):
    with pytest.raises(ValueError, match="coordinates"):
        features({"length": 100, "pfam": {"pfam_A": {"instance": [entry]}}})


def test_duplicate_instance_rejected():
    with pytest.raises(ValueError, match="Duplicate"):
        features({"length": 100, "pfam": {"pfam_A": {"instance": [[1, 10], [1, 10]]}}})


def test_missing_namespace_rejected():
    with pytest.raises(ValueError, match="namespace"):
        features({"length": 100})
