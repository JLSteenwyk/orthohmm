import copy

import pytest

from benchmark_tools.consolidate_qfo_orthofinder_provenance import bind_conversion


@pytest.fixture
def evidence():
    row = dict(key="sequence_only", prediction_semantics="checkpoint cliques", retained_pairs=10,
               submitted_pairs=10, scores={"GO": 0.5}, secondary_mean=0.5)
    conversion = dict(semantics="checkpoint cliques", retained_pairs=10, total_pairs=10,
        removed_mapping_pairs=0, pairs={"sha256": "one"}, filtered_pairs={"sha256": "one"},
        finished_epoch=125.0, started_epoch=100.0, admission={"sha256": "native"}, job_id="123")
    return row, conversion


def test_conversion_time_is_not_inference_time(evidence):
    result = bind_conversion(*evidence)
    assert result["conversion_wall_seconds"] == 25
    assert result["standalone_sequence_only_inference_seconds"] is None
    assert result["scores"] == evidence[0]["scores"]


@pytest.mark.parametrize("key,value", [("semantics", "root-HOG pairs"), ("retained_pairs", 9),
    ("total_pairs", 11), ("removed_mapping_pairs", 1), ("finished_epoch", 99)])
def test_wrong_binding_rejected(evidence, key, value):
    row, conversion = evidence
    conversion[key] = value
    with pytest.raises(ValueError):
        bind_conversion(row, conversion)


def test_reference_filter_cannot_change_pair_bytes(evidence):
    row, conversion = copy.deepcopy(evidence)
    conversion["filtered_pairs"]["sha256"] = "different"
    with pytest.raises(ValueError):
        bind_conversion(row, conversion)
