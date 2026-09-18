import gzip

import pytest

from benchmark_tools.audit_swiss_missing_input_relations import incident_counts
from benchmark_tools.audit_qfo_swiss_counts import HEADER, read_raw


def fixture(tmp_path, rows):
    path = tmp_path / "raw.gz"
    with gzip.open(path, "wt") as stream:
        stream.write(HEADER + "\n" + "\n".join(rows) + "\n")
    return path


def test_incident_counts_once_and_all_labels(tmp_path):
    path = fixture(tmp_path, ["F\tA\tB\tFN", "F\tA\tC\tTN", "F\tB\tC\tFP", "G\tD\tE\tTP"])
    _, truth, _ = read_raw(path, ["F", "G"])
    result = incident_counts(path, ["F", "G"], {"A", "B", "D"}, truth)
    assert result["affected_totals"] == dict(TP=1, FP=1, FN=1, TN=1)
    assert result["both_endpoints_missing"] == 1
    assert result["affected_relations"] == 4
    assert result["missing_genes_incident"] == ["A", "B", "D"]


def test_no_incident_pairs(tmp_path):
    path = fixture(tmp_path, ["F\tA\tB\tTP"])
    _, truth, _ = read_raw(path, ["F"])
    result = incident_counts(path, ["F"], {"absent"}, truth)
    assert result["affected_relations"] == 0
    assert result["affected_totals"] == dict(TP=0, FP=0, FN=0, TN=0)


@pytest.mark.parametrize("rows", [["F\tA\tB\tTN"], ["F\tA\tC\tTP"]])
def test_truth_or_identity_change_rejected(tmp_path, rows):
    path = fixture(tmp_path, rows)
    with pytest.raises(ValueError, match="truth differ"):
        incident_counts(path, ["F"], {"A"}, {("F", "A", "B"): True})


def test_duplicate_relation_rejected(tmp_path):
    path = fixture(tmp_path, ["F\tA\tB\tTP", "F\tB\tA\tTP"])
    with pytest.raises(ValueError, match="Duplicate"):
        incident_counts(path, ["F"], {"A"}, {})
