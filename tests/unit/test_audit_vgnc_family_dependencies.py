import gzip

import pytest

from benchmark_tools.audit_vgnc_family_dependencies import blocks, summarize_raw


def test_transitive_reference_overlap_not_prediction_based():
    mapping, summary = blocks({(1, 2): "A", (2, 3): "B", (3, 4): "C", (5, 6): "D"})
    assert mapping == {"A": "A", "B": "A", "C": "A", "D": "D"}
    assert summary["shared_proteins"] == 2
    assert summary["reference_blocks"] == 2
    assert summary["merged_label_groups"] == [["A", "B", "C"]]


def raw(tmp_path, text):
    path = tmp_path / "raw.gz"
    with gzip.open(path, "wt") as stream:
        stream.write(text)
    return path


def test_cross_block_fp_retained_and_isolates_counted(tmp_path):
    path = raw(tmp_path, "a\tb\tTP\tA\tB\ts1\ts2\na\tc\tFP\tA\tC\ts1\ts2\n")
    result = summarize_raw(path, {"A": "A", "B": "A", "C": "C", "D": "D"})
    assert result["within_reference_block"]["TP"] == 1
    assert result["between_reference_blocks"]["FP"] == 1
    assert result["prediction_link_components"] == 2
    assert result["largest_prediction_link_component_blocks"] == 2


@pytest.mark.parametrize("text", [
    "a\tb\tTP\tA\tB\ts1\ts2\n", "a\ta\tFP\tA\tB\ts1\ts2\n",
    "a\tb\tFP\tA\tmissing\ts1\ts2\n", "a\tb\tOTHER\tA\tB\ts1\ts2\n",
    "a\tb\tFP\tA\tB\ts1\ts2\na\tb\tFP\tA\tB\ts1\ts2\n",
    "a\tb\tFP\tA\tB\ts1\ts2\na\tc\tFP\tB\tA\ts1\ts2\n",
])
def test_invalid_rows(tmp_path, text):
    with pytest.raises(ValueError):
        summarize_raw(raw(tmp_path, text), {"A": "A", "B": "B"})
