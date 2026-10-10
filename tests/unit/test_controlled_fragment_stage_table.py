import csv
from pathlib import Path

import pytest

from benchmark_tools import readback_controlled_fragment_stage_table as reader


def test_original_failed_reader_is_unchanged():
    assert reader.previous.pin(reader.previous.__file__)["sha256"] == reader.PREVIOUS_SHA


def test_serialization_preserves_unavailable_versus_false():
    row = dict.fromkeys(reader.FIELDS, None)
    row.update(predicted=False, truth=True, seed=20261101)
    observed = reader.serialized(row)
    assert observed["predicted"] == "False" and observed["truth"] == "True"
    assert observed["hit_forward"] == "NA" and observed["seed"] == "20261101"


@pytest.mark.parametrize("replacement", ["NA", "", "False"])
def test_exact_na_cells_not_empty_or_false(tmp_path, replacement):
    row = dict.fromkeys(reader.FIELDS, None)
    path = tmp_path / "stages.tsv"
    serialized = reader.serialized(row)
    serialized["hit_forward"] = replacement
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, reader.FIELDS, delimiter="\t")
        writer.writeheader()
        writer.writerow(serialized)
    if replacement == "NA":
        assert reader.validate_table(path, [row]) == 1
    else:
        with pytest.raises(ValueError, match="TSV value mismatch"):
            reader.validate_table(path, [row])


def test_table_header_and_count_are_strict(tmp_path):
    path = tmp_path / "stages.tsv"
    path.write_text("case_id\nCase0000\n")
    with pytest.raises(ValueError, match="header"):
        reader.validate_table(path, [])


def test_only_independent_kernels_are_reused():
    text = Path(reader.__file__).read_text()
    assert "import trace_controlled_fragment_stages" not in text
    assert "previous.read_context" in text and "previous.observation" in text
