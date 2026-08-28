import gzip
import json

import pytest

from benchmark_tools.qfo_filter_pairs import filter_pairs, load_mapping


def test_filter_pairs_keeps_only_mapped_identifiers(tmp_path):
    mapping_path = tmp_path / "mapping.json.gz"
    with gzip.open(mapping_path, "wt") as handle:
        json.dump({"mapping": {"A": 1, "B": 2, "C": 3}}, handle)

    input_path = tmp_path / "pairs.tsv"
    input_path.write_text("A\tB\nA\tmissing\nC\tB\n")
    output_path = tmp_path / "pairs.qfo.tsv"

    total, retained = filter_pairs(input_path, output_path, load_mapping(mapping_path))

    assert (total, retained) == (3, 2)
    assert output_path.read_text() == "A\tB\nC\tB\n"


def test_filter_pairs_rejects_malformed_rows(tmp_path):
    input_path = tmp_path / "pairs.tsv"
    input_path.write_text("A\tB\textra\n")

    with pytest.raises(ValueError, match="Expected two tab-separated IDs"):
        filter_pairs(input_path, tmp_path / "out.tsv", {"A": 1, "B": 2})
