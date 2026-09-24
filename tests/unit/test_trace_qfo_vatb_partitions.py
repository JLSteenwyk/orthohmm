import pytest

from benchmark_tools.trace_qfo_vatb_partitions import trace


def test_selected_membership_and_pair_opportunity(tmp_path):
    path = tmp_path / "partition"
    path.write_text("sp|A|NAME tr|B|NAME unrelated\nC\nD other\n")
    result = trace(path, ["A", "B", "C", "D"])
    assert result["represented_groups"] == 3
    assert result["reference_genes"] == 4
    assert result["singleton_reference_genes"] == 1
    assert result["co_grouped_reference_pairs"] == 1
    assert [g["total_genes"] for g in result["groups"]] == [3, 1, 2]


@pytest.mark.parametrize("text,accessions", [
    ("sp|A|NAME\ntr|A|OTHER\n", ["A"]),
    ("A A\n", ["A"]),
    ("A\n", ["A", "B"]),
    ("A\n", ["A", "A"]),
    ("A\n", []),
    ("unrecognized|A|NAME\n", ["A"]),
])
def test_invalid_or_ambiguous_inventory(tmp_path, text, accessions):
    path = tmp_path / "partition"
    path.write_text(text)
    with pytest.raises(ValueError):
        trace(path, accessions)
