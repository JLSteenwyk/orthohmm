import pytest

from benchmark_tools.audit_three_kingdoms_orthomcl_inputs import compare, mapping, sequence_inventory


def test_exact_content_ignores_mapping_order():
    assert compare({"a": 1, "b": 2}, {"b": 2, "a": 1})["content_equal"]


def test_missing_extra_and_changed_sequences_retained():
    result = compare({"a": 1, "b": 2}, {"a": 3, "c": 4})
    assert result["missing_ids"] == ["b"]
    assert result["extra_ids"] == ["c"]
    assert result["changed_sequences"] == [dict(gene="a", expected=1, observed=3)]
    assert not result["content_equal"]


def test_wrapped_fasta_identity(tmp_path):
    a, b = tmp_path / "a.fa", tmp_path / "b.fa"
    a.write_text(">gene description\nACD\nEF\n")
    b.write_text(">gene\nACDEF\n")
    assert sequence_inventory(a) == sequence_inventory(b)


@pytest.mark.parametrize("text", ["", ">a\nAC\n>a\nAD\n", ">a\n"])
def test_reject_bad_fasta(tmp_path, text):
    path = tmp_path / "input.fa"
    path.write_text(text)
    with pytest.raises(ValueError):
        sequence_inventory(path)


@pytest.mark.parametrize("text", ["A: a a\nB: b\n", "A: a\n", "A: a b\nB: b\n"])
def test_reject_bad_species_map(tmp_path, text):
    path = tmp_path / "all.gg"
    path.write_text(text)
    with pytest.raises(ValueError):
        mapping(path, {"A": {"a"}, "B": {"b"}})
