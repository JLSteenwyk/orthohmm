import pytest

from benchmark_tools.audit_orthomcl_input_parity import compare_inputs


def inputs(tmp_path):
    a, b, combined, gg = [tmp_path / name for name in ("a.fasta", "b.fasta", "all.fa", "all.gg")]
    a.write_text(">A description\nABZ\n")
    b.write_text(">B\nACD\n")
    combined.write_text(">B\nAC\nD\n>A\nABZ\n")
    gg.write_text("b: B\na: A\n")
    return [a, b], combined, gg


def test_exact_and_changed(tmp_path):
    sources, combined, gg = inputs(tmp_path)
    result = compare_inputs(sources, combined, gg)
    assert result["identical_sequences"] == 2
    assert result["genome_map_complete_and_correct"]
    combined.write_text(">A\nAXX\n>B\nacd\n")
    result = compare_inputs(sources, combined, gg)
    assert result["identical_sequences"] == 0
    assert result["sequence_difference_count"] == 2


@pytest.mark.parametrize("text", [">A\nABZ\n", ">A\nABZ\n>A\nABZ\n", ">C\nACD\n"])
def test_invalid_combined(tmp_path, text):
    sources, combined, gg = inputs(tmp_path)
    combined.write_text(text)
    with pytest.raises(ValueError):
        compare_inputs(sources, combined, gg)


@pytest.mark.parametrize("text", ["a: A\n", "a: A\na: B\n", "a: B\nb: A\n",
                                  "a: A A\nb: B\n", "a: C\nb: B\n", "a A\nb: B\n"])
def test_invalid_genome_map(tmp_path, text):
    sources, combined, gg = inputs(tmp_path)
    gg.write_text(text)
    with pytest.raises(ValueError):
        compare_inputs(sources, combined, gg)


def test_duplicate_original_id(tmp_path):
    sources, combined, gg = inputs(tmp_path)
    sources[1].write_text(">A\nABZ\n")
    with pytest.raises(ValueError, match="Duplicate original"):
        compare_inputs(sources, combined, gg)


def test_empty_original(tmp_path):
    sources, combined, gg = inputs(tmp_path)
    sources[0].write_text("")
    with pytest.raises(ValueError, match="Empty original"):
        compare_inputs(sources, combined, gg)
