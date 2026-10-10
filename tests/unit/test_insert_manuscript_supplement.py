import json

import pytest

from benchmark_tools import insert_manuscript_supplement as current
from benchmark_tools.prepare_ob_candidate_neighborhood import record


def test_exact_parent_restoration_and_stable_newlines():
    parent = "# Title\n\nOld results\n\n### Later\nPreserved body\n"
    revised, piece = current.insert(parent, "### New\nNew table\n", "### Later\n")
    assert piece == "### New\nNew table\n\n"
    assert revised.replace(piece, "", 1) == parent


@pytest.mark.parametrize("parent,supplement,anchor", [
    ("absent", "text", "anchor"), ("anchor anchor", "text", "anchor"),
    ("anchor", "", "anchor"), ("text\n\nanchor", "text\n", "anchor"),
])
def test_ambiguous_empty_or_duplicate_insertion_rejected(parent, supplement, anchor):
    with pytest.raises(ValueError):
        current.insert(parent, supplement, anchor)


def inputs(tmp_path):
    parent = tmp_path / "parent.md"
    parent.write_text("# Parent\n\n### Later\nOld evidence\n")
    supplement = tmp_path / "supplement.md"
    supplement.write_text("### New\nNew table\n")
    return record(parent), record(supplement)


def test_generation_bound_to_inputs_and_parent_byte_restoration(tmp_path):
    parent, supplement = inputs(tmp_path)
    output, receipt = tmp_path / "new.md", tmp_path / "generation.json"
    result = current.run(parent, supplement, output, receipt, "### Later\n")
    assert output.read_bytes().replace(result["inserted_section"].encode(), b"", 1) == (tmp_path / "parent.md").read_bytes()
    assert json.loads(receipt.read_text()) == result
    assert result["publication_ready"] is result["manuscript_rendered"] is False


def test_pins_changed_or_output_occupied_do_not_overwrite(tmp_path):
    parent, supplement = inputs(tmp_path)
    output, receipt = tmp_path / "new.md", tmp_path / "generation.json"
    bad = dict(supplement, sha256="0" * 64)
    with pytest.raises(ValueError):
        current.run(parent, bad, output, receipt, "### Later\n")
    assert not output.exists()
    output.write_text("retain")
    with pytest.raises(ValueError):
        current.run(parent, supplement, output, receipt, "### Later\n")
    assert output.read_text() == "retain" and not receipt.exists()


def test_relocation_that_would_break_relative_links_rejected(tmp_path):
    parent, supplement = inputs(tmp_path)
    other = tmp_path / "other"
    other.mkdir()
    with pytest.raises(ValueError):
        current.run(parent, supplement, other / "new.md", other / "generation.json", "### Later\n")


def test_dangling_output_symlink_is_preserved(tmp_path):
    parent, supplement = inputs(tmp_path)
    output = tmp_path / "new.md"
    target = tmp_path / "missing.md"
    output.symlink_to(target)
    with pytest.raises(ValueError, match="Outputs must be fresh"):
        current.run(parent, supplement, output, tmp_path / "generation.json", "### Later\n")
    assert output.is_symlink() and not target.exists()
