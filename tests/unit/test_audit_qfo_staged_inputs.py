import pytest

from benchmark_tools.audit_qfo_staged_inputs import compare_stage
from benchmark_tools.prepare_ob_candidate_neighborhood import record


def setup(tmp_path):
    original = tmp_path / "species.fasta"
    original.write_text(">A description\nACD\n")
    stage = tmp_path / "stage"
    stage.mkdir()
    staged = stage / "species.fa"
    staged.write_bytes(original.read_bytes())
    return original, stage, staged


def test_identical_and_changed(tmp_path):
    original, stage, staged = setup(tmp_path)
    result = compare_stage([record(original)], stage, ".fa")
    assert result["all_files_identical"] and result["identical_files"] == 1
    assert not result["files"][0]["same_resolved_path"]
    staged.write_text(">A\nACD\n")
    assert not compare_stage([record(original)], stage, ".fa")["all_files_identical"]


def test_symlink_is_explicit(tmp_path):
    original, stage, staged = setup(tmp_path)
    staged.unlink()
    staged.symlink_to(original)
    row = compare_stage([record(original)], stage, ".fa")["files"][0]
    assert row["is_symlink"] and row["same_resolved_path"]


@pytest.mark.parametrize("extra", ["extra.fa", "extra.fasta", "extra.faa"])
def test_unexpected_input(tmp_path, extra):
    original, stage, staged = setup(tmp_path)
    (stage / extra).write_text(">B\nAAA\n")
    with pytest.raises(ValueError, match="inventory"):
        compare_stage([record(original)], stage, ".fa")


def test_missing_input(tmp_path):
    original, stage, staged = setup(tmp_path)
    staged.unlink()
    with pytest.raises(ValueError, match="inventory"):
        compare_stage([record(original)], stage, ".fa")


def test_duplicate_inventory(tmp_path):
    original, stage, staged = setup(tmp_path)
    with pytest.raises(ValueError, match="Duplicate"):
        compare_stage([record(original), record(original)], stage, ".fa")


def test_original_hash_guard(tmp_path):
    original, stage, staged = setup(tmp_path)
    identity = record(original)
    original.write_text(">A\nDDD\n")
    with pytest.raises(ValueError, match="identity changed"):
        compare_stage([identity], stage, ".fa")
