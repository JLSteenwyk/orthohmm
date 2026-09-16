from pathlib import Path

import pytest

from benchmark_tools.run_qfo_publication_replay import command_for, compare_partition, write_json


def test_frozen_label_blind_command():
    argv = command_for(Path("/data"), Path("/launcher"), Path("/output"))
    assert "--official-benchmark" not in argv
    assert argv[argv.index("--cpu") + 1] == "32"
    assert argv[argv.index("--profile-iterations") + 1] == "1"
    assert argv[argv.index("--leiden-seed") + 1] == "4"
    assert "--accuracy-checkpoint" in argv
    assert "--checkpoint-sha256" in argv
    assert not any("jackknife" in arg for arg in argv)


@pytest.mark.parametrize("observed,equal", [("c\nb a\n", True), ("a\nb c\n", False)])
def test_partition_comparison_ignores_order_but_not_membership(tmp_path, observed, equal):
    a, b = tmp_path / "a", tmp_path / "b"
    a.write_text("a b\nc\n")
    b.write_text(observed)
    result = compare_partition(a, b, {"a", "b", "c"})
    assert result["partition_equal"] is equal
    assert result["byte_equal"] is False


@pytest.mark.parametrize("observed", ["a b\n", "a b\nb c\n", "a b\nc d\n"])
def test_incomplete_duplicate_or_foreign_genes_rejected(tmp_path, observed):
    a, b = tmp_path / "a", tmp_path / "b"
    a.write_text("a b\nc\n")
    b.write_text(observed)
    with pytest.raises(ValueError):
        compare_partition(a, b, {"a", "b", "c"})


def test_provenance_cannot_be_overwritten(tmp_path):
    path = tmp_path / "report.json"
    write_json(path, {"status": "original"})
    with pytest.raises(FileExistsError):
        write_json(path, {"status": "replacement"})
