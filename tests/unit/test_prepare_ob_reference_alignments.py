import json
from pathlib import Path
from types import SimpleNamespace

import pytest

from benchmark_tools import prepare_ob_reference_alignments as aligner


def test_explicit_stop_removal_retains_unknown_residue():
    assert aligner.normalize("a*X*c*") == "AXC"


@pytest.mark.parametrize("sequence", ["", "***", "AC-D", "AC?D", "ACUD"])
def test_unexpected_reference_alphabet_rejected(sequence):
    with pytest.raises(ValueError):
        aligner.normalize(sequence)


def test_identity_excludes_ambiguous_and_gap_positions():
    result = aligner.mean_canonical_identity({"a": "ACX-", "b": "ADXQ"})
    assert result == {"mean_pairwise_identity": .5, "pairs": 1, "pairs_without_canonical_overlap": 0}


def test_any_incomparable_pair_makes_family_missing():
    result = aligner.mean_canonical_identity({"a": "A-X", "b": "A-X", "c": "-AX"})
    assert result["mean_pairwise_identity"] is None
    assert result["pairs"] == 3
    assert result["pairs_without_canonical_overlap"] == 2


@pytest.mark.parametrize("contents", [">a\nAC\n>a\nAC\n", ">a\nAC\n>b\nA\n", ">a\nAD\n>b\nAC\n"])
def test_alignment_inventory_length_and_residues_checked(tmp_path, contents):
    path = tmp_path / "alignment.fa"
    path.write_text(contents)
    with pytest.raises(ValueError):
        aligner.validate_alignment(path, {"a": "AC", "b": "AC"})


@pytest.mark.parametrize("failed", [False, True])
def test_alignment_command_normalization_and_failure_record(tmp_path, monkeypatch, failed):
    def run(command, **kwargs):
        assert command[1:5] == ["--amino", "--auto", "--thread", "1"]
        assert Path(command[-1]).read_text() == ">a\nAAA\n>b\nAA\n"
        if not failed:
            kwargs["stdout"].write(">a\nAAA\n>b\nA-A\n")
        return SimpleNamespace(returncode=int(failed))
    monkeypatch.setattr(aligner.subprocess, "run", run)
    result = aligner.align_family("RefOG001.txt", {"a", "b"}, {"a": "AAA", "b": "A*A"}, tmp_path, Path("/mafft"), {})
    assert result["removed_stop_symbols"] == {"b": 1}
    assert result["status"] == ("failed" if failed else "alignment_validated")
    assert json.loads((tmp_path / "RefOG001/status.json").read_text())["status"] == result["status"]
    if not failed:
        assert result["identity"]["mean_pairwise_identity"] == 1.


def test_tool_inventory_includes_companion_executables(tmp_path):
    entry = tmp_path / "bin/mafft"
    entry.parent.mkdir()
    entry.write_text("entry")
    prefix = tmp_path / "libexec/mafft"
    prefix.mkdir(parents=True)
    (prefix / "aligner").write_bytes(b"binary")
    observed_prefix, inventory = aligner.tool_inventory(entry)
    assert observed_prefix == prefix
    assert {Path(item["path"]) for item in inventory} == {entry, prefix / "aligner"}


def test_no_overwrite_before_input_access(tmp_path, monkeypatch):
    monkeypatch.setattr(aligner, "read_frozen", lambda *a: pytest.fail("Unexpected input access"))
    with pytest.raises(FileExistsError):
        aligner.prepare(tmp_path, tmp_path)
