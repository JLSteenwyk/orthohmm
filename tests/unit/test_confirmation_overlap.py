from pathlib import Path
from zipfile import ZipFile

import pytest

from benchmark_tools import audit_confirmation_overlap as auditor


def archive(tmp_path, members, name="candidate.zip"):
    path = tmp_path / name
    with ZipFile(path, "w") as output:
        for label, text in members:
            output.writestr(label, text)
    return path


def test_identifier_renaming_and_duplicate_sequences_cannot_hide_overlap(tmp_path):
    candidate = archive(tmp_path, [("A.faa", ">renamed\nacde*\n>copy\nACDE\n"),
                                   ("B.faa", ">other\nFGHI\n")])
    old = tmp_path / "development.faa"
    old.write_text(">old_id\nACDE\n>second\nACDE\n>unmatched\nKLMN\n")
    result = auditor.audit([candidate], {"exposed": [old]})
    assert result["candidate_proteins"] == 3
    assert result["candidate_unique_normalized_sequences"] == 2
    overlap = result["overlaps"]["exposed"]
    assert overlap["matching_candidate_proteins"] == 2
    assert overlap["matching_development_records"] == 2
    assert overlap["matching_candidate_fraction"] == pytest.approx(2 / 3)
    assert result["strict_sequence_disjoint_confirmation_ruled_out"]
    assert not result["scoring_performed"]
    assert not result["independence_established"]


def test_no_exact_matches_does_not_establish_independence(tmp_path):
    candidate = archive(tmp_path, [("A.fa", ">x\nAC*DE\n")])
    old = tmp_path / "old.fa"
    old.write_text(">old\nACDE\n")
    result = auditor.audit([candidate], {"exposed": [old]})
    assert result["overlaps"]["exposed"]["matching_candidate_proteins"] == 0
    assert not result["strict_sequence_disjoint_confirmation_ruled_out"]
    assert not result["independence_established"]


def test_truth_and_hidden_metadata_are_not_read_or_extracted(tmp_path):
    candidate = archive(tmp_path, [("A.faa", ">a\nACDE\n"),
        ("reference_classification.txt", "not a FASTA and never parsed"),
        ("__MACOSX/._A.faa", "invalid"), (".hidden.faa", "invalid")])
    old = tmp_path / "old.fa"
    old.write_text(">old\nFGHI\n")
    result = auditor.audit([candidate], {"exposed": [old]})
    assert len(result["candidate_fasta_members"]) == 1
    assert not result["reference_classifications_read"]
    assert not (tmp_path / "A.faa").exists()


@pytest.mark.parametrize("name", ["../A.fa", "/A.fa", "A\\B.fa", "C:/A.fa"])
def test_unsafe_zip_paths_rejected(tmp_path, name):
    candidate = archive(tmp_path, [(name, ">a\nACDE\n")])
    with ZipFile(candidate) as contents, pytest.raises(ValueError, match="Unsafe"):
        auditor.fasta_members(contents)


def test_duplicate_member_paths_and_uncompressed_bounds(tmp_path, monkeypatch):
    with pytest.warns(UserWarning, match="Duplicate"):
        candidate = archive(tmp_path, [("A.fa", ">a\nACDE\n"), ("A.fa", ">b\nFGHI\n")])
    with ZipFile(candidate) as contents, pytest.raises(ValueError, match="duplicate"):
        auditor.fasta_members(contents)
    candidate = archive(tmp_path, [("A.fa", ">a\nACDE\n")], "bounded.zip")
    monkeypatch.setattr(auditor, "MAX_BYTES", 1)
    with ZipFile(candidate) as contents, pytest.raises(ValueError, match="size"):
        auditor.fasta_members(contents)


@pytest.mark.parametrize("text", ["", ">a\n", ">a\n*\n", ">a\nACDE\n>a\nFGHI\n"])
def test_empty_and_duplicate_fasta_records_rejected(tmp_path, text):
    candidate = archive(tmp_path, [("A.fa", text)])
    old = tmp_path / "old.fa"
    old.write_text(">old\nACDE\n")
    with pytest.raises(ValueError, match="Empty|duplicate"):
        auditor.audit([candidate], {"exposed": [old]})


def test_same_member_in_two_archives_is_ambiguous(tmp_path):
    a = archive(tmp_path, [("A.fa", ">a\nACDE\n")])
    b = archive(tmp_path, [("A.fa", ">b\nFGHI\n")], "second.zip")
    old = tmp_path / "old.fa"
    old.write_text(">old\nACDE\n")
    with pytest.raises(ValueError, match="across archives"):
        auditor.audit([a, b], {"exposed": [old]})


def test_input_mutation_invalidates_audit(tmp_path, monkeypatch):
    candidate = archive(tmp_path, [("A.fa", ">a\nACDE\n")])
    old = tmp_path / "old.fa"
    old.write_text(">old\nACDE\n")
    original = auditor.proteins
    calls = 0
    def mutated(handle):
        nonlocal calls
        yield from original(handle)
        calls += 1
        if calls == 2:
            old.write_text(">old\nFGHI\n")
    monkeypatch.setattr(auditor, "proteins", mutated)
    with pytest.raises(ValueError, match="changed"):
        auditor.audit([candidate], {"exposed": [old]})


def test_cross_member_id_collisions_are_disclosed_not_merged(tmp_path):
    candidate = archive(tmp_path, [("A.fa", ">same\nACDE\n"), ("B.fa", ">same\nFGHI\n")])
    old = tmp_path / "old.fa"
    old.write_text(">old\nFGHI\n")
    result = auditor.audit([candidate], {"exposed": [old]})
    assert result["candidate_identifier_collisions_across_members"] == 1
    assert result["overlaps"]["exposed"]["matching_candidate_proteins"] == 1
