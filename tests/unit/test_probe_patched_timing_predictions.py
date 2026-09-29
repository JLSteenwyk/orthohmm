import pytest

from benchmark_tools.build_private_timing_environment import record
from benchmark_tools.probe_patched_timing_predictions import parity


def test_prediction_bytes_checked(tmp_path):
    old = tmp_path / "old"
    new = tmp_path / "new"
    old.mkdir()
    new.mkdir()
    (old / "orthohmm_orthogroups.txt").write_text("OG0: a b\n")
    (new / "orthohmm_orthogroups.txt").write_text("OG0: a b\n")
    prior = {"current": {"evidence": [record(old / "orthohmm_orthogroups.txt")]}}
    assert parity(prior, new)[0]["identical"]
    (new / "orthohmm_orthogroups.txt").write_text("OG0: a\n")
    with pytest.raises(ValueError, match="differ"):
        parity(prior, new)


def test_empty_parity_rejected(tmp_path):
    with pytest.raises(ValueError, match="no predictions"):
        parity({"current": {"evidence": []}}, tmp_path)


def test_changed_prior_evidence_rejected(tmp_path):
    path = tmp_path / "S1.fa"
    path.write_text(">a\nAAA\n")
    prior = {"current": {"evidence": [record(path)]}}
    path.write_text(">a\nBBB\n")
    with pytest.raises(ValueError, match="Evidence changed"):
        parity(prior, tmp_path)


def test_phylogeny_output_subdirectory(tmp_path):
    old = tmp_path / "old/orthohmm_phylogeny"
    new = tmp_path / "new/orthohmm_phylogeny"
    old.mkdir(parents=True)
    new.mkdir(parents=True)
    name = "orthohmm_pairwise_orthologs.tsv"
    (old / name).write_text("a\tb\n")
    (new / name).write_text("a\tb\n")
    prior = {"current": {"evidence": [record(old / name)]}}
    assert parity(prior, new.parent)[0]["identical"]
