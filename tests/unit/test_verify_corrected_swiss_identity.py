import math

import pytest

from benchmark_tools.prepare_ob_candidate_neighborhood import record
from benchmark_tools.verify_corrected_swiss_identity import fasta, scalar_identity, verify_family


def test_scalar_uses_pair_mean_and_canonical_overlap_only():
    aligned = dict(a="ACDX-", b="ACEX-", c="A---X")
    result = scalar_identity(aligned)
    assert result["pairs"] == 3
    assert result["pairs_without_canonical_overlap"] == 0
    assert math.isclose(result["mean_pairwise_identity"], (2/3 + 1 + 1)/3)


def test_one_missing_pair_makes_family_missing():
    result = scalar_identity(dict(a="A-", b="-C", c="AC"))
    assert result == dict(pairs=3, pairs_without_canonical_overlap=1, mean_pairwise_identity=None)


@pytest.mark.parametrize("aligned", [{}, {"a": "A"}, {"a": "A", "b": "AA"}, {"a": "", "b": ""}])
def test_invalid_dimensions(aligned):
    with pytest.raises(ValueError, match="dimensions"):
        scalar_identity(aligned)


@pytest.mark.parametrize("problem", [None, "residue", "membership", "identity", "pairs", "stops", "dimensions", "status", "hash"])
def test_family_validation(tmp_path, problem):
    input_path, aligned = tmp_path / "input.faa", tmp_path / "aligned.faa"
    input_path.write_text(">a\nACX\n>b\nAEX\n")
    aligned.write_text(">a\nACX-\n>b\nAEX-\n")
    if problem == "residue":
        aligned.write_text(">a\nACX-\n>b\nADX-\n")
    elif problem == "membership":
        aligned.write_text(">a\nACX-\n>c\nAEX-\n")
    run = dict(status="alignment_validated", exit_code=0, accuracy_evaluated=False,
        input=record(input_path), alignment=record(aligned), genes=2, columns=4,
        identity=dict(pairs=1, pairs_without_canonical_overlap=0, mean_pairwise_identity=.5),
        removed_stop_symbols={"a": 1})
    if problem == "identity":
        run["identity"]["mean_pairwise_identity"] = .9
    elif problem == "pairs":
        run["identity"]["pairs"] = 2
    elif problem == "stops":
        run["removed_stop_symbols"] = {}
    elif problem == "dimensions":
        run["columns"] = 5
    elif problem == "status":
        run["exit_code"] = 1
    elif problem == "hash":
        aligned.write_text(aligned.read_text() + "\n")
    if problem:
        with pytest.raises(ValueError):
            verify_family(run, dict(a="ACX*", b="AEX"))
    else:
        assert verify_family(run, dict(a="ACX*", b="AEX")) == run["identity"]


def test_duplicate_fasta_rejected(tmp_path):
    path = tmp_path / "duplicate.faa"
    path.write_text(">a\nA\n>a\nA\n")
    with pytest.raises(ValueError, match="duplicate"):
        fasta(path)
