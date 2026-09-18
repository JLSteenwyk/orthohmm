from pathlib import Path

import pytest

from benchmark_tools import prepare_qfo_sequence_pairs as module


def admitted():
    prediction = dict(path="/groups", sha256="groups", bytes=1)
    return dict(variant="all_hits", status="corrected_sequence_graph_admitted", accuracy_evaluated=False,
        publication_ready=False, numeric=dict(status="numeric_checkpoint_verified", summary=dict(genes=984137, species=78)),
        coverage=[dict(label="multipass", output={}), dict(label="multipass_refined", output=prediction)],
        prediction=prediction, checked_records=[prediction])


@pytest.mark.parametrize("problem", [None, "variant", "status", "accuracy", "genes", "species", "stages", "prediction", "inventory"])
def test_admission_selection(problem):
    value = admitted()
    if problem == "variant":
        value["variant"] = "top100"
    elif problem == "status":
        value["status"] = "failed"
    elif problem == "accuracy":
        value["accuracy_evaluated"] = True
    elif problem in ("genes", "species"):
        value["numeric"]["summary"][problem] = 1
    elif problem == "stages":
        value["coverage"].reverse()
    elif problem == "prediction":
        value["prediction"] = {}
    elif problem == "inventory":
        value["checked_records"] = []
    if problem:
        with pytest.raises(ValueError):
            module.validate_admission(value, "all_hits")
    else:
        assert module.validate_admission(value, "all_hits") == value["prediction"]


def fixture(tmp_path, groups):
    fasta = tmp_path / "input"
    fasta.mkdir()
    (fasta / "s1.fasta").write_text(">sp|A|A_S1\nMAAA\n>sp|B|B_S1\nMAAA\n")
    (fasta / "s2.fasta").write_text(">sp|C|C_S2\nMAAA\n")
    (fasta / "s3.fasta").write_text(">sp|D|D_S3\nMAAA\n")
    partition = tmp_path / "groups.txt"
    partition.write_text(groups)
    output = tmp_path / "output"
    output.mkdir()
    converter = Path(__file__).resolve().parents[2] / "qfo_benchmark/og_to_pairwise.py"
    return partition, fasta, converter, {"A", "B", "C", "D"}, output


def test_real_converter_cross_species(tmp_path):
    args = fixture(tmp_path, "sp|A|A_S1 sp|B|B_S1 sp|C|C_S2\nsp|D|D_S3\n")
    result = module.convert(*args)
    assert result["total_pairs"] == result["expected_pairs"] == result["retained_pairs"] == 2
    assert (args[-1] / "pairs.qfo.partial.tsv").read_text() == "A\tC\nB\tC\n"


def test_valid_zero_pairs_retained(tmp_path):
    args = fixture(tmp_path, "A B\nC\nD\n")
    assert module.convert(*args)["total_pairs"] == 0
    assert (args[-1] / "pairs.qfo.partial.tsv").read_bytes() == b""


@pytest.mark.parametrize("problem", ["mapping_loss", "missing", "duplicate", "foreign", "pair_count"])
def test_conversion_failure(tmp_path, monkeypatch, problem):
    groups = {"missing": "A B C\n", "duplicate": "A B C\nD A\n", "foreign": "A B C\nX\n"}.get(problem, "A B C\nD\n")
    args = list(fixture(tmp_path, groups))
    if problem == "mapping_loss":
        args[3] = {"A", "B", "D"}
    if problem == "pair_count":
        monkeypatch.setattr(module, "expected_pairs", lambda *args: 3)
    with pytest.raises(ValueError):
        module.convert(*args)
    assert not (args[-1] / "pairs.tsv").exists()


def test_unscheduled_conversion_rejected(tmp_path, monkeypatch):
    monkeypatch.delenv("SLURM_JOB_ID", raising=False)
    with pytest.raises(ValueError, match="scheduled"):
        module.prepare(tmp_path, "all_hits", tmp_path / "absent", "hash", "job")
    assert not (tmp_path / "benchmarks").exists()
