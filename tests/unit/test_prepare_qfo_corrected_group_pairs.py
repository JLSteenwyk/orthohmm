import gzip
import json

import pytest

import benchmark_tools.prepare_qfo_corrected_group_pairs as module
from benchmark_tools.prepare_ob_candidate_neighborhood import record, check


def test_combinatorial_pair_count(tmp_path):
    path = tmp_path / "groups.txt"
    path.write_text("a b c\nd e\nf\n")
    owners = {"a": "s1", "b": "s1", "c": "s2", "d": "s2", "e": "s3", "f": "s1"}
    assert module.expected_pairs(path, owners) == 3


@pytest.mark.parametrize("body", ["a a b\n", "a b\na\n", "a unknown\n", "a\n"])
def test_bad_partition(tmp_path, body):
    path = tmp_path / "groups.txt"
    path.write_text(body)
    with pytest.raises(ValueError):
        module.expected_pairs(path, {"a": "s1", "b": "s2"})


@pytest.mark.parametrize("index", [1, 3, 5, 7, -1, 8, True, "0"])
def test_native_pair_cells_cannot_use_group_conversion(index):
    with pytest.raises(ValueError):
        module.select({}, index)


def setup(tmp_path, monkeypatch, missing_mapping=False):
    fasta = tmp_path / "inputs"
    fasta.mkdir()
    (fasta / "one.fasta").write_text(">sp|A|one\nAAA\n>sp|B|one\nAAA\n")
    (fasta / "two.fasta").write_text(">sp|C|two\nAAA\n")
    partition = tmp_path / "partition.txt"
    partition.write_text("sp|A|one sp|B|one sp|C|two\n")
    mapping = tmp_path / "mapping.json.gz"
    mapped_ids = ("A", "B") if missing_mapping else ("A", "B", "C")
    with gzip.open(mapping, "wt") as stream:
        json.dump({"mapping": dict.fromkeys(mapped_ids, 1)}, stream)
    cells = [{"label": label, "reconciliation": bool(i % 2), "candidate_partition": str(partition)}
             for i, label in enumerate(module.CELLS)]
    manifest = {"cells": cells, "input_fastas": [record(p) for p in sorted(fasta.glob("*.fasta"))],
        "candidate_arms": {label: {"candidate_partition": record(partition)}
                           for label in ("p0_c0", "p0_c1", "p1_c0", "p1_c1")}}
    prepared = tmp_path / "manifest.json"
    prepared.write_text(json.dumps(manifest))
    admission_path = tmp_path / "admission.json"
    admission = {"prepared_manifest": record(prepared)}
    admission_path.write_text(json.dumps(admission))
    results = tmp_path / "benchmark_tools/results"
    results.mkdir(parents=True)
    (results / "qfo_assessment_environment_20260917.json").write_text(json.dumps({"reference_files": [record(mapping)]}))
    monkeypatch.setenv("SLURM_JOB_ID", "123")
    monkeypatch.setenv("SLURM_CPUS_PER_TASK", "2")
    monkeypatch.setattr(module, "verify_admission", lambda *a: (admission, manifest, cells[1], tmp_path, tmp_path, {"JobIDRaw": "122"}))
    monkeypatch.setattr(module, "read_frozen", lambda path, sha: json.loads(path.read_text()))
    return admission_path


def test_real_converter_and_filter_with_mocked_upstream_admission(tmp_path, monkeypatch):
    path = setup(tmp_path, monkeypatch)
    result = module.prepare(tmp_path, 0, path, record(path)["sha256"], "122")
    assert result["status"] == "corrected_factorial_group_pairs_prepared_unscored"
    assert result["total_pairs"] == result["retained_pairs"] == result["expected_pairs"] == 2
    assert result["removed_mapping_pairs"] == 0
    assert result["accuracy_evaluated"] is False
    assert result["participant"] == "ohmm_qfo_corrected_factorial_p0_c0_r0"
    assert open(result["filtered_pairs"]["path"]).read() == "A\tC\nB\tC\n"
    for item in [result["pairs"], result["filtered_pairs"], *result["checked_records"]]:
        check(item)
    with pytest.raises(FileExistsError):
        module.prepare(tmp_path, 0, path, record(path)["sha256"], "122")


def test_mapping_loss_is_preserved_failure(tmp_path, monkeypatch):
    path = setup(tmp_path, monkeypatch, missing_mapping=True)
    with pytest.raises(ValueError, match="mapping loss"):
        module.prepare(tmp_path, 0, path, record(path)["sha256"], "122")
    output = tmp_path / "benchmarks/results/qfo_corrected_factorial_pairs_v1/p0_c0_r0"
    result = json.loads((output / "results.json").read_text())
    assert result["status"] == "failed"
    assert (output / "pairs.partial.tsv").is_file()
    assert not (output / "pairs.tsv").exists()


def test_unscheduled_conversion_rejected(tmp_path, monkeypatch):
    monkeypatch.delenv("SLURM_JOB_ID", raising=False)
    with pytest.raises(ValueError, match="scheduled"):
        module.prepare(tmp_path, 0, tmp_path / "missing", "0" * 64, "122")
