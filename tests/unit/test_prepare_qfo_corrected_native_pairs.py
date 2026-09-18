import gzip
import json
from pathlib import Path

import pytest

import benchmark_tools.prepare_qfo_corrected_native_pairs as module
from benchmark_tools.prepare_ob_candidate_neighborhood import record


def native_fixture():
    pair = {"path": "/native/pairs", "sha256": "pairs", "bytes": 1}
    candidate, prepared, source = ({"path": "/" + key, "sha256": key, "bytes": 1}
                                   for key in ("candidate", "prepared", "source"))
    return {"status": "corrected_qfo_native_pair_output_verified", "cell": "p0_c0_r1", "index": 0,
        "accuracy_evaluated": False, "scoring_admitted": False, "publication_ready": False,
        "candidate_admission": candidate, "prepared": prepared, "source": source,
        "native_pairs": pair, "native_pair_count": 2, "checked_records": [pair]}, candidate, prepared, source


def test_admission_identity():
    native, candidate, prepared, source = native_fixture()
    module.validate_native(native, 1, candidate, prepared, source)


@pytest.mark.parametrize("key,value", [("status", "historical"), ("cell", "p1_c0_r1"), ("index", 1),
    ("accuracy_evaluated", True), ("scoring_admitted", True), ("publication_ready", True),
    ("candidate_admission", {}), ("prepared", {}), ("source", {}), ("native_pair_count", True),
    ("native_pair_count", 0), ("checked_records", [])])
def test_reject_misbound_admission(key, value):
    native, candidate, prepared, source = native_fixture()
    native[key] = value
    with pytest.raises(ValueError):
        module.validate_native(native, 1, candidate, prepared, source)


@pytest.mark.parametrize("index", [0, 2, 4, 6, True, -1, 8])
def test_wrong_semantics_index(index):
    with pytest.raises(ValueError):
        module.validate_native({}, index, {}, {}, {})


@pytest.mark.parametrize("mismatch", [False, True])
def test_real_pair_conversion_with_mocked_native_recheck(tmp_path, monkeypatch, mismatch):
    candidate = tmp_path / "candidate.json"
    candidate.write_text("{}")
    prepared = tmp_path / "prepared.json"
    prepared.write_text("{}")
    partition = tmp_path / "groups.txt"
    partition.write_text("A B C\n")
    pairs = tmp_path / "native_pairs.tsv"
    pairs.write_text("gene_a\tspecies_a\tgene_b\tspecies_b\nA\ts1\tB\ts2\nA\ts1\tC\ts3\n")
    metadata = tmp_path / "native_manifest.json"
    metadata.write_text("{}")
    executor = tmp_path / "benchmarks/work/publication_qfo_corrected_reconcile_admission_v1"
    checker = executor / "benchmark_tools/admit_qfo_corrected_factorial_cell.py"
    checker.parent.mkdir(parents=True)
    checker.write_text("# frozen checker fixture\n")
    native, _, _, _ = native_fixture()
    native.update(candidate_admission=record(candidate), prepared=record(prepared), source=record(checker),
        native_pairs=record(pairs), checked_records=[record(pairs)],
        scheduler={"JobIDRaw": "121"}, native_group_integrity={"native_manifest": record(metadata)})
    native_path = tmp_path / "native_admission.json"
    native_path.write_text(json.dumps(native))
    environment = tmp_path / "benchmark_tools/results/qfo_assessment_environment_20260917.json"
    environment.parent.mkdir(parents=True)
    mapping = tmp_path / "mapping.json.gz"
    with gzip.open(mapping, "wt") as stream:
        json.dump({"mapping": dict.fromkeys("ABC", 1)}, stream)
    environment.write_text(json.dumps({"reference_files": [record(mapping)]}))
    cell = {"label": "p0_c0_r1", "reconciliation": True, "candidate_partition": str(partition)}
    monkeypatch.setenv("SLURM_JOB_ID", "123")
    monkeypatch.setenv("SLURM_CPUS_PER_TASK", "2")
    monkeypatch.setattr(module, "verify_admission", lambda *a: (
        {"prepared_manifest": record(prepared)}, {"input_fastas": []}, cell, tmp_path, tmp_path, {}))
    monkeypatch.setattr(module, "read_frozen", lambda path, sha: json.loads(path.read_text()))
    monkeypatch.setattr(module, "gene_ownership", lambda *a: ({"A": "s1", "B": "s2", "C": "s3"}, dict.fromkeys("ABC", 0)))
    monkeypatch.setattr(module.subprocess, "check_output", lambda *a, **k: module.ADMITTER_COMMIT + "\n")
    calls = []
    def run(argv, **kwargs):
        calls.append(argv)
        if argv[0] == "git":
            return
        assert argv[1] == str(checker)
        assert argv[argv.index("--job") + 1] == "121"
        assert argv[argv.index("--index") + 1] == "0"
        Path(argv[argv.index("--output") + 1]).write_text(json.dumps({} if mismatch else native))
    monkeypatch.setattr(module.subprocess, "run", run)
    args = (tmp_path, 1, candidate, record(candidate)["sha256"], "120", native_path, record(native_path)["sha256"])
    output = tmp_path / "benchmarks/results/qfo_corrected_factorial_pairs_v1/p0_c0_r1"
    if mismatch:
        with pytest.raises(ValueError, match="Fresh independent"):
            module.prepare(*args)
        assert json.loads((output / "results.json").read_text())["status"] == "failed"
        assert not (output / "pairs.tsv").exists()
    else:
        result = module.prepare(*args)
        assert result["status"] == "corrected_factorial_native_pairs_prepared_unscored"
        assert result["total_pairs"] == result["retained_pairs"] == 2
        assert result["removed_mapping_pairs"] == 0
        assert (output / "pairs.qfo.tsv").read_text() == "A\tB\nA\tC\n"
        assert "B\tC" not in (output / "pairs.qfo.tsv").read_text()
        with pytest.raises(FileExistsError):
            module.prepare(*args)
    assert any("--output" in argv for argv in calls)
