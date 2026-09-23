import gzip
import json
from pathlib import Path

import pytest

from benchmark_tools import prepare_qfo_cpm_pairs as module
from benchmark_tools.prepare_ob_candidate_neighborhood import record


def fixture():
    pairs = {"path": "native.tsv", "sha256": "abc", "bytes": 10}
    verified = {"context": {"arm": "cpm_low"}, "admission_record": {"path": "candidate.json"}}
    native = {"status": "cpm_native_pairs_verified_unscored", "arm": "cpm_low", "index": 0,
        "cell": {"label": "candidate_cpm_low"}, "source": {"path": "validator"},
        "context": verified["context"], "candidate_admission": verified["admission_record"],
        "accuracy_evaluated": False, "scoring_admitted": False, "publication_ready": False,
        "native_pair_count": 2, "native_pairs": pairs, "checked_records": [pairs]}
    return native, verified


@pytest.mark.parametrize("change", [None, ("status", "prepared"), ("arm", "control"), ("index", 1),
    ("index", False), ("cell", {"label": "p1_c1_r1"}), ("source", {}), ("context", {}),
    ("candidate_admission", {}), ("accuracy_evaluated", True), ("scoring_admitted", True),
    ("publication_ready", True), ("native_pair_count", True), ("native_pair_count", 0), ("checked_records", [])])
def test_native_identity(change):
    native, verified = fixture()
    if change:
        native[change[0]] = change[1]
        with pytest.raises(ValueError):
            module.validate_native(native, 0, {"path": "validator"}, verified)
    else:
        module.validate_native(native, 0, {"path": "validator"}, verified)


@pytest.mark.parametrize("index", [-1, 2, True, "0", 0.0])
def test_invalid_index(index):
    with pytest.raises(ValueError):
        module.validate_native({}, index, {}, {})
    with pytest.raises(ValueError):
        module.completed_admission("", index)


def accounting(state="COMPLETED", cpus="2", exit_code="0:0", node="bizon"):
    return ("JobID|JobIDRaw|State|ExitCode|NodeList|AllocCPUS\n"
            f"22090_0|22091|{state}|{exit_code}|{node}|{cpus}\n")


@pytest.mark.parametrize("change", [{"state": "RUNNING"}, {"state": "FAILED"}, {"state": "PENDING"},
    {"cpus": "32"}, {"exit_code": "1:0"}, {"node": "spark-7ff0"}])
def test_scheduler_gate(change):
    with pytest.raises(ValueError):
        module.completed_admission(accounting(**change), 0)


def test_scheduler_unique_array_identity():
    assert module.completed_admission(accounting(), 0)["JobIDRaw"] == "22091"
    with pytest.raises(ValueError):
        module.completed_admission(accounting(), 1)
    with pytest.raises(ValueError):
        module.completed_admission(accounting() + accounting().splitlines()[1] + "\n", 0)


def test_allocation_precedes_input_reads(tmp_path, monkeypatch):
    monkeypatch.delenv("SLURM_JOB_ID", raising=False)
    monkeypatch.setattr(module.subprocess, "check_output", lambda *a, **k: pytest.fail("Read scheduler"))
    with pytest.raises(ValueError, match="scheduled"):
        module.prepare(tmp_path, 0)


@pytest.mark.parametrize("problem", [None, "fresh", "mapping", "changed", "source", "revision", "live"])
def test_prepare_real_conversion_with_mocked_native_recheck(tmp_path, monkeypatch, problem):
    executor = tmp_path / "benchmarks/work/publication_qfo_cpm_native_admission_v3"
    checker = executor / "benchmark_tools/admit_qfo_cpm_phylogeny.py"
    checker.parent.mkdir(parents=True)
    checker.write_text("# frozen fixture\n")
    monkeypatch.setattr(module, "ADMITTER_SHA", "wrong" if problem == "source" else record(checker)["sha256"])
    pairs = tmp_path / "native.tsv"
    pairs.write_text("gene_a\tspecies_a\tgene_b\tspecies_b\nA\ts1\tB\ts2\nA\ts1\tC\ts3\n")
    metadata = tmp_path / "metadata.json"
    metadata.write_text("{}")
    native, verified = fixture()
    native.update(source=record(checker), native_pairs=record(pairs), checked_records=[record(pairs)],
        helpers=[], native_group_integrity={"native_manifest": record(metadata)})
    verified.update(arm={"partition": {"path": "CPM_candidate"}}, manifest={"input_fastas": []}, checked_records=[])
    admission = tmp_path / "benchmarks/work/qfo_cpm_native_admission_22090_0.json"
    admission.write_text(json.dumps(native))
    mapping = tmp_path / "mapping.json.gz"
    with gzip.open(mapping, "wt") as stream:
        json.dump({"mapping": dict.fromkeys("AB" if problem == "mapping" else "ABC", 1)}, stream)
    env_path = tmp_path / "benchmark_tools/results/qfo_assessment_environment_20260917.json"
    env_path.parent.mkdir(parents=True)
    env_path.write_text(json.dumps({"reference_files": [record(mapping)]}))
    monkeypatch.setattr(module, "ENV_SHA", record(env_path)["sha256"])
    checks = []
    def verify(*args):
        checks.append(args)
        return {} if problem == "changed" and len(checks) > 1 else verified
    monkeypatch.setattr(module, "verify_sources", verify)
    monkeypatch.setattr(module, "gene_ownership", lambda *a: ({"A": "s1", "B": "s2", "C": "s3"}, {}))
    for key, value in {"SLURM_JOB_ID": "23100", "SLURM_CPUS_PER_TASK": "2",
                       "SLURM_JOB_NODELIST": "bizon", "SLURM_ARRAY_TASK_ID": "0"}.items():
        monkeypatch.setenv(key, value)
    monkeypatch.setattr(module.subprocess, "check_output", lambda command, **k:
        accounting(state="RUNNING" if problem == "live" else "COMPLETED") if command[0] == "sacct"
        else "wrong" if problem == "revision" else module.ADMITTER_COMMIT)
    def run(command, **kwargs):
        if command[0] != "git":
            assert command[1:3] == ["-B", str(checker)]
            assert command[command.index("--index") + 1] == "0"
            Path(command[command.index("--output") + 1]).write_text(json.dumps({} if problem == "fresh" else native))
    monkeypatch.setattr(module.subprocess, "run", run)
    output = tmp_path / "benchmarks/results/qfo_cpm_pairs_v1/cpm_low"
    if problem:
        with pytest.raises(ValueError):
            module.prepare(tmp_path, 0)
        if problem in {"source", "revision", "live"}:
            assert not output.exists()
            return
        report = json.loads((output / "results.json").read_text())
        assert report["status"] == "failed" and not (output / "pairs.tsv").exists()
        if problem == "mapping":
            counts = json.loads((output / "conversion_counts.json").read_text())
            assert counts["removed_mapping_pairs"] == 1
    else:
        report = module.prepare(tmp_path, 0)
        assert report["status"] == "cpm_native_pairs_prepared_unscored"
        assert report["native_admission_recheck"] in report["checked_records"]
        assert (output / "pairs.qfo.tsv").read_text() == "A\tB\nA\tC\n"
        assert report["written_pairs"] == report["total_pairs"] == report["retained_pairs"] == 2
        assert report["removed_mapping_pairs"] == 0 and report["accuracy_evaluated"] is False
        assert len(checks) == 2
    monkeypatch.setattr(module, "verify_sources", lambda *a: verified)
    with pytest.raises(FileExistsError):
        module.prepare(tmp_path, 0)
