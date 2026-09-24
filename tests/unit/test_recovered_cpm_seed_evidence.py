import json
from pathlib import Path
import subprocess
import sys

import pytest

from benchmark_tools import recovered_cpm_seed_evidence as module
from benchmark_tools import prepare_recovered_cpm_candidates as candidates
from benchmark_tools.prepare_ob_candidate_neighborhood import record


@pytest.mark.parametrize("problem", [None, "pending", "failed", "cpu", "memory", "revision", "source",
    "status", "admitted", "accuracy", "native_job", "failure", "comparison", "seed", "stages",
    "records", "changed_seed"])
def test_recovered_seed_gate(tmp_path, monkeypatch, problem):
    executor = tmp_path / "benchmarks/work/cpm_checkpoint_admission_v1_20260923"
    source = executor / "benchmark_tools/admit_cpm_checkpoint_recovery.py"
    source.parent.mkdir(parents=True)
    source.write_text("fixture")
    monkeypatch.setattr(module, "SOURCE_SHA", "wrong" if problem == "source" else record(source)["sha256"])
    directory = tmp_path / "benchmarks/results/qfo_cpm_checkpoint_recovery_v1"
    original = tmp_path / "benchmarks/results/qfo_parameter_cpm_replay_v3/cpm_high/replay"
    directory.mkdir(parents=True)
    original.mkdir(parents=True)
    paths = [original / "orthogroups_multipass.txt", original / "orthogroups_multipass_refined.txt",
             directory / "orthogroups_profiles.txt", directory / "orthogroups_profiles_refined.txt"]
    for path in paths:
        path.write_text("a b\n")
    (directory / "status.json").write_text("{}")
    stages = [dict(label=label, origin=origin, output=record(path)) for path, label, origin in zip(paths,
        ("multipass", "multipass_refined", "strict_profiles", "strict_profiles_refined"),
        ("reused", "reused", "recovered", "recovered"))]
    report = dict(status="cpm_checkpoint_recovered_seed_admitted_unscored", source=record(source),
        seed_admitted=True, accuracy_evaluated=False, publication_ready=False, stages=stages,
        scheduler=dict(JobID="22154", State="COMPLETED", ExitCode="0:0", AllocCPUS="1", ReqMem="64G", NodeList="bizon"),
        comparison=dict(partition_equal=True, genes=984137), original_failure=dict(JobID="22081_1", State="FAILED"),
        seed_partition=record(paths[-1]), source_report=record(directory / "status.json"),
        checked_records=[record(source), record(directory / "status.json"), *[record(path) for path in paths]])
    if problem in ("status", "admitted", "accuracy", "seed", "records"):
        key, value = {"status": ("status", "running"), "admitted": ("seed_admitted", False),
            "accuracy": ("accuracy_evaluated", True), "seed": ("seed_partition", record(paths[2])),
            "records": ("checked_records", [])}[problem]
        report[key] = value
    elif problem == "native_job":
        report["scheduler"]["JobID"] = "22081"
    elif problem == "failure":
        report["original_failure"]["State"] = "COMPLETED"
    elif problem == "comparison":
        report["comparison"]["partition_equal"] = False
    elif problem == "stages":
        report["stages"].reverse()
    path = tmp_path / "benchmarks/results/qfo_cpm_checkpoint_recovery_admission_v1/status.json"
    path.parent.mkdir()
    path.write_text(json.dumps(report))
    if problem == "changed_seed":
        paths[-1].write_text("changed")
    def check_output(command, **kwargs):
        if command[0] == "git":
            return "changed" if problem == "revision" else module.COMMIT
        state = {"pending": "PENDING", "failed": "FAILED"}.get(problem, "COMPLETED")
        cpu = "1" if problem == "cpu" else "2"
        memory = "32G" if problem == "memory" else "64G"
        return f"JobID|State|ExitCode|AllocCPUS|ReqMem|NodeList\n22155|{state}|0:0|{cpu}|{memory}|bizon\n"
    monkeypatch.setattr(module.subprocess, "check_output", check_output)
    monkeypatch.setattr(module.subprocess, "run", lambda *a, **k: None)
    if problem:
        with pytest.raises(ValueError):
            module.evidence(tmp_path)
    else:
        result = module.evidence(tmp_path)
        assert result["seed_partition"] == record(paths[-1])
        assert result["report"] == report


def test_candidate_preparation_rejects_pending_before_scientific_imports(tmp_path, monkeypatch):
    for name, value in dict(SLURM_JOB_ID="fixture", SLURM_CPUS_PER_TASK="2", SLURM_MEM_PER_NODE="65536",
        SLURM_JOB_NODELIST="bizon", PYTHONHASHSEED="0", OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1").items():
        monkeypatch.setenv(name, value)
    def pending(root):
        raise ValueError("admission pending")
    monkeypatch.setattr(candidates, "evidence", pending)
    with pytest.raises(ValueError, match="pending"):
        candidates.prepare(tmp_path)
    assert not (tmp_path / "benchmarks").exists()


def test_candidate_cli_help(tmp_path):
    done = subprocess.run([sys.executable, "-I", str(Path(candidates.__file__).resolve()), "--help"],
                          cwd=tmp_path, capture_output=True, text=True)
    assert done.returncode == 0, done.stderr


def test_existing_builder_bytes_are_pinned():
    assert record(Path(candidates.__file__).with_name("prepare_qfo_cpm_candidates.py"))["sha256"] == candidates.BUILD_SHA
