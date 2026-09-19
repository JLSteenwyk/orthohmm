import json
from pathlib import Path
import subprocess
import sys
from types import SimpleNamespace

import pytest

from benchmark_tools import run_qfo_cpm_variant as module
from benchmark_tools.prepare_ob_candidate_neighborhood import record
from tests.unit.test_run_qfo_corrected_replay import fixture, STAGES


@pytest.mark.parametrize("index", [-1, 2, True, "0"])
def test_unknown_index_rejected(tmp_path, index):
    with pytest.raises(ValueError, match="index"):
        module.run(tmp_path, index)


def test_allocation_required(tmp_path, monkeypatch):
    monkeypatch.delenv("SLURM_JOB_ID", raising=False)
    with pytest.raises(ValueError, match="scheduled"):
        module.run(tmp_path, 0)


@pytest.mark.parametrize("state", ["RUNNING", "PENDING", "FAILED"])
def test_control_scheduler_gate_precedes_files(tmp_path, monkeypatch, state):
    monkeypatch.setattr(module.subprocess, "check_output", lambda *a, **k:
        f"JobIDRaw|State|ExitCode|Elapsed|NodeList|AllocCPUS\n21958|{state}|0:0|00:01:00|bizon|2\n")
    with pytest.raises(ValueError, match="COMPLETED"):
        module.control_evidence(tmp_path, {})


@pytest.mark.parametrize("problem", [None, "cpus", "revision", "source", "context", "authorization", "input"])
def test_control_evidence(tmp_path, monkeypatch, problem):
    executor = tmp_path / "benchmarks/work/publication_qfo_cpm_control_admission_v1"
    source = executor / "benchmark_tools/admit_qfo_cpm_control.py"
    source.parent.mkdir(parents=True)
    source.write_text("fixture")
    monkeypatch.setattr(module, "ADMISSION_SHA", record(source)["sha256"] if problem != "source" else "wrong")
    data = tmp_path / "input"
    data.write_text("fixture")
    context = {"arm": "control"}
    report = {"status": "cpm_control_reproduced_and_admitted", "source": record(source),
        "context": {} if problem == "context" else context,
        "changed_arms_authorized": [] if problem == "authorization" else list(module.ARMS),
        "accuracy_evaluated": False, "publication_ready": False, "checked_records": [record(data)]}
    path = tmp_path / "benchmarks/work/qfo_cpm_control_admission_21958.json"
    path.write_text(json.dumps(report))
    if problem == "input":
        data.write_text("changed")
    def check_output(command, **kwargs):
        if command[0] == "sacct":
            cpu = "32" if problem == "cpus" else "2"
            return f"JobIDRaw|State|ExitCode|Elapsed|NodeList|AllocCPUS\n21958|COMPLETED|0:0|00:01:00|bizon|{cpu}\n"
        return "changed" if problem == "revision" else module.ADMISSION_COMMIT
    monkeypatch.setattr(module.subprocess, "check_output", check_output)
    monkeypatch.setattr(module.subprocess, "run", lambda *a, **k: None)
    if problem:
        with pytest.raises(ValueError):
            module.control_evidence(tmp_path, context)
    else:
        result = module.control_evidence(tmp_path, context)
        assert result["report"] == report
        assert result["report_record"] == record(path)


@pytest.mark.parametrize("index", [0, 1])
def test_worker_entry_does_not_preload_scientific_package(tmp_path, index):
    directory = tmp_path / "benchmark_tools"
    directory.mkdir()
    script = directory / "run_qfo_cpm_variant.py"
    script.write_text(Path(module.__file__).read_text())
    (directory / "run_qfo_cpm_control.py").write_text(
        "import sys\ndef worker(root, arm):\n"
        "    assert 'benchmark_tools' not in sys.modules\n"
        "    (root / 'selected_arm').write_text(arm)\n")
    subprocess.run([sys.executable, "-I", str(script), "--root", str(tmp_path),
                    "--index", str(index), "--worker"], check=True)
    assert (tmp_path / "selected_arm").read_text() == module.ARMS[index]


@pytest.mark.parametrize("outcome", ["success", "authorization", "fresh", "failed", "context", "runtime", "coverage"])
def test_run_preserves_failures_and_requires_fresh_control(tmp_path, monkeypatch, outcome):
    for key, value in {"SLURM_JOB_ID": "test", "SLURM_CPUS_PER_TASK": "32", "SLURM_JOB_NODELIST": "bizon",
                       "SLURM_ARRAY_TASK_ID": "0"}.items():
        monkeypatch.setenv(key, value)
    output = tmp_path / "output"
    names = tmp_path / "names.txt"
    names.write_text("".join(f"g{i}\n" for i in range(984137)))
    plan = {"runtime": {}}
    context = {"output_root": str(output), "cwd": str(tmp_path), "checked_records": [],
               "environment_overrides": {"OMP_NUM_THREADS": "1"}, "expected_stages": STAGES}
    monkeypatch.setattr("benchmark_tools.checked_replay_payload_worker.corrected_evidence",
                        lambda *a: (plan, {}, {}, record(names)))
    monkeypatch.setattr("benchmark_tools.cpm_replay_context.evidence", lambda *a: context)
    authorization = {"executor": str(tmp_path), "report": {"authorized": True}, "checked_records": []}
    def control(*args):
        if outcome == "authorization":
            raise ValueError("Missing authorization")
        return authorization
    monkeypatch.setattr(module, "control_evidence", control)
    calls = []
    def verify(*args):
        calls.append(args)
        return {"changed": True} if outcome == "runtime" and len(calls) > 1 else {}
    monkeypatch.setattr("benchmark_tools.verify_qfo_replay_launcher.verify", verify)
    def read_partition(path, universe):
        assert len(universe) == 984137
        if outcome == "coverage":
            raise ValueError("Incomplete partition")
        return [universe]
    monkeypatch.setattr("benchmark_tools.audit_historical_profile_ablation.read_partition", read_partition)
    monkeypatch.setattr(module.subprocess, "check_output", lambda *a, **k: "frozen\n")
    launches = []
    def execute(command, **kwargs):
        launches.append(command)
        if "--output" in command:
            Path(command[-1]).write_text(json.dumps({} if outcome == "fresh" else authorization["report"]))
            return SimpleNamespace(returncode=0)
        assert (output / "fresh_control_admission.json").exists()
        assert kwargs["env"]["OMP_NUM_THREADS"] == "1"
        if outcome == "failed":
            return SimpleNamespace(returncode=1)
        worker, replay = fixture()
        worker["context"] = {} if outcome == "context" else context
        for i, stage in enumerate(replay["stages"]):
            path = output / f"partition{i}.txt"
            path.write_text("fixture")
            stage["output"] = record(path)
        (output / "checked_worker.json").write_text(json.dumps(worker))
        (output / "replay.json").write_text(json.dumps(replay))
        return SimpleNamespace(returncode=0)
    monkeypatch.setattr(module.subprocess, "run", execute)
    if outcome == "success":
        result = module.run(tmp_path, 0)
        assert result["status"] == "cpm_variant_replay_complete_pending_admission"
        assert result["accuracy_evaluated"] is False and result["publication_ready"] is False
        with pytest.raises(FileExistsError):
            module.run(tmp_path, 0)
    else:
        with pytest.raises((ValueError, RuntimeError)):
            module.run(tmp_path, 0)
        if outcome == "authorization":
            assert not output.exists() and not launches
        else:
            result = json.loads((output / "results.json").read_text())
            assert result["status"] == "failed" and result["accuracy_evaluated"] is False
            assert len(launches) == (1 if outcome == "fresh" else 2)
