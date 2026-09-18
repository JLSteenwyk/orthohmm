from pathlib import Path
import subprocess

import pytest

from benchmark_tools import admit_qfo_sequence_search_control as module
from benchmark_tools.prepare_ob_candidate_neighborhood import record


def fixture(tmp_path):
    executor = tmp_path / "executor"
    (executor / "benchmark_tools").mkdir(parents=True)
    for name in ("run_qfo_sequence_search_control.py", "run_sequence_search_control.py"):
        (executor / "benchmark_tools" / name).write_text("# test source\n")
    manifest = {"path": "/manifest", "bytes": 1, "sha256": "fixture"}
    plan, targets = {"searches": []}, []
    for i in range(78):
        directory = Path(f"/targets/{i}")
        expected = {"index": i, "output": str(directory / "hits.tsv"),
                    "makedb": ["diamond", "makedb", str(i)], "search": ["diamond", "blastp", str(i)]}
        plan["searches"].append(expected)
        phases = {p: {"argv": expected[p], "exit_code": 0, "wall_s": 1.,
                      "log": {"path": str(directory / (p + ".log"))},
                      "gnu_time": {"path": str(directory / (p + ".time.log"))}} for p in ("makedb", "search")}
        targets.append({"index": i, "status": "complete_pending_numeric_validation", "phases": phases,
                        "hits": {"path": expected["output"]}, "database": {"path": str(directory / "target.dmnd")}})
    report = {"status": "complete_pending_numeric_validation", "job_id": "123", "node": "bizon",
        "manifest": manifest, "executor": record(executor / "benchmark_tools/run_qfo_sequence_search_control.py"),
        "phase_helper": record(executor / "benchmark_tools/run_sequence_search_control.py"),
        "time_binary": record(Path("/usr/bin/time")), "targets": targets,
        "environment_overrides": {"PYTHONHASHSEED": "0", "OMP_NUM_THREADS": "1", "OPENBLAS_NUM_THREADS": "1", "MKL_NUM_THREADS": "1"},
        "accuracy_evaluated": False, "numeric_validated": False, "publication_ready": False}
    scheduler = {"State": "COMPLETED", "ExitCode": "0:0", "NodeList": "bizon", "AllocCPUS": "32",
                 "ReqMem": "192G", "JobIDRaw": "123"}
    return report, plan, scheduler, manifest, executor


def test_complete_panel(tmp_path):
    records = module.validate_execution(*fixture(tmp_path))
    assert len(records) == 4 + 78 * 6


@pytest.mark.parametrize("problem", ["missing_target", "order", "failed_target", "phase_exit", "phase_boolean",
    "command", "log", "hits", "database", "nan", "negative", "environment", "executor", "job", "status",
    "numeric", "scheduler", "memory", "cpu"])
def test_invalid_panel(tmp_path, problem):
    report, plan, scheduler, manifest, executor = fixture(tmp_path)
    target = report["targets"][0]
    phase = target["phases"]["search"]
    if problem == "missing_target":
        report["targets"].pop()
    elif problem == "order":
        report["targets"].reverse()
    elif problem == "failed_target":
        target["status"] = "failed"
    elif problem in ("phase_exit", "phase_boolean"):
        phase["exit_code"] = 1 if problem == "phase_exit" else False
    elif problem == "command":
        phase["argv"] = ["different"]
    elif problem == "log":
        phase["log"]["path"] = "/wrong"
    elif problem in ("hits", "database"):
        target[problem]["path"] = "/wrong"
    elif problem in ("nan", "negative"):
        phase["wall_s"] = float("nan") if problem == "nan" else -1
    elif problem == "environment":
        report["environment_overrides"]["OMP_NUM_THREADS"] = "32"
    elif problem == "executor":
        report["executor"]["sha256"] = "wrong"
    elif problem == "job":
        report["job_id"] = "999"
    elif problem == "status":
        report["status"] = "running"
    elif problem == "numeric":
        report["numeric_validated"] = True
    else:
        scheduler[{"scheduler": "State", "memory": "ReqMem", "cpu": "AllocCPUS"}[problem]] = "wrong"
    with pytest.raises(ValueError):
        module.validate_execution(report, plan, scheduler, manifest, executor)


def test_live_job_refused_before_reading_partial_output(tmp_path, monkeypatch):
    monkeypatch.setattr(module.subprocess, "check_output", lambda *a, **k:
        "JobIDRaw|State|ExitCode|Elapsed|NodeList|AllocCPUS|ReqMem\n123|RUNNING|0:0|01:00|bizon|32|192G\n")
    monkeypatch.setattr(module, "verify_plan", lambda *a: pytest.fail("Read partial evidence"))
    with pytest.raises(ValueError):
        module.admit(tmp_path, "123", tmp_path / "admission.json")
    assert not (tmp_path / "admission.json").exists()


def test_batch_syntax_resources():
    path = Path(module.__file__).parent / "results/qfo_sequence_search_admit_batch_20260918.sh"
    subprocess.run(["bash", "-n", str(path)], check=True)
    for setting in ("--cpus-per-task=2", "--mem=64G", "--time=24:00:00", "--no-requeue", "--nodelist=bizon"):
        assert "#SBATCH " + setting in path.read_text()
