from copy import deepcopy

import pytest

from benchmark_tools.prepare_qfo_corrected_factorial import prepare, validate_replay_admission


def fixture():
    plan, source, report = {"plan": True}, {"source": True}, {"report": True}
    stages = [{"label": label, "output": {"path": label, "sha256": str(i)}} for i, label in enumerate(
        ("multipass", "multipass_refined", "strict_profiles", "strict_profiles_refined"))]
    admission = {"status": "corrected_checked_replay_admitted", "accuracy_evaluated": False,
        "publication_ready": False, "source": source, "plan": plan, "coverage": stages,
        "source_report": report, "checked_records": [report, *[s["output"] for s in stages]],
        "clustering": [{"genes": 984137} for _ in range(4)], "native_partition_comparison": {"partition_equal": False},
        "scheduler": {"State": "COMPLETED", "ExitCode": "0:0", "NodeList": "bizon", "AllocCPUS": "32"}}
    return admission, plan, source


@pytest.mark.parametrize("problem", [None, "equal", "status", "accuracy", "source", "plan", "scheduler", "cpus",
    "stages", "order", "unchecked", "report", "missing_cluster", "old_genes", "comparison"])
def test_corrected_replay_binding(problem):
    admission, plan, source = deepcopy(fixture())
    if problem == "equal":
        admission["native_partition_comparison"]["partition_equal"] = True
    elif problem == "status":
        admission["status"] = "checked_full_replay_verified"
    elif problem == "accuracy":
        admission["accuracy_evaluated"] = True
    elif problem == "source":
        admission["source"] = {}
    elif problem == "plan":
        admission["plan"] = {}
    elif problem == "scheduler":
        admission["scheduler"]["State"] = "RUNNING"
    elif problem == "cpus":
        admission["scheduler"]["AllocCPUS"] = "16"
    elif problem == "stages":
        admission["coverage"].pop()
    elif problem == "order":
        admission["coverage"].reverse()
    elif problem == "unchecked":
        admission["checked_records"].pop()
    elif problem == "report":
        admission["checked_records"].pop(0)
    elif problem == "missing_cluster":
        admission["clustering"].pop()
    elif problem == "old_genes":
        admission["clustering"][0]["genes"] = 976504
    elif problem == "comparison":
        admission["native_partition_comparison"]["partition_equal"] = "unknown"
    if problem not in (None, "equal"):
        with pytest.raises(ValueError):
            validate_replay_admission(admission, plan, source)
    else:
        seeds = validate_replay_admission(admission, plan, source)
        assert seeds == [(False, admission["coverage"][1]["output"]), (True, admission["coverage"][3]["output"])]


@pytest.mark.parametrize("problem", ["output", "environment", "unscheduled", "cpus"])
def test_prepare_rejects_before_input_access(tmp_path, monkeypatch, problem):
    for key in ("PYTHONHASHSEED", "OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
        monkeypatch.setenv(key, "0" if key == "PYTHONHASHSEED" else "1")
    monkeypatch.setenv("SLURM_JOB_ID", "fixture")
    monkeypatch.setenv("SLURM_CPUS_PER_TASK", "2")
    output = tmp_path / "out"
    if problem == "output":
        output.mkdir()
    elif problem == "environment":
        monkeypatch.setenv("OMP_NUM_THREADS", "32")
    elif problem == "unscheduled":
        monkeypatch.delenv("SLURM_JOB_ID")
    else:
        monkeypatch.setenv("SLURM_CPUS_PER_TASK", "32")
    with pytest.raises((ValueError, FileExistsError)):
        prepare(tmp_path, tmp_path / "missing.json", "missing", "1", output)
    assert not (output / "manifest.json").exists()
