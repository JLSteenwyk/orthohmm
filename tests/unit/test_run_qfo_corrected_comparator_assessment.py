import pytest

from benchmark_tools.run_qfo_corrected_comparator_assessment import validate_stage


def fixture():
    stage = {"status": "corrected_comparator_pairs_prepared_unscored", "accuracy_evaluated": False,
             "method": "proteinortho", "participant": "qfo_corrected_proteinortho", "job_id": "1",
             "total_pairs": 2, "retained_pairs": 2, "removed_mapping_pairs": 0,
             "pairs": {"bytes": 8, "sha256": "a"}, "filtered_pairs": {"bytes": 8, "sha256": "a"}}
    scheduler = {"State": "COMPLETED", "ExitCode": "0:0", "NodeList": "bizon", "AllocCPUS": "2", "JobIDRaw": "1"}
    return stage, scheduler


def test_valid():
    stage, scheduler = fixture()
    validate_stage(stage, "proteinortho", scheduler)


@pytest.mark.parametrize("key,value", [("status", "failed"), ("accuracy_evaluated", True),
    ("method", "sonic"), ("participant", "historical"), ("job_id", "2"),
    ("total_pairs", 0), ("total_pairs", True), ("retained_pairs", 1), ("removed_mapping_pairs", 1),
    ("filtered_pairs", {"bytes": 8, "sha256": "changed"})])
def test_reject_changed_stage(key, value):
    stage, scheduler = fixture()
    stage[key] = value
    with pytest.raises(ValueError):
        validate_stage(stage, "proteinortho", scheduler)


@pytest.mark.parametrize("key,value", [("State", "RUNNING"), ("ExitCode", "1:0"),
                                      ("NodeList", "other"), ("AllocCPUS", "1")])
def test_reject_scheduler(key, value):
    stage, scheduler = fixture()
    scheduler[key] = value
    with pytest.raises(ValueError, match="scheduler"):
        validate_stage(stage, "proteinortho", scheduler)


@pytest.mark.parametrize("exit_code", [0, 1])
def test_run_preserves_native_status(monkeypatch, tmp_path, exit_code):
    import json
    from types import SimpleNamespace
    from benchmark_tools import run_qfo_corrected_comparator_assessment as module
    source = tmp_path / "source"
    source.write_text("fixture")
    output = tmp_path / "execution"
    report = {"cwd": str(output), "results": str(tmp_path / "scores"), "command": ["fixture"],
              "source": module.record(source), "verified_records": [], "environment_overrides": {},
              "accuracy_admitted": False}
    monkeypatch.setattr(module, "prepare", lambda *args: report)
    monkeypatch.setenv("SLURM_CPUS_PER_TASK", "8")
    monkeypatch.setenv("SLURM_JOB_ID", "123")
    monkeypatch.setattr(module.subprocess, "run", lambda *args, **kwargs: SimpleNamespace(returncode=exit_code))
    if exit_code:
        with pytest.raises(RuntimeError, match="scoring failed"):
            module.run(tmp_path, "proteinortho", "fixture", 1)
    else:
        assert module.run(tmp_path, "proteinortho", "fixture", 1)["status"] == "process_succeeded_pending_independent_admission"
    final = json.loads((output / "results.json").read_text())
    assert final["exit_code"] == exit_code
    assert final["accuracy_admitted"] is False
    assert (output / "preflight.json").exists()
    assert final["status"] == ("failed" if exit_code else "process_succeeded_pending_independent_admission")
