import pytest

from benchmark_tools import admit_qfo_corrected_orthofinder as module
from tests.unit.test_admit_qfo_corrected_high_sensitivity import fixture as hmm_fixture


def fixture():
    plan, execution, scheduler = hmm_fixture()
    config = plan["methods"].pop("orthohmm_high_sensitivity")
    config["native_argv"] = ["/orthofinder", "-f", "/input", "-t", "32", "-a", "32", "-S", "diamond"]
    plan["methods"][module.METHOD] = config
    execution.update(method=module.METHOD, index=1, array_task_id="1", native_argv=config["native_argv"][:])
    for key, filename in (("log", "native.log"), ("timing", "time.txt")):
        execution[key]["path"] = f"/output/execution/{module.METHOD}/{filename}"
    return plan, execution, scheduler


def validate(plan, execution, scheduler):
    module.validate_execution(plan, execution, scheduler, {"sha256": "plan"}, {"sha256": "runner"})


def test_completed_execution():
    validate(*fixture())


@pytest.mark.parametrize("key,value", [
    ("status", "running"), ("exit_code", 1), ("method", "orthohmm_high_sensitivity"),
    ("index", 0), ("index", True), ("array_task_id", "0"), ("array_job_id", ""),
    ("job_id", "102"), ("node", "other"), ("source", {}), ("manifest", {}),
    ("native_argv", []), ("cwd", "/wrong"), ("started_epoch", 5.),
    ("finished_epoch", float("nan")), ("accuracy_admitted", True), ("native_outputs_validated", True),
    ("log", {"path": "/wrong", "bytes": 10}),
])
def test_reject_execution_change(key, value):
    plan, execution, scheduler = fixture()
    execution[key] = value
    with pytest.raises(ValueError):
        validate(plan, execution, scheduler)


@pytest.mark.parametrize("key,value", [("State", "RUNNING"), ("ExitCode", "1:0"),
    ("NodeList", "other"), ("AllocCPUS", "16"), ("JobIDRaw", "102")])
def test_reject_scheduler_change(key, value):
    plan, execution, scheduler = fixture()
    scheduler[key] = value
    with pytest.raises(ValueError):
        validate(plan, execution, scheduler)


def test_pending_job_rejected_before_file_reads(tmp_path, monkeypatch):
    monkeypatch.setattr(module.subprocess, "check_output", lambda *a, **k:
        "JobIDRaw|State|ExitCode|Elapsed|NodeList|AllocCPUS\n101|PENDING|0:0|0:00|bizon|32\n")
    output = tmp_path / "admission.json"
    with pytest.raises(ValueError, match="COMPLETED"):
        module.admit(tmp_path, 101, output)
    assert not output.exists()
