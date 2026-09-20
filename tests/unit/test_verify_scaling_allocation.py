import json
from pathlib import Path

import pytest

from benchmark_tools import verify_scaling_allocation as module
from benchmark_tools.capture_array_scheduler import TERMINAL

PLAN = Path(__file__).resolve().parents[2] / "benchmark_tools/results/dgx_root_context_scaling_plan_v2_20260920.json"


def raw(**changes):
    fields = dict(JobId="12345", Partition="spark", JobState="RUNNING", ExitCode="0:0",
        Restarts="0", Requeue="0", NodeList="spark-7ff0", OverSubscribe="NO",
        MinMemoryNode="96G", NumNodes="1", NumCPUs="20", NumTasks="1",
        TimeLimit="1-00:00:00", Command=str(module.SUBMISSION_SCRIPT),
        WorkDir=str(module.RECIPE_ROOT))
    fields["CPUs/Task"] = "20"
    fields.update(changes)
    return " ".join(f"{k}={v}" for k, v in fields.items()) + "\n"


@pytest.mark.parametrize("state", sorted(TERMINAL))
def test_terminal_failure_states_retained(state):
    result = module.validate(raw(JobState=state, ExitCode="1:9"), 12345, "terminal")
    assert result["fields"]["JobState"] == state
    assert result["fields"]["ExitCode"] == "1:9"
    assert result["scheduler_terminal_verified"] is True
    for key in ("execution_authorized", "next_submission_authorized", "automatic_retry", "scientific_timings_admitted"):
        assert result[key] is False


@pytest.mark.parametrize("limit", ["1-00:00:00", "24:00:00"])
def test_running_is_not_terminal_or_authorized(limit):
    result = module.validate(raw(TimeLimit=limit), 12345, "running")
    assert result["scheduler_terminal_verified"] is False
    assert result["execution_authorized"] is False


@pytest.mark.parametrize("key,value", [
    ("JobId", "12346"), ("Partition", "other"), ("Restarts", "1"), ("Requeue", "1"),
    ("NodeList", "bizon"), ("OverSubscribe", "YES"), ("MinMemoryNode", "128G"),
    ("NumNodes", "2"), ("NumCPUs", "32"), ("NumTasks", "2"), ("CPUs/Task", "10"),
    ("TimeLimit", "01:00:00"), ("TimeLimit", "UNLIMITED"), ("Command", "/tmp/other.sh"),
    ("WorkDir", "/tmp"), ("ArrayJobId", "12300"), ("ArrayTaskId", "0"),
    ("HetJobId", "12300"), ("HetJobOffset", "0"), ("ExitCode", "bad"),
    ("JobState", "PENDING"), ("JobState", "COMPLETING"), ("JobState", "COMPLETED"),
])
def test_different_allocation_rejected(key, value):
    with pytest.raises(ValueError):
        module.validate(raw(**{key: value}), 12345, "running")


@pytest.mark.parametrize("job", [True, 0, -1, "12345", 12345.0])
def test_job_type_rejected(job):
    with pytest.raises(ValueError):
        module.validate(raw(), job, "running")


@pytest.mark.parametrize("content", ["", raw() + raw(), raw().strip() + " NumCPUs=20\n"])
def test_ambiguous_record_rejected(content):
    with pytest.raises(ValueError):
        module.validate(content, 12345, "running")


def test_running_cannot_be_terminal():
    with pytest.raises(ValueError, match="terminal"):
        module.validate(raw(), 12345, "terminal")


@pytest.mark.parametrize("phase", ["", "complete", "observation_timeout", None])
def test_invalid_phase(phase):
    with pytest.raises(ValueError, match="phase"):
        module.validate(raw(), 12345, phase)


@pytest.mark.parametrize("key", ["JobId", "Partition", "JobState", "ExitCode", "Restarts",
    "Requeue", "NodeList", "OverSubscribe", "MinMemoryNode", "NumNodes", "NumCPUs",
    "NumTasks", "CPUs/Task", "TimeLimit", "Command", "WorkDir"])
def test_required_field_missing(key):
    text = " ".join(item for item in raw().split() if not item.startswith(key + "="))
    with pytest.raises(ValueError):
        module.validate(text, 12345, "running")


@pytest.mark.parametrize("index", range(27))
def test_all_frozen_slots_without_claiming_task_binding(tmp_path, index):
    scheduler = tmp_path / "scheduler.txt"
    scheduler.write_text(raw())
    result = module.audit(PLAN, index, scheduler, 12345, "running")
    assert result["requested_index"] == index
    assert result["task_identity_bound"] is False
    assert result["execution_authorized"] is False


def test_cli_and_no_overwrite(tmp_path):
    scheduler, output = tmp_path / "scheduler.txt", tmp_path / "result.json"
    scheduler.write_text(raw(JobState="FAILED", ExitCode="1:0"))
    argv = ["--plan", str(PLAN), "--index", "0", "--scheduler", str(scheduler),
            "--job", "12345", "--phase", "terminal", "--output", str(output)]
    assert module.main(argv) == 0
    contents = output.read_bytes()
    assert json.loads(contents)["fields"]["JobState"] == "FAILED"
    with pytest.raises(FileExistsError):
        module.main(argv)
    assert output.read_bytes() == contents


def test_changed_evidence_rejected(tmp_path, monkeypatch):
    scheduler = tmp_path / "scheduler.txt"
    scheduler.write_text(raw())
    validate = module.validate

    def changed(*args):
        result = validate(*args)
        scheduler.write_text(raw(NumCPUs="32"))
        return result

    monkeypatch.setattr(module, "validate", changed)
    with pytest.raises(ValueError):
        module.audit(PLAN, 0, scheduler, 12345, "running")
