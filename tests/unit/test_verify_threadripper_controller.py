import json
from types import SimpleNamespace

import pytest

from benchmark_tools.verify_threadripper_controller import validate, remaining_budget, ReleaseBudgetGuard


RAW = ("JobId=42 JobState=RUNNING Partition=gpu NodeList=bizon NumNodes=1 "
       "NumCPUs=192 NumTasks=1 CPUs/Task=64 OverSubscribe=NO MinMemoryNode=128G "
       "Requeue=0 Restarts=0 Command=/recipe/run.sh WorkDir=/recipe "
       "TimeLimit=1-00:00:00 ExitCode=0:0")


def check(raw=RAW, phase="running"):
    return validate(raw, 42, phase, command="/recipe/run.sh", cwd="/recipe")


def test_running_allocation_separates_node_and_task_cpus():
    r = check()
    assert r["fields"]["NumCPUs"] == "192"
    assert r["fields"]["CPUs/Task"] == "64"
    assert r["scientific_execution_authorized"] is False


@pytest.mark.parametrize("old,new", [
    ("NumCPUs=192", "NumCPUs=64"), ("CPUs/Task=64", "CPUs/Task=32"),
    ("MinMemoryNode=128G", "MinMemoryNode=96G"), ("OverSubscribe=NO", "OverSubscribe=YES"),
    ("NodeList=bizon", "NodeList=other"), ("Requeue=0", "Requeue=1"),
    ("Restarts=0", "Restarts=1"), ("Command=/recipe/run.sh", "Command=/other"),
    ("WorkDir=/recipe", "WorkDir=/other"), ("TimeLimit=1-00:00:00", "TimeLimit=00:02:00"),
    ("JobId=42", "JobId=43"), ("Partition=gpu", "Partition=other")])
def test_policy_drift_rejected(old, new):
    with pytest.raises(ValueError):
        check(RAW.replace(old, new))


@pytest.mark.parametrize("suffix", [" JobId=42", " ArrayJobId=42", " HetJobId=42", "\n"+RAW])
def test_duplicate_or_composite_record_rejected(suffix):
    with pytest.raises(ValueError):
        check(RAW+suffix)


@pytest.mark.parametrize("state", ["PENDING", "COMPLETING", "COMPLETED"])
def test_running_phase_requires_running(state):
    with pytest.raises(ValueError):
        check(RAW.replace("RUNNING", state))


def test_running_record_cannot_prove_terminal():
    with pytest.raises(ValueError):
        check(phase="terminal")


def test_invalid_exit_status():
    with pytest.raises(ValueError, match="exit"):
        check(RAW.replace("ExitCode=0:0", "ExitCode=unknown"))


@pytest.mark.parametrize("state,exit_code", [("COMPLETED", "0:0"), ("FAILED", "1:0"),
                                            ("TIMEOUT", "0:15"), ("CANCELLED", "0:15")])
def test_terminal_outcome_is_preserved_not_admitted(state, exit_code):
    raw = RAW.replace("RUNNING", state).replace("ExitCode=0:0", "ExitCode="+exit_code)
    r = check(raw, "terminal")
    assert r["scheduler_terminal_verified"] is True
    assert r["scheduler_state"] == state
    assert r["scheduler_exit_code"] == exit_code
    assert r["scientific_timings_admitted"] is False


@pytest.mark.parametrize("runtime,accepted", [("00:00:00", True), ("00:59:28", True),
    ("00:59:29", False), ("01:00:00", False), ("1-00:00:00", False),
    ("UNKNOWN", False), ("00:60:00", False), ("24:00:00", False)])
def test_release_budget(runtime, accepted):
    raw = RAW.replace("1-00:00:00", "1-02:00:00") + " RunTime=" + runtime
    def invoke():
        return remaining_budget(raw, 42, command="/recipe/run.sh", cwd="/recipe", query_elapsed_s=0.2)
    if accepted:
        result = invoke()
        assert result["conservative_available_s"] >= 90000
        assert result["native_timeout_s"] == 85800
    else:
        with pytest.raises(ValueError):
            invoke()


@pytest.mark.parametrize("age", [-1, 5.01, float("nan"), float("inf"), True])
def test_release_budget_rejects_invalid_freshness(age):
    raw = RAW.replace("1-00:00:00", "1-02:00:00") + " RunTime=00:00:10"
    with pytest.raises(ValueError, match="freshness"):
        remaining_budget(raw, 42, command="/recipe/run.sh", cwd="/recipe", query_elapsed_s=age)


@pytest.mark.parametrize("failure", [None, "budget", "query", "timeout"])
def test_guard_retains_query_and_failure(tmp_path, failure):
    def runner(argv, **kwargs):
        assert argv == ["scontrol", "show", "job", "42", "--oneliner"]
        assert kwargs["timeout"] == 5
        if failure == "timeout":
            raise TimeoutError("query stalled")
        raw = RAW.replace("1-00:00:00", "1-02:00:00")
        raw += " RunTime=" + ("01:00:00" if failure == "budget" else "00:00:10")
        return SimpleNamespace(returncode=1 if failure == "query" else 0, stdout=raw, stderr="")
    times = iter([1., 1.2])
    guard = ReleaseBudgetGuard(42, command="/recipe/run.sh", cwd="/recipe", runner=runner,
                               clock=lambda: next(times))
    if failure:
        with pytest.raises((ValueError, TimeoutError)):
            guard(tmp_path)
    else:
        assert guard(tmp_path)["status"] == "threadripper_release_budget_checked"
    result = json.loads((tmp_path / "release_budget.json").read_text())
    assert result["status"] == ("release_budget_check_failed" if failure else "release_budget_check_passed")
