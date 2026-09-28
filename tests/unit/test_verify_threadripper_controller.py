import pytest

from benchmark_tools.verify_threadripper_controller import validate


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
