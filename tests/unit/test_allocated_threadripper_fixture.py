import json

import pytest

from benchmark_tools import run_allocated_threadripper_fixture as fixture


def controller(elapsed="00:00:10"):
    return (f"JobId=42 Partition=gpu NodeList=bizon NumNodes=1 NumCPUs=64 NumTasks=1 "
        f"OverSubscribe=OK MinMemoryNode=128G Requeue=0 Restarts=0 CPUs/Task=64 "
        f"Command={fixture.SCRIPT} WorkDir={fixture.ROOT} TimeLimit=00:20:00 "
        f"ExitCode=0:0 JobState=RUNNING RunTime={elapsed}")


@pytest.mark.parametrize("elapsed", ["00:00:00", "00:14:59", "00:15:00"])
def test_fixture_remaining_allowance(elapsed):
    assert fixture.budget(controller(elapsed), 42)["allocation_mode"] == "shared"


@pytest.mark.parametrize("elapsed", ["00:15:01", "00:19:59", "01:00:00", "00:99:00", "1-00:00:00", "bad"])
def test_fixture_insufficient_or_malformed_allowance(elapsed):
    with pytest.raises(ValueError):
        fixture.budget(controller(elapsed), 42)


@pytest.mark.parametrize("replacement", ["NumCPUs=192", "Requeue=1", "Restarts=1", "JobState=COMPLETED",
    "MinMemoryNode=64G", "TimeLimit=1-02:00:00", "NodeList=other"])
def test_fixture_wrong_allocation(replacement):
    key = replacement.split("=", 1)[0]
    raw = " ".join(replacement if token.startswith(key+"=") else token for token in controller().split())
    with pytest.raises(ValueError):
        fixture.budget(raw, 42)


def test_fixture_command_is_isolated_engineering_work():
    command = fixture.command()
    assert command[:5] == [str(fixture.PYTHON), "-I", "-S", "-B", "-c"]
    assert command[5] == fixture.PROGRAM
    assert "max_workers=4" in command[5] and "< 6" in command[5]
    assert "orthohmm" not in command[5]


@pytest.mark.parametrize("change", ["schema", "root", "command", "inference", "retry", "sources", "plan", "digest"])
def test_unprepared_fixture_refused_before_output(tmp_path, change):
    output = tmp_path / "unused"
    value = dict(schema="allocated_threadripper_fixture_prepared_v1", root=str(fixture.ROOT),
        command=fixture.command(), new_sources=fixture.sources(), historical_plan=fixture.record(fixture.PLAN),
        output_directory=str(output), native_inference_authorized=False, automatic_retry=False)
    if change == "schema":
        value["schema"] = "native"
    elif change == "root":
        value["root"] += "/other"
    elif change == "command":
        value["command"] = ["/bin/true"]
    elif change == "inference":
        value["native_inference_authorized"] = True
    elif change == "retry":
        value["automatic_retry"] = True
    elif change == "sources":
        value["new_sources"] = []
    elif change == "plan":
        value["historical_plan"]["sha256"] = "0"*64
    path = tmp_path / "prepared.json"
    path.write_text(json.dumps(value))
    digest = "0"*64 if change == "digest" else fixture.record(path)["sha256"]
    with pytest.raises(ValueError):
        fixture.execute(path, digest)
    assert not output.exists()
