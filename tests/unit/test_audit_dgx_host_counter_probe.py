import json
from pathlib import Path

import pytest

from benchmark_tools.audit_dgx_host_counter_probe import validate


def fixture():
    path = Path(__file__).resolve().parents[2] / "benchmark_tools/results/dgx_slurm_host_counter_probe_21799.json"
    data = json.loads(path.read_text())
    scheduler = "JobId=21799 JobState=COMPLETED ExitCode=0:0 NodeList=spark-7ff0 Partition=spark CPUs/Task=1 NumCPUs=20 MinMemoryNode=1G OverSubscribe=NO Restarts=0"
    return data, scheduler


def test_retained_scheduled_probe():
    data, scheduler = fixture()
    fields, scope = validate(data, scheduler, 21799)
    assert fields["NumCPUs"] == "20"
    assert "/job_21799/step_batch/" in scope


@pytest.mark.parametrize("problem", ["live", "wrong_node", "duplicate", "wrong_job", "source", "errors", "scope", "missing", "replay"])
def test_corruption(problem):
    data, scheduler = fixture()
    if problem == "live":
        scheduler = scheduler.replace("COMPLETED", "RUNNING")
    elif problem == "wrong_node":
        scheduler = scheduler.replace("spark-7ff0", "bizon")
    elif problem == "duplicate":
        scheduler += " JobId=21799"
    elif problem == "wrong_job":
        scheduler = scheduler.replace("21799", "21800")
    elif problem == "source":
        data["source"]["sha256"] = "changed"
    elif problem == "errors":
        data["snapshots"][0]["errors"] = [{"type": "OSError"}]
    elif problem == "scope":
        data["snapshots"][0]["raw"]["cgroup_membership"] = "0::/user.slice/session.scope"
    elif problem == "missing":
        data["snapshots"][0]["optional"].pop("host_cpu_pressure")
    else:
        data["summary"]["accounted_host_busy_cpu_s"] += 1
    with pytest.raises(ValueError):
        validate(data, scheduler, 21799)
