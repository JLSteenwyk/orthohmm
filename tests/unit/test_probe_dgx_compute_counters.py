import copy
import hashlib
import json
from pathlib import Path

import pytest

from benchmark_tools import probe_dgx_compute_counters as module


@pytest.fixture
def trial():
    path = Path(__file__).resolve().parents[2] / "benchmark_tools/results/dgx_step_separation_probe_21800.json"
    retained = json.loads(path.read_text())["trials"][0]
    before, after = copy.deepcopy(retained["native_snapshots"])
    before["optional"]["cgroup_cpu.stat"] = "usage_usec 1000000\n"
    after["optional"]["cgroup_cpu.stat"] = "usage_usec 6000000\n"
    after["optional"]["cgroup_memory.peak"] = str(module.BUFFER_BYTES * 2)
    return dict(monitor=False, native_snapshots=[before, after], snapshots=[retained["snapshots"][0]],
        clock_ticks_per_second=100, native_exit_code=0,
        work_started_ns=before["finished_monotonic_ns"] + 1,
        work_finished_ns=after["started_monotonic_ns"] - 1,
        children=[dict(checksum=module.expected_checksum(), cpu_s=1., cgroup=before["raw"]["cgroup_membership"])
                  for _ in range(module.CHILDREN)])


def test_small_fixed_work():
    assert module.checksum(3, 10) == hashlib.sha256(hashlib.sha256(b"x" * 10).digest() * 3).hexdigest()


def test_positive_accounting(trial):
    result = module.validate_trial(trial, 21800)
    assert result["native_cgroup_cpu_s"] == 5
    assert result["child_cpu_s"] == 4


@pytest.mark.parametrize("failure", ["checksum", "scope", "cpu", "peak", "bracket", "missing", "reads", "monitor", "exit"])
def test_reject_invalid_evidence(trial, failure):
    if failure == "checksum":
        trial["children"][0]["checksum"] = "bad"
    elif failure == "scope":
        trial["children"][0]["cgroup"] = "0::/job_21800/step_batch\n"
    elif failure == "cpu":
        trial["native_snapshots"][1]["optional"]["cgroup_cpu.stat"] = "usage_usec 1000001\n"
    elif failure == "peak":
        trial["native_snapshots"][1]["optional"]["cgroup_memory.peak"] = "1"
    elif failure == "bracket":
        trial["work_started_ns"] = 0
    elif failure == "missing":
        trial["children"].pop()
    elif failure == "reads":
        trial["snapshots"][0]["errors"] = [dict(field="host_cpu_pressure")]
    elif failure == "monitor":
        trial["monitor"] = True
    else:
        trial["native_exit_code"] = 1
    with pytest.raises(ValueError):
        module.validate_trial(trial, 21800)


def test_unscheduled_rejected(tmp_path, monkeypatch):
    monkeypatch.delenv("SLURM_JOB_ID", raising=False)
    with pytest.raises(KeyError):
        module.run(tmp_path / "out")
    assert not (tmp_path / "out").exists()
