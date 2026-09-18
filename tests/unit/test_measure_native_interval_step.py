from copy import deepcopy
import json
from pathlib import Path
import sys

import pytest

from benchmark_tools.measure_native_interval_step import evaluate, run_command
from benchmark_tools.run_dgx_interval_native_smoke import ROOT, relocate


@pytest.fixture
def evidence():
    path = Path(__file__).resolve().parents[2] / "benchmark_tools/results/dgx_interval_controls_21806.json"
    trial = json.loads(path.read_text())["trials"][0]
    points = trial["points"]
    before, after = deepcopy(points[0]["host"][1]), deepcopy(points[-1]["host"][0])
    for sample, point, shift in ((before, points[0], 1000), (after, points[-1], -1000)):
        sample["started_monotonic_ns"] += shift
        sample["finished_monotonic_ns"] += shift
        sample["raw"]["cgroup_membership"] = point["native_membership"]
        sample["optional"]["cgroup_cpu.stat"] = point["native_cpu_stat"]
    done = dict(exit_code=0, timed_out=False, snapshots=[before, after],
        started_ns=before["finished_monotonic_ns"] + 1000,
        finished_ns=after["started_monotonic_ns"] - 1000)
    return points, done


def test_complete_boundaries_and_interval_replay(evidence):
    points, done = evidence
    result = evaluate(points, done, 21806)
    assert len(result["intervals"]) == 10
    assert result["whole_command_screen"]["screen_passed"] is True
    assert result["scientific_timings_admitted"] is False


@pytest.mark.parametrize("failure", ["missing", "start", "finish", "inner", "gap"])
def test_boundary_or_missingness_rejection(evidence, failure):
    points, done = evidence
    if failure == "missing":
        points = points[:1]
    elif failure == "start":
        done["started_ns"] = points[0]["host"][0]["started_monotonic_ns"]
    elif failure == "finish":
        done["finished_ns"] = points[-1]["host"][1]["finished_monotonic_ns"]
    elif failure == "inner":
        done["snapshots"][0]["finished_monotonic_ns"] = done["started_ns"] + 1
    else:
        del points[4]
    with pytest.raises(ValueError):
        evaluate(points, done, 21806)


@pytest.mark.parametrize("code", [0, 7])
def test_native_exit_preserved(tmp_path, code):
    with (tmp_path / "log").open("w") as log:
        assert run_command([sys.executable, "-c", f"raise SystemExit({code})"], log, 5) == (code, False)


def test_timeout_records_failure(tmp_path):
    with (tmp_path / "log").open("w") as log:
        assert run_command([sys.executable, "-c", "import time; time.sleep(60)"], log, .05) == (124, True)


def test_only_output_prefix_is_relocated():
    value = dict(path=str(ROOT / "launcher_smoke_v1/run_00/output"),
                 command=[str(ROOT / "envs/orthohmm/bin/python"), "--cpu", "20"],
                 other=str(ROOT / "launcher_smoke_v1_other"))
    result = relocate(value)
    assert result["path"] == str(ROOT / "interval_native_smoke_v2/run_00/output")
    assert result["command"] == value["command"]
    assert result["other"] == value["other"]


def test_packaged_worker_imports_recipe_not_working_checkout(tmp_path):
    import shutil
    import subprocess

    root = Path(__file__).resolve().parents[2]
    package = tmp_path / "benchmark_tools"
    package.mkdir()
    for name in ("__init__.py", "measure_native_interval_step.py", "measure_native_counter_step.py",
                 "probe_host_counters.py", "probe_dgx_step_separation.py", "probe_interval_cpu.py",
                 "screen_bracketed_cpu.py"):
        shutil.copyfile(root / "benchmark_tools" / name, package / name)
    code = ("import runpy, pathlib; "
            f"runpy.run_path({str(package / 'measure_native_interval_step.py')!r}, run_name='import_probe'); "
            "import benchmark_tools, benchmark_tools.measure_native_counter_step as dependency; "
            f"assert pathlib.Path(benchmark_tools.__file__).parent == pathlib.Path({str(package)!r}); "
            f"assert pathlib.Path(dependency.__file__).parent == pathlib.Path({str(package)!r})")
    completed = subprocess.run([sys.executable, "-B", "-c", code], cwd=root, capture_output=True, text=True)
    assert completed.returncode == 0, completed.stderr
