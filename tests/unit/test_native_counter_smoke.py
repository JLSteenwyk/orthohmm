import copy
import json
from pathlib import Path

import pytest

from benchmark_tools.measure_native_counter_step import validate
from benchmark_tools.run_dgx_counter_native_smoke import ROOT, SPEC_SHA, relocate, read_pinned


def test_relocated_commands_change_only_output_paths():
    path = Path(__file__).resolve().parents[2] / "benchmark_tools/results/dgx_launcher_smoke_spec_20260917.json"
    spec = read_pinned(path, SPEC_SHA)
    original = copy.deepcopy(spec)
    for row in spec["runs"]:
        moved = relocate(row)
        assert moved["dataset"] == row["dataset"]
        assert moved["cwd"] == row["cwd"]
        old = str(ROOT / "launcher_smoke_v1")
        new = str(ROOT / "counter_native_smoke_v1")
        assert json.loads(json.dumps(moved).replace(new + "/", old + "/")) == row
        assert moved["native_argv"] != row["native_argv"]
        validate(moved["native_argv"], 20, 900, 1.)
    assert spec == original


@pytest.mark.parametrize("command,cpus,timeout,interval", [
    ([],20,900,1), (["python"],20,900,1), (["/bin/true", ""],20,900,1),
    (["/bin/true"],1,900,1), (["/bin/true"],20,901,1),
    (["/bin/true"],20,0,1), (["/bin/true"],20,900,.2)])
def test_frozen_limits(command, cpus, timeout, interval):
    with pytest.raises(ValueError):
        validate(command, cpus, timeout, interval)


def test_non_output_prefix_is_unchanged():
    other = str(ROOT / "launcher_smoke_v10/input")
    assert relocate(other) == other


def test_retained_three_native_smokes():
    from benchmark_tools.probe_host_counters import summarize
    root = Path(__file__).resolve().parents[2]
    report = json.loads((root / "benchmark_tools/results/dgx_counter_native_smokes_21802.json").read_text())
    assert report["scientific_timings_admitted"] == 0
    assert report["publication_ready"] is False
    assert report["controlled_workload_verified"] is False
    assert [r["index"] for r in report["runs"]] == [0, 1, 2]
    for row in report["runs"]:
        proof = row["verification"]
        assert proof["status"] == "command_exited_zero"
        assert proof["before"] == {k: v for k, v in proof["after"].items() if k != "copied_inputs"}
        m = proof["measurement"]
        assert m["native"]["exit_code"] == 0
        assert m["counter_read_errors"] == 0
        assert summarize(*m["native_snapshots"], 100) == m["host_summary"]
        assert row["output_validation"]["input_genes"] == 645
        assert row["output_validation"]["status"] == "native_scaling_outputs_checked"
        assert row["output_validation"]["resource_measurements_admitted"] is False
        assert "JobState=COMPLETED " in row["scheduler"]
        assert "ExitCode=0:0" in row["scheduler"]
