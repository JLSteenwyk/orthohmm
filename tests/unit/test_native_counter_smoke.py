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
