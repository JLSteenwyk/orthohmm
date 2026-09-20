from copy import deepcopy
import json
from pathlib import Path

import pytest

from benchmark_tools import prepare_root_context_overhead as module

RESULTS = Path(module.__file__).parent / "results"
PARENT = RESULTS / "dgx_root_context_native_plan_20260919.json"
PROTOCOL = RESULTS / "ROOT_CONTEXT_OVERHEAD_PROTOCOL_20260919.md"


def test_design_and_native_commands_preserved():
    plan = module.build(PARENT, PROTOCOL)
    parent = json.loads(PARENT.read_text())
    assert len(plan["runs"]) == 18
    assert [r["index"] for r in plan["runs"]] == list(range(18))
    for pair in range(9):
        tasks = plan["runs"][2*pair:2*pair+2]
        assert [r["pair"] for r in tasks] == [pair, pair]
        assert [r["arm"] for r in tasks] == (["root_context", "lineage"] if pair//3 == 1 else ["lineage", "root_context"])
        assert tasks[0]["method"] == tasks[1]["method"]
        for task in tasks:
            original = parent["runs"][task["native_parent_index"]]
            recovered = module.relocate(task["run"], str(module.OUTPUT_ROOT / f"run_{task['index']:02d}"),
                str(module.ROOT / "root_context_native_v1" / f"run_{task['native_parent_index']:02d}"))
            assert recovered == original["run"]
    assert [r["native_parent_index"] for r in plan["runs"]][::2] == [0, 1, 2, 1, 2, 0, 2, 0, 1]
    for method in module.METHODS:
        assert len({r["pair"] for r in plan["runs"] if r["method"] == method}) == 3
    for key in ("baseline_sha256", "core_commit", "enumerator", "environment_overrides", "environment_paths",
                "interval_s", "launcher_python", "native_timeout_s", "order", "runtime_manifests", "unset_environment"):
        assert plan[key] == parent[key]


def test_periodic_arms_and_allocation_bounds():
    plan = module.build(PARENT, PROTOCOL)
    assert plan["collectors"]["lineage"]["module"] == "benchmark_tools.measure_native_lineage_step"
    assert plan["collectors"]["root_context"]["module"] == "benchmark_tools.measure_native_root_context"
    assert plan["interval_s"] == 1 and plan["native_timeout_s"] == 900
    assert plan["allocation"]["time_limit_s"] > 18*plan["native_timeout_s"]
    assert plan["waiting_session"]["local_timeout_s"] > plan["waiting_session"]["remote_timeout_s"] + 10 > 18000
    assert plan["execution_authorized"] is False and plan["scientific_timings_admitted"] is False
    assert plan["engineering_budget"] == dict(every_pair_max=.1, per_method_median_max=.05, required_pairs_per_method=3)


@pytest.mark.parametrize("source", [PARENT, PROTOCOL])
def test_pinned_input_byte_drift(tmp_path, source):
    altered = tmp_path / source.name
    altered.write_bytes(source.read_bytes() + b" ")
    with pytest.raises(ValueError):
        module.build(altered if source == PARENT else PARENT, altered if source == PROTOCOL else PROTOCOL)


def test_retained_plan_reproduces():
    assert module.build(PARENT, PROTOCOL) == json.loads((RESULTS / "dgx_root_context_overhead_plan_20260919.json").read_text())
