from collections import Counter
from copy import deepcopy
import json
from pathlib import Path

import pytest

from benchmark_tools import prepare_frontier_overhead_panel as module

ROOT = Path(__file__).resolve().parents[2] / "benchmark_tools/results"


@pytest.fixture
def panel():
    return module.build(ROOT / "dgx_scaling_commands_20260917.json",
                        ROOT / "dgx_scientific_execution_20260917.json")


def test_all_methods_have_complete_preordered_pairs(panel):
    rows = panel["runs"]
    assert len(rows) == 18
    assert [row["index"] for row in rows] == list(range(18))
    assert Counter((r["method"], r["mode"]) for r in rows) == Counter(
        {(method, mode): 3 for method in module.METHODS for mode in ("boundary", "periodic")})
    for i in range(0, 18, 2):
        left, right = rows[i:i+2]
        assert left["method"] == right["method"]
        assert left["pair"] == right["pair"]
        assert [left["mode"], right["mode"]] == (["periodic", "boundary"] if left["pair"] == 1 else ["boundary", "periodic"])
    assert len({r["run"]["measurement_directory"] for r in rows}) == 18
    assert panel["execution_authorized"] is False
    assert panel["scientific_timings_admitted"] is False


def test_retained_plan_reproduces_from_pinned_sources(panel):
    assert json.loads((ROOT / "dgx_frontier_overhead_plan_20260918.json").read_text()) == panel


def test_commands_inputs_and_parameters_round_trip_to_original(panel):
    original = module.read_pinned(ROOT / "dgx_scaling_commands_20260917.json", module.PLAN_SHA)
    for row in panel["runs"]:
        template = original["runs"][row["original_index"]]
        restored = module.relocate(row["run"],
            f"{module.ROOT}/frontier_overhead_v1/run_{row['index']:02d}",
            f"{module.ROOT}/scaling_native_v1/run_{row['original_index']:02d}")
        restored.update(index=template["index"], repeat=template["repeat"])
        assert restored == template
    assert panel["native_timeout_s"] == 900
    assert panel["engineering_budget"] == {"per_method_median_max": .05, "every_pair_max": .10}


def test_relocation_preserves_input_and_neighbor_prefixes():
    value = {"paths": ["/old/file", "/older/file", "/old", 2]}
    saved = deepcopy(value)
    assert module.relocate(value, "/old", "/new") == {"paths": ["/new/file", "/older/file", "/new", 2]}
    assert value == saved


def test_changed_source_rejected(tmp_path):
    path = tmp_path / "changed.json"
    path.write_text("{}")
    with pytest.raises(ValueError, match="digest differs"):
        module.build(path, ROOT / "dgx_scientific_execution_20260917.json")
