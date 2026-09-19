import hashlib
import json
from pathlib import Path

import pytest

from benchmark_tools import run_lineage_native_diagnostic as module
from benchmark_tools.prepare_ob_candidate_neighborhood import record

FOLDER = Path(module.__file__).parent / "results"


def build():
    return module.build(FOLDER / "dgx_pressure_overhead_plan_v2_20260919.json",
                        FOLDER / "LINEAGE_NATIVE_PROTOCOL_20260919.md")


def setup(tmp_path):
    plan = tmp_path / "plan.json"
    plan.write_text(json.dumps(build()))
    files = [Path(module.__file__).resolve(), Path(module.measure.__code__.co_filename).resolve(), plan,
             FOLDER / "LINEAGE_NATIVE_PROTOCOL_20260919.md",
             FOLDER / "dgx_pressure_overhead_plan_v2_20260919.json"]
    recipe = tmp_path / "recipe.json"
    recipe.write_text(json.dumps(dict(records=[dict(record(p), kind="file") for p in files])))
    return plan, recipe


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def test_frozen_plan_and_relocated_outputs():
    plan = build()
    assert json.loads((FOLDER / "dgx_lineage_native_plan_20260919.json").read_text()) == plan
    assert [r["method"] for r in plan["runs"]] == list(module.METHODS)
    base = json.loads((FOLDER / "dgx_pressure_overhead_plan_v2_20260919.json").read_text())
    for i, row in enumerate(plan["runs"]):
        assert row["run"]["dataset"] == base["runs"][module.INDICES[i]]["run"]["dataset"]
        assert "lineage_native_v1" in row["run"]["measurement_directory"]
        assert "pressure_frontier_overhead_v2" not in json.dumps(row)
    assert not plan["scientific_timings_admitted"]


@pytest.mark.parametrize("index", [0, 1, 2])
def test_exact_selection(tmp_path, index):
    plan, recipe = setup(tmp_path)
    selected, row = module.select(plan, sha(plan), recipe, sha(recipe), index)
    assert row == selected["runs"][index]


@pytest.mark.parametrize("index", [-1, 3, True, 1.0])
def test_bad_indices(tmp_path, index):
    plan, recipe = setup(tmp_path)
    with pytest.raises(ValueError):
        module.select(plan, sha(plan), recipe, sha(recipe), index)


@pytest.mark.parametrize("fault", ["plan_sha", "recipe_sha", "changed_plan", "missing_collector", "changed_collector"])
def test_binding_failures(tmp_path, fault):
    plan, recipe = setup(tmp_path)
    if fault == "changed_plan":
        data = json.loads(plan.read_text())
        data["native_timeout_s"] = 1000
        plan.write_text(json.dumps(data))
    elif fault in ("missing_collector", "changed_collector"):
        data = json.loads(recipe.read_text())
        if fault == "missing_collector":
            data["records"].pop(1)
        else:
            data["records"][1]["sha256"] = "changed"
        recipe.write_text(json.dumps(data))
    with pytest.raises(ValueError):
        module.select(plan, "bad" if fault == "plan_sha" else sha(plan),
                      recipe, "bad" if fault == "recipe_sha" else sha(recipe), 0)
