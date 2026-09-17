from copy import deepcopy
import json
from pathlib import Path

import pytest

from benchmark_tools.admit_dgx_launcher_smoke import PROJECT, RECIPE_SHA, verify_run


def fixture():
    results = Path(__file__).resolve().parents[2] / "benchmark_tools/results"
    spec = json.loads((results / "dgx_launcher_smoke_spec_20260917.json").read_text())
    run, order = spec["runs"][0], spec["orders"][0]
    expected = deepcopy(run)
    expected["gnu_time"] = {"executable": "/usr/bin/time", "output": str(Path(run["measurement_directory"]).parent / "native.time.tsv")}
    prepared = {"run": expected, "status": "fresh_native_inputs_prepared"}
    measured = {"requested_cpus": 20, "requested_memory_bytes": 96 * 1024 ** 3,
                "timeout_s": 900, "interval_s": 1, "host_interval_s": 30}
    runtimes = [dict(row, records=count, scientific_execution_authorized=False,
                     status="runtime_tree_identity_matches")
                for row, count in zip(spec["runtime_manifests"], [26673, 10066])]
    runtimes.append({"path": str(PROJECT / "runtime_inventory_v1/native_launcher_recipe_v1.json"),
                     "sha256": RECIPE_SHA, "records": 12, "scientific_execution_authorized": False,
                     "status": "runtime_tree_identity_matches"})
    observed = {"runtime": runtimes, "original_inputs": order["inputs_in_native_order"],
                "native_order": order["native_order"]}
    verified = {"status": "command_exited_zero", "measurement": measured,
                "before": observed, "after": deepcopy(observed),
                "source_sha256": "36243fcdfe6f49e73dbe1dd1b627bf4a330f19f3cd52dcf94a9618ddf0415cf8"}
    return run, order, prepared, verified, measured


def test_exact_verification():
    args = fixture()
    assert verify_run(*args) == args[2]["run"]


@pytest.mark.parametrize("problem", ["command", "order", "recipe", "allocation", "status"])
def test_mutations_rejected(problem):
    run, order, prepared, verified, measured = fixture()
    if problem == "command":
        prepared["run"]["native_argv"].append("--changed")
    elif problem == "order":
        verified["after"]["native_order"] = []
    elif problem == "recipe":
        verified["after"]["runtime"].pop()
    elif problem == "allocation":
        measured["requested_cpus"] = 4
    else:
        verified["status"] = "command_failed"
    with pytest.raises(ValueError):
        verify_run(run, order, prepared, verified, measured)


def test_smoke_keeps_all_scientific_nonpath_arguments():
    results = Path(__file__).resolve().parents[2] / "benchmark_tools/results"
    smoke = json.loads((results / "dgx_launcher_smoke_spec_20260917.json").read_text())
    plan = json.loads((results / "dgx_scaling_commands_20260917.json").read_text())
    assert not smoke["scientific_timing_runs_authorized"]
    for before, after in zip(plan["runs"][:3], smoke["runs"]):
        assert before["native_method"] == after["native_method"]
        assert len(before["native_argv"]) == len(after["native_argv"])
        for old, new in zip(before["native_argv"], after["native_argv"]):
            if not old.startswith("/"):
                assert old == new
            else:
                expected = old.replace(before["dataset"]["input_directory"], after["dataset"]["input_directory"])
                expected = expected.replace("/scaling_native_v1/", "/launcher_smoke_v1/")
                assert new == expected
