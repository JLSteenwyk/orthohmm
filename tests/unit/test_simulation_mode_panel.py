from copy import deepcopy
from types import SimpleNamespace

import pytest

from benchmark_tools.prepare_simulation_mode_panel import panel_rows, CONDITIONS
from benchmark_tools.run_simulation_tree_mode_control import METHODS, parser_checks, run


def fixture():
    datasets = [{"condition": c, "seed": s, "label": f"{c}_{s}"}
                for c in CONDITIONS for s in range(20261101, 20261111)]
    records = [{**d, "method": m, "status": "complete"} for d in datasets for m in METHODS]
    return {"datasets": datasets}, {"records": records}


def test_full_inventory_reuses_pilot_without_omitting_failed_tool():
    manifest, evidence = fixture()
    failed = next(r for r in evidence["records"] if r["condition"] == "divergent" and r["method"] == METHODS[1])
    failed.update(status="failed", failure_stage="native_output", reason="nonfinite graph")
    rows = panel_rows(manifest, evidence)
    assert len(rows) == 70 and sum(r["reuse_pilot"] for r in rows) == 1
    row = next(r for r in rows if r["label"] == failed["label"])
    assert row["methods"] == [METHODS[0]]
    assert row["unavailable"][METHODS[1]]["reason"] == "nonfinite graph"
    assert sum(len(r["methods"]) + len(r["unavailable"]) for r in rows) == 140


@pytest.mark.parametrize("problem", ["missing_dataset", "duplicate_dataset", "missing_result", "duplicate_result", "running", "pilot_failed"])
def test_invalid_baseline_inventory_rejected(problem):
    manifest, evidence = fixture()
    if problem == "missing_dataset":
        manifest["datasets"].pop()
    elif problem == "duplicate_dataset":
        manifest["datasets"][-1] = deepcopy(manifest["datasets"][0])
    elif problem == "missing_result":
        evidence["records"].pop()
    elif problem == "duplicate_result":
        evidence["records"].append(deepcopy(evidence["records"][0]))
    elif problem == "running":
        evidence["records"][-1]["status"] = "running"
    else:
        evidence["records"][0]["status"] = "failed"
    with pytest.raises(ValueError):
        panel_rows(manifest, evidence)


@pytest.mark.parametrize("method", METHODS)
def test_single_method_native_parser_does_not_require_other_baseline(monkeypatch, method):
    calls = []
    def invoke(argv, **kwargs):
        calls.append(argv)
        return SimpleNamespace(stdout="accepted\n", stderr="", returncode=0)
    monkeypatch.setattr("benchmark_tools.run_simulation_tree_mode_control.subprocess.run", invoke)
    manifest = {"core_root": "/frozen", "tool_entrypoints": {
        "orthohmm_python": {"absolute_path": "/oh/python"}, "orthofinder": {"absolute_path": "/of/orthofinder"}}}
    baseline = {method: {"tree": {"path": "/tree"}}}
    assert len(parser_checks(manifest, baseline, ["a", "b"], {})) == 1
    assert len(calls) == 1 and calls[0][0] == ("/oh/python" if method == METHODS[0] else "/of/python")


@pytest.mark.parametrize("methods", [(), ("unknown",), METHODS[::-1], (METHODS[0], METHODS[0])])
def test_invalid_method_subset_rejected_before_loading_manifests(tmp_path, monkeypatch, methods):
    monkeypatch.chdir(tmp_path)
    with pytest.raises(ValueError, match="canonical subset"):
        run(tmp_path, "baseline_20261101", tmp_path / "new", methods)
