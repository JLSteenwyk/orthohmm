import json
from pathlib import Path

import pytest

from benchmark_tools.prepare_dgx_scaling_plan import plan, remap, ROOT


@pytest.fixture
def manifests():
    root = Path(__file__).resolve().parents[2] / "benchmark_tools/results"
    return [json.loads((root / name).read_text()) for name in ("dgx_scaling_inputs_20260917.json",
            "publication_scaling_commands_20260917.json", "publication_variable_native_methods_20260916.json")]


def test_exact_panel_and_environments(manifests):
    report = plan(*manifests)
    assert not report["execution_authorized"] and not report["inference_started"]
    assert len(report["runs"]) == 27
    assert len({r["measurement_directory"] for r in report["runs"]}) == 27
    assert report["resource_plan"]["cpus"] == 20
    assert report["resource_plan"]["memory_gib"] == 96
    assert report["environment_overrides"]["PYTHONHASHSEED"] == "0"
    assert report["environment_overrides"]["PYTHONPATH"] == str(ROOT / "core_arm_v2")
    assert any("fasttree-2.2.0" in p for p in report["environment_paths"]["orthohmm"])
    assert any("fasttree-2.1.11" in p for p in report["environment_paths"]["orthofinder"])
    for run in report["runs"]:
        assert run["cwd"] == str(ROOT / "core_arm_v2")
        if run["environment_role"] == "orthohmm":
            argv = run["native_argv"]
            assert argv[argv.index("--threads_per_worker") + 1] == "4"
        else:
            assert "-og" not in run["native_argv"]


@pytest.mark.parametrize("change", ["input", "duplicate", "count", "method", "cpu", "unmapped"])
def test_reject_changes(manifests, change):
    inputs, original, baseline = manifests
    if change == "input":
        inputs["datasets"][0]["inputs"][0]["sha256"] = "0" * 64
    elif change == "duplicate":
        inputs["datasets"][0]["inputs"][1] = inputs["datasets"][0]["inputs"][0]
    elif change == "count":
        inputs["datasets"][0]["proteins"] += 1
    elif change == "method":
        original["runs"][0]["native_argv"].append("--unexpected")
    elif change == "cpu":
        argv = original["runs"][0]["native_argv"]
        argv[argv.index("-c") + 1] = "16"
    else:
        original["runs"][0]["native_argv"][0] = "/unexpected/python"
    with pytest.raises(ValueError):
        plan(inputs, original, baseline)


def test_path_mapping_uses_component_boundaries():
    assert remap("/a/core/x", {Path("/a/core"): Path("/b/core")}) == "/b/core/x"
    with pytest.raises(ValueError):
        remap("/a/core_other/x", {Path("/a/core"): Path("/b/core")})


def test_metadata_listing_order_is_not_native_enumeration(manifests):
    manifests[0]["datasets"][0]["inputs"].reverse()
    result = plan(*manifests)
    assert any("enumeration" in gate for gate in result["remaining_gates"])
