import json
from pathlib import Path

import pytest

from benchmark_tools import prepare_root_context_native as module

FOLDER = Path(module.__file__).parent / "results"


def build():
    return module.build(FOLDER / "dgx_lineage_native_plan_20260919.json", FOLDER / "ROOT_CONTEXT_NATIVE_PROTOCOL_20260919.md")


def test_only_output_paths_and_diagnostic_metadata_change():
    plan = build()
    parent = json.loads((FOLDER / "dgx_lineage_native_plan_20260919.json").read_text())
    for index, (task, original) in enumerate(zip(plan["runs"], parent["runs"])):
        restored = module.relocate(task, str(module.ROOT / "root_context_native_v1" / f"run_{index:02d}"),
                                   str(module.ROOT / "lineage_native_v1" / f"run_{index:02d}"))
        assert restored == original
    for key in ("runtime_manifests", "enumerator", "environment_paths", "environment_overrides",
                "unset_environment", "order", "launcher_python", "core_commit", "native_timeout_s", "interval_s"):
        assert plan[key] == parent[key]
    assert plan["allocation"]["time_limit_s"] == 3600
    assert plan["failure_policy"] == "stop_after_failure_retain_unrun"
    assert plan["execution_authorized"] is False and plan["scientific_timings_admitted"] is False


def test_retained_plan_reproduces():
    assert build() == json.loads((FOLDER / "dgx_root_context_native_plan_20260919.json").read_text())


@pytest.mark.parametrize("fault", ["parent", "protocol"])
def test_modified_source_rejected(tmp_path, fault):
    parent = FOLDER / "dgx_lineage_native_plan_20260919.json"
    protocol = FOLDER / "ROOT_CONTEXT_NATIVE_PROTOCOL_20260919.md"
    changed = tmp_path / "changed"
    changed.write_bytes((parent if fault == "parent" else protocol).read_bytes() + b"\n")
    with pytest.raises(ValueError):
        module.build(changed if fault == "parent" else parent, changed if fault == "protocol" else protocol)
