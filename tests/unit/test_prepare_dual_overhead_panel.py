from pathlib import Path

from benchmark_tools import prepare_dual_overhead_panel as module


def test_actual_frozen_native_work_and_order_preserved():
    base = Path(module.__file__).parent / "results"
    parent = base / "dgx_pressure_overhead_plan_v2_20260919.json"
    original = module.read_pinned(parent, module.PARENT_SHA)
    result = module.build(parent, base / "DUAL_COLLECTOR_OVERHEAD_PROTOCOL_20260919.md")
    relocated = module.relocate(original, module.ROOT + "/pressure_frontier_overhead_v2",
                               module.ROOT + "/dual_collector_overhead_v1")
    assert result["runs"] == relocated["runs"]
    for key in ("runtime_manifests", "order", "native_timeout_s", "launcher_python"):
        assert result[key] == original[key]
    assert len(result["runs"]) == 18
    assert result["execution_authorized"] is False
    assert result["scientific_timings_admitted"] is False
