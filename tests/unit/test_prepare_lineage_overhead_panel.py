import json
from pathlib import Path

import pytest

from benchmark_tools import prepare_lineage_overhead_panel as module

BASE = Path(module.__file__).parent / "results"
PARENT = BASE / "dgx_pressure_overhead_plan_v2_20260919.json"
PROTOCOL = BASE / "LINEAGE_COLLECTOR_OVERHEAD_PROTOCOL_20260919.md"


def test_retained_plan_reproduces_exactly():
    retained = json.loads((BASE / "dgx_lineage_overhead_plan_20260919.json").read_text())
    assert retained == module.build(PARENT, PROTOCOL)


def test_all_frozen_workloads_and_settings_preserved():
    original = module.read_pinned(PARENT, module.PARENT_SHA)
    result = module.build(PARENT, PROTOCOL)
    expected = module.relocate(original, module.ROOT + "/pressure_frontier_overhead_v2",
                               module.OUTPUT_ROOT)
    assert result["runs"] == expected["runs"]
    for key in ("runtime_manifests", "order", "native_timeout_s", "launcher_python",
                "resource_plan", "engineering_budget", "interval_s", "eligibility_delay_s",
                "environment_paths", "environment_overrides", "unset_environment"):
        assert result[key] == original[key]
    pairs = {}
    for row in result["runs"]:
        pairs.setdefault((row["method"], row["pair"]), []).append(row["mode"])
        assert row["run"]["measurement_directory"].startswith(module.OUTPUT_ROOT + "/")
    assert len(pairs) == 9
    assert all(sorted(modes) == ["boundary", "periodic"] for modes in pairs.values())
    assert result["periodic_collector"] == "benchmark_tools.measure_native_lineage_step.measure"
    assert result["boundary_collector"] == "benchmark_tools.measure_lineage_boundary_step.measure"
    assert result["native_pressure"] is True
    for key in ("execution_authorized", "scientific_timings_admitted", "publication_ready"):
        assert result[key] is False
    assert result["protocol_sha256"] == module.hashlib.sha256(PROTOCOL.read_bytes()).hexdigest()


def test_changed_parent_bytes_rejected(tmp_path):
    changed = tmp_path / "parent.json"
    changed.write_bytes(PARENT.read_bytes() + b"\n")
    with pytest.raises(ValueError):
        module.build(changed, PROTOCOL)


@pytest.mark.parametrize("indices", [list(range(17)), list(reversed(range(18))), [0] * 18])
def test_incomplete_or_reordered_parent_rejected(monkeypatch, indices):
    original = json.loads(PARENT.read_text())
    original["runs"] = [{"index": index} for index in indices]
    monkeypatch.setattr(module, "read_pinned", lambda *args: original)
    with pytest.raises(ValueError, match="complete frozen"):
        module.build(PARENT, PROTOCOL)
