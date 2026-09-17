import json
from pathlib import Path

import pytest

from benchmark_tools import authorize_dgx_scaling as module

RESULTS = Path(__file__).resolve().parents[2] / "benchmark_tools/results"


def test_authorization_preserves_all27_runs():
    spec = module.authorize(RESULTS)
    plan = json.loads((RESULTS / "dgx_scaling_commands_20260917.json").read_text())
    assert spec["runs"] == plan["runs"]
    assert spec["purpose"] == "scientific_scaling" and spec["execution_authorized"]
    assert spec["native_timeout_s"] == 85800
    assert [row["proteomes"] for row in spec["orders"]] == [4, 8, 12]


@pytest.mark.parametrize("problem", ["smoke", "overhead", "environment", "input"])
def test_gate_or_identity_changes_refuse_authorization(monkeypatch, problem):
    original = module.read_pinned
    def read(path, sha):
        data = original(path, sha)
        if problem == "smoke" and path.name == "dgx_native_launcher_smokes_admitted_20260917.json":
            data["status"] = "failed"
        elif problem == "overhead" and path.name == "dgx_verified_overhead_20260917.json":
            data["all_protocol_gates_met"] = False
        elif problem == "environment" and path.name == "dgx_launcher_smoke_spec_20260917.json":
            data["environment_overrides"]["OMP_NUM_THREADS"] = "99"
        elif problem == "input" and path.name == "dgx_native_input_order_20260917.json":
            data["datasets"][0]["inputs_in_native_order"].pop()
        return data
    monkeypatch.setattr(module, "read_pinned", read)
    with pytest.raises(ValueError):
        module.authorize(RESULTS)
