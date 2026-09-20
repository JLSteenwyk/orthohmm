from copy import deepcopy
import gzip
import json
from pathlib import Path

import pytest

from benchmark_tools import prepare_root_context_scaling as module

RESULTS = Path(__file__).resolve().parents[2] / "benchmark_tools/results"


@pytest.fixture(scope="module")
def inputs():
    names = ("dgx_scaling_commands_20260917.json", "dgx_native_input_order_20260917.json",
             "dgx_root_context_overhead_plan_20260919.json")
    data = [json.loads((RESULTS / name).read_text()) for name in names]
    with gzip.open(RESULTS / "root_context_overhead_audit_22022_20260920.json.gz", "rt") as stream:
        data.append(json.load(stream))
    return data


def test_preserves_all_original_work_and_declares_unresolved_gates(inputs):
    original, order, overhead, audit = inputs
    plan = module.assemble(*inputs)
    assert len(plan["runs"]) == 27
    for task, before in zip(plan["runs"], original["runs"]):
        restored = module.relocate(task["run"], str(module.OUTPUT_ROOT), str(module.ROOT / "scaling_native_v1"))
        assert restored == before
        for key in ("index", "method", "proteomes", "repeat"):
            assert task[key] == before[key]
        assert task["run"]["dataset"] == before["dataset"]
        assert task["run"]["cwd"] == before["cwd"]
        assert task["run"]["measurement_directory"] != before["measurement_directory"]
    assert plan["orders"] == order["datasets"]
    assert plan["collector"] == overhead["collectors"]["root_context"]
    assert plan["allocation"]["time_limit_s"] == 86400
    assert plan["native_timeout_s"] == 85800
    assert plan["waiting_session"]["remote_timeout_s"] == 86520
    assert plan["waiting_session"]["local_timeout_s"] == 86550
    assert plan["execution_authorized"] is False
    assert plan["environment_policy"]["service_change_authorized"] is False
    assert plan["scientific_timings_admitted"] is False
    assert plan["failure_policy"]["automatic_retry"] is False
    assert plan["reporting"]["partial_medians"] is False


@pytest.mark.parametrize("case", ["missing", "reorder", "duplicate", "bool_index", "bool_repeat",
    "environment", "root", "order_bytes", "audit_incomplete", "audit_issue", "budget", "admission"])
def test_rejects_identity_or_evidence_drift(inputs, case):
    original, order, overhead = deepcopy(inputs[:3])
    audit = {k:deepcopy(inputs[3][k]) for k in ("validated_tasks", "issues", "comparison", "scientific_timings_admitted")}
    if case == "missing": original["runs"].pop()
    elif case == "reorder": original["runs"].reverse()
    elif case == "duplicate": original["runs"][0]["repeat"] = 1
    elif case == "bool_index": original["runs"][0]["index"] = False
    elif case == "bool_repeat": original["runs"][0]["repeat"] = False
    elif case == "environment": original["environment_overrides"]["OMP_NUM_THREADS"] = "20"
    elif case == "root": original["output_root"] += "other"
    elif case == "order_bytes": order["datasets"][0]["inputs_in_native_order"][0]["sha256"] = "changed"
    elif case == "audit_incomplete": audit["validated_tasks"] = 17
    elif case == "audit_issue": audit["issues"] = ["missing session"]
    elif case == "budget": audit["comparison"]["engineering_budget_passed"] = None
    elif case == "admission": audit["scientific_timings_admitted"] = True
    with pytest.raises(ValueError):
        module.assemble(original, order, overhead, audit)


def test_build_reproduces_committed_plan():
    assert module.build(RESULTS) == json.loads((RESULTS / "dgx_root_context_scaling_plan_20260920.json").read_text())


def test_protocol_hash_checked_before_other_sources(tmp_path):
    (tmp_path / "DGX_SCALING_REPLACEMENT_PROTOCOL_20260920.md").write_text("changed")
    with pytest.raises(ValueError, match="protocol"):
        module.build(tmp_path)


def test_long_run_amendment_changes_only_declared_collector_metadata():
    original = module.build(RESULTS)
    updated = module.build_long_run(RESULTS)
    retained = json.loads((RESULTS / "dgx_root_context_scaling_plan_v2_20260920.json").read_text())
    assert updated == retained
    changed = {k for k in original if original[k] != updated[k]}
    assert changed == {"status", "collector", "remaining_gates"}
    assert updated["runs"] == original["runs"]
    assert updated["collector"]["module"] == "benchmark_tools.measure_scaling_root_context"
    assert updated["execution_authorized"] is False
