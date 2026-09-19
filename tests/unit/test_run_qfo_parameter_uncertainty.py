import json
from pathlib import Path

import pytest

from benchmark_tools import run_qfo_parameter_uncertainty as module
from benchmark_tools.prepare_ob_candidate_neighborhood import record
from tests.unit.test_audit_qfo_parameter_swiss import fixture
from benchmark_tools.audit_qfo_parameter_swiss import assemble


def setup(tmp_path, monkeypatch, missing=()):
    entries, baseline = fixture(tmp_path)
    for index in missing:
        entries[index] = {"arm": entries[index]["arm"], "status": "not_admitted", "reason": "pending prerequisite"}
    counts = assemble(entries, baseline)
    inventory = tmp_path / "inventory.json"
    inventory.write_text("{}")
    baseline_path = tmp_path / "baseline.json"
    baseline_path.write_text("{}")
    plan = tmp_path / "plan.json"
    plan.write_text("{}")
    protocol = tmp_path / "protocol.md"
    protocol.write_text("Frozen test protocol\n")
    monkeypatch.setattr(module, "PLAN_SHA", record(plan)["sha256"])
    monkeypatch.setattr(module, "PROTOCOL_SHA", record(protocol)["sha256"])
    counts.update(source=record(Path(module.__file__).with_name("audit_qfo_parameter_swiss.py")),
                  checked_inputs=[record(inventory), record(baseline_path)])
    monkeypatch.setattr(module, "audit", lambda *args: counts)
    return (inventory, record(inventory)["sha256"], baseline_path, plan, protocol, tmp_path / "result.json"), counts


@pytest.mark.parametrize("missing", [(), (1, 2), (1, 2, 3, 4, 5, 6), (0,)])
def test_real_kernel_complete_partial_and_unestimable_reports(tmp_path, monkeypatch, missing):
    args, _ = setup(tmp_path, monkeypatch, missing)
    result = module.run(*args)
    assert result == json.loads(args[-1].read_text())
    expected = 0 if 0 in missing else 6 - len(missing)
    assert result["estimated_contrasts"] == expected
    assert result["complete_panel"] == (expected == 6)
    assert result["uncertainty_admitted"] == (expected > 0)
    assert result["publication_ready"] is False
    assert result["replicates"] == 100000 and result["seed"] == 20260925
    assert result["multiplicity_endpoints"] == 18
    assert len(result["comparisons"]) == 6
    for row in result["comparisons"]:
        if row["status"] == "not_estimable":
            assert row["metrics"] is None and row["family_differences"] is None


@pytest.mark.parametrize("change", ["source", "protocol", "plan", "inventory", "audit_source",
                                   "mutated_input", "controls", "replicates", "seed", "multiplicity"])
def test_provenance_and_controls_fail_closed(tmp_path, monkeypatch, change):
    args, counts = setup(tmp_path, monkeypatch)
    if change == "source":
        monkeypatch.setattr(module, "SOURCES", {"audit_qfo_parameter_swiss.py": "wrong"})
    elif change in ("protocol", "plan", "inventory"):
        args[{"protocol": 4, "plan": 3, "inventory": 0}[change]].write_text("changed")
    elif change == "audit_source":
        counts["source"] = {"path": "other"}
    elif change == "mutated_input":
        def audit(*unused):
            args[0].write_text("changed during calculation")
            return counts
        monkeypatch.setattr(module, "audit", audit)
    else:
        calculate = module.calculate
        def changed(report):
            result = calculate(report)
            key = {"controls": "protocol_controls_match", "replicates": "replicates",
                   "seed": "seed", "multiplicity": "multiplicity_endpoints"}[change]
            result[key] = False if change == "controls" else -1
            return result
        monkeypatch.setattr(module, "calculate", changed)
    with pytest.raises((ValueError, json.JSONDecodeError)):
        module.run(*args)
    assert not args[-1].exists()


@pytest.mark.parametrize("symlink", [False, True])
def test_existing_output_never_overwritten(tmp_path, monkeypatch, symlink):
    args, _ = setup(tmp_path, monkeypatch)
    if symlink:
        args[-1].symlink_to(tmp_path / "missing")
    else:
        args[-1].write_text("preserved")
    monkeypatch.setattr(module, "audit", lambda *args: pytest.fail("Audited before overwrite gate"))
    with pytest.raises(FileExistsError):
        module.run(*args)


def test_raw_audit_failure_cannot_reach_bootstrap(tmp_path, monkeypatch):
    args, _ = setup(tmp_path, monkeypatch)
    def failed(*unused):
        raise ValueError("Truth differs")
    monkeypatch.setattr(module, "audit", failed)
    monkeypatch.setattr(module, "calculate", lambda *args: pytest.fail("Calculated unadmitted input"))
    with pytest.raises(ValueError, match="Truth differs"):
        module.run(*args)
    assert not args[-1].exists()
