from copy import deepcopy

import pytest

from benchmark_tools import admit_simulation_mode_panel as audit
from benchmark_tools.prepare_simulation_methods import commands
from benchmark_tools.simulation_supplied_commands import fresh_supplied_method


def test_independent_command_reconstruction_matches_frozen_builder(tmp_path):
    tree = tmp_path / "tree"
    tree.write_text("((a,b),(c,d));")
    dataset = {"label": "x", "input": str(tmp_path / "inputs")}
    dataset["methods"] = commands(dataset, tmp_path / "old", tmp_path / "frozen", tmp_path / "python", tmp_path / "of")
    original = deepcopy(dataset)
    baseline = {m: {"tree": {"path": str(tree)}} for m in audit.METHODS}
    expected = audit.expected_configuration(dataset, baseline, tmp_path / "new", audit.METHODS)
    assert expected["methods"] == {m: fresh_supplied_method(m, dataset["methods"][m], tree, tmp_path / "new" / m)
                                   for m in audit.METHODS}
    assert dataset == original
    subset = audit.expected_configuration(dataset, baseline, tmp_path / "new", audit.METHODS[:1])
    assert tuple(subset["methods"]) == audit.METHODS[:1]


def comparison_fixture(tmp_path, monkeypatch, changed=False):
    method = audit.METHODS[0]
    before, after = tmp_path / "before", tmp_path / "after"
    for directory in (before, after):
        directory.mkdir()
        (directory / "tree").write_text("((a:1,b:1):1,(c:1,d:1):1):0;")
    dataset, configured = [{"methods": {method: {"output": str(path)}}} for path in (before, after)]
    old, new = [("a", "b")], ([("a", "c")] if changed else [("b", "a")])
    admission = {"status": "admitted", "native_validation": {}}
    monkeypatch.setattr(audit, "admit_method", lambda *a: admission)
    monkeypatch.setattr(audit, "load_predictions", lambda m, p, *a: (old if p == before else new, []))
    monkeypatch.setattr(audit, "native_tree", lambda m, p: p / "tree")
    inventory = {"x": {"path": "x", "bytes": 1, "sha256": "aa"}}
    monkeypatch.setattr(audit, "artifact_inventory", lambda *a: inventory)
    observed = {"admission": admission, "pairs": audit.compare_pairs(old, new), "rooted_tree_identical": True,
                "retained_artifacts": inventory, "artifact_comparison": audit.compare_inventory(inventory, inventory)}
    args = [method, dataset, configured, {}, {}, {}, {method: {"retained_artifacts": inventory}}, observed, {}, []]
    return args, observed


@pytest.mark.parametrize("changed", [False, True])
def test_valid_non_equivalence_is_retained_not_rejected(tmp_path, monkeypatch, changed):
    args, _ = comparison_fixture(tmp_path, monkeypatch, changed)
    result = audit.compare_method(*args)
    assert result["status"] == ("not_equivalent" if changed else "equivalent")
    assert result["pairs"]["identical"] is not changed


def test_report_comparison_disagreement_is_invalid_evidence(tmp_path, monkeypatch):
    args, observed = comparison_fixture(tmp_path, monkeypatch)
    observed["pairs"]["after_pairs"] += 1
    with pytest.raises(ValueError, match="disagrees"):
        audit.compare_method(*args)


def test_native_failure_retained_without_attempting_prediction_conversion(monkeypatch):
    failure = {"status": "failed", "failure_stage": "native_output", "reason": "incomplete pairs"}
    monkeypatch.setattr(audit, "admit_method", lambda *a: failure)
    def forbidden(*args):
        pytest.fail("Failed native output must not be converted")
    monkeypatch.setattr(audit, "load_predictions", forbidden)
    result = audit.compare_method("orthofinder_full", {}, {}, {}, {}, {}, {}, {"admission": failure}, {}, [])
    assert result == {"status": "failed", "admission": failure}


def test_partial_scheduler_panel_rejected_before_outputs_created(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    accounting = "JobID|JobIDRaw|State|ExitCode|Elapsed\n" + "\n".join(
        f"21334_{i}|{30000+i}|{'RUNNING' if i == 69 else 'COMPLETED'}|0:0|00:00:01" for i in range(70))
    monkeypatch.setattr(audit.subprocess, "check_output", lambda *a, **k: accounting)
    output = tmp_path / "new"
    with pytest.raises(ValueError, match="not uniquely terminal"):
        audit.admit(tmp_path, output)
    assert not output.exists()


def test_existing_admission_refused(tmp_path):
    with pytest.raises(FileExistsError):
        audit.admit(tmp_path, tmp_path)


@pytest.mark.parametrize("change", [None, "source", "evidence"])
def test_fresh_pilot_allows_only_identical_source_relocation(tmp_path, change):
    old, new = tmp_path / "old.py", tmp_path / "new.py"
    old.write_text("same")
    new.write_text("different" if change == "source" else "same")
    prior = {"source": audit.record(old), "methods": {"x": "equivalent"}}
    fresh = {**prior, "source": audit.record(new)}
    if change == "evidence":
        fresh["methods"] = {"x": "different"}
    if change:
        with pytest.raises(ValueError, match="Fresh pilot"):
            audit.verify_fresh_pilot(fresh, prior)
    else:
        audit.verify_fresh_pilot(fresh, prior)
