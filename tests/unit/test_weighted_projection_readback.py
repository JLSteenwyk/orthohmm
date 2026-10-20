from copy import deepcopy
import json
from pathlib import Path

import pytest

from benchmark_tools import report_weighted_projection as reader
from benchmark_tools import report_profile_evidence_pilot as base
from benchmark_tools import weighted_projection_experiment as experiment
from benchmark_tools.run_profile_weighting_pilot import advancement
from benchmark_tools.prepare_ob_candidate_neighborhood import record
from tests.unit.test_profile_evidence_readback import arm_fixture, completed, family_cells


def cells_and_rows():
    spec = experiment.specification()
    cells, rows = [], []
    for i, label in enumerate(experiment.LABELS):
        arm = experiment.configuration(label)
        point = spec["control_points"].get(label, {"f1": .745, "precision": .82, "recall": .68})
        active = {"arm": label, "mode": arm["profile_sequence_weighting"],
            "projection_policy": arm["profile_projection_policy"], "role": arm["role"],
            "profiles_built": 2, "joint_profiles": 2 if i >= 4 else 0}
        cells.append({"label": label, "dataset": "OrthoBench", "status": "independent_content_readback_verified",
            "scores": {base.ENDPOINTS[0]: point["f1"]}, "details": {k: point[k] for k in ("precision", "recall")},
            "activation": active, "arm_parameters": arm})
        rows.append({"label": label, "arm_parameters": arm, "graph_control_matches": True})
    controls = {c["label"]: experiment.point(c) for c in cells[:4]}
    for cell, row in zip(cells[4:], rows[4:]):
        row["advancement"] = experiment.matched_gate(cell["label"], experiment.point(cell), controls, advancement)
    return cells, rows


def test_both_passing_combinations_are_preserved_not_only_apparent_best():
    cells, rows = cells_and_rows()
    gates = reader.verify_panel(cells, rows)
    assert list(gates) == list(experiment.LABELS[4:])
    assert all(gate["eligible"] for gate in gates.values())
    assert all("advancement" not in c for c in cells[:4])


@pytest.mark.parametrize("kind", ["missing", "order", "status", "dataset", "arm", "requested", "mode", "policy",
    "role", "inactive", "joint_missing", "control_joint", "control_point", "control_promoted", "control_graph", "gate"])
def test_reader_rejects_partial_or_changed_panel_and_claimed_gates(kind):
    cells, rows = cells_and_rows()
    if kind == "missing": cells.pop()
    elif kind == "order": rows.reverse()
    elif kind == "status": cells[5]["status"] = "failed"
    elif kind == "dataset": cells[5]["dataset"] = "QfO"
    elif kind == "arm": cells[5]["activation"]["arm"] = "repair_henikoff"
    elif kind == "requested": cells[5]["arm_parameters"] = {}
    elif kind == "mode": cells[5]["activation"]["mode"] = "none"
    elif kind == "policy": cells[5]["activation"]["projection_policy"] = "legacy"
    elif kind == "role": cells[5]["activation"]["role"] = "constituent_control"
    elif kind == "inactive": cells[5]["activation"]["profiles_built"] = 0
    elif kind == "joint_missing": cells[5]["activation"]["joint_profiles"] = 0
    elif kind == "control_joint": cells[2]["activation"]["joint_profiles"] = 2
    elif kind == "control_point": cells[2]["details"]["recall"] += .001
    elif kind == "control_promoted": rows[2]["advancement"] = {"eligible": True}
    elif kind == "control_graph": rows[2]["graph_control_matches"] = False
    else: rows[5]["advancement"]["eligible"] = False
    with pytest.raises(ValueError):
        reader.verify_panel(cells, rows)


def test_beating_canonical_but_not_stronger_constituent_remains_stopped():
    cells, rows = cells_and_rows()
    cells[4]["scores"][base.ENDPOINTS[0]] = .7415
    controls = {c["label"]: experiment.point(c) for c in cells[:4]}
    rows[4]["advancement"] = experiment.matched_gate(cells[4]["label"], experiment.point(cells[4]), controls, advancement)
    gates = reader.verify_panel(cells, rows)
    assert not gates[experiment.LABELS[4]]["eligible"]
    assert gates[experiment.LABELS[5]]["eligible"]


def test_six_arm_terminal_completion_contract_is_reused_without_weakening():
    report, plan, arms = completed("weighted_projection")
    base.require_completed(report, plan, arms, "weighted_projection", "24477")
    report["arms"][-1]["status"] = "not_run_failed_constituent_control_gate"
    with pytest.raises(ValueError, match="complete unchanged"):
        base.require_completed(report, plan, arms, "weighted_projection", "24477")


def test_actual_paired_family_bootstrap_covers_all18_registered_contrasts():
    original = family_cells()
    cells = [deepcopy(original[i >= 4]) for i in range(6)]
    for label, cell in zip(experiment.LABELS, cells): cell["label"] = label
    result = reader.uncertainty(cells)
    assert result["multiplicity_endpoints"] == 18 and len(result["contrasts"]) == 6
    assert [[r["candidate"], r["control"]] for r in result["contrasts"]] == experiment.specification()["uncertainty"]["contrasts"]
    for row in result["contrasts"]:
        value = row["result"]
        assert value["replicates"] == 100000 and value["seed"] == 20261019 and len(value["families"]) == 70
        assert "18 planned endpoints" in value["multiplicity"]
        metrics = value["comparisons"][row["candidate"]]["metrics"]
        assert set(metrics) == {"f_score", "precision", "recall"}
        assert all(len(v["paired_percentile_ci"]) == len(v["bonferroni_percentile_ci"]) == 2 for v in metrics.values())


@pytest.mark.parametrize("kind", ["missing", "nonfinite", "family_count", "point"])
def test_bootstrap_cannot_use_partial_or_wrong_family_statistics(kind):
    original = family_cells()
    cells = [deepcopy(original[0]) for _ in range(6)]
    for label, cell in zip(experiment.LABELS, cells): cell["label"] = label
    if kind == "missing": cells.pop()
    elif kind == "nonfinite": cells[0]["details"]["family_details"]["refog_records"][0]["true_positive"] = float("nan")
    elif kind == "family_count": cells[0]["details"]["family_details"]["refog_records"].pop()
    else: cells[0]["details"]["family_details"]["f_score"] += .001
    with pytest.raises(ValueError): reader.uncertainty(cells)


@pytest.mark.parametrize("revision", [None, "HEAD", "abc1234", "A" * 40, "f" * 39, "f" * 41])
def test_mutable_launch_revision_rejected_before_any_lookup(tmp_path, monkeypatch, revision):
    monkeypatch.setattr(reader.subprocess, "check_output", lambda *a, **kw: pytest.fail("No lookup"))
    with pytest.raises(ValueError, match="immutable full"):
        reader.launch(tmp_path, revision)


def launch_entry():
    return {"arms": list(experiment.LABELS), "runner_commit": "a" * 40, "scientific_source": "science",
        "execution_plan": {"path": "pilot/plan.json", "bytes": 1, "sha256": "pin"},
        "inference": {"submitted": True, "job_id": "24600", "planned_cpus": 16,
                      "planned_memory": "128G", "planned_wall_limit": "1-00:00:00"},
        "independent_reader": {"runner_commit": "b" * 40, "scientific_reader_source": "reader"}}


@pytest.mark.parametrize("kind", [None, "unsubmitted", "job", "arms", "cpus", "memory", "wall", "reader", "escape", "science"])
def test_launch_binds_exact_resources_and_two_frozen_sources(tmp_path, monkeypatch, kind):
    entry = launch_entry()
    if kind == "unsubmitted": entry["inference"]["submitted"] = False
    elif kind == "job": entry["inference"]["job_id"] = "unknown"
    elif kind == "arms": entry["arms"].reverse()
    elif kind in ("cpus", "memory", "wall"):
        entry["inference"][{"cpus": "planned_cpus", "memory": "planned_memory", "wall": "planned_wall_limit"}[kind]] = "wrong"
    elif kind == "reader": entry["independent_reader"]["runner_commit"] = "HEAD"
    elif kind == "escape": entry["independent_reader"]["scientific_reader_source"] = "../escape"
    elif kind == "science": entry["scientific_source"] = "/absolute"
    payload = json.dumps({reader.ENTRY: entry}).encode()
    monkeypatch.setattr(reader.subprocess, "check_output", lambda command, **kw:
        "commit" if command[1] == "cat-file" else payload)
    if kind:
        with pytest.raises(ValueError): reader.launch(tmp_path, "f" * 40)
    else:
        observed, pin, binding = reader.launch(tmp_path, "f" * 40)
        assert observed == entry and pin["path"] == str(tmp_path / "pilot/plan.json")
        assert binding["entry_key"] == reader.ENTRY and len(binding["git_blob_sha256"]) == 64


def test_previous_readback_is_never_overwritten(tmp_path):
    with pytest.raises(FileExistsError, match="Preserve every"):
        reader.run(tmp_path, "f" * 40, tmp_path)


@pytest.mark.parametrize("label", experiment.LABELS)
def test_raw_arm_readback_revalidates_all_six_structured_roles_and_products(tmp_path, monkeypatch, label):
    plan, row, _, _, _ = arm_fixture(tmp_path, monkeypatch, label, experiment="weighted_projection")
    actual = base.admit_arm(tmp_path, plan, row, "weighted_projection", ["a", "b"])
    assert actual["arm_parameters"] == experiment.configuration(label)
    assert actual["rooted_species_tree"] == row["rooted_species_tree"]
    assert actual["graph_stages"] and actual["checked"]
    assert actual["activation"]["joint_profiles"] == int(label in experiment.LABELS[4:])


@pytest.mark.parametrize("kind", ["graph", "point", "coverage", "tree", "candidate", "role", "turnover"])
def test_noncanonical_control_raw_readback_cannot_bypass_constituent_validation(tmp_path, monkeypatch, kind):
    label = "legacy_henikoff"
    plan, row, _, _, _ = arm_fixture(tmp_path, monkeypatch, label, experiment="weighted_projection")
    reference = plan["constituent_controls"][label]
    if kind == "graph": reference["graph_stages"][0]["sha256"] = "changed"
    elif kind == "point": reference["expected"]["f1"] += .001
    elif kind == "coverage": reference["coverage"]["pairs"] += 1
    elif kind == "tree": row["rooted_species_tree"] = {**row["rooted_species_tree"], "sha256": "changed"}
    elif kind == "candidate": reference["candidate"] = {**reference["candidate"], "sha256": "changed"}
    elif kind == "role": row["arm_parameters"] = {}
    else: monkeypatch.setattr(base, "native_pair_turnover", lambda *a: {"files": [], "added_pairs": 1, "lost_pairs": 0})
    with pytest.raises(ValueError): base.admit_arm(tmp_path, plan, row, "weighted_projection", ["a", "b"])


@pytest.mark.parametrize("kind", [None, "partial", "failed_arm", "references", "runtime", "reader"])
def test_reader_orchestration_requires_complete_success_and_preserves_failed_readbacks(tmp_path, monkeypatch, kind):
    source, root = tmp_path / "science", tmp_path / "producer"
    (source / "benchmark_tools").mkdir(parents=True)
    slurm = source / "benchmark_tools/run_profile_weighting_pilot.slurm"
    slurm.write_text("fixture\n")
    root.mkdir()
    controls = {label: {"reference": label} for label in experiment.LABELS[:4]}
    plan = {"source_root": str(source), "source_revision": "a" * 40, "repo": str(tmp_path),
        "output": str(root), "cpu": 16, "experiment": "weighted_projection", "arms": list(experiment.LABELS),
        "arm_parameters": experiment.specification()["arms"], "constituent_controls": controls,
        "sources": [record(slurm)], "inputs": [], "runtime": {"files": []}}
    plan_path = root / "execution_plan.json"
    plan_path.write_text(json.dumps(plan))
    pin = record(plan_path)
    cells, rows = cells_and_rows()
    for cell in cells: cell["checked"] = []
    for row in rows: row.update(dataset="OrthoBench", status="independent_content_readback_verified")
    report = {"execution_plan": pin, "status": "running" if kind == "partial" else "finished_all_arms_verified",
        "source_revision": plan["source_revision"], "experiment": "weighted_projection", "job_id": "7",
        "arms": rows, "defaults_promoted": False, "goal_completed": False}
    (root / "results.json").write_text(json.dumps(report))
    entry = {"scientific_source": "science", "runner_commit": plan["source_revision"], "inference": {"job_id": "7"},
        "independent_reader": {"runner_commit": "c" * 40 if kind == "reader" else "b" * 40,
            "scientific_reader_source": str(Path(reader.__file__).resolve().parents[1])}}
    monkeypatch.setattr(reader, "launch", lambda *a: (entry, pin, {"revision": "f" * 40}))
    monkeypatch.setattr(reader, "sources", lambda _: [])
    monkeypatch.setattr(reader.subprocess, "check_output", lambda *a, **kw: plan["source_revision"])
    monkeypatch.setattr(base, "terminal", lambda _: {"accounting": {"State": "COMPLETED"}})
    monkeypatch.setattr(experiment, "control_references", lambda *a: ({} if kind == "references" else controls, []))
    monkeypatch.setattr(reader.recovery, "ownership", lambda *a: ["a", "b"])
    monkeypatch.setattr(reader, "reader_inventory", lambda _: ("b" * 40, []))
    monkeypatch.setattr(base, "retained_context", lambda _: ([], []))
    monkeypatch.setattr(reader, "uncertainty", lambda _: {"scope": "orchestration only; bootstrap tested separately"})
    runtime_calls = []
    def runtime():
        runtime_calls.append(True)
        return {"changed": True} if kind == "runtime" and len(runtime_calls) > 1 else plan["runtime"]
    monkeypatch.setattr(reader, "runtime", runtime)
    calls = []
    def admit(repo, plan, row, method, names):
        calls.append(row["label"])
        if kind == "failed_arm" and row["label"] == experiment.LABELS[4]:
            raise ValueError("Retained native readback failure")
        return deepcopy(cells[experiment.LABELS.index(row["label"])])
    monkeypatch.setattr(base, "admit_arm", admit)
    output = tmp_path / "readback"
    if kind:
        with pytest.raises(ValueError): reader.run(tmp_path, "f" * 40, output)
        if kind in ("partial", "references", "reader"):
            assert not output.exists() and not calls
        else:
            result = json.loads((output / "summary.json").read_text())
            assert result["status"] == "failed" and not result["accuracy_result_admitted"]
            with pytest.raises(FileExistsError): reader.run(tmp_path, "f" * 40, output)
    else:
        result = reader.run(tmp_path, "f" * 40, output)
        assert result["status"] == "weighted_projection_independently_read_back" and result["accuracy_result_admitted"]
        assert calls == list(experiment.LABELS) and result["advancing_labels"] == list(experiment.LABELS[4:])
        assert not result["defaults_promoted"] and not result["goal_completed"]


@pytest.mark.parametrize("changed", [False, True])
def test_reader_inventory_verifies_committed_bytes_including_its_launcher(tmp_path, monkeypatch, changed):
    folder = tmp_path / "benchmark_tools"
    folder.mkdir()
    paths = [folder / "reader.py", *[folder / p for p in (
        experiment.PROTOCOL, base.PROTOCOL, "report_weighted_projection.slurm")]]
    for path in paths: path.write_text("fixture\n")
    monkeypatch.setattr(reader, "sources", lambda _: [record(paths[0])])
    monkeypatch.setattr(reader.subprocess, "check_output", lambda command, **kw:
        "f" * 40 if command[1] == "rev-parse" else b"changed" if changed else b"fixture\n")
    if changed:
        with pytest.raises(ValueError, match="committed source"): reader.reader_inventory(tmp_path)
    else:
        revision, pins = reader.reader_inventory(tmp_path)
        assert revision == "f" * 40 and len(pins) == 4
