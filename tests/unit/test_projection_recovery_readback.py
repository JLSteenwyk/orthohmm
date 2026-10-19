from copy import deepcopy
import json
from pathlib import Path

import pytest
import numpy as np

from benchmark_tools import report_projection_recovery as reader
from benchmark_tools import report_profile_evidence_pilot as base
from benchmark_tools.prepare_ob_candidate_neighborhood import record
from tests.unit.test_profile_evidence_readback import arm_fixture, completed, family_cells


def cells_and_rows():
    original = family_cells()
    cells = [deepcopy(original[0]), *[deepcopy(original[1]) for _ in range(4)]]
    rows = [{"label": label} for label in reader.ARMS]
    for i, cell in enumerate(cells):
        label = reader.ARMS[i]
        family = cell["details"]["family_details"]
        cell.update(label=label, dataset="OrthoBench", status="independent_content_readback_verified",
            scores={base.ENDPOINTS[0]: family["f_score"] / 100})
        cell["details"].update(precision=family["precision"] / 100, recall=family["recall"] / 100)
        identity = "legacy" if i == 0 else "projection"
        pin = {"path": "/" + label, "bytes": 10, "sha256": identity}
        cell.update(products={"product.tsv": dict(pin)}, rooted_species_tree=dict(pin),
            graph_stages=[{"label": "profiles", "bytes": 10, "sha256": identity}],
            candidate={"candidate_partition": dict(pin), "membership_constraints": dict(pin), "resources": {"wall_s": i}},
            coverage={"genes": 251378, "species": 12, "groups": 30 if i else 31, "pairs": 100 if i else 99},
            activation={"arm": label, "projection_policy": "legacy" if i == 0 else "residue_preserving",
                "candidate_digest": "same", "completely_measured_queries": 5, "applicability": "completely_measured"})
    rows[1]["advancement"] = reader.advancement(reader.point(cells[0]), reader.point(cells[1]))
    rows[3]["measurement_control_matches"] = True
    for i in (2, 4):
        rows[i].update(advancement=base.selective_gate(reader.point(cells[0]), cells[i]),
            matched_projection_advancement=reader.advancement(reader.point(cells[1]), reader.point(cells[i])),
            accuracy_interpretation="completely_measured")
    return cells, rows


def test_repair_only_gain_cannot_be_reported_as_a_backend_or_decision_gain():
    cells, rows = cells_and_rows()
    contrast = reader.verify_panel(cells, rows)
    assert contrast["f1"] == contrast["precision"] == contrast["recall"] == 0
    assert all(cells[i]["advancement"]["eligible"] for i in (1, 2, 4))
    assert all(not cells[i]["matched_projection_advancement"]["eligible"] for i in (2, 4))
    assert "advancement" not in cells[0] and "advancement" not in cells[3]


@pytest.mark.parametrize("kind", ["missing", "order", "digest", "graph", "native", "tree", "candidate",
    "constraints", "scores", "families", "coverage", "policy", "status", "dataset", "control_flag", "gate",
    "matched_gate", "applicability", "control_advancement"])
def test_panel_rejects_changed_or_selectively_reported_evidence(kind):
    cells, rows = cells_and_rows()
    if kind == "missing": cells.pop()
    elif kind == "order": cells.reverse()
    elif kind == "digest": cells[4]["activation"]["candidate_digest"] = "changed"
    elif kind == "graph": cells[3]["graph_stages"][0]["sha256"] = "changed"
    elif kind == "native": cells[3]["products"]["product.tsv"]["sha256"] = "changed"
    elif kind == "tree": cells[3]["rooted_species_tree"]["sha256"] = "changed"
    elif kind in ("candidate", "constraints"):
        cells[3]["candidate"]["candidate_partition" if kind == "candidate" else "membership_constraints"]["sha256"] = "changed"
    elif kind == "scores": cells[3]["scores"][base.ENDPOINTS[0]] += .001
    elif kind == "families": cells[3]["details"]["family_details"]["refog_records"][0]["genes"] += 1
    elif kind == "coverage": cells[3]["coverage"]["pairs"] += 1
    elif kind == "policy": cells[4]["activation"]["projection_policy"] = "legacy"
    elif kind == "status": cells[4]["status"] = "running"
    elif kind == "dataset": cells[4]["dataset"] = "QfO"
    elif kind == "control_flag": rows[3]["measurement_control_matches"] = False
    elif kind == "gate": rows[4]["advancement"]["eligible"] = False
    elif kind == "matched_gate": rows[4]["matched_projection_advancement"]["eligible"] = True
    elif kind == "applicability": rows[4]["accuracy_interpretation"] = "no_complete_measurements"
    else: rows[3]["advancement"] = {"eligible": True}
    with pytest.raises(ValueError): reader.verify_panel(cells, rows)


def test_empty_active_backend_requires_applicability_fallback_not_accuracy_advancement():
    cells, rows = cells_and_rows()
    cells[4]["activation"].update(completely_measured_queries=0, applicability="no_complete_measurements")
    rows[4].update(advancement=base.selective_gate(reader.point(cells[0]), cells[4]),
                   accuracy_interpretation="no_complete_measurements")
    reader.verify_panel(cells, rows)
    assert not cells[4]["advancement"]["eligible"]
    assert "applicability fallback" in reader.next_actions(cells)[2]["action"]


def test_real_paired_uncertainty_reports_all_24_prespecified_metric_contrasts():
    cells, _ = cells_and_rows()
    result = reader.uncertainty(cells)
    for name, count in [("versus_legacy", 4), ("versus_projection", 3), ("hmmer_versus_builtin", 1)]:
        assert len(result[name]["comparisons"]) == count
        assert result[name]["replicates"] == 100000 and result[name]["seed"] == 20261018
        assert "24 planned endpoints" in result[name]["multiplicity"]
    assert result["hmmer_versus_builtin"]["comparisons"][reader.ARMS[4]]["metrics"]["f_score"]["difference_percentage_points"] == 0


@pytest.mark.parametrize("revision", [None, "HEAD", "d3fdca0d", "A" * 40, "f" * 39, "f" * 41])
def test_launch_requires_a_full_immutable_commit_before_git_lookup(tmp_path, monkeypatch, revision):
    monkeypatch.setattr(reader.subprocess, "check_output", lambda *a, **k: pytest.fail("No Git lookup for mutable revision"))
    with pytest.raises(ValueError, match="immutable full"):
        reader.launch(tmp_path, revision)


@pytest.mark.parametrize("change", [None, "unsubmitted", "job", "arms", "cpu", "memory", "wall", "blob",
                                    "reader_revision", "reader_absolute", "reader_escape"])
def test_original_launch_snapshot_and_resources_are_bound(tmp_path, monkeypatch, change):
    revision = "f" * 40
    entry = {"execution_plan": {"path": "plan.json", "bytes": 3, "sha256": "fixed", "input_records": 42},
        "arms": list(reader.ARMS), "inference": {"submitted": True, "job_id": "24513", "planned_cpus": 16,
        "planned_memory": "128G", "planned_wall_limit": "1-00:00:00"},
        "independent_reader": {"runner_commit": "b" * 40, "scientific_reader_source": "reader"}}
    if change == "unsubmitted": entry["inference"]["submitted"] = False
    elif change == "job": entry["inference"]["job_id"] = None
    elif change == "arms": entry["arms"].reverse()
    elif change == "cpu": entry["inference"]["planned_cpus"] = 32
    elif change == "memory": entry["inference"]["planned_memory"] = "64G"
    elif change == "wall": entry["inference"]["planned_wall_limit"] = "2-00:00:00"
    elif change == "reader_revision": entry["independent_reader"]["runner_commit"] = "HEAD"
    elif change == "reader_absolute": entry["independent_reader"]["scientific_reader_source"] = "/reader"
    elif change == "reader_escape": entry["independent_reader"]["scientific_reader_source"] = "../reader"
    calls = []

    def git(command, **kwargs):
        calls.append(command)
        return ("blob" if change == "blob" else "commit") if command[1] == "cat-file" else json.dumps({"projection_recovery_pilot": entry}).encode()

    monkeypatch.setattr(reader.subprocess, "check_output", git)
    if change:
        with pytest.raises(ValueError): reader.launch(tmp_path, revision)
    else:
        _, pin, binding = reader.launch(tmp_path, revision)
        assert pin == {"path": str(tmp_path / "plan.json"), "bytes": 3, "sha256": "fixed"}
        assert binding["revision"] == revision
        assert calls[-1] == ["git", "show", revision + ":" + reader.LEDGER]


@pytest.mark.parametrize("kind", [None, "partial", "order", "source", "failure", "job"])
def test_original_complete_successful_recovery_producer_is_required(kind):
    report, plan, arms = completed("projection_recovery")
    if kind == "partial": report["status"] = "running"
    elif kind == "order": report["arms"].reverse()
    elif kind == "source": report["source_revision"] = "other"
    elif kind == "failure": report["arms"][3]["status"] = "failed"
    elif kind == "job": report["job_id"] = "other"
    if kind:
        with pytest.raises(ValueError): base.require_completed(report, plan, arms, "projection_recovery", "24477")
    else:
        base.require_completed(report, plan, arms, "projection_recovery", "24477")


@pytest.mark.parametrize("kind", [None, "wrong_tree", "no_reference", "wrong_reference_score"])
def test_measurement_admission_uses_original_repaired_tree_and_scores(tmp_path, monkeypatch, kind):
    activation = {"arm": "projection_hmmer_score_only", "projection_policy": "residue_preserving"}
    plan, row, native, details, folder = arm_fixture(tmp_path, monkeypatch,
        "projection_hmmer_score_only", experiment="projection_recovery", activation=activation)
    monkeypatch.setattr(base, "validate_activation", lambda *a, **k: activation)
    tree = folder / "native/orthohmm_phylogeny/species_tree.rooted.nwk"
    tree.write_text("(a,b);\n")
    row["rooted_species_tree"] = record(tree)
    row["native_pair_turnover"] = {"added_pairs": 0, "files": []}
    monkeypatch.setattr(base, "native_pair_turnover", lambda *a: row["native_pair_turnover"])
    reference = {"label": "projection_off", "candidate": row["candidate"],
        "products": plan["case"]["golden"], "rooted_species_tree": record(tree),
        "scores": {base.ENDPOINTS[0]: .7}, "details": details}
    row["golden_matches"]["species_tree.rooted.nwk"] = True
    plan["case"]["expected"]["f1"] = .69  # Must not be used for repaired measurement.
    if kind == "wrong_tree": reference["rooted_species_tree"]["sha256"] = "wrong"
    elif kind == "no_reference": reference = None
    elif kind == "wrong_reference_score": reference["scores"][base.ENDPOINTS[0]] = .68
    if kind:
        with pytest.raises(ValueError):
            base.admit_arm(tmp_path, plan, row, "projection_recovery", ["a", "b"], projection_control=reference)
    else:
        cell = base.admit_arm(tmp_path, plan, row, "projection_recovery", ["a", "b"], projection_control=reference)
        assert cell["rooted_species_tree"] == record(tree) and cell["graph_stages"]


def test_existing_readback_is_preserved_before_launch_lookup(tmp_path):
    with pytest.raises(FileExistsError, match="Preserve"):
        reader.run(tmp_path, "f" * 40, tmp_path)


@pytest.mark.parametrize("changed", [False, True])
def test_reader_inventory_requires_actual_committed_bytes(tmp_path, monkeypatch, changed):
    folder = tmp_path / "benchmark_tools"
    folder.mkdir()
    files = [folder / "fixture.py", *[folder / name for name in (
        reader.PROTOCOL, base.PROTOCOL, "hmm_projection_integrity_recovery_20261018.json")]]
    for path in files: path.write_text("fixture\n")
    monkeypatch.setattr(reader, "sources", lambda _: [record(files[0])])

    def git(command, **kwargs):
        if command[1] == "rev-parse": return "f" * 40
        return b"old\n" if changed else b"fixture\n"

    monkeypatch.setattr(reader.subprocess, "check_output", git)
    if changed:
        with pytest.raises(ValueError, match="committed source"):
            reader.reporter_inventory(tmp_path)
    else:
        revision, pins = reader.reporter_inventory(tmp_path)
        assert revision == "f" * 40 and len(pins) == 4


@pytest.mark.parametrize("kind", [None, "duplicate", "empty", "wrong_species", "wrong_report"])
def test_complete_fasta_checkpoint_ownership_is_repeated(tmp_path, monkeypatch, kind):
    from orthohmm.accuracy import write_accuracy_checkpoint
    fasta = tmp_path / "a.fa"
    fasta.write_text(">a\nACDE\n>b\nACDE\n>c\nACDE\n")
    species = [0, 0, 0] if kind != "wrong_species" else [0, 0, 1]
    checkpoint = write_accuracy_checkpoint(str(tmp_path / "cache"), ["a", "b", "c"], species, [0], [1], [10.])
    actual = json.loads(json.dumps(reader.validate_ownership(["a", "b", "c"], np.zeros(3, dtype=np.int32),
                                                            dict.fromkeys(["a", "b", "c"], "a.fa"))))
    if kind == "duplicate": fasta.write_text(">a\nACDE\n>a\nACDE\n>c\nACDE\n")
    elif kind == "empty": fasta.write_text(">a\nACDE\n>b\n>c\nACDE\n")
    elif kind == "wrong_report": actual["genes"] = 4
    plan = {"case": {"checkpoint": str(checkpoint), "fastas": [record(fasta)]}}
    monkeypatch.setattr(reader, "UNIVERSE", (3, 1))
    if kind:
        with pytest.raises(ValueError): reader.ownership(plan, {"ownership": actual})
    else:
        assert reader.ownership(plan, {"ownership": actual}) == ["a", "b", "c"]


@pytest.mark.parametrize("kind", [None, "incomplete_producer", "failed_arm", "changed_runtime", "wrong_reader"])
def test_full_reader_orchestration_preserves_failure_and_only_admits_complete_success(tmp_path, monkeypatch, kind):
    repo = tmp_path
    source = repo / "science"
    (source / "benchmark_tools").mkdir(parents=True)
    slurm = source / "benchmark_tools/run_profile_weighting_pilot.slurm"
    slurm.write_text("fixture\n")
    root = repo / "producer"
    root.mkdir()
    cells, rows = cells_and_rows()
    for cell in cells: cell["checked"] = []
    for row in rows: row.update(dataset="OrthoBench", status="independent_content_readback_verified")
    plan = {"source_root": str(source), "source_revision": "a" * 40, "repo": str(repo),
        "output": str(root), "cpu": 16, "experiment": "projection_recovery", "arms": list(reader.ARMS),
        "sources": [record(slurm)], "inputs": [], "runtime": {"files": []}}
    path = root / "execution_plan.json"
    path.write_text(json.dumps(plan))
    pin = record(path)
    entry = {"scientific_source": "science", "runner_commit": plan["source_revision"], "inference": {"job_id": "7"},
             "independent_reader": {"runner_commit": "c" * 40 if kind == "wrong_reader" else "b" * 40,
                 "scientific_reader_source": str(Path(reader.__file__).resolve().parents[1])}}
    report = {"execution_plan": pin, "status": "running" if kind == "incomplete_producer" else "finished_all_arms_verified",
        "source_revision": plan["source_revision"], "experiment": "projection_recovery", "job_id": "7",
        "arms": rows, "defaults_promoted": False, "goal_completed": False}
    (root / "results.json").write_text(json.dumps(report))
    monkeypatch.setattr(reader, "launch", lambda *a: (entry, pin, {"revision": "b" * 40}))
    monkeypatch.setattr(reader, "sources", lambda _: [])
    monkeypatch.setattr(base, "terminal", lambda _: {"accounting": {"State": "COMPLETED"}})
    monkeypatch.setattr(reader.subprocess, "check_output", lambda *a, **k: plan["source_revision"])
    monkeypatch.setattr(reader, "ownership", lambda *a: ["a", "b", "c"])
    monkeypatch.setattr(reader, "reporter_inventory", lambda _: ("b" * 40, []))
    monkeypatch.setattr(base, "retained_context", lambda _: ([], []))
    monkeypatch.setattr(reader, "uncertainty", lambda _: {"scope": "orchestration fixture; actual bootstrap tested separately"})
    calls = []

    def admit(repo, plan, row, experiment, names, **options):
        calls.append(row["label"])
        if kind == "failed_arm" and row["label"] == reader.ARMS[3]:
            raise ValueError("Retained native/measurement readback failure")
        if row["label"] == reader.ARMS[3]:
            assert options["projection_control"]["label"] == "projection_off"
        return deepcopy(cells[reader.ARMS.index(row["label"])])

    monkeypatch.setattr(base, "admit_arm", admit)
    runtime_calls = []

    def runtime():
        runtime_calls.append(True)
        return {"changed": True} if kind == "changed_runtime" and len(runtime_calls) > 1 else plan["runtime"]

    monkeypatch.setattr(reader, "runtime", runtime)
    output = repo / "readback"
    if kind:
        with pytest.raises(ValueError): reader.run(repo, "b" * 40, output)
        if kind in ("incomplete_producer", "wrong_reader"):
            assert not output.exists() and calls == []
        else:
            result = json.loads((output / "summary.json").read_text())
            assert result["status"] == "failed" and not result["accuracy_result_admitted"]
            with pytest.raises(FileExistsError): reader.run(repo, "b" * 40, output)
    else:
        result = reader.run(repo, "b" * 40, output)
        assert result["status"] == "projection_recovery_independently_read_back" and result["accuracy_result_admitted"]
        assert calls == list(reader.ARMS) and len(result["cells"]) == 5
        assert result["advancing_labels"] == [reader.ARMS[i] for i in (1, 2, 4)]
        assert not result["defaults_promoted"] and not result["goal_completed"]
        assert len(result["required_next_actions"]) == 4
        for checked in result["checked"]: reader.check(checked)
