from copy import deepcopy
import csv
import json
from pathlib import Path

import numpy as np
import pytest

from benchmark_tools import report_weighting_followup as reader
from benchmark_tools.prepare_ob_candidate_neighborhood import record
from benchmark_tools.report_hillclimb import METRIC_NAMES, endpoint_vector
from orthohmm.accuracy import write_accuracy_checkpoint
from tests.unit.test_profile_evidence_readback import arm_fixture, timing


def fixture(tmp_path, monkeypatch, dataset="ThreeKingdoms", mode="henikoff"):
    old, row, native, _, folder = arm_fixture(tmp_path, monkeypatch, mode)
    case = old["case"]
    scores = {metric: .7 for metric in METRIC_NAMES[dataset].values()}
    vector = endpoint_vector(dataset, scores)
    details = {"reference_scope": "mocked truth; real partition/activation/command/timing checks"}
    nested = {"dataset": dataset, "case": case, "checkpoint": case["checkpoint"],
        "checkpoint_sha": case["checkpoint_manifest"]["sha256"], "control": {"scores": vector},
        "candidate": row["candidate"], "golden": case["golden"]}
    plan = {"scientific_source": old["source_root"], "output": old["output"], "cpu": 16,
        "scientific_revision": reader.followup.REVISION, "dataset": dataset, "case": nested,
        "scoring_ids": {m: reader.followup.scoring_id(Path(old["output"]), m) for m in reader.followup.ARMS[1:]}, "inputs": []}
    seed = folder / "graph/orthogroups_profiles_refined.txt"
    seed.parent.mkdir()
    seed.write_text("a b\n")
    row["candidate"]["seed_partition"] = record(seed)
    graph_path = folder / "graph_metrics.json"
    graph = json.loads(graph_path.read_text())
    graph["stages"][0]["output"] = record(seed)
    graph_path.write_text(json.dumps(graph))
    graph_command, command = reader.commands(Path(plan["scientific_source"]), plan, case, folder,
        Path(nested["checkpoint"]), nested["checkpoint_sha"])
    if mode != "none": graph_command += ["--profile-sequence-weighting", mode]
    (folder / "preflight.json").write_text(json.dumps({"graph": graph_command, "native": command}))
    row.update(dataset=dataset, graph_execution=timing(folder, "graph", graph_command),
        native_execution=timing(folder, "native", command), scores=scores, details=details,
        metric_vector=vector, checked=[], graph_report=record(graph_path))
    native.update(genes=reader.followup.UNIVERSES[dataset][0], species=reader.followup.UNIVERSES[dataset][1])
    monkeypatch.setattr(reader, "replay_native", lambda *a: native)
    calls = []
    def assess(repo, observed_dataset, observed_row, observed_native, scoring_folder, participant):
        calls.append((scoring_folder, participant))
        assert observed_dataset == dataset and observed_native is native
        return deepcopy(scores), deepcopy(details), []
    monkeypatch.setattr(reader, "readback_scores", assess)
    if dataset == "QfO" and mode == "none":
        recovery = tmp_path / "benchmarks/results/canonical_qfo_scoring_recovery_r1_v1_20261014"
        recovery.mkdir(parents=True)
        recovery_plan = recovery / "execution_plan.json"
        recovery_plan.write_text(json.dumps({"scoring_id": "canonical_fixture"}))
        report = {"status": "canonical_qfo_scoring_recovered_independently_verified", "execution_plan": record(recovery_plan),
            "scores": scores, "conversion": {"fixture": True}, "metric_files": [], "references": []}
        (recovery / "results.json").write_text(json.dumps(report))
        row.update(conversion=report["conversion"], metric_files=[], references=[],
            control_scoring_reuse=record(recovery / "results.json"))
        plan["inputs"] = [record(recovery / name) for name in ("results.json", "execution_plan.json")]
    return plan, row, native, folder, calls


@pytest.mark.parametrize("dataset", reader.DATASETS)
@pytest.mark.parametrize("mode", reader.followup.ARMS)
def test_original_commands_activation_partitions_and_control_scoring_reuse_are_checked(tmp_path, monkeypatch, dataset, mode):
    plan, row, _, folder, calls = fixture(tmp_path, monkeypatch, dataset, mode)
    cell = reader.admit_arm(tmp_path, plan, row, ["a", "b"])
    assert cell["scores"] == row["metric_vector"] and cell["coverage"]["genes"] == reader.followup.UNIVERSES[dataset][0]
    assert cell["activation"] == row["activation"] and cell["checked"]
    if dataset == "QfO" and mode == "none":
        assert cell["control_assessment_reuse"] == row["control_scoring_reuse"]
        assert calls[0] == (tmp_path / "benchmarks/results/canonical_qfo_scoring_recovery_r1_v1_20261014",
                            "ohmm_hillclimb_r1_canonical_fixture")
    elif dataset == "QfO":
        assert calls[0] == (folder, "ohmm_hillclimb_r1_" + plan["scoring_ids"][mode])
        assert cell["control_assessment_reuse"] is None
    else: assert calls == [(folder, None)]


@pytest.mark.parametrize("kind", ["command", "activation", "graph", "candidate", "seed", "coverage", "universe", "summary",
                                  "scores", "details", "vector", "scoring_id", "reuse"])
def test_changed_followup_receipts_cannot_be_admitted(tmp_path, monkeypatch, kind):
    plan, original, native, _, _ = fixture(tmp_path, monkeypatch, "QfO")
    row = deepcopy(original)
    if kind == "command": row["graph_execution"]["command"].append("--unexpected")
    elif kind == "activation": row["activation"]["unexpected"] = True
    elif kind == "graph": row["graph_report"]["sha256"] = "wrong"
    elif kind == "candidate": row["candidate"]["candidate_partition"]["sha256"] = "wrong"
    elif kind == "seed": row["candidate"]["seed_partition"] = record(Path(row["graph_report"]["path"]))
    elif kind == "coverage": row["coverage"]["all_candidate_genes_preserved"] = False
    elif kind == "universe": native["genes"] -= 1
    elif kind == "summary": row["native_summary"]["unexpected"] = True
    elif kind == "scores": row["scores"]["VGNC F"] += .01
    elif kind == "details": row["details"]["reference_scope"] = "changed"
    elif kind == "vector": row["metric_vector"]["QfO: GO"] += .01
    elif kind == "scoring_id": plan["scoring_ids"]["henikoff"] = "wrong"
    else: row["control_scoring_reuse"] = {"unexpected": True}
    with pytest.raises(ValueError): reader.admit_arm(tmp_path, plan, row, ["a", "b"])


@pytest.mark.parametrize("kind", ["unbound", "changed_recovery", "status", "assessment"])
def test_qfo_control_reuse_requires_exact_recovered_assessment(tmp_path, monkeypatch, kind):
    plan, row, _, _, _ = fixture(tmp_path, monkeypatch, "QfO", "none")
    if kind == "unbound": plan["inputs"] = []
    elif kind == "changed_recovery": row["control_scoring_reuse"]["sha256"] = "wrong"
    elif kind == "assessment": row["conversion"] = {"changed": True}
    else:
        path = Path(row["control_scoring_reuse"]["path"])
        report = json.loads(path.read_text())
        report["status"] = "running"
        path.write_text(json.dumps(report))
        row["control_scoring_reuse"] = record(path)
        plan["inputs"][0] = record(path)
    with pytest.raises(ValueError, match="assessment reuse"):
        reader.admit_arm(tmp_path, plan, row, ["a", "b"])


@pytest.mark.parametrize("kind", ["valid", "duplicate", "empty", "ownership", "universe"])
def test_real_checkpoint_and_fasta_ownership_with_json_class_keys(tmp_path, monkeypatch, kind):
    fastas = []
    for i, name in enumerate(("a", "b")):
        path = tmp_path / f"sp{i}.faa"
        path.write_text(f">{name}\nAAAA\n")
        fastas.append(record(path))
    if kind in ("duplicate", "empty"):
        path = Path(fastas[1]["path"])
        path.write_text(">a\nAAAA\n" if kind == "duplicate" else ">b\n")
        fastas[1] = record(path)
    names = ["a", "b"]
    checkpoint = write_accuracy_checkpoint(str(tmp_path / "checkpoint"), names, np.array([0, 1], dtype=np.int32),
        np.array([0], dtype=np.int32), np.array([1], dtype=np.int32), np.array([1.], dtype=np.float64))
    plan = {"dataset": "ThreeKingdoms", "case": {"checkpoint": str(checkpoint), "case": {"fastas": fastas}}}
    report = {"ownership": {"genes": 2, "species": 2, "classes": {"0": "sp0.faa", "1": "sp1.faa"}}}
    monkeypatch.setattr(reader.followup, "UNIVERSES", {"ThreeKingdoms": (2, 2)})
    if kind == "ownership": report["ownership"]["classes"]["0"] = "other.faa"
    elif kind == "universe": reader.followup.UNIVERSES["ThreeKingdoms"] = (3, 2)
    if kind == "valid": assert reader.ownership(plan, report) == names
    else:
        with pytest.raises(ValueError): reader.ownership(plan, report)


def vector_fixture():
    pilot = {"cells": [{"label": mode, "dataset": "OrthoBench", "scores": {reader.ENDPOINTS[0]: .7 + (i % 2) * .001}}
                       for i, mode in enumerate((*reader.followup.ARMS, "species_henikoff_cluster"))]}
    cells = [{"label": mode, "dataset": dataset, "scores": {key: .7 for key in reader.ENDPOINTS if key.startswith(dataset + ":")}}
             for dataset in reader.DATASETS for mode in reader.followup.ARMS]
    return pilot, cells


def test_eight_endpoint_vectors_preserve_both_variants_and_stopped_missing_slots(tmp_path):
    pilot, cells = vector_fixture()
    methods = reader.assemble(pilot, cells)
    assert len(methods) == 4 and len(methods[0]["scores"]) == 8
    assert methods[1]["full_vector_screen"]["eligible"]
    assert not methods[2]["full_vector_screen"]["eligible"]
    assert len(methods[3]["missing_endpoints"]) == 7 and "full_vector_screen" not in methods[3]
    assert reader.dataset_deltas(cells[:3])["henikoff"] == {key: 0 for key in cells[0]["scores"]}
    table = tmp_path / "scores.tsv"
    reader.write_table(table, methods, [])
    rows = list(csv.DictReader(table.open(), delimiter="\t"))
    assert len(rows) == 4 and rows[-1][reader.ENDPOINTS[-1]] == "NA"
    assert all(len(row) == 10 for row in rows)


@pytest.mark.parametrize("kind", ["missing", "duplicate", "functional_regression", "supplementary_regression", "invalid"])
def test_no_cherry_picking_or_mixed_average_can_hide_missingness_or_regressions(kind):
    pilot, cells = vector_fixture()
    if kind == "missing": cells.pop()
    elif kind == "duplicate": cells.append(deepcopy(cells[0]))
    elif kind == "functional_regression": cells[1]["scores"]["QfO: GO"] -= .006
    elif kind == "supplementary_regression": cells[4]["scores"][reader.ENDPOINTS[-1]] -= .006
    else: cells[0]["scores"]["QfO: GO"] = float("nan")
    if kind in ("missing", "duplicate", "invalid"):
        with pytest.raises(ValueError): reader.assemble(pilot, cells)
    else:
        methods = reader.assemble(pilot, cells)
        assert not methods[1]["full_vector_screen"]["eligible"]
        assert methods[1]["full_vector_screen"]["declines_exceeding_0_005"]


@pytest.mark.parametrize("kind", ["valid", "partial", "order", "job", "source", "driver", "error", "defaults"])
def test_complete_panels_must_match_original_jobs_and_sources(monkeypatch, kind):
    monkeypatch.setattr(reader.followup, "validate_plan", lambda p: None)  # Dedicated follow-up tests cover this helper.
    plan = {"dataset": "QfO", "driver_revision": "6e684fb2d050563d03578d636e3de9325ecec73c"}
    report = {"status": "all_weighting_followup_arms_scored", "dataset": "QfO", "job_id": "24490",
        "source_revision": reader.followup.REVISION, "driver_revision": plan["driver_revision"],
        "defaults_promoted": False, "goal_completed": False,
        "arms": [{"dataset": "QfO", "label": mode, "status": "native_and_score_independently_verified"} for mode in reader.followup.ARMS]}
    entry = {"dataset": "QfO", "inference": {"job_id": "24490"}}
    if kind == "partial": report["arms"].pop()
    elif kind == "order": report["arms"].reverse()
    elif kind == "job": report["job_id"] = "other"
    elif kind == "source": report["source_revision"] = "other"
    elif kind == "driver": report["driver_revision"] = "other"
    elif kind == "error": report["arms"][1]["error"] = "failed"
    elif kind == "defaults": report["defaults_promoted"] = True
    if kind == "valid": reader.require_completed(plan, report, entry)
    else:
        with pytest.raises(ValueError, match="complete successful ordered"):
            reader.require_completed(plan, report, entry)


def test_existing_namespace_is_preserved(tmp_path):
    with pytest.raises(FileExistsError, match="Preserve every"):
        reader.run(tmp_path, tmp_path)


def test_exact_point_screen_boundary_matches_original_subtraction_rule():
    pilot, cells = vector_fixture()
    cells[1]["scores"]["QfO: GO"] = .7 - .005
    assert reader.assemble(pilot, cells)[1]["full_vector_screen"]["eligible"]
    cells[1]["scores"]["QfO: GO"] -= 1e-12
    assert not reader.assemble(pilot, cells)[1]["full_vector_screen"]["eligible"]


@pytest.mark.parametrize("kind", ["partial", "mixed", "endpoints"])
def test_per_dataset_deltas_require_complete_matching_endpoint_sets(kind):
    _, cells = vector_fixture()
    cells = cells[:3]
    if kind == "partial": cells.pop()
    elif kind == "mixed": cells[1]["dataset"] = "ThreeKingdoms"
    else: cells[1]["scores"].pop("QfO: GO")
    with pytest.raises(ValueError): reader.dataset_deltas(cells)


def test_launch_is_bound_to_original_immutable_commit_not_live_ledger(tmp_path, monkeypatch):
    entry = {"scientific_revision": reader.followup.REVISION,
             "datasets": [{"dataset": name} for name in reader.DATASETS]}
    def git(command, **kwargs):
        assert command == ["git", "show", reader.LAUNCH + ":" + reader.LEDGER]
        return json.dumps({"profile_weighting_followup_execution": entry}).encode()
    monkeypatch.setattr(reader.subprocess, "check_output", git)
    entries, binding = reader.launch(tmp_path)
    assert entries == entry["datasets"] and binding["revision"] == reader.LAUNCH
