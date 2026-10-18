from copy import deepcopy
import csv
import json
from pathlib import Path

import pytest

from benchmark_tools import report_hmm_evidence_campaign as campaign
from benchmark_tools import report_profile_evidence_pilot as profile
from benchmark_tools.prepare_ob_candidate_neighborhood import record
from tests.unit.test_selective_profile_readback import cells_and_rows


def evidence():
    hmmer_cells, rows = cells_and_rows()
    for cell in hmmer_cells:
        cell.update(dataset="OrthoBench", status="independent_content_readback_verified")
    backend = profile.verify_selective_panel(hmmer_cells, rows)
    pilot = []
    for i, label in enumerate(("none", "henikoff", "henikoff_cluster", "species_henikoff_cluster")):
        cell = deepcopy(hmmer_cells[0 if i in (0, 3) else 1])
        cell["label"] = label
        if i == 3: cell["scores"][campaign.ENDPOINTS[0]] -= .001
        pilot.append(cell)
    following = []
    for dataset in ("QfO", "ThreeKingdoms"):
        for label in ("none", "henikoff", "henikoff_cluster"):
            scores = {key: .7 for key in campaign.ENDPOINTS if key.startswith(dataset + ":")}
            following.append({"label": label, "dataset": dataset, "scores": scores,
                "status": "independent_content_readback_verified",
                "details": {name: {"assessed_relations": 100 if label == "none" else 95, "similarity": .7}
                            for name in ("EC", "GO", "FAS")} if dataset == "QfO" else {"reference": "fixture"}})
    weighting = {"cells": following, "methods": campaign.assemble({"cells": pilot}, following),
                 "orthobench_uncertainty": {"fixture": "not biological evidence"}}
    guarded = deepcopy(pilot[3])
    guarded["label"] = "strict_span_guard"
    weighting["methods"].append({"key": "alignment:strict_span_guard", "label": "OrthoHMM alignment guard",
        "status": "stopped_after_orthobench_point_gate", "scores": guarded["scores"], "cells": [guarded],
        "missing_endpoints": list(campaign.ENDPOINTS[1:])})
    context = [{"key": "canonical_collection:canonical", "label": "Canonical OrthoHMM", "status": "retained_context",
                "scores": deepcopy(weighting["methods"][0]["scores"])},
               {"key": "orthofinder_3_1_5_full", "label": "OrthoFinder 3.1.5 full", "status": "retained_context",
                "scores": dict.fromkeys(campaign.ENDPOINTS, .8)}]
    weighting["context"] = deepcopy(context)
    local_cells = deepcopy([hmmer_cells[0], hmmer_cells[1]])
    local_cells[1]["label"] = "local_joint_v1"
    for cell in local_cells: cell.pop("activation")
    local_cells[1]["advancement"] = campaign.advancement(campaign.point(local_cells[0]), campaign.point(local_cells[1]))
    local = {"context": deepcopy(context), "cells": local_cells, "uncertainty": {"fixture": "local"}}
    hmmer = {"context": deepcopy(context), "cells": hmmer_cells, "backend_contrast": backend, "uncertainty": {"fixture": "selective"}}
    return weighting, local, hmmer, {"arms": rows}


def test_collation_preserves_all_mechanisms_controls_populations_and_separate_next_actions(tmp_path):
    result = campaign.collate(*evidence())
    assert len(result["methods"]) == 11
    roles = {method["key"]: method for method in result["methods"]}
    assert roles["selective:hmmer_score_only"]["role"] == "measurement_control"
    assert roles["selective:hmmer_score_only"]["point_gate"] is None
    assert roles["local:off"]["role"] == "unchanged_control"
    assert len(result["required_next_actions"]) == 5
    assert not any(action["method"].endswith("hmmer_score_only") for action in result["required_next_actions"])
    assert roles["profile_weighting:species_henikoff_cluster"]["missing_endpoints"] == list(campaign.ENDPOINTS[1:])
    assert roles["alignment:strict_span_guard"]["missing_endpoints"] == list(campaign.ENDPOINTS[1:])
    assert "local:local_joint_v1" in result["frontier"]["excluded_incomplete"]
    assert len(result["functional_population_comparisons"]) == 6
    assert result["functional_population_comparisons"][0]["candidate_assessed_relations"] == 95
    assert "defaults unchanged" in result["recommendation"]
    path = tmp_path / "scores.tsv"
    campaign.write_table(path, result["methods"], result["context"])
    table = list(csv.DictReader(path.open(), delimiter="\t"))
    assert len(table) == 13 and all(len(row) == 10 for row in table)
    assert table[-1][campaign.ENDPOINTS[-1]] == "NA"


@pytest.mark.parametrize("kind", ["context", "canonical", "missing_weight", "missing_alignment", "weight_vector", "weight_screen",
    "missing_local", "missing_hmmer", "local_control", "hmmer_control", "local_gate", "hmmer_gate", "backend", "dataset", "status"])
def test_no_changed_control_gate_population_snapshot_or_selected_subset_can_enter_campaign(kind):
    weight, local, hmmer, producer = evidence()
    if kind == "context": local["context"][0]["scores"][campaign.ENDPOINTS[0]] += .01
    elif kind == "canonical":
        for report in (weight, local, hmmer): report["context"][0]["scores"][campaign.ENDPOINTS[0]] += .01
    elif kind == "missing_weight": weight["methods"].pop(2)
    elif kind == "missing_alignment": weight["methods"].pop()
    elif kind == "weight_vector": weight["methods"][1]["scores"]["QfO: GO"] += .01
    elif kind == "weight_screen": weight["methods"][1]["full_vector_screen"]["eligible"] = False
    elif kind == "missing_local": local["cells"].pop()
    elif kind == "missing_hmmer": hmmer["cells"].pop(2)
    elif kind == "local_control": local["cells"][0]["details"]["precision"] -= .01
    elif kind == "hmmer_control": hmmer["cells"][0]["details"]["recall"] -= .01
    elif kind == "local_gate": local["cells"][1]["advancement"]["eligible"] = False
    elif kind == "hmmer_gate": producer["arms"][3]["advancement"]["eligible"] = False
    elif kind == "backend": hmmer["backend_contrast"]["f1"] = 1.
    elif kind == "dataset": local["cells"][1]["dataset"] = "QfO"
    else: local["cells"][1]["status"] = "running"
    with pytest.raises(ValueError): campaign.collate(weight, local, hmmer, producer)


def test_unmeasured_hmmer_requires_fixed_applicability_fallback_not_negative_accuracy():
    weight, local, hmmer, producer = evidence()
    cell, row = hmmer["cells"][3], producer["arms"][3]
    cell["activation"].update(completely_measured_queries=0, applicability="no_complete_measurements")
    row["advancement"] = profile.selective_gate(campaign.point(hmmer["cells"][0]), cell)
    row["accuracy_interpretation"] = "no_complete_measurements"
    result = campaign.collate(weight, local, hmmer, producer)
    method, = [m for m in result["methods"] if m["key"] == "selective:hmmer_reselect"]
    assert method["status"] == "applicability_fallback_required" and not method["point_gate"]["eligible"]
    action, = [a for a in result["required_next_actions"] if a["method"] == method["key"]]
    assert "QfO then Three Kingdoms" in action["action"]


def test_full_vector_regression_stops_weighted_branch_without_changing_other_methods():
    weight, local, hmmer, producer = evidence()
    weight["cells"][1]["scores"]["QfO: GO"] -= .006
    pilot = {"cells": [method["cells"][0] for method in weight["methods"][:4]]}
    weight["methods"][:4] = campaign.assemble(pilot, weight["cells"])
    result = campaign.collate(weight, local, hmmer, producer)
    candidate = result["methods"][1]
    assert candidate["campaign_decision"] == "fixed_full_vector_branch_stopped_without_retuning"
    assert not any(a["method"] == candidate["key"] for a in result["required_next_actions"])
    assert any(a["method"] == "profile_weighting:henikoff_cluster" for a in result["required_next_actions"])


def test_equal_decision_control_does_not_become_hmmer_backend_specific_benefit():
    result = campaign.collate(*evidence())
    assert result["backend_contrast"]["f1"] == 0
    assert "Off-relative gain alone" in result["backend_contrast"]["scope"]


@pytest.mark.parametrize("kind", ["valid", "partial", "revision", "defaults", "empty_checked", "inventory"])
def test_reader_admission_requires_own_terminal_resource_check_and_frozen_source(tmp_path, monkeypatch, kind):
    source = tmp_path / "source"
    source.mkdir()
    helper = source / "helper.py"
    helper.write_text("original\n")
    path, table = tmp_path / "summary.json", tmp_path / "scores.tsv"
    table.write_text("fixture\n")
    report = {"status": "admitted", "accuracy_result_admitted": True, "reporter_revision": "revision",
        "defaults_promoted": False, "goal_completed": False, "checked": [record(table)],
        "reporter_files": [record(helper)], "score_table": record(table)}
    if kind == "partial": report["status"] = "running"
    elif kind == "revision": report["reporter_revision"] = "other"
    elif kind == "defaults": report["defaults_promoted"] = True
    elif kind == "empty_checked": report["checked"] = []
    elif kind == "inventory": report["reporter_files"][0]["sha256"] = "wrong"
    path.write_text(json.dumps(report))
    entry = {"scientific_reader_source": "source", "runner_commit": "revision", "job": {"job_id": "job", "summary": "summary.json"}}
    calls = []
    monkeypatch.setattr(campaign, "successful_job", lambda *a: calls.append(a) or {"State": "COMPLETED"})
    monkeypatch.setattr(campaign.subprocess, "check_output", lambda *a, **kw: "revision\n")
    monkeypatch.setattr(campaign, "sources", lambda *a: [record(helper)])
    if kind == "valid":
        value, pins, provenance = campaign.read_admission(tmp_path, entry, "admitted")
        assert value == report and record(path) in pins and provenance["revision"] == "revision"
    else:
        with pytest.raises(ValueError): campaign.read_admission(tmp_path, entry, "admitted")
    assert calls == [("job", 1, "16G")]


def test_duplicate_content_pins_are_deduplicated_but_conflicts_are_rejected(tmp_path):
    path = tmp_path / "data"
    path.write_text("original\n")
    pin = record(path)
    assert campaign.unique_pins([pin, {**pin, "scope": "extra metadata"}]) == [pin]
    with pytest.raises(ValueError, match="Conflicting"):
        campaign.unique_pins([pin, {**pin, "sha256": "wrong"}])


def test_real_output_table_and_summary_preserve_sources_and_require_fresh_namespace(tmp_path, monkeypatch):
    weight, local, hmmer, producer = evidence()
    source = tmp_path / "producer.json"
    source.write_text(json.dumps(producer))
    hmmer["producer"] = record(source)
    reports = {"weighting": weight, "local": local, "hmmer": hmmer}
    monkeypatch.setattr(campaign, "launch", lambda repo: ({key: key for key in reports}, {"revision": campaign.LAUNCH}))
    monkeypatch.setattr(campaign, "read_admission", lambda repo, key, status: (reports[key], [record(source)], {"fixture": key}))
    monkeypatch.setattr(campaign, "sources", lambda *a: [])
    monkeypatch.setattr(campaign.subprocess, "check_output", lambda *a, **kw: "test_source_revision\n")
    before = record(source)
    output = tmp_path / "report"
    result = campaign.run(tmp_path, output)
    assert result["status"] == "admitted_campaign_results_collated_followup_decisions_explicit"
    assert result["defaults_promoted"] is result["goal_completed"] is False
    assert result["score_table"] == record(output / "scores.tsv")
    assert json.loads((output / "summary.json").read_text()) == result
    assert before == record(source)
    with pytest.raises(FileExistsError): campaign.run(tmp_path, output)


def test_launch_uses_original_committed_reader_submission_ledger(tmp_path, monkeypatch):
    ledger = {name: {"fixture": key} for key, (name, _) in campaign.READERS.items()}
    def git(command, **kwargs):
        assert command == ["git", "show", campaign.LAUNCH + ":" + campaign.LEDGER]
        return json.dumps(ledger).encode()
    monkeypatch.setattr(campaign.subprocess, "check_output", git)
    entries, binding = campaign.launch(tmp_path)
    assert set(entries) == set(campaign.READERS) and binding["revision"] == campaign.LAUNCH
