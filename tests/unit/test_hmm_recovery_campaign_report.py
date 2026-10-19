from copy import deepcopy
import csv
import json
from pathlib import Path

import pytest

from benchmark_tools import report_hmm_recovery_campaign as campaign
from benchmark_tools import report_hmm_evidence_campaign as base
from benchmark_tools import report_projection_recovery as recovery_reader
from benchmark_tools.prepare_ob_candidate_neighborhood import record
from tests.unit.test_hmm_evidence_campaign_report import evidence as historical_evidence
from tests.unit.test_projection_recovery_readback import cells_and_rows


def evidence():
    weighting, local, _, _ = historical_evidence()
    cells, rows = cells_and_rows()
    recovery = {"context": deepcopy(weighting["context"]), "cells": cells,
        "backend_contrast": recovery_reader.verify_panel(cells, rows),
        "uncertainty": {"fixture": "recovery24"}}
    return weighting, local, recovery, {"job_id": "24514", "arms": rows}


def test_complete_recovery_collation_retains_controls_all_methods_and_original_failure(tmp_path):
    result = campaign.collate(*evidence())
    assert len(result["methods"]) == 13
    methods = {m["key"]: m for m in result["methods"]}
    assert methods["projection_recovery:legacy_off"]["role"] == "unchanged_control"
    assert methods["projection_recovery:projection_off"]["construction_role"] == "projection_only"
    assert methods["projection_recovery:projection_hmmer_score_only"]["role"] == "measurement_control"
    hmmer = methods["projection_recovery:projection_hmmer_reselect"]
    assert hmmer["point_gate"]["eligible"] and not hmmer["matched_projection_gate"]["eligible"]
    assert "repair-only" in hmmer["interpretation"] and result["backend_contrast"]["f1"] == 0
    failed = methods["selective_original:failed_panel"]
    assert failed["scores"] == {} and failed["missing_endpoints"] == list(base.ENDPOINTS)
    assert failed["key"] in result["frontier"]["excluded_incomplete"]
    assert result["uncertainty"]["projection_recovery"] == {"fixture": "recovery24"}
    assert len(result["functional_population_comparisons"]) == 6
    fas = [p for p in result["functional_population_comparisons"] if p["endpoint"] == "QfO: FAS"]
    assert len(fas) == 2 and all("not the numeric scored-sample count" in p["population_count_semantics"] for p in fas)
    path = tmp_path / "scores.tsv"
    base.write_table(path, result["methods"], result["context"])
    rows = list(csv.DictReader(path.open(), delimiter="\t"))
    assert len(rows) == 15 and all(len(row) == 10 for row in rows)
    assert all(rows[-1][key] == "NA" for key in base.ENDPOINTS)
    assert not any(a["method"].endswith("projection_hmmer_score_only") for a in result["required_next_actions"])
    assert "defaults unchanged" in result["recommendation"]


@pytest.mark.parametrize("kind", ["context", "canonical", "missing_weight", "missing_alignment", "weight_screen",
    "missing_local", "local_control", "local_gate", "missing_recovery", "recovery_control", "digest", "tree",
    "measurement_score", "legacy_gate", "matched_gate", "backend", "status"])
def test_selected_panels_or_changed_controls_and_gates_cannot_enter_recovery_report(kind):
    weight, local, recovery, producer = evidence()
    if kind == "context": local["context"][0]["scores"][base.ENDPOINTS[0]] += .01
    elif kind == "canonical":
        for report in (weight, local, recovery): report["context"][0]["scores"][base.ENDPOINTS[0]] += .01
    elif kind == "missing_weight": weight["methods"].pop(1)
    elif kind == "missing_alignment": weight["methods"].pop()
    elif kind == "weight_screen": weight["methods"][1]["full_vector_screen"]["eligible"] = False
    elif kind == "missing_local": local["cells"].pop()
    elif kind == "local_control": local["cells"][0]["details"]["precision"] += .01
    elif kind == "local_gate": local["cells"][1]["advancement"]["eligible"] = False
    elif kind == "missing_recovery": recovery["cells"].pop(3)
    elif kind == "recovery_control": recovery["cells"][0]["details"]["recall"] += .01
    elif kind == "digest": recovery["cells"][4]["activation"]["candidate_digest"] = "changed"
    elif kind == "tree": recovery["cells"][3]["rooted_species_tree"]["sha256"] = "changed"
    elif kind == "measurement_score": recovery["cells"][3]["scores"][base.ENDPOINTS[0]] += .01
    elif kind == "legacy_gate": producer["arms"][4]["advancement"]["eligible"] = False
    elif kind == "matched_gate": producer["arms"][4]["matched_projection_advancement"]["eligible"] = True
    elif kind == "backend": recovery["backend_contrast"]["f1"] = 1.
    else: recovery["cells"][4]["status"] = "running"
    with pytest.raises(ValueError): campaign.collate(weight, local, recovery, producer)


def test_neutral_local_stops_and_inactive_hmmer_keeps_original_fallback():
    weight, local, recovery, producer = evidence()
    local["cells"][1] = deepcopy(local["cells"][0])
    local["cells"][1]["label"] = "local_joint_v1"
    local["cells"][1]["advancement"] = base.advancement(base.point(local["cells"][0]), base.point(local["cells"][1]))
    cell, row = recovery["cells"][4], producer["arms"][4]
    cell["activation"].update(completely_measured_queries=0, applicability="no_complete_measurements")
    row.update(advancement=recovery_reader.base.selective_gate(base.point(recovery["cells"][0]), cell),
               accuracy_interpretation="no_complete_measurements")
    result = campaign.collate(weight, local, recovery, producer)
    methods = {m["key"]: m for m in result["methods"]}
    assert methods["local:local_joint_v1"]["status"] == "stopped_after_orthobench_point_gate"
    assert methods["projection_recovery:projection_hmmer_reselect"]["status"] == "applicability_fallback_required"
    assert not any(a["method"] == "local:local_joint_v1" for a in result["required_next_actions"])


@pytest.mark.parametrize("change", [None, "unsubmitted", "revision", "source", "dependency"])
def test_launch_uses_exact_original_successful_readers_and_separate_recovery_binding(tmp_path, monkeypatch, change):
    old = {key: {"fixture": key} for key in ("weighting", "local", "hmmer")}
    entry = {"runner_commit": "reader", "scientific_reader_source": "reader_source",
             "job": {"submitted": True, "dependency": "afterany:24514"}}
    science = {"inference": {"submitted": True, "job_id": "24514"},
               "independent_reader": {"runner_commit": "reader", "scientific_reader_source": "reader_source"}}
    if change == "unsubmitted": entry["job"]["submitted"] = False
    elif change == "revision": entry["runner_commit"] = "other"
    elif change == "source": entry["scientific_reader_source"] = "other"
    elif change == "dependency": entry["job"]["dependency"] = "afterany:other"
    monkeypatch.setattr(base, "launch", lambda _: (old, {"revision": base.LAUNCH}))
    calls = []

    def git(command, **kwargs):
        calls.append(command)
        return json.dumps({"projection_recovery_postrun_readback": entry, "projection_recovery_pilot": science}).encode()

    monkeypatch.setattr(campaign.subprocess, "check_output", git)
    if change:
        with pytest.raises(ValueError): campaign.launch(tmp_path)
    else:
        entries, actual, bindings = campaign.launch(tmp_path)
        assert set(entries) == {"weighting", "local", "recovery"} and actual == science
        assert entries["weighting"] == old["weighting"] and entries["recovery"] == entry
        assert bindings["historical"]["revision"] == base.LAUNCH
    assert calls == [["git", "show", campaign.LAUNCH + ":" + base.LEDGER]]


@pytest.mark.parametrize("kind", [None, "changed_bytes", "producer_state", "job_state"])
def test_original_failures_are_pinned_and_never_become_accuracy(kind, tmp_path, monkeypatch):
    paths = [tmp_path / name for name in ("failed.json", "reader.log", "collator.log")]
    paths[0].write_text(json.dumps({"job_id": "24493", "status": "wrong" if kind == "producer_state" else "finished_with_explicit_gaps"}))
    for path in paths[1:]: path.write_text("retained failure\n")
    spec = {"old_failure_records": [record(p) for p in paths]}
    if kind == "changed_bytes": paths[1].write_text("changed\n")

    def accounting(command, **kwargs):
        job = command[command.index("-j") + 1]
        return "JobIDRaw|State|ExitCode\n" + job + ("|COMPLETED|0:0\n" if kind == "job_state" else "|FAILED|1:0\n")

    monkeypatch.setattr(campaign.subprocess, "check_output", accounting)
    if kind:
        with pytest.raises(ValueError): campaign.failure_history(tmp_path, spec)
    else:
        history, pins = campaign.failure_history(tmp_path, spec)
        assert len(history["jobs"]) == len(pins) == 3
        assert history["status"] == "original_failures_preserved_not_accuracy_evidence"


@pytest.mark.parametrize("kind", [None, "failed_reader", "binding", "science_revision", "producer_job"])
def test_actual_report_writers_require_all_admissions_and_preserve_sources(tmp_path, monkeypatch, kind):
    weight, local, recovery, producer = evidence()
    producer_path = tmp_path / "producer.json"
    if kind == "producer_job": producer["job_id"] = "other"
    producer_path.write_text(json.dumps(producer))
    recovery.update(producer=record(producer_path), launch_ledger={"revision": "other" if kind == "binding" else "7f82"},
                    source_revision="other" if kind == "science_revision" else "science")
    reports = {"weighting": weight, "local": local, "recovery": recovery}
    entries = {key: {"fixture": key} for key in reports}
    entries["recovery"]["launch_binding"] = {"revision": "7f82"}
    science = {"runner_commit": "science", "inference": {"job_id": "24514"}}
    monkeypatch.setattr(campaign, "reporter_inventory", lambda _: ("reporter", []))
    monkeypatch.setattr(campaign, "launch", lambda _: (entries, science, {"fixture": "bound"}))
    monkeypatch.setattr(campaign, "failure_history", lambda *a: ({"fixture": "explicit failure history"}, []))
    monkeypatch.setattr(campaign, "retained_fas_samples", lambda *a: ({"fixture": "audited numeric populations"}, []))
    seen = []

    def read(repo, entry, status):
        key = entry["fixture"]
        seen.append((key, status))
        if kind == "failed_reader" and key == "recovery": raise ValueError("Unadmitted failed reader")
        return reports[key], [record(producer_path)], {"fixture": key}

    monkeypatch.setattr(base, "read_admission", read)
    before = record(producer_path)
    output = tmp_path / "report"
    if kind:
        with pytest.raises(ValueError): campaign.run(tmp_path, output)
        assert not output.exists()
    else:
        result = campaign.run(tmp_path, output)
        assert result["status"] == "admitted_recovery_campaign_collated_followups_explicit"
        assert result["defaults_promoted"] is result["goal_completed"] is False
        assert json.loads((output / "summary.json").read_text()) == result
        assert result["score_table"] == record(output / "scores.tsv")
        assert seen == list(campaign.REQUIRED.items())
        with pytest.raises(FileExistsError): campaign.run(tmp_path, output)
    assert record(producer_path) == before


def test_existing_report_is_rejected_before_any_admission(tmp_path):
    with pytest.raises(FileExistsError): campaign.run(tmp_path, tmp_path)


@pytest.mark.parametrize("kind", [None, "source", "pool", "sample", "modes"])
def test_fas_numeric_counts_are_bound_to_prior_audit_and_current_panel(tmp_path, monkeypatch, kind):
    weight, _, _, _ = evidence()
    populations = [{"mode": mode, "FAS_reported_eligible_pool": 100 if mode == "none" else 95,
                    "FAS_numeric_sample_pairs": 40} for mode in ("none", "henikoff", "henikoff_cluster")]
    if kind == "pool": populations[1]["FAS_reported_eligible_pool"] += 1
    elif kind == "sample": populations[1]["FAS_numeric_sample_pairs"] = 96
    elif kind == "modes": populations.reverse()
    blob = json.dumps({"weighting": {"functional_populations": populations}}).encode()
    path = tmp_path / campaign.FAS_EVIDENCE
    path.parent.mkdir(parents=True)
    path.write_bytes(blob + (b"\n" if kind == "source" else b""))
    monkeypatch.setattr(campaign.subprocess, "check_output", lambda *a, **k: blob)
    if kind:
        with pytest.raises(ValueError): campaign.retained_fas_samples(tmp_path, weight)
    else:
        result, pins = campaign.retained_fas_samples(tmp_path, weight)
        assert result["populations"] == populations and pins == [record(path)]
