from copy import deepcopy
import csv
import hashlib
import json
from pathlib import Path

import pytest

from benchmark_tools import report_hmm_transfer_campaign as campaign
from benchmark_tools.prepare_ob_candidate_neighborhood import record
from tests.unit.test_completed_hmm_resources import cell_fixture
from tests.unit.test_hmm_evidence_campaign_report import evidence as historical_evidence
from tests.unit.test_stockholm_projection_followup_plan import stockholm_evidence


def evidence():
    weighting, local, _, _ = historical_evidence()
    pilot, producer = stockholm_evidence(passing=["projection_off"])
    pilot.update(context=deepcopy(weighting["context"]), uncertainty={"fixture": "original paired family uncertainty"})
    selected = campaign.planner.selection(pilot, producer, scientific_revision=campaign.SCIENCE,
                                          reader_revision=campaign.SCIENCE)
    cells = []
    for dataset in campaign.transfer.DATASETS:
        for label in selected["execution_arms"]:
            cell = {"label": label, "dataset": dataset, "role": selected["roles"][label],
                    "scores": {key: .7 for key in campaign.base.ENDPOINTS if key.startswith(dataset + ":")}}
            if dataset == "QfO":
                cell["details"] = {key: {"assessed_relations": 100, "similarity": .7} for key in ("EC", "GO", "FAS")}
                cell["functional_populations"] = {key: {"sample_pairs" if key == "FAS" else "numeric_scored_pairs": 40}
                                                   for key in ("EC", "GO", "FAS")}
            cells.append(cell)
    methods, contrasts = campaign.transfer.assemble(pilot, cells, selected)
    following = {"status": "projection_full_vectors_independently_read_back", "accuracy_result_admitted": True,
        "defaults_promoted": False, "goal_completed": False, "context": deepcopy(pilot["context"]),
        "orthobench_uncertainty": deepcopy(pilot["uncertainty"]), "selection": selected, "cells": cells,
        "methods": methods, "matched_contrasts": contrasts,
        "eligible_labels": [m["cells"][0]["label"] for m in methods if m.get("full_vector_screen", {}).get("eligible")]}
    return weighting, local, pilot, producer, following


def test_complete_merge_preserves_original_cells_gates_controls_and_stopped_branches(tmp_path):
    weight, local, pilot, producer, following = evidence()
    original = campaign.recovery.collate(weight, local, pilot, producer)
    before = deepcopy(original)
    result = campaign.merge_transfer(original, pilot, producer, following)
    methods = {m["key"]: m for m in result["methods"]}
    for old in before["methods"]:
        new = methods[old["key"]]
        if not old["key"].startswith("projection_recovery:"):
            assert new == old
            continue
        assert new["cells"][0] == old["cells"][0]
        assert new["point_gate"] == old["point_gate"]
        assert new["construction_role"] == old["construction_role"]
        assert new["pilot_role"] == old["role"]
        if "matched_projection_gate" in old: assert new["matched_projection_gate"] == old["matched_projection_gate"]
    candidate = methods["projection_recovery:projection_off"]
    assert candidate["role"] == "advancing_candidate" and not candidate["missing_endpoints"]
    assert candidate["full_vector_screen"]["eligible"]
    for label in ("projection_builtin_reselect", "projection_hmmer_reselect"):
        stopped = methods["projection_recovery:" + label]
        assert stopped["role"] == "stopped_candidate"
        assert stopped["missing_endpoints"] == list(campaign.base.ENDPOINTS[1:])
        assert not any(a["method"] == stopped["key"] for a in result["required_next_actions"])
    assert methods["projection_recovery:projection_hmmer_score_only"]["role"] == "measurement_control"
    assert methods["selective_original:failed_panel"]["scores"] == {}
    assert methods["projection_recovery:legacy_off"]["scores"] == methods["profile_weighting:none"]["scores"]
    assert result["uncertainty"]["projection_recovery"] == pilot["uncertainty"]
    assert result["matched_transfer_contrasts"] == following["matched_contrasts"]
    assert len(result["functional_population_comparisons"]) == 9
    assert all(p["candidate_numeric_scored_pairs"] == 40 for p in result["functional_population_comparisons"][-3:])
    action, = [a for a in result["required_next_actions"] if a["method"] == candidate["key"]]
    assert "justified combinations" in action["action"] and "Run unchanged corrected" not in action["action"]
    assert "defaults unchanged" in result["recommendation"]
    path = tmp_path / "scores.tsv"
    campaign.base.write_table(path, result["methods"], result["context"])
    rows = list(csv.DictReader(path.open(), delimiter="\t"))
    assert len(rows) == 15 and all(len(r) == 10 for r in rows)
    assert all(rows[-1][key] == "NA" for key in campaign.base.ENDPOINTS)


@pytest.mark.parametrize("kind", ["status", "admission", "defaults", "closed", "context", "uncertainty", "selection",
    "method", "contrast", "eligible", "missing_transfer", "missing_endpoint", "canonical", "pilot_cell"])
def test_changed_partial_or_reselected_transfer_is_rejected(kind):
    weight, local, pilot, producer, following = evidence()
    result = campaign.recovery.collate(weight, local, pilot, producer)
    if kind == "status": following["status"] = "running"
    elif kind == "admission": following["accuracy_result_admitted"] = False
    elif kind == "defaults": following["defaults_promoted"] = True
    elif kind == "closed": following["goal_completed"] = True
    elif kind == "context": following["context"].reverse()
    elif kind == "uncertainty": following["orthobench_uncertainty"] = {}
    elif kind == "selection": following["selection"]["roles"]["projection_off"] = "measurement_control"
    elif kind == "method": following["methods"][1]["scores"] = {}
    elif kind == "contrast": following["matched_contrasts"]["versus_projection"] = {}
    elif kind == "eligible": following["eligible_labels"] = []
    elif kind == "missing_transfer": following["cells"].pop()
    elif kind == "missing_endpoint": following["cells"][1]["scores"].pop("QfO: GO")
    elif kind == "canonical":
        following["cells"][0]["scores"]["QfO: GO"] += .001
        following["methods"], following["matched_contrasts"] = campaign.transfer.assemble(pilot, following["cells"], following["selection"])
    else:
        result["methods"][next(i for i, m in enumerate(result["methods"]) if m["key"] == "projection_recovery:projection_off")]["cells"][0] = {"changed": True}
    with pytest.raises(ValueError): campaign.merge_transfer(result, pilot, producer, following)


def test_original_full_vector_failure_stops_branch_without_retuning():
    weight, local, pilot, producer, following = evidence()
    following["cells"][1]["scores"]["QfO: GO"] -= .006
    following["methods"], following["matched_contrasts"] = campaign.transfer.assemble(pilot, following["cells"], following["selection"])
    following["eligible_labels"] = []
    result = campaign.merge_transfer(campaign.recovery.collate(weight, local, pilot, producer), pilot, producer, following)
    action, = [a for a in result["required_next_actions"] if a["method"] == "projection_recovery:projection_off"]
    assert "Stop this fixed branch without retuning" in action["action"]
    candidate, = [m for m in result["methods"] if m["key"] == action["method"]]
    assert not candidate["full_vector_screen"]["eligible"]


@pytest.mark.parametrize("kind", [None, "pending", "revision", "source", "submitted", "cpu", "memory", "time",
    "campaign", "report_binding", "binding", "dependency", "reader_revision", "reader_source"])
def test_transfer_admission_uses_immutable_reader_and_scientific_bindings(tmp_path, monkeypatch, kind):
    source = "benchmarks/work/stockholm_projection_followup_readback_v1_20261019"
    binding = {"revision": "83ac", "git_blob_sha256": "science"}
    entry = {"runner_commit": campaign.READER, "scientific_reader_source": source, "launch_binding": deepcopy(binding),
        "job": {"submitted": True, "job_id": "24525", "cpus": 1, "memory_allocation": "16G",
                "wall_limit": "02:00:00", "dependency": "afterany:24524"}}
    report = {"campaign": "stockholm_annotation", "launch_ledger": deepcopy(binding)}
    science = {"datasets": [{"inference": {"job_id": "24524"}}],
        "independent_reader": {"runner_commit": campaign.READER, "scientific_reader_source": source}}
    if kind == "revision": entry["runner_commit"] = "other"
    elif kind == "source": entry["scientific_reader_source"] = "other"
    elif kind == "submitted": entry["job"]["submitted"] = False
    elif kind in ("cpu", "memory", "time"):
        key, value = {"cpu": ("cpus", 16), "memory": ("memory_allocation", "128G"), "time": ("wall_limit", "1-00:00:00")}[kind]
        entry["job"][key] = value
    elif kind == "campaign": report["campaign"] = "projection"
    elif kind == "report_binding": report["launch_ledger"] = {}
    elif kind == "binding": entry["launch_binding"] = {"revision": "other"}
    elif kind == "dependency": entry["job"]["dependency"] = "afterany:other"
    elif kind == "reader_revision": science["independent_reader"]["runner_commit"] = "other"
    elif kind == "reader_source": science["independent_reader"]["scientific_reader_source"] = "other"
    blob = json.dumps({"stockholm_projection_followup_postrun_readback": entry}).encode()
    calls = []
    monkeypatch.setattr(campaign.subprocess, "check_output", lambda command, **kw: calls.append(command) or blob)
    seen = []
    def read(repo, actual, status):
        seen.append((actual, status))
        if kind == "pending": raise ValueError("Reader pending")
        return report, [], {"fixture": "readback"}
    monkeypatch.setattr(campaign.base, "read_admission", read)
    monkeypatch.setattr(campaign.transfer, "launch", lambda repo, revision, *, campaign_name: (science, binding))
    if kind:
        with pytest.raises(ValueError): campaign.admit_transfer(tmp_path)
    else:
        actual, pins, provenance = campaign.admit_transfer(tmp_path)
        assert actual is report and pins == []
        assert provenance["submission"]["git_blob_sha256"] == hashlib.sha256(blob).hexdigest()
        assert seen == [(entry, "projection_full_vectors_independently_read_back")]
    assert calls == [["git", "show", campaign.READBACK_LAUNCH + ":" + campaign.base.LEDGER]]


@pytest.mark.parametrize("kind", [None, "bytes", "producer", "job", "state"])
def test_original_projection_failures_are_pinned_and_not_biological_results(tmp_path, monkeypatch, kind):
    files = [tmp_path / name for name in ("producer.json", "reader.log", "report.log", "preparation.log")]
    files[0].write_text(json.dumps({"job_id": "other" if kind == "job" else "24514",
        "status": "running" if kind == "producer" else "finished_with_explicit_gaps"}))
    for path in files[1:]: path.write_text("retained failed execution\n")
    spec = {"projection_failure_records": [record(p) for p in files]}
    if kind == "bytes": files[1].write_text("changed\n")
    def accounting(command, **kw):
        job = command[command.index("-j") + 1]
        return "JobIDRaw|State|ExitCode\n" + job + ("|COMPLETED|0:0\n" if kind == "state" else "|FAILED|1:0\n")
    monkeypatch.setattr(campaign.subprocess, "check_output", accounting)
    if kind:
        with pytest.raises(ValueError): campaign.projection_failures(tmp_path, spec)
    else:
        history, pins = campaign.projection_failures(tmp_path, spec)
        assert len(history["jobs"]) == len(pins) == 4
        assert "never negative accuracy" in history["scope"]


def test_resource_table_retains_real_phase_receipts_local_scope_and_missing_values(tmp_path):
    cell, _ = cell_fixture(tmp_path / "weighted")
    local, _ = cell_fixture(tmp_path / "local", local=True)
    output = tmp_path / "output"
    output.mkdir()
    rows, pins = campaign.resource_table(output, [{"key": "profile_weighting:henikoff", "cells": [cell]},
        {"key": "local:local_joint_v1", "cells": [local]}, {"key": "failed", "cells": []}])
    assert len(rows) == 5 and len(pins) == 6
    assert rows[-1]["phase"] == "native" and "excludes nomination" in rows[-1]["scope"]
    text = (output / "resources.tsv").read_text()
    assert "\tNA\t" in text and "raw_GNU_time" in text


@pytest.mark.parametrize("kind", [None, "dirty_source", "dirty_protocol"])
def test_reporter_inventory_requires_exact_committed_bytes_for_source_and_both_protocols(tmp_path, monkeypatch, kind):
    folder = tmp_path / "benchmark_tools"
    folder.mkdir()
    helper = folder / "report_hmm_transfer_campaign.py"
    helper.write_text("committed helper\n")
    for name in (campaign.PROTOCOL, campaign.recovery.PROTOCOL): (folder / name).write_text("{}\n")
    monkeypatch.setattr(campaign, "sources", lambda _: [record(helper)])
    def git(command, **kw):
        if command[1] == "rev-parse": return "revision\n"
        path = tmp_path / command[-1].split(":", 1)[1]
        if (kind == "dirty_source" and path == helper) or (kind == "dirty_protocol" and path.name == campaign.PROTOCOL):
            return b"different committed bytes\n"
        return path.read_bytes()
    monkeypatch.setattr(campaign.subprocess, "check_output", git)
    if kind:
        with pytest.raises(ValueError): campaign.reporter_inventory(tmp_path)
    else:
        revision, pins = campaign.reporter_inventory(tmp_path)
        assert revision == "revision" and len(pins) == 3


@pytest.mark.parametrize("kind", [None, "pending_transfer", "failed_pilot_reader", "launch_binding", "science", "producer", "changed_pin"])
def test_full_report_lifecycle_requires_all_admissions_before_output(tmp_path, monkeypatch, kind):
    weight, local, pilot, producer, following = evidence()
    producer["job_id"] = "other" if kind == "producer" else "24519"
    producer_path = tmp_path / "producer.json"
    producer_path.write_text(json.dumps(producer))
    pilot.update(producer=record(producer_path), source_revision="other" if kind == "science" else campaign.SCIENCE,
                 launch_ledger={"revision": "other" if kind == "launch_binding" else "8731"})
    if kind == "changed_pin": producer_path.write_text("changed\n")
    reports = {"weighting": weight, "local": local, "recovery": pilot}
    entries = {key: {"fixture": key} for key in reports}
    entries["recovery"]["launch_binding"] = {"revision": "8731"}
    monkeypatch.setattr(campaign, "reporter_inventory", lambda _: ("reporter", []))
    def admitted(repo):
        if kind == "pending_transfer": raise ValueError("Reader24525 pending")
        return following, [], {"fixture": "complete transfer"}
    monkeypatch.setattr(campaign, "admit_transfer", admitted)
    monkeypatch.setattr(campaign.base, "launch", lambda _: (entries, {"fixture": "historical"}))
    monkeypatch.setattr(campaign.planner, "launch", lambda repo, name: (entries,
        {"runner_commit": campaign.SCIENCE, "inference": {"job_id": "24519"}}, {"fixture": name}))
    seen = []
    def read(repo, entry, status):
        key = entry["fixture"]
        seen.append((key, status))
        if kind == "failed_pilot_reader" and key == "recovery": raise ValueError("Reader24520 not admitted")
        return reports[key], [record(producer_path)], {"fixture": key}
    monkeypatch.setattr(campaign.base, "read_admission", read)
    monkeypatch.setattr(campaign.recovery, "retained_fas_samples", lambda *a: ({"fixture": "numeric FAS sample audit"}, []))
    monkeypatch.setattr(campaign.recovery, "failure_history", lambda *a: ({"fixture": "old HMMER failure"}, []))
    monkeypatch.setattr(campaign, "projection_failures", lambda *a: ({"fixture": "old projection failure"}, []))
    def resources(output, methods):
        (output / "resources.tsv").write_text("method\tphase\n")
        return [], []
    monkeypatch.setattr(campaign, "resource_table", resources)
    output = tmp_path / "report"
    before = record(producer_path)
    if kind:
        with pytest.raises(ValueError): campaign.run(tmp_path, output)
        assert not output.exists()
        if kind == "pending_transfer": assert not seen
    else:
        result = campaign.run(tmp_path, output)
        assert result["status"] == "admitted_individual_hmm_methods_and_projection_transfers_collated"
        assert result["accuracy_result_admitted"] is True
        assert result["defaults_promoted"] is result["goal_completed"] is False
        assert json.loads((output / "summary.json").read_text()) == result
        assert result["score_table"] == record(output / "scores.tsv")
        assert result["resource_table"] == record(output / "resources.tsv")
        failed = result["methods"][-1]
        assert failed["role"] == "failed_execution" and failed["scores"] == {}
        assert failed["key"] in result["frontier"]["excluded_incomplete"]
        assert seen == list(campaign.recovery.REQUIRED.items())
        with pytest.raises(FileExistsError): campaign.run(tmp_path, output)
    assert record(producer_path) == before


def test_existing_report_rejected_before_any_reads(tmp_path):
    with pytest.raises(FileExistsError): campaign.run(tmp_path, tmp_path)
