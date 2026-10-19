from copy import deepcopy
import hashlib
import json
from pathlib import Path

import pytest

from benchmark_tools import plan_projection_followup as planner
from tests.unit.test_projection_followup_plan import authenticated_fixture, evidence


def stockholm_evidence(passing=None, inactive=()):
    report, producer = evidence(passing, inactive)
    report.update(source_revision=planner.STOCKHOLM_SCIENCE, reporter_revision=planner.STOCKHOLM_SCIENCE)
    for cell in report["cells"]: cell["source_revision"] = planner.STOCKHOLM_SCIENCE
    return report, producer


def select(report, producer):
    return planner.selection(report, producer, scientific_revision=planner.STOCKHOLM_SCIENCE,
                             reader_revision=planner.STOCKHOLM_SCIENCE)


@pytest.mark.parametrize("passing", [None, [], ["projection_off"], ["projection_builtin_reselect"], ["projection_hmmer_reselect"]])
@pytest.mark.parametrize("inactive", [(), ("projection_hmmer_reselect",)])
def test_stockholm_binding_preserves_every_original_gate_role_and_control_closure(passing, inactive):
    old_report, old_producer = evidence(passing, inactive)
    new_report, new_producer = stockholm_evidence(passing, inactive)
    assert select(new_report, new_producer) == planner.selection(old_report, old_producer)


@pytest.mark.parametrize("kind", ["old_report", "old_cell", "old_reader", "partial", "unadmitted", "gate", "backend", "promoted"])
def test_new_binding_cannot_admit_old_partial_or_changed_evidence(kind):
    report, producer = stockholm_evidence()
    if kind == "old_report": report["source_revision"] = planner.SCIENCE
    elif kind == "old_cell": report["cells"][4]["source_revision"] = planner.SCIENCE
    elif kind == "old_reader": report["reporter_revision"] = planner.READER
    elif kind == "partial": report["cells"].pop()
    elif kind == "unadmitted": report["accuracy_result_admitted"] = False
    elif kind == "gate": producer["arms"][4]["advancement"]["eligible"] = False
    elif kind == "backend": report["backend_contrast"]["f1"] += .001
    else: report["defaults_promoted"] = True
    with pytest.raises(ValueError): select(report, producer)


def launch_fixture(tmp_path, monkeypatch):
    binding = {"revision": planner.STOCKHOLM_SCIENCE_LAUNCH, "entry_key": "stockholm_annotation_recovery_pilot"}
    science = {"runner_commit": planner.STOCKHOLM_SCIENCE, "scientific_source": "science",
        "inference": {"job_id": "24519"}, "execution_plan": {"path": "plan.json", "bytes": 1, "sha256": "pin"},
        "independent_reader": {"runner_commit": planner.STOCKHOLM_SCIENCE, "scientific_reader_source": "reader"}}
    entry = {"runner_commit": planner.STOCKHOLM_SCIENCE, "scientific_reader_source": "reader", "launch_binding": deepcopy(binding),
        "job": {"submitted": True, "job_id": "24520", "dependency": "afterany:24519", "cpus": 1,
                "memory_allocation": "16G", "wall_limit": "02:00:00", "summary": "readback/summary.json"}}
    ledger = {"stockholm_annotation_recovery_pilot": deepcopy(science), "stockholm_annotation_recovery_postrun_readback": entry}
    pin = {**science["execution_plan"], "path": str(tmp_path / "plan.json")}
    monkeypatch.setattr(planner.subprocess, "check_output", lambda *a, **k: json.dumps(ledger).encode())
    calls = []

    def scientific_launch(repo, revision, *, entry_key):
        calls.append((repo, revision, entry_key))
        return science, pin, binding

    monkeypatch.setattr(planner.recovery_reader, "launch", scientific_launch)
    return ledger, science, pin, calls


def test_launch_binds_both_immutable_snapshots_and_separate_reader(tmp_path, monkeypatch):
    ledger, science, pin, calls = launch_fixture(tmp_path, monkeypatch)
    entries, actual, provenance = planner.launch(tmp_path, "stockholm_annotation")
    assert actual == science and entries["recovery"] == ledger["stockholm_annotation_recovery_postrun_readback"]
    assert calls == [(tmp_path, planner.STOCKHOLM_SCIENCE_LAUNCH, "stockholm_annotation_recovery_pilot")]
    assert provenance["recovery"]["git_blob_sha256"] == hashlib.sha256(json.dumps(ledger).encode()).hexdigest()


@pytest.mark.parametrize("kind", ["science", "parent", "reader", "reader_binding", "reader_path", "shared_source",
    "unsubmitted", "job", "dependency", "cpu", "memory", "time", "absolute_path", "escaped_path", "plan"])
def test_launch_rejects_stale_jobs_changed_resources_and_path_or_source_substitution(tmp_path, monkeypatch, kind):
    ledger, science, pin, _ = launch_fixture(tmp_path, monkeypatch)
    entry = ledger["stockholm_annotation_recovery_postrun_readback"]
    job = entry["job"]
    if kind == "science": ledger["stockholm_annotation_recovery_pilot"]["runner_commit"] = planner.SCIENCE
    elif kind == "parent": science["inference"]["job_id"] = "24514"
    elif kind == "reader": entry["runner_commit"] = planner.READER
    elif kind == "reader_binding": entry["launch_binding"]["revision"] = "old"
    elif kind == "reader_path": entry["scientific_reader_source"] = "other_reader"
    elif kind == "shared_source":
        science["scientific_source"] = entry["scientific_reader_source"]
        ledger["stockholm_annotation_recovery_pilot"] = deepcopy(science)
    elif kind == "unsubmitted": job["submitted"] = False
    elif kind == "job": job["job_id"] = "24515"
    elif kind == "dependency": job["dependency"] = "afterany:24514"
    elif kind == "cpu": job["cpus"] = 16
    elif kind == "memory": job["memory_allocation"] = "128G"
    elif kind == "time": job["wall_limit"] = "1-00:00:00"
    elif kind == "absolute_path": job["summary"] = "/other/summary.json"
    elif kind == "escaped_path": job["summary"] = "../other/summary.json"
    else: pin["sha256"] = "other"
    with pytest.raises(ValueError): planner.launch(tmp_path, "stockholm_annotation")


def test_stockholm_authentication_repeats_parent_runtime_source_and_selection(tmp_path, monkeypatch):
    report, producer, pilot, product = authenticated_fixture(tmp_path, monkeypatch)
    report["source_revision"] = report["reporter_revision"] = planner.STOCKHOLM_SCIENCE
    for cell in report["cells"]: cell["source_revision"] = planner.STOCKHOLM_SCIENCE
    producer.update(job_id="24519", source_revision=planner.STOCKHOLM_SCIENCE)
    pilot["source_revision"] = planner.STOCKHOLM_SCIENCE
    plan_path = Path(report["execution_plan"]["path"])
    plan_path.write_text(json.dumps(pilot))
    report["execution_plan"] = producer["execution_plan"] = planner.record(plan_path)
    product.write_text(json.dumps(producer))
    report["producer"] = planner.record(product)
    entry = {"runner_commit": planner.STOCKHOLM_SCIENCE, "launch_binding": report["launch_ledger"]}
    science = {"runner_commit": planner.STOCKHOLM_SCIENCE, "scientific_source": "science", "inference": {"job_id": "24519"}}
    monkeypatch.setattr(planner, "launch", lambda repo, name: ({"recovery": entry}, science, {"campaign": name}))
    monkeypatch.setattr(planner.subprocess, "check_output", lambda *a, **k: planner.STOCKHOLM_SCIENCE)
    actual, plan, chosen, pins, provenance = planner.authenticate(tmp_path, campaign_name="stockholm_annotation")
    assert actual is report and plan == pilot and chosen == select(report, producer)
    assert provenance["science_accounting"]["fixture"] == ["24519", 16, "128G"]
    assert provenance["launches"] == {"campaign": "stockholm_annotation"}
    for pin in pins: planner.check(pin)


@pytest.mark.parametrize("eligible", [True, False])
@pytest.mark.parametrize("change", [None, "launch", "source", "reader", "resources", "parent", "entry"])
def test_preparation_uses_explicit_protocol_and_never_runs_inference(tmp_path, monkeypatch, eligible, change):
    repo, driver, output = tmp_path / "repo", tmp_path / "driver", tmp_path / "prepared"
    repo.mkdir()
    folder = driver / "benchmark_tools"
    folder.mkdir(parents=True)
    spec = json.loads(Path(planner.__file__).with_name(planner.STOCKHOLM_PROTOCOL).read_text())
    if change == "launch": spec["scientific_launch"] = planner.recovery_campaign.LAUNCH
    elif change == "source": spec["scientific_revision"] = planner.SCIENCE
    elif change == "reader": spec["reader_revision"] = planner.READER
    elif change == "resources": spec["resources"]["cpus_per_run"] = 32
    elif change == "parent": spec["parent"] = "other.json"
    elif change == "entry": spec["launch_entry"] = "projection_recovery_pilot"
    (folder / planner.STOCKHOLM_PROTOCOL).write_text(json.dumps(spec))
    for name in (planner.PROTOCOL, "hmm_stockholm_annotation_recovery_20261019.json"):
        (folder / name).write_bytes(Path(planner.__file__).with_name(name).read_bytes())
    helper = folder / "plan_projection_followup.py"
    helper.write_text("frozen planner\n")
    report, producer = stockholm_evidence(None if eligible else [])
    chosen = select(report, producer)
    pin = planner.record(helper)
    report["execution_plan"] = pin
    pilot = {"source_root": "/frozen/stockholm", "runtime": {"files": []}, "tool_path": "/tools"}
    calls = []

    def authenticate(repo, *, campaign_name):
        assert campaign_name == "stockholm_annotation"
        return report, pilot, chosen, [], {"reader": {"reader": pin}}

    def selected(api, repo, dataset):
        calls.append(dataset)
        checkpoint = tmp_path / dataset / "checkpoint"
        checkpoint.mkdir(parents=True)
        return {"dataset": dataset, "checkpoint": str(checkpoint), "case": {"baseline_output": str(tmp_path / dataset / "native")}}, []

    monkeypatch.setattr(planner, "__file__", str(helper))
    monkeypatch.setattr(planner, "authenticate", authenticate)
    monkeypatch.setattr(planner, "selected_case", selected)
    monkeypatch.setattr(planner, "sources", lambda _: [planner.record(helper)])
    monkeypatch.setattr(planner.subprocess, "check_output", lambda command, **k: "planner_revision" if command[1] == "rev-parse"
                        else (driver / command[-1].split(":", 1)[1]).read_bytes())
    if change:
        with pytest.raises(ValueError): planner.prepare(repo, driver, output, campaign_name="stockholm_annotation")
        assert calls == [] and not output.exists()
    else:
        result = planner.prepare(repo, driver, output, campaign_name="stockholm_annotation")
        assert result["scientific_revision"] == planner.STOCKHOLM_SCIENCE and result["campaign"] == "stockholm_annotation"
        assert {pin["path"] for pin in result["preparer_files"]} >= {
            str(folder / name) for name in (planner.PROTOCOL, planner.STOCKHOLM_PROTOCOL,
                                          "hmm_stockholm_annotation_recovery_20261019.json")}
        assert calls == (list(planner.DATASETS) if eligible else [])
        assert result["status"] == ("projection_transfer_inputs_prepared_execution_driver_required" if eligible else "no_projection_followup_required")
        assert json.loads((output / "prepared_inputs.json").read_text()) == result


def test_unknown_campaign_never_falls_back_to_failed_original():
    with pytest.raises(ValueError): planner.configuration("stockholm")
