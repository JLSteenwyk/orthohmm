from copy import deepcopy
import json

import pytest

from benchmark_tools import report_additional_hillclimb as reporter
from benchmark_tools import report_shared_qfo_export as recovery
from benchmark_tools.metric_frontier import ENDPOINTS
from benchmark_tools.prepare_ob_candidate_neighborhood import record
from benchmark_tools.report_hillclimb import endpoint_vector
from benchmark_tools.run_membership_hillclimb import PRODUCTS, write_json
from benchmark_tools.shared_qfo_export import normalize, verify_projection
from benchmark_tools.summarize_matched_resources import LABELS

ACCOUNTING = "JobIDRaw|State|ExitCode|AllocCPUS|NodeList\n24360|COMPLETED|0:0|8|bizon\n"


def test_terminal_resources_required(monkeypatch):
    monkeypatch.setattr(recovery.subprocess, "check_output", lambda *args, **kwargs: ACCOUNTING)
    assert recovery.terminal({"job_id": "24360"})["AllocCPUS"] == "8"
    with pytest.raises(ValueError):
        recovery.terminal({"job_id": "24247"})


@pytest.mark.parametrize("text", [ACCOUNTING.replace("COMPLETED", "RUNNING"),
    ACCOUNTING.replace("COMPLETED", "FAILED"), ACCOUNTING.replace("0:0", "1:0"),
    ACCOUNTING.replace("|8|", "|16|"), ACCOUNTING.replace("bizon", "other"),
    ACCOUNTING + "24360|COMPLETED|0:0|8|bizon\n"])
def test_nonterminal_failed_or_wrong_resources_refused(monkeypatch, text):
    monkeypatch.setattr(recovery.subprocess, "check_output", lambda *args, **kwargs: text)
    with pytest.raises(ValueError):
        recovery.terminal({"job_id": "24360"})


def panel(tmp_path, monkeypatch, status):
    root = tmp_path / "producer"
    root.mkdir()
    write_json(root / "results.json", {"status": status, "job_id": "24360", "scores": {"GO": 1}})
    table = tmp_path / "old.tsv"
    table.write_text("old table\n")
    others = [{"key": "other", "label": "other", "scores": {ENDPOINTS[0]: .5}, "cells": []},
        {"key": "canonical_collection:canonical", "label": "canonical",
         "scores": {key: .8 for key in ENDPOINTS}, "cells": [{"dataset": "QfO", "evidence": "unchanged"}]}]
    skipped = {"dataset": "QfO", "status": "not_run_unverified_control", "scores": {}}
    method = {"key": "shared_hit:shared_hit1", "label": "shared_hit1",
        "scores": {ENDPOINTS[0]: .73, ENDPOINTS[-1]: .86}, "cells": [skipped]}
    parent = tmp_path / "parent.json"
    write_json(parent, {"methods": [*others, method], "score_table": record(table),
                        "round_screening": {}, "limitations": []})
    failed = {"dataset": "QfO", "status": "failed", "scores": {}, "producer_job_id": "24247",
              "error": "Native pairs are not canonical, unique and sorted"}
    failures = [skipped, failed]
    monkeypatch.setattr(recovery, "load_plan", lambda *args: ({"runtime": {"files": []}}, [], deepcopy(failures)))
    monkeypatch.setattr(reporter, "runtime", lambda: {"files": []})
    monkeypatch.setattr(recovery.subprocess, "check_output", lambda *args, **kwargs: ACCOUNTING)
    return root, parent, others, method, failures


@pytest.mark.parametrize("status", ["scoring_derived_shared_qfo_relation_set", "failed_shared_qfo_lossless_export_recovery"])
def test_pending_or_failed_scores_not_admitted_with_history_retained(tmp_path, monkeypatch, status):
    root, parent, others, old, failures = panel(tmp_path, monkeypatch, status)
    monkeypatch.setattr(recovery, "cell", lambda *args: pytest.fail("Partial scores admitted"))
    result = reporter.build(tmp_path, parent, tmp_path / "out", shared_export_root=root)
    assert result["methods"][:-1] == others
    method = result["methods"][-1]
    assert method["scores"] == old["scores"]
    assert all(attempt in method["cells"] for attempt in failures)
    assert "NA" in (tmp_path / "out/scores.tsv").read_text()
    assert result["defaults_promoted"] is result["goal_completed"] is False


def test_completed_updates_only_six_qfo_values_and_keeps_control_record(tmp_path, monkeypatch):
    root, parent, others, old, failures = panel(tmp_path, monkeypatch, recovery.SUCCESS)
    scores = {key: .81 for key in ENDPOINTS if key.startswith("QfO:")}
    good = {"dataset": "QfO", "status": "independent_content_readback_verified", "scores": scores,
        "matched_control": {"scores": {key: .8 for key in scores}},
        "products": {name: {"path": name, "bytes": 1, "sha256": "same"} for name in PRODUCTS}}
    monkeypatch.setattr(recovery, "cell", lambda *args: deepcopy(good))
    result = reporter.build(tmp_path, parent, tmp_path / "out", shared_export_root=root)
    assert result["methods"][:-1] == others
    method = result["methods"][-1]
    assert method["scores"] == {**old["scores"], **scores}
    assert all(attempt in method["cells"] for attempt in failures)
    assert (tmp_path / "out/snapshots/shared_export_reference.json").exists()


@pytest.mark.parametrize("change", ["missing", "different", "incomplete", "duplicate"])
def test_existing_control_must_agree_without_being_rewritten(change):
    scores = {key: .8 for key in ENDPOINTS if key.startswith("QfO:")}
    methods = [{"key": "canonical_collection:canonical", "scores": deepcopy(scores)}]
    if change == "missing":
        methods = []
    elif change == "different":
        methods[0]["scores"][next(iter(scores))] = .7
    elif change == "incomplete":
        scores.pop(next(iter(scores)))
    else:
        methods *= 2
    with pytest.raises(ValueError, match="Existing canonical control"):
        recovery.require_existing_control(methods, {"matched_control": {"scores": scores}})


@pytest.fixture
def completed(tmp_path, monkeypatch):
    root = tmp_path / "producer"
    root.mkdir()
    artifact = tmp_path / "native"
    artifact.write_text("native fixture\n")
    pin = record(artifact)
    raw = {key.split(": ", 1)[1]: .81 for key in ENDPOINTS if key.startswith("QfO:")}
    control = {"report": pin, "cell": {"scores": endpoint_vector("QfO", {key: .8 for key in raw})}}
    original_path = tmp_path / "original/results.json"
    original_path.parent.mkdir()
    failed = {"execution": {"original_native_exit_zero_export_failed": True}, "control": control}
    write_json(original_path, failed)
    original = {"output": str(original_path.parent)}
    converter = tmp_path / "executor/qfo_benchmark/og_to_pairwise.py"
    converter.parent.mkdir(parents=True)
    converter.write_text("# trusted converter fixture\n")
    plan = {"source_root": str(converter.parent.parent), "scoring_id": "fixture", "failure": pin}
    native = {"genes": 984137, "species": 78, "groups": 370016, "pairs": recovery.DISTINCT,
        "products": {name: pin for name in PRODUCTS}, "summary": {"native": True}}
    projection = {"status": "full_byte_projection_verified"}
    write_json(root / "execution_plan.json", plan)
    report = {"status": recovery.SUCCESS, "job_id": "24360", "postflight_pins_verified": True,
        "native_inference_rerun": False, "defaults_promoted": False, "goal_completed": False,
        "execution_plan": record(root / "execution_plan.json"), "native_products": native["products"],
        "native_summary": native["summary"], "normalization": pin, "converter": record(converter),
        "coverage": {key: native[key] for key in ("genes", "species", "groups", "pairs")},
        "scoring": {"scores": deepcopy(raw)}, "scores": deepcopy(raw),
        "original_failed_producer": record(original_path),
        "control": {"report": pin, "scores": control["cell"]["scores"]},
        "deltas": {key: .81 - .8 for key in control["cell"]["scores"]}}
    write_json(root / "results.json", report)
    (root / "recovery.time.txt").write_text("\n".join(f"{label}: {value}" for label, value in
        zip(LABELS, ("0:30.00", "20", "2", "1000", "0"))) + "\n")
    monkeypatch.setattr(recovery, "parent_context", lambda *args: (original, {}, failed, [pin]))
    monkeypatch.setattr(recovery, "export_readback", lambda *args: (deepcopy(native), projection, [pin]))
    def reader(repo, dataset, scored, actual, folder, participant):
        assert dataset == "QfO" and actual == native and folder == root
        assert participant == "ohmm_hillclimb_r1_fixture"
        return deepcopy(raw), {}, [pin]
    monkeypatch.setattr(recovery, "readback_scores", reader)
    monkeypatch.setattr(recovery.subprocess, "check_output", lambda *args, **kwargs: ACCOUNTING)
    return root, plan, report


def test_completed_cell_reads_raw_scores_with_separate_native_and_recovery_costs(tmp_path, completed):
    cell = recovery.cell(tmp_path, *completed)
    assert len(cell["scores"]) == 6 and cell["coverage"]["genes"] == 984137
    assert cell["resources"]["completed_native_stage_in_failed_export_job"] == {
        "original_native_exit_zero_export_failed": True}
    assert cell["resources"]["derived_export_validation_and_scoring"]["gnu_time"]["exit_status"] == 0
    assert cell["source_revision"] == recovery.REVISION


@pytest.mark.parametrize("field,value", [("status", "running"), ("postflight_pins_verified", False),
    ("native_inference_rerun", True), ("defaults_promoted", True), ("goal_completed", True),
    ("error", "failed"), ("postflight_error", "changed")])
def test_unverified_execution_refused(tmp_path, completed, field, value):
    root, plan, report = completed
    report[field] = value
    with pytest.raises(ValueError, match="successfully verified"):
        recovery.cell(tmp_path, root, plan, report)


@pytest.mark.parametrize("field", ["native_products", "native_summary", "coverage",
    "original_failed_producer", "control", "converter"])
def test_changed_native_or_control_provenance_refused(tmp_path, completed, field):
    root, plan, report = completed
    report[field] = {}
    with pytest.raises(ValueError, match="authenticated native exports"):
        recovery.cell(tmp_path, root, plan, report)


@pytest.mark.parametrize("field", ["scores", "deltas"])
def test_changed_score_or_delta_refused(tmp_path, completed, field):
    root, plan, report = completed
    report[field][next(iter(report[field]))] = 1
    with pytest.raises(ValueError, match="matched-control deltas"):
        recovery.cell(tmp_path, root, plan, report)


@pytest.mark.parametrize("text,match", [("Exit status: 0\n", "Incomplete GNU-time"),
    ("failed", "Failed stage")])
def test_incomplete_or_failed_resource_record_refused(tmp_path, completed, text, match):
    root, plan, report = completed
    path = root / "recovery.time.txt"
    path.write_text(path.read_text().replace("Exit status: 0", "Exit status: 1") if text == "failed" else text)
    with pytest.raises(ValueError, match=match):
        recovery.cell(tmp_path, root, plan, report)


@pytest.fixture
def export_fixture(tmp_path, monkeypatch):
    root, folder = tmp_path / "recovery", tmp_path / "native"
    view, phylo = root / "canonical_export", folder / "output/orthohmm_phylogeny"
    view.mkdir(parents=True)
    phylo.mkdir(parents=True)
    header = "gene_a\tspecies_a\tgene_b\tspecies_b"
    fields = ("a", "sp1", "b", "sp2")
    row = "\t".join(fields)
    (phylo / PRODUCTS[1]).write_text(header + "\n" + (row + "\n") * 2)
    (phylo / PRODUCTS[2]).write_text(header + "\tconfidence\n" + (row + "\tmedium\n") * 2)
    proof = normalize(phylo / PRODUCTS[1], phylo / PRODUCTS[2], view / PRODUCTS[1], view / PRODUCTS[2], fields, 2, 1)
    copied = []
    for name in (PRODUCTS[0], "provenance_manifest.json", "reconciliation_summary.json", "species_tree.rooted.nwk"):
        (phylo / name).write_text(name + "\n")
        (view / name).write_bytes((phylo / name).read_bytes())
        copied += [record(phylo / name), record(view / name)]
    (folder / "metrics.json").write_text("{}\n")
    manifest = {"status": "explicit_lossless_canonical_export_view_verified", "normalization": proof,
        "independent_projection": verify_projection(proof), "unchanged_copied_products": copied,
        "original_completed_native_metadata": record(folder / "metrics.json"), "native_inference_rerun": False}
    write_json(root / "normalization.json", manifest)
    report = {"normalization": record(root / "normalization.json")}
    native = {"status": "derived_native_content_verified", "genes": 984137, "species": 78,
              "groups": 370016, "pairs": 1, "products": {}, "summary": {}}
    monkeypatch.setattr(recovery, "DUPLICATE", fields)
    monkeypatch.setattr(recovery, "DUPLICATE_LINE", 3)
    monkeypatch.setattr(recovery, "ROWS", 2)
    monkeypatch.setattr(recovery, "DISTINCT", 1)
    def native_reader(case, parent, pair_export):
        assert pair_export == proof
        return deepcopy(native)
    monkeypatch.setattr(recovery, "shared_native", native_reader)
    return root, report, {"case": {"output": str(folder)}}, manifest, native


def test_export_readback_requires_real_byte_projection_and_unchanged_metadata(export_fixture):
    root, report, original, manifest, _ = export_fixture
    native, projection, pins = recovery.export_readback(root, report, original, {})
    assert projection == manifest["independent_projection"]
    assert len(native["products"]) == 4 and native["phylogeny"] == str(root / "canonical_export")
    assert all(pin["sha256"] for pin in pins)


@pytest.mark.parametrize("field", ["status", "native_inference_rerun", "original_completed_native_metadata",
    "independent_projection", "unchanged_copied_products", "duplicate", "duplicate_confidence",
    "original_rows", "distinct_rows", "line", "original", "derived"])
def test_changed_normalization_manifest_refused(export_fixture, field):
    root, report, original, manifest, _ = export_fixture
    if field in manifest:
        manifest[field] = "changed"
    elif field == "line":
        manifest["normalization"]["skips"]["pairs"]["line"] = 4
    else:
        manifest["normalization"][field] = "changed" if field not in ("original", "derived") else {}
    write_json(root / "normalization.json", manifest)
    report["normalization"] = record(root / "normalization.json")
    with pytest.raises(ValueError):
        recovery.export_readback(root, report, original, {})


def test_changed_copied_native_metadata_refused(export_fixture):
    root, report, original, _, _ = export_fixture
    (root / "canonical_export/provenance_manifest.json").write_text("different\n")
    with pytest.raises(ValueError, match="metadata or root groups"):
        recovery.export_readback(root, report, original, {})


@pytest.mark.parametrize("field,value", [("status", "fresh_native_content_verified"),
    ("genes", 1), ("species", 1), ("groups", 1), ("pairs", 2)])
def test_incomplete_native_coverage_refused(export_fixture, field, value):
    root, report, original, _, native = export_fixture
    native[field] = value
    with pytest.raises(ValueError, match="complete original native coverage"):
        recovery.export_readback(root, report, original, {})


@pytest.mark.parametrize("field,value", [("job_id", "other"), ("source_revision", "other"),
    ("scientific_revision", "other")])
def test_changed_launch_identity_refused_before_parent_readback(tmp_path, monkeypatch, field, value):
    root = tmp_path / "producer"
    root.mkdir()
    source = tmp_path / "executor"
    plan = {"repo": str(tmp_path), "output": str(root), "source_root": str(source),
        "source_revision": "executor", "scientific_revision": recovery.REVISION,
        "native_inference_rerun": False, "scoring_id": "fixture"}
    path = root / "execution_plan.json"
    write_json(path, plan)
    launch = {"job_id": "24360", "execution_plan": record(path), "source_root": str(source),
        "source_revision": "executor", "scientific_revision": recovery.REVISION}
    launch_path = tmp_path / recovery.SUBMISSION
    launch_path.parent.mkdir(parents=True)
    write_json(launch_path, launch)
    report = {"job_id": "24360", "execution_plan": record(path), "source_revision": "executor",
              "scientific_revision": recovery.REVISION}
    report[field] = value
    monkeypatch.setattr(recovery, "short_id", lambda *args: "fixture")
    monkeypatch.setattr(recovery, "parent_context", lambda *args: pytest.fail("Wrong launch admitted"))
    with pytest.raises(ValueError):
        recovery.load_plan(tmp_path, root, report)


@pytest.fixture
def planned(tmp_path, monkeypatch):
    root, parent_root = tmp_path / "recovery", tmp_path / "original"
    root.mkdir()
    parent_root.mkdir()
    artifact = tmp_path / "artifact.json"
    write_json(artifact, {})
    pin = record(artifact)
    plan = {"repo": str(tmp_path), "output": str(root), "source_root": str(tmp_path / "executor"),
        "source_revision": "executor", "scientific_revision": recovery.REVISION,
        "native_inference_rerun": False, "scoring_id": "fixture", "failure": pin,
        "sources": [pin], "runtime": {"files": [pin]}}
    write_json(root / "execution_plan.json", plan)
    plan_pin = record(root / "execution_plan.json")
    launch = {"job_id": "24360", "execution_plan": plan_pin, "source_root": plan["source_root"],
        "source_revision": "executor", "scientific_revision": recovery.REVISION}
    launch_path = tmp_path / recovery.SUBMISSION
    launch_path.parent.mkdir(parents=True)
    write_json(launch_path, launch)
    report = {"job_id": "24360", "execution_plan": plan_pin, "source_revision": "executor",
              "scientific_revision": recovery.REVISION}
    skip = parent_root / "original_skip_snapshot.json"
    write_json(skip, {"job_id": "24242"})
    failed = {"original_skip": record(skip), "status": "failed_or_unverified_shared_qfo_continuation",
        "error": "Native pairs are not canonical, unique and sorted", "error_type": "ValueError",
        "execution": {"returncode": 0}}
    write_json(parent_root / "results.json", failed)
    original = {"output": str(parent_root), "runtime": deepcopy(plan["runtime"]), "parent_plan": pin}
    monkeypatch.setattr(recovery, "short_id", lambda *args: "fixture")
    monkeypatch.setattr(recovery, "parent_context", lambda *args: (original, {}, failed, [pin]))
    monkeypatch.setattr(recovery, "runtime", lambda: deepcopy(plan["runtime"]))
    def helpers(source, context, revision, *, protocol, slurm):
        assert context == original and revision == "executor"
        assert protocol == recovery.PROTOCOL and slurm == recovery.SLURM
        return [pin]
    monkeypatch.setattr(recovery, "verify_helpers", helpers)
    def skipped(archived, parent, identity):
        assert archived == {"job_id": "24242"} and parent == {} and identity == pin
        return {"dataset": "QfO", "status": "not_run_unverified_control", "scores": {}}
    monkeypatch.setattr(recovery, "original_skip", skipped)
    return root, report, original, failed


def test_load_plan_preserves_skip_and_terminal_failed_native_without_copying_huge_report(tmp_path, planned):
    root, report, _, failed = planned
    plan, pins, attempts = recovery.load_plan(tmp_path, root, report)
    assert attempts[0]["status"] == "not_run_unverified_control"
    assert attempts[0]["producer_job_id"] == "24242"
    assert attempts[1]["status"] == "failed" and attempts[1]["producer_job_id"] == "24247"
    assert attempts[1]["resources"] == failed["execution"]
    assert attempts[1]["source_snapshot"]["sha256"]
    assert all(attempt["scores"] == {} for attempt in attempts)
    assert plan["native_inference_rerun"] is False and pins


@pytest.mark.parametrize("change", ["executor", "runtime", "parent_runtime", "skip_location"])
def test_plan_provenance_must_be_reverified(tmp_path, planned, monkeypatch, change):
    root, report, original, failed = planned
    if change == "executor":
        monkeypatch.setattr(recovery, "verify_helpers", lambda *args, **kwargs: [])
    elif change == "runtime":
        monkeypatch.setattr(recovery, "runtime", lambda: {})
    elif change == "parent_runtime":
        original["runtime"] = {}
    else:
        path = tmp_path / "other_skip.json"
        write_json(path, {"job_id": "24242"})
        failed["original_skip"] = record(path)
    with pytest.raises(ValueError):
        recovery.load_plan(tmp_path, root, report)
