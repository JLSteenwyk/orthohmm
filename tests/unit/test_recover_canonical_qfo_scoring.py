from copy import deepcopy
import json
from pathlib import Path
from tempfile import TemporaryDirectory

import pytest

from benchmark_tools.recover_canonical_qfo_scoring import PATH_ERROR, REVISION, failed_row, scoring_id
from benchmark_tools import recover_canonical_qfo_scoring as recovery


def parent():
    return {"job_id": "24240", "source_revision": REVISION, "status": "running", "arms": [
        {"dataset": "QfO", "label": "canonical", "status": "failed", "scores": {},
         "error_type": "ValueError", "error": PATH_ERROR,
         "graph_execution": {"returncode": 0}, "native_execution": {"returncode": 0},
         "outputs": [{"path": "native"}], "native_summary": {"genes": 984137}},
        {"dataset": "ThreeKingdoms", "label": "arrival", "status": "running"}]}


def test_selected_failure_is_stable_while_unrelated_arms_progress():
    report = parent()
    row = failed_row(report)
    report["arms"][1]["status"] = "independent_content_readback_verified"
    assert failed_row(report) == row
    row["status"] = "changed"
    assert report["arms"][0]["status"] == "failed"


@pytest.mark.parametrize("change", [
    {"status": "running"}, {"error": "Segmentation fault"}, {"error_type": "RuntimeError"},
    {"scores": {"GO": .5}}, {"native_execution": {"returncode": -11}},
    {"graph_execution": {"returncode": 1}}, {"outputs": []}, {"native_summary": {}}])
def test_reject_inference_failure_or_other_scoring_failure(change):
    report = parent()
    report["arms"][0].update(change)
    with pytest.raises(ValueError, match="successful native"):
        failed_row(report)


@pytest.mark.parametrize("change", [{"job_id": "other"}, {"source_revision": "other"}])
def test_original_producer_identity_is_required(change):
    report = parent()
    report.update(change)
    with pytest.raises(ValueError, match="unique original"):
        failed_row(report)


def test_ambiguous_parent_rejected():
    report = parent()
    report["arms"].append(deepcopy(report["arms"][0]))
    with pytest.raises(ValueError, match="unique original"):
        failed_row(report)


def test_short_unique_id_and_existing_namespace(tmp_path):
    output = tmp_path / "output"
    with TemporaryDirectory(dir="/tmp", prefix="cq-") as temporary:
        repo = Path(temporary)
        identifier = scoring_id(repo, output)
        assert len(identifier) == 15
        assert scoring_id(repo, output) == identifier
        assert scoring_id(repo, output / "new") != identifier
        (repo / "qfo_benchmark/w" / ("hc1_" + identifier)).mkdir(parents=True)
        with pytest.raises(FileExistsError):
            scoring_id(repo, output)


def test_overlong_host_path_rejected():
    with pytest.raises(ValueError, match="validated length"):
        scoring_id(Path("/" + "x" * 160), Path("/output"))


@pytest.fixture
def validation_case(tmp_path, monkeypatch):
    report = parent()
    row = report["arms"][0]
    folder = tmp_path / "QfO/canonical"
    folder.mkdir(parents=True)
    artifact = tmp_path / "artifact"
    artifact.write_text("unchanged")
    pin = recovery.record(artifact)
    row.update(ordering={"checkpoint": str(tmp_path / "checkpoint"), "manifest": pin},
               candidate={"membership_constraints": pin})
    expected = (["graph"], ["native"])
    row["graph_execution"]["command"] = expected[0]
    row["native_execution"]["command"] = expected[1]
    preflight = {"graph": expected[0], "native": expected[1]}
    (folder / "preflight.json").write_text(json.dumps(preflight))
    plan = {"source_root": str(tmp_path / "source"), "source_revision": REVISION,
            "cpu": 16, "cases": [{"dataset": "QfO"}]}
    native = {"genes": 984137, "species": 78, "summary": row["native_summary"], "checked": [pin]}
    monkeypatch.setattr(recovery, "check_plan", lambda *args: (plan, []))
    monkeypatch.setattr(recovery, "commands", lambda *args: expected)
    monkeypatch.setattr(recovery, "replay_native", lambda *args: native)

    def reject_original_path(work):
        raise ValueError(PATH_ERROR)

    monkeypatch.setattr(recovery, "require_darwin_work_path", reject_original_path)
    return tmp_path, report, row, plan, native, artifact


def test_native_readback_and_artifact_hashes_required(validation_case):
    root, report, row, _, native, artifact = validation_case
    result, pins = recovery.verify_parent(root, root, report, deepcopy(row))
    assert result == native
    assert pins
    artifact.write_text("modified")
    with pytest.raises(ValueError):
        recovery.verify_parent(root, root, report, deepcopy(row))


@pytest.mark.parametrize("field,value", [("genes", 100), ("species", 12), ("summary", {"changed": True})])
def test_native_universe_and_summary_required(validation_case, field, value):
    root, report, row, _, native, _ = validation_case
    native[field] = value
    with pytest.raises(ValueError, match="universe or summary"):
        recovery.verify_parent(root, root, report, deepcopy(row))


def test_actual_native_command_must_match_plan(validation_case):
    root, report, row, _, _, _ = validation_case
    row["native_execution"]["command"] = ["different"]
    with pytest.raises(ValueError, match="Executed commands"):
        recovery.verify_parent(root, root, report, deepcopy(row))


def test_preflight_must_match_frozen_command(validation_case):
    root, report, row, _, _, _ = validation_case
    (root / "QfO/canonical/preflight.json").write_text(json.dumps({"graph": ["different"], "native": ["native"]}))
    with pytest.raises(ValueError, match="frozen plan"):
        recovery.verify_parent(root, root, report, deepcopy(row))


def test_failure_path_must_reproduce(validation_case, monkeypatch):
    root, report, row, _, _, _ = validation_case
    monkeypatch.setattr(recovery, "require_darwin_work_path", lambda *args: None)
    with pytest.raises(ValueError, match="does not reproduce"):
        recovery.verify_parent(root, root, report, deepcopy(row))
