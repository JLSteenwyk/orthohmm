from copy import deepcopy
import json
from pathlib import Path

import pytest

from benchmark_tools import report_local_joint_pilot as reader
from benchmark_tools import run_local_joint_pilot as pilot
from benchmark_tools.prepare_ob_candidate_neighborhood import record
from orthohmm.accuracy import load_accuracy_checkpoint
from tests.unit.test_local_rescue_readback import accepted
from tests.unit.test_run_local_joint_pilot import prepared
from tests.unit.test_profile_evidence_readback import family_cells


@pytest.fixture
def finished(prepared, monkeypatch):
    path, plan, _ = prepared
    report = pilot.run(path)
    monkeypatch.setattr(reader, "UNIVERSE", (16, 4))
    monkeypatch.setattr(reader, "readback_scores", lambda *a:
        ({"OrthoBench F1": .7}, {"precision": .8, "recall": .6}, []))
    return plan, report


def test_actual_full_proposal_recheck_and_final_native_admission(finished, tmp_path, monkeypatch):
    plan, report = finished
    before = [record(p["path"]) for p in plan["inputs"]]
    def forbidden(*a, **kw):
        pytest.fail("Post-run readback must reuse source trees, not infer replacements")
    monkeypatch.setattr("orthohmm.phylogeny_pipeline._run_to_file", forbidden)
    original, repeated, pins = reader.recheck_proposals(plan, report, tmp_path / "recheck")
    assert repeated["ownership"] == original["ownership"]
    assert set(repeated["ownership"]["classes"]) == {"0", "1", "2", "3"}
    assert repeated["contexts"] == original["contexts"]
    names, *_ = load_accuracy_checkpoint(Path(plan["case"]["checkpoint"]), verify=True)
    off, local = [reader.admit_arm(Path(plan["repo"]), plan, row, original, repeated, names)
                  for row in report["arms"]]
    assert off["scores"] == local["scores"] == {"OrthoBench: OrthoBench F1": .7}
    assert off["coverage"]["genes"] == local["coverage"]["genes"] == 16
    assert local["accepted_contexts"] == ["LocalContext000001"]
    assert local["native_pair_turnover"]["added_pairs"] > 0
    assert local["final_effects"]["genes_with_changed_root_og_membership"] > 0
    assert len(pins) > 10 and all(record(p["path"]) == p for p in pins)
    assert before == [record(p["path"]) for p in plan["inputs"]]


@pytest.mark.parametrize("kind", ["contexts", "ownership", "effective_partition", "nomination", "candidate", "cache"])
def test_complete_proposal_comparison_rejects_changed_producer_receipts(finished, tmp_path, kind):
    plan, report = finished
    path = Path(report["proposal_readback"]["path"])
    original = json.loads(path.read_text())
    if kind == "contexts": original["contexts"][0]["accepted"] = False
    elif kind == "ownership": original["ownership"]["classes"]["0"] = "other.faa"
    elif kind == "effective_partition": original["effective_partition"]["unexpected"] = 1
    elif kind == "nomination": original["nomination"]["sha256"] = "wrong"
    elif kind == "candidate": original["candidate"]["sha256"] = "wrong"
    else: original["composite_cache"]["accepted_contexts"] = []
    path.write_text(json.dumps(original))
    report["proposal_readback"] = record(path)
    with pytest.raises(ValueError, match="does not reproduce"):
        reader.recheck_proposals(plan, report, tmp_path / "recheck")


@pytest.mark.parametrize("kind", ["coverage", "summary", "score", "details", "metric", "command", "golden",
                                  "accepted", "turnover", "effects", "species_tree", "candidate"])
def test_final_native_admission_rejects_changed_results(finished, kind):
    plan, report = finished
    original = json.loads(Path(report["proposal_readback"]["path"]).read_text())
    row = deepcopy(report["arms"][0 if kind == "golden" else 1])
    if kind == "coverage": row["coverage"]["all_candidate_genes_preserved"] = False
    elif kind == "summary": row["native_summary"]["unexpected"] = 1
    elif kind == "score": row["scores"]["OrthoBench F1"] += .01
    elif kind == "details": row["details"]["recall"] += .01
    elif kind == "metric": row["metric_vector"]["OrthoBench: OrthoBench F1"] += .01
    elif kind == "command": row["native_execution"]["command"].append("--unexpected")
    elif kind == "golden": row["golden_matches"][next(iter(row["golden_matches"]))] = False
    elif kind == "accepted": row["accepted_contexts"] = []
    elif kind == "turnover": row["native_pair_turnover"]["added_pairs"] += 1
    elif kind == "effects": row["final_effects"]["genes_with_changed_root_og_membership"] += 1
    elif kind == "species_tree": row["final_species_tree_changed"] = not row["final_species_tree_changed"]
    else: row["candidate"]["candidate_partition"]["sha256"] = "wrong"
    names, *_ = load_accuracy_checkpoint(Path(plan["case"]["checkpoint"]), verify=True)
    with pytest.raises(ValueError):
        reader.admit_arm(Path(plan["repo"]), plan, row, original, original, names)


@pytest.mark.parametrize("kind", ["duplicate", "escape", "copy", "source", "status"])
def test_cache_signature_checks_actual_content_and_ownership(tmp_path, kind):
    root = tmp_path / "cache"
    root.mkdir()
    source, target = tmp_path / "source", root / "tree"
    source.write_text("source\n")
    target.write_bytes(source.read_bytes())
    cache = {"status": "authenticated_composite_cache_built", "path": str(root),
             "files": [{"source": record(source), "copy": record(target)}]}
    assert reader.cache_signature(cache)["tree"]["source"] == record(source)
    if kind == "duplicate": cache["files"] *= 2
    elif kind == "escape": cache["files"][0]["copy"] = record(source)
    elif kind == "copy": target.write_text("changed\n")
    elif kind == "source": source.write_text("changed\n")
    else: cache["status"] = "running"
    with pytest.raises(ValueError): reader.cache_signature(cache)


def completed_panel():
    plan = {"cpu": 16, "source_revision": "frozen", "status": "prepared_local_joint_readback_and_native_unrun"}
    report = {"status": "finished_local_readback_and_native_pilot_verified", "job_id": "24495",
        "source_revision": "frozen", "defaults_promoted": False, "goal_completed": False,
        "arms": [{"label": label, "dataset": "OrthoBench", "status": "independent_content_readback_verified"}
                 for label in reader.LABELS]}
    return plan, report


@pytest.mark.parametrize("kind", ["valid", "partial", "reordered", "failed", "wrong_job", "wrong_source", "cpu", "dataset", "defaults"])
def test_complete_successful_ordered_panel_is_required(kind):
    plan, report = completed_panel()
    if kind == "partial": report["arms"].pop()
    elif kind == "reordered": report["arms"].reverse()
    elif kind == "failed": report["arms"][1]["error"] = "failed"
    elif kind == "wrong_job": report["job_id"] = "other"
    elif kind == "wrong_source": report["source_revision"] = "other"
    elif kind == "cpu": plan["cpu"] = 1
    elif kind == "dataset": report["arms"][1]["dataset"] = "other"
    elif kind == "defaults": report["defaults_promoted"] = True
    if kind == "valid": reader.require_completed(report, plan, "24495")
    else:
        with pytest.raises(ValueError, match="complete successful ordered"):
            reader.require_completed(report, plan, "24495")


def test_launch_uses_immutable_commit_and_normalizes_pin(tmp_path, monkeypatch):
    data = {"local_joint_final_pilot": {"execution_plan": {"path": "results/plan.json", "bytes": 12,
        "sha256": "abc", "input_count": 999}}}
    def git(command, **kwargs):
        assert command == ["git", "show", reader.LAUNCH + ":" + reader.LEDGER]
        return json.dumps(data).encode()
    monkeypatch.setattr(reader.subprocess, "check_output", git)
    entry, pin, binding = reader.launch(tmp_path)
    assert pin == {"path": str(tmp_path / "results/plan.json"), "bytes": 12, "sha256": "abc"}
    assert entry == data["local_joint_final_pilot"] and binding["revision"] == reader.LAUNCH


def test_existing_postrun_namespace_is_preserved(tmp_path):
    with pytest.raises(FileExistsError, match="Preserve every"):
        reader.run(tmp_path, tmp_path)


def test_actual_paired_family_uncertainty_uses_single_three_metric_contrast():
    cells = family_cells()[:2]
    for label, cell in zip(reader.LABELS, cells): cell["label"] = label
    result = reader.local_uncertainty(cells)
    assert result["replicates"] == 100000 and result["seed"] == 20261018
    assert set(result["comparisons"]) == {"local_joint_v1"}
    assert "3 planned endpoints" in result["multiplicity"]
    cells.reverse()
    with pytest.raises(ValueError, match="uncertainty protocol or contrast"):
        reader.local_uncertainty(cells)
