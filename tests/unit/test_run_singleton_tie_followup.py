from copy import deepcopy
import importlib.util
import json
from pathlib import Path
import subprocess
import sys
from types import SimpleNamespace

import numpy as np
import pytest

from benchmark_tools.prepare_ob_candidate_neighborhood import check, record
from benchmark_tools.run_membership_hillclimb import PRODUCTS, write_json

PATH = Path(__file__).resolve().parents[2] / "benchmark_tools/run_singleton_tie_followup.py"
spec = importlib.util.spec_from_file_location("singleton_followup_driver", PATH)
runner = importlib.util.module_from_spec(spec)
spec.loader.exec_module(runner)


def prepared():
    return {"status": "prespecified_followup_inputs_prepared_not_submitted", "scientific_revision": runner.REVISION,
            "cpu": 16, "pilot_gate": {"eligible": True, "control_f1": .74, "variant_f1": .736,
                                      "delta": .736 - .74},
            "cases": [{"dataset": dataset} for dataset in runner.UNIVERSES]}


@pytest.mark.parametrize("dataset", list(runner.UNIVERSES))
def test_only_prespecified_dataset_case_selected(dataset):
    plan = prepared()
    assert runner.selected_case(plan, dataset) == {"dataset": dataset}


@pytest.mark.parametrize("key,value", [("cpu", 32), ("scientific_revision", "new"), ("status", "running")])
def test_changed_fixed_contrast_rejected(key, value):
    plan = prepared()
    plan[key] = value
    with pytest.raises(ValueError, match="unchanged prespecified"):
        runner.selected_case(plan, "QfO")


@pytest.mark.parametrize("change", [{"delta": float("nan")}, {"delta": 0.1}, {"eligible": False},
                                    {"variant_f1": float("inf")}, {"variant_f1": -.1}])
def test_pilot_gate_cannot_be_forged(change):
    plan = prepared()
    plan["pilot_gate"].update(change)
    with pytest.raises(ValueError):
        runner.selected_case(plan, "QfO")


def test_unsupported_dataset_and_reordered_inventory_rejected():
    with pytest.raises(ValueError, match="Unsupported"):
        runner.selected_case(prepared(), "OrthoBench")
    plan = prepared()
    plan["cases"].reverse()
    with pytest.raises(ValueError, match="unchanged prespecified"):
        runner.selected_case(plan, "QfO")


def test_scoring_namespace_is_short_stable_and_distinct(tmp_path):
    assert len(runner.scoring_id(tmp_path)) == 15
    assert runner.scoring_id(tmp_path) == runner.scoring_id(tmp_path)
    assert runner.scoring_id(tmp_path) != runner.scoring_id(tmp_path / "other")


def test_mixed_scientific_imports_rejected(tmp_path):
    with pytest.raises(ValueError, match="already imported"):
        runner.science_api(tmp_path)


def test_standalone_driver_loads_existing_frozen_science_in_fresh_process():
    repo = PATH.parent.parent
    source = repo / "benchmarks/work/singleton_tie_orthobench_r1_v1_20261014"
    if not source.exists():
        pytest.skip("Retained local scientific snapshot unavailable")
    code = """
import importlib.util, inspect, json, pathlib, sys
spec = importlib.util.spec_from_file_location('driver', sys.argv[1])
driver = importlib.util.module_from_spec(spec)
spec.loader.exec_module(driver)
root = pathlib.Path(sys.argv[2]).resolve()
api = driver.science_api(root)
paths = [inspect.getfile(getattr(api, name)) for name in
         ('commands', 'prepare_partition', 'load_accuracy_checkpoint', '_expand_phylogeny_candidates', 'readback_scores')]
assert all(pathlib.Path(path).resolve().is_relative_to(root) for path in paths)
print(json.dumps({'all_scientific_functions_from_frozen_source': True}))
"""
    result = subprocess.run([sys.executable, "-c", code, str(PATH), str(source)],
                            capture_output=True, text=True, timeout=60)
    assert result.returncode == 0, result.stderr
    assert json.loads(result.stdout)["all_scientific_functions_from_frozen_source"]


def worker_fixture(tmp_path, monkeypatch, graph_failure=False, wrong_universe=False):
    output = tmp_path / "output"
    output.mkdir()
    source = tmp_path / "science"
    driver = PATH.parent.parent
    artifact = tmp_path / "candidate"
    artifact.write_text("candidate")
    candidate = {key: record(artifact) for key in ("candidate_partition", "membership_constraints")}
    case = {"dataset": "ThreeKingdoms", "checkpoint": str(tmp_path / "checkpoint"), "checkpoint_sha": "fixed",
            "case": {}, "candidate": candidate, "golden": {},
            "control": {"scores": {"ThreeKingdoms: Three Kingdoms F1": .872}}}
    initial = prepared()
    initial.update(scientific_source=str(source), runtime={"files": []}, tool_path="fixed")
    initial["cases"][1] = case
    initial_path = tmp_path / "prepared.json"
    write_json(initial_path, initial)
    plan = {"repo": str(tmp_path), "driver_root": str(driver), "scientific_source": str(source),
            "output": str(output), "inputs": [], "driver_sources": [], "prepared_inputs": record(initial_path),
            "dataset": "ThreeKingdoms", "case": case, "scientific_revision": runner.REVISION,
            "runtime": initial["runtime"], "tool_path": "fixed", "cpu": 16,
            "scoring_id": runner.scoring_id(output), "driver_revision": "fixed", "scope": "fixture"}
    path = output / "execution_plan.json"
    write_json(path, plan)
    calls = []
    def stage(command, folder, name, actual_source, env):
        assert actual_source == source and env["PYTHONPATH"] == str(source)
        assert env["NUMBA_NUM_THREADS"] == "16"
        snapshot = json.loads((output / "results.json").read_text())
        assert snapshot["arms"][-1]["status"] == name + "_running"
        calls.append((folder.name, name))
        if name == "graph":
            write_json(folder / "graph_metrics.json", {"parameters": {"singleton_tie_policy": folder.name,
                "cpm_resolution": .1, "leiden_seed": 4, "profile_iterations": 1,
                "profile_min_species": 1, "matrix": "BLOSUM62"}})
        else:
            phylo = folder / "native/orthohmm_phylogeny"
            phylo.mkdir(parents=True)
            for name in PRODUCTS:
                (phylo / name).write_text("native fixture")
        return {"returncode": int(graph_failure and name == "graph"), "command": command}
    api = SimpleNamespace(check=check, record=record, write_json=write_json, PRODUCTS=PRODUCTS,
        runtime=lambda: initial["runtime"], host_activity=lambda: {},
        load_accuracy_checkpoint=lambda *a, **k: (["x"] * 443217, np.arange(443217) % 12, [], [], []),
        commands=lambda *a: (["graph"], ["native"]), timed_stage=stage,
        prepare_partition=lambda *a: deepcopy(candidate), _expand_phylogeny_candidates=None,
        load_membership_constraints=None, validate_partition=lambda *a: {},
        replay_native=lambda *a: {"genes": 443216 if wrong_universe else 443217, "species": 12,
                                 "prediction": str(artifact), "checked": [], "summary": {}},
        compare_golden=lambda *a: {},
        score_three_kingdoms=lambda repo, prediction, folder: {"scores": {
            "Three Kingdoms F1": .872 if folder.name == "first_hit" else .875}},
        readback_scores=lambda repo, dataset, row, *a: (row["scores"], {}, []),
        endpoint_vector=lambda dataset, scores: {dataset + ": " + key: value for key, value in scores.items()})
    monkeypatch.setattr(runner, "science_api", lambda *a: api)
    monkeypatch.setenv("CUDA_VISIBLE_DEVICES", "")
    return path, calls, output


def test_failed_control_graph_never_runs_variant(tmp_path, monkeypatch):
    path, calls, output = worker_fixture(tmp_path, monkeypatch, graph_failure=True)
    result = runner.run(path)
    assert calls == [("first_hit", "graph")]
    assert [row["status"] for row in result["arms"]] == ["failed", "not_run_failed_control_gate"]
    assert not (output / "defer_ambiguous").exists()


def test_wrong_control_universe_never_runs_variant(tmp_path, monkeypatch):
    path, calls, output = worker_fixture(tmp_path, monkeypatch, wrong_universe=True)
    result = runner.run(path)
    assert result["arms"][0]["error"] == "Native complete gene universe differs"
    assert len(calls) == 2 and not (output / "defer_ambiguous").exists()


def test_complete_matched_worker_records_stage_progress_and_separate_scores(tmp_path, monkeypatch):
    path, calls, _ = worker_fixture(tmp_path, monkeypatch)
    result = runner.run(path)
    assert result["status"] == "both_singleton_policies_scored"
    assert calls == [("first_hit", "graph"), ("first_hit", "native"),
                     ("defer_ambiguous", "graph"), ("defer_ambiguous", "native")]
    assert result["deltas"]["ThreeKingdoms: Three Kingdoms F1"] == pytest.approx(.003)
    assert result["defaults_promoted"] is False and result["goal_completed"] is False
    with pytest.raises(ValueError, match="preserve previous"):
        runner.run(path)


def test_execution_cannot_retarget_scientific_source(tmp_path, monkeypatch):
    path, calls, _ = worker_fixture(tmp_path, monkeypatch)
    plan = json.loads(path.read_text())
    plan["scientific_source"] = str(tmp_path / "another")
    write_json(path, plan)
    with pytest.raises(ValueError, match="fixed scientific contrast"):
        runner.run(path)
    assert calls == []


def test_reporting_evidence_frozen_without_retargeting_scientific_scorers(tmp_path):
    repo, driver = tmp_path / "repo", tmp_path / "driver"
    path = repo / "benchmark_tools/reporter.py"
    path.parent.mkdir(parents=True)
    path.write_text("admitted reporting source")
    scorer = repo / "three_kingdoms/score.py"
    scorer.parent.mkdir()
    scorer.write_text("scientific scorer")
    api = SimpleNamespace(check=check, record=record)
    frozen, origins = runner.freeze_reporting_evidence(api, [record(path), record(scorer)], repo, driver)
    assert frozen[0]["path"] == str(driver / "evidence/benchmark_tools/reporter.py")
    assert frozen[1] == record(scorer)
    assert origins[0]["original"] == record(path)
    path.write_text("later unrelated reporting change")
    for pin in frozen:
        check(pin)
    scorer.write_text("scientific change")
    with pytest.raises(ValueError):
        check(frozen[1])


def test_existing_reporting_capture_must_match_the_admitted_version(tmp_path):
    repo, driver = tmp_path / "repo", tmp_path / "driver"
    path = repo / "benchmark_tools/reporter.py"
    path.parent.mkdir(parents=True)
    path.write_text("admitted")
    target = driver / "evidence/benchmark_tools/reporter.py"
    target.parent.mkdir(parents=True)
    target.write_text("other")
    with pytest.raises(ValueError, match="Frozen reporting evidence"):
        runner.freeze_reporting_evidence(SimpleNamespace(check=check, record=record), [record(path)], repo, driver)


def test_qfo_control_reuses_pinned_assessment_without_resampling_fas(tmp_path):
    recovery = tmp_path / "benchmarks/results/canonical_qfo_scoring_recovery_r1_v1_20261014"
    recovery.mkdir(parents=True)
    write_json(recovery / "execution_plan.json", {"scoring_id": "recovered"})
    path = recovery / "results.json"
    write_json(path, {"status": "canonical_qfo_scoring_recovered_independently_verified",
        "execution_plan": record(recovery / "execution_plan.json"),
        "scores": {"FAS": .7637247098243404}, "conversion": {}, "metric_files": []})
    calls = []
    def readback(repo, dataset, row, native, folder, participant):
        calls.append((folder, participant))
        return row["scores"], {}, []
    api = SimpleNamespace(record=record, readback_scores=readback,
        endpoint_vector=lambda dataset, scores: {dataset + ": " + key: value for key, value in scores.items()},
        score_qfo=lambda *a: pytest.fail("Control FAS resampled"))
    plan = {"repo": str(tmp_path), "dataset": "QfO", "inputs": [record(path), record(recovery / "execution_plan.json")],
            "case": {"control": {"scores": {"QfO: FAS": .7637247098243404}}}}
    row = {}
    scores, _, _ = runner.control_scores(api, plan, row, {}, tmp_path / "control")
    assert scores["FAS"] == .7637247098243404
    assert calls == [(recovery, "ohmm_hillclimb_r1_recovered")]
    assert row["control_scoring_reuse"] == record(path)
    plan["inputs"].pop()
    with pytest.raises(ValueError, match="pinned producer"):
        runner.control_scores(api, plan, {}, {}, tmp_path / "control")
