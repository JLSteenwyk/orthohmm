from copy import deepcopy
import csv
import hashlib
import json

import pytest

from benchmark_tools.assemble_controlled_fragment_results import terminal_job, prediction_scores, fragment_rows, metric_means, tables, assemble
from benchmark_tools.controlled_fragment_observations import SEEDS, METHODS, METRICS, paired_intervals
from benchmark_tools.prepare_controlled_fragment_observations import record


def test_only_unique_terminal_job_is_scored():
    header = "JobIDRaw|State|ExitCode|Elapsed\n"
    success = "1|COMPLETED|0:0|00:01\n"
    assert terminal_job(header + success, 1)["State"] == "COMPLETED"
    assert terminal_job(header + "1|OUT_OF_MEMORY|0:9|00:01\n", 1)["State"] == "OUT_OF_MEMORY"
    for text in (header, header + success + success, header + "1|RUNNING|0:0|00:01\n", header + "1|COMPLETED|1:0|00:01\n"):
        with pytest.raises(ValueError):
            terminal_job(text, 1)


def test_native_failure_is_retained_without_a_score(monkeypatch):
    monkeypatch.setattr("benchmark_tools.assemble_controlled_fragment_results.admit_method", lambda *a: {"status": "failed", "failure_stage": "native_output", "reason": "native stopped"})
    result = prediction_scores(METHODS[0], {}, {}, {}, {}, {}, [], {}, {})
    assert result["status"] == "failed" and "score" not in result and result["reason"] == "native stopped"


def test_prediction_conversion_requires_native_inventory_and_counts_all_false_positives(tmp_path, monkeypatch):
    artifact = tmp_path / "pairs.txt"
    artifact.write_text("a d\n")
    dataset = {"methods": {METHODS[0]: {"output": str(tmp_path)}}}
    execution = {"methods": {METHODS[0]: {"outputs": [record(artifact)]}}}
    monkeypatch.setattr("benchmark_tools.assemble_controlled_fragment_results.admit_method", lambda *a: {"status": "admitted"})
    monkeypatch.setattr("benchmark_tools.assemble_controlled_fragment_results.load_predictions", lambda *a: ([("a", "d")], [artifact]))
    owners = {"a": "A", "b": "B", "c": "A", "d": "B"}
    result = prediction_scores(METHODS[0], dataset, execution, {}, {}, owners, ["A", "B"], {"ortholog_pairs": [["a", "b"], ["c", "d"]]}, {"a": True, "b": False, "c": False, "d": True})
    assert (result["score"]["tp"], result["score"]["fp"], result["score"]["fn"]) == (0, 1, 2)
    assert result["strata"][2]["score"]["fp"] == 1
    execution["methods"][METHODS[0]]["outputs"] = []
    with pytest.raises(ValueError, match="native execution inventory"):
        prediction_scores(METHODS[0], dataset, execution, {}, {}, owners, [], {"ortholog_pairs": []}, {g: False for g in owners})


def test_unreached_or_preflight_failure_keeps_all_four_outcomes_missing():
    dataset = {"seed": SEEDS[0]}
    for panel_row in (None, {"status": "failed", "reason": "preflight failed"}):
        rows = fragment_rows(dataset, panel_row, {}, {})
        assert len(rows) == 4 and {r["method"] for r in rows} == set(METHODS)
        assert all(r["status"] == "unavailable" and "score" not in r for r in rows)


def test_execution_path_cannot_belong_to_another_dataset(tmp_path):
    dataset = {"seed": SEEDS[0], "methods": {METHODS[0]: {"output": str(tmp_path / "output/method")}}}
    with pytest.raises(ValueError, match="different inference identity"):
        fragment_rows(dataset, {"execution": {"absolute_path": str(tmp_path / "other.json")}}, {}, {})


def records():
    score = {"tp": 1, "fp": 1, "fn": 1, "f1": .5, "precision": .5, "recall": .5,
             "predicted_pairs": 2, "eligible_true_pairs": 2, "pair_endpoint_coverage": 1.0}
    return [{"arm": arm, "seed": seed, "method": method, "status": "complete", "score": deepcopy(score),
             "strata": [{"fragment_endpoints": i, "score": deepcopy(score)} for i in range(3)]}
            for arm in ("baseline", "fragment") for seed in SEEDS for method in METHODS]


def test_means_only_defined_successes_with_explicit_seed_lists():
    rows = records()
    target = next(r for r in rows if r["arm"] == "fragment" and r["method"] == METHODS[0] and r["seed"] == SEEDS[0])
    target.update(status="failed")
    target.pop("score")
    second = next(r for r in rows if r["arm"] == "fragment" and r["method"] == METHODS[0] and r["seed"] == SEEDS[1])
    second["score"]["precision"] = None
    actual = next(r for r in metric_means(rows) if r["arm"] == "fragment" and r["method"] == METHODS[0])
    assert actual["failed_or_unavailable_seeds"] == [SEEDS[0]]
    assert actual["metrics"]["f1"]["eligible_seeds"] == list(SEEDS[1:])
    assert actual["metrics"]["precision"]["eligible_seeds"] == list(SEEDS[2:])
    assert actual["metrics"]["precision"]["mean"] == .5
    with pytest.raises(ValueError, match="Incomplete seed"):
        metric_means(rows[:-1])


def test_complete_tsv_inventory_nulls_and_failed_rows_not_zero(tmp_path):
    rows = records()
    rows[-1] = {k: v for k, v in rows[-1].items() if k not in {"score", "strata"}}
    rows[-1].update(status="unavailable", reason="native failure")
    rows[0]["score"]["precision"] = None
    comparisons = paired_intervals(rows)
    tables(tmp_path, rows, comparisons)
    with (tmp_path / "scores.tsv").open() as handle:
        scores = list(csv.DictReader(handle, delimiter="\t"))
    assert len(scores) == 80 and scores[0]["precision"] == "NA"
    assert scores[-1]["tp"] == "NA" and scores[-1]["f1"] == "NA"
    with (tmp_path / "strata.tsv").open() as handle:
        strata = list(csv.DictReader(handle, delimiter="\t"))
    assert len(strata) == 240 and all(r["fp"] == "NA" for r in strata[-3:])
    with (tmp_path / "comparisons.tsv").open() as handle:
        effects = list(csv.DictReader(handle, delimiter="\t"))
    assert len(effects) == 15
    for row, source in zip(effects, comparisons):
        assert int(row["paired_seeds"]) == len(source["eligible_seeds"])
        assert float(row["estimate"]) == source["estimate"]
        assert float(row["adjusted_low"]) == source["adjusted_interval"][0]
    with pytest.raises(FileExistsError):
        tables(tmp_path, rows, comparisons)


def test_existing_scientific_report_never_replaced(tmp_path):
    output = tmp_path / "old_report"
    output.mkdir()
    with pytest.raises(FileExistsError):
        assemble(tmp_path, tmp_path, "sha", tmp_path, "sha", tmp_path, 1, output)


def test_terminal_all_failed_panel_preserves_complete_inventory_and_missing_estimates(tmp_path, monkeypatch):
    result_directory = tmp_path / "benchmark_tools/results"
    result_directory.mkdir(parents=True)
    protocol = result_directory / "CONTROLLED_FRAGMENT_OBSERVATION_PROTOCOL_20261009.md"
    protocol.write_text("test protocol")
    parent = result_directory / "test_parent.json"
    parent.write_text("{}")
    monkeypatch.setattr("benchmark_tools.assemble_controlled_fragment_results.PINS", {parent.name: record(parent)["sha256"]})
    baseline = [row for row in records() if row["arm"] == "baseline"]
    datasets = [{"label": f"fragment20_center60_v1_{seed}", "seed": seed,
                 "baseline_records": [row for row in baseline if row["seed"] == seed]} for seed in SEEDS]
    panel = {"schema": "controlled_fragment_observations_v1", "status": "prepared_unexecuted",
             "datasets": datasets, "inference_identities": 30,
             "protocol": record(protocol), "baseline_pins": [record(parent)]}
    manifest, runtime = tmp_path / "manifest.json", tmp_path / "runtime.json"
    manifest.write_text(json.dumps(panel))
    runtime.write_text("{}")
    execution_root = tmp_path / "execution"
    execution_root.mkdir()
    status = {"schema": "controlled_fragment_panel_execution_v1", "status": "finished_pending_native_validation",
              "accuracy_evaluated": False, "provenance": {"fragment_manifest": record(manifest),
                "runtime_manifest": record(runtime), "allocation": {"job_id": "1"}, "sources": []},
              "datasets": [{"label": d["label"], "status": "failed", "reason": "fixture preflight failure"} for d in datasets]}
    (execution_root / "panel_status.json").write_text(json.dumps(status))
    monkeypatch.setattr("benchmark_tools.assemble_controlled_fragment_results.verify_environment", lambda *a: None)
    monkeypatch.setattr("benchmark_tools.assemble_controlled_fragment_results.observed_inputs", lambda *a: {})
    monkeypatch.setattr("benchmark_tools.assemble_controlled_fragment_results.subprocess.check_output",
                        lambda *a, **k: "JobIDRaw|State|ExitCode|Elapsed\n1|COMPLETED|0:0|00:01\n")
    output = tmp_path / "results"
    result = assemble(tmp_path, manifest, hashlib.sha256(manifest.read_bytes()).hexdigest(), runtime,
                      hashlib.sha256(runtime.read_bytes()).hexdigest(), execution_root, 1, output)
    assert result["status"] == "complete_with_explicit_outcomes"
    assert len(result["records"]) == 80 and result["fragment_successes"] == 0
    assert all(row["estimate"] is None and len(row["excluded_seeds"]) == 10 for row in result["comparisons"])
    assert result["publication_ready"] is False
    assert all("score" not in row for row in result["records"] if row["arm"] == "fragment")
    assert json.loads((output / "report.json").read_text()) == result
