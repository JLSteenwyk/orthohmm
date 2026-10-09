from copy import deepcopy
import hashlib
import json
from pathlib import Path

import numpy as np
import pytest

from benchmark_tools.controlled_fragment_observations import (
    CONDITION, SEEDS, METHODS, INFERENCE_METHODS, METRICS, COMPARISONS,
    observe, verify_observation, score_strata, paired_intervals, fresh_methods,
)
from benchmark_tools.prepare_controlled_fragment_observations import (
    baseline_rows, bind_baseline, read_input_sequences, record, PINS,
)
from benchmark_tools.simulation_conditions import score_pairs


def observations(n=15):
    sequences = {f"g{i}": "ACDEFGHIKLMNPQRSTVWY"[:10 + i % 10] for i in range(n)}
    owners = {gene: "A" if i % 2 else "B" for i, gene in enumerate(sequences)}
    return sequences, owners


def test_exact_label_independent_selection_and_center_substrings():
    sequences, owners = observations()
    original = deepcopy(sequences)
    actual, coordinates = observe(sequences, owners, 20261101)
    independent = sorted(sequences, key=lambda gene: (
        hashlib.sha256(f"20261101:{CONDITION}:{gene}".encode()).hexdigest(), gene))[:3]
    assert [row["gene"] for row in coordinates if row["fragment"]] == sorted(independent)
    for row in coordinates:
        length = len(sequences[row["gene"]])
        if row["fragment"]:
            assert row["observed_length"] == 3 * length // 5
            assert row["start"] == (length - row["observed_length"]) // 2
        else:
            assert actual[row["gene"]] == sequences[row["gene"]]
    assert sequences == original and len(actual) == 15
    assert verify_observation(sequences, actual, owners, coordinates, 20261101)["fragment_genes"] == 3
    reversed_sequences = dict(reversed(list(sequences.items())))
    assert observe(reversed_sequences, owners, 20261101) == (actual, coordinates)
    renamed_owners = {gene: "other" for gene in owners}
    assert observe(sequences, renamed_owners, 20261101)[0] == actual


@pytest.mark.parametrize("n", [1, 4, 5, 9, 10, 14, 100])
def test_exact_floor_fraction(n):
    sequences, owners = observations(n)
    actual, coordinates = observe(sequences, owners, 2)
    assert sum(row["fragment"] for row in coordinates) == n // 5
    assert verify_observation(sequences, actual, owners, coordinates, 2)["genes"] == n


@pytest.mark.parametrize("seed", [True, 0, -1, 1.1, "1"])
def test_invalid_seeds(seed):
    sequences, owners = observations()
    with pytest.raises(ValueError, match="integer seed"):
        observe(sequences, owners, seed)


@pytest.mark.parametrize("fault", ["empty", "missing_owner", "empty_owner", "empty_sequence", "whitespace", "short_selected"])
def test_invalid_parent_not_rescued(fault):
    sequences, owners = observations()
    if fault == "empty":
        sequences, owners = {}, {}
    elif fault == "missing_owner":
        owners.pop("g1")
    elif fault == "empty_owner":
        owners["g1"] = ""
    elif fault == "empty_sequence":
        sequences["g1"] = ""
    elif fault == "whitespace":
        sequences["g1"] = "AA BB"
    else:
        selected = next(row["gene"] for row in observe(sequences, owners, 2)[1] if row["fragment"])
        sequences[selected] = "A"
    with pytest.raises(ValueError):
        observe(sequences, owners, 2)


@pytest.mark.parametrize("fault", ["changed_sequence", "missing_gene", "added_gene", "wrong_owner", "wrong_flag", "wrong_hash", "wrong_coordinate", "duplicate_row", "wrong_order"])
def test_readback_rejects_observation_corruption(fault):
    sequences, owners = observations()
    actual, coordinates = observe(sequences, owners, 2)
    gene = coordinates[0]["gene"]
    if fault == "changed_sequence":
        actual[gene] = "XXXXX"
    elif fault == "missing_gene":
        actual.pop(gene)
    elif fault == "added_gene":
        actual["new"] = "AAA"
    elif fault == "wrong_owner":
        coordinates[0]["species"] = "other"
    elif fault == "wrong_flag":
        coordinates[0]["fragment"] = 1
    elif fault == "wrong_hash":
        coordinates[0]["parent_sha256"] = "0" * 64
    elif fault == "wrong_coordinate":
        coordinates[0]["start"] += 1
    elif fault == "duplicate_row":
        coordinates[-1] = deepcopy(coordinates[0])
    else:
        coordinates.reverse()
    with pytest.raises(ValueError):
        verify_observation(sequences, actual, owners, coordinates, 2)


def test_all_pairs_partition_even_inter_family_false_positives_and_duplicates():
    owners = {"a": "A", "b": "B", "c": "A", "d": "B", "e": "A", "f": "B"}
    flags = {"a": False, "b": False, "c": True, "d": False, "e": True, "f": True}
    truth = [("a", "b"), ("c", "d"), ("e", "f")]
    predictions = [("a", "b"), ("b", "a"), ("c", "b"), ("e", "f"), ("c", "f")]
    result = score_strata(iter(predictions), iter(truth), owners, flags)
    assert result["score"] == score_pairs(predictions, truth, owners)
    assert result["score"]["duplicate_prediction_rows"] == 1
    assert [(r["score"]["tp"], r["score"]["fp"], r["score"]["fn"]) for r in result["strata"]] == [(1, 0, 0), (0, 1, 1), (1, 1, 0)]
    assert result["score"]["tp"] == 2 and result["score"]["fp"] == 2 and result["score"]["fn"] == 1


def test_zero_denominators_are_null_not_perfect_or_imputed():
    owners = {"a": "A", "b": "B"}
    result = score_strata([], [], owners, {"a": False, "b": True})
    for score in [result["score"], *(r["score"] for r in result["strata"])]:
        assert all(score[metric] is None for metric in METRICS)
        assert score["undefined_ratios"] == list(METRICS)
    result = score_strata([], [("a", "b")], owners, {"a": False, "b": True})
    assert result["score"]["precision"] is None
    assert result["score"]["recall"] == 0 and result["score"]["f1"] == 0


@pytest.mark.parametrize("fault", ["unknown", "same_species", "duplicate_truth", "missing_flag", "nonboolean_flag"])
def test_invalid_pair_or_flag_universe(fault):
    owners = {"a": "A", "b": "B", "c": "A"}
    truth, predictions, flags = [("a", "b")], [("a", "b")], {"a": False, "b": True, "c": False}
    if fault == "unknown":
        predictions = [("a", "z")]
    elif fault == "same_species":
        predictions = [("a", "c")]
    elif fault == "duplicate_truth":
        truth.append(("b", "a"))
    elif fault == "missing_flag":
        flags.pop("c")
    else:
        flags["a"] = 1
    with pytest.raises(ValueError):
        score_strata(predictions, truth, owners, flags)


def comparison_records():
    return [{"arm": arm, "seed": seed, "method": method, "status": "complete",
             "score": {metric: .5 + (seed - SEEDS[0]) / 100 + METHODS.index(method) / 20
                       + (.02 * (seed - SEEDS[0]) if arm == "fragment" else 0) for metric in METRICS}}
            for arm in ("baseline", "fragment") for seed in SEEDS for method in METHODS]


def test_bootstrap_matches_independent_seed_resampling_fixed_family():
    records = comparison_records()
    actual = paired_intervals(records)
    assert len(actual) == 15
    for row in actual:
        values = [next(r for r in records if (r["arm"], r["seed"], r["method"]) == (row["target_arm"], s, row["target_method"]))["score"][row["metric"]]
                  - next(r for r in records if (r["arm"], r["seed"], r["method"]) == (row["reference_arm"], s, row["reference_method"]))["score"][row["metric"]] for s in SEEDS]
        rng = np.random.Generator(np.random.PCG64(20261011))
        draws = rng.integers(0, 10, (20000, 10))
        samples = np.asarray(values)[draws].mean(axis=1)
        assert row["estimate"] == np.mean(values)
        assert row["nominal_interval"] == np.quantile(samples, [.025, .975], method="linear").tolist()
        assert row["adjusted_interval"] == np.quantile(samples, [.05 / 30, 1 - .05 / 30], method="linear").tolist()
        assert row["eligible_seeds"] == list(SEEDS) and row["planned_endpoints"] == 15
    for start in range(0, 15, 3):
        assert actual[start]["nominal_interval"] == actual[start + 1]["nominal_interval"] == actual[start + 2]["nominal_interval"]


def test_failed_and_undefined_metrics_use_paired_lists_and_explicit_exclusions():
    records = comparison_records()
    failed = next(r for r in records if r["arm"] == "fragment" and r["seed"] == SEEDS[0] and r["method"] == METHODS[0])
    failed.update(status="failed", reason="native failure")
    failed.pop("score")
    null = next(r for r in records if r["arm"] == "baseline" and r["seed"] == SEEDS[1] and r["method"] == METHODS[0])
    null["score"]["precision"] = None
    result = paired_intervals(records)
    assert result[0]["eligible_seeds"] == list(SEEDS[1:])
    assert result[1]["eligible_seeds"] == list(SEEDS[2:])
    assert result[1]["excluded_seeds"][0]["reasons"][0]["reason"] == "native failure"
    assert result[1]["excluded_seeds"][1]["reasons"][0]["status"] == "undefined_metric"
    assert all(row["planned_endpoints"] == 15 for row in result)


def test_no_successful_pairs_is_unavailable_not_zero():
    records = [{k: v for k, v in row.items() if k != "score"} for row in comparison_records()]
    for row in records:
        row["status"] = "failed"
    for row in paired_intervals(records):
        assert row["status"] == "unavailable" and row["estimate"] is None and row["nominal_interval"] is None
        assert len(row["excluded_seeds"]) == 10


@pytest.mark.parametrize("fault", ["missing", "duplicate", "unknown_method", "unknown_seed", "unknown_arm", "running", "failed_score", "nan", "infinite", "boolean", "negative", "over_one"])
def test_invalid_comparison_inventory_rejected(fault):
    rows = comparison_records()
    if fault == "missing":
        rows.pop()
    elif fault == "duplicate":
        rows.append(deepcopy(rows[0]))
    elif fault.startswith("unknown_"):
        rows[0][fault.removeprefix("unknown_")] = "other"
    elif fault == "running":
        rows[0]["status"] = "running"
    elif fault == "failed_score":
        rows[0]["status"] = "failed"
    else:
        rows[0]["score"]["f1"] = {"nan": float("nan"), "infinite": float("inf"), "boolean": True,
                                   "negative": -.1, "over_one": 1.1}[fault]
    with pytest.raises(ValueError):
        paired_intervals(rows)


@pytest.fixture
def commands(tmp_path):
    old = tmp_path / "old"
    methods = {}
    for name in INFERENCE_METHODS[:2]:
        methods[name] = {"output": str(old / name), "metrics": str(old / (name + ".json")),
                         "argv": ["python", "harness", str(old / "inputs"), str(old / name), str(old / (name + ".json")),
                                  "--cpu", "4", "--threads-per-worker", "4", "--accuracy-profile", "high_sensitivity"]}
    methods[METHODS[1]]["argv"] += ["--phylogeny", "reconcile", "--species-tree-mode", "infer"]
    methods[METHODS[2]] = {"output": str(old / METHODS[2]), "copy_inputs_from": str(old / "inputs"),
                          "copy_inputs_to": str(old / METHODS[2] / "input"),
                          "argv": ["of", "-f", str(old / METHODS[2] / "input"), "-t", "4", "-a", "4", "-S", "diamond"]}
    methods[METHODS[3]] = {"output": methods[METHODS[2]]["output"], "parent_method": METHODS[2]}
    return {"input": str(old / "inputs"), "methods": methods}


def test_new_commands_change_only_paths_and_preserve_sequence_checkpoint(commands, tmp_path):
    original = deepcopy(commands)
    new = fresh_methods(commands, tmp_path / "new_inputs", tmp_path / "new_outputs")
    assert commands == original
    for method in INFERENCE_METHODS[:2]:
        assert new[method]["argv"][:2] == original["methods"][method]["argv"][:2]
        assert new[method]["argv"][5:] == original["methods"][method]["argv"][5:]
        assert new[method]["argv"][2] == str(tmp_path / "new_inputs")
    before, after = original["methods"][METHODS[2]]["argv"], new[METHODS[2]]["argv"]
    assert before[:2] == after[:2] and before[3:] == after[3:]
    assert new[METHODS[3]]["output"] == new[METHODS[2]]["output"]
    assert "argv" not in new[METHODS[3]]


@pytest.mark.parametrize("fault", ["existing", "same_input", "nested_input", "overlap_new", "old_output", "old_metrics", "bad_hmm", "bad_copy", "duplicate_f"])
def test_fresh_commands_refuse_overlap_restart_or_malformed_baselines(commands, tmp_path, fault):
    inputs, output = tmp_path / "new_input", tmp_path / "new_output"
    if fault == "existing":
        output.mkdir()
    elif fault == "same_input":
        inputs = Path(commands["input"])
    elif fault == "nested_input":
        output = Path(commands["input"]) / "new"
    elif fault == "overlap_new":
        inputs = output / "input"
    elif fault == "old_output":
        output = Path(commands["methods"][METHODS[0]]["output"])
    elif fault == "old_metrics":
        inputs = Path(commands["methods"][METHODS[0]]["metrics"])
    elif fault == "bad_hmm":
        commands["methods"][METHODS[0]]["argv"][2] = "other"
    elif fault == "bad_copy":
        commands["methods"][METHODS[2]]["copy_inputs_to"] = "other"
    else:
        commands["methods"][METHODS[2]]["argv"] += ["-f", "other"]
    with pytest.raises((ValueError, FileExistsError)):
        fresh_methods(commands, inputs, output)


def test_complete_baseline_inventory_and_failure_not_selected_away():
    rows = [{"condition": "baseline", "seed": seed, "method": method, "status": "complete"} for seed in SEEDS for method in METHODS]
    assert len(baseline_rows({"records": rows})) == 40
    for changed in (rows[:-1], rows + [rows[0]], [dict(r, status="failed") if i == 0 else r for i, r in enumerate(rows)]):
        with pytest.raises(ValueError, match="40 prespecified"):
            baseline_rows({"records": changed})


@pytest.fixture
def input_fixture(tmp_path):
    inputs = tmp_path / "input"
    inputs.mkdir()
    (inputs / "A.fasta").write_text(">a\nACDEFG\n>c\nHIKLMN\n")
    (inputs / "B.fasta").write_text(">b\nACDEFG\n>d\nHIKLMN\n")
    truth = {"extant_genes": 4, "species": ["A", "B"], "ortholog_pairs": [["a", "b"], ["c", "d"]],
             "ortholog_pair_count": 2, "families": {"1": ["a", "b"], "2": ["c", "d"]},
             "prepared_inputs": [{"path": "input/" + p.name, "bytes": p.stat().st_size,
                                  "sha256": hashlib.sha256(p.read_bytes()).hexdigest()} for p in sorted(inputs.iterdir())]}
    truth_path = tmp_path / "truth.json"
    truth_path.write_text(json.dumps(truth))
    verified = {"status": "ready", "truth": record(truth_path), "inputs": [record(p) for p in sorted(inputs.iterdir())]}
    return inputs, truth, truth_path, verified


def test_inputs_bound_to_truth_families_and_exact_file_inventory(input_fixture):
    inputs, truth, _, verified = input_fixture
    sequences, owners, species = read_input_sequences(verified, truth, inputs)
    assert sequences["a"] == "ACDEFG" and owners == {"a": "A", "c": "A", "b": "B", "d": "B"}
    assert species == ["A", "B"]


@pytest.mark.parametrize("fault", ["extra_file", "changed_bytes", "wrong_pin", "missing_family_gene", "duplicate_family_gene", "unknown_pair", "duplicate_pair"])
def test_input_and_truth_corruption_rejected(input_fixture, fault):
    inputs, truth, _, verified = input_fixture
    if fault == "extra_file":
        (inputs / "note.txt").write_text("unrecorded")
    elif fault == "changed_bytes":
        (inputs / "A.fasta").write_text(">other\nAA\n")
    elif fault == "wrong_pin":
        truth["prepared_inputs"][0]["sha256"] = "wrong"
    elif fault == "missing_family_gene":
        truth["families"]["1"].pop()
    elif fault == "duplicate_family_gene":
        truth["families"]["2"].append("a")
    elif fault == "unknown_pair":
        truth["ortholog_pairs"][0] = ["a", "unknown"]
    else:
        truth["ortholog_pairs"].append(["b", "a"])
        truth["ortholog_pair_count"] += 1
    with pytest.raises(ValueError):
        read_input_sequences(verified, truth, inputs)


@pytest.fixture
def binding_fixture(input_fixture, tmp_path, monkeypatch):
    inputs, truth, truth_path, verified = input_fixture
    seed = SEEDS[0]
    dataset = {"input": str(inputs), "truth": str(truth_path), "label": f"baseline_{seed}", "seed": seed, "methods": {}}
    rows, tasks, executors = {}, {}, {}
    native_hash = PINS["publication_variable_native_methods_20260916.json"]
    comparator_hash = PINS["publication_variable_methods_20260916.json"]
    generation_hash = PINS["publication_variable_simulation_manifest_20260916.json"]
    names = ("run_simulation_methods.py", "verify_simulation_histories.py", "run_simulation_generation.py", "benchmark_production.py")
    predictions = [("a", "b")]
    for kind, method_names, manifest_hash in (("native", METHODS[:2], native_hash), ("comparator", METHODS[2:], comparator_hash)):
        directory = tmp_path / kind
        directory.mkdir()
        executor = tmp_path / (kind + "_executor")
        (executor / "benchmark_tools").mkdir(parents=True)
        for name in names:
            (executor / "benchmark_tools" / name).write_text("pass\n")
        sources = [{k: v for k, v in record(executor / "benchmark_tools" / n).items() if k != "absolute_path"} for n in names]
        task = {"JobID": "10_0", "JobIDRaw": "11", "State": "COMPLETED", "ExitCode": "0:0"}
        tasks[kind], executors[kind] = task, executor
        status = {"dataset": dataset["label"], "verified_inputs": verified, "status": "finished_pending_native_validation",
                  "provenance": {"method_manifest_sha256": manifest_hash, "generation_manifest_sha256": generation_hash,
                                 "slurm_job_id": "11", "slurm_array_task_id": "0", "sources": sources}, "methods": {}}
        for method in method_names:
            parent = METHODS[2] if method == METHODS[3] else method
            output = directory / parent
            output.mkdir(exist_ok=True)
            artifact = output / (method + ".txt")
            artifact.write_text("a b\n")
            dataset["methods"][method] = {"output": str(output), "argv": [parent]}
            status["methods"].setdefault(parent, {"status": "process_succeeded", "argv": [parent], "exit_code": 0, "outputs": []})["outputs"].append(record(artifact))
            rows[seed, method] = {"inference_method_manifest_sha256": manifest_hash, "reused_comparator": kind == "comparator",
                                  "scheduler": task, "truth_sha256": verified["truth"]["sha256"],
                                  "native_validation": {"completion": "test admission"},
                                  "prediction_artifacts": [record(artifact)], "score": score_pairs(predictions, truth["ortholog_pairs"], {"a": "A", "b": "B", "c": "A", "d": "B"})}
        (directory / "execution").mkdir()
        status_path = directory / "execution/status.json"
        status_path.write_text(json.dumps(status))
        for method in method_names:
            rows[seed, method]["execution_evidence"] = record(status_path)
    monkeypatch.setattr("benchmark_tools.prepare_controlled_fragment_observations.admit_method", lambda *args: {"status": "admitted", "native_validation": {"completion": "test admission"}})
    monkeypatch.setattr("benchmark_tools.prepare_controlled_fragment_observations.load_predictions",
                        lambda method, output, *args: (predictions, [output / (method + ".txt")]))
    manifest = {"generation_manifest": {"sha256": generation_hash}}
    return dataset, rows, manifest, manifest, tasks, executors


def test_baseline_binding_rechecks_retained_scores_and_native_artifacts(binding_fixture):
    bound = bind_baseline(*binding_fixture)
    assert len(bound["bindings"]) == 4 and bound["predictions"][METHODS[0]] == [("a", "b")]
    assert bound["owners"] == {"a": "A", "c": "A", "b": "B", "d": "B"}


@pytest.mark.parametrize("fault", ["manifest", "reuse", "scheduler", "execution_path", "truth_sha", "native_admission", "predictions", "score", "wrong_argv", "changed_truth"])
def test_binding_rejects_unmatched_provenance_or_results(binding_fixture, tmp_path, fault):
    dataset, rows, manifest, previous, tasks, executors = binding_fixture
    row = rows[SEEDS[0], METHODS[0]]
    if fault == "manifest":
        row["inference_method_manifest_sha256"] = "other"
    elif fault == "reuse":
        row["reused_comparator"] = True
    elif fault == "scheduler":
        row["scheduler"] = dict(row["scheduler"], JobIDRaw="other")
    elif fault == "execution_path":
        row["execution_evidence"]["absolute_path"] = str(tmp_path / "other.json")
    elif fault == "truth_sha":
        row["truth_sha256"] = "other"
    elif fault == "native_admission":
        row["native_validation"]["completion"] = "other"
    elif fault == "predictions":
        row["prediction_artifacts"][0]["sha256"] = "other"
    elif fault == "score":
        row["score"]["tp"] = 99
    elif fault == "wrong_argv":
        dataset["label"] = "other"
    else:
        Path(dataset["truth"]).write_text("{}")
    with pytest.raises(ValueError):
        bind_baseline(dataset, rows, manifest, previous, tasks, executors)
