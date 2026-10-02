import copy
import json
from pathlib import Path
import subprocess
import sys

import numpy as np
import pytest

from benchmark_tools import reproduce_simulation_panels as replay
from benchmark_tools.summarize_simulation_panel import summarize


RESULTS = Path(__file__).resolve().parents[2] / "benchmark_tools" / "results"
PATHS = [RESULTS / f"simulation_{length}_native_results_20260916.json" for length in ("fixed", "variable")]


@pytest.fixture(scope="module")
def reports():
    return [json.loads(path.read_bytes()) for path in PATHS]


@pytest.mark.parametrize("index,complete,failed,estimated,bounds", [(0, 134, 146, 0, 0), (1, 267, 13, 14, 112)])
def test_original_reports_reproduce(reports, index, complete, failed, estimated, bounds):
    result = replay.verify(reports[index])
    assert result["scored_records"] == complete
    assert result["failed_records"] == failed
    assert result["estimated_contrasts"] == estimated
    assert result["unavailable_contrasts"] == 14 - estimated
    assert result["interval_bound_values"] == bounds
    assert result["method_mean_cells"] == 84


@pytest.mark.parametrize("fault", ["replicates", "bootstrap_seed", "adjustment", "diagnostic_contrast", "duplicate", "missing",
    "wrong_truth", "bad_count", "boolean_count", "nan_metric", "coverage", "pending", "failure_score", "missing_reason",
    "checkpoint_without_parent", "imputed_mean", "imputed_contrast", "outcome_seeds", "conditioning", "bound", "role", "bound_shape"])
def test_corruption_rejected(reports, fault):
    report = copy.deepcopy(reports[1] if fault in ("bound", "role", "bound_shape") else reports[0])
    complete = next(r for r in report["records"] if r["status"] == "complete")
    failed = next(r for r in report["records"] if r["status"] == "failed")
    block = report["conditions"]["baseline"]
    contrast = block["contrasts"][replay.METHODS[0]]
    if fault == "replicates":
        report["bootstrap"]["replicates"] = 1000
    elif fault == "bootstrap_seed":
        report["bootstrap"]["seed"] += 1
    elif fault == "adjustment":
        report["bootstrap"]["f1_multiplicity_count"] = 21
    elif fault == "diagnostic_contrast":
        block["contrasts"][replay.METHODS[3]] = copy.deepcopy(contrast)
    elif fault == "duplicate":
        report["records"].append(copy.deepcopy(complete))
    elif fault == "missing":
        report["records"].pop()
    elif fault == "wrong_truth":
        complete["truth_sha256"] = "b" * 64
    elif fault == "bad_count":
        complete["score"]["tp"] += 1
    elif fault == "boolean_count":
        complete["score"]["tp"] = True
    elif fault == "nan_metric":
        complete["score"]["f1"] = float("nan")
    elif fault == "coverage":
        complete["score"]["genes_in_predicted_pairs"] += 1
    elif fault == "pending":
        failed["status"] = "pending"
    elif fault == "failure_score":
        failed["score"] = complete["score"]
    elif fault == "missing_reason":
        failed.pop("reason")
    elif fault == "checkpoint_without_parent":
        checkpoint = next(r for r in report["records"] if r["method"] == replay.METHODS[3] and r["condition"] == complete["condition"] and r["seed"] == complete["seed"])
        checkpoint.update(status="complete", truth_sha256=complete["truth_sha256"], score=complete["score"])
    elif fault == "imputed_mean":
        block["methods"][replay.METHODS[2]]["available_case_means"]["f1"] = 0
    elif fault == "imputed_contrast":
        contrast["metrics"] = {}
    elif fault == "outcome_seeds":
        block["methods"][replay.METHODS[0]]["complete_seeds"].pop()
    elif fault == "conditioning":
        contrast["conditional_on_success"] = False
    elif fault == "bound":
        contrast["metrics"]["f1"]["bonferroni_14_ci"][0] += .01
    elif fault == "role":
        contrast["metrics"]["precision"]["inference_role"] = "primary"
    else:
        contrast["metrics"]["f1"]["paired_95_percent_ci"] = [contrast["metrics"]["f1"]["paired_95_percent_ci"]]
    with pytest.raises(ValueError):
        replay.verify(report)


def synthetic_report(pair_count):
    rows = []
    for condition in replay.CONDITIONS:
        for index, seed in enumerate(replay.PANELS["fixed_length_v1"]["seeds"]):
            for method in replay.METHODS:
                row = dict(condition=condition, seed=seed, method=method)
                if index >= pair_count and method in replay.METHODS[2:]:
                    row.update(status="failed", reason="synthetic terminal failure")
                else:
                    tp, fp = (1, 1) if index == 0 and method == replay.METHODS[0] else (100, 0)
                    row.update(status="complete", truth_sha256="a" * 64, score=dict(
                        tp=tp, fp=fp, fn=100 - tp, input_genes=20, eligible_true_pairs=100, predicted_pairs=tp + fp,
                        genes_in_predicted_pairs=20, duplicate_prediction_rows=0, pair_endpoint_coverage=1.,
                        coverage_definition="input genes occurring in at least one predicted cross-species pair",
                        f1=2 * tp / (tp + fp + 100), precision=tp / (tp + fp), recall=tp / 100, undefined_ratios=[]))
                rows.append(row)
    return summarize(rows)


@pytest.mark.parametrize("pairs", [0, 1, 2])
def test_unavailable_single_pair_and_seed_means_not_pooled(pairs):
    report = synthetic_report(pairs)
    result = replay.verify(report)
    assert result["unavailable_contrasts"] == (14 if pairs == 0 else 0)
    assert result["insufficient_seed_contrasts"] == (14 if pairs == 1 else 0)
    if pairs == 2:
        point = report["conditions"]["baseline"]["contrasts"][replay.METHODS[0]]["metrics"]["f1"]["comparator_mean"]
        assert point == pytest.approx((2 / 102 + 1) / 2)
        assert point != pytest.approx(202 / 302)


def test_explicit_draws_match_original_seed_resampling():
    differences = np.array([[-20., -3., 2.], [1., 4., -5.]])
    nominal, adjusted = replay.paired_intervals(differences, 20261130)
    weights = np.random.Generator(np.random.PCG64(20261130)).multinomial(2, [.5, .5], size=20000)
    draws = (weights[:, :1] * differences[0] + weights[:, 1:] * differences[1]) / 2
    assert nominal == pytest.approx(np.quantile(draws, [.025, .975], axis=0))
    assert adjusted == pytest.approx(np.quantile(draws[:, 0], [.025 / 14, 1 - .025 / 14]))


@pytest.mark.parametrize("actual", [True, float("nan"), [1], [[1, 1]], "1"])
def test_no_broadcast_or_nonfinite_arithmetic(actual):
    with pytest.raises(ValueError):
        replay.close(actual, 1.)


def test_undefined_zero_ratios_are_explicit():
    score = dict(tp=0, fp=0, fn=0, input_genes=0, eligible_true_pairs=0, predicted_pairs=0,
        genes_in_predicted_pairs=0, duplicate_prediction_rows=0, pair_endpoint_coverage=0.,
        coverage_definition="input genes occurring in at least one predicted cross-species pair",
        f1=0., precision=0., recall=0., undefined_ratios=list(replay.METRICS))
    assert replay.score_values(score) == [0., 0., 0.]
    score["undefined_ratios"] = []
    with pytest.raises(ValueError):
        replay.score_values(score)


@pytest.mark.parametrize("fault", ["one_input", "duplicate_panels", "altered_bytes", "symlink", "numpy_version"])
def test_execution_guards_retain_failure(tmp_path, monkeypatch, fault):
    paths = list(PATHS)
    if fault == "one_input":
        paths.pop()
    elif fault == "duplicate_panels":
        paths[1] = paths[0]
    elif fault == "altered_bytes":
        changed = tmp_path / "changed.json"
        changed.write_bytes(paths[0].read_bytes() + b"\n")
        paths[0] = changed
    elif fault == "symlink":
        link = tmp_path / "linked.json"
        link.symlink_to(paths[0])
        paths[0] = link
    else:
        monkeypatch.setattr(np, "__version__", "wrong")
    output = tmp_path / "failure.json"
    with pytest.raises(ValueError):
        replay.reproduce(paths, output)
    receipt = json.loads(output.read_bytes())
    assert receipt["status"] == "validation_failed"
    assert receipt["native_inference_repeated"] is False
    assert receipt["publication_ready"] is False


def test_no_overwrite(tmp_path):
    output = tmp_path / "existing.json"
    output.write_bytes(b"retained")
    with pytest.raises(FileExistsError):
        replay.reproduce(PATHS, output)
    assert output.read_bytes() == b"retained"


def test_standalone_cli_from_copied_directory(tmp_path):
    script = tmp_path / "reproduce_simulation_panels.py"
    script.write_bytes(Path(replay.__file__).read_bytes())
    paths = []
    for source in PATHS:
        target = tmp_path / source.name
        target.write_bytes(source.read_bytes())
        paths += ["--results", str(target)]
    output = tmp_path / "reproduced.json"
    child = subprocess.run([sys.executable, "-I", "-B", str(script), *paths, "--output", str(output)],
                           cwd=tmp_path, capture_output=True, text=True, timeout=30)
    assert child.returncode == 0, child.stderr
    result = json.loads(output.read_bytes())
    assert result["status"] == "frozen_simulation_panels_arithmetically_reproduced"
    assert result["panels_pooled"] is False
    assert result["raw_scoring_repeated"] is False
    assert result["native_inference_repeated"] is False
    assert result["publication_ready"] is False
    assert result == json.loads(child.stdout)


def test_guide_keeps_panels_and_inferential_scope_separate():
    guide = (RESULTS.parent / "SIMULATION_ARITHMETIC_REPLAY.md").read_text()
    assert "| Fixed length | 134 | 146 | 0 | 14 |" in guide
    assert "| Variable length | 267 | 13 | 14 | 0 |" in guide
    assert "42 paired metric" in guide and "112 scalar interval bounds" in guide
    assert "does not validate tree-perturbation or matched-recall panels" in guide
