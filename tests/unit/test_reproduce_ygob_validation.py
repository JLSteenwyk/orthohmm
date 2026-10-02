import copy
import gzip
import hashlib
import json
from pathlib import Path

import numpy as np
import pytest

from benchmark_tools import reproduce_ygob_validation as replay


ROOT = Path(__file__).resolve().parents[2]
SUMMARY = ROOT / "benchmark_tools/results/ygob_frozen_results_20260916.json"
SNAPSHOT = ROOT / "benchmark_tools/results/ygob_sufficient_counts_20261002.json.gz"
SNAPSHOT_SHA = "62885bb85eb357d97e55d053134d262e1ebb1e2b088a37df96986f98f8620d66"


@pytest.fixture
def inputs():
    return json.loads(SUMMARY.read_bytes()), json.loads(gzip.decompress(SNAPSHOT.read_bytes()))


def test_retained_snapshot_identity_and_all_point_metrics(inputs):
    summary, snapshot = inputs
    assert hashlib.sha256(SNAPSHOT.read_bytes()).hexdigest() == SNAPSHOT_SHA
    assert hashlib.sha256(SUMMARY.read_bytes()).hexdigest() == replay.SUMMARY_SHA
    assert hashlib.sha256(gzip.decompress(SNAPSHOT.read_bytes())).hexdigest() == "2b931b0d11dba223f5f4ee3617eb88ae9a019f1fecb71809472d2ea812d6b89a"
    arrays = replay.verify_points(summary, snapshot)
    assert set(arrays) == set(replay.METHODS)
    assert all(array.shape == (10250, 3) and array.dtype == np.int64 for array in arrays.values())
    assert sum(len(row) for method in snapshot["methods"].values() for row in method) == 205000
    assert all(type(v) is int for row in snapshot["methods"][replay.METHODS[0]] for v in row)


@pytest.mark.parametrize("problem", ["native_scope", "extra_identifiers", "method", "columns", "missing_row",
                                   "genes", "float_count", "bool_count", "negative", "truth_count", "half_count",
                                   "covered", "exact", "point", "independent_counts", "coverage", "projected",
                                   "seed", "multiplicity", "diagnostic_contrast", "interval", "summary_binding"])
def test_corrupt_counts_controls_and_scopes_are_rejected(problem, inputs):
    summary, snapshot = inputs
    row = snapshot["methods"][replay.METHODS[0]][0]
    if problem == "native_scope":
        snapshot["raw_reference_or_prediction_memberships_included"] = True
    elif problem == "extra_identifiers":
        snapshot["gene_ids"] = ["not permitted"]
    elif problem == "method":
        del snapshot["methods"][replay.METHODS[3]]
    elif problem == "columns":
        snapshot["columns"][1] = "fp"
    elif problem == "missing_row":
        snapshot["methods"][replay.METHODS[0]].pop()
    elif problem == "genes":
        snapshot["pillar_sizes"][0] += 1
    elif problem == "float_count":
        row[0] = float(row[0])
    elif problem == "bool_count":
        row[4] = bool(row[4])
    elif problem == "negative":
        row[1] = -1
    elif problem == "truth_count":
        row[2] += 1
    elif problem == "half_count":
        row[1] += 1
    elif problem == "covered":
        row[3] = snapshot["pillar_sizes"][0] + 1
    elif problem == "exact":
        row[4] = 1 - row[4]
    elif problem == "point":
        summary["scores"][replay.METHODS[0]]["metrics"]["f1"] += .001
    elif problem == "independent_counts":
        summary["independent_enumerated_counts"][replay.METHODS[0]]["tp"] += 1
    elif problem == "coverage":
        summary["scores"][replay.METHODS[0]]["exact_reference_groups"] += 1
    elif problem == "projected":
        summary["scores"][replay.METHODS[0]]["predicted_groups_with_scored_genes"] += 1
    elif problem == "seed":
        summary["uncertainty"]["seed"] += 1
    elif problem == "multiplicity":
        summary["uncertainty"]["multiplicity_count"] = 3
    elif problem == "diagnostic_contrast":
        summary["uncertainty"]["comparisons"][replay.METHODS[3]] = {}
    elif problem == "interval":
        summary["uncertainty"]["comparisons"][replay.METHODS[0]]["f1"]["bonferroni_ci"] = [1, -1]
    else:
        snapshot["summary"]["sha256"] = "0" * 64
    with pytest.raises(ValueError):
        replay.verify_points(summary, snapshot)


def test_integer_half_counts_and_zero_ratios():
    assert replay.statistics([0, 0, 0]).tolist() == [0, 0, 0]
    assert replay.statistics([2, 1, 2]).tolist() == pytest.approx([4 / 7, 2 / 3, 1 / 2])


def test_paired_draws_match_explicit_toy_resampling_and_batch_invariance():
    arrays = {"base": np.array([[2, 1, 0], [0, 1, 2], [0, 0, 0]]),
              "other": np.array([[2, 0, 0], [2, 0, 0], [0, 0, 0]])}
    sample = replay.paired_draws(arrays, 257, 31, 17)
    other_batch = replay.paired_draws(arrays, 257, 31, 128)
    weights = np.random.Generator(np.random.PCG64(31)).multinomial(3, [1 / 3] * 3, size=257)
    for method, rows in arrays.items():
        expected = []
        for w in weights:
            tp, fp, fn = [sum(int(w[i]) * int(rows[i, j]) for i in range(3)) for j in range(3)]
            expected.append([2 * tp / (2 * tp + fp + fn) if 2 * tp + fp + fn else 0,
                             tp / (tp + fp) if tp + fp else 0,
                             tp / (tp + fn) if tp + fn else 0])
        np.testing.assert_array_equal(sample[method], expected)
        np.testing.assert_array_equal(sample[method], other_batch[method])


@pytest.mark.parametrize("batch", [0, 129, True])
def test_unbounded_or_invalid_batch_rejected(batch):
    with pytest.raises(ValueError, match="bounded batch"):
        replay.paired_draws({"a": np.array([[0, 0, 0]])}, 100, 1, batch)


def test_all_six_interval_fields_are_checked_without_repeating_full_panel(inputs, monkeypatch):
    summary, snapshot = inputs
    arrays = replay.verify_points(summary, snapshot)
    draws = {m: np.repeat(replay.statistics(a.sum(axis=0))[None, :], 10, axis=0) for m, a in arrays.items() if m in replay.INFERENTIAL}
    monkeypatch.setattr(replay, "paired_draws", lambda *a: draws)
    fixture = copy.deepcopy(summary)
    for method in replay.METHODS[:2]:
        for i, metric in enumerate(replay.METRICS):
            delta = float(100 * (draws[method][0, i] - draws[replay.METHODS[2]][0, i]))
            fixture["uncertainty"]["comparisons"][method][metric] = dict(
                difference_percentage_points=delta, paired_95_percent_ci=[delta, delta], bonferroni_ci=[delta, delta])
    assert replay.verify_intervals(fixture, arrays) == 6
    for method in replay.METHODS[:2]:
        for metric in replay.METRICS:
            for field in ("difference_percentage_points", "paired_95_percent_ci", "bonferroni_ci"):
                changed = copy.deepcopy(fixture)
                item = changed["uncertainty"]["comparisons"][method][metric]
                if field == "difference_percentage_points":
                    item[field] += .01
                else:
                    item[field][0] -= .01
                with pytest.raises(ValueError, match="arithmetic differs"):
                    replay.verify_intervals(changed, arrays)


def test_wrong_digest_retains_failure_before_resampling(tmp_path, monkeypatch):
    output = tmp_path / "failed.json"
    monkeypatch.setattr(replay, "paired_draws", lambda *a: pytest.fail("must not resample"))
    with pytest.raises(ValueError, match="digest differs"):
        replay.reproduce(SNAPSHOT, "0" * 64, SUMMARY, output)
    assert json.loads(output.read_text())["status"] == "validation_failed"


@pytest.mark.parametrize("kind", ["export", "reproduce"])
def test_no_overwrite_before_reading(tmp_path, monkeypatch, kind):
    output = tmp_path / "existing"
    output.write_text("preserved")
    monkeypatch.setattr(replay, "read_pinned", lambda *a: pytest.fail("must not read"))
    with pytest.raises(FileExistsError):
        if kind == "export":
            replay.export_counts(Path("missing"), Path("missing"), output)
        else:
            replay.reproduce(Path("missing"), "0" * 64, Path("missing"), output)
    assert output.read_text() == "preserved"


@pytest.mark.parametrize("name", ["PUBLICATION_MAIN_TEXT_20260927.md", "PUBLICATION_MANUSCRIPT_DRAFT_20260916.md"])
def test_manuscript_availability_links_replay_with_its_limits(name):
    text = (ROOT / "benchmark_tools/results" / name).read_text()
    assert "[standalone YGOB arithmetic replay](YGOB_ARITHMETIC_REPLAY_RESULT_20261002.md)" in text
    assert "This is not native re-admission, proof of pillar exchangeability or independent" in text
    assert "family validation; original transfer scores and overlap limitations remain." in text
