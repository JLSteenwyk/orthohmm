"""Native-count resampling tests, not production native benchmark outcomes."""

from copy import deepcopy
import sys

import numpy as np
import pytest

from benchmark_tools import bootstrap_composed_native_qfo_swiss as module
from benchmark_tools.audit_qfo_swiss_counts import statistics
from benchmark_tools.bootstrap_qfo_factorial import CELLS
from benchmark_tools.bootstrap_qfo_swiss_stages import aggregate, METRICS
from tests.unit import test_composed_native_qfo_swiss_uncertainty as handoff
from tests.unit.test_bind_native_qfo_swiss_uncertainty import retained_bootstrap
from tests.unit.test_composed_native_qfo_swiss_uncertainty import frozen


def rows(counts, indexes):
    return {CELLS[i]: deepcopy(counts["cells"][i]) for i in indexes}


def test_all_observed_counts_reproduce_original_complete_bootstrap(retained_bootstrap):
    counts, original = retained_bootstrap
    result = module.analyze(rows(counts, range(8)), counts)
    assert result["point_estimates"] == original["point_estimates"]
    for fresh, old in zip(result["comparisons"], original["comparisons"]):
        assert fresh["metrics"] == old["metrics"]
        assert fresh["family_differences"] == old["family_differences"]
    assert result["new_bootstrap_draws"] == 100000
    assert not result["retained_intervals_reused"] and not result["unobserved_cells_imputed"]
    assert result["multiplicity_endpoints"] == 42 and result["seed"] == 20260922


def test_two_observed_cells_never_fill_missing_contrasts(retained_bootstrap):
    counts, _ = retained_bootstrap
    result = module.analyze(rows(counts, (5, 7)), counts)
    assert result["observed_cells"] == ["p1_c0_r1", "p1_c1_r1"]
    assert set(result["point_estimates"]) == set(result["observed_cells"])
    available = [row for row in result["comparisons"] if row["metrics"] is not None]
    assert len(available) == 1 and available[0]["name"] == "C_at_P1_R1"
    assert len(result["comparisons"]) == 14
    for row in result["comparisons"]:
        if row["metrics"] is None:
            assert row["missing_cells"] and row["family_differences"] is None


def test_changed_actual_counts_recompute_benchmark_statistic_not_cached_f1(retained_bootstrap):
    counts, original = retained_bootstrap
    native = rows(counts, (5, 7))
    last = native["p1_c1_r1"]
    for i, family in enumerate(last["families"]):
        if i < 9:
            family["counts_without_prior"]["TP"] += i % 3 + 1
            family["counts_without_prior"]["FN"] -= i % 3 + 1
            family["statistics_with_prior"] = statistics(family["counts_without_prior"])
    means = np.mean([[r["statistics_with_prior"][m] for m in ("PPV", "TPR")]
        for r in last["families"]], axis=0)
    last["aggregate"] = dict(zip(METRICS, aggregate(means).tolist()))
    result = module.analyze(native, counts)
    observed = result["point_estimates"]["p1_c1_r1"]["F1"]
    assert observed == last["aggregate"]["F1"]
    assert abs(observed - original["point_estimates"]["p1_c1_r1"]["F1"]) > .01
    assert abs(observed - np.mean([r["statistics_with_prior"]["F1"] for r in last["families"]])) > 1e-6
    effect = next(row for row in result["comparisons"] if row["name"] == "C_at_P1_R1")
    assert effect["metrics"]["F1"]["paired_percentile_ci"][0] < effect["metrics"]["F1"]["paired_percentile_ci"][1]
    assert effect["metrics"]["F1"]["bonferroni_percentile_ci"][0] <= effect["metrics"]["F1"]["paired_percentile_ci"][0]


@pytest.mark.parametrize("problem", ["members", "truth", "negative", "boolean", "statistics", "aggregate", "family", "missing"])
def test_new_intervals_require_consistent_real_native_counts(retained_bootstrap, problem):
    counts, _ = retained_bootstrap
    native = rows(counts, (5, 7))
    row = native["p1_c1_r1"]
    family = row["families"][0]
    if problem == "members": family["represented_genes"][0] = "not_reference_member"
    elif problem == "truth": family["counts_without_prior"]["FN"] += 1
    elif problem == "negative": family["counts_without_prior"]["TP"] = -1
    elif problem == "boolean": family["counts_without_prior"]["TP"] = True
    elif problem == "statistics": family["statistics_with_prior"]["PPV"] = .123
    elif problem == "aggregate": row["aggregate"]["F1"] = .123
    elif problem == "family": family["family"] = "unknown"
    else: row["families"].pop()
    with pytest.raises(ValueError):
        module.analyze(native, counts)


@pytest.mark.parametrize("indexes", [(), (7,), (0, 7)])
def test_no_estimable_contrast_does_not_generate_draws(retained_bootstrap, indexes, monkeypatch):
    counts, _ = retained_bootstrap
    monkeypatch.setattr(module.np.random, "Generator", lambda *args: pytest.fail("unneeded draws"))
    with pytest.raises(ValueError):
        module.analyze(rows(counts, indexes), counts)


def test_checked_different_raw_counts_enable_new_intervals(tmp_path, monkeypatch, frozen):
    refs = handoff.mixed_fixture(tmp_path, monkeypatch, frozen, matching=False)
    snapshot, counts, composed, allocated, old_bootstrap, _ = refs
    result = module.run(snapshot["path"], snapshot["sha256"],
        [(r["path"], r["sha256"]) for r in (composed, allocated)], counts["path"], old_bootstrap["path"])
    assert result["schema"] == "composed_native_qfo_swiss_bootstrap_v1"
    assert result["bound_cells"]["p1_c1_r1"]["status"] == "native_records_differ"
    effect = next(row for row in result["comparisons"] if row["name"] == "C_at_P1_R1")
    assert effect["metrics"] is not None and result["new_bootstrap_draws"] == 100000
    assert not result["publication_ready"] and not result["independent_confirmation"]
    assert not result["new_accuracy_or_resource_admission"]


def test_inconsistent_reference_cannot_be_resampled_even_after_guarded_binding(tmp_path, monkeypatch, frozen):
    refs = handoff.mixed_fixture(tmp_path, monkeypatch, frozen, changed=True)
    snapshot, counts, composed, allocated, old_bootstrap, _ = refs
    with pytest.raises(ValueError, match="reference members"):
        module.run(snapshot["path"], snapshot["sha256"],
            [(r["path"], r["sha256"]) for r in (composed, allocated)], counts["path"], old_bootstrap["path"])


def test_cli_retains_existing_output_before_work(tmp_path, monkeypatch):
    output = tmp_path / "existing.json"
    output.write_text("retain me")
    monkeypatch.setattr(sys, "argv", ["test", "--snapshot", "absent", "--snapshot-sha256", "unused",
        "--retained-counts", "absent", "--bootstrap", "absent", "--counts-audit", "absent", "unused",
        "--output", str(output)])
    with pytest.raises(ValueError, match="Output already exists"):
        module.main()
    assert output.read_text() == "retain me"
