from copy import deepcopy
import json

import pytest

from benchmark_tools import reproduce_qfo_parameter_uncertainty as module
from benchmark_tools.bootstrap_qfo_parameter_neighborhood import calculate
from tests.unit.test_bootstrap_qfo_parameter_neighborhood import fixture


def result_for(counts):
    report = calculate(counts)
    estimated = sum(r["status"] == "estimated" for r in report["comparisons"])
    report.update(status="corrected_qfo_parameter_uncertainty_audited", scientific_inputs_admitted=True,
                  uncertainty_admitted=estimated > 0, complete_panel=estimated == 6,
                  estimated_contrasts=estimated, reconstructed_counts=counts)
    return report


@pytest.fixture(scope="module")
def result():
    return result_for(fixture())


def test_complete_numerical_reproduction(result):
    assert module.verify(result) == 18


@pytest.mark.parametrize("missing, expected", [((1, 2), 12), (tuple(range(1, 7)), 0), ((0,), 0)])
def test_missing_arms_preserved(missing, expected):
    counts = fixture()
    for index in missing:
        counts["arms"][index] = {"arm": module.ARMS[index], "status": "not_admitted", "reason": "pending"}
    report = result_for(counts)
    assert module.verify(report) == expected
    report["comparisons"][0]["metrics"] = {}
    with pytest.raises(ValueError, match="Unavailable contrast"):
        module.verify(report)


@pytest.mark.parametrize("problem", ["seed", "multiplicity", "order", "point", "interval", "shape", "nan",
    "family", "wins", "complete", "count", "raw", "overlap", "truth", "summary"])
def test_corruption_rejected(result, problem):
    report = deepcopy(result)
    item = report["comparisons"][0]["metrics"]["F1"]
    if problem in ("seed", "multiplicity"):
        report["seed" if problem == "seed" else "multiplicity_endpoints"] += 1
    elif problem == "order":
        report["comparisons"].reverse()
    elif problem == "point":
        item["difference"] += .01
    elif problem == "interval":
        item["bonferroni_percentile_ci"][0] += .01
    elif problem == "shape":
        item["paired_percentile_ci"] = item["paired_percentile_ci"][:1]
    elif problem == "nan":
        item["difference"] = float("nan")
    elif problem == "family":
        report["comparisons"][0]["family_differences"][0]["F1"] += .01
    elif problem == "wins":
        item["family_wins"] += 1
    elif problem == "complete":
        report["complete_panel"] = False
    elif problem == "count":
        report["estimated_contrasts"] = 5
    elif problem == "raw":
        report["reconstructed_counts"]["arms"][0]["families"][0]["counts_without_prior"]["TP"] = -1
    elif problem == "overlap":
        rows = report["reconstructed_counts"]["arms"][0]["families"]
        rows[1]["represented_genes"] = rows[0]["represented_genes"]
    elif problem == "truth":
        report["reconstructed_counts"]["arms"][1]["families"][0]["counts_without_prior"]["TN"] += 1
    else:
        report["arms"][0]["status"] = "not_admitted"
    with pytest.raises(ValueError):
        module.verify(report)


def test_runner_hash_checks_and_fresh_output(result, tmp_path):
    helper = tmp_path / "helper.txt"
    helper.write_text("frozen evidence")
    report = deepcopy(result)
    report.update(source=module.record(helper), helpers=[], checked_inputs=[module.record(helper)])
    source, output = tmp_path / "input.json", tmp_path / "output.json"
    source.write_text(json.dumps(report))
    digest = module.record(source)["sha256"]
    actual = module.run(source, digest, output)
    assert actual["endpoints"] == 18 and actual["publication_ready"] is False
    with pytest.raises(FileExistsError):
        module.run(source, digest, output)
    helper.write_text("changed")
    with pytest.raises(ValueError):
        module.run(source, digest, tmp_path / "changed.json")
    assert not (tmp_path / "changed.json").exists()
