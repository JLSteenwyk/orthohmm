import copy
from fractions import Fraction

import pytest

from benchmark_tools import audit_fas_stratum_weights as audit
from benchmark_tools.audit_fas_stratum_weights import METHODS, decompose, panel, render


def fixture(p=4, m=2, k=4, c=2, r=1, pre_mean=.75, new_mean=.25):
    sample = dict(method="example", strata={
        "precomputed": dict(population=p, requested=k, saved=k, mean=pre_mean),
        "missing": dict(population=m, requested=c, saved=r, mean=new_mean)},
        omitted_requested_new_scores=c-r, intended_sample_size=k+c,
        saved_mean=(k*pre_mean+r*new_mean)/(k+r))
    population = dict(method="example", precomputed=p, missing=m, eligible_pairs=p+m,
                      native_counts_match=True, saved_lookup_strata_and_values_match=True)
    return sample, population


def test_exact_mixture_and_attrition_decomposition():
    value = decompose(*fixture())
    assert value["eligible_precomputed_weight"] == 2/3
    assert value["saved_precomputed_weight"] == .8
    assert value["rounding_component"] == 0
    assert value["attrition_weight_shift"] == float(Fraction(2, 15))
    assert value["attrition_component"] == pytest.approx(1/15)
    assert value["reweighted_saved_strata_diagnostic"] == pytest.approx(7/12)
    assert value["native_saved_mean"] == .65
    assert value["rounding_half_unit_condition_met"] is True


def test_complete_census_without_omissions_has_no_weight_shift():
    value = decompose(*fixture(r=2))
    assert value["native_minus_diagnostic"] == 0
    assert value["weight_only_absolute_difference_bound"] == 0


def test_rounding_and_attrition_have_separate_signs():
    value = decompose(*fixture(p=2, m=30000, k=1, c=9000, r=8999))
    assert value["rounding_component"] > 0
    assert value["attrition_component"] > 0
    assert value["native_minus_diagnostic"] == pytest.approx(
        value["rounding_component"] + value["attrition_component"])


def test_rounding_down_and_reversed_stratum_gap():
    value = decompose(*fixture(p=8, m=30000, k=2, c=9000, r=8999, pre_mean=.1, new_mean=.9))
    assert value["rounding_component"] > 0
    assert value["attrition_component"] < 0
    assert abs(value["native_minus_diagnostic"]) <= value["weight_only_absolute_difference_bound"]


@pytest.mark.parametrize("stratum,key,value", [
    ("precomputed", "saved", 3), ("precomputed", "requested", 5),
    ("precomputed", "population", True), ("missing", "saved", 0),
    ("missing", "saved", 3), ("missing", "requested", 1),
    ("missing", "mean", float("nan")), ("missing", "mean", float("inf")),
    ("missing", "mean", -.1), ("missing", "mean", 1.1), ("missing", "mean", True),
])
def test_invalid_stratum_rejected(stratum, key, value):
    sample, population = fixture()
    sample["strata"][stratum][key] = value
    with pytest.raises(ValueError):
        decompose(sample, population)


@pytest.mark.parametrize("key,value", [
    ("method", "other"), ("precomputed", 5), ("eligible_pairs", 7),
    ("native_counts_match", False), ("saved_lookup_strata_and_values_match", False),
])
def test_population_disagreement_rejected(key, value):
    sample, population = fixture()
    population[key] = value
    with pytest.raises(ValueError):
        decompose(sample, population)


@pytest.mark.parametrize("key,value", [
    ("saved_mean", .7), ("omitted_requested_new_scores", 2),
    ("omitted_requested_new_scores", True), ("intended_sample_size", 7),
])
def test_saved_endpoint_or_totals_disagreement_rejected(key, value):
    sample, population = fixture()
    sample[key] = value
    with pytest.raises(ValueError):
        decompose(sample, population)


def reports():
    rows = [fixture() for _ in METHODS]
    for name, (sample, population) in zip(METHODS, rows):
        sample["method"] = population["method"] = name
    attrition = dict(status="fas_requested_sample_attrition_bounded",
        methods=[a for a, _ in rows], benchmark_scores_changed=False, uncertainty_admitted=False)
    population = dict(status="retained_fas_eligible_populations_completed_with_reuse",
        methods=[p for _, p in rows], benchmark_scores_changed=False, uncertainty_admitted=False,
        publication_ready=False)
    return attrition, population


def test_complete_panel_and_qualified_table():
    values = panel(*reports())
    assert tuple(row["method"] for row in values) == METHODS
    text = render(dict(methods=values))
    assert "not a corrected benchmark" in text
    assert "confidence interval" in text
    assert "+6.666667e-02" in text


@pytest.mark.parametrize("index,change", [
    (0, "missing"), (1, "missing"), (0, "reversed"), (1, "reversed"),
    (0, "duplicate"), (1, "duplicate"), (0, "promoted"), (1, "promoted"),
])
def test_partial_reordered_duplicate_or_promoted_panel_rejected(index, change):
    values = list(copy.deepcopy(reports()))
    report = values[index]
    if change == "missing":
        report["methods"].pop()
    elif change == "reversed":
        report["methods"].reverse()
    elif change == "duplicate":
        report["methods"][1] = report["methods"][0]
    else:
        report["uncertainty_admitted"] = True
    with pytest.raises(ValueError):
        panel(*values)


def test_identical_stratum_means_can_hide_weight_shift():
    value = decompose(*fixture(pre_mean=.5, new_mean=.5))
    assert value["attrition_weight_shift"] > 0
    assert value["native_minus_diagnostic"] == 0


def test_no_mutation_of_retained_inputs():
    values = reports()
    before = copy.deepcopy(values)
    panel(*values)
    assert values == before


def test_changed_retained_report_pin_rejected(tmp_path):
    path = tmp_path / audit.INPUTS[0][0]
    path.parent.mkdir(parents=True)
    path.write_text("{}")
    with pytest.raises(ValueError, match="identity changed"):
        audit.run(tmp_path)


def test_aggregate_run_checks_inputs_twice_without_large_inputs(tmp_path, monkeypatch):
    import hashlib
    import json

    inputs = []
    for i, report in enumerate(reports()):
        path = tmp_path / (str(i) + ".json")
        path.write_text(json.dumps(report))
        inputs.append((path.name, hashlib.sha256(path.read_bytes()).hexdigest()))
    monkeypatch.setattr(audit, "INPUTS", inputs)
    checked = []
    original = audit.check
    def check(ref):
        original(ref)
        checked.append(ref)
    monkeypatch.setattr(audit, "check", check)
    result = audit.run(tmp_path)
    assert result["status"] == "retained_fas_stratum_weight_decomposition_verified"
    assert result["benchmark_scores_changed"] is False
    assert result["uncertainty_admitted"] is False
    assert result["publication_ready"] is False
    assert checked == result["checked_records"]
    assert len(checked) == 4
    assert not any(ref["path"].endswith(".db") for ref in checked)


@pytest.mark.parametrize("same", [False, True])
def test_cli_requires_fresh_distinct_outputs(tmp_path, monkeypatch, same):
    import sys

    output, table = tmp_path / "report.json", tmp_path / "table.md"
    output.write_text("retained")
    if same:
        table = output
    monkeypatch.setattr(sys, "argv", ["audit", "--output", str(output), "--table", str(table)])
    monkeypatch.setattr(audit, "run", lambda *_: pytest.fail("Read inputs before output gate"))
    with pytest.raises((ValueError, FileExistsError)):
        audit.main()
    assert output.read_text() == "retained"


def test_cli_leaves_no_report_when_input_checks_fail(tmp_path, monkeypatch):
    import sys

    output, table = tmp_path / "report.json", tmp_path / "table.md"
    monkeypatch.setattr(sys, "argv", ["audit", "--repo", str(tmp_path),
                                     "--output", str(output), "--table", str(table)])
    with pytest.raises(FileNotFoundError):
        audit.main()
    assert not output.exists() and not table.exists()
