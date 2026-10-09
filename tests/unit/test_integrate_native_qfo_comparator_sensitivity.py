"""All endpoints and unchanged parent manuscript survive statistical integration."""

import copy
import csv
import io
import json
from pathlib import Path
import shutil

import pytest

from benchmark_tools import integrate_native_qfo_comparator_sensitivity as current


ROOT = Path(__file__).resolve().parents[2]
DIRECTORY = ROOT / "benchmark_tools/results"


def result():
    return json.loads((DIRECTORY / current.PINS["result"][0]).read_text())


def test_all_planned_endpoints_and_exact_retained_values():
    data = result()
    before = copy.deepcopy(data)
    rows = current.endpoints(data)
    assert data == before
    assert len(rows) == 48 and sum(r["difference"] is None for r in rows) == 24
    for row in rows:
        comparison = next(c for c in data["comparisons"]
                          if (c["candidate"], c["reference"]) == (row["candidate"], row["reference"]))
        if row["difference"] is not None:
            metric = comparison["metrics"][row["metric"]]
            assert row["difference"] == metric["difference"]
            assert [row["adjusted_low"], row["adjusted_high"]] == metric["bonferroni_percentile_ci"]
        else:
            assert row["missing_reason"] == comparison["missing_reason"]


@pytest.mark.parametrize("change", ["schema", "scope", "seed", "draws", "quantile", "observed",
    "ready", "independent", "duplicate", "missing_imputation", "metric", "nan", "interval", "outcomes"])
def test_changed_scope_and_plot_numbers_refused(change):
    data = result()
    value = data["comparisons"][0]["metrics"]["F1"]
    if change == "schema": data["schema"] = "old_scope"
    elif change == "scope": data["multiplicity_endpoints"] = 24
    elif change == "seed": data["seed"] += 1
    elif change == "draws": data["new_bootstrap_draws"] = 0
    elif change == "quantile": data["quantile_method"] = "nearest"
    elif change == "observed": data["observed_cells"].pop()
    elif change == "ready": data["publication_ready"] = True
    elif change == "independent": data["independent_confirmation"] = True
    elif change == "duplicate": data["comparisons"][1] = copy.deepcopy(data["comparisons"][0])
    elif change == "missing_imputation": data["comparisons"][6]["metrics"] = {}
    elif change == "metric": data["comparisons"][0]["metrics"].pop("TPR")
    elif change == "nan": value["difference"] = float("nan")
    elif change == "interval": value["bonferroni_percentile_ci"].reverse()
    elif change == "outcomes": value["family_wins"] += 1
    with pytest.raises(ValueError): current.endpoints(data)


def test_manuscript_insertions_restore_exact_parent_and_scoped_statements():
    parent = (DIRECTORY / current.PINS["parent"][0]).read_text()
    revised, methods, results = current.manuscript(parent, result(), "comparison_figure")
    assert revised.replace(methods, "", 1).replace(results, "", 1) == parent
    assert revised.count("### Retrospective Native-Cell Comparator Protocol") == 1
    assert revised.count("### Retrospective Native-Cell Comparator Sensitivity") == 1
    assert "not this separately specified retrospective calculation" in results
    assert "563 disjoint represented proteins" in methods
    assert "0.05/96" in methods and "24 of the\n48 endpoints estimable" in methods
    assert "not globally\nmissing" in results
    f1_rows = [line for line in results.splitlines() if line.startswith("| P")]
    assert len(f1_rows) == 8
    for comparison in result()["comparisons"]:
        if comparison["metrics"] is not None:
            metric = comparison["metrics"]["F1"]
            assert f"{metric['difference']:+.6f}" in results
            assert f"[{metric['bonferroni_percentile_ci'][0]:.6f}, {metric['bonferroni_percentile_ci'][1]:.6f}]" in results


def test_ambiguous_parent_anchor_refused():
    parent = (DIRECTORY / current.PINS["parent"][0]).read_text()
    with pytest.raises(ValueError): current.manuscript(parent + current.METHOD_ANCHOR, result(), "figure")


def test_output_existing_refused_before_loading(tmp_path):
    with pytest.raises(FileExistsError): current.run(ROOT, tmp_path, tmp_path / "new.md")


def test_scoped_fresh_generation_with_real_figure_not_new_bootstrap(tmp_path):
    directory = tmp_path / "benchmark_tools/results"
    directory.mkdir(parents=True)
    for name, _ in current.PINS.values():
        target = directory / name
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(DIRECTORY / name, target)
    output = directory / "figure"
    revised_path = directory / "revised.md"
    manifest = current.run(tmp_path, output, revised_path)
    assert manifest["planned_endpoints"] == 48 and manifest["estimated_endpoints"] == 24
    assert manifest["no_new_bootstrap"] is True
    assert manifest["visual_reviewed"] is manifest["manuscript_rendered"] is False
    assert (output / "native_comparator_sensitivity.png").stat().st_size > 10000
    assert (output / "native_comparator_sensitivity.pdf").read_bytes().startswith(b"%PDF")
    rows = list(csv.DictReader(io.StringIO((output / "endpoints.tsv").read_text()), delimiter="\t"))
    assert len(rows) == 48 and sum(r["difference"] == "" for r in rows) == 24
    methods = (output / "methods_section.md").read_text()
    results = (output / "results_section.md").read_text()
    assert revised_path.read_text().replace(methods, "", 1).replace(results, "", 1) == (directory / current.PINS["parent"][0]).read_text()
