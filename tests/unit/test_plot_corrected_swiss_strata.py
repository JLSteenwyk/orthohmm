import csv
import json

import matplotlib.pyplot as plt
import pytest

from benchmark_tools import plot_corrected_swiss_strata as module
from benchmark_tools.bootstrap_corrected_swiss_strata import bootstrap
from tests.unit.test_bootstrap_corrected_swiss_strata import fixture


def synthetic_report(missing=False):
    counts, membership = fixture()
    if missing:
        membership = {f: "missing" for f in membership}
    report = bootstrap(counts, membership, replicates=100)
    # Synthetic schema fixture only; not a scientific admission or resampling result.
    report.update(status="corrected_swiss_primary_stratified_intervals", replicates=100000,
                  scientific_inputs_admitted=True, uncertainty_admitted=True)
    return report


def test_all_endpoints_and_interactions_preserved():
    report = synthetic_report()
    rows = module.endpoints(report)
    assert len(rows) == 27
    for row in rows:
        comparisons = (report["interactions"] if row["stratum"] == "interaction"
                       else report["bins"][row["stratum"]]["comparisons"])
        original = next(c for c in comparisons if (c["candidate"], c["reference"]) ==
                        (row["candidate"], row["reference"]))["metrics"][row["metric"]]
        assert row["difference"] == original["difference"]
        assert row["nominal"] == original["paired_percentile_ci"]
        assert row["adjusted"] == original["bonferroni_percentile_ci"]


@pytest.mark.parametrize("missing", [False, True])
def test_plot_preserves_units_and_unavailable_intervals(missing):
    report = synthetic_report(missing)
    fig, rows = module.plot(report)
    try:
        assert len(fig.axes) == 3
        for ax, metric in zip(fig.axes, module.METRICS):
            selected = [r for r in rows if r["metric"] == metric]
            dots = [l for l in ax.lines if l.get_marker() == "o"]
            if missing:
                assert not dots
                assert sum(t.get_text() == "CI unavailable" for t in ax.texts) == 9
            else:
                assert [float(l.get_xdata()[0]) for l in dots] == [100*r["difference"] for r in selected]
                low, high = ax.get_xlim()
                assert all(low <= 100*v <= high for r in selected for v in r["adjusted"])
        fig.canvas.draw()
    finally:
        plt.close(fig)


@pytest.mark.parametrize("fault", ["status", "admission", "seed", "inventory", "overlap", "contrast",
                                   "metric", "nan", "bounds", "order", "missing_ci", "eligibility"])
def test_changed_or_invalid_results_rejected(fault):
    report = synthetic_report()
    value = report["bins"]["lower"]["comparisons"][0]["metrics"]["F1"]
    if fault == "status":
        report["status"] = "historical"
    elif fault == "admission":
        report["scientific_inputs_admitted"] = False
    elif fault == "seed":
        report["seed"] += 1
    elif fault == "inventory":
        report["bins"]["extra"] = report["bins"]["missing"]
    elif fault == "overlap":
        report["bins"]["higher"]["families"][0] = report["bins"]["lower"]["families"][0]
    elif fault == "contrast":
        report["interactions"].reverse()
    elif fault == "metric":
        del report["interactions"][0]["metrics"]["F1"]
    elif fault == "nan":
        value["difference"] = float("nan")
    elif fault == "bounds":
        value["difference"] = 1.1
    elif fault == "order":
        value["bonferroni_percentile_ci"] = [1, -1]
    elif fault == "missing_ci":
        value["paired_percentile_ci"] = None
    else:
        report["bins"]["lower"]["interval_eligible"] = False
    with pytest.raises(ValueError):
        module.endpoints(report)


@pytest.mark.parametrize("missing", [False, True])
def test_render_binds_source_and_exports_lossless_table(tmp_path, missing):
    source = tmp_path / "synthetic.json"
    source.write_text(json.dumps(synthetic_report(missing)))
    digest = module.record(source)["sha256"]
    output = tmp_path / "figure"
    with pytest.raises(ValueError):
        module.render(source, "0"*64, output)
    assert not output.exists()
    manifest = module.render(source, digest, output)
    assert manifest["results"]["sha256"] == digest
    assert len(manifest["outputs"]) == 4
    for item in manifest["outputs"]:
        module.check(item)
    with (output / "endpoints.tsv").open() as stream:
        rows = list(csv.DictReader(stream, delimiter="\t"))
    assert len(rows) == 27
    expected = module.endpoints(json.loads(source.read_text()))
    assert [float(r["difference"]) if r["difference"] else None for r in rows] == [r["difference"] for r in expected]
    if missing:
        assert all(not r[k] for r in rows for k in ("difference", "nominal_lower", "adjusted_lower"))
    with pytest.raises(FileExistsError):
        module.render(source, digest, output)


def test_nonempty_bin_cannot_hide_point_estimate():
    report = synthetic_report()
    report["bins"]["lower"]["comparisons"][0]["metrics"]["F1"]["difference"] = None
    with pytest.raises(ValueError, match="Point availability"):
        module.endpoints(report)
