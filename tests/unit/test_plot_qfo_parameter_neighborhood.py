import copy
import csv
import json
from pathlib import Path

import matplotlib.pyplot as plt
import pytest

from benchmark_tools import plot_qfo_parameter_neighborhood as module
from benchmark_tools import bootstrap_qfo_parameter_neighborhood as kernel
from benchmark_tools.prepare_ob_candidate_neighborhood import record
from benchmark_tools.validate_qfo_native_assessment import AXES, validate_records
from tests.unit.test_audit_qfo_parameter_swiss import fixture as family_fixture


def fixture(tmp_path, missing=()):
    entries, baseline = family_fixture(tmp_path)
    for index in missing:
        entries[index] = dict(arm=module.ARMS[index], status="not_admitted", reason="failed or pending; retained")
    counts = module.legacy.assemble(entries, baseline)
    report = kernel.calculate(counts)
    estimated = sum(row["status"] == "estimated" for row in report["comparisons"])
    root = Path(module.__file__).resolve().parents[1]
    report.update(status="corrected_qfo_parameter_uncertainty_audited", reconstructed_counts=counts,
        scientific_inputs_admitted=True, uncertainty_admitted=estimated > 0, estimated_contrasts=estimated,
        complete_panel=estimated == 6, protocol=record(root / "benchmark_tools/results/QFO_PARAMETER_NEIGHBORHOOD_PROTOCOL_20260919.md"),
        plan=record(root / "benchmark_tools/results/qfo_parameter_neighborhood_plan_20260919.json"),
        source=record(root / "benchmark_tools/run_qfo_parameter_uncertainty.py"), helpers=[], checked_inputs=[])
    reproduced = dict(status="qfo_parameter_uncertainty_numerically_reproduced", endpoints=estimated * 3,
        planned_endpoints=18, absolute_tolerance=1e-12, publication_ready=False,
        source=record(root / "benchmark_tools/reproduce_qfo_parameter_uncertainty.py"))
    return report, reproduced, entries


@pytest.mark.parametrize("missing", [(), (2,), (0,), (1, 2, 3, 4, 5, 6)])
def test_all_planned_rows_and_missing_effects_visible(tmp_path, missing):
    report, reproduced, _ = fixture(tmp_path, missing)
    figure = module.plot(report, reproduced)
    figure.canvas.draw()
    assert len(figure.axes) == 4
    assert len(figure.axes[0].tables[0].get_celld()) == 32
    unavailable = sum(row["status"] == "not_estimable" for row in report["comparisons"])
    for ax in figure.axes[1:]:
        assert len(ax.lines) == 1 + 3 * (6 - unavailable)
        assert len(ax.texts) == unavailable
    assert all(len(ax.get_yticks()) == 6 for ax in figure.axes[1:])
    assert any("18 planned endpoints" in text.get_text() for text in figure.texts)
    assert len(figure.canvas.buffer_rgba()) > 0
    plt.close(figure)


@pytest.mark.parametrize("problem", ["status", "scientific", "publication", "protocol", "plan", "source", "seed",
    "replicates", "multiplicity", "quantile", "point", "interval", "effect", "directions", "order", "complete",
    "reproduction", "reproduction_source", "reproduction_endpoints", "missing_imputed"])
def test_invalid_or_promoted_plot_evidence(tmp_path, problem):
    report, reproduced, _ = fixture(tmp_path, (2,) if problem == "missing_imputed" else ())
    if problem in ("status", "scientific", "publication", "seed", "replicates", "multiplicity", "quantile", "complete"):
        key, value = {"status": ("status", "unadmitted"), "scientific": ("scientific_inputs_admitted", False),
            "publication": ("publication_ready", True), "seed": ("seed", 1), "replicates": ("replicates", 1000),
            "multiplicity": ("multiplicity_endpoints", 15), "quantile": ("quantile_method", "nearest"),
            "complete": ("complete_panel", False)}[problem]
        report[key] = value
    elif problem in ("protocol", "plan", "source"):
        report[problem]["sha256"] = "wrong"
    elif problem == "point":
        report["point_estimates"]["control"]["F1"] = .99
    elif problem == "interval":
        report["comparisons"][0]["metrics"]["F1"]["bonferroni_percentile_ci"] = [1, -1]
    elif problem == "effect":
        report["comparisons"][0]["metrics"]["F1"]["difference"] = .9
    elif problem == "directions":
        report["comparisons"][0]["metrics"]["F1"]["family_wins"] = True
    elif problem == "order":
        report["comparisons"].reverse()
    elif problem == "reproduction":
        reproduced["status"] = "pending"
    elif problem == "reproduction_source":
        reproduced["source"]["sha256"] = "wrong"
    elif problem == "reproduction_endpoints":
        reproduced["endpoints"] = 15
    else:
        report["point_estimates"]["cpm_high"] = {key: 0. for key in module.METRICS}
    with pytest.raises(ValueError):
        module.validate(report, reproduced)


def setup_export(tmp_path, monkeypatch):
    report, reproduced, entries = fixture(tmp_path, (2,))
    inventory = dict(status="qfo_parameter_score_admission_inventory", arms=[])
    for entry in entries:
        if entry["status"] == "not_admitted":
            inventory["arms"].append(copy.deepcopy(entry))
            continue
        assessment = entry["assessment"]
        native = validate_records(assessment["native_assessments"], assessment["participant"], set(report["families"]))
        assessment["endpoints"] = {}
        for challenge in module.ENDPOINTS:
            x, y = [native[challenge, axis]["metrics"]["value"] for axis in AXES[challenge]]
            score = 2 * x * y / (x + y) if challenge in ("VGNC", "SwissTrees", "TreeFam-A") else y
            assessment["endpoints"][challenge] = dict(score=score,
                axes=dict(x_axis=AXES[challenge][0], y_axis=AXES[challenge][1]),
                native_participant=dict(participant_id=assessment["participant"], metric_x=x, metric_y=y))
        assessment["secondary_six_metric_mean"] = sum(row["score"] for row in assessment["endpoints"].values()) / 6
        pair = tmp_path / (entry["arm"] + "_pairs.json")
        pair.write_text("{}")
        admission = tmp_path / (entry["arm"] + "_admission.json")
        admission.write_text(json.dumps(dict(status="fixture_admitted", pairs_manifest=record(pair), assessment=assessment)))
        inventory["arms"].append(dict(arm=entry["arm"], status="admitted", admission=record(admission)))
    inv_path = tmp_path / "inventory.json"
    inv_path.write_text(json.dumps(inventory))
    report["admission_inventory"] = record(inv_path)
    path = tmp_path / "uncertainty.json"
    path.write_text(json.dumps(report))
    reproduced["input"] = record(path)
    repro_path = tmp_path / "reproduction.json"
    repro_path.write_text(json.dumps(reproduced))
    seen = []
    monkeypatch.setattr(module.legacy, "validate_admission", lambda arm, *args: seen.append(arm))
    return (path, record(path)["sha256"], repro_path, record(repro_path)["sha256"], tmp_path / "export"), seen


def test_tables_native_arithmetic_and_six_file_figure_export(tmp_path, monkeypatch):
    args, seen = setup_export(tmp_path, monkeypatch)
    result = module.export(*args)
    assert seen == [arm for arm in module.ARMS if arm != "cpm_high"]
    assert result["estimated_endpoints"] == 15 and result["planned_endpoints"] == 18
    assert result["complete_panel"] is False and result["publication_ready"] is False
    assert len(result["outputs"]) == 6 and all(Path(ref["path"]).stat().st_size > 100 for ref in result["outputs"])
    missing = result["rows"][2]
    assert all(missing[key] is None for key in (*module.ENDPOINTS, "secondary_mean"))
    assert len(result["intervals"]) == 18
    assert all(row["difference_pp"] is None for row in result["intervals"] if row["candidate"] == "cpm_high")
    with (args[-1] / "intervals.tsv").open() as stream:
        assert len(list(csv.DictReader(stream, delimiter="\t"))) == 18
    assert "GO/EC similarity and FAS are not F1" in (args[-1] / "scores.md").read_text()
    assert json.loads((args[-1] / "manifest.json").read_text()) == result


@pytest.mark.parametrize("problem", ["wrong_reproduction", "input_pin", "output", "symlink", "source_changed"])
def test_export_binding_and_overwrite_contract(tmp_path, monkeypatch, problem):
    args, _ = setup_export(tmp_path, monkeypatch)
    args = list(args)
    if problem == "wrong_reproduction":
        reproduced = json.loads(args[2].read_text())
        reproduced["input"]["sha256"] = "wrong"
        args[2].write_text(json.dumps(reproduced))
        args[3] = record(args[2])["sha256"]
    elif problem == "input_pin":
        args[1] = "wrong"
    elif problem == "output":
        args[-1].mkdir()
    elif problem == "symlink":
        args[-1].symlink_to(tmp_path / "missing")
    else:
        original = module.native_rows
        def changed(report, checked):
            result = original(report, checked)
            args[0].write_text("changed during export")
            return result
        monkeypatch.setattr(module, "native_rows", changed)
    with pytest.raises((ValueError, FileExistsError)):
        module.export(*args)
    assert not (args[-1] / "manifest.json").exists()


@pytest.mark.parametrize("problem", ["inventory_order", "availability", "native_score", "native_metric", "private_wrong_arm"])
def test_native_table_cannot_diverge_from_admitted_uncertainty(tmp_path, monkeypatch, problem):
    args, _ = setup_export(tmp_path, monkeypatch)
    args = list(args)
    report = json.loads(args[0].read_text())
    inventory_path = Path(report["admission_inventory"]["path"])
    inventory = json.loads(inventory_path.read_text())
    if problem == "inventory_order":
        inventory["arms"].reverse()
    elif problem == "availability":
        inventory["arms"][2] = dict(arm="cpm_high", status="not_admitted", reason="different reason")
    else:
        entry = inventory["arms"][0]
        path = Path(entry["admission"]["path"])
        admitted = json.loads(path.read_text())
        if problem == "native_score":
            admitted["assessment"]["endpoints"]["SwissTrees"]["score"] = .99
        elif problem == "native_metric":
            admitted["assessment"]["endpoints"]["GO"]["native_participant"]["metric_y"] = .99
        else:
            admitted["status"] = "private_recovered_cpm_assessment_admitted"
        path.write_text(json.dumps(admitted))
        entry["admission"] = record(path)
    inventory_path.write_text(json.dumps(inventory))
    report["admission_inventory"] = record(inventory_path)
    args[0].write_text(json.dumps(report))
    args[1] = record(args[0])["sha256"]
    reproduced = json.loads(args[2].read_text())
    reproduced["input"] = record(args[0])
    args[2].write_text(json.dumps(reproduced))
    args[3] = record(args[2])["sha256"]
    with pytest.raises(ValueError):
        module.export(*args)
    assert not args[-1].exists()
