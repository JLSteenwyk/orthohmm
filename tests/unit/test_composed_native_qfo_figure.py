"""Synthetic final-cell handoff; real count, render and asset decoding kernels."""

from copy import deepcopy
import csv
import json
from pathlib import Path
import shutil
import sys

import pytest

from benchmark_tools import bootstrap_composed_native_qfo_swiss as bootstrap
from benchmark_tools import plot_composed_native_qfo_scores as plot
from benchmark_tools import review_composed_native_qfo_figure as reviewer
from benchmark_tools.prepare_ob_candidate_neighborhood import record
from tests.unit.test_export_native_qfo_factorial_scores import write

ROOT = Path(__file__).resolve().parents[2]


@pytest.fixture(scope="module")
def data():
    snapshot = json.loads((ROOT / "benchmark_tools/results/native_qfo_scientific_scores_20261007_v3/report.json").read_text())
    failed = json.loads((ROOT / "benchmark_tools/results/native11_qfo_scoring_failure_addendum_20261008_v1/report.json").read_text())["row"]
    snapshot["rows"][5].update(status="retained_composed_scoring_failure", scoring_status="OUT_OF_MEMORY",
        input_accessions=failed["input_accessions"], relation_accessions=failed["relation_accessions"],
        submitted_pairs=failed["submitted_pairs"], relation_coverage=failed["relation_coverage"])
    last = deepcopy(snapshot["rows"][4])
    last.update(index=12, cell="p1_c1_r1", native_job_id=24036, status="supplied_composed_native_admission")
    snapshot["rows"][6] = last
    snapshot.update(schema="composed_native_qfo_scientific_reporting_v1", supplied_admissions=5,
        historical_four_cell_rows_preserved=True)
    counts = json.loads((ROOT / "benchmark_tools/results/qfo_corrected_factorial_complete_20260919/swiss_counts.json").read_text())
    native = {row["cell"]: deepcopy(row) for row in counts["cells"] if row["cell"] in plot.CELLS}
    # This invented final profile reproduces the prior profile's counts. It is
    # neither an actual native12 result nor an independently valid admission.
    native["p1_c1_r1"].update(families=deepcopy(native["p1_c0_r1"]["families"]),
        aggregate=deepcopy(native["p1_c0_r1"]["aggregate"]))
    uncertainty = bootstrap.analyze(native, counts)
    uncertainty.update(schema="composed_native_qfo_swiss_bootstrap_v1", bound_cells={row["cell"]:
        dict(admission=row["admission"], index=row["index"], native_job_id=row["native_job_id"])
        for row in snapshot["rows"] if row["accuracy_admitted"]})
    return snapshot, uncertainty


def test_all_cells_and_native_statistics_retained(data):
    rows, scores, coverage, statuses, intervals = plot.figure_data(*deepcopy(data))
    assert len(rows) == len(coverage) == len(statuses) == 7
    assert len(scores) == 42 and sum(row["value"] is not None for row in scores) == 30
    assert len(intervals) == 12
    assert [r["cell"] for r in coverage if r["relation_coverage"] is None] == ["p0_c1_r1"]
    assert all(r["value"] is None for r in scores if r["cell"] == "p1_c1_r0")
    assert any(r["cell"] == "p1_c1_r0" and not r["accuracy_admitted"] for r in coverage)


def test_positive_interval_not_forced_to_include_zero(data):
    snapshot, intervals = deepcopy(data)
    effect = next(r for r in intervals["comparisons"] if r["name"] == "R_at_P0_C0")
    value = effect["metrics"]["F1"]
    assert value["difference"] > .05
    value.update(paired_percentile_ci=[.05, .2], bonferroni_percentile_ci=[.04, .21])
    plot.figure_data(snapshot, intervals)


@pytest.mark.parametrize("change", ["schema", "extra_admission", "missing_score", "failed_score", "failed_mean",
    "f1", "denominator", "timing", "multiplicity", "seed", "draws", "cached", "imputed",
    "unknown_contrast", "interval_order", "point", "admission", "unavailable_interval"])
def test_partial_or_inconsistent_scope_rejected(data, change):
    snapshot, uncertainty = deepcopy(data)
    final = snapshot["rows"][6]
    effect = next(r for r in uncertainty["comparisons"] if r["metrics"] is not None)
    if change == "schema": snapshot["schema"] = "allocated_native_qfo_scientific_reporting_snapshot_v1"
    elif change == "extra_admission": snapshot["rows"][3]["accuracy_admitted"] = True
    elif change == "missing_score": final["scores"].pop("VGNC")
    elif change == "failed_score": snapshot["rows"][5]["scores"]["GO"] = .5
    elif change == "failed_mean": snapshot["rows"][5]["secondary_mean"] = 0
    elif change == "f1": final["scores"]["SwissTrees"] = .123
    elif change == "denominator": final["input_accessions"] += 1
    elif change == "timing": snapshot["rows"][1]["timing_eligible"] = True
    elif change == "multiplicity": uncertainty["multiplicity_endpoints"] = 12
    elif change == "seed": uncertainty["seed"] = 1
    elif change == "draws": uncertainty["new_bootstrap_draws"] = 0
    elif change == "cached": uncertainty["retained_intervals_reused"] = True
    elif change == "imputed": uncertainty["unobserved_cells_imputed"] = True
    elif change == "unknown_contrast": effect["weights"][0] += 1
    elif change == "interval_order": effect["metrics"]["F1"]["paired_percentile_ci"] = [1, 0]
    elif change == "point": uncertainty["point_estimates"]["p1_c1_r1"]["F1"] = .123
    elif change == "admission": uncertainty["bound_cells"]["p1_c1_r1"]["admission"] = {}
    else:
        missing = next(r for r in uncertainty["comparisons"] if r["metrics"] is None)
        missing["metrics"] = effect["metrics"]
    with pytest.raises(ValueError): plot.figure_data(snapshot, uncertainty)


@pytest.fixture(scope="module")
def rendered(tmp_path_factory, data):
    base = tmp_path_factory.mktemp("synthetic_composed_figure")
    snapshot, uncertainty = deepcopy(data)
    snapshot.update(source=record(ROOT / "benchmark_tools/export_composed_native_qfo_scientific_scores.py"),
        evidence=[], outputs=[])
    snapshot_ref = write(base / "snapshot.json", snapshot)
    uncertainty.update(source=record(bootstrap.__file__), snapshot=snapshot_ref, evidence=[], helpers=[])
    interval_ref = write(base / "intervals.json", uncertainty)
    output = base / "figure"
    # The production replay is separately tested. This fixture's admission is
    # deliberately synthetic; keep only the render/table/asset kernels real.
    with pytest.MonkeyPatch.context() as patch:
        patch.setattr(plot, "replay_snapshot", lambda *args: dict(python=record(sys.executable),
            replay_helper=record(ROOT / "benchmark_tools/audit_composed_native_qfo_swiss_counts.py"),
            result=dict(exact_snapshot_replay=True, admitted_cells=5)))
        plot.run(Path(snapshot_ref["path"]), snapshot_ref["sha256"], Path(interval_ref["path"]),
            interval_ref["sha256"], output, Path(sys.executable))
    return output


def copied(tmp_path, rendered):
    output = tmp_path / "figure"
    shutil.copytree(rendered, output)
    manifest = json.loads((output / "manifest.json").read_text())
    manifest["outputs"] = [record(output / Path(ref["path"]).name) for ref in manifest["outputs"]]
    write(output / "manifest.json", manifest)
    return output


def test_real_render_and_independent_asset_tables_readback(tmp_path, rendered):
    result = reviewer.review(copied(tmp_path, rendered), tmp_path / "preview.png")
    assert result["score_rows"] == 42 and result["interval_rows"] == 12
    assert result["assets_decoded"] and result["scientific_tables_matched"]
    assert result["png"]["size"] == [3000, 2400]
    assert not result["visual_review_complete"] and not result["publication_ready"]


@pytest.mark.parametrize("change", ["score", "null", "coverage", "status", "interval", "source", "scope", "png", "svg"])
def test_resealed_asset_or_table_mutations_rejected(tmp_path, rendered, change):
    output = copied(tmp_path, rendered)
    manifest = json.loads((output / "manifest.json").read_text())
    if change in ("source", "scope"):
        manifest["source" if change == "source" else "new_bootstrap_draws"] = {} if change == "source" else 100000
    elif change == "png": (output / f"{plot.STEM}.png").write_bytes(b"not an image")
    elif change == "svg":
        path = output / f"{plot.STEM}.svg"
        path.write_text(path.read_text().replace("Coverage is not accuracy", "Incorrect coverage"))
    else:
        name = {"score": "scores.tsv", "null": "scores.tsv", "coverage": "coverage.tsv",
            "status": "status.tsv", "interval": "swiss_intervals.tsv"}[change]
        path = output / name
        with path.open(newline="") as stream:
            rows = list(csv.DictReader(stream, delimiter="\t"))
        if change == "score": rows[0]["value"] = "0.123"
        elif change == "null": next(r for r in rows if r["value"] == "")["value"] = "0"
        elif change == "coverage": rows[0]["input_accessions"] = "1"
        elif change == "status": rows[5]["scoring_status"] = "COMPLETED"
        else: rows[0]["difference_pp"] = "123"
        with path.open("w", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
            writer.writeheader()
            writer.writerows(rows)
    manifest["outputs"] = [record(output / Path(ref["path"]).name) for ref in manifest["outputs"]]
    write(output / "manifest.json", manifest)
    with pytest.raises((ValueError, reviewer.Image.UnidentifiedImageError)):
        reviewer.review(output, tmp_path / "preview.png")


def test_existing_output_or_preview_refused_before_reads(tmp_path):
    existing = tmp_path / "existing"
    existing.mkdir()
    with pytest.raises(ValueError, match="Output already exists"):
        plot.run(Path("absent"), "unused", Path("absent"), "unused", existing, Path("absent"))
    with pytest.raises(ValueError, match="Preview already exists"):
        reviewer.review(Path("absent"), existing)


@pytest.mark.parametrize("change", ["snapshot_source", "interval_source", "cross_snapshot"])
def test_changed_input_binding_refused_before_render(tmp_path, data, change, monkeypatch):
    snapshot, uncertainty = deepcopy(data)
    snapshot["source"] = record(ROOT / "benchmark_tools/export_composed_native_qfo_scientific_scores.py")
    if change == "snapshot_source": snapshot["source"] = {}
    snapshot_ref = write(tmp_path / "snapshot.json", snapshot)
    uncertainty.update(source=record(bootstrap.__file__), snapshot=snapshot_ref)
    if change == "interval_source": uncertainty["source"] = {}
    if change == "cross_snapshot": uncertainty["snapshot"] = {}
    ref = write(tmp_path / "intervals.json", uncertainty)
    monkeypatch.setattr(plot, "render", lambda *args: pytest.fail("unexpected render"))
    with pytest.raises(ValueError, match="source or cross-snapshot"):
        plot.run(Path(snapshot_ref["path"]), snapshot_ref["sha256"], Path(ref["path"]), ref["sha256"],
            tmp_path / "figure", Path("absent"))
    assert not (tmp_path / "figure").exists()


@pytest.mark.parametrize("status,stdout", [(1, ""), (0, "not JSON"), (0, '{"admitted_cells":4}')])
def test_failed_or_partial_scientific_replay_refused(tmp_path, monkeypatch, status, stdout):
    from types import SimpleNamespace
    monkeypatch.setattr(plot.subprocess, "run", lambda *args, **kwargs:
        SimpleNamespace(returncode=status, stdout=stdout, stderr="retained failure"))
    with pytest.raises(ValueError): plot.replay_snapshot(tmp_path / "python", tmp_path / "snapshot", "unused")


def test_scientific_replay_sanitizes_environment_and_checks_five_admissions(tmp_path, monkeypatch):
    from types import SimpleNamespace
    python = tmp_path / "python"
    python.write_bytes(b"synthetic executable record; subprocess stubbed\n")
    calls = []
    for key in ("PYTHONPATH", "PYTHONHOME", "PYTHONUSERBASE", "LD_PRELOAD", "LD_LIBRARY_PATH", "LD_AUDIT"):
        monkeypatch.setenv(key, "injected")
    def observed(command, **kwargs):
        calls.append((command, kwargs))
        return SimpleNamespace(returncode=0, stdout=json.dumps(dict(exact_snapshot_replay=True, admitted_cells=5)), stderr="")
    monkeypatch.setattr(plot.subprocess, "run", observed)
    value = plot.replay_snapshot(python, tmp_path / "snapshot", "unused")
    command, kwargs = calls[0]
    assert command[1:3] == ["-I", "-B"]
    assert all(key not in kwargs["env"] for key in ("PYTHONPATH", "PYTHONHOME", "PYTHONUSERBASE",
        "LD_PRELOAD", "LD_LIBRARY_PATH", "LD_AUDIT"))
    assert kwargs["env"]["OPENBLAS_NUM_THREADS"] == kwargs["env"]["PYTHONNOUSERSITE"] == "1"
    assert value["result"]["admitted_cells"] == 5


def test_genuine_zero_similarity_is_not_missing(data):
    snapshot, uncertainty = deepcopy(data)
    snapshot["rows"][6]["scores"]["GO"] = 0
    snapshot["rows"][6]["secondary_mean"] = sum(snapshot["rows"][6]["scores"].values()) / 6
    _, scores, _, _, _ = plot.figure_data(snapshot, uncertainty)
    row = next(r for r in scores if r["cell"] == "p1_c1_r1" and r["endpoint"] == "GO")
    assert row["value"] == 0 and row["accuracy_admitted"]
