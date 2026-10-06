import copy
import json
from pathlib import Path
from types import SimpleNamespace
import xml.etree.ElementTree as ET

import numpy as np
from PIL import Image
import pytest

from benchmark_tools import plot_native_qfo_scientific_scores as module
from benchmark_tools.prepare_ob_candidate_neighborhood import record
from tests.unit.test_export_native_qfo_factorial_scores import write

ROOT = Path(__file__).resolve().parents[2] / "benchmark_tools/results"
VALIDATOR = ROOT.parents[1] / "benchmarks/work/native_factorial_review_py310_20261004/bin/python"


def inputs():
    snapshot = json.loads((ROOT / "native_qfo_scientific_scores_20261006_v1/report.json").read_text())
    binding = json.loads((ROOT / "native_qfo_swiss_uncertainty_binding_22449_20261006.json").read_text())
    return snapshot, binding


def test_native_points_and_non_f1_endpoints_remain_distinct():
    snapshot, binding = inputs()
    rows, scores, intervals, admitted = module.figure_data(snapshot, binding)
    assert len(scores) == 12 and len(intervals) == 3 and admitted == 2
    assert sum(row["statistic"] == "F1" for row in scores) == 6
    assert all(row["endpoint"] != "Secondary mean" for row in scores)
    assert rows["p0_c0_r1"]["resources"] is None
    assert intervals[0]["adjusted_low_pp"] < 0 < intervals[0]["adjusted_high_pp"]


@pytest.mark.parametrize("value", [True, float("nan"), float("inf"), "0.5", -.1, 1.1])
def test_invalid_plotted_score_refused(value):
    snapshot, binding = inputs()
    snapshot["rows"][1]["scores"]["FAS"] = value
    with pytest.raises(ValueError):
        module.figure_data(snapshot, binding)


@pytest.mark.parametrize("key,value", [("accuracy_admitted", False), ("resources", {}),
    ("timing_admitted", True), ("timing_eligible", True)])
def test_recovered_accuracy_does_not_relabel_failed_timing(key, value):
    snapshot, binding = inputs()
    snapshot["rows"][1][key] = value
    with pytest.raises(ValueError):
        module.figure_data(snapshot, binding)


@pytest.mark.parametrize("change", ["missing_endpoint", "bad_f1", "duplicate_cell", "publication_ready",
    "unmatched_contrast", "wrong_reference", "wrong_adjustment", "missing_metric", "invalid_interval",
    "nonnested_interval", "wrong_difference"])
def test_plot_scope_and_uncertainty_contracts(change):
    snapshot, binding = inputs()
    effect = next(row for row in binding["contrasts"] if row["name"] == "R_at_P0_C0")
    if change == "missing_endpoint": snapshot["rows"][1]["scores"].pop("GO")
    elif change == "bad_f1": snapshot["rows"][1]["scores"]["SwissTrees"] = .1
    elif change == "duplicate_cell": snapshot["rows"][2]["cell"] = "p0_c0_r1"
    elif change == "publication_ready": snapshot["publication_ready"] = True
    elif change == "unmatched_contrast": effect["status"] = "native_records_unavailable"
    elif change == "wrong_reference": effect["reference"] = "p1_c0_r0"
    elif change == "wrong_adjustment": binding["multiplicity_endpoints"] = 3
    elif change == "missing_metric": effect["metrics"].pop("TPR")
    elif change == "invalid_interval": effect["metrics"]["F1"]["bonferroni_percentile_ci"] = [float("nan"), .1]
    elif change == "nonnested_interval": effect["metrics"]["F1"]["bonferroni_percentile_ci"] = [.1, .2]
    else: effect["metrics"]["F1"]["difference"] = .9
    with pytest.raises(ValueError):
        module.figure_data(snapshot, binding)


def test_existing_output_preserved_before_input_reads(tmp_path):
    output = tmp_path / "existing"
    output.mkdir()
    (output / "retain").write_text("retain me")
    with pytest.raises(ValueError, match="Output already exists"):
        module.run(Path("absent"), "absent", Path("absent"), "absent", output, Path("absent"))
    assert (output / "retain").read_text() == "retain me"


def test_cross_snapshot_binding_rejected_before_output_creation(tmp_path):
    snapshot, binding = inputs()
    snapshot_ref = write(tmp_path / "snapshot.json", snapshot)
    binding_ref = write(tmp_path / "binding.json", binding)
    output = tmp_path / "plot"
    with pytest.raises(ValueError, match="different snapshots"):
        module.run(Path(snapshot_ref["path"]), snapshot_ref["sha256"],
                   Path(binding_ref["path"]), binding_ref["sha256"], output, Path("absent"))
    assert not output.exists()


def test_render_fixture_exports_pixels_vector_labels_and_bound_endpoints(tmp_path, monkeypatch):
    snapshot, binding = inputs()
    snapshot_ref = write(tmp_path / "snapshot.json", snapshot)
    binding["snapshot"] = snapshot_ref
    binding_ref = write(tmp_path / "binding.json", binding)
    monkeypatch.setattr(module, "replay_binding", lambda *args: dict(python=record(VALIDATOR), fixture=True))
    result = module.run(Path(snapshot_ref["path"]), snapshot_ref["sha256"],
        Path(binding_ref["path"]), binding_ref["sha256"], tmp_path / "plot", VALIDATOR)
    assert result["new_bootstrap_draws"] == 0
    assert result["new_scoring_or_admission"] is result["scientific_timings_admitted"] is False
    assert result["publication_ready"] is result["visual_review_complete"] is False
    assert len(result["outputs"]) == 5
    assert all(record(ref["path"]) == ref for ref in result["outputs"])
    with Image.open(tmp_path / "plot/native_qfo_p0c0.png") as picture:
        assert picture.size == (2500, 1600)
        pixels = np.asarray(picture.convert("RGB"))
        assert np.count_nonzero(np.any(pixels < 240, axis=2)) > 5000
        for color in module.COLORS:
            rgb = [int(color[i:i+2], 16) for i in (1, 3, 5)]
            assert np.count_nonzero(np.all(pixels == rgb, axis=2)) > 100
    svg = ET.parse(tmp_path / "plot/native_qfo_p0c0.svg")
    labels = " ".join(node.text or "" for node in svg.iter() if node.tag.endswith("}text"))
    for text in ("Orthology F1", "Similarity endpoints (not F1)", "precision-recall",
                 "42-endpoint adjusted interval", "2/7 fresh cells admitted", "F1 adjusted interval includes zero"):
        assert text in labels


def fake_worker(monkeypatch, tmp_path, binding):
    python = tmp_path / "venv/bin/python"
    ref = dict(path=str(python), bytes=1, sha256=module.VALIDATION_PYTHON_SHA)
    response = dict(binding=copy.deepcopy(binding), versions=dict(module.VALIDATION_VERSIONS),
                    executable=str(python), prefix=str(python.parent.parent))
    process = SimpleNamespace(returncode=0, stderr="", stdout=json.dumps(response))
    calls = []
    monkeypatch.setattr(module, "record", lambda *args: ref)
    monkeypatch.setattr(module, "check", lambda *args: None)
    monkeypatch.setattr(module.subprocess, "run", lambda *args, **kwargs: (calls.append((args, kwargs)) or process))
    return python, ref, response, process, calls


def test_worker_preserves_venv_invocation_and_sanitizes_environment(tmp_path, monkeypatch):
    _, binding = inputs()
    python, _, _, _, calls = fake_worker(monkeypatch, tmp_path, binding)
    for name in ("PYTHONHOME", "PYTHONPATH", "LD_PRELOAD", "LD_LIBRARY_PATH", "LD_AUDIT"):
        monkeypatch.setenv(name, "untrusted-injection")
    result = module.replay_binding(python, tmp_path / "snapshot", "digest", binding)
    command = calls[0][0][0]
    assert command[:3] == [str(python), "-B", "-c"]
    assert command[3] == module.WORKER
    request = json.loads(command[4])
    assert request["versions"] == module.VALIDATION_VERSIONS
    assert request["snapshot"] == str(tmp_path / "snapshot")
    kwargs = calls[0][1]
    assert kwargs["timeout"] == 180 and kwargs["check"] is False
    assert all(name not in kwargs["env"] for name in
               ("PYTHONHOME", "PYTHONPATH", "LD_PRELOAD", "LD_LIBRARY_PATH", "LD_AUDIT"))
    assert all(kwargs["env"][name] == "1" for name in
               ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"))
    assert result["exact_binding_replay"] is True
    assert result["new_scoring_or_admission"] is False


@pytest.mark.parametrize("change", ["binary", "python_version", "package", "prefix", "executable",
                                   "binding", "worker_failed", "stderr", "malformed_json"])
def test_changed_worker_identity_or_exact_replay_refused(tmp_path, monkeypatch, change):
    _, binding = inputs()
    python, ref, response, process, calls = fake_worker(monkeypatch, tmp_path, binding)
    if change == "binary": ref["sha256"] = "0" * 64
    elif change == "python_version": response["versions"]["python"] = "3.12.3"
    elif change == "package": response["versions"]["numpy"] = "2.0.0"
    elif change == "prefix": response["prefix"] = "/other/environment"
    elif change == "executable": response["executable"] = "/symlink-resolved/python"
    elif change == "binding": response["binding"]["alpha"] = .05000000000000001
    elif change == "worker_failed": process.returncode = 1
    elif change == "stderr": process.stderr = "unexpected warning"
    process.stdout = "not-json" if change == "malformed_json" else json.dumps(response)
    with pytest.raises(ValueError):
        module.replay_binding(python, tmp_path / "snapshot", "digest", binding)
    if change == "binary": assert not calls


def test_actual_retained_binding_replays_exactly_in_original_environment():
    if not VALIDATOR.exists():
        pytest.skip("Original scientific environment not installed on this host")
    snapshot_path = ROOT / "native_qfo_scientific_scores_20261006_v1/report.json"
    snapshot, binding = inputs()
    values = [snapshot["rows"][1]["scores"][name] for name in module.reporter.original.ENDPOINTS]
    sequential = 0.0
    for value in values:
        sequential += value
    assert sequential / 6 == snapshot["rows"][1]["secondary_mean"]
    result = module.replay_binding(VALIDATOR, snapshot_path, record(snapshot_path)["sha256"], binding)
    assert result["versions"] == module.VALIDATION_VERSIONS
    assert result["exact_binding_replay"] is True
    assert result["invocation"] == str(VALIDATOR)
