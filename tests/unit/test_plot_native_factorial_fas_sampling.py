import copy
import json
from pathlib import Path

from PIL import Image
import numpy as np
import pytest

from benchmark_tools import plot_native_factorial_fas_sampling as current


def report():
    path = Path(__file__).resolve().parents[2] / "benchmark_tools/results/native_factorial_fas_sampling_20261010_v1/report.json"
    return json.loads(path.read_text())


def test_all_actual_rows_contrasts_and_target_are_preserved():
    data = report()
    methods, contrasts = current.validate(data)
    assert len(methods) == 4 and len(contrasts) == 6
    assert sum(r["zero_included"] for r in contrasts) == 4


@pytest.mark.parametrize("change", ["target", "biological", "allocation", "inventory", "bound", "observed", "zero"])
def test_changed_scope_or_contrast_rejected(change):
    data = copy.deepcopy(report())
    if change == "target": data["target"] = "biological_accuracy"
    elif change == "biological": data["biological_generalization_intervals"] = True
    elif change == "allocation": data["component_error"] = .05
    elif change == "inventory": data["contrasts"].pop()
    elif change == "bound": data["contrasts"][0]["conditional_expected_difference_bounds"][0] += .001
    elif change == "observed": data["contrasts"][0]["observed_difference"] += .001
    elif change == "zero": data["contrasts"][0]["zero_included"] = False
    with pytest.raises(ValueError): current.validate(data)


def test_actual_render_is_nonblank_with_all_vector_labels(tmp_path):
    current.render(report(), tmp_path)
    with Image.open(tmp_path / "native_factorial_fas_sampling.png") as im:
        assert im.size == (2900, 1360)
        pixels = np.asarray(im.convert("RGB"))
        assert np.count_nonzero(np.any(pixels < 240, axis=2)) > 30000
    svg = (tmp_path / "native_factorial_fas_sampling.svg").read_text()
    for text in ("P0/C0/R0", "P1/C0/R1", "Observed native Z (descriptive)",
                 "not biological/generalization error bars", "All six expected-ratio differences"):
        assert text in svg
    assert (tmp_path / "native_factorial_fas_sampling.pdf").read_bytes().startswith(b"%PDF")


def test_occupied_output_rejected_before_reading_input(tmp_path):
    with pytest.raises(ValueError, match="Output already exists"):
        current.run(tmp_path / "missing_repo", tmp_path)
