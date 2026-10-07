"""Fixed selection, descriptive scope and fixture rendering of stratum panels."""

from copy import deepcopy
import json
from pathlib import Path

import pytest

from benchmark_tools import plot_native_qfo_three_cell_strata as plot

ROOT = Path(__file__).resolve().parents[2]
RESULTS = ROOT / "benchmark_tools/results"


def inputs():
    return (json.loads((RESULTS / name).read_text()) for name, _ in plot.PINS.values())


def test_all_nonempty_nonredundant_bins_and_exact_values_are_selected():
    report, readback = inputs()
    points, excluded = plot.figure_data(report, readback)
    assert len(points) == 66 and len(excluded) == 9
    assert sum(r["reason"] == "empty_bin" for r in excluded) == 5
    assert sum(r["reason"] == "redundant_all_family_copy" for r in excluded) == 4
    assert {r["families"] for r in points} == {3, 6, 7, 9, 11, 12, 15, 18}
    for point in points:
        row = next(r for r in report["differences"] if all(
            point[k] == r[k] for k in ("suite", "stratum", "contrast")))
        assert point["difference_pp"] == 100 * row[point["metric"]]
    assert "not confidence intervals" in plot.CAPTION


@pytest.mark.parametrize("change", ["scope", "missing", "duplicate", "nan", "unit", "unplotted", "empty"])
def test_invalid_selection_and_scope_rejected(change):
    report, readback = inputs()
    report = deepcopy(report)
    if change == "scope":
        report["new_uncertainty"] = True
    elif change == "missing":
        report["differences"].pop()
    elif change == "duplicate":
        report["differences"][1] = deepcopy(report["differences"][0])
    elif change == "nan":
        next(r for r in report["differences"] if r["stratum"] == "higher_entropy")["PPV"] = float("nan")
    elif change == "unit":
        report["differences"][0]["reference"] = "wrong"
    elif change == "unplotted":
        report["bins"]["sequence"]["missing"] = ["APP"]
    else:
        readback["score_rows_checked"] = 0
    with pytest.raises((ValueError, KeyError)):
        plot.figure_data(report, readback)


def test_fixture_render_is_nonblank_and_has_all_formats(tmp_path):
    from PIL import Image
    import fitz
    report, readback = inputs()
    # Render invented fixture differences, not the selected numeric figure.
    for row in report["differences"]:
        for k in plot.METRICS:
            if row[k] is not None:
                row[k] = .01 if k == "TPR" else -.01
    points, _ = plot.figure_data(report, readback)
    plot.render(points, tmp_path)
    with Image.open(tmp_path / "fixed_stratum_contrasts.png") as image:
        assert image.size == (2240, 1280)
        assert len(image.convert("RGB").getcolors(image.width * image.height)) > 50
    with fitz.open(tmp_path / "fixed_stratum_contrasts.pdf") as doc:
        assert len(doc) == 1
        text = doc[0].get_text()
        assert "Candidate expansion minus baseline" in text
        assert "not confidence intervals" in text
    assert (tmp_path / "fixed_stratum_contrasts.svg").stat().st_size > 10000


def test_existing_namespace_refused_before_loading(tmp_path):
    with pytest.raises(ValueError, match="already exists"):
        plot.build(tmp_path / "missing_repo", tmp_path)
