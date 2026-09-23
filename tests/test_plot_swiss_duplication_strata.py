import csv
from pathlib import Path

import numpy as np
from PIL import Image
import pytest

from benchmark_tools.plot_swiss_duplication_strata import render
from benchmark_tools.prepare_ob_candidate_neighborhood import check


def test_render_all_endpoints_and_missing(tmp_path):
    root = Path(__file__).resolve().parents[1]
    output = tmp_path / "figure"
    report = render(root, output)
    assert report["endpoints"] == 48 and report["unavailable_endpoints"] == 6
    assert report["new_inferential_claims"] is False
    for item in report["outputs"]:
        check(item)
    with (output / "endpoints.tsv").open() as stream:
        rows = list(csv.DictReader(stream, delimiter="\t"))
    assert sum(r["difference"] == "NA" for r in rows) == 6
    assert {r["method"] for r in rows if r["difference"] == "NA"} == {"orthomcl_1_4"}
    pixels = np.asarray(Image.open(output / "swiss_duplication_descriptive.png").convert("RGB"))
    assert pixels.shape == (1620, 2700, 3) and pixels.std() > 5
    with pytest.raises(FileExistsError):
        render(root, output)


def test_changed_source_rejected(tmp_path):
    path = tmp_path / "benchmark_tools/results/swiss_duplication_strata_20260923/scores.tsv"
    path.parent.mkdir(parents=True)
    path.write_text("changed\n")
    with pytest.raises(ValueError):
        render(tmp_path, tmp_path / "figure")
    assert not (tmp_path / "figure").exists()
