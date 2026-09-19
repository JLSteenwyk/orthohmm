import csv
import json

from PIL import Image
import numpy as np
import pytest

from benchmark_tools.plot_corrected_swiss_comparison import endpoints, run
from benchmark_tools.prepare_ob_candidate_neighborhood import record
from tests.unit.test_reproduce_corrected_swiss_comparison import result


def test_tsv_retains_all_planned_endpoints_and_missing_rows(result):
    rows = endpoints(result)
    assert len(rows) == 24
    assert sum(row["difference"] is None for row in rows) == 6
    assert rows[0]["difference"] == result["comparisons"][0]["metrics"]["F1"]["difference"]


def test_rendered_assets_and_raw_scale_table(tmp_path, result):
    path = tmp_path / "input.json"
    path.write_text(json.dumps(result))
    output = tmp_path / "figure"
    manifest = run(path, record(path)["sha256"], output)
    assert len(manifest["outputs"]) == 5
    assert manifest["endpoints"] == 24
    assert manifest["display_scale"] == 100 and manifest["tsv_scale"] == 1
    with (output / "endpoints.tsv").open() as stream:
        rows = list(csv.DictReader(stream, delimiter="\t"))
    assert len(rows) == 24
    assert float(rows[0]["difference"]) == result["comparisons"][0]["metrics"]["F1"]["difference"]
    assert sum(row["difference"] == "" for row in rows) == 6
    pixels = np.asarray(Image.open(output / "corrected_swiss_comparison.png").convert("RGB"))
    assert pixels.shape == (1260, 2700, 3)
    assert np.mean(np.min(pixels, axis=2) < 200) > .025
    with pytest.raises(FileExistsError):
        run(path, record(path)["sha256"], output)
