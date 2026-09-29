import json
from pathlib import Path

import pytest

from benchmark_tools.plot_ob_complete_uncertainty import endpoints, run
from benchmark_tools.prepare_ob_candidate_neighborhood import record

RESULT = Path(__file__).resolve().parents[2] / "benchmark_tools/results/ob_complete_uncertainty_20260928.json"


def test_all_endpoints_and_existing_units():
    result = json.loads(RESULT.read_text())
    rows = endpoints(result)
    assert len(rows) == 21
    row = next(r for r in rows if r["method"] == "orthohmm_phylogeny_satellite_v2" and r["metric"] == "f_score")
    assert row["difference_pp"] == pytest.approx(1.369593514512374)
    assert row["adjusted_low_pp"] < 0 < row["adjusted_high_pp"]


@pytest.mark.parametrize("change", ["missing", "baseline", "multiplicity", "nan", "shape", "reversed", "difference"])
def test_invalid_input_rejected(change):
    result = json.loads(RESULT.read_text())
    values = result["comparisons"]["orthohmm_high_sensitivity"]["metrics"]["f_score"]
    if change == "missing": result["comparisons"].pop("fastoma_0_3_5")
    elif change == "baseline": result["baseline"] = "other"
    elif change == "multiplicity": result["multiplicity"] = "six endpoints"
    elif change == "nan": values["difference_percentage_points"] = float("nan")
    elif change == "shape": values["paired_percentile_ci"] = [0]
    elif change == "reversed": values["bonferroni_percentile_ci"].reverse()
    else: values["difference_percentage_points"] += 1
    with pytest.raises(ValueError):
        endpoints(result)


def test_render_and_refuse_overwrite(tmp_path):
    output = tmp_path / "figure"
    manifest = run(RESULT, record(RESULT)["sha256"], output)
    assert manifest["endpoints"] == 21
    assert len((output / "endpoints.tsv").read_text().splitlines()) == 22
    assert len(manifest["outputs"]) == 4
    assert manifest["new_statistics"] is False
    with pytest.raises(FileExistsError):
        run(RESULT, record(RESULT)["sha256"], output)


def test_checksum_rejected_before_output(tmp_path):
    output = tmp_path / "figure"
    with pytest.raises(ValueError, match="checksum"):
        run(RESULT, "0" * 64, output)
    assert not output.exists()
