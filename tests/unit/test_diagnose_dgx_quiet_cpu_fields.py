import copy
import json
from pathlib import Path

import pytest

from benchmark_tools.diagnose_dgx_quiet_cpu_fields import diagnose, run

RESULTS = Path(__file__).resolve().parents[2] / "benchmark_tools/results"
SOURCE = RESULTS / "dgx_hierarchy_quiet_smokes_21820.json"


def measurement(index=1):
    return json.loads(SOURCE.read_text())["runs"][index]["verification"]["measurement"]


def test_retained_flag_and_field_identity():
    result = diagnose(measurement())
    assert result["screening"]["original_threshold_screen"]["flagged_intervals"] == [3]
    flagged = result["intervals"][3]
    assert flagged["diagnostic_field_residual_cpu_s"] == pytest.approx(dict(
        user_and_nice=.251996, system=.095421, irq_and_softirq=.02))
    assert flagged["host_outer_minus_job_outer_cpu_s"] == pytest.approx(.367417)
    assert flagged["read_window_busy_cpu_s"] == pytest.approx(.01)
    for row in result["intervals"]:
        assert sum(row["diagnostic_field_residual_cpu_s"].values()) == pytest.approx(
            row["host_outer_minus_job_outer_cpu_s"] + row["job_usage_minus_user_system_cpu_s"])


def test_reject_changed_screen():
    changed = copy.deepcopy(measurement())
    changed["screening"]["original_threshold_screen"]["flagged_intervals"] = []
    with pytest.raises(ValueError, match="counter replay"):
        diagnose(changed)


def test_all_methods_retained_and_result_reproducible(tmp_path):
    result = run(SOURCE, tmp_path / "result.json")
    retained = json.loads((RESULTS / "dgx_quiet_cpu_fields_20260918.json").read_text())
    assert result == retained
    assert len(result["runs"]) == 3
    assert result["scientific_timings_admitted"] is False
    assert any(row["host_outer_minus_job_outer_cpu_s"] < 0
               for method in result["runs"] for row in method["intervals"])
    with pytest.raises(FileExistsError):
        run(SOURCE, tmp_path / "result.json")


def test_reject_changed_source(tmp_path):
    path = tmp_path / "source.json"
    path.write_text(SOURCE.read_text() + "\n")
    with pytest.raises(ValueError, match="Changed retained"):
        run(path, tmp_path / "output.json")
    assert not (tmp_path / "output.json").exists()
