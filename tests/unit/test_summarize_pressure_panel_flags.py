import pytest

import gzip
import json

from benchmark_tools.prepare_ob_candidate_neighborhood import record
from benchmark_tools.summarize_pressure_panel_flags import distribution, summarize, unique_report, report


@pytest.mark.parametrize("values", [[], [float("nan")], [True], [float("inf")]])
def test_invalid_distribution(values):
    with pytest.raises(ValueError):
        distribution(values)


def test_signed_distribution_preserved():
    assert distribution([-3, 1, 2]) == dict(minimum=-3, median=1, maximum=2)


def screening():
    return dict(original_threshold_screen=dict(
        intervals=[dict(screen_passed=False, reasons=["negative"], signed_unassigned_average_cores=-1,
                        outer_read_overhang_s=.02)],
        flagged_intervals=[0], whole_command_screen={}),
        frontier_intervals=[dict(outside_target_frontier_cpu_s=.01, root_minus_frontier_cpu_s=-.2)],
        native_pressure_whole_command={})


def test_reasons_and_signs_retained():
    result = summarize(screening())
    assert result["reasons"] == {"negative": 1}
    assert result["root_minus_frontier_cpu_s"]["median"] == -.2
    assert result["outer_read_overhang_s"]["median"] == .02


def test_changed_flags_rejected():
    data = screening()
    data["original_threshold_screen"]["flagged_intervals"] = []
    with pytest.raises(ValueError, match="flags differ"):
        summarize(data)


def test_mismatched_intervals_rejected():
    data = screening()
    data["frontier_intervals"] = []
    with pytest.raises(ValueError, match="inventories differ"):
        summarize(data)


def test_duplicate_identical_report_allowed():
    item = dict(path="/a/frontier_report.json", bytes=10, sha256="a")
    assert unique_report([item, dict(item)]) == item


def test_conflicting_duplicate_report_rejected():
    item = dict(path="/a/frontier_report.json", bytes=10, sha256="a")
    with pytest.raises(ValueError, match="Conflicting"):
        unique_report([item, dict(item, sha256="b")])


@pytest.mark.parametrize("paths", [[], ["/a/frontier_report.json", "/b/frontier_report.json"]])
def test_missing_or_multiple_reports_rejected(paths):
    with pytest.raises(ValueError, match="unique"):
        unique_report([dict(path=p) for p in paths])


def test_report_end_to_end_and_tamper(tmp_path):
    path = tmp_path / "frontier_report.json"
    path.write_text(json.dumps(dict(screening=screening())))
    evidence = record(path)
    runs = [dict(index=i, method="test", mode="boundary", status="failed") for i in range(18)]
    runs[1].update(mode="periodic", status="validated", evidence=[evidence, dict(evidence)])
    audit = tmp_path / "audit.json.gz"
    audit.write_bytes(gzip.compress(json.dumps(dict(runs=runs)).encode()))
    checksum = record(audit)["sha256"]
    result = report(audit, checksum)
    assert len(result["runs"]) == 18
    assert result["runs"][0]["status"] == "failed"
    assert result["runs"][1]["diagnostic"]["flagged_intervals"] == 1
    assert result["scientific_timings_admitted"] is False
    with pytest.raises(ValueError, match="checksum"):
        report(audit, "incorrect")
    path.write_text("changed")
    with pytest.raises(ValueError):
        report(audit, checksum)
