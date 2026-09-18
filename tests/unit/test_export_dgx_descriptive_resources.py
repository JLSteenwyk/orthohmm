import pytest

from benchmark_tools.export_dgx_descriptive_resources import assemble


def fixture():
    native = {"status": "native_panel_validation_complete_not_timing_admission", "failures": 0, "runs": []}
    resource = {"failures": 0, "runs": []}
    host = {"review_failures": 0, "runs": []}
    plan = {"runs": []}
    for method in ("hmm", "phylo", "of"):
        for size in (4, 8, 12):
            for repeat in range(3):
                index = len(native["runs"])
                common = {"index": index, "method": method, "proteomes": size, "repeat": repeat}
                native["runs"].append({**common, "status": "native_outputs_validated"})
                host["runs"].append({**common, "host_review": {"job_id": index,
                    "retained_host_summary": {"status": "inconclusive"}}})
                plan["runs"].append({**common, "native_method": method,
                    "dataset": {"proteins": size * 100, "sequence_characters": size * 1000}})
                resource["runs"].append({"index": index, "job_id": index,
                    "status": "resource_accounting_reproduced_not_timing_admission",
                    "gnu_time": {"elapsed_seconds": [10, 30, 20][repeat], "user_seconds": 10,
                        "system_seconds": 2, "max_process_rss_kib": 100},
                    "summary": {"maximum_reported_cgroup_peak_bytes": 300000,
                        "maximum_sampled_sum_rss_bytes": 200000, "observations": 10,
                        "samples_with_process_errors": repeat}})
    return native, resource, host, plan


def test_all_runs_and_median_range_keep_uncertainty():
    report = assemble(*fixture())
    assert len(report["runs"]) == 27 and len(report["summaries"]) == 9
    assert report["summaries"][0]["native_elapsed_seconds"] == {"median": 20, "minimum": 10, "maximum": 30}
    assert all(r["host_inconclusive_runs"] == 3 and r["rss_incomplete_runs"] == 2 for r in report["summaries"])
    assert report["scientific_timings_admitted"] == 0


def test_failed_native_run_retained_without_survivor_only_average():
    args = fixture()
    args[0]["runs"][0].update(status="native_validation_failed", error="Missing gene")
    args[0]["failures"] = 1
    report = assemble(*args)
    assert len(report["runs"]) == 27
    assert report["runs"][0]["native_validation_error"] == "Missing gene"
    cell = report["summaries"][0]
    assert cell["native_valid_runs"] == 2
    assert cell["native_elapsed_seconds"] is None


@pytest.mark.parametrize("problem", ["missing", "resource_failure", "job", "method", "index", "repeat", "failure_count"])
def test_mismatched_or_incomplete_evidence_rejected(problem):
    n, r, h, p = fixture()
    if problem == "missing":
        r["runs"].pop()
    elif problem == "resource_failure":
        r["failures"] = 1
    elif problem == "job":
        r["runs"][0]["job_id"] = 999
    elif problem == "method":
        n["runs"][0]["method"] = "other"
    elif problem == "index":
        h["runs"][1]["index"] = 0
    elif problem == "repeat":
        for x in (n, h, p):
            x["runs"][1]["repeat"] = 0
    else:
        n["failures"] = 1
    with pytest.raises(ValueError):
        assemble(n, r, h, p)
