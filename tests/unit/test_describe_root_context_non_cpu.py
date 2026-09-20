from copy import deepcopy
import json
from pathlib import Path
import shutil

import pytest

from benchmark_tools.describe_root_context_non_cpu import describe, run


def fixture():
    row = dict(status="validated", method="example", arm="lineage", scientific_timings_admitted=False,
        original_flagged_intervals=[1], narrow_flagged_intervals=[1],
        pressure_observation_window=dict(status="native_step_pressure_diagnostic",
            scientific_timings_admitted=False, controlled_workload_verified=False,
            native_stall_usec={r: dict(some=0, full=0) for r in ("cpu", "memory", "io")}),
        memory=dict(errors=[], raw={"memory.current": "10\n", "memory.peak": "20\n",
            "memory.events": "low 0\nhigh 0\nmax 0\noom 0\noom_kill 0\n"}))
    return dict(validated_tasks=18, issues=[], scientific_timings_admitted=False,
                runs=[dict(deepcopy(row), index=i) for i in range(18)])


def test_all_outcomes_and_positive_events_preserved_without_admission():
    audit = fixture()
    audit["runs"][2]["pressure_observation_window"]["native_stall_usec"]["memory"]["some"] = 100
    audit["runs"][4]["pressure_observation_window"]["native_stall_usec"]["io"]["full"] = 200
    audit["runs"][7]["memory"]["raw"]["memory.events"] += "oom_group_kill 1\n"
    result = describe(audit)
    assert len(result["runs"]) == 18
    assert result["observed_memory_stall_tasks"] == [2]
    assert result["observed_io_stall_tasks"] == [4]
    assert result["nonzero_memory_event_tasks"] == [7]
    assert result["scientific_timings_admitted"] is False
    assert result["environmental_validity_established"] is False
    assert all(r["original_flags"] == r["narrow_flags"] == 1 for r in result["runs"])


@pytest.mark.parametrize("case", ["missing", "order", "bool_index", "invalid", "admission", "issues",
    "negative", "bool_counter", "resource", "memory_error", "events_missing", "events_duplicate", "peak"])
def test_invalid_observations_never_become_zero_or_quiet(case):
    audit = fixture()
    row = audit["runs"][0]
    if case == "missing": audit["runs"].pop()
    elif case == "order": row["index"] = 3
    elif case == "bool_index": row["index"] = False
    elif case == "invalid": row["status"] = "failed_or_invalid"
    elif case == "admission": audit["scientific_timings_admitted"] = True
    elif case == "issues": audit["issues"] = ["missing"]
    elif case == "negative": row["pressure_observation_window"]["native_stall_usec"]["io"]["some"] = -1
    elif case == "bool_counter": row["pressure_observation_window"]["native_stall_usec"]["io"]["some"] = False
    elif case == "resource": del row["pressure_observation_window"]["native_stall_usec"]["memory"]
    elif case == "memory_error": row["memory"]["errors"] = ["missing"]
    elif case == "events_missing": row["memory"]["raw"]["memory.events"] = "oom 0\n"
    elif case == "events_duplicate": row["memory"]["raw"]["memory.events"] += "oom 0\n"
    elif case == "peak": row["memory"]["raw"]["memory.peak"] = "9\n"
    with pytest.raises(ValueError):
        describe(audit)


def test_real_audit_relocation_reproduces_retained_description(tmp_path):
    results = Path(__file__).resolve().parents[2] / "benchmark_tools/results"
    path = tmp_path / "audit.json.gz"
    shutil.copyfile(results / "root_context_overhead_audit_22022_20260920.json.gz", path)
    observed = run(path)
    retained = json.loads((results / "root_context_non_cpu_22022_20260920.json").read_text())
    assert {k:v for k,v in observed.items() if k != "sources"} == {
        k:v for k,v in retained.items() if k != "sources"}
    assert observed["observed_memory_stall_tasks"] == [2, 3, 4, 5, 6, 7, 8, 9, 12, 13, 16, 17]
    assert observed["observed_io_stall_tasks"] == list(range(18))
    assert observed["nonzero_memory_event_tasks"] == []


def test_changed_audit_rejected_before_parsing(tmp_path):
    path = tmp_path / "audit.json.gz"
    path.write_bytes(b"not the pinned audit")
    with pytest.raises(ValueError, match="hash"):
        run(path)
