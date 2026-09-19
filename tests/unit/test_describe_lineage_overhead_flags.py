import copy

import pytest

from benchmark_tools.describe_lineage_overhead_flags import PERIODIC, describe, report


@pytest.fixture
def audit():
    scopes = ["/", "/system.slice", "/system.slice/spark-7ff0_slurmstepd.scope",
              "/system.slice/spark-7ff0_slurmstepd.scope/job_1"]
    interval = dict(narrow=dict(screen_passed=False, reasons=["excess_unassigned_cpu"],
                    signed_unassigned_average_cores=.3, native_cpu_s=10., wall_s=1.,
                    outer_read_overhang_s=.002),
                    lineage=dict(root_minus_target_cpu_usec=20,
                        signed_complements=[dict(ancestor=a, excluded_child=b, signed_cpu_usec=v)
                            for a, b, v in zip(scopes, scopes[1:], [30, -15, 5])]))
    runs = []
    for index in range(18):
        periodic = index in PERIODIC
        runs.append(dict(index=index, status="validated", method="fixture",
            mode="periodic" if periodic else "boundary", interval_screening_available=periodic,
            observation_points=2, narrow_flagged_intervals=[0] if periodic else None,
            lineage_screening=dict(intervals=[copy.deepcopy(interval)], narrow_flagged_intervals=[0])))
    return dict(validated_tasks=18, runs=runs)


def test_complete_description_keeps_signed_complements_and_boundary_unknown(audit):
    result = describe(audit)
    assert result["all"]["count"] == result["flagged"]["count"] == 9
    assert result["unflagged"]["count"] == 0
    assert len(result["boundary_tasks_without_interval_coverage"]) == 9
    assert result["flagged_intervals"][0]["system_minus_slurm_cpu_usec"] == -15
    assert result["scientific_timings_admitted"] is False
    assert result["environmental_validity_established"] is False


def test_unflagged_intervals_are_included(audit):
    for run in audit["runs"]:
        if run["index"] in PERIODIC:
            run["narrow_flagged_intervals"] = []
            run["lineage_screening"]["narrow_flagged_intervals"] = []
            run["lineage_screening"]["intervals"][0]["narrow"].update(screen_passed=True, reasons=[])
    result = describe(audit)
    assert result["all"]["count"] == result["unflagged"]["count"] == 9
    assert result["flagged"]["count"] == 0


@pytest.mark.parametrize("fault", ["missing", "failed", "boundary", "coverage", "flags", "path", "sum", "nan"])
def test_incomplete_or_inconsistent_evidence_rejected(audit, fault):
    run = audit["runs"][1]
    interval = run["lineage_screening"]["intervals"][0]
    if fault == "missing":
        audit["runs"].pop()
    elif fault == "failed":
        run["status"] = "failed"
    elif fault == "boundary":
        audit["runs"][0]["narrow_flagged_intervals"] = []
    elif fault == "coverage":
        run["observation_points"] = 3
    elif fault == "flags":
        run["narrow_flagged_intervals"] = []
    elif fault == "path":
        interval["lineage"]["signed_complements"][0]["excluded_child"] = "/other"
    elif fault == "sum":
        interval["lineage"]["root_minus_target_cpu_usec"] = 21
    else:
        interval["narrow"]["native_cpu_s"] = float("nan")
    with pytest.raises(ValueError):
        describe(audit)


def test_wrong_audit_hash_rejected(tmp_path):
    path = tmp_path / "audit.gz"
    path.write_bytes(b"not the frozen audit")
    with pytest.raises(ValueError, match="Wrong retained"):
        report(path)
