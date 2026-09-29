import copy
import json
import os
from pathlib import Path
import subprocess
import sys

import pytest

from benchmark_tools import calibrate_threadripper_observer as workload
from benchmark_tools import audit_threadripper_observer as audit


@pytest.mark.parametrize("values", [(0, 1, 1, 1), (1, 0, 1, 1), (1, 1, 0, 1),
    (1, 1, 1, 0), (33, 1, 1, 1), (1, 5, 1, 1), (1, 1, 65, 1), (1, 1, 1, 31),
    (True, 1, 1, 1), (1, 1, 1, 1.0)])
def test_dimensions_bounded(values):
    with pytest.raises(ValueError):
        workload.validate(*values)


def test_maximum_dimensions_allowed():
    workload.validate(32, 4, 64, 30)


@pytest.mark.parametrize("index", [-1, 32, True])
def test_bad_worker_index_rejected(tmp_path, index):
    with pytest.raises(ValueError, match="index"):
        workload.worker(tmp_path, index, 1, 1, 1)


def test_real_small_workload_and_refused_retry(tmp_path):
    if not set(range(32)) <= os.sched_getaffinity(0):
        pytest.skip("This diagnostic requires CPUs 0-31")
    path = tmp_path / "witnesses"
    command = ["taskset", "-c", "0-31", sys.executable, "-B", str(Path(workload.__file__).resolve()),
               "--workload", "--output", str(path), "--workers", "2", "--threads", "2",
               "--arena-mib", "1", "--duration", "1"]
    result = subprocess.run(command, capture_output=True, text=True, timeout=25)
    assert result.returncode == 0, result.stderr
    witness = json.loads((path / "complete.json").read_text())
    assert witness["workers"] == 2
    for row in witness["witnesses"]:
        assert len(row["identity"]["threads"]) == 2
        assert row["identity"]["arena_bytes"] == 1024**2
        assert row["user_seconds"] + row["system_seconds"] > .05
        assert not Path(f'/proc/{row["identity"]["pid"]}').exists()
    result = subprocess.run(command, capture_output=True, text=True, timeout=10)
    assert result.returncode != 0
    assert "FileExistsError" in result.stderr


def test_failed_worker_cleans_up_owned_children(tmp_path, monkeypatch):
    children = []
    real_popen = subprocess.Popen
    def failing_popen(command, **kwargs):
        child = real_popen([sys.executable, "-c", "raise SystemExit(4)"] if not children
                           else [sys.executable, "-c", "import time; time.sleep(60)"], **kwargs)
        children.append(child)
        return child
    monkeypatch.setattr(workload.subprocess, "Popen", failing_popen)
    with pytest.raises(RuntimeError, match="readiness"):
        workload.workload(tmp_path / "failed", 2, 1, 1, 1)
    assert all(child.poll() is not None for child in children)
    assert not (tmp_path / "failed/complete.json").exists()


def fixture():
    scope = "/job_1/step_0/user/task_0"
    membership = f"0::{scope}\n"
    identity = dict(pid=10, index=0, start_ticks=123, affinity=list(range(32)),
        membership=membership, arena_bytes=1024**2,
        threads=[dict(tid=11, start_ticks=124, affinity=list(range(32)), membership=membership)])
    witness = dict(workers=1, threads_per_worker=1, arena_mib_per_worker=1, duration_seconds=1,
        scientific_timings_admitted=False, witnesses=[dict(identity=identity, started_ns=10**9,
        finished_ns=8 * 10**9, user_seconds=5., system_seconds=.1,
        final_membership=membership, final_affinity=list(range(32)))])
    scoped = dict(cpu=dict(scope=scope), native_outcome="exited_zero", native_exit_code=0,
                  primary=dict(cpu_seconds=6., peak_memory_bytes=2 * 1024**2, wall_seconds=9.))
    points = [dict(host=[dict(started_monotonic_ns=t * 10**9)], thread_affinity=dict(
        started_ns=t * 10**9, finished_ns=t * 10**9 + 100,
        status="observed_within_affinity", errors=[], violating_tids=[],
        threads=[dict(tid=10, start_ticks=123, cgroup=membership, affinity=list(range(32)), outside_cpus=[]),
                 dict(tid=11, start_ticks=124, cgroup=membership, affinity=list(range(32)), outside_cpus=[])]))
        for t in range(2, 7)]
    done = dict(started_ns=1, finished_ns=9 * 10**9)
    return witness, scoped, points, done


def evaluate(values):
    return audit.evaluate(*values, dimensions=dict(workers=1, threads_per_worker=1,
                          arena_mib_per_worker=1, duration_seconds=1))


def test_complete_calibration_never_admits_production():
    result = evaluate(fixture())
    assert result["status"] == "calibration_checks_passed"
    assert result["expected_threads"] == 2
    assert result["complete_interior_samples"] == 5
    for key in ("controlled_timing_admitted", "full_run_containment_verified", "slowdown_overhead_validated",
                "publication_ready"):
        assert result[key] is False


@pytest.mark.parametrize("defect", ["missing_thread", "wrong_ticks", "wrong_membership", "wide_affinity",
    "racing_sample", "duplicate_thread", "point_cost", "sample_gap", "few_samples", "low_cpu",
    "high_cpu", "low_peak", "failed_native"])
def test_failed_checks_preserved(defect):
    values = copy.deepcopy(fixture())
    witness, scoped, points, done = values
    sample = points[0]["thread_affinity"]
    if defect == "missing_thread":
        sample["threads"].pop()
    elif defect == "wrong_ticks":
        sample["threads"][1]["start_ticks"] += 1
    elif defect == "wrong_membership":
        sample["threads"][1]["cgroup"] = "0::/wrong\n"
    elif defect == "wide_affinity":
        sample["threads"][1]["affinity"].append(32)
    elif defect == "racing_sample":
        sample["errors"] = [dict(error="raced")]
    elif defect == "duplicate_thread":
        sample["threads"].append(copy.deepcopy(sample["threads"][0]))
    elif defect == "point_cost":
        sample["finished_ns"] += 2 * 10**9
    elif defect == "sample_gap":
        points.pop(1)
    elif defect == "few_samples":
        points.pop()
    elif defect == "low_cpu":
        scoped["primary"]["cpu_seconds"] = 1.
    elif defect == "high_cpu":
        scoped["primary"]["cpu_seconds"] = 40.
    elif defect == "low_peak":
        scoped["primary"]["peak_memory_bytes"] = 1
    else:
        scoped.update(native_outcome="exited_nonzero", native_exit_code=1)
    assert evaluate(values)["status"] == "calibration_checks_failed"


@pytest.mark.parametrize("defect", ["worker_scope", "final_scope", "allocation", "threads", "cpu_nan",
                                        "cpu_negative", "interval", "duplicate_witness"])
def test_invalid_witness_rejected(defect):
    values = copy.deepcopy(fixture())
    row = values[0]["witnesses"][0]
    if defect == "worker_scope":
        row["identity"]["membership"] = "wrong"
    elif defect == "final_scope":
        row["final_membership"] = "wrong"
    elif defect == "allocation":
        row["identity"]["arena_bytes"] = 1
    elif defect == "threads":
        row["identity"]["threads"].clear()
    elif defect == "cpu_nan":
        row["user_seconds"] = float("nan")
    elif defect == "cpu_negative":
        row["system_seconds"] = -1
    elif defect == "interval":
        row["finished_ns"] = row["started_ns"]
    else:
        row["identity"]["threads"][0]["tid"] = row["identity"]["pid"]
    with pytest.raises(ValueError):
        evaluate(values)


def test_existing_audit_refused(tmp_path):
    with pytest.raises(FileExistsError):
        audit.audit(tmp_path, tmp_path / "missing", "unused", tmp_path)


def test_protocol_digest_checked_first(tmp_path):
    path = tmp_path / "protocol.json"
    path.write_text("not JSON")
    with pytest.raises(ValueError, match="digest"):
        audit.verify_protocol(path, "wrong")


def test_protocol_checks_sources_and_prespecified_envelope(tmp_path):
    endpoint = tmp_path / "endpoint.json"
    endpoint.write_text("{}")
    value = dict(schema="threadripper_observer_calibration_v1", dimensions=audit.DIMENSIONS,
                 limits=audit.LIMITS, production_timing=False,
                 resource_protocol=audit.record(endpoint), interpreter=audit.record(sys.executable),
                 sources=[audit.record(audit.__file__), audit.record(workload.__file__)])
    path = tmp_path / "protocol.json"
    path.write_text(json.dumps(value))
    assert audit.verify_protocol(path, audit.record(path)["sha256"]) == value
    value["limits"] = {**audit.LIMITS, "min_complete_samples": 1}
    path.write_text(json.dumps(value))
    with pytest.raises(ValueError, match="differs"):
        audit.verify_protocol(path, audit.record(path)["sha256"])
    value["limits"] = audit.LIMITS
    value["sources"].pop()
    path.write_text(json.dumps(value))
    with pytest.raises(ValueError, match="omitted"):
        audit.verify_protocol(path, audit.record(path)["sha256"])
