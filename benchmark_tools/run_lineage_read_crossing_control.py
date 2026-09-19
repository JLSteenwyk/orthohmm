"""Finite owned-service read crossing; never a scientific timing experiment."""

import argparse
import hashlib
import os
from pathlib import Path, PurePosixPath
import shutil
import subprocess
import sys
import time

from benchmark_tools.probe_cgroup_lineage import snapshot, compare, LineageSnapshotError
from benchmark_tools.probe_dgx_step_separation import save
from benchmark_tools.run_lineage_lifecycle_control import membership, trial as lifecycle_trial


def assess(before, crossing, after, event, nested):
    spans = dict(before_to_crossing=compare(before, crossing),
                 crossing_to_after=compare(crossing, after), full=compare(before, after))
    rows = crossing["rows"]
    if (any(type(event[k]) is not int for k in ("started_ns", "finished_ns"))
            or rows[0]["scope"] != "/" or len(rows) < 2
            or not rows[0]["finished_ns"] <= event["started_ns"] < event["finished_ns"] <= rows[1]["started_ns"]):
        raise ValueError("Lifecycle event did not cross the intended counter reads")
    if (nested["status"] != "lifecycle_control_evaluated"
            or nested["result"]["descendant_cpu_retention_response"] is not True
            or nested["result"]["outside_target_response"] is not True):
        raise ValueError("Nested lifecycle response failed")
    response = spans["full"]["root_minus_target_cpu_usec"] >= 500000
    return dict(status="read_crossing_control_evaluated", event=event, spans=spans,
                full_span_response=response, control_met=response,
                scientific_timings_admitted=False, environmental_validity_established=False)


def trial(directory, index, job, manager, target, cpu):
    directory.mkdir()
    root = Path("/sys/fs/cgroup")
    before = snapshot(root, target)
    save(directory / "before.json", before)
    event, nested = {}, None

    def inject(scope):
        nonlocal nested
        if scope != "/":
            return
        event["started_ns"] = time.monotonic_ns()
        try:
            nested = lifecycle_trial(directory / "lifecycle", index, job, manager, target, cpu)
            save(directory / "lifecycle_result.json", nested)
        except Exception as error:
            raise ValueError("Owned lifecycle failed: " + str(error)) from error
        finally:
            event["finished_ns"] = time.monotonic_ns()
            save(directory / "event.json", event)

    crossing = snapshot(root, target, after_read=inject)
    save(directory / "crossing.json", crossing)
    after = snapshot(root, target)
    save(directory / "after.json", after)
    if membership() != target:
        raise ValueError("Observer membership changed")
    return assess(before, crossing, after, event, nested)


def run(output):
    if (os.uname().nodename != "spark-7ff0" or os.environ.get("SLURM_CPUS_PER_TASK") != "2"
            or len(os.sched_getaffinity(0)) != 2):
        raise ValueError("Require exclusive DGX allocation with two-CPU observer step")
    job = int(os.environ["SLURM_JOB_ID"])
    target = membership()
    if f"job_{job}" not in PurePosixPath(target).parts:
        raise ValueError("Observer outside assigned job")
    manager = subprocess.check_output(["systemctl", "--user", "show", "--property=ControlGroup",
                                       "--value"], text=True, timeout=10).strip()
    output.mkdir(exist_ok=False)
    sources = list(Path(__file__).parent.glob("*.py"))
    sources += [Path(__file__).parent / "results/LINEAGE_READ_CROSSING_PROTOCOL_20260919.md"]
    hashes = {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in sources}
    save(output / "identity.json", dict(job=job, target=target, manager=manager, sources=hashes,
        affinity=sorted(os.sched_getaffinity(0)), uname=list(os.uname()), python=sys.version,
        executables={name: dict(path=str(Path(shutil.which(name)).resolve()),
            sha256=hashlib.sha256(Path(shutil.which(name)).read_bytes()).hexdigest())
            for name in (sys.executable, "systemd-run", "systemctl", "taskset")}))
    rows = []
    for index in range(3):
        directory = output / f"trial_{index}"
        try:
            result = trial(directory, index, job, manager, target, min(os.sched_getaffinity(0)))
        except Exception as error:
            result = dict(status="failed", error=str(error), error_type=type(error).__name__)
            if isinstance(error, LineageSnapshotError):
                result["partial_snapshot"] = error.evidence
        save(directory / "result.json", result)
        rows.append(result)
    unchanged = hashes == {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in sources}
    result = dict(status="read_crossing_controls_finished", trials=rows, sources_unchanged=unchanged,
        all_controls_met=unchanged and all(r.get("control_met") is True for r in rows),
        scientific_timings_admitted=False, environmental_validity_established=False,
        limitations=["Deliberate read-window crossing, not production monitor overhead.",
                     "No process attribution, timing correction or scientific inclusion."])
    save(output / "report.json", result)
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = run(args.output.absolute())
    raise SystemExit(0 if result["all_controls_met"] else 1)
