"""Bounded sleeping-thread controls; deliberately widen only our own child thread."""

import argparse
import json
import os
from pathlib import Path
import subprocess
import sys
import threading
import time

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.measure_threadripper_scaling import measure, TIMEOUT
from benchmark_tools.replay_threadripper_scaling import replay
from benchmark_tools.probe_dgx_step_separation import save
from benchmark_tools.prepare_ob_candidate_neighborhood import record


def workload(mode, child=False):
    if sorted(os.sched_getaffinity(0)) != list(range(32)):
        raise ValueError("Control requires inherited CPUs 0-31")
    ready = threading.Barrier(3, timeout=10)
    release = threading.Event()
    rows, errors = [], []

    def target(index):
        try:
            if index == 0:
                os.sched_setaffinity(0, {0, 96} if child and mode == "widen" else {0})
            rows.append(dict(tid=threading.get_native_id(), affinity=sorted(os.sched_getaffinity(0))))
            ready.wait()
            release.wait(12)
        except Exception as error:
            errors.append(str(error))
            ready.abort()

    threads = [threading.Thread(target=target, args=(i,)) for i in range(2)]
    process = None
    try:
        for thread in threads:
            thread.start()
        ready.wait()
        print(json.dumps(dict(pid=os.getpid(), child=child, mode=mode,
                              leader_affinity=sorted(os.sched_getaffinity(0)), threads=rows)), flush=True)
        if not child:
            process = subprocess.Popen([sys.executable, "-B", str(Path(__file__).resolve()),
                                        "--workload", mode, "--child"])
        time.sleep(8)
        if process is not None and process.wait(timeout=15) != 0:
            raise ValueError("Child control failed")
        if errors:
            raise ValueError(str(errors))
    finally:
        release.set()
        for thread in threads:
            if thread.ident is not None:
                thread.join(timeout=15)
        if process is not None and process.poll() is None:
            # Only this control's own child, never an unrelated workload.
            process.terminate()
            process.wait(timeout=15)


def validate_control(mode, identities, points):
    if mode not in ("clean", "widen") or len(identities) != 2:
        raise ValueError("Require two process identities and a known mode")
    if sorted(row["child"] for row in identities) != [False, True]:
        raise ValueError("Missing parent or child identity")
    expected, offenders = {}, []
    for row in identities:
        if row["mode"] != mode or row["leader_affinity"] != list(range(32)) or len(row["threads"]) != 2:
            raise ValueError("Control identity differs")
        expected[row["pid"]] = list(range(32))
        masks = sorted(t["affinity"] for t in row["threads"])
        narrow = [0, 96] if mode == "widen" and row["child"] else [0]
        if masks != sorted([narrow, list(range(32))]):
            raise ValueError("Control thread masks differ")
        for thread in row["threads"]:
            if thread["tid"] in expected:
                raise ValueError("Duplicate identity")
            expected[thread["tid"]] = thread["affinity"]
            if 96 in thread["affinity"]:
                offenders.append(thread["tid"])
    matches = []
    for i, point in enumerate(points):
        value = point["thread_affinity"]
        observed = {row["tid"]: row["affinity"] for row in value["threads"]}
        if all(observed.get(tid) == mask for tid, mask in expected.items()) and not value["errors"]:
            status = "violation" if offenders else "observed_within_affinity"
            if value["status"] != status or sorted(value["violating_tids"]) != sorted(offenders):
                raise ValueError("Monitor did not detect the expected thread-level outcome")
            matches.append(i)
    if not matches:
        raise ValueError("No complete simultaneous observation of all control threads")
    return dict(mode=mode, expected_masks=expected, expected_violating_tids=offenders,
                complete_matching_points=matches, scientific_timings_admitted=False)


def run(output):
    output.mkdir(exist_ok=False)
    job = int(os.environ["SLURM_JOB_ID"])
    sources = [record(Path(__file__).resolve())]
    controls = []
    for mode in ("clean", "widen"):
        directory = output / mode
        command = [sys.executable, "-B", str(Path(__file__).resolve()), "--workload", mode]
        measured = measure(command, directory, job, 32, 128*1024**3, TIMEOUT, 1.)
        audited = replay(directory, job, command)
        save(output / f"{mode}_replay.json", audited)
        if audited["native_outcome"] != "exited_zero":
            raise ValueError("Native control did not exit zero")
        identities = [json.loads(line) for line in (directory / "native.log").read_text().splitlines()]
        controls.append(validate_control(mode, identities, measured["points"]))
    save(output / "controls.json", dict(job_id=job, sources=sources, controls=controls,
        status="thread_affinity_controls_observed", scientific_timings_admitted=False,
        limitations=["Sleeping threads, not a compute-heavy native engine or overhead assessment.",
                     "Deliberate widening changes only the control child thread, within its reserved cpuset."]))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--workload", choices=("clean", "widen"))
    group.add_argument("--output", type=Path)
    parser.add_argument("--child", action="store_true")
    args = parser.parse_args()
    if args.workload:
        workload(args.workload, args.child)
    else:
        run(args.output.resolve())
