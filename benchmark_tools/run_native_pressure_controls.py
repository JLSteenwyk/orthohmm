"""Frozen same-core CPU-pressure controls, not native inference timing."""

import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import time

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.probe_native_pressure import read_point, compare
from benchmark_tools.probe_dgx_step_separation import burn, save, wait_file

ORDER = (("quiet", "native-only", "contended"), ("native-only", "contended", "quiet"),
         ("contended", "quiet", "native-only"))
PROTOCOL_SHA = "812bedfc6b9127053a1de426c17f36c864005164d6928da6f396db43c46109a7"


def work(cpu, quiet=False):
    os.sched_setaffinity(0, {cpu})
    start = time.monotonic_ns()
    if quiet:
        cpu_start = time.process_time()
        time.sleep(1.5)
        result = dict(cpu_seconds=time.process_time()-cpu_start,
                      cgroup=Path("/proc/self/cgroup").read_text(), pid=os.getpid())
    else:
        result = burn()
    result.update(started_ns=start, finished_ns=time.monotonic_ns(), affinity=sorted(os.sched_getaffinity(0)))
    return result


def worker(directory):
    allowed = sorted(os.sched_getaffinity(0))
    cpu = allowed[0]
    os.sched_setaffinity(0, {cpu})
    save(directory / "ready.json", dict(pid=os.getpid(), cpu=cpu, allowed=allowed,
         membership=Path("/proc/self/cgroup").read_text()))
    go = wait_file(directory / "go.json")
    save(directory / "done.json", work(cpu, quiet=go["mode"] == "quiet"))
    wait_file(directory / "release.json")


def validate_trial(row, job):
    pressure = compare(*row["points"], job)
    if pressure != row["pressure"] or row["worker_exit_code"] != 0:
        raise ValueError("Pressure replay or worker completion differs")
    cpu = row["ready"]["cpu"]
    if cpu not in row["ready"]["allowed"] or cpu not in row["observer_allowed"]:
        raise ValueError("Worker CPU outside recorded allocation affinity")
    if row["observer_cpu"] == cpu or row["observer_cpu"] not in row["observer_allowed"]:
        raise ValueError("Observer must use another permitted CPU")
    native = row["native"]
    if (native["cgroup"] != row["points"][0]["native_membership"] or native["affinity"] != [cpu]
            or native["pid"] != row["ready"]["pid"]
            or row["ready"]["membership"] != native["cgroup"]):
        raise ValueError("Native work has wrong membership or affinity")
    if not (row["points"][0]["host"][1]["finished_monotonic_ns"] < native["started_ns"]
            < native["finished_ns"] < row["points"][1]["host"][0]["started_monotonic_ns"]):
        raise ValueError("Pressure observations do not enclose native work")
    mode = row["mode"]
    if mode not in ("quiet", "native-only", "contended"):
        raise ValueError("Unknown control mode")
    if mode != "quiet" and not .75 <= native["cpu_seconds"] <= .90:
        raise ValueError("Native CPU dose outside frozen range")
    if mode == "quiet" and native["finished_ns"] - native["started_ns"] < 1_500_000_000:
        raise ValueError("Quiet observation shorter than frozen duration")
    competitor = row["competitor"]
    overlap = None
    if mode == "contended":
        if (competitor is None or competitor["affinity"] != [cpu]
                or competitor["cgroup"] != row["points"][0]["host"][0]["raw"]["cgroup_membership"]
                or not .75 <= competitor["cpu_seconds"] <= .90):
            raise ValueError("Wrong competitor scope, affinity or CPU dose")
        overlap = (min(native["finished_ns"], competitor["finished_ns"])
                   - max(native["started_ns"], competitor["started_ns"])) / 1e9
        if overlap < .25:
            raise ValueError("Insufficient recorded work overlap")
    elif competitor is not None:
        raise ValueError("Unexpected competitor in control")
    return dict(status="control_injection_validated", overlap_s=overlap)


def summarize(rows):
    expected = [(block, mode) for block, modes in enumerate(ORDER) for mode in modes]
    if [(r["block"], r["mode"]) for r in rows] != expected:
        raise ValueError("Control inventory/order differs from frozen design")
    blocks = []
    for block in range(3):
        subset = {r["mode"]: r for r in rows if r["block"] == block}
        valid = all(r["status"] == "validated" for r in subset.values())
        delta = (subset["contended"]["pressure"]["native_stall_usec"]["cpu"]["some"]
                 - subset["native-only"]["pressure"]["native_stall_usec"]["cpu"]["some"]) if valid else None
        blocks.append(dict(block=block, injection_valid=valid, difference_usec=delta,
                           response_check_passed=(delta >= 100000 if valid else None)))
    return dict(blocks=blocks, all_response_checks_passed=all(b["response_check_passed"] is True for b in blocks),
                scientific_timings_admitted=False, controlled_workload_verified=False)


def trial(directory, block, mode, job):
    directory.mkdir()
    original = set(os.sched_getaffinity(0))
    source = Path(__file__).resolve()
    argv = ["srun", "--exclusive", "--exact", "--nodes=1", "--ntasks=1", "--cpus-per-task=1",
            sys.executable, "-B", str(source), "--worker", str(directory)]
    competitor = None
    with (directory / "step.log").open("x") as log:
        process = subprocess.Popen(argv, stdout=log, stderr=subprocess.STDOUT)
        try:
            ready = wait_file(directory / "ready.json")
            candidates = sorted(original - {ready["cpu"]})
            if not candidates or ready["cpu"] not in original:
                raise ValueError("Cannot separate observer and native affinities")
            observer_cpu = candidates[0]
            os.sched_setaffinity(0, {observer_cpu})
            before = read_point(ready["pid"], ready["membership"], job)
            save(directory / "before.json", before)
            save(directory / "go.json", dict(mode=mode))
            if mode == "contended":
                competitor = subprocess.Popen([sys.executable, "-B", str(source), "--burn-on", str(ready["cpu"])],
                                              stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
            native = wait_file(directory / "done.json")
            load = None
            if competitor is not None:
                stdout, stderr = competitor.communicate(timeout=15)
                save(directory / "competitor.json", dict(stdout=stdout, stderr=stderr, exit_code=competitor.returncode))
                if competitor.returncode:
                    raise ValueError("Competitor failed")
                load = json.loads(stdout)
            time.sleep(.2)
            after = read_point(ready["pid"], ready["membership"], job)
            save(directory / "after.json", after)
        finally:
            if competitor is not None and competitor.poll() is None:
                competitor.kill()
                competitor.wait()
            for name, value in (("go.json", dict(mode="quiet")), ("release.json", dict(release=True))):
                if not (directory / name).exists():
                    save(directory / name, value)
            try:
                code = process.wait(timeout=60)
            finally:
                os.sched_setaffinity(0, original)
    row = dict(block=block, mode=mode, ready=ready, observer_allowed=sorted(original), observer_cpu=observer_cpu,
               native=native, competitor=load, points=[before, after], worker_exit_code=code, argv=argv,
               pressure=compare(before, after, job))
    row["validation"] = validate_trial(row, job)
    row["status"] = "validated"
    return row


def run(output):
    if (os.uname().nodename != "spark-7ff0" or os.environ.get("SLURM_CPUS_PER_TASK") != "2"
            or os.environ.get("SLURM_MEM_PER_NODE") != "256"):
        raise ValueError("Require DGX two-CPU/256MiB control allocation")
    base = Path(__file__).resolve().parent
    paths = [base / name for name in ("run_native_pressure_controls.py", "probe_native_pressure.py",
             "audit_dgx_pressure.py", "probe_host_counters.py", "probe_dgx_step_separation.py",
             "summarize_frontier_overhead.py", "results/NATIVE_PRESSURE_CONTROL_PROTOCOL_20260918.md")]
    sources = {str(p.relative_to(base)): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
    if sources["results/NATIVE_PRESSURE_CONTROL_PROTOCOL_20260918.md"] != PROTOCOL_SHA:
        raise ValueError("Frozen control protocol differs")
    output.mkdir(exist_ok=False)
    job, rows = int(os.environ["SLURM_JOB_ID"]), []
    for block, modes in enumerate(ORDER):
        for mode in modes:
            directory = output / f"trial_{len(rows):02d}"
            try:
                row = trial(directory, block, mode, job)
            except Exception as error:
                row = dict(block=block, mode=mode, status="failed", error_type=type(error).__name__, error=str(error))
            save(directory / "trial.json", row)
            rows.append(row)
    for name, digest in sources.items():
        if hashlib.sha256((base / name).read_bytes()).hexdigest() != digest:
            raise ValueError("Control source changed")
    result = dict(job_id=job, sources=sources, trials=rows, summary=summarize(rows), publication_ready=False,
                  limitations=["CPU-only injected-demand controls, not observer overhead or scientific timing admission.",
                               "All failed controls retained; no selective repetition or threshold adjustment."])
    save(output / "result.json", result)
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--worker", type=Path)
    group.add_argument("--burn-on", type=int)
    group.add_argument("--output", type=Path)
    args = parser.parse_args()
    if args.worker is not None:
        worker(args.worker)
    elif args.burn_on is not None:
        print(json.dumps(work(args.burn_on)))
    else:
        result = run(args.output.resolve())
        raise SystemExit(0 if result["summary"]["all_response_checks_passed"] else 1)
