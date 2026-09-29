"""Bounded worker/thread calibration for the local collector, with witnesses."""

import argparse
import json
import os
from pathlib import Path
import resource
import subprocess
import sys
import threading
import time

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.probe_dgx_step_separation import save, wait_file
from benchmark_tools.run_integrated_full_job import record
from benchmark_tools.observe_thread_affinity import identity as proc_identity


def validate(workers, threads, arena_mib, duration):
    if any(type(v) is not int or v < 1 for v in (workers, threads, arena_mib, duration)):
        raise ValueError("Require positive integer workload dimensions")
    if workers > 32 or threads > 4 or arena_mib > 64 or duration > 30:
        raise ValueError("Workload exceeds the fixed calibration envelope")


def worker(directory, index, threads, arena_mib, duration):
    validate(1, threads, arena_mib, duration)
    if type(index) is not int or not 0 <= index < 32:
        raise ValueError("Invalid worker index")
    affinity = sorted(os.sched_getaffinity(0))
    if affinity != list(range(32)):
        raise ValueError("Require inherited CPUs 0-31")
    arena = bytearray(arena_mib * 1024**2)
    arena[::4096] = b"\x01" * len(arena[::4096])
    barrier = threading.Barrier(threads + 1, timeout=20)
    stop = threading.Event()
    identities, errors = [], []

    def target():
        try:
            tid = threading.get_native_id()
            identities.append(dict(tid=tid, affinity=sorted(os.sched_getaffinity(0)),
                                   start_ticks=proc_identity(Path(f"/proc/{tid}/stat"))[1],
                                   membership=Path(f"/proc/{tid}/cgroup").read_text()))
            barrier.wait()
            deadline = time.monotonic() + 45
            gate = directory / "workers_go.json"
            while not gate.exists():
                if stop.wait(.02):
                    return
                if time.monotonic() >= deadline:
                    raise TimeoutError(str(gate))
            if json.loads(gate.read_text()) != {"go": True}:
                raise ValueError("Invalid thread release")
            value = 1
            while not stop.is_set():
                for _ in range(1000):
                    value = (1664525 * value + 1013904223) & 0xffffffff
        except BaseException as error:
            errors.append(str(error))
            barrier.abort()
            stop.set()

    running = [threading.Thread(target=target, daemon=True) for _ in range(threads)]
    try:
        for thread in running:
            thread.start()
        barrier.wait()
        identity = dict(pid=os.getpid(), index=index, affinity=affinity,
                        start_ticks=proc_identity(Path("/proc/self/stat"))[1], threads=identities,
                        membership=Path("/proc/self/cgroup").read_text(), arena_bytes=len(arena))
        save(directory / f"ready_{index:02d}.json", identity)
        if wait_file(directory / "workers_go.json", seconds=45) != {"go": True}:
            raise ValueError("Invalid calibration release")
        before = resource.getrusage(resource.RUSAGE_SELF)
        started = time.monotonic_ns()
        stop.wait(duration)
        stop.set()
        for thread in running:
            thread.join(timeout=10)
        if errors or any(thread.is_alive() for thread in running):
            raise RuntimeError("Worker thread failed or did not finish: " + repr(errors))
        finished = time.monotonic_ns()
        after = resource.getrusage(resource.RUSAGE_SELF)
        save(directory / f"done_{index:02d}.json", dict(identity=identity,
            started_ns=started, finished_ns=finished,
            user_seconds=after.ru_utime - before.ru_utime,
            system_seconds=after.ru_stime - before.ru_stime,
            final_membership=Path("/proc/self/cgroup").read_text(),
            final_affinity=sorted(os.sched_getaffinity(0))))
    finally:
        stop.set()
        for thread in running:
            if thread.ident is not None:
                thread.join(timeout=10)


def workload(directory, workers, threads, arena_mib, duration):
    validate(workers, threads, arena_mib, duration)
    directory.mkdir(parents=True, exist_ok=False)
    children = []
    try:
        for index in range(workers):
            command = [sys.executable, "-B", str(Path(__file__).resolve()), "--worker-index", str(index),
                       "--output", str(directory), "--threads", str(threads),
                       "--arena-mib", str(arena_mib), "--duration", str(duration)]
            with (directory / f"worker_{index:02d}.log").open("x") as log:
                children.append(subprocess.Popen(command, stdout=log, stderr=subprocess.STDOUT))
        deadline = time.monotonic() + 35
        ready = []
        for index in range(workers):
            path = directory / f"ready_{index:02d}.json"
            while not path.exists():
                if any(child.poll() is not None for child in children):
                    raise RuntimeError("A worker exited before readiness")
                if time.monotonic() >= deadline:
                    raise TimeoutError(str(path))
                time.sleep(.02)
            ready.append(json.loads(path.read_text()))
        if any(child.poll() is not None for child in children):
            raise RuntimeError("A worker exited before release")
        save(directory / "workers_go.json", {"go": True})
        deadline = time.monotonic() + duration + 45
        for child in children:
            if child.wait(timeout=max(.01, deadline - time.monotonic())) != 0:
                raise RuntimeError("Calibration worker failed")
        witnesses = [json.loads((directory / f"done_{i:02d}.json").read_text()) for i in range(workers)]
        if [w["identity"] for w in witnesses] != ready:
            raise ValueError("Worker identities changed")
        result = dict(workers=workers, threads_per_worker=threads, arena_mib_per_worker=arena_mib,
                      duration_seconds=duration, witnesses=witnesses,
                      parent_pid=os.getpid(), source=record(__file__), scientific_timings_admitted=False)
        save(directory / "complete.json", result)
        return result
    finally:
        for child in children:
            if child.poll() is None:
                child.terminate()
        for child in children:
            try:
                child.wait(timeout=10)
            except subprocess.TimeoutExpired:
                child.kill()
                child.wait(timeout=10)


def run(output, protocol_path, protocol_sha):
    from benchmark_tools.measure_threadripper_scaling import measure, TIMEOUT
    from benchmark_tools.derive_threadripper_resources import derive
    from benchmark_tools.audit_threadripper_observer import verify_protocol
    if os.environ.get("SLURM_CPUS_PER_TASK") != "64" or os.uname().nodename != "bizon":
        raise ValueError("Require the local 64-slot diagnostic allocation")
    protocol = verify_protocol(protocol_path, protocol_sha)
    if record(sys.executable) != protocol["interpreter"]:
        raise ValueError("Diagnostic interpreter differs from the frozen protocol")
    output.mkdir(parents=True, exist_ok=False)
    job = int(os.environ["SLURM_JOB_ID"])
    command = [sys.executable, "-B", str(Path(__file__).resolve()), "--workload",
               "--output", str(output / "witnesses")]
    save(output / "started.json", dict(job_id=job, source=record(__file__), command=command,
                                      protocol=record(protocol_path),
                                      production_timing=False))
    measure(command, output / "measurement", job, 32, 128 * 1024**3, TIMEOUT, 1.)
    derive(output / "measurement", job, Path(protocol["resource_protocol"]["path"]),
           protocol["resource_protocol"]["sha256"],
           output / "scoped_resources.json")
    verify_protocol(protocol_path, protocol_sha)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--workload", action="store_true")
    parser.add_argument("--worker-index", type=int)
    parser.add_argument("--workers", type=int, default=32)
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument("--arena-mib", type=int, default=64)
    parser.add_argument("--duration", type=int, default=30)
    parser.add_argument("--protocol", type=Path)
    parser.add_argument("--protocol-sha256")
    args = parser.parse_args()
    if args.worker_index is not None:
        worker(args.output.resolve(), args.worker_index, args.threads, args.arena_mib, args.duration)
    elif args.workload:
        workload(args.output.resolve(), args.workers, args.threads, args.arena_mib, args.duration)
    else:
        if args.protocol is None or args.protocol_sha256 is None:
            parser.error("Diagnostic launch requires a pinned prospective protocol")
        run(args.output.resolve(), args.protocol.resolve(), args.protocol_sha256)
