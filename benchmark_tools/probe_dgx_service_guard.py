"""Exercise the service guard before submission; never launch a benchmark/job."""

import argparse
import hashlib
import json
import os
from pathlib import Path
import signal
import time

from benchmark_tools import dgx_service_guard
from benchmark_tools.probe_dgx_step_separation import save


def run(output, interrupt=False):
    guard = dgx_service_guard.ServiceGuard(output)
    paths = [Path(__file__).resolve(), Path(dgx_service_guard.__file__).resolve(),
             Path(__file__).resolve().with_name("probe_dgx_step_separation.py"),
             Path(__file__).resolve().with_name("probe_host_counters.py")]
    sources = {str(path): hashlib.sha256(path.read_bytes()).hexdigest() for path in paths}
    result = {"status": "probe_started", "sources": sources, "interrupt_requested": interrupt,
              "benchmark_submitted": False, "scientific_timings_admitted": False}

    def stop(signum, frame):
        raise InterruptedError(f"Service probe signal {signum}")

    previous = {sig: signal.signal(sig, stop) for sig in (signal.SIGHUP, signal.SIGINT, signal.SIGTERM)}
    try:
        guard.begin()
        for index in range(3):
            guard.check()
            save(output / f"sample_{index}.json", {"suppression_verified": True,
                 "observed_unix_ns": time.time_ns()})
            if interrupt and index == 1:
                os.kill(os.getpid(), signal.SIGTERM)
            if index != 2:
                time.sleep(2)
        result["status"] = "normal_probe_completed"
    except InterruptedError as error:
        result.update(status="probe_interrupted", error=str(error))
    except BaseException as error:
        result.update(status="probe_failed", error_type=type(error).__name__, error=str(error))
        raise
    finally:
        try:
            # No submission attempt exists in this probe, so prelaunch cleanup
            # does not need or fabricate a terminal Slurm receipt.
            if guard.stopped:
                guard.restore()
            result["restoration_verified"] = guard.restored
            if sources != {str(path): hashlib.sha256(path.read_bytes()).hexdigest() for path in paths}:
                raise ValueError("Probe source changed during execution")
        except BaseException as error:
            result.update(restoration_verified=False, cleanup_error_type=type(error).__name__,
                          cleanup_error=str(error))
            raise
        finally:
            for sig, handler in previous.items():
                signal.signal(sig, handler)
            if output.is_dir():
                save(output / "probe.json", result)
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--interrupt", action="store_true")
    args = parser.parse_args()
    report = run(args.output.resolve(), args.interrupt)
    print(json.dumps(report, sort_keys=True))
