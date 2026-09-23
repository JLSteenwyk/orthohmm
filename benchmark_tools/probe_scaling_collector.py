"""Short native command under the full long-run collector, without timing admission."""

import argparse
import os
from pathlib import Path

from benchmark_tools.measure_scaling_root_context import measure, TIMEOUT
from benchmark_tools.replay_scaling_root_context import replay
from benchmark_tools.probe_dgx_step_separation import save


def run(output):
    output.mkdir(exist_ok=False)
    job = int(os.environ["SLURM_JOB_ID"])
    command = ["/bin/sleep", "35"]
    directory = output / "measurement"
    measured = measure(command, directory, job, 20, 96 * 1024**3, TIMEOUT, 1.)
    audited = replay(directory, job, command)
    save(output / "replay.json", audited)
    host = audited["host_process_replay"]
    if (measured["native"]["exit_code"] != 0 or measured["native"]["timed_out"]
            or host is None or host["successful_snapshots"] < 3
            or not host["command_bracketed_by_samples"]):
        raise ValueError("Short collector diagnostic failed or lacks bracketing/periodic evidence")
    result = dict(status="collector_probe_completed", job_id=job,
        native_wall_s=measured["native_wall_s"], host_process_replay=host,
        scientific_timings_admitted=False, environmental_validity_established=False)
    save(output / "probe.json", result)
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    run(parser.parse_args().output.resolve())
