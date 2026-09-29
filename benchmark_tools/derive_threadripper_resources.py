"""Derive prospectively selected resource endpoints from a checked raw replay."""

import argparse
import json
from pathlib import Path

from benchmark_tools.run_integrated_full_job import record, save
from benchmark_tools.run_integrated_publication_workflow import check
from benchmark_tools.replay_threadripper_scaling import replay
from benchmark_tools.summarize_fixture_memory_scopes import native_cpu, summarize

SCOPES = {
    "wall_seconds": "native_command_monotonic_interval",
    "cpu_seconds": "native_task_subtree_cpu_stat_bracket_including_wrapper",
    "peak_memory_bytes": "native_step_lifetime_memory_peak_including_launcher",
}


def protocol(path, digest):
    pin = record(path)
    if pin["sha256"] != digest:
        raise ValueError("Resource protocol digest differs")
    value = json.loads(path.read_text())
    if (value["schema"] != "threadripper_resource_endpoints_v1"
            or value["primary_scopes"] != SCOPES
            or value["whole_job_teardown_required"] is not False
            or value["production_execution_authorized"] is not False):
        raise ValueError("Resource endpoint protocol differs")
    for item in value["sources"]:
        check(item)
    current = {record(Path(fn.__code__.co_filename))["sha256"] for fn in (replay, native_cpu, summarize)}
    if not current <= {p["sha256"] for p in value["sources"]}:
        raise ValueError("Resource protocol omits current derivation helpers")
    return pin, value


def endpoints(replayed, done, job_id):
    if replayed["measured"]["job_id"] != job_id or replayed["measured"]["native"] != done:
        raise ValueError("Replay and native job/command receipts differ")
    memory = summarize(replayed)
    cpu = native_cpu(done, job_id, replayed["native_completion"])
    if cpu["command_wall_seconds"] != replayed["native_wall_s"]:
        raise ValueError("CPU bracket and command wall differ")
    return dict(primary=dict(wall_seconds=replayed["native_wall_s"], cpu_seconds=cpu["cpu_seconds"],
                            peak_memory_bytes=memory["measurements"]["native_step"]["bytes"]),
                primary_scopes=SCOPES, cpu=cpu, memory=memory,
                native_outcome=replayed["native_outcome"], native_exit_code=replayed["native_exit_code"],
                controlled_timing_admitted=False, publication_ready=False)


def derive(directory, job_id, protocol_path, protocol_sha, output):
    if output.exists():
        raise FileExistsError(output)
    pin, _ = protocol(protocol_path, protocol_sha)
    command_path = directory / "command.json"
    command_pin = record(command_path)
    command = json.loads(command_path.read_text())["command"]
    replayed = replay(directory, job_id, command)
    done_path = directory / "done.json"
    done_pin = record(done_path)
    if done_pin not in replayed["evidence"]:
        raise ValueError("Native receipt is not bound to replay")
    result = endpoints(replayed, json.loads(done_path.read_text()), job_id)
    for item in [pin, command_pin, done_pin, *replayed["evidence"]]:
        check(item)
    protocol(protocol_path, protocol_sha)
    result.update(status="scoped_resources_derived_pending_timing_admission", job_id=job_id,
                  protocol=pin, source=record(__file__), evidence=replayed["evidence"],
                  final_whole_job_cpu_seconds=None, final_whole_job_peak_bytes=None)
    save(output, result)
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--directory", type=Path, required=True)
    parser.add_argument("--job", type=int, required=True)
    parser.add_argument("--protocol", type=Path, required=True)
    parser.add_argument("--protocol-sha256", required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    print(derive(args.directory.resolve(), args.job, args.protocol.resolve(),
                 args.protocol_sha256, args.output.resolve())["status"])
