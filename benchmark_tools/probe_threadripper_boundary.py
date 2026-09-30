"""One pinned /usr/bin/true boundary integration fixture, never tool timing."""

import argparse
import json
import os
from pathlib import Path
import sys

from benchmark_tools import measure_threadripper_scaling as periodic
from benchmark_tools.measure_threadripper_boundary import measure
from benchmark_tools.probe_dgx_step_separation import save
from benchmark_tools.replay_threadripper_boundary import replay
from benchmark_tools.run_integrated_publication_workflow import check, record
from benchmark_tools.summarize_fixture_memory_scopes import native_cpu, peak
from benchmark_tools.verify_lineage_native_provenance import same

SETTINGS = dict(cpus=32, memory_bytes=128 * 1024**3, timeout_s=85800,
                interval_s=1., monitor_host=True, host_interval_s=30.)
SCHEDULING = dict(node="bizon", partition="gpu", cpus_per_task=64,
                  memory_gib=128, exclusive=True, time_limit="00:05:00",
                  requeue=False, attempts=1, automatic_retry=False)
ADMISSIONS = ("scientific_timings_admitted", "controlled_workload_verified",
              "native_outputs_validated", "publication_ready")


def protocol(path, digest):
    path = Path(path).absolute()
    if path.is_symlink() or path.resolve() != path or not path.is_file():
        raise ValueError("Require direct prospective protocol")
    pin = record(path)
    if pin["sha256"] != digest:
        raise ValueError("Boundary fixture protocol digest differs")
    value = json.loads(path.read_text())
    if (value.get("schema") != "threadripper_boundary_fixture_v1"
            or value.get("command") != ["/usr/bin/true"]
            or not same(value.get("settings"), SETTINGS)
            or not same(value.get("scheduling"), SCHEDULING)
            or value.get("release_guard") is not None
            or value.get("production_execution_authorized") is not False
            or any(value.get(k) is not False for k in ADMISSIONS)):
        raise ValueError("Boundary fixture policy differs")
    root = Path(__file__).resolve().parent.parent
    if value.get("cwd") != str(root):
        raise ValueError("Fixture working directory differs")
    directory = Path(value["run_directory"])
    if (not directory.is_absolute() or directory.resolve() != directory
            or directory.parent != root / "benchmarks/work"):
        raise ValueError("Require direct fixture directory under local work root")
    current = [record(p) for p in sorted((root / "benchmark_tools").glob("*.py"))]
    if not same(value.get("sources"), current):
        raise ValueError("Complete current helper inventory differs")
    if value["native_executable"]["path"] != "/usr/bin/true":
        raise ValueError("Native executable binding differs")
    for item in (pin, value["interpreter"], value["native_executable"],
                 value["submission_script"], value["controller_provenance"], value["policy"]):
        check(item)
    return pin, value


def audit(directory, job_id, protocol_path, digest, output):
    directory, output = Path(directory).absolute(), Path(output).absolute()
    if output.exists() or output.is_symlink():
        raise FileExistsError(output)
    pin, value = protocol(protocol_path, digest)
    if directory != Path(value["run_directory"]):
        raise ValueError("Fixture directory differs from protocol")
    if type(job_id) is not int or job_id <= 0:
        raise ValueError("Require external positive expected job identity")
    if (directory / "fixture_failure.json").exists():
        raise ValueError("Failed fixture cannot pass integration audit")
    started_path = directory / "started.json"
    if started_path.is_symlink() or not started_path.is_file():
        raise ValueError("Missing direct fixture launch receipt")
    started_pin = record(started_path)
    started = json.loads(started_path.read_text())
    expected = dict(schema="threadripper_boundary_fixture_started_v1", job_id=job_id,
        command=value["command"], protocol=pin, source=record(__file__),
        interpreter=value["interpreter"], release_guard=None, production_timing=False)
    if not same(started, expected):
        raise ValueError("Fixture launch receipt differs from external bindings")
    raw = replay(directory / "measurement", job_id, value["command"],
        expected_launcher=value["interpreter"]["path"],
        expected_worker=str(Path(periodic.__file__).resolve()))
    if (raw["native_outcome"] != "exited_zero" or raw["native_exit_code"] != 0
            or any(raw[k] is not False for k in ADMISSIONS)):
        raise ValueError("Fixture native execution failed or timing was admitted")
    measured = raw["measured"]
    cpu = native_cpu(measured["native"], job_id, raw["native_completion"])
    memory = peak(measured["step_memory"])
    result = dict(schema="threadripper_boundary_fixture_audit_v1",
        status="native_boundary_component_replayed_pending_terminal_audit",
        job_id=job_id, protocol=pin, source=record(__file__),
        native_exit_code=0, native_timed_out=False, native_points=len(measured["point_records"]),
        native_wall_s=raw["native_wall_s"], cpu=cpu, step_memory=memory,
        host_process_observation=measured["host_process_observation"],
        evidence=[started_pin, *raw["evidence"]],
        **{k: False for k in ADMISSIONS},
        limitations=["Only /usr/bin/true integration; not OrthoHMM/OrthoFinder inference or slowdown validation.",
            "Synthetic tests are complemented by native receipts; terminal scheduler provenance remains separate.",
            "No environmental release guard or prospective contamination review is tested.",
            "Five-minute diagnostic allocation does not validate the full native timeout or reporting reserve.",
            "CPU bracket includes wrapper work; step lifetime peak includes launcher; no subtraction or addition.",
            "Interpreter/helper pins and retained controller provenance are not complete loaded-runtime closure.",
            "Two boundary points do not prove continuous affinity, containment, isolation or total monitor cost."])
    for item in result["evidence"]:
        check(item)
    protocol(protocol_path, digest)
    save(output, result)
    return result


def run(output, protocol_path, digest):
    pin, value = protocol(protocol_path, digest)
    output = Path(output).absolute()
    if (output != Path(value["run_directory"]) or Path.cwd() != Path(value["cwd"])
            or record(sys.executable) != value["interpreter"]
            or any(os.environ.get(k) for k in ("LD_PRELOAD", "LD_LIBRARY_PATH", "LD_AUDIT"))
            or os.uname().nodename != "bizon"
            or os.environ.get("SLURM_CPUS_PER_TASK") != "64"
            or os.environ.get("SLURM_MEM_PER_NODE") != "131072"):
        raise ValueError("Require pinned local runtime, cwd and 64-slot/128GiB allocation")
    job_id = int(os.environ["SLURM_JOB_ID"])
    if job_id <= 0:
        raise ValueError("Invalid allocation identity")
    output.mkdir(exist_ok=False)
    save(output / "started.json", dict(schema="threadripper_boundary_fixture_started_v1",
        job_id=job_id, command=value["command"], protocol=pin, source=record(__file__),
        interpreter=value["interpreter"], release_guard=None, production_timing=False))
    try:
        measure(value["command"], output / "measurement", job_id, **SETTINGS)
        return audit(output, job_id, protocol_path, digest, output / "component_audit.json")
    except BaseException as error:
        save(output / "fixture_failure.json", dict(status="boundary_fixture_failed",
            job_id=job_id, error_type=type(error).__name__, error=str(error),
            scientific_timings_admitted=False, publication_ready=False))
        raise


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--protocol", type=Path, required=True)
    parser.add_argument("--protocol-sha256", required=True)
    parser.add_argument("--audit-directory", type=Path)
    parser.add_argument("--job", type=int)
    args = parser.parse_args()
    if args.audit_directory is not None:
        result = audit(args.audit_directory, args.job, args.protocol, args.protocol_sha256, args.output)
    else:
        if args.job is not None:
            parser.error("External --job is only valid with --audit-directory")
        result = run(args.output, args.protocol, args.protocol_sha256)
    print(result["status"])
