"""Validate retained counter-native smoke evidence without admitting timings."""

import argparse
import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.run_dgx_counter_native_smoke import ROOT, SPEC_SHA, relocate, read_pinned
from benchmark_tools.validate_scaling_outputs import validate
from benchmark_tools.prepare_ob_candidate_neighborhood import record, check
from benchmark_tools.probe_host_counters import summarize
from benchmark_tools.probe_dgx_step_separation import validate_scopes


def audit(archive, spec_path):
    spec = read_pinned(spec_path, SPEC_SHA)
    inventory = [record(p) for p in sorted(archive.rglob("*")) if p.is_file()]
    runs = []
    for i, original in enumerate(spec["runs"]):
        directory = archive / "counter_native_smoke_v1" / f"run_{i:02d}"
        prepared = json.loads((directory / "preparation.json").read_text())
        frozen = relocate(original)
        run = prepared["run"]
        if {k: v for k, v in run.items() if k != "gnu_time"} != frozen:
            raise ValueError("Native run differs from frozen commands")
        proof = json.loads((directory / "verification.json").read_text())
        measured = json.loads((directory / "measurement/counter_report.json").read_text())
        command = json.loads((directory / "measurement/command.json").read_text())
        # OrthoFinder copies are created after preflight and validated separately below.
        after_identity = {k: v for k, v in proof["after"].items() if k != "copied_inputs"}
        if proof["status"] != "command_exited_zero" or proof["before"] != after_identity or proof["measurement"] != measured:
            raise ValueError("Execution/runtime validation failed")
        if measured["counter_read_errors"] != 0 or measured["native"]["exit_code"] != 0:
            raise ValueError("Counter read or native process failure")
        if measured["controlled_workload_verified"] or measured["publication_ready"]:
            raise ValueError("Unexpected admission claim")
        before, after = measured["native_snapshots"]
        samples = [json.loads(line) for line in (directory / "measurement/observer.jsonl").read_text().splitlines()]
        if len(samples) != measured["observer_snapshots"] or any(s["errors"] for s in [before, after, *samples]):
            raise ValueError("Incomplete/error-bearing samples")
        if summarize(before, after, 100) != measured["host_summary"]:
            raise ValueError("Host summary replay differs")
        for sample in samples:
            validate_scopes(sample, before, measured["job_id"])
        done = measured["native"]
        if not before["finished_monotonic_ns"] < done["started_ns"] < done["finished_ns"] < after["started_monotonic_ns"]:
            raise ValueError("Work not bracketed")
        scheduler = (archive / f"scheduler_{i}.txt").read_text()
        for field in (f"JobId={measured['job_id']} ", f"ArrayTaskId={i} ", "JobState=COMPLETED ",
                      "ExitCode=0:0", "Restarts=0", "NodeList=spark-7ff0", "CPUs/Task=20", "MinMemoryNode=96G"):
            if field not in scheduler:
                raise ValueError("Scheduler evidence differs: " + field)
        # Adapt recorded native exit/argv to the existing output validator's schema.
        adapted = dict(command=command["command"], cwd=frozen["cwd"], exit_code=done["exit_code"], timed_out=False)
        checked = validate(run, adapted, {ROOT: archive})
        runs.append(dict(index=i, method=run["native_method"], scheduler=scheduler, verification=proof,
                         output_validation=checked))
    for item in inventory:
        check(item)
    return dict(status="three_counter_native_smokes_validated_not_timing_admission", source=record(__file__),
        specification=record(spec_path), archived_files=len(inventory), archive_bytes=sum(r["bytes"] for r in inventory),
        runs=runs, publication_ready=False, controlled_workload_verified=False, scientific_timings_admitted=0,
        limitations=["One small development-exposed fixture per method, not independent accuracy or scaling.",
            "Native outputs checked using existing group/pair/graph validators, not a new correctness oracle.",
            "No general observer overhead, workload isolation or exact tree equivalence is established.",
            "Input/runtime checks were outside the timed native command; raw archives remain local."])


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--archive", type=Path, required=True)
    parser.add_argument("--spec", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    result = audit(args.archive.resolve(), args.spec.resolve())
    with args.output.open("x") as stream:
        json.dump(result, stream, indent=2, sort_keys=True, allow_nan=False)
        stream.write("\n")
    print(result["status"])
