"""Replay complete-command hierarchy smokes and validate native outputs."""

import argparse
import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.run_dgx_hierarchy_native_smoke import ROOT, SPEC_SHA, relocate, read_pinned
from benchmark_tools.measure_native_hierarchy_step import evaluate, interval_point
from benchmark_tools.validate_scaling_outputs import validate
from benchmark_tools.prepare_ob_candidate_neighborhood import record, check


def audit(archive, spec_path):
    spec = read_pinned(spec_path, SPEC_SHA)
    inventory = [record(p) for p in sorted(archive.rglob("*")) if p.is_file()]
    runs = []
    for i, original in enumerate(spec["runs"]):
        directory = archive / "hierarchy_native_smoke_v1" / f"run_{i:02d}"
        prepared = json.loads((directory / "preparation.json").read_text())
        frozen = relocate(original)
        run = prepared["run"]
        if {k: v for k, v in run.items() if k != "gnu_time"} != frozen:
            raise ValueError("Native run differs from frozen commands")
        proof = json.loads((directory / "verification.json").read_text())
        measured = json.loads((directory / "measurement/hierarchy_report.json").read_text())
        command = json.loads((directory / "measurement/command.json").read_text())
        after_identity = {k: v for k, v in proof["after"].items() if k != "copied_inputs"}
        if proof["status"] != "command_exited_zero" or proof["before"] != after_identity or proof["measurement"] != measured:
            raise ValueError("Execution/runtime validation failed")
        done = measured["native"]
        if done["exit_code"] != 0 or done["timed_out"]:
            raise ValueError("Native process failure")
        if any(measured[k] for k in ("controlled_workload_verified", "publication_ready", "scientific_timings_admitted")):
            raise ValueError("Unexpected timing admission")
        if evaluate(measured["points"], done, measured["job_id"]) != measured["screening"]:
            raise ValueError("Interval or command-boundary replay differs")
        raw_points = [json.loads(p.read_text()) for p in sorted((directory / "measurement").glob("point_*.json"))]
        if raw_points != measured["points"]:
            raise ValueError("Raw observation inventory differs")
        for name, expected in (("done.json", done), ("step_memory.json", measured["step_memory"])):
            if json.loads((directory / "measurement" / name).read_text()) != expected:
                raise ValueError("Raw worker or memory evidence differs")
        memory = measured["step_memory"]
        if memory["errors"] or memory["scope"] != interval_point(measured["points"][-1], measured["job_id"])["native_cpu_scope"]:
            raise ValueError("Memory scope or read failure")
        current, peak = (int(memory["raw"][key]) for key in ("memory.current", "memory.peak"))
        if not 0 <= current <= peak or not done["finished_ns"] < memory["started_ns"] <= memory["finished_ns"]:
            raise ValueError("Invalid final memory evidence")
        scheduler = (archive / f"scheduler_{i}.txt").read_text()
        for field in (f"JobId={measured['job_id']} ", f"ArrayTaskId={i} ", "JobState=COMPLETED ",
                      "ExitCode=0:0", "Restarts=0", "NodeList=spark-7ff0", "CPUs/Task=20", "MinMemoryNode=96G"):
            if field not in scheduler:
                raise ValueError("Scheduler evidence differs: " + field)
        adapted = dict(command=command["command"], cwd=frozen["cwd"], exit_code=done["exit_code"], timed_out=False)
        checked = validate(run, adapted, {ROOT: archive})
        runs.append(dict(index=i, method=run["native_method"], scheduler=scheduler, verification=proof,
                         output_validation=checked))
    for item in inventory:
        check(item)
    return dict(status="three_hierarchy_native_smokes_validated_not_timing_admission", source=record(__file__),
        specification=record(spec_path), archived_files=len(inventory), archive_bytes=sum(r["bytes"] for r in inventory),
        runs=runs, publication_ready=False, controlled_workload_verified=False, scientific_timings_admitted=0,
        limitations=["One small development-exposed fixture per method, not accuracy or scaling.",
            "Native outputs use existing partition/pair/graph validators, not a new correctness oracle.",
            "Diagnostic screening results are retained, not used to selectively suppress completed runs.",
            "No general observer overhead, accounting bound, non-CPU isolation or scientific inclusion policy.",
            "Final step memory includes wrapper and cache and is not process RSS."])


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
