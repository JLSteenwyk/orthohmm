"""Trace three native fixture commands; tracing invalidates comparative timings."""

import argparse
import json
import os
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.probe_threadripper_native_fixture import fixture_run
from benchmark_tools.run_simulation_methods import read_frozen, execution_environment
from benchmark_tools.measure_threadripper_run import measure_run
from benchmark_tools.measure_threadripper_scaling import measure
from benchmark_tools.snapshot_runtime_trees import inventory, digest
from benchmark_tools.snapshot_orthohmm_input_order import record
from benchmark_tools.probe_dgx_step_separation import save
from benchmark_tools.validate_threadripper_outputs import validate


def traced_command(command, trace):
    return ["/usr/bin/strace", "-f", "-yy", "-s", "4096", "-e",
            "trace=%file,%process", "-o", str(trace), "--", *command]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("plan", "baseline", "input", "output", "tmpfs"):
        parser.add_argument("--" + name, type=Path, required=True)
    parser.add_argument("--plan-sha256", required=True)
    parser.add_argument("--baseline-sha256", required=True)
    args = parser.parse_args()
    plan = read_frozen(args.plan, args.plan_sha256)
    baseline = read_frozen(args.baseline, args.baseline_sha256)
    job = int(os.environ["SLURM_JOB_ID"])
    args.output.mkdir(parents=True, exist_ok=False)
    env, _ = execution_environment(baseline)
    env.update(PYTHONDONTWRITEBYTECODE="1", PYTHONPYCACHEPREFIX=str(args.output / "python_cache"))
    os.environ.update(env)
    os.chdir(baseline["core_root"])
    runtime = args.output / "core_tree.json"
    save(runtime, inventory([Path(baseline["core_root"]) / "orthohmm"]))
    result = dict(status="running", job_id=job, source=record(__file__),
        plan=record(args.plan), baseline=record(args.baseline), tracer=record("/usr/bin/strace"),
        scientific_timings_admitted=False, runs=[],
        limitations=["Traced 16-gene fixture only; strace changes runtime and resource use.",
                     "File/process syscalls, including failed lookups, not full syscall or loader isolation.",
                     "Output validation removes only the recorded exact tracer prefix."])
    save(args.output / "started.json", result)
    try:
        for index, method in enumerate(("orthohmm_high_sensitivity", "orthohmm_satellite_v2", "orthofinder_full")):
            root = args.output / f"run_{index:02d}"
            trace = root / "native.strace"
            template = next(r for r in plan["runs"] if r["native_method"] == method)
            run = fixture_run(template, args.input, root, args.tmpfs)
            save(args.output / f"command_{index:02d}.json", run)

            def collector(command, *arguments, **keywords):
                return measure(traced_command(command, trace), *arguments, **keywords)

            wrapped = measure_run(run, baseline, [(runtime, digest(runtime))], collector, job)
            row = dict(method=method, wrapper=record(root / "verification.json"), status=wrapped["status"])
            result["runs"].append(row)
            directory = Path(run["measurement_directory"])
            actual = json.loads((directory / "command.json").read_text())["command"]
            if actual != traced_command(run["native_argv"], trace):
                raise ValueError("Recorded traced command differs")
            done = json.loads((directory / "done.json").read_text())
            row.update(trace=record(trace), traced_command=actual)
            row["outputs"] = validate(run, dict(done, command=run["native_argv"], cwd=os.getcwd()), baseline)
            save(args.output / f"result_{index:02d}.json", row)
            if wrapped["status"] != "command_exited_zero":
                raise RuntimeError("Traced native run or runtime verification failed")
        result["status"] = "traced_native_fixture_completed"
    except Exception as error:
        result.update(status="traced_native_fixture_failed", error_type=type(error).__name__, error=str(error))
        raise
    finally:
        save(args.output / "result.json", result)


if __name__ == "__main__":
    main()
