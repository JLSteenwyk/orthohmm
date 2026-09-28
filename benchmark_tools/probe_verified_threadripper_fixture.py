"""Run one installation fixture with full tree and lookup checks around inference."""

import argparse
import json
import os
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.check_threadripper_runtime import RuntimeChecker
from benchmark_tools.probe_threadripper_native_fixture import fixture_run
from benchmark_tools.run_simulation_methods import read_frozen, execution_environment
from benchmark_tools.measure_threadripper_run import measure_run
from benchmark_tools.measure_threadripper_scaling import measure
from benchmark_tools.validate_threadripper_outputs import validate
from benchmark_tools.prepare_ob_candidate_neighborhood import record
from benchmark_tools.probe_dgx_step_separation import save
from benchmark_tools.verify_threadripper_controller import ReleaseBudgetGuard

LOOKUP_SHA = "d5f26d4c31346f6d710ca911a5220a0448a0a58e1997d708ceb821a8e2a18e34"
PLAN_SHA = "c384e27730e3802b39ba14a42f7f50e84da5ce6deb9de9b2c32a74a745aed296"
METHODS = ("orthohmm_high_sensitivity", "orthohmm_satellite_v2", "orthofinder_full")


def select_fixture(plan, method):
    if method not in METHODS:
        raise ValueError("Unknown diagnostic method")
    rows = plan["runs"][:3]
    if ([r["native_method"] for r in rows] != list(METHODS)
            or [r["index"] for r in rows] != [0, 1, 2]):
        raise ValueError("First frozen method block differs")
    return rows[METHODS.index(method)]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("input", "output", "tmpfs"):
        parser.add_argument("--" + name, type=Path, required=True)
    parser.add_argument("--method", choices=METHODS, default=METHODS[0])
    parser.add_argument("--scheduler-command", required=True)
    args = parser.parse_args()
    allocation_cwd = str(Path.cwd())
    results = Path(__file__).resolve().parent / "results"
    lookup_path = results / "threadripper_python_lookup_v3_20260928.json"
    lookup = read_frozen(lookup_path, LOOKUP_SHA)
    baseline = read_frozen(Path(lookup["baseline"]["path"]), lookup["baseline"]["sha256"])
    binding = read_frozen(Path(lookup["binding"]["path"]), lookup["binding"]["sha256"])
    plan_path = results / "threadripper_scaling_commands_20260928.json"
    plan = read_frozen(plan_path, PLAN_SHA)
    args.output.mkdir(parents=True, exist_ok=False)
    env, _ = execution_environment(baseline)
    env.update(PYTHONDONTWRITEBYTECODE="1", PYTHONPYCACHEPREFIX=str(args.output / "python_cache"))
    os.environ.update(env)
    os.chdir(baseline["core_root"])
    checker = RuntimeChecker(lookup_path, LOOKUP_SHA, args.output / "lookup_checks")
    run = fixture_run(select_fixture(plan, args.method), args.input, args.output / "run_00", args.tmpfs)
    save(args.output / "command.json", run)
    save(args.output / "started.json", dict(source=record(__file__),
        checker_source=record(Path(__file__).with_name("check_threadripper_runtime.py")),
        lookup=record(lookup_path), plan=record(plan_path), method=args.method,
        scheduler_command=args.scheduler_command, allocation_cwd=allocation_cwd,
        release_budget_required=True,
        scientific_timings_admitted=False))
    job = int(os.environ["SLURM_JOB_ID"])
    guard = ReleaseBudgetGuard(job, command=args.scheduler_command, cwd=allocation_cwd)
    wrapped = measure_run(run, baseline, binding["runtime_specs"], measure,
                          job, runtime_checker=checker, release_guard=guard)
    result = dict(wrapper=wrapped, scientific_timings_admitted=False,
                  limitations=["One 16-gene fixture, not a production timing identity or executor authorization.",
                               "Full tree and declared-lookup checks do not establish quiet-host eligibility."])
    if wrapped["status"] == "command_exited_zero":
        directory = Path(run["measurement_directory"])
        done = json.loads((directory / "done.json").read_text())
        command = json.loads((directory / "command.json").read_text())["command"]
        result["outputs"] = validate(run, dict(done, command=command, cwd=os.getcwd()), baseline)
    save(args.output / "result.json", result)
    return 0 if wrapped["status"] == "command_exited_zero" else 1


if __name__ == "__main__":
    raise SystemExit(main())
