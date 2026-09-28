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

LOOKUP_SHA = "07eecb538551c72b742e89cdec3691bd95c33d4caa3bad09e7c9f549baaaac09"
PLAN_SHA = "c384e27730e3802b39ba14a42f7f50e84da5ce6deb9de9b2c32a74a745aed296"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("input", "output", "tmpfs"):
        parser.add_argument("--" + name, type=Path, required=True)
    args = parser.parse_args()
    results = Path(__file__).resolve().parent / "results"
    lookup_path = results / "threadripper_python_lookup_20260928.json"
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
    run = fixture_run(plan["runs"][0], args.input, args.output / "run_00", args.tmpfs)
    if run["native_method"] != "orthohmm_high_sensitivity":
        raise ValueError("Unexpected first frozen method")
    save(args.output / "command.json", run)
    save(args.output / "started.json", dict(source=record(__file__),
        checker_source=record(Path(__file__).with_name("check_threadripper_runtime.py")),
        lookup=record(lookup_path), plan=record(plan_path), scientific_timings_admitted=False))
    wrapped = measure_run(run, baseline, binding["runtime_specs"], measure,
                          int(os.environ["SLURM_JOB_ID"]), runtime_checker=checker)
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
