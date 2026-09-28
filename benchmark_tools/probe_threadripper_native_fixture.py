"""Exercise all three native timing paths on the existing installation fixture.

Diagnostic only: shared-host times and cumulative job peaks are not benchmarks.
"""

import argparse
from copy import deepcopy
import json
import os
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from Bio import SeqIO
from benchmark_tools.measure_threadripper_run import measure_run
from benchmark_tools.measure_threadripper_scaling import measure
from benchmark_tools.prepare_threadripper_run import paths
from benchmark_tools.probe_dgx_step_separation import save
from benchmark_tools.run_simulation_methods import read_frozen, execution_environment
from benchmark_tools.snapshot_orthohmm_input_order import record
from benchmark_tools.snapshot_runtime_trees import inventory, digest
from benchmark_tools.validate_threadripper_outputs import validate


def fixture_run(template, source, root, tmpfs):
    files = sorted(source.glob("*.fa"))
    sequences = [s for p in files for s in SeqIO.parse(p, "fasta")]
    if len(files) != 4 or len(sequences) != 16 or len({s.id for s in sequences}) != 16:
        raise ValueError("Require the four-species, 16-gene installation fixture")
    run = deepcopy(template)
    old_root = str(Path(run["measurement_directory"]).parent)
    old_input = run["prepared_input_directory"]
    target = str(tmpfs / root.name / "input")

    def relocate(value):
        if value == old_input:
            return target
        if value.startswith(old_root + "/"):
            return str(root / Path(value).relative_to(old_root))
        return value

    run["native_argv"] = list(map(relocate, run["native_argv"]))
    config = run["configuration"]
    for key in ("output", "metrics"):
        if key in config:
            config[key] = relocate(config[key])
    config["argv"] = list(map(relocate, config["argv"]))
    config.update(copy_inputs_from=str(source), copy_inputs_to=target)
    run.update(measurement_directory=str(root / "measurement"),
               prepared_input_directory=target, proteomes=4,
               dataset=dict(input_directory=str(source), inputs=[record(p) for p in files],
                            proteins=16, proteomes=4),
               input_creation_order=[p.name for p in files],
               expected_native_order=[p.name for p in files],
               diagnostic_only=True)
    paths(run)
    return run


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
    env, resolved = execution_environment(baseline)
    # Excluded bytecode must not be loaded from historical runtime caches.
    cache = args.output / "python_cache"
    cache.mkdir()
    env.update(PYTHONDONTWRITEBYTECODE="1", PYTHONPYCACHEPREFIX=str(cache))
    os.environ.update(env)
    os.chdir(baseline["core_root"])
    runtime = args.output / "core_tree.json"
    save(runtime, inventory([Path(baseline["core_root"]) / "orthohmm"]))
    result = dict(status="running", job_id=job, source=record(__file__),
        plan=record(args.plan), baseline=record(args.baseline), resolved=resolved,
        scientific_timings_admitted=False, accuracy_evaluated=False, runs=[],
        limitations=["Diagnostic fixture, not any of the 27 timing identities.",
                     "Core tree and existing baseline checks, not a complete transitive runtime freeze.",
                     "Multiple controls share one job; later job peaks include earlier controls.",
                     "Shared-host collection does not establish isolation."])
    save(args.output / "started.json", result)
    try:
        for index, method in enumerate(("orthohmm_high_sensitivity", "orthohmm_satellite_v2", "orthofinder_full")):
            template = next(r for r in plan["runs"] if r["native_method"] == method)
            run = fixture_run(template, args.input, args.output / f"run_{index:02d}", args.tmpfs)
            save(args.output / f"command_{index:02d}.json", run)
            wrapped = measure_run(run, baseline, [(runtime, digest(runtime))], measure, job)
            row = dict(method=method, wrapper=wrapped)
            result["runs"].append(row)
            directory = Path(run["measurement_directory"])
            try:
                measurement = json.loads((directory / "done.json").read_text())
                command = json.loads((directory / "command.json").read_text())
                measurement.update(command=command["command"], cwd=os.getcwd())
                row["outputs"] = validate(run, measurement, baseline)
            except Exception as error:
                row["output_error"] = dict(type=type(error).__name__, message=str(error))
            save(args.output / f"result_{index:02d}.json", row)
            if wrapped["status"] in {"verified_wrapper_failed", "runtime_changed_or_unverifiable"}:
                raise RuntimeError("Runtime or infrastructure verification failed")
        result["status"] = "diagnostic_finished"
    except Exception as error:
        result.update(status="diagnostic_failed", error_type=type(error).__name__, error=str(error))
        raise
    finally:
        save(args.output / "result.json", result)


if __name__ == "__main__":
    main()
