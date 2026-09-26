"""Two fixed import-context controls before frozen high-CPM partition parsing."""

import argparse
import importlib.util
import json
import os
from pathlib import Path
import resource
import subprocess
import sys
import time

# Load the stdlib-only helper without caching the development benchmark_tools package.
spec = importlib.util.spec_from_file_location("parser_control_helper", Path(__file__).with_name("probe_cpm_partition_parser.py"))
helper = importlib.util.module_from_spec(spec)
spec.loader.exec_module(helper)
PRIOR_SHA = "ef73bd0badf9752345d6e7e8a02cb4fa8a1568bac554bd99f30eeefcba589a93"
PROTOCOL = "benchmark_tools/results/QFO_CPM_PARSER_IMPORT_PROTOCOL_20260926.md"


def scientific_imports(launcher):
    sys.path.insert(0, str(launcher))
    from benchmark_tools import replay_high_sensitivity as replay
    from benchmark_tools.audit_historical_profile_ablation import read_partition
    import orthohmm.accuracy
    import orthohmm.refinement
    modules = [replay, orthohmm.accuracy, orthohmm.refinement,
               sys.modules[read_partition.__module__], sys.modules[replay.audit_numeric_checkpoint.__module__]]
    if any(not Path(m.__file__).resolve().is_relative_to(launcher) for m in modules):
        raise ValueError("Nonfrozen scientific import")
    return [helper.record(m.__file__) for m in modules]


def worker(root, mode):
    resource.setrlimit(resource.RLIMIT_AS, (4 * 1024**3, 4 * 1024**3))
    resource.setrlimit(resource.RLIMIT_CPU, (120, 120))
    resource.setrlimit(resource.RLIMIT_CORE, (0, 0))
    inputs = helper.pinned_inputs(root)
    print("before_imports", file=sys.stderr, flush=True)
    modules = scientific_imports(root / "benchmarks/work/publication_qfo_replay_native_v1") if mode == "frozen_imports" else []
    print("before_parser", file=sys.stderr, flush=True)
    result = helper.parse_only(*(Path(r["path"]) for r in inputs))
    print("after_parser", file=sys.stderr, flush=True)
    helper.check([*inputs, *modules])
    if (result["genes"], result["groups"], result["memberships"]) != (984137, 390845, 984137):
        raise ValueError("Wrong parser coverage")
    result.update(mode=mode, scientific_sources=modules)
    print(json.dumps(result, sort_keys=True))


def run(root, output):
    if output.exists() or output.is_symlink():
        raise FileExistsError(output)
    prior_path = root / "benchmark_tools/results/qfo_cpm_parser_controls_20260926.json"
    prior_record = helper.record(prior_path)
    if prior_record["sha256"] != PRIOR_SHA:
        raise ValueError("Changed parser-only controls")
    prior = json.loads(prior_path.read_text())
    original = json.loads((root / helper.STATUS).read_text())
    launcher = root / "benchmarks/work/publication_qfo_replay_native_v1"
    python = original["scientific_child_command"][0]
    inputs = [prior_record, *prior["checked_records"], helper.record(__file__), helper.record(root / PROTOCOL)]
    # Recheck the original import sources and installed native-library identities.
    inputs.extend(row["launcher"] for row in original["runtime_before"]["files"])
    inputs.extend(r for r in original["checked_records"] if "/site-packages/" in r["path"])
    helper.check(inputs)
    output.mkdir(parents=True)
    report = dict(status="import_controls_running", checked_records=inputs, arms=[],
        accuracy_admitted=False, publication_ready=False,
        limits=dict(cpu_seconds=120, address_space_bytes=4 * 1024**3, wall_seconds=120),
        limitations=["Import and parser controls omit numerical checkpoint loading and native refinement.",
            "AST parser helper omits original preceding allocation history; no memory-safety or causal claim.",
            "No optimizer, scoring, scientific admission or comparative timing; no retries."])
    try:
        for mode in ("site_only", "frozen_imports"):
            env = {k: v for k, v in os.environ.items() if k not in ("PYTHONHOME", "LD_PRELOAD", "LD_LIBRARY_PATH")}
            overrides = dict(PYTHONPATH=str(launcher), PYTHONMALLOC="debug", PYTHONHASHSEED="0",
                PYTHONNOUSERSITE="1", PYTHONFAULTHANDLER="1", OMP_NUM_THREADS="1",
                OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1")
            env.update(overrides)
            command = [python, "-B", str(Path(__file__).resolve()), "--root", str(root), "--worker", mode]
            arm = dict(mode=mode, command=command, cwd=str(launcher), attempts=1, environment_overrides=overrides,
                       removed_environment_keys=["PYTHONHOME", "LD_PRELOAD", "LD_LIBRARY_PATH"], status="running")
            report["arms"].append(arm)
            started = time.monotonic()
            try:
                done = subprocess.run(command, cwd=launcher, env=env, capture_output=True, timeout=120)
                stdout, stderr = done.stdout, done.stderr
                arm.update(returncode=done.returncode, status="completed" if done.returncode == 0 else "failed")
            except subprocess.TimeoutExpired as error:
                stdout, stderr = error.stdout or b"", error.stderr or b""
                arm.update(returncode=None, status="timed_out")
            arm["wall_seconds_descriptive_only"] = time.monotonic() - started
            for label, content in (("stdout", stdout), ("stderr", stderr)):
                path = output / f"{mode}.{label}"
                path.write_bytes(content)
                arm[label] = helper.record(path)
            if arm["status"] == "completed":
                arm["result"] = json.loads(stdout)
            helper.check(inputs)
        report["status"] = "import_controls_observed_not_admitted"
    except Exception as error:
        report.update(status="import_controls_failed", error_type=type(error).__name__, error=str(error))
        raise
    finally:
        with (output / "report.json").open("x") as stream:
            json.dump(report, stream, indent=2, sort_keys=True)
            stream.write("\n")
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--worker", choices=("site_only", "frozen_imports"))
    args = parser.parse_args()
    if args.worker:
        worker(args.root.resolve(), args.worker)
    else:
        if args.output is None:
            parser.error("--output is required")
        run(args.root.resolve(), args.output.absolute())
