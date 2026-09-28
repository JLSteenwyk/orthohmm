"""Isolate the frozen partition reader; never admit refinement or accuracy."""

import argparse
import ast
import gc
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import time

FAILURE_SHA = "6261e52aa55b2d4f114d63c6387b79e818aa56199a1a2e6f90c2a883aef0cd8f"
SOURCE_SHA = "ffdafda2b55c9580eccc7497881f0e160450af6c9029e3540009bebba2081b31"


def record(path):
    path = Path(path)
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return dict(path=str(path), bytes=path.stat().st_size, sha256=digest.hexdigest())


def check(records):
    for item in records:
        if record(item["path"]) != item:
            raise ValueError("Changed diagnostic input: " + item["path"])


def reader(source):
    tree = ast.parse(source.read_text(), filename=str(source))
    functions = [node for node in tree.body
                 if isinstance(node, ast.FunctionDef) and node.name == "read_partition"]
    if len(functions) != 1 or functions[0].decorator_list:
        raise ValueError("Expected one undecorated frozen partition reader")
    # Compile the original function body and line locations, without module imports.
    namespace = {}
    exec(compile(ast.Module(body=functions, type_ignores=[]), str(source), "exec"), namespace)
    return namespace["read_partition"]


def worker(config):
    if not sys.flags.no_site or not gc.isenabled():
        raise ValueError("Require disabled site initialization and enabled GC")
    check(config["checked_records"])
    read_partition = reader(Path(config["source"]))
    names = Path(config["names"]).read_text().splitlines()
    universe = set(names)
    if len(names) != config["genes"] or len(universe) != len(names):
        raise ValueError("Wrong gene universe")
    before = gc.get_stats()
    read_partition(Path(config["initial"]), universe)
    groups = read_partition(Path(config["refined"]), universe)
    after = gc.get_stats()
    if len(groups) != config["groups"]:
        raise ValueError("Wrong refined group count")
    scientific = sorted(name for name in sys.modules if name.split(".")[0] in
                        {"Bio", "numpy", "scipy", "igraph", "leidenalg", "orthohmm", "numba", "pyhmmer"})
    if scientific:
        raise ValueError("Unexpected scientific imports")
    check(config["checked_records"])
    return dict(status="isolated_readback_completed_not_admitted", genes=len(names),
                groups=len(groups), gc_before=before, gc_after=after,
                gc_enabled=gc.isenabled(), gc_threshold=gc.get_threshold(),
                scientific_modules=scientific, loaded_modules=sorted(sys.modules),
                accuracy_admitted=False)


def run(root, output):
    failure_path = root / "benchmark_tools/results/qfo_cpm_allocator_failure_22158.json"
    if record(failure_path)["sha256"] != FAILURE_SHA:
        raise ValueError("Changed retained failure audit")
    failure = json.loads(failure_path.read_text())
    status_record = failure["checked_records"][0]
    check([status_record])
    original = json.loads(Path(status_record["path"]).read_text())
    python = original["command"][0]
    executable = next(item for item in original["checked_records"]
                      if item["path"] == str(Path(python).resolve()))
    launcher = root / "benchmarks/work/publication_qfo_replay_native_v1"
    source = launcher / "benchmark_tools/audit_historical_profile_ablation.py"
    if record(source)["sha256"] != SOURCE_SHA:
        raise ValueError("Changed frozen reader")
    native = root / "benchmarks/results/qfo_cpm_checkpoint_recovery_v1"
    names = native / "payload/gene_names.txt"
    initial = native / "orthogroups_profiles.txt"
    refined = Path(failure["checked_records"][1]["path"])
    inputs = [next(item for item in original["checked_records"] if item["path"] == str(path))
              for path in (names, initial)]
    records = [record(__file__), record(failure_path), status_record, executable,
               record(source), *inputs, failure["checked_records"][1]]
    check(records)
    output.mkdir(parents=True, exist_ok=False)
    config = dict(checked_records=records, source=str(source), names=str(names),
                  initial=str(initial), refined=str(refined), **failure["coverage"])
    config_path = output / "config.json"
    config_path.write_text(json.dumps(config, indent=2, sort_keys=True) + "\n")
    command = [python, "-B", "-S", str(Path(__file__).resolve()), "--worker", str(config_path)]
    env = os.environ.copy()
    for key in ("PYTHONHOME", "LD_PRELOAD", "LD_LIBRARY_PATH"):
        env.pop(key, None)
    overrides = dict(PYTHONPATH=str(launcher), PYTHONHASHSEED="0", PYTHONNOUSERSITE="1",
                     PYTHONMALLOC="debug", PYTHONFAULTHANDLER="1",
                     OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1")
    env.update(overrides)
    report = dict(status="isolated_readback_running", command=command, cwd=str(launcher),
                  environment_overrides=overrides, checked_records=records,
                  config=record(config_path), attempts=1, timeout_seconds=300,
                  job_id=os.environ.get("SLURM_JOB_ID"), accuracy_admitted=False)
    start = time.monotonic()
    try:
        with (output / "child.json").open("xb") as out, (output / "child.log").open("xb") as err:
            child = subprocess.run(command, cwd=launcher, env=env, timeout=300,
                                   stdout=out, stderr=err)
        report["returncode"] = child.returncode
        if child.returncode:
            raise RuntimeError("Isolated reader failed; no retry")
        report["child"] = json.loads((output / "child.json").read_text())
        check([*records, report["config"]])
        report["status"] = "isolated_readback_observed_not_admitted"
    except BaseException as error:
        report.update(status="isolated_readback_failed", error_type=type(error).__name__, error=str(error))
        raise
    finally:
        report["wall_seconds_descriptive_only"] = time.monotonic() - start
        report["artifacts"] = [record(output / name) for name in ("child.json", "child.log")
                               if (output / name).is_file()]
        with (output / "report.json").open("x") as stream:
            json.dump(report, stream, indent=2, sort_keys=True)
            stream.write("\n")
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--worker", type=Path)
    args = parser.parse_args()
    if args.worker:
        print(json.dumps(worker(json.loads(args.worker.read_text())), sort_keys=True))
    elif args.root and args.output:
        run(args.root.resolve(), args.output.absolute())
    else:
        parser.error("Supply --worker or both --root and --output")
