"""Bounded parser-only controls for the post-refinement high-CPM crash."""

import argparse
import ast
import gc
import hashlib
import json
import os
from pathlib import Path
import resource
import subprocess
import sys
import time

PINS = {
    "benchmarks/work/publication_qfo_replay_native_v1/benchmark_tools/audit_historical_profile_ablation.py":
        "ffdafda2b55c9580eccc7497881f0e160450af6c9029e3540009bebba2081b31",
    "benchmarks/results/qfo_cpm_checkpoint_recovery_v1/payload/gene_names.txt":
        "ad6ee208adc2c3e90bab86d68e40c32857c0f93a89989bd1a70d76385ef0edd7",
    "benchmarks/results/qfo_cpm_checkpoint_recovery_v1/orthogroups_profiles_refined.txt":
        "f4c6f1973bc9636828baf1e6d9be3f416a18fc8ada5302fb502c082495fc1811",
}
STATUS = "benchmarks/results/qfo_cpm_refinement_memcheck_diagnostic_v1/status.json"
STATUS_SHA = "33d543f8d47b1df8175030aa209a1ea19406545f37401ba1b7c3f467ebea3d3a"
PROTOCOL = "benchmark_tools/results/QFO_CPM_PARSER_CONTROL_PROTOCOL_20260926.md"


def record(path):
    path = Path(path).resolve()
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return dict(path=str(path), bytes=path.stat().st_size, sha256=digest.hexdigest())


def check(records):
    if any(record(item["path"]) != item for item in records):
        raise ValueError("Changed parser control input")


def parse_only(source, names_path, partition):
    tree = ast.parse(source.read_text(), filename=str(source))
    functions = [node for node in tree.body if isinstance(node, ast.FunctionDef) and node.name == "read_partition"]
    if len(functions) != 1 or functions[0].decorator_list:
        raise ValueError("Expected one undecorated frozen parser")
    namespace = {}
    exec(compile(ast.Module(body=functions, type_ignores=[]), str(source), "exec"), namespace)
    names = names_path.read_text().splitlines()
    universe = set(names)
    if len(names) != len(universe):
        raise ValueError("Duplicate parser gene universe")
    before = dict(enabled=gc.isenabled(), thresholds=gc.get_threshold(), stats=gc.get_stats())
    groups = namespace["read_partition"](partition, universe)
    return dict(genes=len(names), groups=len(groups), memberships=sum(map(len, groups)),
                gc_before=before, gc_after=dict(enabled=gc.isenabled(), thresholds=gc.get_threshold(),
                                               stats=gc.get_stats()),
                imported_modules=sorted(sys.modules), python_version=sys.version)


def pinned_inputs(root):
    records = [record(root / name) for name in PINS]
    if [r["sha256"] for r in records] != list(PINS.values()):
        raise ValueError("Changed frozen parser/data")
    return records


def worker(root):
    resource.setrlimit(resource.RLIMIT_AS, (4 * 1024**3, 4 * 1024**3))
    resource.setrlimit(resource.RLIMIT_CPU, (120, 120))
    resource.setrlimit(resource.RLIMIT_CORE, (0, 0))
    inputs = pinned_inputs(root)
    result = parse_only(*(Path(r["path"]) for r in inputs))
    check(inputs)
    if (result["genes"], result["groups"], result["memberships"]) != (984137, 390845, 984137):
        raise ValueError("Wrong parser control coverage")
    if any(name.split(".")[0] in {"numpy", "Bio", "igraph", "leidenalg", "orthohmm"}
           for name in result["imported_modules"]):
        raise ValueError("Scientific module imported into parser control")
    result.update(cpu_affinity=sorted(os.sched_getaffinity(0)), inputs=inputs)
    print(json.dumps(result, sort_keys=True))


def run(root, output):
    if output.exists() or output.is_symlink():
        raise FileExistsError(output)
    inputs = pinned_inputs(root)
    status_record = record(root / STATUS)
    if status_record["sha256"] != STATUS_SHA:
        raise ValueError("Changed retained diagnostic")
    status = json.loads(Path(status_record["path"]).read_text())
    python = Path(status["scientific_child_command"][0])
    wanted = {str(python.resolve()), "/usr/lib/x86_64-linux-gnu/libc.so.6"}
    runtime = {r["path"]: r for r in status["checked_records"] if r["path"] in wanted}
    if set(runtime) != wanted:
        raise ValueError("Missing retained interpreter/libc")
    inputs.extend([status_record, *runtime.values(), record(__file__), record(root / PROTOCOL)])
    check(inputs)
    output.mkdir(parents=True)
    report = dict(status="parser_controls_running", checked_records=inputs, arms=[],
                  accuracy_admitted=False, publication_ready=False,
                  limits=dict(cpu_seconds=120, address_space_bytes=4 * 1024**3, wall_seconds=120),
                  limitations=["Parser-only allocation history, not refinement or scientific-runtime replay.",
                      "AST extracts the exact frozen function but omits module imports and preceding native activity.",
                      "Two preplanned allocator arms, one execution each, no retry or scientific admission.",
                      "Shared-host observations are not controlled comparative timing or proof of memory safety."])
    try:
        for allocator in ("default", "debug"):
            env = {k: v for k, v in os.environ.items() if k not in (
                "PYTHONHOME", "PYTHONPATH", "PYTHONMALLOC", "LD_PRELOAD", "LD_LIBRARY_PATH")}
            env.update(PYTHONHASHSEED="0", PYTHONNOUSERSITE="1", PYTHONFAULTHANDLER="1")
            if allocator == "debug":
                env["PYTHONMALLOC"] = "debug"
            command = [str(python), "-B", "-S", str(Path(__file__).resolve()), "--root", str(root), "--worker"]
            arm = dict(allocator=allocator, command=command, cwd=str(output), attempts=1,
                       environment_overrides={k: env[k] for k in (
                           "PYTHONHASHSEED", "PYTHONNOUSERSITE", "PYTHONFAULTHANDLER")},
                       removed_environment_keys=["PYTHONHOME", "PYTHONPATH", "LD_PRELOAD", "LD_LIBRARY_PATH"],
                       status="running")
            arm["environment_overrides"]["PYTHONMALLOC"] = env.get("PYTHONMALLOC")
            report["arms"].append(arm)
            started = time.monotonic()
            try:
                done = subprocess.run(command, cwd=output, env=env, capture_output=True, timeout=120)
                stdout, stderr = done.stdout, done.stderr
                arm.update(returncode=done.returncode, status="completed" if done.returncode == 0 else "failed")
            except subprocess.TimeoutExpired as error:
                stdout, stderr = error.stdout or b"", error.stderr or b""
                arm.update(returncode=None, status="timed_out")
            arm["wall_seconds_descriptive_only"] = time.monotonic() - started
            for label, content in (("stdout", stdout), ("stderr", stderr)):
                path = output / f"{allocator}.{label}"
                path.write_bytes(content)
                arm[label] = record(path)
            if arm["status"] == "completed":
                arm["result"] = json.loads(stdout)
            check(inputs)
        report["status"] = "parser_controls_observed_not_admitted"
    except Exception as error:
        report.update(status="parser_controls_failed", error_type=type(error).__name__, error=str(error))
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
    parser.add_argument("--worker", action="store_true")
    args = parser.parse_args()
    if args.worker:
        worker(args.root.resolve())
    else:
        if args.output is None:
            parser.error("--output is required")
        run(args.root.resolve(), args.output.absolute())
