"""Prepare paired native collector controls without preparing or launching work."""

import argparse
from copy import deepcopy
import hashlib
import json
from pathlib import Path

from benchmark_tools.prepare_scaling_inputs import METHODS, SIZES, planned_runs
from benchmark_tools.prepare_threadripper_run import paths
from benchmark_tools.probe_dgx_step_separation import save
from benchmark_tools.run_simulation_methods import read_frozen
from benchmark_tools.snapshot_orthohmm_input_order import record

SEED = "threadripper-native-observer-20260930-v1"
ARMS = ("boundary", "periodic")
RESOURCES = dict(host="bizon", native_workers=32, native_affinity=list(range(32)),
                 scheduler_cpu_slots=64, memory_bytes=128 * 1024**3,
                 native_timeout_s=85800)
BUDGET = dict(every_pair_max=0.10, per_method_size_median_max=0.05,
              required_pairs_per_method_size=3)


def same(left, right):
    return json.dumps(left, sort_keys=True, allow_nan=False) == json.dumps(
        right, sort_keys=True, allow_nan=False)


def arm_order(method, size, repeat):
    if method not in METHODS or type(size) is not int or size not in SIZES:
        raise ValueError("Unknown method or size")
    if type(repeat) is not int or repeat not in range(3):
        raise ValueError("Unknown repeat")
    token = f"{SEED}\n{method}\n{size}\n{'tie' if repeat == 2 else 'base'}"
    reverse = bool(hashlib.sha256(token.encode("ascii")).digest()[0] & 1)
    if repeat == 1:
        reverse = not reverse
    return list(reversed(ARMS)) if reverse else list(ARMS)


def fresh_root(value):
    path = Path(value)
    if (not path.is_absolute() or path.resolve() != path or path.exists()
            or path.is_symlink() or ".." in path.parts):
        raise ValueError("Require fresh absolute roots without symlink traversal")
    return path


def move_path(value, old_root, new_root, old_input, new_input):
    path = Path(value)
    if path == old_input:
        return str(new_input)
    if path.is_absolute() and path.is_relative_to(old_root):
        return str(new_root / path.relative_to(old_root))
    return value


def relocate(run, output, inputs, index):
    old_root, old_input = paths(run)
    new_root = output / f"run_{index:02d}"
    new_input = inputs / new_root.name / "input"
    result = deepcopy(run)

    def move(value):
        return move_path(value, old_root, new_root, old_input, new_input)

    result["index"] = index
    result["native_argv"] = [move(v) for v in run["native_argv"]]
    config = result["configuration"]
    config["argv"] = [move(v) for v in config["argv"]]
    for key in ("output", "metrics", "copy_inputs_to"):
        if key in config:
            config[key] = move(config[key])
    result["measurement_directory"] = str(new_root / "measurement")
    result["prepared_input_directory"] = str(new_input)
    paths(result)

    # Check the full object, not merely selected CLI flags, after undoing paths.
    restored = deepcopy(result)
    restored["index"] = run["index"]
    back = lambda v: move_path(v, new_root, old_root, new_input, old_input)
    restored["native_argv"] = [back(v) for v in restored["native_argv"]]
    restored["configuration"]["argv"] = [back(v) for v in config["argv"]]
    for key in ("output", "metrics", "copy_inputs_to"):
        if key in config:
            restored["configuration"][key] = back(config[key])
    restored["measurement_directory"] = run["measurement_directory"]
    restored["prepared_input_directory"] = run["prepared_input_directory"]
    if not same(restored, run):
        raise ValueError("Relocation changed native work beyond fresh run/input paths")
    return result


def build(parent, output_root, input_root):
    identities = [{k: r[k] for k in ("index", "method", "proteomes", "repeat")}
                  for r in parent["runs"]]
    if (not same(identities, planned_runs()) or not same(parent["resources"], RESOURCES)
            or parent["status"] != "threadripper_commands_prepared_unrun"
            or parent["scientific_execution_authorized"] is not False):
        raise ValueError("Require the unrun local frozen 27-identity command plan")
    output, inputs = fresh_root(output_root), fresh_root(input_root)
    if (output.is_relative_to("/dev/shm") or not inputs.is_relative_to("/dev/shm")
            or inputs == Path("/dev/shm")):
        raise ValueError("Require persistent outputs and a private tmpfs input root")
    tasks = []
    cells = {}
    for original in parent["runs"]:
        old_root, old_input = paths(original)
        for new, old in ((output, old_root.parent), (inputs, old_input.parent.parent)):
            if new.is_relative_to(old) or old.is_relative_to(new):
                raise ValueError("Overhead roots overlap the frozen production roots")
        expected = "orthofinder_full" if original["method"] == METHODS[2] else original["method"]
        if original["native_method"] != expected:
            raise ValueError("Method and native method differ")
        size = original["proteomes"]
        dataset = original["dataset"]
        names = original["input_creation_order"]
        if (dataset["proteomes"] != size or len(dataset["inputs"]) != size
                or len(set(names)) != size or len(names) != size
                or set(names) != {Path(r["path"]).name for r in dataset["inputs"]}):
            raise ValueError("Input membership or enumeration differs")
        expected_order = sorted(names) if expected == "orthofinder_full" else names
        if original["expected_native_order"] != expected_order:
            raise ValueError("Expected native enumeration differs")
        key = (original["method"], size)
        invariant = dict(dataset=dataset, input_creation_order=names,
                         expected_native_order=original["expected_native_order"])
        if key in cells and not same(invariant, cells[key]):
            raise ValueError("Input bytes or order differ across same-cell repeats")
        cells[key] = deepcopy(invariant)
        for arm in arm_order(original["method"], size, original["repeat"]):
            index = len(tasks)
            tasks.append(dict(index=index, pair=original["index"], arm=arm,
                method=original["method"], proteomes=size, repeat=original["repeat"],
                run=relocate(original, output, inputs, index)))
    return dict(schema="threadripper_native_overhead_plan_v1",
        status="prospective_native_overhead_unrun", purpose="native_point_collector_incremental_slowdown",
        ordering_seed=SEED, resources=deepcopy(RESOURCES), runs=tasks,
        parent_identities=identities, output_root=str(output), input_root=str(inputs),
        environment_overrides=deepcopy(parent["environment_overrides"]),
        prepend_path=deepcopy(parent["prepend_path"]),
        resolved_executables=deepcopy(parent["resolved_executables"]),
        engineering_budget=deepcopy(BUDGET),
        statistic="periodic_native_wall_s / boundary_native_wall_s - 1",
        collector_arms=dict(boundary="boundary_only_implementation_and_native_validation_required",
            periodic="benchmark_tools.measure_threadripper_scaling.measure"),
        common_host_observation=dict(interval_s=30, required_in_both_arms=True),
        failure_policy="retain_all_attempts_stop_after_failure_no_selective_retry",
        execution_authorized=False, scientific_timings_admitted=False, publication_ready=False,
        limitations=["54 engineering tasks are separate from the unchanged 27 production identities.",
            "Preparation only; no input/output directories, scheduler job or passing eligibility review created.",
            "Boundary control, source/runtime freeze, environment review, native handoff and independent pair audit remain required.",
            "Identical outputs do not prove identical internal computational work.",
            "Common whole-host observation remains in both arms; its cost is not isolated by this contrast.",
            "Cache preparation is symmetric but filesystem caches are not cold or evicted.",
            "Engineering budgets are not confidence intervals, timing corrections or production admission."])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--parent", type=Path, required=True)
    parser.add_argument("--parent-sha256", required=True)
    parser.add_argument("--protocol", type=Path, required=True)
    parser.add_argument("--protocol-sha256", required=True)
    parser.add_argument("--output-root", type=Path, required=True)
    parser.add_argument("--input-root", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, required=True)
    args = parser.parse_args()
    if args.manifest.exists() or args.manifest.is_symlink():
        raise FileExistsError(args.manifest)
    protocol = record(args.protocol)
    if protocol["sha256"] != args.protocol_sha256:
        raise ValueError("Prospective protocol differs")
    parent = read_frozen(args.parent, args.parent_sha256)
    result = build(parent, args.output_root, args.input_root)
    result["sources"] = [record(args.parent), protocol, record(__file__)]
    if result["sources"][0]["sha256"] != args.parent_sha256:
        raise ValueError("Parent changed during preparation")
    if record(args.protocol) != protocol:
        raise ValueError("Protocol changed during preparation")
    result["helpers"] = [record(Path(__file__).with_name(name)) for name in (
        "prepare_scaling_inputs.py", "prepare_threadripper_run.py",
        "prepare_native_scaling_run.py", "run_simulation_methods.py",
        "snapshot_orthohmm_input_order.py", "probe_dgx_step_separation.py")]
    save(args.manifest, result)


if __name__ == "__main__":
    main()
