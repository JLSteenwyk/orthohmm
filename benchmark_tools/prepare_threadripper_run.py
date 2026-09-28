"""Fresh per-run tmpfs inputs and persistent outputs, outside inference timing."""

from copy import deepcopy
from pathlib import Path
import shutil
import time

from benchmark_tools.prepare_native_scaling_run import run_directory, check_original_inputs
from benchmark_tools.snapshot_orthohmm_input_order import record, snapshot
from benchmark_tools.probe_dgx_step_separation import save


def paths(run):
    adapted = deepcopy(run)
    adapted["configuration"].pop("copy_inputs_to", None)
    root = run_directory(adapted)
    if root.is_relative_to("/dev/shm"):
        raise ValueError("Require persistent output directory")
    target = Path(run["prepared_input_directory"])
    if (not target.is_absolute() or not target.is_relative_to("/dev/shm")
            or target.name != "input" or target.parent.name != root.name
            or target.resolve() != target or target.exists() or target.is_symlink()):
        raise ValueError("Require fresh private tmpfs input path matching the run identity")
    config = run["configuration"]
    if config["copy_inputs_to"] != str(target) or config["copy_inputs_from"] != run["dataset"]["input_directory"]:
        raise ValueError("Input-copy configuration differs")
    method = run["native_method"]
    argv = run["native_argv"]
    if method == "orthofinder_full":
        if (argv.count("-f") != 1 or argv[argv.index("-f")+1] != str(target)
                or argv.count("-o") != 1 or argv[argv.index("-o")+1] != config["output"] or "-op" in argv):
            raise ValueError("Unexpected OrthoFinder input/output command")
    elif method in ("orthohmm_high_sensitivity", "orthohmm_satellite_v2"):
        if (argv[1:4] != ["-m", "orthohmm", str(target)] or argv.count("-o") != 1
                or argv[argv.index("-o")+1] != config["output"]):
            raise ValueError("Unexpected OrthoHMM input command")
    else:
        raise ValueError("Unknown method")
    return root, target


def check_prepared(run, baseline):
    target = Path(run["prepared_input_directory"])
    expected = {Path(r["path"]).name: r for r in run["dataset"]["inputs"]}
    entries = list(target.iterdir())
    if (len(entries) != len(expected) or {p.name for p in entries} != set(expected)
            or any(p.is_symlink() or not p.is_file() for p in entries)):
        raise ValueError("Prepared input inventory changed")
    rows = [record(p) for p in entries]
    if any(any(row[k] != expected[Path(row["path"]).name][k] for k in ("bytes", "sha256")) for row in rows):
        raise ValueError("Prepared input bytes changed")
    runtime = {"records": [dict(path=r["absolute_path"], bytes=r["bytes"], sha256=r["sha256"])
                           for r in baseline["core_sources"]]}
    observed = snapshot(baseline["core_root"], {"datasets": [dict(proteomes=run["proteomes"],
        input_directory=str(target), inputs=rows)]}, runtime)
    if observed["datasets"][0]["native_order"] != run["input_creation_order"]:
        raise ValueError("Actual frozen enumeration differs from retained order")
    names = run["input_creation_order"]
    expected_order = sorted(names) if run["native_method"] == "orthofinder_full" else names
    if expected_order != run["expected_native_order"]:
        raise ValueError("Expected native mapping differs")
    return observed


def prepare(run, baseline):
    started = time.monotonic_ns()
    root, target = paths(run)
    if not root.is_dir() or any(root.iterdir()):
        raise ValueError("Caller must create a fresh empty persistent run directory")
    expected = {Path(row["path"]).name: row for row in run["dataset"]["inputs"]}
    names = run["input_creation_order"]
    if len(names) != run["proteomes"] or len(set(names)) != len(names) or set(names) != set(expected):
        raise ValueError("Wrong or duplicate input names")
    order = dict(proteomes=run["proteomes"], input_directory=run["dataset"]["input_directory"],
                 native_order=names, inputs_in_native_order=[expected[name] for name in names])
    try:
        originals = check_original_inputs(run, order)
        target.mkdir(parents=True, exist_ok=False)
        for row in originals:
            source = Path(row["path"])
            with source.open("rb") as src, (target / source.name).open("xb") as dst:
                shutil.copyfileobj(src, dst)
        observed = check_prepared(run, baseline)
        check_original_inputs(run, order)
        if run["native_method"] != "orthofinder_full":
            Path(run["configuration"]["output"]).mkdir(parents=True, exist_ok=False)
        result = dict(status="threadripper_inputs_prepared", input_snapshot=observed,
            original_inputs=originals, input_bytes=sum(r["bytes"] for r in originals),
            expected_native_basename_order=run["expected_native_order"],
            started_ns=started, finished_ns=time.monotonic_ns(),
            inference_started=False, scientific_execution_authorized=False,
            limitations=["Must run inside the allocation's preparation context for intended memory charging.",
                         "Caller must recheck before/after inference and validate actual OrthoFinder SpeciesIDs.",
                         "No runtime, allocation, quiet-host or native-output admission performed here."])
        save(root / "preparation.json", result)
        return result
    except Exception as error:
        save(root / "preparation_failed.json", dict(error_type=type(error).__name__, error=str(error),
             started_ns=started, finished_ns=time.monotonic_ns(), inference_started=False))
        raise
