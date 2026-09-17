"""Prepare fresh native scaling inputs outside the measured inference boundary."""

from copy import deepcopy
from pathlib import Path
import shutil

from benchmark_tools.gnu_time_companion import command as time_command
from benchmark_tools.snapshot_orthohmm_input_order import record


def run_directory(run):
    measurement = Path(run["measurement_directory"])
    if not measurement.is_absolute() or measurement.name != "measurement":
        raise ValueError("Unexpected measurement path")
    root = measurement.parent
    if root.is_symlink() or root.resolve() != root:
        raise ValueError("Run directory must not traverse symlinks")
    for key in ("output", "metrics", "copy_inputs_to"):
        if key not in run["configuration"]:
            continue
        path = Path(run["configuration"][key])
        if not path.is_absolute() or not path.is_relative_to(root) or path == root or ".." in path.parts:
            raise ValueError("Native path escapes fresh run directory")
        if path.is_relative_to(measurement):
            raise ValueError("Native output overlaps collector output")
    return root


def check_original_inputs(run, order):
    dataset = run["dataset"]
    if order["input_directory"] != dataset["input_directory"] or order["proteomes"] != dataset["proteomes"]:
        raise ValueError("Native order describes a different dataset")
    expected = {Path(row["path"]).name: row for row in dataset["inputs"]}
    if len(expected) != len(dataset["inputs"]) or len(expected) != dataset["proteomes"]:
        raise ValueError("Ambiguous input membership")
    names = order["native_order"]
    if len(names) != len(expected) or set(names) != set(expected):
        raise ValueError("Frozen native order has different membership")
    if order["inputs_in_native_order"] != [expected[name] for name in names]:
        raise ValueError("Native order input identities differ")
    observed = [record(Path(dataset["input_directory"]) / name) for name in names]
    if observed != order["inputs_in_native_order"]:
        raise ValueError("Original input bytes or paths changed")
    return observed


def check_copies(run):
    directory = Path(run["configuration"]["copy_inputs_to"])
    expected = {Path(row["path"]).name: (row["bytes"], row["sha256"]) for row in run["dataset"]["inputs"]}
    copied = [record(path) for path in sorted(directory.iterdir()) if path.is_file()]
    actual = {Path(row["path"]).name: (row["bytes"], row["sha256"]) for row in copied}
    if len(copied) != len(expected) or actual != expected:
        raise ValueError("Copied native inputs differ from original basenames/bytes")
    return copied


def prepare(run, order):
    root = run_directory(run)
    if not root.is_dir() or any(root.iterdir()):
        raise ValueError("Caller must create a fresh empty run directory")
    originals = check_original_inputs(run, order)
    config = run["configuration"]
    method = run["native_method"]
    if method == "orthofinder_full":
        destination = Path(config["copy_inputs_to"])
        if Path(config["copy_inputs_from"]) != Path(run["dataset"]["input_directory"]):
            raise ValueError("Copy source differs from frozen dataset")
        destination.mkdir(parents=True, exist_ok=False)
        for item in originals:
            source = Path(item["path"])
            with source.open("rb") as src, (destination / source.name).open("xb") as dst:
                shutil.copyfileobj(src, dst)
        copies = check_copies(run)
        native_order = sorted(Path(item["path"]).name for item in copies)
    elif method in {"orthohmm_high_sensitivity", "orthohmm_satellite_v2"}:
        Path(config["output"]).mkdir(parents=True, exist_ok=False)
        copies = []
        native_order = order["native_order"]
    else:
        raise ValueError("Unknown native method")
    # Recheck original bytes after copying; partial preparation is never erased.
    check_original_inputs(run, order)
    measured_run = deepcopy(run)
    measured_run["gnu_time"] = {"executable": "/usr/bin/time", "output": str(root / "native.time.tsv")}
    return {"run": measured_run,
            "measured_argv": time_command(run["native_argv"], measured_run["gnu_time"]["output"]),
            "original_inputs": originals, "copied_inputs": copies,
            "expected_native_basename_order": native_order,
            "status": "fresh_native_inputs_prepared", "inference_started": False,
            "limitations": ["No inference or output admission is performed here.",
                            "Caller must invoke the actual OrthoHMM enumerator immediately before/after inference.",
                            "Copies and input hashing warm file caches and are outside the inference timer."]}
