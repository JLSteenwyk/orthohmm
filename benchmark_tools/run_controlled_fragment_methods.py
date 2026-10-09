"""Run the 30 new fragment identities sequentially; never retry an identity."""

import argparse
import json
import os
from pathlib import Path
import subprocess
import sys
import time

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from Bio import SeqIO
from benchmark_tools.assemble_simulation_results import input_universe
from benchmark_tools.controlled_fragment_observations import CONDITION, SEEDS, METHODS, INFERENCE_METHODS, fresh_methods, verify_observation
from benchmark_tools.controlled_fragment_runtime import SCHEMA, PARENT, adopt_inventory
from benchmark_tools.prepare_controlled_fragment_observations import PINS, PROTOCOL, read_input_sequences, record
from benchmark_tools.run_simulation_generation import verify_file
from benchmark_tools.run_simulation_methods import read_frozen, verify_environment, execution_environment, execute

SOURCES = ("run_controlled_fragment_methods.py", "controlled_fragment_observations.py",
           "controlled_fragment_runtime.py", "prepare_controlled_fragment_observations.py",
           "run_simulation_methods.py", "run_simulation_generation.py", "benchmark_production.py",
           "build_publication_runtime.py", "validate_profile_runtime.py", "assemble_simulation_results.py",
           "validate_simulation_outputs.py", "simulation_method_outputs.py", "simulation_conditions.py")


def observed_inputs(dataset):
    verified = dataset["verified_inputs"]
    if verified["status"] != "ready" or Path(verified["truth"]["absolute_path"]).resolve() != Path(dataset["truth"]).resolve():
        raise ValueError("Fragment truth path or readiness differs")
    verify_file(Path(dataset["truth"]), verified["truth"])
    truth = json.loads(Path(dataset["truth"]).read_text())
    owners, species = input_universe(verified["inputs"], truth)
    paths = {Path(item["absolute_path"]).resolve() for item in verified["inputs"]}
    inputs = Path(dataset["input"]).resolve()
    if paths != set(inputs.iterdir()) or any(p.parent != inputs or p.suffix != ".fasta" for p in paths):
        raise ValueError("Fragment input file inventory changed")
    observed = {sequence.id: str(sequence.seq) for path in sorted(paths) for sequence in SeqIO.parse(path, "fasta")}
    parent = dataset["parent_inputs"]
    verify_file(Path(parent["truth"]["absolute_path"]), parent["truth"])
    if (parent["truth"]["bytes"], parent["truth"]["sha256"]) != (verified["truth"]["bytes"], verified["truth"]["sha256"]):
        raise ValueError("Observation changed evolutionary truth bytes")
    directories = {Path(item["absolute_path"]).parent.resolve() for item in parent["inputs"]}
    if len(directories) != 1:
        raise ValueError("Baseline input directories are inconsistent")
    full, old_owners, old_species = read_input_sequences(parent, truth, directories.pop())
    if owners != old_owners or set(species) != set(old_species):
        raise ValueError("Observation changed species ownership")
    coordinate_record = dataset["coordinates"]
    coordinate_path = Path(coordinate_record["absolute_path"])
    verify_file(coordinate_path, coordinate_record)
    coordinates = json.loads(coordinate_path.read_text())
    checked = verify_observation(full, observed, owners, coordinates, dataset["seed"])
    if checked != dataset["verification"]:
        raise ValueError("Saved transformation verification differs from disk readback")
    return verified


def load_panel(root, path, sha, runtime_path, runtime_sha, destination):
    panel = read_frozen(path, sha)
    if (panel.get("schema") != "controlled_fragment_observations_v1" or panel["status"] != "prepared_unexecuted"
            or panel["condition"] != CONDITION or panel["inference_order"] != list(INFERENCE_METHODS)
            or panel["inference_identities"] != 30 or panel["inference_launched"] is not False
            or panel["new_accuracy_evaluated"] is not False
            or [d["seed"] for d in panel["datasets"]] != list(SEEDS)):
        raise ValueError("Wrong fragment panel, identity inventory or execution order")
    expected_pins = [record(root / "benchmark_tools/results" / name) for name in PINS]
    if panel["baseline_pins"] != expected_pins:
        raise ValueError("Fragment parent artifact identity differs")
    for name, item in zip(PINS, expected_pins):
        if item["sha256"] != PINS[name]:
            raise ValueError("Changed pinned baseline artifact")
    protocol = root / "benchmark_tools/results" / PROTOCOL
    if panel["protocol"] != record(protocol):
        raise ValueError("Fragment protocol changed")
    for item in (panel["source"], panel["transform_source"]):
        verify_file(Path(item["absolute_path"]), item)
    runtime = read_frozen(runtime_path, runtime_sha)
    parent_path = root / "benchmark_tools/results" / PARENT
    parent = read_frozen(parent_path, PINS[PARENT])
    expected = adopt_inventory(parent, runtime["environments"], record(parent_path))
    if runtime.get("schema") != SCHEMA or {k: v for k, v in runtime.items() if k not in {"adoption_source", "path_resolved_executables"}} != expected:
        raise ValueError("Prospective runtime changed outside the explicit inventory adoption")
    verify_file(Path(runtime["adoption_source"]["absolute_path"]), runtime["adoption_source"])
    verify_environment(runtime)
    env, resolved = execution_environment(runtime)
    if resolved != runtime["path_resolved_executables"]:
        raise ValueError("Prospective tool resolution changed")
    baselines = {d["seed"]: d for d in parent["datasets"] if d["condition"] == "baseline"}
    for dataset in panel["datasets"]:
        label = f"{CONDITION}_{dataset['seed']}"
        if dataset["label"] != label or dataset["condition"] != CONDITION or dataset["parent"] != f"baseline_{dataset['seed']}":
            raise ValueError("Fragment dataset identity changed")
        if dataset["methods"] != fresh_methods(baselines[dataset["seed"]], dataset["input"], destination / label):
            raise ValueError("Fragment scientific command or output binding differs")
        observed_inputs(dataset)
    return panel, runtime, env


def allocation():
    job = os.environ.get("SLURM_JOB_ID")
    cpus = int(os.environ.get("SLURM_CPUS_PER_TASK", "0"))
    memory = int(os.environ.get("SLURM_MEM_PER_NODE", "0"))
    if not job or cpus < 16 or memory < 16384:
        raise ValueError("Require this project's scheduler allocation: at least 16 CPUs and 16 GiB")
    raw = subprocess.check_output(["scontrol", "show", "job", job], text=True)
    return {"job_id": job, "cpus_per_task": cpus, "memory_mib": memory, "scontrol": raw,
            "scope": "Shared-host enforced job limits, not node isolation"}


def run(root, path, sha, runtime_path, runtime_sha, destination):
    if destination.exists():
        raise FileExistsError("Existing fragment inference panel; no automatic restart")
    panel, runtime, env = load_panel(root, path, sha, runtime_path, runtime_sha, destination)
    resources = allocation()
    sources = [record(Path(__file__).with_name(name)) for name in SOURCES]
    provenance = {"fragment_manifest": record(path), "runtime_manifest": record(runtime_path),
                  "executor_commit_at_start": subprocess.check_output(["git", "-C", str(root), "rev-parse", "HEAD"], text=True).strip(),
                  "sources": sources, "allocation": resources,
                  "timing_limitation": "Observed shared-Threadripper inference resources; contention unknown and tool-dependent; no isolated performance claim"}
    destination.mkdir(parents=True, exist_ok=False)
    status = {"schema": "controlled_fragment_panel_execution_v1", "status": "running",
              "started_epoch": time.time(), "provenance": provenance, "datasets": [],
              "accuracy_evaluated": False, "native_outputs_validated": False}
    status_path = destination / "panel_status.json"

    def save():
        temporary = destination / "panel_status.tmp"
        temporary.write_text(json.dumps(status, indent=2, sort_keys=True) + "\n")
        temporary.replace(status_path)

    save()
    try:
        for dataset in panel["datasets"]:
            row = {"label": dataset["label"], "seed": dataset["seed"], "status": "preflight"}
            status["datasets"].append(row)
            save()
            try:
                for item in sources:
                    verify_file(Path(item["absolute_path"]), item)
                read_frozen(path, sha)
                read_frozen(runtime_path, runtime_sha)
                verify_environment(runtime)
                verified = observed_inputs(dataset)
                evidence = destination / dataset["label"] / "execution"
                result = execute(dataset, INFERENCE_METHODS, env, evidence, verified, provenance, runtime)
                row.update(status="finished_pending_native_validation", execution=record(evidence / "status.json"),
                           failed_methods=result["failed_methods"])
            except Exception as error:
                row.update(status="failed", failure_stage="preflight_or_execution", error_type=type(error).__name__, reason=str(error))
            save()
        for item in sources:
            verify_file(Path(item["absolute_path"]), item)
        status["status"] = "finished_pending_native_validation"
    except BaseException as error:
        status.update(status="interrupted_or_integrity_failure", error_type=type(error).__name__, reason=str(error))
        raise
    finally:
        status["finished_epoch"] = time.time()
        save()
    return status


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--manifest-sha256", required=True)
    parser.add_argument("--runtime", type=Path, required=True)
    parser.add_argument("--runtime-sha256", required=True)
    parser.add_argument("--destination", type=Path, required=True)
    parser.add_argument("--check-only", action="store_true")
    args = parser.parse_args()
    values = (args.root.resolve(), args.manifest.resolve(), args.manifest_sha256,
              args.runtime.resolve(), args.runtime_sha256, args.destination.resolve())
    if args.check_only:
        panel, _, _ = load_panel(*values)
        print(json.dumps({"status": "fragment_launch_preflight_verified", "datasets": len(panel["datasets"]), "new_inference": False}))
    else:
        result = run(*values)
        print(json.dumps({"status": result["status"], "datasets": len(result["datasets"]), "accuracy_evaluated": False}))
