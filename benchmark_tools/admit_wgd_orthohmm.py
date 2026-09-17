"""Admit completed OrthoHMM application outputs without computing endpoints."""

import argparse
import json
from pathlib import Path
import shlex
import subprocess
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.run_wgd_application import DIRECTORIES, pinned
from benchmark_tools.snapshot_orthohmm_input_order import record
from benchmark_tools.validate_scaling_outputs import input_universe, validate_orthohmm


def check_receipt(receipt, selected, root, spec, accounting, task):
    if (receipt["status"] != "native_exited_zero" or receipt["native"]["exit_code"] != 0
            or receipt["native"]["timed_out"] or receipt["method"] != selected
            or receipt["hostname"] != "bizon" or receipt["environment"] != spec["environment"]):
        raise ValueError("Execution receipt does not establish frozen native completion")
    if receipt["cwd"] != str(root) or receipt["executed_argv"] != [
            "/usr/bin/time", "-v", "-o", str(root / "native.time.log"), *selected["argv"]]:
        raise ValueError("Changed executed command or working directory")
    expected = [{"records": len(pinned(row)["records"]), "scientific_execution_authorized": False,
                 "status": "runtime_tree_identity_matches"} for row in spec["runtime_manifests"]]
    if receipt["runtime_before"] != expected or receipt["runtime_after"] != expected:
        raise ValueError("Missing or changed before/after runtime identities")
    matching = [line.split("|") for line in accounting.splitlines()
                if len(line.split("|")) == 4 and line.split("|")[1] == task]
    if matching != [[str(receipt["job_id"]), task, "COMPLETED", "0:0"]]:
        raise ValueError("Scheduler does not establish this task completed successfully")


def expected_harness(selected):
    argv = selected["argv"]
    command = [argv[0], "-m", "orthohmm", argv[2], "-o", argv[3], "-c", "32",
               "--threads_per_worker", "8", "-x", "BLOSUM62", "-e", "0.0001",
               "--clustering", "leiden", "--cpm_resolution", "0.1", "--refinement_profile", "default",
               "--accuracy_profile", "high_sensitivity", "--metrics_json", argv[4], "--stop", "infer"]
    if selected["method"] == "orthohmm_satellite_v2":
        command += ["--phylogeny", "reconcile", "--species_tree_mode", "infer", "--aligner", "mafft",
                    "--tree_builder", "FastTree", "--phylogeny_candidates", "satellite_v2",
                    "--phylogeny_root_rule", "species_overlap", "--phylogeny_pair_rule", "positive_paralogy",
                    "--species_tree_rooting", "min_variance"]
    elif selected["method"] != "orthohmm_high_sensitivity":
        raise ValueError("Not an OrthoHMM application method")
    return command


def admit(repo, index, accounting):
    if index not in (0, 1):
        raise ValueError("This validator covers only OrthoHMM application tasks")
    spec_path = repo / "benchmark_tools/results/biological_wgd_execution_20260917.json"
    spec = pinned({"path": str(spec_path), "sha256": "c43704020c56ead316678461f5be3e8d4efc43a56ff5ced4f0e3cfd3c189025c"})
    plan = pinned(spec["command_plan"])
    selected = plan["runs"][index]
    root = Path(plan["output_root"]) / DIRECTORIES[index]
    receipt_path = root / "execution.json"
    receipt = json.loads(receipt_path.read_text())
    check_receipt(receipt, selected, root, spec, accounting, f"21661_{index}")
    lines = (root / "native.time.log").read_text().splitlines()
    timed = [shlex.split(line.strip().removeprefix('Command being timed: "').removesuffix('"'))
             for line in lines if line.strip().startswith("Command being timed:")]
    if timed != [selected["argv"]] or [line.strip() for line in lines if "Exit status:" in line] != ["Exit status: 0"]:
        raise ValueError("GNU-time does not confirm exact native harness command and success")
    inputs = pinned(plan["inputs"])
    owners, species = input_universe({**inputs, "proteomes": 4})
    command = expected_harness(selected)
    config = {"metrics": selected["argv"][4], "output": selected["argv"][3]}
    metrics = json.loads(Path(config["metrics"]).read_text())
    harness = metrics["harness"]
    if (harness["command"] != command or harness["exit_code"] != 0
            or harness["git_commit"] != plan["core_commit"] or harness["git_dirty"]):
        raise ValueError("Changed frozen harness identity")
    expected_inputs = [{**r, "path": Path(r["path"]).name} for r in inputs["inputs"]]
    if harness["input_manifest"] != expected_inputs:
        raise ValueError("Harness input manifest differs")
    artifacts = []
    for key, base in (("source_manifest", Path(plan["core_root"])), ("output_manifest", Path(config["output"]))):
        for row in harness[key]:
            relative = Path(row["path"])
            if relative.is_absolute() or ".." in relative.parts:
                raise ValueError("Harness manifest path escapes its root")
            actual = record(base / relative)
            if actual["sha256"] != row["sha256"] or actual["bytes"] != row["bytes"]:
                raise ValueError("Harness source/output checksum differs")
            artifacts.append(actual)
    native = validate_orthohmm({"native_method": selected["method"], "native_argv": command,
                               "cwd": plan["core_root"], "configuration": config}, owners, species)
    return {"status": "orthohmm_application_native_outputs_admitted", "accuracy_evaluated": False,
            "method": selected["method"], "validator": record(__file__), "spec": record(spec_path),
            "receipt": record(receipt_path), "accounting": accounting, "native": native,
            "verified_artifacts": artifacts, "native_command": metrics["command"],
            "interpretation": "Native completion, hashes, counts and exact partitions verified; biological endpoints not evaluated."}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", type=Path, required=True)
    parser.add_argument("--index", type=int, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    accounting = subprocess.check_output(["sacct", "-j", "21661", "--noheader", "--parsable2",
                                         "--format=JobIDRaw,JobID,State,ExitCode", "-X"], text=True)
    result = admit(args.repo.resolve(), args.index, accounting)
    with args.output.open("x") as handle:
        json.dump(result, handle, indent=2, sort_keys=True)
        handle.write("\n")
