"""Freeze corrected-input replay commands after native HMM evidence admission."""

import argparse
import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.prepare_ob_candidate_neighborhood import record, check
from benchmark_tools.run_simulation_methods import read_frozen
from benchmark_tools.admit_qfo_corrected_high_sensitivity import PLAN_SHA, METHOD
from benchmark_tools.verify_qfo_replay_launcher import verify


def validate_admission(admission, primary, inputs):
    if (admission.get("status") != "corrected_high_sensitivity_native_evidence_admitted"
            or admission.get("accuracy_evaluated") is not False):
        raise ValueError("Corrected native HMM evidence not admitted")
    scheduler = admission["scheduler"]
    if (scheduler["State"] != "COMPLETED" or scheduler["ExitCode"] != "0:0"
            or scheduler["NodeList"] != "bizon" or scheduler["AllocCPUS"] != "32"):
        raise ValueError("Unsuccessful or incompatible inference accounting")
    content = admission["content"]
    if (content["status"] != "high_sensitivity_output_content_verified"
            or type(content["genes"]) is not int or content["genes"] != 984137
            or content["accuracy_evaluated"] is not False):
        raise ValueError("Wrong corrected gene universe or content status")
    if len(inputs) != 78 or len({r["path"] for r in inputs}) != 78:
        raise ValueError("Require 78 distinct corrected FASTAs")
    if set(content["species_ownership"]) != {r["path"] for r in inputs}:
        raise ValueError("Input inventory differs from admitted species ownership")
    if len(set(content["species_ownership"].values())) != 78:
        raise ValueError("Admitted species codes are not distinct")
    if {Path(r["path"]).parent for r in inputs} != {Path(primary["input_directory"])}:
        raise ValueError("FASTA paths outside corrected directory")
    checkpoint = Path(primary["methods"][METHOD]["output"]) / "orthohmm_working_res/high_sensitivity_checkpoint"
    manifest = admission["checkpoint_manifest"]
    if manifest["path"] != str(checkpoint / "manifest.json") or manifest != content["numeric_checkpoint"]["manifest"]:
        raise ValueError("Checkpoint manifest differs from admitted corrected checkpoint")
    checked = admission["checked_records"]
    if any(item not in checked for item in inputs + [manifest]):
        raise ValueError("Inputs/checkpoint absent from admission's checked records")
    return checkpoint, manifest["sha256"]


def command_for(python, launcher, output, fasta, checkpoint, checkpoint_sha):
    return [str(python), str(launcher / "benchmark_tools/replay_high_sensitivity.py"),
            "--accuracy-checkpoint", str(checkpoint), "--checkpoint-sha256", checkpoint_sha,
            "--fasta-directory", str(fasta), "--output-directory", str(output / "replay"),
            "--json", str(output / "replay.json"), "--cpu", "32", "--matrix", "BLOSUM62",
            "--cpm-resolution", "0.1", "--leiden-seed", "4", "--profile-iterations", "1",
            "--profile-min-species", "1"]


def prepare(root, admission_path, admission_sha, output, destination):
    if output.exists() or destination.exists():
        raise FileExistsError("Require fresh replay output and manifest")
    primary_path = root / "benchmark_tools/results/qfo_corrected_primary_commands_20260918.json"
    primary = read_frozen(primary_path, PLAN_SHA)
    admission = read_frozen(admission_path, admission_sha)
    stage_path = Path(primary["input_directory"]) / "staging_manifest.json"
    stage_record = record(stage_path)
    if stage_record not in primary["inputs"]:
        raise ValueError("Stage manifest is not pinned by corrected primary plan")
    inputs = json.loads(stage_path.read_text())["input_fastas"]
    checkpoint, checkpoint_sha = validate_admission(admission, primary, inputs)
    checked = [record(admission_path), record(primary_path), stage_record,
               *primary["inputs"], *admission["checked_records"], *admission["content"]["checked_records"]]
    for item in checked:
        check(item)
    core = root / "benchmarks/work/publication_method_native_v2"
    launcher = root / "benchmarks/work/publication_qfo_replay_native_v1"
    runtime = verify(core, launcher, root / "benchmark_tools/results/publication_native_runtime_20260916.json")
    python = Path(primary["methods"][METHOD]["native_argv"][0])
    checked.append(record(python))
    command = command_for(python, launcher, output, Path(primary["input_directory"]), checkpoint, checkpoint_sha)
    report = {"status": "corrected_replay_command_frozen_unrun", "source": record(__file__),
        "admission": record(admission_path), "primary_plan": record(primary_path),
        "checked_records": checked, "input_fastas": inputs, "checkpoint_manifest": admission["checkpoint_manifest"],
        "native_command": command, "cwd": str(launcher), "output_root": str(output), "runtime": runtime,
        "environment_overrides": {"PYTHONPATH": str(launcher), "PYTHONHASHSEED": "0",
            "OMP_NUM_THREADS": "1", "OPENBLAS_NUM_THREADS": "1", "MKL_NUM_THREADS": "1"},
        "expected_stages": ["multipass", "multipass_refined", "strict_profiles", "strict_profiles_refined"],
        "execution_authorized": False, "accuracy_evaluated": False,
        "remaining_gates": [
            "Freeze corrected checked-clustering worker and driver against this plan; do not execute native_command alone.",
            "Require all four checked native clustering boundaries and complete corrected gene coverage.",
            "Compare final replay partition to fresh native high-sensitivity output; report nonequivalence, never transfer scores.",
            "Admit replay outputs before building the four candidate arms and eight P/C/R cells.",
            "Incremental replay on the shared host is not dedicated end-to-end timing."]}
    for item in checked:
        check(item)
    with destination.open("x") as stream:
        json.dump(report, stream, indent=2, sort_keys=True)
        stream.write("\n")
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("root", "admission", "output-root", "manifest"):
        parser.add_argument("--" + name, type=Path, required=True)
    parser.add_argument("--admission-sha256", required=True)
    args = parser.parse_args()
    prepare(args.root.resolve(), args.admission.resolve(), args.admission_sha256,
            args.output_root.resolve(), args.manifest.resolve())
