"""Independently admit corrected QfO native reconciliation outputs before scoring."""

import argparse
import json
from pathlib import Path
import subprocess
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from Bio import SeqIO
from benchmark_tools.admit_qfo_factorial_cell import check_pairs
from benchmark_tools.prepare_ob_candidate_neighborhood import record, check
from benchmark_tools.run_qfo_corrected_factorial_cell import verify_admission
from benchmark_tools.run_qfo_factorial_cell import native_command, ENVIRONMENT_SHA
from benchmark_tools.run_simulation_methods import read_frozen, verify_environment, execution_environment
from benchmark_tools.validate_factorial_native import validate_native_cell
from benchmark_tools.validate_simulation_outputs import verify_process
from benchmark_tools.verify_qfo_replay_launcher import LAUNCHER_COMMIT
from benchmark_tools.verify_ygob_validation import require_completed_job

EXECUTOR_COMMIT = "2bf70fb27cc63edc7c49d38d3c7f7d09ddccca9b"


def require_success(scheduler, status, postflight, index):
    if (scheduler["State"], scheduler["ExitCode"], scheduler["NodeList"], scheduler["AllocCPUS"]) != (
            "COMPLETED", "0:0", "bizon", "32"):
        raise ValueError("Require successful terminal 32-CPU workstation task")
    provenance = status["provenance"]
    if (provenance["slurm_job_id"] != scheduler["JobIDRaw"]
            or provenance["slurm_array_task_id"] != str(index) or provenance["cell_index"] != index):
        raise ValueError("Execution scheduler/cell identity differs")
    if (status.get("status") != "finished_pending_native_validation" or status.get("failed_methods") != []
            or status.get("accuracy_evaluated") is not False or status.get("native_outputs_validated") is not False):
        raise ValueError("Incomplete or already-scored execution")
    if postflight != {"status": "corrected_inputs_runtime_sources_reverified", "failed_methods": [],
                      "native_outputs_validated": False, "accuracy_evaluated": False}:
        raise ValueError("Successful corrected postflight missing")


def validate_provenance(status, expected, manifest, cell):
    provenance = status["provenance"]
    if provenance != expected:
        raise ValueError("Execution provenance differs from reconstructed corrected command")
    inputs = {"status": "ready", "inputs": [{**r, "absolute_path": r["path"]} for r in manifest["input_fastas"]]}
    if (status["verified_inputs"] != inputs or set(status["methods"]) != {cell["label"]}
            or status["dataset"] != cell["label"]):
        raise ValueError("Execution input/method inventory mismatch")


def gene_ownership(manifest, native, candidate_path):
    taxa = {r["filename"]: r["taxon"] for r in native["input_proteomes"]}
    if len(taxa) != 78 or len(set(taxa.values())) != 78:
        raise ValueError("Require 78 distinct corrected species")
    owners = {}
    for item in manifest["input_fastas"]:
        for sequence in SeqIO.parse(item["path"], "fasta"):
            if sequence.id in owners:
                raise ValueError("Duplicate FASTA gene")
            owners[sequence.id] = taxa[Path(item["path"]).name]
    candidates = {}
    with candidate_path.open() as stream:
        for family, line in enumerate(stream):
            for gene in line.split():
                if gene in candidates:
                    raise ValueError("Duplicate candidate gene")
                candidates[gene] = family
    if len(owners) != 984137 or set(candidates) != set(owners):
        raise ValueError("Corrected candidate/FASTA universe differs")
    return owners, candidates


def admit(root, index, job, admission_path, admission_sha, admission_job):
    if type(index) is not int or not 0 <= index < 4:
        raise ValueError("Unknown corrected reconciliation cell")
    accounting = subprocess.check_output(["sacct", "-j", str(job), "--parsable2",
        "--format=JobIDRaw,State,ExitCode,Elapsed,NodeList,AllocCPUS"], text=True)
    scheduler = require_completed_job(accounting, job)
    if scheduler["NodeList"] != "bizon" or scheduler["AllocCPUS"] != "32":
        raise ValueError("Unexpected reconciliation resources")
    admission, manifest, cell, output, prepared_executor, admission_scheduler = verify_admission(
        root, admission_path, admission_sha, admission_job, index)
    executor = root / "benchmarks/work/publication_qfo_corrected_reconcile_v1"
    if subprocess.check_output(["git", "-C", str(executor), "rev-parse", "HEAD"], text=True).strip() != EXECUTOR_COMMIT:
        raise ValueError("Reconciliation executor revision changed")
    subprocess.run(["git", "-C", str(executor), "diff", "--exit-code", "HEAD", "--", "benchmark_tools", "orthohmm"], check=True)
    launcher = Path(manifest["launcher_root"])
    argv, sources = native_command(cell, launcher, prepared_executor)
    environment_path = root / "benchmark_tools/results/publication_variable_native_methods_20260916.json"
    environment = read_frozen(environment_path, ENVIRONMENT_SHA)
    verify_environment(environment)
    _, resolved = execution_environment(environment)
    evidence = output / "execution" / cell["label"]
    status_path, post_path = evidence / "status.json", evidence / "postflight.json"
    checked = [record(status_path), record(post_path)]
    status, postflight = json.loads(status_path.read_text()), json.loads(post_path.read_text())
    require_success(scheduler, status, postflight, index)
    expected = {"candidate_admission": record(admission_path), "admission_scheduler": admission_scheduler,
        "prepared_manifest": admission["prepared_manifest"], "environment": record(environment_path),
        "resolved_tools": resolved, "planned_argv": cell["argv"], "executed_argv": argv,
        "launcher_source_equivalence": sources, "executor": record(executor / "benchmark_tools/run_qfo_corrected_factorial_cell.py"),
        "execution_helper": record(executor / "benchmark_tools/run_simulation_methods.py"),
        "command_helper": record(executor / "benchmark_tools/run_qfo_factorial_cell.py"),
        "slurm_job_id": scheduler["JobIDRaw"], "slurm_array_task_id": str(index),
        "cell_index": index, "cwd": str(launcher), "scope": "corrected-release incremental reconciliation; no scoring"}
    validate_provenance(status, expected, manifest, cell)
    config = {"argv": argv, "output": argv[argv.index("--output-directory") + 1],
              "metrics": argv[argv.index("--json") + 1]}
    artifacts = verify_process(config, status["methods"][cell["label"]])
    integrity = {"scheduler": scheduler, "accounting": accounting, "execution_status": checked[0],
                 "postflight": checked[1], "artifact_count": len(artifacts), "executor_commit": EXECUTOR_COMMIT}
    adapted = {**manifest, "fasta_inputs": manifest["input_fastas"], "launcher": sources[0]["executed"]}
    native = validate_native_cell(adapted, environment, {**cell, "argv": argv}, output, launcher, integrity,
                                  expected_revision=LAUNCHER_COMMIT)
    directory = Path(config["output"]) / "orthohmm_phylogeny"
    native_manifest = json.loads((directory / "provenance_manifest.json").read_text())
    owners, candidates = gene_ownership(manifest, native_manifest, Path(cell["candidate_partition"]))
    pairs = directory / "orthohmm_pairwise_orthologs.tsv"
    pair_record = record(pairs)
    count = check_pairs(pairs, owners, candidates, native_manifest["results"]["ortholog_pairs"])
    if type(count) is not int or count <= 0:
        raise ValueError("No native ortholog predictions")
    checked.extend([pair_record, *[expected[k] for k in (
        "candidate_admission", "prepared_manifest", "environment", "executor", "execution_helper", "command_helper")]])
    for pair in sources:
        checked.extend([pair["prepared"], pair["executed"]])
    for item in checked:
        check(item)
    verify_process(config, status["methods"][cell["label"]])
    verify_admission(root, admission_path, admission_sha, admission_job, index)
    verify_environment(environment)
    return {"status": "corrected_qfo_native_pair_output_verified", "cell": cell["label"], "index": index,
            "source": record(__file__), "candidate_admission": record(admission_path),
            "prepared": admission["prepared_manifest"], "environment": record(environment_path),
            "native_group_integrity": native, "native_pairs": pair_record, "native_pair_count": count,
            "scheduler": scheduler, "accounting": accounting, "checked_records": checked,
            "accuracy_evaluated": False, "scoring_admitted": False, "publication_ready": False,
            "helpers": [record(Path(__file__).with_name(n)) for n in (
                "admit_qfo_factorial_cell.py", "validate_factorial_native.py", "validate_factorial_partition.py",
                "run_qfo_corrected_factorial_cell.py", "run_qfo_factorial_cell.py", "simulation_method_outputs.py")],
            "limitations": ["Native integrity, not independent tree correctness or event-pair reconstruction.",
                "QfO R-on requires native pairs, not RootHOG clique conversion.",
                "Reference mapping and independent assessment admission remain required."]}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("root", "admission", "output"):
        parser.add_argument("--" + name, type=Path, required=True)
    for name in ("job", "admission-sha256", "admission-job"):
        parser.add_argument("--" + name, required=True)
    parser.add_argument("--index", type=int, choices=range(4), required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    result = admit(args.root.resolve(), args.index, args.job, args.admission.resolve(), args.admission_sha256, args.admission_job)
    with args.output.open("x") as stream:
        json.dump(result, stream, indent=2, sort_keys=True)
        stream.write("\n")
