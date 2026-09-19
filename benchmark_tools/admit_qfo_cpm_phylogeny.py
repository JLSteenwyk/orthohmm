"""Independently admit CPM-specific native phylogeny and pairs, without scoring."""

import argparse
from copy import deepcopy
import csv
import io
import json
from pathlib import Path
import subprocess
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.admit_qfo_corrected_factorial_cell import gene_ownership
from benchmark_tools.admit_qfo_factorial_cell import check_pairs
from benchmark_tools.admit_qfo_parameter_phylogeny import check_execution as check_native_execution
from benchmark_tools.prepare_ob_candidate_neighborhood import check, record
from benchmark_tools.run_ob_candidate_neighborhood import variant_cell
from benchmark_tools.run_qfo_cpm_phylogeny import verify_sources
from benchmark_tools.run_qfo_cpm_variant import ARMS
from benchmark_tools.run_qfo_factorial_cell import native_command
from benchmark_tools.run_simulation_methods import execution_environment, read_frozen
from benchmark_tools.validate_factorial_native import validate_native_cell
from benchmark_tools.validate_simulation_outputs import verify_process
from benchmark_tools.verify_qfo_replay_launcher import LAUNCHER_COMMIT

JOB = "21972"
EXECUTOR_COMMIT = "f1a2bb875f5a56a827e5669780b2f1c8718dd1cd"
EXECUTOR_SHA = "f33e02be25e91728ea83edcf64eb021b990e324a9a9036657ca7fbc2511b63d6"


def completed_task(accounting, index):
    if type(index) is not int or index not in range(2):
        raise ValueError("Unknown CPM phylogeny index")
    rows = [r for r in csv.DictReader(io.StringIO(accounting), delimiter="|")
            if r["JobID"] == f"{JOB}_{index}"]
    if len(rows) != 1 or tuple(rows[0][k] for k in ("State", "ExitCode", "NodeList", "AllocCPUS")) != (
            "COMPLETED", "0:0", "bizon", "32"):
        raise ValueError("Require successfully completed 32-CPU CPM phylogeny task")
    return rows[0]


def adapted_manifest(verified, launcher_record):
    result = deepcopy(verified["manifest"])
    arm = verified["arm"]
    candidate = verified["admission"]["candidate_arm"]
    for field, key in (("seed_partition", "seed_partition"), ("candidate_partition", "partition"),
                       ("membership_constraints", "constraints")):
        if candidate[field] != arm[key]:
            raise ValueError("CPM seed/candidate admission disagrees")
    result["fasta_inputs"] = result["input_fastas"]
    result["launcher"] = launcher_record
    result["candidate_arms"]["p1_c1"] = deepcopy(candidate)
    return result


def check_execution(status, preflight, postflight, expected, manifest, cell, command, fresh):
    native_post = {"status": "complete_pending_native_validation", "cell": cell,
                   "accuracy_evaluated": False, "native_outputs_validated": False}
    if postflight != {**native_post, "admission_command": command, "fresh_candidate_admission": fresh}:
        raise ValueError("Successful postflight with exact fresh candidate admission required")
    check_native_execution(status, preflight, native_post, expected, manifest, cell)


def admit(root, index):
    if type(index) is not int or index not in range(2):
        raise ValueError("Unknown CPM phylogeny index")
    accounting = subprocess.check_output(["sacct", "-j", JOB, "--parsable2",
        "--format=JobID,JobIDRaw,State,ExitCode,Elapsed,NodeList,AllocCPUS"], text=True)
    scheduler = completed_task(accounting, index)
    source = record(__file__)
    helpers = [record(m.__file__) for n, m in sorted(sys.modules.items())
               if n.startswith("benchmark_tools.") and getattr(m, "__file__", None)]
    verified = verify_sources(root, index)
    executor = root / "benchmarks/work/publication_qfo_cpm_phylogeny_v1"
    if subprocess.check_output(["git", "-C", str(executor), "rev-parse", "HEAD"], text=True).strip() != EXECUTOR_COMMIT:
        raise ValueError("Changed CPM phylogeny executor revision")
    subprocess.run(["git", "-C", str(executor), "diff", "--exit-code", "HEAD", "--", "benchmark_tools", "orthohmm"], check=True)
    producer = record(executor / "benchmark_tools/run_qfo_cpm_phylogeny.py")
    if producer["sha256"] != EXECUTOR_SHA:
        raise ValueError("Changed CPM phylogeny producer source")
    execution_helpers = [record(p) for p in sorted((executor / "benchmark_tools").glob("*.py"))]
    output = root / "benchmarks/results/qfo_cpm_phylogeny_v1" / ARMS[index]
    launcher, prepared = Path(verified["launcher"]), Path(verified["prepared"])
    cell = variant_cell(verified["original"], verified["arm"], output)
    argv, equivalence = native_command(cell, launcher, prepared)
    _, resolved = execution_environment(verified["environment"])
    records = [record(output / name) for name in ("preflight.json", "postflight.json", "execution/status.json")]
    preflight, postflight, status = [read_frozen(Path(r["path"]), r["sha256"]) for r in records]
    expected = {"source": producer, "helpers": execution_helpers, "verified": verified, "cell": cell,
        "executed_argv": argv, "launcher_source_equivalence": equivalence, "resolved_tools": resolved,
        "cwd": str(launcher), "job_id": scheduler["JobIDRaw"], "array_task_id": str(index),
        "scope": "CPM-specific seed/candidate families; independently inferred phylogeny; no accuracy scoring"}
    fresh_path = output / "fresh_candidate_admission.json"
    fresh = record(fresh_path)
    command = [argv[0], "-B", str(Path(verified["admission_executor"]) / "benchmark_tools/admit_qfo_cpm_candidates.py"),
               "--root", str(root), "--index", str(index), "--output", str(fresh_path)]
    check_execution(status, preflight, postflight, expected, verified["manifest"], cell, command, fresh)
    if read_frozen(fresh_path, fresh["sha256"]) != verified["admission"]:
        raise ValueError("Fresh candidate admission disagrees")
    config = {"argv": argv, "output": str(output / "output"), "metrics": str(output / "metrics.json")}
    artifacts = verify_process(config, status["methods"][cell["label"]])
    integrity = {"scheduler": scheduler, "accounting": accounting, "execution_status": records[2],
                 "postflight": records[1], "artifact_count": len(artifacts), "executor_commit": EXECUTOR_COMMIT}
    native = validate_native_cell(adapted_manifest(verified, equivalence[0]["executed"]),
        verified["environment"], {**cell, "argv": argv}, output, launcher, integrity, expected_revision=LAUNCHER_COMMIT)
    directory = output / "output/orthohmm_phylogeny"
    native_manifest = json.loads((directory / "provenance_manifest.json").read_text())
    owners, candidates = gene_ownership(verified["manifest"], native_manifest, Path(verified["arm"]["partition"]["path"]))
    pairs = record(directory / "orthohmm_pairwise_orthologs.tsv")
    count = check_pairs(Path(pairs["path"]), owners, candidates, native_manifest["results"]["ortholog_pairs"])
    if type(count) is not int or count <= 0:
        raise ValueError("No native pairs")
    checked = [*records, producer, *execution_helpers, *verified["checked_records"], fresh,
               record(output / "admission.log"), pairs]
    for pair in equivalence:
        checked.extend([pair["prepared"], pair["executed"]])
    for item in checked:
        check(item)
    verify_process(config, status["methods"][cell["label"]])
    if verify_sources(root, index) != verified:
        raise ValueError("CPM phylogeny inputs changed during validation")
    for item in [source, *helpers]:
        check(item)
    return {"status": "cpm_native_pairs_verified_unscored", "arm": ARMS[index], "index": index,
        "context": verified["context"], "cell": cell, "scheduler": scheduler, "accounting": accounting,
        "candidate_admission": verified["admission_record"], "native_group_integrity": native,
        "native_pairs": pairs, "native_pair_count": count, "checked_records": checked,
        "source": source, "helpers": helpers, "accuracy_evaluated": False, "scoring_admitted": False,
        "publication_ready": False, "limitations": [
            "Integrity validation, not independent reconstruction of trees or orthology events.",
            "Native pairs, not RootHOG cliques; reference mapping and scoring admission remain required.",
            "Incremental shared-host runtime is not controlled comparative timing."]}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--index", type=int, choices=range(2), required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists() or args.output.is_symlink():
        raise FileExistsError(args.output)
    result = admit(args.root.resolve(), args.index)
    with args.output.open("x") as stream:
        json.dump(result, stream, indent=2, sort_keys=True)
        stream.write("\n")
