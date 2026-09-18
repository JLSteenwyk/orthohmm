"""Validate a terminal QfO factorial run and its native pairwise predictions."""

import argparse
import csv
import io
import json
from pathlib import Path
import subprocess
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from Bio import SeqIO
from benchmark_tools.prepare_ob_candidate_neighborhood import record, check
from benchmark_tools.run_qfo_factorial_cell import select_cell, native_command, verify_prepared, ENVIRONMENT_SHA
from benchmark_tools.run_simulation_methods import read_frozen, verify_environment, execution_environment
from benchmark_tools.validate_simulation_outputs import verify_process
from benchmark_tools.validate_factorial_native import validate_native_cell
from benchmark_tools.verify_qfo_replay_launcher import LAUNCHER_COMMIT
from benchmark_tools.simulation_method_outputs import orthohmm_pairs

PREPARED_SHA = "706b07c91e9a130dae229837641a7ad62d7d36a09679e5e0daa0959e182b7d64"
EXECUTOR_COMMIT = "de3202f1c4b9a7b57ccdd903f5a362194cdd5371"
ARRAY_JOB = 21671


def require_success(accounting, index, status, postflight):
    rows = list(csv.DictReader(io.StringIO(accounting), delimiter="|"))
    matches = [r for r in rows if r["JobID"] == f"{ARRAY_JOB}_{index}"]
    if len(matches) != 1:
        raise ValueError("Missing or ambiguous scheduler task")
    row = matches[0]
    if (row["State"], row["ExitCode"], row["NodeList"], row["AllocCPUS"]) != ("COMPLETED", "0:0", "bizon", "32"):
        raise ValueError("Require successful terminal matched-resource workstation task")
    provenance = status["provenance"]
    if provenance["slurm_job_id"] != row["JobIDRaw"] or provenance["slurm_array_task_id"] != str(index):
        raise ValueError("Execution scheduler identity differs")
    if (status.get("status") != "finished_pending_native_validation" or status.get("failed_methods") != []
            or status.get("accuracy_evaluated") is not False or status.get("native_outputs_validated") is not False):
        raise ValueError("Incomplete or already-scored execution")
    if postflight != {"status": "inputs_runtime_sources_reverified", "failed_methods": [],
                      "native_outputs_validated": False, "accuracy_evaluated": False}:
        raise ValueError("Successful postflight missing")
    return row


def check_pairs(path, owners, candidates, expected_count):
    previous, count = None, 0
    for pair in orthohmm_pairs(path, owners):
        if pair[0] >= pair[1] or (previous is not None and pair <= previous):
            raise ValueError("Native pairs are not canonical, unique and sorted")
        if any(g not in candidates for g in pair) or candidates[pair[0]] != candidates[pair[1]]:
            raise ValueError("Native pair crosses candidate families")
        previous, count = pair, count + 1
    if count != expected_count:
        raise ValueError("Native pair count differs from completion summary")
    return count


def admit(root, index):
    if Path.cwd().resolve() != root:
        raise ValueError("Run from original repository verification directory")
    results = root / "benchmark_tools/results"
    prepared_path = results / "qfo_factorial_prepared_20260917.json"
    environment_path = results / "publication_variable_native_methods_20260916.json"
    prepared = read_frozen(prepared_path, PREPARED_SHA)
    environment = read_frozen(environment_path, ENVIRONMENT_SHA)
    cell, output, prepared_executor = select_cell(prepared, index)
    launcher = Path(prepared["launcher_root"])
    argv, sources = native_command(cell, launcher, prepared_executor)
    evidence = output / "execution" / cell["label"]
    status_path, post_path = evidence / "status.json", evidence / "postflight.json"
    checked = [record(status_path), record(post_path)]
    status, postflight = json.loads(status_path.read_text()), json.loads(post_path.read_text())
    accounting = subprocess.check_output(["sacct", "-j", str(ARRAY_JOB), "--parsable2",
        "--format=JobID,JobIDRaw,State,ExitCode,Elapsed,NodeList,AllocCPUS"], text=True)
    scheduler = require_success(accounting, index, status, postflight)
    executor = root / "benchmarks/work/publication_qfo_factorial_reconcile_v1"
    if subprocess.check_output(["git", "-C", str(executor), "rev-parse", "HEAD"], text=True).strip() != EXECUTOR_COMMIT:
        raise ValueError("Executor revision changed")
    subprocess.run(["git", "-C", str(executor), "diff", "--exit-code", "HEAD", "--", "benchmark_tools"], check=True)
    provenance = status["provenance"]
    for key, path in {"prepared_manifest": prepared_path, "environment": environment_path,
                      "executor": executor / "benchmark_tools/run_qfo_factorial_cell.py",
                      "execution_helper": executor / "benchmark_tools/run_simulation_methods.py"}.items():
        if provenance[key] != record(path):
            raise ValueError(f"Changed execution provenance: {key}")
    if (provenance["planned_argv"] != cell["argv"] or provenance["executed_argv"] != argv or
            provenance["launcher_source_equivalence"] != sources or provenance["cwd"] != str(launcher)):
        raise ValueError("Wrong launcher relocation or command")
    inputs = {"status": "ready", "inputs": [{**r, "absolute_path": r["path"]} for r in prepared["input_fastas"]]}
    if status["verified_inputs"] != inputs or set(status["methods"]) != {cell["label"]} or status["dataset"] != cell["label"]:
        raise ValueError("Execution input/method inventory mismatch")
    if provenance["preparation_scheduler"] != verify_prepared(prepared, cell, 21670):
        raise ValueError("Preparation scheduler changed")
    verify_environment(environment)
    _, resolved = execution_environment(environment)
    if provenance["resolved_tools"] != resolved:
        raise ValueError("Executable resolution changed")
    config = {"argv": argv, "output": argv[argv.index("--output-directory") + 1],
              "metrics": argv[argv.index("--json") + 1]}
    artifacts = verify_process(config, status["methods"][cell["label"]])
    integrity = {"scheduler": scheduler, "accounting": accounting, "execution_status": checked[0],
                 "postflight": checked[1], "artifact_count": len(artifacts), "executor_commit": EXECUTOR_COMMIT}
    adapted = {**prepared, "fasta_inputs": prepared["input_fastas"], "launcher": sources[0]["executed"]}
    native_cell = {**cell, "argv": argv}
    native = validate_native_cell(adapted, environment, native_cell, output, launcher, integrity,
                                  expected_revision=LAUNCHER_COMMIT)
    directory = Path(config["output"]) / "orthohmm_phylogeny"
    manifest = json.loads((directory / "provenance_manifest.json").read_text())
    owners = {}
    taxa = {r["filename"]: r["taxon"] for r in manifest["input_proteomes"]}
    for item in prepared["input_fastas"]:
        for sequence in SeqIO.parse(item["path"], "fasta"):
            if sequence.id in owners:
                raise ValueError("Duplicate FASTA gene")
            owners[sequence.id] = taxa[Path(item["path"]).name]
    candidates = {}
    with Path(cell["candidate_partition"]).open() as stream:
        for family, line in enumerate(stream):
            for gene in line.split():
                if gene in candidates:
                    raise ValueError("Duplicate candidate gene")
                candidates[gene] = family
    if set(candidates) != set(owners):
        raise ValueError("Candidate gene universe differs")
    pairs = directory / "orthohmm_pairwise_orthologs.tsv"
    checked.append(record(pairs))
    count = check_pairs(pairs, owners, candidates, manifest["results"]["ortholog_pairs"])
    for item in checked:
        check(item)
    return {"status": "qfo_native_pair_output_verified", "cell": cell["label"], "source": record(__file__),
            "prepared": record(prepared_path), "environment": record(environment_path), "native_group_integrity": native,
            "native_pairs": checked[-1], "native_pair_count": count, "accuracy_evaluated": False,
            "scoring_admitted": False, "helpers": [record(Path(__file__).with_name(n)) for n in
                ("validate_factorial_native.py", "validate_factorial_partition.py", "simulation_method_outputs.py", "run_qfo_factorial_cell.py")],
            "limitations": ["Native pair integrity and complete group coverage; not independent gene-tree correctness or pair-event reconstruction.",
                "RootHOGs validate grouping only; QfO R-on uses native pairs, not RootHOG clique conversion.",
                "Reference-mapping conversion and native assessment admission remain required."]}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--index", type=int, choices=range(4), required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    result = admit(args.root.resolve(), args.index)
    with args.output.open("x") as stream:
        json.dump(result, stream, indent=2, sort_keys=True)
        stream.write("\n")
