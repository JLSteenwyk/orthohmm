"""Admit corrected-QfO candidate-variant native outputs without scoring."""

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
from benchmark_tools.prepare_ob_candidate_neighborhood import check, record
from benchmark_tools.run_ob_candidate_neighborhood import variant_cell
from benchmark_tools.run_qfo_factorial_cell import native_command
from benchmark_tools.run_qfo_parameter_phylogeny import VARIANTS, verify_sources
from benchmark_tools.run_simulation_methods import execution_environment
from benchmark_tools.validate_factorial_native import validate_native_cell
from benchmark_tools.validate_simulation_outputs import verify_process
from benchmark_tools.verify_qfo_replay_launcher import LAUNCHER_COMMIT

EXECUTOR_COMMIT = "aa8c0e1937b898a9da83c69bf36ec342a4e04b89"
# Replacement after the original array failed before inference; executor unchanged.
JOB = "22034"


def completed_task(accounting, index):
    if type(index) is not int or index not in range(4):
        raise ValueError("Unknown candidate variant")
    rows = [r for r in csv.DictReader(io.StringIO(accounting), delimiter="|")
            if r["JobID"] == f"{JOB}_{index}"]
    if len(rows) != 1 or tuple(rows[0][k] for k in ("State", "ExitCode", "NodeList", "AllocCPUS")) != (
            "COMPLETED", "0:0", "bizon", "32"):
        raise ValueError("Require successfully completed 32-CPU candidate task")
    return rows[0]


def check_execution(status, preflight, postflight, expected, manifest, cell):
    inputs = {"status": "ready", "inputs": [{**r, "absolute_path": r["path"]} for r in manifest["input_fastas"]]}
    if (preflight != expected or status["provenance"] != expected
            or status["verified_inputs"] != inputs or status["dataset"] != cell["label"]
            or set(status["methods"]) != {cell["label"]}
            or status["status"] != "finished_pending_native_validation" or status["failed_methods"] != []
            or status["accuracy_evaluated"] is not False or status["native_outputs_validated"] is not False):
        raise ValueError("Execution identity, provenance or success differs")
    if postflight != {"status": "complete_pending_native_validation", "cell": cell,
                      "accuracy_evaluated": False, "native_outputs_validated": False}:
        raise ValueError("Successful unscored postflight required")


def adapted_manifest(manifest, arm, launcher_record):
    result = deepcopy(manifest)
    result["fasta_inputs"] = result["input_fastas"]
    result["launcher"] = launcher_record
    target = result["candidate_arms"]["p1_c1"]
    target["candidate_partition"] = arm["partition"]
    target["membership_constraints"] = arm["constraints"]
    return result


def admit(root, index):
    if type(index) is not int or index not in range(4):
        raise ValueError("Unknown candidate variant")
    accounting = subprocess.check_output(["sacct", "-j", JOB, "--parsable2",
        "--format=JobID,JobIDRaw,State,ExitCode,Elapsed,NodeList,AllocCPUS"], text=True)
    scheduler = completed_task(accounting, index)
    validator_source = record(__file__)
    validator_helpers = [record(m.__file__) for n, m in sorted(sys.modules.items())
                         if n.startswith("benchmark_tools.") and getattr(m, "__file__", None)]
    arm, manifest, original, launcher, prepared, environment, inputs, admission_scheduler = verify_sources(root, index)
    executor = root / "benchmarks/work/publication_qfo_parameter_phylogeny_v1"
    if subprocess.check_output(["git", "-C", str(executor), "rev-parse", "HEAD"], text=True).strip() != EXECUTOR_COMMIT:
        raise ValueError("Changed executor revision")
    subprocess.run(["git", "-C", str(executor), "diff", "--exit-code", "HEAD", "--", "benchmark_tools", "orthohmm"], check=True)
    output = root / "benchmarks/results/qfo_parameter_phylogeny_v1" / VARIANTS[index]
    cell = variant_cell(original, arm, output)
    argv, equivalence = native_command(cell, launcher, prepared)
    _, resolved = execution_environment(environment)
    records = [record(output / name) for name in ("preflight.json", "postflight.json", "execution/status.json")]
    preflight, postflight, status = [json.loads(Path(r["path"]).read_text()) for r in records]
    helpers = preflight["helpers"]
    paths = [Path(r["path"]) for r in helpers]
    required = {"run_simulation_methods.py", "run_ob_candidate_neighborhood.py", "run_qfo_factorial_cell.py",
                "run_qfo_corrected_factorial_cell.py", "validate_simulation_outputs.py"}
    if (not required.issubset({p.name for p in paths}) or len(paths) != len(set(paths))
            or any(p.parent != executor / "benchmark_tools" for p in paths)):
        raise ValueError("Missing or foreign execution helpers")
    for item in helpers:
        check(item)
    expected = {"source": record(executor / "benchmark_tools/run_qfo_parameter_phylogeny.py"),
        "helpers": helpers, "inputs": inputs, "candidate_admission_scheduler": admission_scheduler,
        "arm": arm, "cell": cell, "executed_argv": argv, "launcher_source_equivalence": equivalence,
        "resolved_tools": resolved, "cwd": str(launcher), "job_id": scheduler["JobIDRaw"],
        "array_task_id": str(index),
        "scope": "Corrected candidate variant; independently inferred phylogeny; unscored incremental execution"}
    check_execution(status, preflight, postflight, expected, manifest, cell)
    config = {"argv": argv, "output": str(output / "output"), "metrics": str(output / "metrics.json")}
    artifacts = verify_process(config, status["methods"][cell["label"]])
    integrity = {"scheduler": scheduler, "accounting": accounting, "execution_status": records[2],
                 "postflight": records[1], "artifact_count": len(artifacts), "executor_commit": EXECUTOR_COMMIT}
    native = validate_native_cell(adapted_manifest(manifest, arm, equivalence[0]["executed"]),
        environment, {**cell, "argv": argv}, output, launcher, integrity, expected_revision=LAUNCHER_COMMIT)
    directory = output / "output/orthohmm_phylogeny"
    native_manifest = json.loads((directory / "provenance_manifest.json").read_text())
    owners, candidates = gene_ownership(manifest, native_manifest, Path(arm["partition"]["path"]))
    pair_record = record(directory / "orthohmm_pairwise_orthologs.tsv")
    count = check_pairs(Path(pair_record["path"]), owners, candidates, native_manifest["results"]["ortholog_pairs"])
    if type(count) is not int or count <= 0:
        raise ValueError("No native pairs")
    checked = [*records, *inputs, *helpers, expected["source"], pair_record]
    for pair in equivalence:
        checked.extend([pair["prepared"], pair["executed"]])
    for item in checked:
        check(item)
    verify_process(config, status["methods"][cell["label"]])
    verify_sources(root, index)
    for item in [validator_source, *validator_helpers]:
        check(item)
    return {"status": "corrected_qfo_parameter_native_pairs_verified", "variant": arm["label"], "index": index,
        "cell": cell, "scheduler": scheduler, "accounting": accounting, "native_group_integrity": native,
        "native_pairs": pair_record, "native_pair_count": count, "checked_records": checked,
        "source": validator_source, "helpers": validator_helpers,
        "accuracy_evaluated": False, "scoring_admitted": False, "publication_ready": False,
        "limitations": ["Integrity validation, not independent reconstruction of trees or orthology events.",
                        "Native pairs, not RootHOG cliques; reference mapping and scoring admission remain required.",
                        "Incremental shared-host runtime is not controlled comparative timing."]}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", required=True, type=Path)
    parser.add_argument("--index", required=True, type=int, choices=range(4))
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    result = admit(args.root.resolve(), args.index)
    with args.output.open("x") as stream:
        json.dump(result, stream, indent=2, sort_keys=True)
        stream.write("\n")
