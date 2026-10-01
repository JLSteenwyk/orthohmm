"""Independently admit the private full-QfO deployment control, without scoring."""

import argparse
import csv
import importlib.util
import io
import json
from pathlib import Path
import subprocess
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.prepare_ob_candidate_neighborhood import check, record

JOB = "22387"
COMMIT = "bf3cc9e0a1f754492280f74d39c763cab49c4006"
EXECUTOR = "benchmarks/work/qfo_private_phylogeny_control_executor_20261001"
OUTPUT = "benchmarks/results/qfo_private_phylogeny_control_v1"
SUBMISSION = "benchmark_tools/results/qfo_private_phylogeny_control_submission_22387.json"
PROTOCOL = "benchmark_tools/results/QFO_PRIVATE_PHYLOGENY_ADMISSION_PROTOCOL_20261001.md"
DRIVER_SHA = "2fcee25ba93964544df53690aa546dc8fff3fe36953c360a88595ead7bb13c96"
ENVIRONMENT_SHA = "ee61bc0cd524e0e85aa8ac42991d00f3c64baf07ee0a0a22d66e1120449a331e"
CONTROL_PROTOCOL_SHA = "da7bf35d089f97671804090c953af27c8559b6430ad60b27be6dcee5a97c1756"
LABEL = "private_phylogeny_control"


def completed(accounting):
    rows = list(csv.DictReader(io.StringIO(accounting), delimiter="|"))
    parents = [r for r in rows if r["JobID"] == JOB]
    if len(parents) != 1 or tuple(parents[0][k] for k in
            ("JobIDRaw", "State", "ExitCode", "AllocCPUS", "ReqMem", "NodeList")) != (
                JOB, "COMPLETED", "0:0", "32", "192G", "bizon"):
        raise ValueError("Require completed private QfO control " + JOB)
    steps = [r for r in rows if r["JobID"].startswith(JOB + ".")]
    if not any(r["JobID"] == JOB + ".batch" for r in steps) or any(
            r["State"] != "COMPLETED" or r["ExitCode"] != "0:0" for r in steps):
        raise ValueError("Require successful terminal private QfO control steps")
    return parents[0]


def accounting():
    return subprocess.check_output(["sacct", "-j", JOB, "--parsable2",
        "--format=JobID,JobIDRaw,State,ExitCode,Elapsed,AllocCPUS,ReqMem,NodeList"], text=True)


def read(path):
    return json.loads(Path(path).read_bytes())


def execution_identity(preflight, postflight, status, expected, inputs):
    if (preflight != expected or status["provenance"] != expected
            or status["verified_inputs"] != inputs or status["dataset"] != LABEL
            or set(status["methods"]) != {LABEL} or status["failed_methods"] != []
            or status["status"] != "finished_pending_native_validation"
            or status["accuracy_evaluated"] is not False or status["native_outputs_validated"] is not False):
        raise ValueError("Private control execution/provenance differs")
    if (set(postflight) != {"status", "native_comparison", "accuracy_evaluated", "publication_ready",
                           "recovered_cpm_inference_authorized", "controlled_timing"}
            or postflight["status"] != "private_qfo_baseline_parity_complete_pending_admission"
            or any(postflight[k] is not False for k in
                ("accuracy_evaluated", "publication_ready", "recovered_cpm_inference_authorized", "controlled_timing"))):
        raise ValueError("Private control lacks successful unscored parity postflight")


def cache_summary(summary):
    keys = ("candidate_families", "reconciled_families", "bypassed_families", "checkpoint_hits",
            "remapped_checkpoint_hits", "species_tree_families")
    if (any(type(summary[k]) is not int or summary[k] < 0 for k in keys)
            or summary["candidate_families"] != 351739
            or summary["reconciled_families"] + summary["bypassed_families"] != summary["candidate_families"]
            or not 0 <= summary["remapped_checkpoint_hits"] <= summary["checkpoint_hits"] <= summary["reconciled_families"]
            or type(summary["species_tree_checkpoint_hit"]) is not bool):
        raise ValueError("Invalid full-baseline checkpoint accounting")
    return {k: summary[k] for k in (*keys, "species_tree_checkpoint_hit")}


def control_environment(root, executor):
    path = executor / "benchmark_tools/qfo_private_phylogeny_environment.py"
    if record(path)["sha256"] != ENVIRONMENT_SHA:
        raise ValueError("Control environment verifier changed")
    spec = importlib.util.spec_from_file_location("verified_private_qfo_control_environment", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module, module.verify_baseline(root)


def admit(root, submission_sha, protocol_sha, destination):
    if not destination.is_absolute() or destination.resolve() != destination:
        raise ValueError("Require direct absolute private QfO admission output")
    if destination.exists() or destination.is_symlink():
        raise FileExistsError(destination)
    observed = accounting()
    scheduler = completed(observed)
    protocol = record(root / PROTOCOL)
    if protocol["sha256"] != protocol_sha:
        raise ValueError("Unreviewed private QfO control admission protocol")
    submission_pin = record(root / SUBMISSION)
    if submission_pin["sha256"] != submission_sha:
        raise ValueError("Wrong private QfO control submission identity")
    submission = read(submission_pin["path"])
    executor = root / EXECUTOR
    if (submission["job_id"] != JOB or submission["status"] != "private_full_qfo_phylogeny_control_submitted"
            or submission["executor_commit"] != COMMIT or submission["executor"] != str(executor)
            or any(submission[k] is not False for k in
                ("accuracy_evaluated", "recovered_cpm_inference_authorized", "controlled_timing", "publication_ready"))):
        raise ValueError("Wrong private QfO control submission")
    if subprocess.check_output(["git", "-C", str(executor), "rev-parse", "HEAD"], text=True).strip() != COMMIT:
        raise ValueError("Private QfO control executor revision changed")
    subprocess.run(["git", "-C", str(executor), "diff", "--exit-code", "HEAD", "--", "benchmark_tools", "orthohmm"],
                   check=True, capture_output=True)
    for item in submission["source_records"]:
        check(item)
    source = record(executor / "benchmark_tools/run_qfo_private_phylogeny_control.py")
    if source["sha256"] != DRIVER_SHA or source not in submission["source_records"]:
        raise ValueError("Wrong private QfO control driver binding")
    module, verified = control_environment(root, executor)
    output = root / OUTPUT
    provenance_pins = [record(output / name) for name in ("preflight.json", "postflight.json", "execution/status.json")]
    preflight, postflight, status = [read(item["path"]) for item in provenance_pins]
    helpers = [record(path) for path in sorted((executor / "benchmark_tools").glob("*.py"))]
    if preflight["helpers"] != helpers:
        raise ValueError("Private QfO execution helper inventory differs")
    for item in helpers:
        current = record(Path(__file__).parent / Path(item["path"]).name)
        if (current["bytes"], current["sha256"]) != (item["bytes"], item["sha256"]):
            raise ValueError("Admission helper differs from frozen control: " + Path(item["path"]).name)
    launcher, prepared = Path(verified["launcher"]), Path(verified["prepared"])
    argv, equivalence = module.control_command(verified["original"], launcher, prepared, verified["environment"], output)
    from benchmark_tools.run_simulation_methods import execution_environment
    from benchmark_tools.validate_simulation_outputs import verify_process
    from benchmark_tools.validate_factorial_native import validate_native_cell
    from benchmark_tools.admit_qfo_corrected_factorial_cell import gene_ownership
    from benchmark_tools.admit_qfo_factorial_cell import check_pairs
    from benchmark_tools.verify_qfo_replay_launcher import LAUNCHER_COMMIT
    from benchmark_tools.inspect_native_python_lookup import PROBE, scientific_origins
    from benchmark_tools.run_qfo_private_phylogeny_control import compare_outputs

    _, resolved = execution_environment(verified["environment"])
    control_protocol = record(root / "benchmark_tools/results/QFO_PRIVATE_PHYLOGENY_CONTROL_PROTOCOL_20261001.md")
    if control_protocol["sha256"] != CONTROL_PROTOCOL_SHA:
        raise ValueError("Private control protocol changed")
    expected = {"source": source, "helpers": helpers, "protocol": control_protocol, "verified": verified,
        "executed_argv": argv, "launcher_equivalence": equivalence, "resolved_tools": resolved,
        "cwd": str(launcher), "job_id": JOB,
        "scope": "Entire frozen QfO p1_c1_r1 private deployment parity; validated checkpoint reuse; unscored incremental shared-host control"}
    inputs = {"status": "ready", "inputs": [{**r, "absolute_path": r["path"]} for r in verified["manifest"]["input_fastas"]]}
    execution_identity(preflight, postflight, status, expected, inputs)
    config = {"argv": argv, "output": str(output / "output"), "metrics": str(output / "metrics.json")}
    artifacts = verify_process(config, status["methods"][LABEL])
    lookup_pin, process_pin = record(output / "lookup.json"), record(output / "lookup_process.json")
    lookup, process = read(lookup_pin["path"]), read(process_pin["path"])
    report = lookup["report"]
    lookup_command = [argv[0], "-B", "-c", PROBE, json.dumps(["orthohmm.phylogeny_pipeline", "benchmark_tools.replay_phylogeny"])]
    if (process["returncode"] != 0 or process["command"] != lookup_command or json.loads(process["stdout"]) != report
            or report["cwd"] != str(launcher) or report["executable"] != str(Path(argv[0]).resolve())
            or report["requested"] != ["orthohmm.phylogeny_pipeline", "benchmark_tools.replay_phylogeny"]
            or report["dont_write_bytecode"] is not True or report["pycache_prefix"] != str(output / "bytecode_cache")
            or lookup["continuous_enforcement"] is not False
            or report["modules"].get("benchmark_tools.replay_phylogeny") != str(launcher / "benchmark_tools/replay_phylogeny.py")
            or lookup["origin"] != scientific_origins(report, "orthohmm", launcher / "orthohmm")):
        raise ValueError("Private QfO launcher lookup provenance differs")
    paths = sorted(set(report["modules"].values()) | set(report["mapped_files"]) | {report["executable"]})
    if (any(p.startswith(("/home/bizon/anaconda3/", "/home/bizon/.local/")) for p in paths)
            or lookup["checked_records"] != [record(path) for path in paths]):
        raise ValueError("Private QfO observed runtime changed or uses retired roots")
    cell = {**verified["original"], "label": LABEL, "argv": argv,
            "checkpoint_source": argv[argv.index("--checkpoint-source") + 1]}
    manifest = verified["manifest"]
    adapted = {**manifest, "fasta_inputs": manifest["input_fastas"], "launcher": equivalence[0]["executed"]}
    integrity = {"scheduler": scheduler, "accounting": observed, "execution_status": provenance_pins[2],
        "postflight": provenance_pins[1], "artifact_count": len(artifacts), "executor_commit": COMMIT}
    native = validate_native_cell(adapted, verified["environment"], cell, output, launcher, integrity,
                                  expected_revision=LAUNCHER_COMMIT)
    directory = output / "output/orthohmm_phylogeny"
    native_manifest = read(directory / "provenance_manifest.json")
    cache = cache_summary(native_manifest["results"])
    owners, candidates = gene_ownership(manifest, native_manifest, Path(cell["candidate_partition"]))
    pair_pin = record(directory / "orthohmm_pairwise_orthologs.tsv")
    count = check_pairs(Path(pair_pin["path"]), owners, candidates, native_manifest["results"]["ortholog_pairs"])
    if type(count) is not int or count <= 0:
        raise ValueError("Private control has no native pairs")
    comparisons = compare_outputs(Path(cell["checkpoint_source"]) / "orthohmm_phylogeny", directory)
    if comparisons != postflight["native_comparison"] or not all(item["byte_equal"] for item in comparisons.values()):
        raise ValueError("Independent full private QfO output comparison failed")
    records = [record(__file__), protocol, submission_pin, *submission["source_records"],
        submission["deployment"], submission["original_native_admission"], submission["preserved_shared_rejection"],
        *provenance_pins, *helpers, *verified["checked_records"], lookup_pin, process_pin, *lookup["checked_records"], pair_pin,
        native["native_manifest"], native["native_metrics"], native["species_tree"],
        *[record(Path(__file__).parent / name) for name in ("validate_factorial_native.py", "validate_factorial_partition.py",
            "admit_qfo_corrected_factorial_cell.py", "admit_qfo_factorial_cell.py", "simulation_method_outputs.py", "run_qfo_private_phylogeny_control.py")]]
    for item in comparisons.values():
        records.extend([item["historical"], item["current"]])
    for pair in equivalence:
        records.extend([pair["prepared"], pair["executed"]])
    unique = {}
    for item in records:
        if item["path"] in unique and unique[item["path"]] != item:
            raise ValueError("Conflicting private QfO admission provenance")
        unique[item["path"]] = item
        check(item)
    if module.verify_baseline(root) != verified or completed(accounting()) != scheduler:
        raise ValueError("Private QfO runtime/baseline/completion changed during admission")
    verify_process(config, status["methods"][LABEL])
    result = {"status": "private_qfo_phylogeny_deployment_admitted_unscored", "job_id": JOB,
        "source": record(__file__), "protocol": protocol, "submission": submission_pin,
        "scheduler": scheduler, "accounting": observed, "native_group_integrity": native,
        "native_pairs": pair_pin, "native_pair_count": count, "cache_use": cache,
        "native_comparison": comparisons, "checked_records": list(unique.values()),
        "recovered_cpm_inference_authorized": True, "accuracy_evaluated": False, "scoring_admitted": False,
        "controlled_timing": False, "publication_ready": False,
        "limitations": ["Deployment admission under validated checkpoint reuse, not fresh all-tree/search equivalence.",
            "Native integrity and byte parity, not independent gene-tree/event correctness or new reference accuracy.",
            "Import-chain/file identity checks are not continuous access enforcement or a hermetic OS snapshot.",
            "Recovered candidates require their own fixed admission and explicit private command handoff.",
            "Incremental shared-host resources do not establish controlled comparative timing."]}
    with destination.open("x") as stream:
        json.dump(result, stream, indent=2, sort_keys=True, allow_nan=False)
        stream.write("\n")
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--submission-sha256", required=True)
    parser.add_argument("--protocol-sha256", required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    admit(args.root.resolve(), args.submission_sha256, args.protocol_sha256, args.output.absolute())
