"""Independently admit recovered private high-CPM native groups/pairs, unscored."""

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

JOB = "22390"
COMMIT = "8952ebed34921f873eac83ce63319ee3d7d8dbf7"
EXECUTOR = "benchmarks/work/qfo_private_cpm_phylogeny_executor_20261001"
SUBMISSION = "benchmark_tools/results/qfo_private_cpm_phylogeny_submission_22390.json"
PROTOCOL = "benchmark_tools/results/QFO_PRIVATE_CPM_NATIVE_ADMISSION_PROTOCOL_20261001.md"
PRODUCER_SHA = "ea4e6c2e8ab74d4bcc9ecc2f09d625e32ed2f8b32005e8be90812eacfe3192fe"
PRODUCER_PROTOCOL = "benchmark_tools/results/QFO_PRIVATE_CPM_PHYLOGENY_PROTOCOL_20261001.md"
PRODUCER_PROTOCOL_SHA = "d35db243dfed5f5724cfd0610e87c1642c8b9d7eee36d6ef6715ccbbcdf82aa4"


def accounting():
    return subprocess.check_output(["sacct", "-j", JOB, "-P",
        "--format=JobID,JobIDRaw,State,ExitCode,Elapsed,AllocCPUS,ReqMem,NodeList"], text=True)


def completed(text):
    rows = list(csv.DictReader(io.StringIO(text), delimiter="|"))
    selected = [row for row in rows if row["JobID"] == JOB]
    if len(selected) != 1 or tuple(selected[0][k] for k in
            ("JobIDRaw", "State", "ExitCode", "AllocCPUS", "ReqMem", "NodeList")) != (
                JOB, "COMPLETED", "0:0", "32", "192G", "bizon"):
        raise ValueError("Require successfully completed private recovered native job " + JOB)
    steps = [row for row in rows if row["JobID"].startswith(JOB + ".")]
    if not any(row["JobID"] == JOB + ".batch" for row in steps) or any(
            row["State"] != "COMPLETED" or row["ExitCode"] != "0:0" for row in steps) or len(
            {row["JobID"] for row in steps}) != len(steps):
        raise ValueError("Require successful terminal private recovered native steps")
    return selected[0]


def read(path):
    return json.loads(Path(path).read_bytes())


def producer(root, executor):
    path = executor / "benchmark_tools/run_private_helper_cpm_phylogeny.py"
    if record(path)["sha256"] != PRODUCER_SHA:
        raise ValueError("Private recovered producer source changed")
    spec = importlib.util.spec_from_file_location("verified_private_cpm_producer", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module, module.verify_sources(root)


def execution_identity(preflight, postflight, status, expected, manifest, cell, command, fresh):
    from benchmark_tools.admit_qfo_parameter_phylogeny import check_execution

    native_post = dict(status="complete_pending_native_validation", cell=cell,
                       accuracy_evaluated=False, native_outputs_validated=False)
    if (any(postflight.get(k) is not False for k in
            ("accuracy_evaluated", "native_outputs_validated", "publication_ready"))
            or postflight != {**native_post, "publication_ready": False,
                             "admission_command": command, "fresh_candidate_admission": fresh}):
        raise ValueError("Require exact successful unscored fresh-candidate postflight")
    check_execution(status, preflight, native_post, expected, manifest, cell)


def cache_accounting(summary):
    names = ("candidate_families", "reconciled_families", "bypassed_families", "checkpoint_hits",
             "remapped_checkpoint_hits", "species_tree_families")
    if (any(type(summary[k]) is not int or summary[k] < 0 for k in names)
            or summary["candidate_families"] != 346866
            or summary["reconciled_families"] + summary["bypassed_families"] != summary["candidate_families"]
            or not 0 <= summary["remapped_checkpoint_hits"] <= summary["checkpoint_hits"] <= summary["reconciled_families"]
            or summary["species_tree_families"] > summary["candidate_families"]
            or type(summary["species_tree_checkpoint_hit"]) is not bool):
        raise ValueError("Invalid recovered high-CPM checkpoint accounting")
    return {k: summary[k] for k in (*names, "species_tree_checkpoint_hit")}


def lookup_identity(output, argv, launcher):
    from benchmark_tools.inspect_native_python_lookup import PROBE, scientific_origins

    lookup_pin, process_pin = record(output / "lookup.json"), record(output / "lookup_process.json")
    lookup, process = read(lookup_pin["path"]), read(process_pin["path"])
    report = lookup["report"]
    requested = ["orthohmm.phylogeny_pipeline", "benchmark_tools.replay_phylogeny"]
    if (process["returncode"] != 0 or process["command"] != [argv[0], "-B", "-c", PROBE, json.dumps(requested)]
            or json.loads(process["stdout"]) != report or report["requested"] != requested
            or report["executable"] != str(Path(argv[0]).resolve()) or report["cwd"] != str(launcher)
            or report["dont_write_bytecode"] is not True or report["pycache_prefix"] != str(output / "bytecode_cache")
            or lookup["continuous_enforcement"] is not False
            or report["modules"].get("orthohmm.phylogeny_pipeline") != str(launcher / "orthohmm/phylogeny_pipeline.py")
            or report["modules"].get("benchmark_tools.replay_phylogeny") != str(launcher / "benchmark_tools/replay_phylogeny.py")
            or lookup["origin"] != scientific_origins(report, "orthohmm", launcher / "orthohmm")):
        raise ValueError("Private recovered replay lookup differs")
    paths = sorted(set(report["modules"].values()) | set(report["mapped_files"]) | {report["executable"]})
    if (any(path.startswith(("/home/bizon/anaconda3/", "/home/bizon/.local/")) for path in paths)
            or lookup["checked_records"] != [record(path) for path in paths]):
        raise ValueError("Private observed runtime changed or uses retired roots")
    return [lookup_pin, process_pin, *lookup["checked_records"]]


def admit(root, submission_sha, protocol_sha, destination):
    if not destination.is_absolute() or destination.resolve() != destination:
        raise ValueError("Require direct absolute private recovered admission output")
    if destination.exists() or destination.is_symlink():
        raise FileExistsError(destination)
    observed = accounting()
    scheduler = completed(observed)
    protocol, submission_pin = record(root / PROTOCOL), record(root / SUBMISSION)
    if protocol["sha256"] != protocol_sha or submission_pin["sha256"] != submission_sha:
        raise ValueError("Unreviewed private recovered admission protocol/submission")
    submission = read(submission_pin["path"])
    executor = root / EXECUTOR
    if (submission["status"] != "private_recovered_qfo_high_cpm_phylogeny_submitted"
            or submission["job_id"] != JOB or submission["executor_commit"] != COMMIT
            or submission["executor"] != str(executor) or submission["executor_clean"] is not True
            or any(submission[k] is not False for k in ("native_completion_observed", "native_outputs_validated",
                "accuracy_evaluated", "controlled_timing", "publication_ready"))):
        raise ValueError("Wrong private recovered native submission identity")
    if subprocess.check_output(["git", "-C", str(executor), "rev-parse", "HEAD"], text=True).strip() != COMMIT:
        raise ValueError("Private recovered executor revision changed")
    subprocess.run(["git", "-C", str(executor), "diff", "--exit-code", "HEAD", "--", "benchmark_tools", "orthohmm"],
                   check=True, capture_output=True)
    for item in submission["source_records"]:
        check(item)
    source = record(executor / "benchmark_tools/run_private_helper_cpm_phylogeny.py")
    if source["sha256"] != PRODUCER_SHA or source not in submission["source_records"]:
        raise ValueError("Private recovered submitted source changed")
    module, verified = producer(root, executor)
    baseline, candidates = verified["baseline"], verified["candidates"]
    output = root / module.OUTPUT
    records = [record(output / name) for name in ("preflight.json", "postflight.json", "execution/status.json")]
    preflight, postflight, status = [read(pin["path"]) for pin in records]
    helpers = [record(path) for path in sorted((executor / "benchmark_tools").glob("*.py"))]
    if preflight["helpers"] != helpers:
        raise ValueError("Private recovered execution helper inventory differs")
    for item in helpers:
        current = record(Path(__file__).parent / Path(item["path"]).name)
        if (current["bytes"], current["sha256"]) != (item["bytes"], item["sha256"]):
            raise ValueError("Native admission helper differs from frozen producer: " + Path(item["path"]).name)
    cell, planned, argv, equivalence = module.native_command(verified, output)
    from benchmark_tools.run_simulation_methods import execution_environment
    from benchmark_tools.validate_simulation_outputs import verify_process
    from benchmark_tools.admit_qfo_cpm_phylogeny import adapted_manifest
    from benchmark_tools.validate_factorial_native import validate_native_cell
    from benchmark_tools.admit_qfo_corrected_factorial_cell import gene_ownership
    from benchmark_tools.admit_qfo_factorial_cell import check_pairs
    from benchmark_tools.verify_qfo_replay_launcher import LAUNCHER_COMMIT

    _, resolved = execution_environment(baseline["environment"])
    producer_protocol = record(root / PRODUCER_PROTOCOL)
    if producer_protocol["sha256"] != PRODUCER_PROTOCOL_SHA:
        raise ValueError("Private recovered native protocol changed")
    expected = dict(source=source, protocol=producer_protocol, helpers=helpers, verified=verified,
        cell=cell, planned_cell=planned, executed_argv=argv, launcher_source_equivalence=equivalence,
        resolved_tools=resolved, cwd=baseline["launcher"], job_id=JOB,
        seed_handoff="explicit_helper_runtime_seed_amendment",
        native_handoff="explicit_admitted_private_phylogeny_deployment",
        scope="Unscored fixed recovered high-CPM arm; inferred phylogeny with validated raw-tree checkpoint reuse; incremental shared-host execution")
    fresh = record(output / "fresh_candidate_admission.json")
    command = [module.PYTHON, "-B", str(Path(candidates["admission_executor"]) / "benchmark_tools/admit_helper_cpm_candidates.py"),
        "--root", str(root), "--preparation-sha256", candidates["admission"]["preparation"]["sha256"],
        "--protocol-sha256", candidates["admission"]["protocol"]["sha256"], "--output", fresh["path"]]
    execution_identity(preflight, postflight, status, expected, baseline["manifest"], cell, command, fresh)
    if read(fresh["path"]) != candidates["admission"]:
        raise ValueError("Private recovered fresh candidate admission disagrees")
    config = dict(argv=argv, output=str(output / "output"), metrics=str(output / "metrics.json"))
    artifacts = verify_process(config, status["methods"][cell["label"]])
    launcher = Path(baseline["launcher"])
    lookup = lookup_identity(output, argv, launcher)
    adapted = adapted_manifest(dict(manifest=baseline["manifest"], arm=candidates["arm"], admission=candidates["admission"]),
                               equivalence[0]["executed"])
    integrity = dict(scheduler=scheduler, accounting=observed, execution_status=records[2], postflight=records[1],
                     artifact_count=len(artifacts), executor_commit=COMMIT)
    native = validate_native_cell(adapted, baseline["environment"], cell, output, launcher, integrity,
                                  expected_revision=LAUNCHER_COMMIT)
    directory = output / "output/orthohmm_phylogeny"
    native_manifest = read(directory / "provenance_manifest.json")
    cache = cache_accounting(native_manifest["results"])
    owners, families = gene_ownership(baseline["manifest"], native_manifest, Path(candidates["arm"]["partition"]["path"]))
    pairs = record(directory / "orthohmm_pairwise_orthologs.tsv")
    count = check_pairs(Path(pairs["path"]), owners, families, native_manifest["results"]["ortholog_pairs"])
    if type(count) is not int or count <= 0:
        raise ValueError("Private recovered native arm has no pairs")
    checked = [record(__file__), protocol, submission_pin, *submission["source_records"], source,
        submission["source_preflight"], submission["private_admission"], submission["private_readback"], submission["candidate_admission"],
        producer_protocol, *records, *helpers, *verified["checked_records"], fresh, record(output / "admission.log"),
        *lookup, pairs, native["native_manifest"], native["native_metrics"], native["species_tree"]]
    for pair in equivalence:
        checked.extend([pair["prepared"], pair["executed"]])
    unique = {}
    for item in checked:
        if item["path"] in unique and unique[item["path"]] != item:
            raise ValueError("Conflicting private recovered native identity")
        unique[item["path"]] = item
        check(item)
    verify_process(config, status["methods"][cell["label"]])
    if module.verify_sources(root) != verified or completed(accounting()) != scheduler:
        raise ValueError("Private recovered evidence/runtime/completion changed during admission")
    for item in unique.values():
        check(item)
    result = dict(status="private_recovered_cpm_native_pairs_verified_unscored", arm="cpm_high", index=1,
        cell=cell, scheduler=scheduler, accounting=observed, candidate_admission=candidates["admission_record"],
        private_native_admission=verified["private_control"]["admission"], source=record(__file__),
        protocol=protocol, submission=submission_pin, native_group_integrity=native, native_pairs=pairs,
        native_pair_count=count, cache_use=cache, checked_records=list(unique.values()),
        accuracy_evaluated=False, scoring_admitted=False, controlled_timing=False, publication_ready=False,
        limitations=["Native integrity, not independent gene-tree/event correctness or new reference accuracy.",
            "Recovered helper seeds/private native deployment are explicitly amended; original failed gates remain unchanged.",
            "Native pairs, not RootHOG cliques; lossless conversion and independent assessment admission remain required.",
            "Validated checkpoint reuse/shared-host incremental resources do not establish controlled scaling."])
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
