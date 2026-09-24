"""Assemble fresh preflight evidence for the single high-CPM continuation."""

import csv
import io
import json
from pathlib import Path
import subprocess

from benchmark_tools.prepare_ob_candidate_neighborhood import check, record
from benchmark_tools.run_simulation_methods import read_frozen
from benchmark_tools.run_blast_recovery_batch import save_status
from benchmark_tools.validate_cpm_refinement_reconstruction import validate as validate_refinement
from benchmark_tools.audit_cpm_recovery_prefix import run as audit_predecessors
from benchmark_tools.audit_cpm_high_failed_payload import audit as audit_failure
from benchmark_tools.probe_cpm_high_constructor import AUDIT_SHA, HELPERS, require_audit
from benchmark_tools.probe_cpm_worker_boundary import require_stop
from benchmark_tools.verify_blast_recovery_panel import unique_records

PROTOCOL_SHA = "67615f1d7d66dd208a5e5ad6119fa0bdb282c4e3a5f412f8b95f5a86246aa0ad"
BOUNDARY_SHA = "376fc46a435f1988683630717c945f0f51751ab7f5fd5594178973bcad5e8be1"
ADAPTER_SHA = "443f2cb6b1f6f1fcc4d69a4410f0eeddaaa3cfee83fd0013f475221f7a39684f"


def failed_job(accounting):
    rows = [row for row in csv.DictReader(io.StringIO(accounting), delimiter="|") if row["JobID"] == "22081_1"]
    if len(rows) != 1 or tuple(rows[0][k] for k in ("State", "ExitCode", "AllocCPUS", "NodeList")) != (
            "FAILED", "1:0", "32", "bizon"):
        raise ValueError("Original high-CPM failure identity differs")
    return rows[0]


def failure_identity(parent, execution, manifest, worker, context, scheduler, executor):
    if (parent["status"] != "failed" or parent["arm"] != "cpm_high"
            or type(parent["index"]) is not int or parent["index"] != 1
            or parent["job_id"] != scheduler["JobIDRaw"]
            or parent["executor_commit"] != "a2486a9ef39afc2035c3fd20942696be48c5647d"
            or parent["context"] != context or parent["runtime_before"] != context["baseline_plan"]["runtime"]
            or parent["source"] != record(executor / "benchmark_tools/run_qfo_cpm_variant.py")
            or parent["accuracy_evaluated"] is not False or parent["publication_ready"] is not False):
        raise ValueError("Original failed parent differs")
    if (execution != worker["calls"][3] or execution["status"] != "failed"
            or type(execution["index"]) is not int or execution["index"] != 3
            or execution["stage"] != "profile_expanded" or "SIGSEGV" not in execution["error"]
            or execution["manifest"] != manifest):
        raise ValueError("Original fourth-stage failure differs")


def prepare(root, output):
    if output.exists() or output.is_symlink():
        raise FileExistsError(output)
    # A pending/failed reconstruction must not create a recovery directory.
    refinement = validate_refinement(root)
    if (refinement["status"] != "original_cpm_refinement_independently_verified"
            or refinement["comparison"]["partition_equal"] is not True
            or refinement["comparison"]["genes"] != 984137):
        raise ValueError("Unverified original refinement")
    results = root / "benchmark_tools/results"
    protocol = record(results / "QFO_CPM_CHECKPOINT_RECOVERY_PROTOCOL_20260923.md")
    if protocol["sha256"] != PROTOCOL_SHA:
        raise ValueError("Recovery protocol changed")
    prior_path = results / "qfo_cpm_high_failed_payload_audit_20260923.json"
    prior = read_frozen(prior_path, AUDIT_SHA)
    require_audit(prior)
    boundary_path = results / "qfo_cpm_worker_boundary_result_22152.json"
    boundary = read_frozen(boundary_path, BOUNDARY_SHA)
    if (boundary["status"] != "frozen_worker_stopped_before_optimizer_unscored"
            or boundary["job_id"] != "22152" or boundary["returncode"] != 0
            or boundary["optimizer_called"] is not False):
        raise ValueError("Invalid worker-boundary diagnostic")
    adapters = [r for r in boundary["observations"] if Path(r["path"]).name == "constructor_adapter.json"]
    if len(adapters) != 1:
        raise ValueError("Ambiguous diagnostic constructor observation")
    adapter = read_frozen(Path(adapters[0]["path"]), adapters[0]["sha256"])
    require_stop(boundary["result"], adapter, boundary["result"]["saved"])
    helpers = []
    for name, digest in {**HELPERS, "checked_python_pair_worker.py": ADAPTER_SHA}.items():
        item = record(Path(__file__).with_name(name))
        if item["sha256"] != digest:
            raise ValueError("Changed recovery native helper")
        helpers.append(item)
    accounting = subprocess.check_output(["sacct", "-j", "22081", "--parsable2",
        "--format=JobID,JobIDRaw,State,ExitCode,Elapsed,AllocCPUS,NodeList"], text=True)
    scheduler = failed_job(accounting)
    records = unique_records([record(__file__), protocol, record(prior_path), record(boundary_path),
        *helpers, *prior["checked_records"], *boundary["checked_records"], *boundary["observations"],
        boundary["worker_log"], *refinement["checked_records"]])
    for item in records:
        check(item)
    output.mkdir()
    report = dict(status="cpm_checkpoint_preflight_running", source=record(__file__),
                  checked_records=records, original_scheduler=scheduler, accounting=accounting,
                  optimizer_executed=False, accuracy_evaluated=False, publication_ready=False)
    save_status(output / "status.json", report)
    try:
        save_status(output / "refinement_validation.json", refinement)
        prefix = audit_predecessors(root, output / "predecessors.json")
        fresh_failure = audit_failure(root, output / "failed_payload.json")
        require_audit(fresh_failure)
        if fresh_failure["observed"] != prior["observed"]:
            raise ValueError("Failed graph differs from reviewed diagnostic")
        directory = Path(prefix["context"]["output_root"])
        stage = directory / "clustering/cluster_3_profile_expanded"
        paths = [directory / "results.json", directory / "checked_worker.json", stage / "execution.json",
                 stage / "payload_manifest.json"]
        original_records = [record(path) for path in paths]
        parent, worker, execution, manifest = [read_frozen(Path(r["path"]), r["sha256"]) for r in original_records]
        failure_identity(parent, execution, original_records[-1], worker, prefix["context"], scheduler,
                         Path(prefix["original_executor"]))
        if manifest["inputs"] != [record(stage / "payload" / name) for name in (
                "gene_names.txt", "sources.npy", "targets.npy", "weights.npy", "metadata.json")]:
            raise ValueError("Failed graph manifest changed")
        from benchmark_tools.probe_leiden_boundary import saved_fingerprint
        fingerprint = saved_fingerprint(stage / "payload")
        if fingerprint != boundary["result"]["saved"]:
            raise ValueError("Failed graph differs from observed full-size boundary")
        from benchmark_tools.verify_qfo_replay_launcher import verify
        runtime = verify(root / "benchmarks/work/publication_method_native_v2",
                         Path(prefix["context"]["cwd"]), results / "publication_native_runtime_20260916.json")
        if runtime != prefix["context"]["baseline_plan"]["runtime"]:
            raise ValueError("Recovery runtime changed")
        records = unique_records([*records, *prefix["checked_records"], *fresh_failure["checked_records"],
            *original_records, *[record(output / name) for name in (
                "refinement_validation.json", "predecessors.json", "failed_payload.json")]])
        for item in records:
            check(item)
        report.update(status="cpm_checkpoint_preflight_verified_unscored", checked_records=records,
            preflight_passed=True, saved_payload=str(stage / "payload"), saved_graph=fingerprint,
            context=prefix["context"], runtime=runtime, downstream_admitted=False,
            limitations=["No optimizer or accuracy evaluation has run; this is a fresh checkpoint preflight.",
                         "The original failure remains failed; lost profile statistics and timings remain unavailable."])
    except BaseException as error:
        report.update(status="cpm_checkpoint_preflight_failed", preflight_passed=False,
                      error_type=type(error).__name__, error=str(error))
        raise
    finally:
        save_status(output / "status.json", report)
    return report
