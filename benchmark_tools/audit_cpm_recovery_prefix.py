"""Revalidate retained high-CPM predecessors without admitting a recovery."""

import argparse
import json
from pathlib import Path
import subprocess
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.checked_replay_payload_worker import FILES, STAGES
from benchmark_tools.prepare_ob_candidate_neighborhood import check, record
from benchmark_tools.validate_checked_replay_payload import validate


def require_failed_worker(worker):
    calls = worker["calls"]
    if (worker["status"] != "failed" or worker["accuracy_evaluated"] is not False
            or len(calls) != 4
            or [(c["index"], c["stage"]) for c in calls] != list(enumerate(STAGES))
            or any(type(c["index"]) is not int or c["accuracy_evaluated"] is not False for c in calls)
            or any(c["status"] != "checked" or type(c["exit_code"]) is not int
                   or c["exit_code"] != 0 for c in calls[:3])
            or calls[3]["status"] != "failed"):
        raise ValueError("Require three checked predecessors and failed fourth stage")
    return calls


def audit_prefix(root, executor, directory, worker, plan_record):
    """Check native predecessors; caller must separately bind job/runtime/context."""
    calls = require_failed_worker(worker)
    check(plan_record)
    records, summaries = [plan_record], []
    for row in calls[:3]:
        stage = directory / "clustering" / f"cluster_{row['index']}_{row['stage']}"
        payload = stage / "payload"
        execution_record = record(stage / "execution.json")
        if json.loads(Path(execution_record["path"]).read_text()) != row:
            raise ValueError("Predecessor execution differs from parent")
        manifest_record = record(stage / "payload_manifest.json")
        if row["manifest"] != manifest_record:
            raise ValueError("Predecessor manifest changed")
        manifest = json.loads(Path(manifest_record["path"]).read_text())
        inputs = [record(payload / name) for name in FILES]
        if (type(manifest["index"]) is not int or manifest["index"] != row["index"]
                or manifest["stage"] != row["stage"] or manifest["accuracy_evaluated"] is not False
                or manifest["output_directory"] != str(directory / "replay")
                or manifest["inputs"] != inputs or len(manifest["original_payload"]) != len(FILES)):
            raise ValueError("Predecessor payload inventory differs")
        original_dirs = {Path(r["path"]).parent for r in manifest["original_payload"]}
        if len(original_dirs) != 1:
            raise ValueError("Ambiguous original payload directory")
        if manifest["original_command"] != [sys.executable, "-m", "orthohmm.leiden_worker",
                                            str(next(iter(original_dirs)))]:
            raise ValueError("Original predecessor command differs")
        for name, prior, copied in zip(FILES, manifest["original_payload"], inputs):
            if (Path(prior["path"]).name != name
                    or any(prior[k] != copied[k] for k in ("bytes", "sha256"))):
                raise ValueError("Original/copied predecessor payload differs")
        expected = [sys.executable, str(executor / "benchmark_tools/checked_replay_payload_worker.py"),
            "--root", str(root), "--payload", str(payload), "--manifest", manifest_record["path"],
            "--manifest-sha256", manifest_record["sha256"], "--corrected-plan", plan_record["path"],
            "--corrected-plan-sha256", plan_record["sha256"], "--cpm-arm", "cpm_high"]
        overrides = {k: "1" for k in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS")}
        threads = {"child_overrides": overrides,
                   "inherited": {**overrides, "OMP_NUM_THREADS": "32" if row["index"] == 2 else "1"}}
        if row["command"] != expected or row["thread_environment"] != threads:
            raise ValueError("Predecessor command or threads differ")
        partition = record(stage / "partition.txt")
        if partition != row["partition"]:
            raise ValueError("Retained predecessor partition changed")
        fresh = validate(payload, manifest, root, executor, corrected_plan=plan_record,
                         retained_partition=partition, cpm_arm="cpm_high")
        # Native validation occurred before its partition was copied for retention.
        original_validation = {**fresh, "provenance_checked": fresh["provenance_checked"][:-1]}
        if (fresh["provenance_checked"][-1] != partition or row["validation"] != original_validation
                or fresh["status"] != "payload_checked" or fresh["accuracy_evaluated"] is not False
                or fresh["genes"] != 984137 or fresh["saved_graph"]["vertices"] != 984137):
            raise ValueError("Fresh predecessor validation differs")
        summaries.append({"index": row["index"], "stage": row["stage"], "partition": partition,
                          "genes": fresh["genes"], "groups": fresh["groups"], "saved_graph": fresh["saved_graph"]})
        records.extend([execution_record, manifest_record, partition, record(stage / "worker.log"),
                        *fresh["provenance_checked"]])
    multipass = record(directory / "replay/orthogroups_multipass.txt")
    if any(multipass[k] != summaries[1]["partition"][k] for k in ("bytes", "sha256")):
        raise ValueError("Multipass output differs from retained predecessor")
    records.append(multipass)
    unique = {}
    for item in records:
        if item["path"] in unique and unique[item["path"]] != item:
            raise ValueError("Conflicting predecessor provenance")
        unique[item["path"]] = item
    for item in unique.values():
        check(item)
    return {"status": "cpm_recovery_predecessors_verified_not_admitted", "clustering": summaries,
            "checked_records": list(unique.values()), "recovery_authorized": False,
            "accuracy_evaluated": False, "publication_ready": False,
            "limitations": ["Only stages 0-2: original job, runtime, failed graph and refinement require separate gates.",
                            "No optimizer continuation, recovered partition or accuracy result is admitted."]}


def run(root, output):
    from benchmark_tools.cpm_replay_context import evidence, REPLAY_SHA
    from benchmark_tools.checked_replay_payload_worker import corrected_evidence
    from benchmark_tools.run_blast_recovery_batch import save_status

    if output.exists() or output.is_symlink():
        raise FileExistsError(output)
    executor = root / "benchmarks/work/publication_qfo_cpm_variant_v3"
    revision = subprocess.check_output(["git", "-C", str(executor), "rev-parse", "HEAD"], text=True).strip()
    if revision != "a2486a9ef39afc2035c3fd20942696be48c5647d":
        raise ValueError("Changed original CPM executor")
    subprocess.run(["git", "-C", str(executor), "diff", "--exit-code", "HEAD", "--",
                    "benchmark_tools", "orthohmm"], check=True)
    plan, plan_record, admission_record, names_record = corrected_evidence(
        root / "benchmarks/work/qfo_corrected_replay_commands_20260918.json", REPLAY_SHA)
    context = evidence(root, plan, plan_record, "cpm_high")
    directory = Path(context["output_root"])
    worker_record = record(directory / "checked_worker.json")
    worker = json.loads(Path(worker_record["path"]).read_text())
    if (worker["context"] != context
            or worker["source"] != record(executor / "benchmark_tools/run_qfo_cpm_control.py")
            or worker["replay_source"] != record(Path(context["cwd"]) / "benchmark_tools/replay_high_sensitivity.py")):
        raise ValueError("Changed failed worker context or sources")
    helpers = [record(Path(__file__)), record(Path(__file__).with_name("validate_checked_replay_payload.py"))]
    result = audit_prefix(root, executor, directory, worker, plan_record)
    result.update(source=helpers[0], original_executor=str(executor), original_executor_commit=revision,
                  failed_worker=worker_record, context=context)
    result["checked_records"].extend([worker_record, worker["source"], worker["replay_source"],
        admission_record, names_record, *helpers, *context["checked_records"]])
    for item in result["checked_records"]:
        check(item)
    save_status(output, result)
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    run(args.root.resolve(), args.output.absolute())
