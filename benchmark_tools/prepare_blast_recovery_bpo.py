"""Convert an independently admitted recovered search into a fresh BPO checkpoint."""

import argparse
import csv
import io
import os
from pathlib import Path
import subprocess
import sys
import time

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from benchmark_tools.prepare_ob_candidate_neighborhood import check, record
from benchmark_tools.prepare_orthomcl_bpo_checkpoint import prepare as checkpoint
from benchmark_tools.run_blast_recovery_batch import save_status
from benchmark_tools.run_simulation_methods import read_frozen
from benchmark_tools.verify_orthomcl_python_runtime import verify_runtime

ADMITTER = "a21c65f2b449d828e3fefc86148cbf4fd87b6ade"
REPLACEMENT_ADMITTER = "198014bbaab12465620f1c6572a03459cdb1759e"
NATIVE_ADMITTER = "42e6dcff5170d17c490273ab93229833b91a6b23"


def admission_contract(root, path):
    original = root / "benchmarks/results/qfo_blast_recovery_search_admission_v1/report.json"
    replacement = root / "benchmarks/results/qfo_blast_replacement_search_admission_v1/report.json"
    native = root / "benchmarks/results/qfo_blast_native_representation_admission_v1/report.json"
    if path == native:
        return True, "22166", NATIVE_ADMITTER, root / "benchmarks/work/blast_native_representation_admission_v1_20260925"
    if path == original:
        return False, "22151", ADMITTER, root / "benchmarks/work/blast_recovery_search_admission_v1_20260923"
    if path == replacement:
        return True, "22163", REPLACEMENT_ADMITTER, root / "benchmarks/work/blast_replacement_search_admission_v1_20260925"
    raise ValueError("Unexpected recovered admission path")


def validate_admission(admission, root, replacement=False, native_representation=False):
    expected_status = ("recovered_search_native_representation_verified" if native_representation
                       else "recovered_orthomcl_search_evidence_verified")
    if (admission["status"] != expected_status
            or admission["search_admitted"] is not True
            or any(admission[k] is not False for k in (
                "accuracy_admitted", "publication_ready", "downstream_execution_authorized"))):
        raise ValueError("Require independently admitted recovered search evidence")
    expected = dict(JobID="22162" if replacement else "22150", State="COMPLETED", ExitCode="0:0", NodeList="bizon", AllocCPUS="2", ReqMem="64G")
    if any(admission["scheduler"].get(k) != v for k, v in expected.items()):
        raise ValueError("Wrong admitted merge scheduler identity")
    if (admission["query_coverage"]["input_proteins"] != 984137
            or admission["database_content"]["input_sequences"] != 984137):
        raise ValueError("Wrong recovered search universe or database parity")
    if native_representation:
        if not replacement:
            raise ValueError("Native representation requires replacement search")
        from benchmark_tools.reviewed_legacy_database import verify as verify_representation
        expected_representation = verify_representation(root, dict(
            status="database_sequence_differences_require_review", content=admission["database_content"]))
        actual = admission.get("database_representation", {})
        # Source paths differ across frozen executors; checked identities remain bound below.
        if {k: v for k, v in actual.items() if k != "checked_records"} != {
                k: v for k, v in expected_representation.items() if k != "checked_records"}:
            raise ValueError("Native representation evidence differs")
        for item in actual["checked_records"]:
            if item not in admission["checked_records"]:
                raise ValueError("Native representation record missing from admission")
        expected_helper = root / "benchmarks/work/blast_native_representation_admission_v1_20260925/benchmark_tools/reviewed_legacy_database.py"
        if record(expected_helper) not in actual["checked_records"]:
            raise ValueError("Native representation not bound to frozen helper")
    elif admission["database_content"]["exact_sequence_parity"] is not True:
        raise ValueError("Wrong recovered search universe or database parity")
    merge = "qfo_blast_replacement_merge_v1" if replacement else "qfo_blast_recovery_merge_v1"
    paths = [root / f"benchmarks/results/{merge}/table/all.blast.candidate",
             root / "benchmarks/results/qfo_corrected_orthomcl_v1/work/all.fa"]
    inputs = []
    for path in paths:
        matches = [r for r in admission["checked_records"] if r["path"] == str(path)]
        if not matches or any(r != matches[0] for r in matches) or matches[0]["bytes"] <= 0:
            raise ValueError("Missing/conflicting recovered BPO input")
        inputs.append(matches[0])
    if admission["candidate"] != inputs[0]:
        raise ValueError("Candidate differs from independently checked input")
    return inputs


def verify_admission(root, path, digest):
    replacement, job, commit, executor = admission_contract(root, path)
    admission = read_frozen(path, digest)
    inputs = validate_admission(admission, root, replacement, job == "22166")
    accounting = subprocess.check_output(["sacct", "-j", job, "--parsable2",
        "--format=JobID,State,ExitCode,NodeList,AllocCPUS,ReqMem,Elapsed"], text=True)
    rows = [r for r in csv.DictReader(io.StringIO(accounting), delimiter="|") if r["JobID"] == job]
    if len(rows) != 1 or tuple(rows[0][k] for k in (
            "State", "ExitCode", "NodeList", "AllocCPUS", "ReqMem")) != (
            "COMPLETED", "0:0", "bizon", "2", "64G"):
        raise ValueError("Require successfully completed recovery validator " + job)
    if subprocess.check_output(["git", "-C", str(executor), "rev-parse", "HEAD"], text=True).strip() != commit:
        raise ValueError("Changed recovery validator executor")
    subprocess.run(["git", "-C", str(executor), "diff", "--exit-code", "HEAD", "--", "benchmark_tools"], check=True)
    source = record(executor / "benchmark_tools/admit_blast_recovery_search.py")
    if source not in admission["checked_records"]:
        raise ValueError("Recovery admission is not bound to frozen validator source")
    checked = [record(path), *admission["checked_records"]]
    for item in checked:
        check(item)
    return admission, inputs, checked, rows[0], accounting


def prepare(root, path, digest):
    if (not os.environ.get("SLURM_JOB_ID") or os.environ.get("SLURM_CPUS_PER_TASK") != "2"
            or os.environ.get("SLURM_MEM_PER_NODE") != "65536" or os.uname().nodename != "bizon"):
        raise ValueError("Require scheduled two-CPU 64-GiB allocation on bizon")
    runtime_before = verify_runtime(root)
    admission, inputs, checked, scheduler, accounting = verify_admission(root, path, digest)
    output = root / "benchmarks/results/qfo_blast_recovery_bpo_v1"
    output.mkdir(exist_ok=False)
    helpers = Path(__file__).parent
    checked.extend(record(p) for p in (Path(__file__), Path(sys.executable),
        helpers / "prepare_orthomcl_bpo_checkpoint.py", helpers / "verify_orthomcl_python_runtime.py"))
    report = dict(status="preparing_recovered_bpo", source=record(__file__), checked_records=checked,
        admission_scheduler=scheduler, admission_accounting=accounting,
        query_coverage=admission["query_coverage"], recovered_search=record(path),
        job_id=os.environ["SLURM_JOB_ID"], started_epoch=time.time(), runtime_before=runtime_before,
        accuracy_admitted=False, publication_ready=False, downstream_execution_authorized=False)
    if "database_representation" in admission:
        report["database_representation"] = admission["database_representation"]
    save_status(output / "report.json", report)
    try:
        result = checkpoint(root, Path(inputs[0]["path"]), Path(inputs[1]["path"]), output / "checkpoint")
        if (result["status"] != "bpo_checkpoint_content_and_indexes_verified"
                or result["content"]["input_proteins"] != 984137):
            raise ValueError("Wrong recovered BPO checkpoint or input universe")
        for actual, expected in (("source_hsp_rows", "hsp_rows"), ("source_pair_blocks", "distinct_directed_pairs")):
            if result["content"][actual] != admission["query_coverage"][expected]:
                raise ValueError("BPO source counts differ from admitted recovered search")
        for item in [*checked, *result["checked_records"], *result["outputs"]]:
            check(item)
        report["runtime_after"] = verify_runtime(root)
        report.update(status="recovered_bpo_prepared_pending_admission",
            checkpoint=record(output / "checkpoint/report.json"), content=result["content"],
            index_validation=result["index_validation"], limitations=[
                "Recovered search failures remain explicit; conversion does not repair outgoing hits.",
                "Independent terminal checkpoint admission is required before native inference.",
                "No original search output is replaced and no held downstream job is released.",
                "Shared-host conversion accounting is not controlled end-to-end timing."])
    except BaseException as error:
        report.update(status="recovered_bpo_preparation_failed", error_type=type(error).__name__, error=str(error))
        raise
    finally:
        report["finished_epoch"] = time.time()
        save_status(output / "report.json", report)
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--admission", type=Path, required=True)
    parser.add_argument("--admission-sha256", required=True)
    args = parser.parse_args()
    prepare(args.root.resolve(), args.admission.resolve(), args.admission_sha256)
