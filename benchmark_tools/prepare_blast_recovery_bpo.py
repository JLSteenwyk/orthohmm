"""Convert an independently admitted recovered search into a fresh BPO checkpoint."""

import argparse
import csv
import io
import os
from pathlib import Path
import subprocess
import sys
import time

from benchmark_tools.prepare_ob_candidate_neighborhood import check, record
from benchmark_tools.prepare_orthomcl_bpo_checkpoint import prepare as checkpoint
from benchmark_tools.run_blast_recovery_batch import save_status
from benchmark_tools.run_simulation_methods import read_frozen
from benchmark_tools.verify_orthomcl_python_runtime import verify_runtime

ADMITTER = "a21c65f2b449d828e3fefc86148cbf4fd87b6ade"


def validate_admission(admission, root):
    if (admission["status"] != "recovered_orthomcl_search_evidence_verified"
            or admission["search_admitted"] is not True
            or any(admission[k] is not False for k in (
                "accuracy_admitted", "publication_ready", "downstream_execution_authorized"))):
        raise ValueError("Require independently admitted recovered search evidence")
    expected = dict(JobID="22150", State="COMPLETED", ExitCode="0:0", NodeList="bizon", AllocCPUS="2", ReqMem="64G")
    if any(admission["scheduler"].get(k) != v for k, v in expected.items()):
        raise ValueError("Wrong admitted merge scheduler identity")
    if (admission["query_coverage"]["input_proteins"] != 984137
            or admission["database_content"]["input_sequences"] != 984137
            or admission["database_content"]["exact_sequence_parity"] is not True):
        raise ValueError("Wrong recovered search universe or database parity")
    paths = [root / "benchmarks/results/qfo_blast_recovery_merge_v1/table/all.blast.candidate",
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
    expected_path = root / "benchmarks/results/qfo_blast_recovery_search_admission_v1/report.json"
    if path != expected_path:
        raise ValueError("Unexpected recovered admission path")
    admission = read_frozen(path, digest)
    inputs = validate_admission(admission, root)
    accounting = subprocess.check_output(["sacct", "-j", "22151", "--parsable2",
        "--format=JobID,State,ExitCode,NodeList,AllocCPUS,ReqMem,Elapsed"], text=True)
    rows = [r for r in csv.DictReader(io.StringIO(accounting), delimiter="|") if r["JobID"] == "22151"]
    if len(rows) != 1 or tuple(rows[0][k] for k in (
            "State", "ExitCode", "NodeList", "AllocCPUS", "ReqMem")) != (
            "COMPLETED", "0:0", "bizon", "2", "64G"):
        raise ValueError("Require successfully completed recovery validator 22151")
    executor = root / "benchmarks/work/blast_recovery_search_admission_v1_20260923"
    if subprocess.check_output(["git", "-C", str(executor), "rev-parse", "HEAD"], text=True).strip() != ADMITTER:
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
