"""Independently admit a recovered BPO checkpoint without releasing native jobs."""

import argparse
import json
import math
import os
from pathlib import Path
import subprocess
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.admit_qfo_corrected_bpo import validate_checkpoint, recheck_content
from benchmark_tools.prepare_blast_recovery_bpo import verify_admission
from benchmark_tools.prepare_ob_candidate_neighborhood import check, record
from benchmark_tools.prepare_qfo_corrected_orthomcl_native import RUNTIME_SHA, HELPERS_SHA
from benchmark_tools.run_blast_recovery_batch import save_status
from benchmark_tools.run_simulation_methods import read_frozen
from benchmark_tools.snapshot_runtime_trees import verify
from benchmark_tools.verify_orthomcl_python_runtime import verify_runtime
from benchmark_tools.verify_ygob_validation import require_completed_job

PREPARER_SHA = "86d5ebd70f73611aee47a84567996c3cc2ae77ff3db15a68b25df50345036cf4"
VALIDATOR_SHA = "848fcfe249b39bd9965a41c7aafd9338ffd18f9ebd6cae4f530ded2e632b805c"


def execution_identity():
    job = os.environ.get("SLURM_JOB_ID", "")
    if (not job.isascii() or not job.isdigit() or int(job) < 1
            or os.environ.get("SLURM_CPUS_PER_TASK") != "2"
            or os.environ.get("SLURM_MEM_PER_NODE") != "65536"
            or os.uname().nodename != "bizon"):
        raise ValueError("Require scheduled two-CPU 64-GiB recovered BPO admission on bizon")
    return dict(admission_job_id=job, node="bizon", allocated_cpus=2, memory_mib=65536)


def preparation_contract(report, scheduler, source, checkpoint, runtime_record):
    expected = dict(State="COMPLETED", ExitCode="0:0", NodeList="bizon", AllocCPUS="2", ReqMem="64G")
    if any(scheduler.get(k) != v for k, v in expected.items()):
        raise ValueError("Require completed two-CPU recovered BPO preparation")
    if (report["status"] != "recovered_bpo_prepared_pending_admission"
            or report["job_id"] != scheduler["JobIDRaw"] or report["source"] != source
            or report["checkpoint"] != checkpoint
            or any(report[k] is not False for k in ("accuracy_admitted", "publication_ready", "downstream_execution_authorized"))):
        raise ValueError("Recovered preparation identity/status differs")
    times = [report["started_epoch"], report["finished_epoch"]]
    if any(type(t) not in (int, float) or not math.isfinite(t) for t in times) or not 0 < times[0] <= times[1]:
        raise ValueError("Invalid recovered preparation timestamps")
    for key in ("runtime_before", "runtime_after"):
        runtime = report[key]
        if runtime["status"] != "dedicated_bpo_python_runtime_verified" or runtime["manifest"] != runtime_record:
            raise ValueError("Wrong recovered preparation runtime")
        for item in runtime["mapped_files"]:
            check(item)


def admit(root, job, executor, commit, output):
    if output.exists() or output.is_symlink():
        raise FileExistsError(output)
    identity = execution_identity()
    runtime = verify_runtime(root)
    accounting = subprocess.check_output(["sacct", "-j", str(job), "--parsable2",
        "--format=JobIDRaw,State,ExitCode,Elapsed,NodeList,AllocCPUS,ReqMem"], text=True)
    scheduler = require_completed_job(accounting, job)
    if not executor.is_relative_to(root / "benchmarks/work"):
        raise ValueError("Require retained preparation executor")
    if subprocess.check_output(["git", "-C", str(executor), "rev-parse", "HEAD"], text=True).strip() != commit:
        raise ValueError("Changed recovered preparation executor")
    subprocess.run(["git", "-C", str(executor), "diff", "--exit-code", "HEAD", "--", "benchmark_tools", "orthohmm"], check=True)
    source = record(executor / "benchmark_tools/prepare_blast_recovery_bpo.py")
    helper = record(Path(__file__).with_name("admit_qfo_corrected_bpo.py"))
    if source["sha256"] != PREPARER_SHA or helper["sha256"] != VALIDATOR_SHA:
        raise ValueError("Unreviewed recovered preparation or content validator")
    base = root / "benchmarks/results/qfo_blast_recovery_bpo_v1"
    parent_record, checkpoint_record = record(base / "report.json"), record(base / "checkpoint/report.json")
    parent = read_frozen(Path(parent_record["path"]), parent_record["sha256"])
    checkpoint = read_frozen(Path(checkpoint_record["path"]), checkpoint_record["sha256"])
    preparation_contract(parent, scheduler, source, checkpoint_record, runtime["manifest"])
    validate_checkpoint(checkpoint, base / "checkpoint", executor)
    admission_record = parent["recovered_search"]
    admission, inputs, search_checked, admission_scheduler, _ = verify_admission(
        root, Path(admission_record["path"]), admission_record["sha256"])
    if (admission_record not in parent["checked_records"] or admission_scheduler != parent["admission_scheduler"]
            or parent["query_coverage"] != admission["query_coverage"]
            or parent["content"] != checkpoint["content"] or parent["index_validation"] != checkpoint["index_validation"]
            or checkpoint["content"]["input_proteins"] != 984137):
        raise ValueError("Recovered search scope or checkpoint summaries differ")
    for actual, expected in (("source_hsp_rows", "hsp_rows"), ("source_pair_blocks", "distinct_directed_pairs")):
        if checkpoint["content"][actual] != admission["query_coverage"][expected]:
            raise ValueError("Recovered source row/pair counts differ")
    manifests = [read_frozen(root / "benchmark_tools/results" / name, sha) for name, sha in (
        ("qfo_corrected_orthomcl_perl_runtime_20260918.json", RUNTIME_SHA),
        ("qfo_corrected_orthomcl_system_helpers_20260918.json", HELPERS_SHA))]
    for manifest in manifests:
        verify(manifest)
    checked = [record(__file__), source, helper, parent_record, checkpoint_record,
        *parent["checked_records"], *search_checked, *checkpoint["checked_records"], *checkpoint["outputs"],
        *[record(path) for path in sorted(Path(__file__).parent.glob("*.py"))],
        *[record(Path(__file__).with_name(name)) for name in ("run_orthomcl_perl_script.pl", "validate_orthomcl_bpo_indexes.pl")]]
    for item in checked:
        check(item)
    output.mkdir(parents=True, exist_ok=False)
    result = dict(status="recovered_bpo_validation_running", **identity, source=record(__file__), scheduler=scheduler,
        accounting=accounting, preparation=parent_record, checkpoint=checkpoint_record, checked_records=checked,
        recovered_search=admission_record, query_coverage=admission["query_coverage"], runtime_before=runtime,
        checkpoint_admitted=False, accuracy_admitted=False, publication_ready=False, downstream_execution_authorized=False)
    save_status(output / "report.json", result)
    try:
        validation = recheck_content(Path(inputs[0]["path"]), Path(inputs[1]["path"]), base / "checkpoint",
                                     output / "recheck", checkpoint["content"], checkpoint["index_validation"])
        for manifest in manifests:
            verify(manifest)
        for item in [*checked, *validation["outputs"]]:
            check(item)
        result.update(status="recovered_orthomcl_bpo_checkpoint_admitted", checkpoint_admitted=True,
            validation=validation, runtime_after=verify_runtime(root), native_inputs=[record(base / "checkpoint" / name)
                for name in ("all.bpo", "indexes/all_bpo.idx", "indexes/all_bpo.se")], limitations=[
                "Independent content/index recheck; not biological accuracy validation.",
                "Original and recovered BLAST failures remain in query coverage; conversion does not repair missing hits.",
                "Native inference requires a separate provenance-aware handoff; old held jobs remain untouched."])
    except BaseException as error:
        result.update(status="recovered_bpo_validation_failed", error_type=type(error).__name__, error=str(error))
        raise
    finally:
        save_status(output / "report.json", result)
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for flag in ("root", "executor", "output"):
        parser.add_argument("--" + flag, type=Path, required=True)
    parser.add_argument("--job", type=int, required=True)
    parser.add_argument("--commit", required=True)
    args = parser.parse_args()
    admit(args.root.resolve(), args.job, args.executor.resolve(), args.commit, args.output.absolute())
