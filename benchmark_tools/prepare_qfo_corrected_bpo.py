"""Bind corrected BLAST admission to a separately audited native BPO checkpoint."""

import argparse
import json
import os
from pathlib import Path
import subprocess
import sys
import time

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.prepare_ob_candidate_neighborhood import record, check
from benchmark_tools.prepare_orthomcl_bpo_checkpoint import prepare as checkpoint
from benchmark_tools.run_simulation_methods import read_frozen
from benchmark_tools.verify_ygob_validation import require_completed_job

ADMITTER = "ae7de310cf7fbcba59107f2b054fc186f206eec4"


def validate_admission(admission, root):
    if (admission["status"] != "corrected_orthomcl_search_evidence_verified"
            or admission["search_admitted"] is not True
            or admission["accuracy_admitted"] is not False
            or admission["publication_ready"] is not False
            or admission["downstream_execution_authorized"] is not False):
        raise ValueError("Require independently admitted corrected BLAST evidence")
    scheduler = admission["scheduler"]
    expected = {"State": "COMPLETED", "ExitCode": "0:0", "NodeList": "bizon",
                "AllocCPUS": "180", "ReqMem": "900G"}
    if any(scheduler.get(k) != v for k, v in expected.items()):
        raise ValueError("Wrong native BLAST scheduler allocation")
    if (admission["query_coverage"]["input_proteins"] != 984137
            or admission["database_content"]["input_sequences"] != 984137):
        raise ValueError("Wrong corrected BLAST input universe")
    work = root / "benchmarks/results/qfo_corrected_orthomcl_v1/work"
    inputs = []
    for name in ("all.blast", "all.fa"):
        matches = [r for r in admission["checked_records"] if r["path"] == str(work / name)]
        if not matches or any(r != matches[0] for r in matches) or matches[0]["bytes"] <= 0:
            raise ValueError("Missing/conflicting admitted BLAST input: " + name)
        inputs.append(matches[0])
    return inputs


def verify_admission(root, path, sha, job):
    admission = read_frozen(path, sha)
    inputs = validate_admission(admission, root)
    accounting = subprocess.check_output(["sacct", "-j", str(job), "--parsable2",
        "--format=JobIDRaw,State,ExitCode,Elapsed,NodeList,AllocCPUS,ReqMem"], text=True)
    scheduler = require_completed_job(accounting, job)
    if any(scheduler.get(k) != v for k, v in {"NodeList": "bizon", "AllocCPUS": "2", "ReqMem": "64G"}.items()):
        raise ValueError("Wrong BLAST admission allocation")
    executor = root / "benchmarks/work/publication_qfo_corrected_blast_admission_v1"
    if subprocess.check_output(["git", "-C", str(executor), "rev-parse", "HEAD"], text=True).strip() != ADMITTER:
        raise ValueError("Changed frozen search admission executor")
    subprocess.run(["git", "-C", str(executor), "diff", "--exit-code", "HEAD", "--", "benchmark_tools"], check=True)
    if admission["source"] != record(executor / "benchmark_tools/admit_qfo_corrected_blast.py"):
        raise ValueError("Wrong search admission source")
    checked = [record(path), admission["source"], *admission["checked_records"],
               admission["database_audit"], admission["table_audit"]]
    for item in checked:
        check(item)
    return admission, inputs, checked, scheduler, accounting


def prepare(root, admission_path, admission_sha, admission_job):
    if (not os.environ.get("SLURM_JOB_ID") or os.environ.get("SLURM_CPUS_PER_TASK") != "2"
            or os.environ.get("SLURM_MEM_PER_NODE") != "65536" or os.uname().nodename != "bizon"):
        raise ValueError("Require scheduled two-CPU 64-GiB allocation on bizon")
    admission, inputs, checked, scheduler, accounting = verify_admission(
        root, admission_path, admission_sha, admission_job)
    output = root / "benchmarks/results/qfo_corrected_orthomcl_v1/bpo_preparation"
    output.mkdir(parents=True, exist_ok=False)
    helpers = Path(__file__).resolve().parent
    checked.extend(record(p) for p in (Path(__file__), Path(sys.executable),
        helpers / "prepare_orthomcl_bpo_checkpoint.py", helpers / "build_orthomcl_bpo_indexes.pl"))
    report = {"status": "preparing", "source": record(__file__), "checked_records": checked,
              "admission_scheduler": scheduler, "admission_accounting": accounting,
              "admitted_search_scheduler": admission["scheduler"],
              "query_coverage": admission["query_coverage"], "job_id": os.environ["SLURM_JOB_ID"],
              "started_epoch": time.time(), "python_version": sys.version,
              "accuracy_admitted": False, "publication_ready": False}
    try:
        result = checkpoint(root, Path(inputs[0]["path"]), Path(inputs[1]["path"]), output / "checkpoint")
        if (result["status"] != "bpo_checkpoint_content_and_indexes_verified"
                or result["content"]["input_proteins"] != 984137):
            raise ValueError("Wrong BPO input universe or checkpoint status")
        for actual, expected in (("source_hsp_rows", "hsp_rows"),
                                 ("source_pair_blocks", "distinct_directed_pairs")):
            if result["content"][actual] != admission["query_coverage"][expected]:
                raise ValueError("BPO source counts differ from admitted search")
        for item in [*checked, *result["checked_records"], *result["outputs"]]:
            check(item)
        report.update(status="corrected_bpo_checkpoint_prepared_pending_admission",
                      checkpoint=record(output / "checkpoint/report.json"),
                      content=result["content"], index_validation=result["index_validation"])
        report["limitations"] = [
            "Retains source query failures; BPO conversion does not repair absent outgoing hits.",
            "No final orthogroup inference, scoring or matched-resource timing is performed.",
            "Recorded Python binary/version is not a complete frozen Python runtime; an execution freeze is required before production launch.",
            "Terminal scheduler and independent checkpoint admission remain required before inference."]
    except BaseException as error:
        report.update(status="failed", error_type=type(error).__name__, error=str(error))
        raise
    finally:
        report["finished_epoch"] = time.time()
        with (output / "report.json").open("x") as stream:
            json.dump(report, stream, indent=2, sort_keys=True)
            stream.write("\n")
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("root", "admission"):
        parser.add_argument("--" + name, type=Path, required=True)
    parser.add_argument("--admission-sha256", required=True)
    parser.add_argument("--admission-job", required=True, type=int)
    args = parser.parse_args()
    result = prepare(args.root.resolve(), args.admission.resolve(), args.admission_sha256, args.admission_job)
    print(json.dumps({"status": result["status"]}))
