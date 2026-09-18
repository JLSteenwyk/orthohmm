"""Independently compare corrected search checkpoints against all source hits."""

import argparse
import json
from pathlib import Path
import sqlite3
import subprocess
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.admit_qfo_sequence_search_control import validate_execution
from benchmark_tools.convert_qfo_sequence_search_control import validate_metadata, ADMITTER
from benchmark_tools.prepare_ob_candidate_neighborhood import record, check
from benchmark_tools.run_qfo_sequence_search_control import verify_plan
from benchmark_tools.run_simulation_methods import read_frozen
from benchmark_tools.verify_sequence_checkpoint_equivalence import reconstruct, verify
from benchmark_tools.verify_ygob_validation import require_completed_job

CONVERTER = "f5c23b81e3cf72573fafd5d1bb087fd2e61392d2"


def completed(job, cpus, memory):
    text = subprocess.check_output(["sacct", "-j", str(job), "--parsable2",
        "--format=JobIDRaw,State,ExitCode,Elapsed,NodeList,AllocCPUS,ReqMem"], text=True)
    row = require_completed_job(text, job)
    if (row["NodeList"], row["AllocCPUS"]) != ("bizon", str(cpus)) or row["ReqMem"] not in (memory, memory + "n"):
        raise ValueError("Unexpected completed allocation")
    return row


def validate_conversion(report, scheduler, source, output_root):
    if (report["status"] != "corrected_qfo_sequence_numeric_checkpoints_verified"
            or report["source"] != source or report["job_id"] != scheduler["JobIDRaw"]
            or report["numeric_validated"] is not True or report["accuracy_evaluated"] is not False
            or report["publication_ready"] is not False or report["genes"] != 984137
            or report["proteomes"] != 78 or type(report["hits"]) is not int or report["hits"] <= 0
            or set(report["variants"]) != {"all_hits", "top100"}):
        raise ValueError("Wrong corrected numeric conversion identity")
    for label, cap in (("all_hits", None), ("top100", 100)):
        variant = report["variants"][label]
        expected = output_root / label / "orthohmm_working_res/high_sensitivity_checkpoint"
        if (variant["checkpoint"] != str(expected) or variant["cap"] != cap
                or type(variant["cap"]) is not type(cap)
                or variant["manifest"]["path"] != str(expected / "manifest.json")
                or variant["audit"]["status"] != "numeric_checkpoint_verified"
                or variant["audit"]["manifest"] != variant["manifest"]):
            raise ValueError("Numeric variant identity differs")


def admit(root, job, output):
    if output.exists() or output.is_symlink():
        raise FileExistsError(output)
    scheduler = completed(job, 2, "192G")
    executor = root / "benchmarks/work/publication_qfo_sequence_numeric_v1"
    if subprocess.check_output(["git", "-C", str(executor), "rev-parse", "HEAD"], text=True).strip() != CONVERTER:
        raise ValueError("Frozen numeric converter changed")
    subprocess.run(["git", "-C", str(executor), "diff", "--exit-code", "HEAD", "--", "benchmark_tools", "orthohmm"], check=True)
    numeric_root = root / "benchmarks/results/qfo_sequence_numeric_v1"
    manifest_path = numeric_root / "manifest.json"
    manifest_record = record(manifest_path)
    report = read_frozen(manifest_path, manifest_record["sha256"])
    validate_conversion(report, scheduler, record(executor / "benchmark_tools/convert_qfo_sequence_search_control.py"), numeric_root)
    expected_helpers = [record(executor / "benchmark_tools" / n) for n in (
        "convert_sequence_search_control.py", "audit_accuracy_checkpoint.py", "admit_qfo_sequence_search_control.py")]
    if (report["helpers"] != expected_helpers or report["checkpoint_writer"] != record(executor / "orthohmm/accuracy.py")
            or report["admission"] not in report["checked_records"]):
        raise ValueError("Conversion helper or admission binding differs")
    prior = completed(report["admission_scheduler"]["JobIDRaw"], 2, "64G")
    if prior != report["admission_scheduler"]:
        raise ValueError("Search-admission accounting changed")
    checked = [manifest_record, *report["checked_records"], report["source"], *report["helpers"], report["checkpoint_writer"]]
    for item in checked:
        check(item)
    admission_record = report["admission"]
    check(admission_record)
    admission = read_frozen(Path(admission_record["path"]), admission_record["sha256"])
    auditor = root / "benchmarks/work/publication_qfo_sequence_search_admission_v1"
    if (subprocess.check_output(["git", "-C", str(auditor), "rev-parse", "HEAD"], text=True).strip() != ADMITTER
            or admission["source"] != record(auditor / "benchmark_tools/admit_qfo_sequence_search_control.py")
            or admission["status"] != "corrected_qfo_search_panel_admitted_pending_numeric_validation"):
        raise ValueError("Wrong independent search admission")
    subprocess.run(["git", "-C", str(auditor), "diff", "--exit-code", "HEAD", "--", "benchmark_tools", "orthohmm"], check=True)
    plan_path = root / "benchmarks/work/qfo_sequence_search_control_v1/manifest.json"
    plan = verify_plan(plan_path)
    execution_record = admission["execution"]
    check(execution_record)
    execution = read_frozen(Path(execution_record["path"]), execution_record["sha256"])
    search_scheduler = completed(admission["scheduler"]["JobIDRaw"], 32, "192G")
    checked.extend(validate_execution(execution, plan, search_scheduler, record(plan_path),
                                     root / "benchmarks/work/publication_qfo_sequence_search_v1"))
    checked.extend([admission_record, execution_record])
    metadata_record = plan["gene_metadata"]
    metadata = read_frozen(Path(metadata_record["path"]), metadata_record["sha256"])
    validate_metadata(plan["inputs"], metadata)
    checked.append(metadata_record)
    for variant in report["variants"].values():
        check(variant["manifest"])
        checked.append(variant["manifest"])
    for item in checked:
        check(item)
    output.mkdir(parents=True, exist_ok=False)
    result = {"status": "verifying", "source": record(__file__), "conversion": manifest_record,
        "scheduler": scheduler, "search_admission": admission_record, "checked_records": checked,
        "accuracy_evaluated": False, "publication_ready": False,
        "helpers": [record(Path(__file__).with_name(n)) for n in (
            "verify_sequence_checkpoint_equivalence.py", "convert_qfo_sequence_search_control.py",
            "admit_qfo_sequence_search_control.py", "run_qfo_sequence_search_control.py")]}
    try:
        with sqlite3.connect(output / "expected.sqlite") as database:
            database.execute("PRAGMA temp_store=FILE")
            database.execute("PRAGMA cache_size=-131072")
            sources = [(Path(t["output"]), Path(t["target_fasta"]["path"]).name) for t in plan["searches"]]
            count = reconstruct(database, sources, metadata, sorted(metadata))
            if count != report["hits"]:
                raise ValueError("Source reconstruction and converter counts differ")
            variants = {label: verify(database, Path(value["checkpoint"]), metadata, value["cap"])
                        for label, value in report["variants"].items()}
        for item in [*checked, result["source"], *result["helpers"]]:
            check(item)
        result.update(status="corrected_qfo_numeric_source_equivalence_admitted", variants=variants,
                      hits=count, genes=len(metadata), proteomes=78, numeric_equivalence=True,
                      limitations=["Exact equivalence to emitted search hits, not proof of biological search completeness.",
                          "Graph replay, native output validation and benchmark scoring remain separate gates.",
                          "Does not establish equal HMM/DIAMOND sensitivity or computational cost."])
    except BaseException as error:
        result.update(status="failed", numeric_equivalence=False, error_type=type(error).__name__, error=str(error))
        raise
    finally:
        with (output / "report.json").open("x") as stream:
            json.dump(result, stream, indent=2, sort_keys=True)
            stream.write("\n")
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for flag in ("root", "output"):
        parser.add_argument("--" + flag, type=Path, required=True)
    parser.add_argument("--job", required=True)
    args = parser.parse_args()
    admit(args.root.resolve(), args.job, args.output.absolute())
