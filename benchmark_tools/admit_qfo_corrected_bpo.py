"""Independently recheck terminal corrected BPO preparation and native indexes."""

import argparse
import json
import math
from pathlib import Path
import subprocess
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.audit_orthomcl_bpo_content import audit
from benchmark_tools.prepare_ob_candidate_neighborhood import record, check
from benchmark_tools.prepare_orthomcl_bpo_checkpoint import native_step
from benchmark_tools.prepare_qfo_corrected_bpo import verify_admission
from benchmark_tools.prepare_qfo_corrected_orthomcl_native import TOOL, RUNTIME_SHA, HELPERS_SHA
from benchmark_tools.run_qfo_corrected_blast import environment
from benchmark_tools.run_simulation_methods import read_frozen
from benchmark_tools.snapshot_runtime_trees import verify
from benchmark_tools.verify_orthomcl_python_runtime import verify_runtime
from benchmark_tools.verify_ygob_validation import require_completed_job

EXECUTOR = "826e963cf06ed414609f0609b96c6bfd283fc2d8"
OUTPUT_NAMES = {"all.bpo", "content.json", "build.stdout", "build.stderr", "validate.stdout",
                "validate.stderr", "indexes/all_bpo.idx", "indexes/all_bpo.se"}


def validate_preparation(report, scheduler, source, checkpoint_record, runtime_record):
    if (report["status"] != "corrected_bpo_checkpoint_prepared_pending_admission"
            or report["accuracy_admitted"] is not False or report["publication_ready"] is not False
            or report["job_id"] != scheduler["JobIDRaw"]):
        raise ValueError("Wrong preparation status or job identity")
    required = {"State": "COMPLETED", "ExitCode": "0:0", "NodeList": "bizon", "AllocCPUS": "2", "ReqMem": "64G"}
    if any(scheduler.get(key) != value for key, value in required.items()):
        raise ValueError("Wrong preparation scheduler evidence")
    if report["source"] != source or report["checkpoint"] != checkpoint_record:
        raise ValueError("Changed preparation source/checkpoint identity")
    times = [report["started_epoch"], report["finished_epoch"]]
    if any(type(v) not in (int, float) or not math.isfinite(v) for v in times) or not 0 < times[0] <= times[1]:
        raise ValueError("Invalid preparation timestamps")
    for stage in ("runtime_before", "runtime_after"):
        if (report[stage]["status"] != "dedicated_bpo_python_runtime_verified"
                or report[stage]["manifest"] != runtime_record):
            raise ValueError("Missing pinned preparation runtime check")
        for item in report[stage]["mapped_files"]:
            check(item)


def validate_checkpoint(report, directory, executor):
    if (report["status"] != "bpo_checkpoint_content_and_indexes_verified"
            or any(report[key] is not False for key in ("search_admitted", "accuracy_admitted", "publication_ready"))
            or report["environment"] != environment()):
        raise ValueError("Wrong checkpoint status, semantics or environment")
    paths = [r["path"] for r in report["outputs"]]
    if len(paths) != len(set(paths)) or set(paths) != {str(directory / name) for name in OUTPUT_NAMES}:
        raise ValueError("Wrong native checkpoint output inventory")
    prefix = [str(TOOL / "venv_orthomcl/bin/perl"), "-I" + str(TOOL),
              str(executor / "benchmark_tools/run_orthomcl_perl_script.pl")]
    bpo, indexes = directory / "all.bpo", directory / "indexes"
    expected = [
        {"stage": "build", "cwd": str(directory), "argv": prefix + [str(executor / "benchmark_tools/build_orthomcl_bpo_indexes.pl"), str(bpo), str(indexes)]},
        {"stage": "validate", "cwd": str(directory), "argv": prefix + [str(executor / "benchmark_tools/validate_orthomcl_bpo_indexes.pl"), str(bpo), str(indexes / "all_bpo.idx"), str(indexes / "all_bpo.se")]}]
    if report["commands"] != expected:
        raise ValueError("Changed native checkpoint commands")
    for item in [*report["checked_records"], *report["outputs"]]:
        check(item)
    if any((directory / name).stat().st_size for name in ("build.stderr", "validate.stderr")):
        raise ValueError("Unreviewed native diagnostics")


def recheck_content(blast, fasta, directory, output, expected_content, expected_indexes):
    output.mkdir(parents=True, exist_ok=False)
    content = audit(blast, fasta, directory / "all.bpo", output / "content.json")
    if content["content"] != expected_content:
        raise ValueError("Independent BPO content counts differ")
    helpers = Path(__file__).resolve().parent
    argv = [str(TOOL / "venv_orthomcl/bin/perl"), str(helpers / "run_orthomcl_perl_script.pl"),
            str(helpers / "validate_orthomcl_bpo_indexes.pl"), str(directory / "all.bpo"),
            str(directory / "indexes/all_bpo.idx"), str(directory / "indexes/all_bpo.se")]
    native_step(argv, output, environment(), "validate")
    indexes = json.loads((output / "validate.stdout").read_text())
    if indexes != expected_indexes:
        raise ValueError("Independent native index summary differs")
    return {"content": content["content"], "index_validation": indexes, "index_command": argv,
            "outputs": [record(p) for p in sorted(output.iterdir()) if p.is_file()]}


def admit(root, job, output):
    if output.exists() or output.is_symlink():
        raise FileExistsError(output)
    runtime_before = verify_runtime(root)
    accounting = subprocess.check_output(["sacct", "-j", str(job), "--parsable2",
        "--format=JobIDRaw,State,ExitCode,Elapsed,NodeList,AllocCPUS,ReqMem"], text=True)
    scheduler = require_completed_job(accounting, job)
    executor = root / "benchmarks/work/publication_qfo_corrected_bpo_v1"
    if subprocess.check_output(["git", "-C", str(executor), "rev-parse", "HEAD"], text=True).strip() != EXECUTOR:
        raise ValueError("Changed frozen preparation executor")
    subprocess.run(["git", "-C", str(executor), "diff", "--exit-code", "HEAD", "--", "benchmark_tools", "orthohmm"], check=True)
    base = root / "benchmarks/results/qfo_corrected_orthomcl_v1/bpo_preparation"
    preparation_record, checkpoint_record = record(base / "report.json"), record(base / "checkpoint/report.json")
    preparation = json.loads((base / "report.json").read_text())
    checkpoint = json.loads((base / "checkpoint/report.json").read_text())
    validate_preparation(preparation, scheduler, record(executor / "benchmark_tools/prepare_qfo_corrected_bpo.py"),
                         checkpoint_record, runtime_before["manifest"])
    validate_checkpoint(checkpoint, base / "checkpoint", executor)
    admission_path = root / "benchmarks/work/qfo_corrected_blast_admission_20260918/report.json"
    records = [r for r in preparation["checked_records"] if r["path"] == str(admission_path)]
    if len(records) != 1:
        raise ValueError("Missing unique source search admission record")
    admission, inputs, search_checked, admission_scheduler, _ = verify_admission(
        root, admission_path, records[0]["sha256"], int(preparation["admission_scheduler"]["JobIDRaw"]))
    if (admission_scheduler != preparation["admission_scheduler"]
            or admission["scheduler"] != preparation["admitted_search_scheduler"]
            or admission["query_coverage"] != preparation["query_coverage"]
            or preparation["content"] != checkpoint["content"]
            or preparation["index_validation"] != checkpoint["index_validation"]
            or checkpoint["content"]["input_proteins"] != 984137):
        raise ValueError("Changed admitted scope or preparation summaries")
    manifests = [read_frozen(root / "benchmark_tools/results" / name, sha) for name, sha in (
        ("qfo_corrected_orthomcl_perl_runtime_20260918.json", RUNTIME_SHA),
        ("qfo_corrected_orthomcl_system_helpers_20260918.json", HELPERS_SHA))]
    for manifest in manifests:
        verify(manifest)
    checked = [preparation_record, checkpoint_record, *preparation["checked_records"], *search_checked,
               *checkpoint["checked_records"], *checkpoint["outputs"], record(__file__),
               *[record(Path(__file__).with_name(name)) for name in (
                   "audit_orthomcl_bpo_content.py", "validate_orthomcl_bpo_indexes.pl", "run_orthomcl_perl_script.pl")]]
    for item in checked:
        check(item)
    output.mkdir(parents=True, exist_ok=False)
    result = {"status": "validating", "source": record(__file__), "scheduler": scheduler,
              "accounting": accounting, "checked_records": checked, "runtime_before": runtime_before,
              "query_coverage": admission["query_coverage"], "accuracy_admitted": False, "publication_ready": False}
    try:
        validation = recheck_content(Path(inputs[0]["path"]), Path(inputs[1]["path"]), base / "checkpoint",
                                     output / "recheck", checkpoint["content"], checkpoint["index_validation"])
        for manifest in manifests:
            verify(manifest)
        for item in [*checked, *validation["outputs"]]:
            check(item)
        result.update(status="corrected_orthomcl_bpo_checkpoint_admitted", validation=validation,
                      runtime_after=verify_runtime(root), native_inputs=[record(base / "checkpoint" / name)
                          for name in ("all.bpo", "indexes/all_bpo.idx", "indexes/all_bpo.se")])
        result["limitations"] = ["Complete content/index recheck, not an independent biological truth assessment.",
                                 "Source failed-query diagnostics are retained; conversion does not repair absent hits.",
                                 "Native inference, final-group validation and scoring remain required."]
    except BaseException as error:
        result.update(status="failed", error_type=type(error).__name__, error=str(error))
        raise
    finally:
        with (output / "report.json").open("x") as stream:
            json.dump(result, stream, indent=2, sort_keys=True)
            stream.write("\n")
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("root", "output"):
        parser.add_argument("--" + name, type=Path, required=True)
    parser.add_argument("--job", type=int, required=True)
    args = parser.parse_args()
    print(json.dumps({"status": admit(args.root.resolve(), args.job, args.output.resolve())["status"]}))
