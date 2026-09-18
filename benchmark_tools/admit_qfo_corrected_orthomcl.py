"""Independently admit completed corrected native OrthoMCL outputs for conversion."""

import argparse
import json
import math
from pathlib import Path
import subprocess
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.audit_orthomcl_native_groups import audit
from benchmark_tools.prepare_ob_candidate_neighborhood import record, check
from benchmark_tools.prepare_orthomcl_bpo_checkpoint import native_step
from benchmark_tools.prepare_qfo_corrected_orthomcl_native import TOOL, RUNTIME_SHA, HELPERS_SHA
from benchmark_tools.run_qfo_corrected_blast import environment
from benchmark_tools.run_qfo_corrected_orthomcl import admitted_inputs, pair_cache_records, ADMITTER, SOURCES_SHA
from benchmark_tools.run_simulation_methods import read_frozen
from benchmark_tools.snapshot_runtime_trees import verify
from benchmark_tools.stage_orthomcl_native_inputs import NAMES
from benchmark_tools.verify_orthomcl_python_runtime import verify_runtime
from benchmark_tools.verify_ygob_validation import require_completed_job

EXECUTOR = "22f4eec4314304a0cda71e02595adfb20827c34b"


def completed(job, cpus, memory):
    accounting = subprocess.check_output(["sacct", "-j", str(job), "--parsable2",
        "--format=JobIDRaw,State,ExitCode,Elapsed,NodeList,AllocCPUS,ReqMem"], text=True)
    scheduler = require_completed_job(accounting, job)
    if any(scheduler.get(k) != v for k, v in {"NodeList": "bizon", "AllocCPUS": str(cpus), "ReqMem": memory}.items()):
        raise ValueError("Wrong native/admission scheduler allocation")
    return scheduler, accounting


def frozen(directory, revision):
    if subprocess.check_output(["git", "-C", str(directory), "rev-parse", "HEAD"], text=True).strip() != revision:
        raise ValueError("Changed frozen executor")
    subprocess.run(["git", "-C", str(directory), "diff", "--exit-code", "HEAD", "--", "benchmark_tools", "orthohmm"], check=True)


def unique_records(records):
    result = {}
    for item in records:
        if item["path"] in result and item != result[item["path"]]:
            raise ValueError("Conflicting provenance records")
        result[item["path"]] = item
    return list(result.values())


def validate_execution(report, scheduler, executor, base, runtime):
    execution, tool = base / "inference_execution", base / "native_tool"
    expected_scheduler = {"State": "COMPLETED", "ExitCode": "0:0", "NodeList": "bizon", "AllocCPUS": "180", "ReqMem": "900G"}
    if (any(scheduler.get(k) != v for k, v in expected_scheduler.items())
            or report["status"] != "corrected_orthomcl_native_exited_zero_pending_admission"
            or report["job_id"] != scheduler["JobIDRaw"] or report["node"] != "bizon"
            or type(report["native_exit_code"]) is not int or report["native_exit_code"] != 0
            or report["accuracy_admitted"] is not False or report["publication_ready"] is not False):
        raise ValueError("Invalid native completion evidence")
    if report["source"] != record(executor / "benchmark_tools/run_qfo_corrected_orthomcl.py"):
        raise ValueError("Changed native runner source")
    expected = [str(TOOL / "venv_orthomcl/bin/perl"), "-I" + str(tool),
                str(executor / "benchmark_tools/run_orthomcl_perl_script.pl"), str(tool / "orthomcl.pl"),
                "--mode", "4", "--bpo_file", str(execution / "inputs/all.bpo"),
                "--gg_file", str(execution / "inputs/all.gg")]
    if (report["command"] != expected or report["cwd"] != str(execution)
            or report["environment"] != {**environment(), "ORTHOMCL_PAIR_WORKERS": "64"}):
        raise ValueError("Changed native command/environment")
    times = [report[k] for k in ("started_epoch", "native_started_epoch", "native_finished_epoch", "finished_epoch")]
    if any(type(v) not in (int, float) or not math.isfinite(v) or v <= 0 for v in times) or times != sorted(times):
        raise ValueError("Invalid native execution timestamps")
    for key in ("runtime_before", "runtime_after"):
        if report[key]["status"] != "dedicated_bpo_python_runtime_verified" or report[key]["manifest"] != runtime:
            raise ValueError("Missing native Python runtime identity")
        for item in report[key]["mapped_files"]:
            check(item)
    for key, name in (("native_log", "native.log"), ("native_timing", "native.time.txt"), ("staging", "inputs/staging.json")):
        if report[key] != record(execution / name):
            raise ValueError("Changed native execution artifact: " + key)


def validate_staging(staging, report, inputs, indexes, directory, executor):
    if (staging["status"] != "native_inputs_staged_and_indexes_verified"
            or staging["accuracy_admitted"] is not False or staging["publication_ready"] is not False
            or staging["sources"] != inputs or staging["staged"] != report["staged_inputs"]
            or staging["input_proteins"] != 984137 or staging["species"] != 78
            or staging["index_validation"] != indexes or set(staging["staged"]) != set(NAMES)
            or set(report["input_mtimes_ns"]) != set(NAMES)):
        raise ValueError("Changed staged input scope or identity")
    for key, name in NAMES.items():
        path = directory / name
        item = staging["staged"][key]
        if (item["path"] != str(path) or path.is_symlink() or path.stat().st_nlink != 1
                or any(item[field] != inputs[key][field] for field in ("bytes", "sha256"))
                or path.stat().st_mtime_ns != report["input_mtimes_ns"][key]):
            raise ValueError("Changed/shared staged native input")
        check(item)
    argv = [str(TOOL / "venv_orthomcl/bin/perl"), str(executor / "benchmark_tools/run_orthomcl_perl_script.pl"),
            str(executor / "benchmark_tools/validate_orthomcl_bpo_indexes.pl"), str(directory / "all.bpo"),
            str(directory / "all_bpo.idx"), str(directory / "all_bpo.se")]
    if staging["index_validation_command"] != argv:
        raise ValueError("Changed staged index validation command")
    if ((directory / "validation/indexes.stderr").stat().st_size
            or json.loads((directory / "validation/indexes.stdout").read_text()) != indexes):
        raise ValueError("Changed staged index validation result")


def validate_outputs(report, base):
    tool = base / "native_tool"
    groups = list(tool.glob("*/all_orthomcl.out"))
    if len(groups) != 1 or report["native_groups"] != record(groups[0]):
        raise ValueError("Missing/changed unique native final groups")
    native = groups[0].parent
    paths = list(native.rglob("*"))
    if native.is_symlink() or any(p.is_symlink() for p in paths):
        raise ValueError("Symlink in native output inventory")
    observed = [record(p) for p in sorted(paths) if p.is_file()]
    if observed != report["native_outputs"]:
        raise ValueError("Changed native output inventory")
    for name in ("tmp/all_ortho.idx", "tmp/all_ortho.mtx", "tmp/all_ortho.mcl"):
        if not (native / name).is_file() or (native / name).stat().st_size == 0:
            raise ValueError("Missing native graph/partition/index")
    if pair_cache_records(base / "inference_execution/inputs") != report["pair_caches"]:
        raise ValueError("Changed native pair caches")
    return native


def admit(root, job, output):
    if output.exists() or output.is_symlink():
        raise FileExistsError(output)
    runtime = verify_runtime(root)
    scheduler, accounting = completed(job, 180, "900G")
    executor = root / "benchmarks/work/publication_qfo_corrected_orthomcl_v1"
    frozen(executor, EXECUTOR)
    base = root / "benchmarks/results/qfo_corrected_orthomcl_v1"
    report_path = base / "inference_execution/status.json"
    report_record = record(report_path)
    report = json.loads(report_path.read_text())
    validate_execution(report, scheduler, executor, base, runtime["manifest"])
    admission_path = root / "benchmarks/work/qfo_corrected_bpo_admission_20260918/report.json"
    admission_record = [r for r in unique_records(report["checked_records"]) if r["path"] == str(admission_path)]
    if len(admission_record) != 1:
        raise ValueError("Missing source BPO admission identity")
    admission = read_frozen(admission_path, admission_record[0]["sha256"])
    inputs = admitted_inputs(admission, root)
    prior_scheduler, _ = completed(int(report["admission_scheduler"]["JobIDRaw"]), 2, "64G")
    prior_executor = root / "benchmarks/work/publication_qfo_corrected_bpo_admission_v1"
    frozen(prior_executor, ADMITTER)
    if (prior_scheduler != report["admission_scheduler"]
            or admission["source"] != record(prior_executor / "benchmark_tools/admit_qfo_corrected_bpo.py")
            or admission["query_coverage"] != report["query_coverage"]):
        raise ValueError("Changed source BPO admission or query diagnostics")
    staging = json.loads(Path(report["staging"]["path"]).read_text())
    validate_staging(staging, report, inputs, admission["validation"]["index_validation"],
                     base / "inference_execution/inputs", executor)
    results = root / "benchmark_tools/results"
    source_path = results / "qfo_corrected_orthomcl_native_sources_20260918.json"
    sources = read_frozen(source_path, SOURCES_SHA)
    manifests = [read_frozen(results / name, sha) for name, sha in (
        ("qfo_corrected_orthomcl_perl_runtime_20260918.json", RUNTIME_SHA),
        ("qfo_corrected_orthomcl_system_helpers_20260918.json", HELPERS_SHA))]
    for manifest in manifests:
        verify(manifest)
    checked = unique_records([report_record, *report["checked_records"], report["staging"],
        report["native_log"], report["native_timing"], *admission["checked_records"],
        *admission["validation"]["outputs"], *inputs.values(), *report["staged_inputs"].values(),
        record(source_path), *sources["checked_records"], *sources["originals"], *sources["configured_sources"],
        *[record(Path(__file__).with_name(name)) for name in (
            "admit_qfo_corrected_orthomcl.py", "audit_orthomcl_native_groups.py",
            "run_qfo_corrected_orthomcl.py", "validate_orthomcl_bpo_indexes.pl", "run_orthomcl_perl_script.pl")]])
    for item in checked:
        check(item)
    native = validate_outputs(report, base)
    output.mkdir(parents=True, exist_ok=False)
    result = {"status": "validating", "source": record(__file__), "scheduler": scheduler,
              "accounting": accounting, "checked_records": checked, "runtime_before": runtime,
              "query_coverage": report["query_coverage"], "accuracy_admitted": False, "publication_ready": False}
    try:
        validation = audit(native / "all_orthomcl.out", native / "tmp/all_ortho.mcl", native / "tmp/all_ortho.idx",
                           base / "inference_execution/inputs/all.gg", output / "groups.json")
        content = validation["content"]
        if (content["input_proteins"] != 984137 or content["input_species"] != 78
                or report["content"] != {"groups": content["final_groups"], "grouped_proteins": content["grouped_proteins"],
                                         "ungrouped_proteins": content["ungrouped_input_proteins"]}):
            raise ValueError("Independent final-group counts differ")
        helpers, directory = Path(__file__).resolve().parent, base / "inference_execution/inputs"
        argv = [str(TOOL / "venv_orthomcl/bin/perl"), str(helpers / "run_orthomcl_perl_script.pl"),
                str(helpers / "validate_orthomcl_bpo_indexes.pl"), *[str(directory / NAMES[k]) for k in ("bpo", "offsets", "ranges")]]
        native_step(argv, output, environment(), "indexes")
        if json.loads((output / "indexes.stdout").read_text()) != admission["validation"]["index_validation"]:
            raise ValueError("Independent staged index check differs")
        for item in checked:
            check(item)
        validate_staging(staging, report, inputs, admission["validation"]["index_validation"], directory, executor)
        validate_outputs(report, base)
        for manifest in manifests:
            verify(manifest)
        result.update(status="corrected_orthomcl_native_outputs_admitted", content=content,
                      runtime_after=verify_runtime(root), native_groups=report["native_groups"],
                      pair_semantics="cross_species_final_group_cliques", index_validation_command=argv,
                      outputs=[record(output / name) for name in ("groups.json", "indexes.stdout", "indexes.stderr")])
        result["limitations"] = ["Native final-output admission authorizes conversion, not accuracy or publication claims.",
                                  "Failed-query diagnostics remain retained; inference does not repair missing hits.",
                                  "Shared-host resource observations are not matched end-to-end timings."]
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
    print(json.dumps({"status": admit(args.root.resolve(), args.job, args.output.absolute())["status"]}))
