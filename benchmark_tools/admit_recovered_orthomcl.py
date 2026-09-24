"""Independently validate recovered native groups before pair conversion."""

import argparse
import json
import math
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.admit_blast_recovery_bpo import execution_identity
from benchmark_tools.admit_qfo_corrected_orthomcl import completed, frozen, unique_records, validate_staging
from benchmark_tools.audit_orthomcl_native_groups import audit
from benchmark_tools.configure_orthomcl_1_4 import replace_assignment
from benchmark_tools.parallelize_orthomcl_pairs import parallelize_source
from benchmark_tools.prepare_ob_candidate_neighborhood import check, record
from benchmark_tools.prepare_orthomcl_bpo_checkpoint import native_step
from benchmark_tools.prepare_qfo_corrected_orthomcl_native import TOOL, RUNTIME_SHA, HELPERS_SHA
from benchmark_tools.recovered_orthomcl_inputs import verify_inputs
from benchmark_tools.run_blast_recovery_batch import save_status
from benchmark_tools.run_qfo_corrected_blast import environment
from benchmark_tools.run_qfo_corrected_orthomcl import pair_cache_records, SOURCES_SHA
from benchmark_tools.run_simulation_methods import read_frozen
from benchmark_tools.snapshot_runtime_trees import verify
from benchmark_tools.stage_orthomcl_native_inputs import NAMES
from benchmark_tools.verify_orthomcl_python_runtime import verify_runtime

RUNNER_SHA = "90ff4e5ba00393e82f73a14576de21ede003ad5614fe3457e70f3a794597bf3e"


def validate_execution(report, scheduler, executor, base, runtime):
    expected = dict(State="COMPLETED", ExitCode="0:0", NodeList="bizon", AllocCPUS="180", ReqMem="900G")
    source = record(executor / "benchmark_tools/run_recovered_orthomcl.py")
    if (any(scheduler.get(k) != v for k, v in expected.items())
            or report["status"] != "recovered_native_exited_zero_pending_admission"
            or report["job_id"] != scheduler["JobIDRaw"] or report["node"] != "bizon"
            or type(report["native_exit_code"]) is not int or report["native_exit_code"] != 0
            or source["sha256"] != RUNNER_SHA or report["source"] != source
            or any(report[k] is not False for k in (
                "accuracy_admitted", "publication_ready", "downstream_execution_authorized"))):
        raise ValueError("Invalid recovered native completion evidence")
    tool, inputs = base / "native_tool", base / "inputs"
    expected_command = [str(TOOL / "venv_orthomcl/bin/perl"), "-I" + str(tool),
        str(executor / "benchmark_tools/run_orthomcl_perl_script.pl"), str(tool / "orthomcl.pl"),
        "--mode", "4", "--bpo_file", str(inputs / "all.bpo"), "--gg_file", str(inputs / "all.gg")]
    if (report["command"] != expected_command or report["cwd"] != str(base)
            or report["environment"] != {**environment(), "ORTHOMCL_PAIR_WORKERS": "64"}):
        raise ValueError("Changed recovered native command/environment")
    times = [report[k] for k in ("started_epoch", "native_started_epoch", "native_finished_epoch", "finished_epoch")]
    if any(type(v) not in (int, float) or not math.isfinite(v) or v <= 0 for v in times) or times != sorted(times):
        raise ValueError("Invalid recovered native timestamps")
    for key in ("runtime_before", "runtime_after"):
        if report[key]["status"] != "dedicated_bpo_python_runtime_verified" or report[key]["manifest"] != runtime:
            raise ValueError("Changed recovered native Python identity")
        for item in report[key]["mapped_files"]:
            check(item)
    for key, name in (("native_log", "native.log"), ("native_timing", "native.time.txt"), ("staging", "inputs/staging.json")):
        if report[key] != record(base / name):
            raise ValueError("Changed recovered native artifact")


def validate_sources(report, originals, base):
    sources = report["native_sources"]
    tool, inputs = base / "native_tool", base / "inputs"
    if (sources["originals"] != originals or sources["threads"] != 180
            or sources["pair_parallel_patch"] is not True or sources["tool_directory"] != str(tool)
            or sources["data_directory"] != str(inputs)):
        raise ValueError("Changed recovered native source configuration")
    expected = [tool / name for name in ("orthomcl.pl", "orthomcl_module.pm")]
    if sources["configured_sources"] != [record(p) for p in expected]:
        raise ValueError("Changed recovered configured source records")
    for original, destination in zip(originals, expected):
        check(original)
        text = Path(original["path"]).read_text()
        if destination.name == "orthomcl.pl":
            text = parallelize_source(text)
        else:
            for name, value in (("PATH_TO_ORTHOMCL", f'"{tool}/"'),
                                ("ORTHOMCL_DATA_DIR", f'"{inputs}/"'), ("BLAST_NOCPU", "180")):
                text = replace_assignment(text, name, value)
        if destination.is_symlink() or destination.read_text() != text:
            raise ValueError("Native source differs from pinned source plus approved configuration")


def validate_outputs(report, base):
    groups = list((base / "native_tool").glob("*/all_orthomcl.out"))
    if len(groups) != 1 or report["native_groups"] != record(groups[0]):
        raise ValueError("Missing/changed recovered final groups")
    native = groups[0].parent
    paths = list(native.rglob("*"))
    if native.is_symlink() or any(p.is_symlink() for p in paths):
        raise ValueError("Symlink in recovered native output inventory")
    if [record(p) for p in sorted(paths) if p.is_file()] != report["native_outputs"]:
        raise ValueError("Changed recovered output inventory")
    for name in ("tmp/all_ortho.idx", "tmp/all_ortho.mtx", "tmp/all_ortho.mcl"):
        if not (native / name).is_file() or (native / name).stat().st_size == 0:
            raise ValueError("Missing native graph/partition/index")
    if pair_cache_records(base / "inputs") != report["pair_caches"]:
        raise ValueError("Changed recovered pair caches")
    return native


def admit(root, job, executor, commit, bpo_executor, bpo_commit, output):
    if output.exists() or output.is_symlink():
        raise FileExistsError(output)
    identity, runtime = execution_identity(), verify_runtime(root)
    scheduler, accounting = completed(job, 180, "900G")
    if not executor.resolve().is_relative_to(root / "benchmarks/work"):
        raise ValueError("Require retained native executor")
    frozen(executor, commit)
    base = root / "benchmarks/results/qfo_blast_recovery_native_v1"
    report_record = record(base / "report.json")
    report = read_frozen(Path(report_record["path"]), report_record["sha256"])
    validate_execution(report, scheduler, executor, base, runtime["manifest"])
    prior = report["admission"]
    evidence = verify_inputs(root, Path(prior["path"]), prior["sha256"], int(report["admission_job"]),
                             bpo_executor, bpo_commit)
    if (evidence["scheduler"] != report["admission_scheduler"]
            or evidence["query_coverage"] != report["query_coverage"]
            or evidence["admission"] != prior or prior not in report["checked_records"]):
        raise ValueError("Changed recovered BPO admission or failure coverage")
    staging = read_frozen(Path(report["staging"]["path"]), report["staging"]["sha256"])
    validate_staging(staging, report, evidence["inputs"], evidence["index_validation"], base / "inputs", executor)
    results = root / "benchmark_tools/results"
    source_path = results / "qfo_corrected_orthomcl_native_sources_20260918.json"
    sources = read_frozen(source_path, SOURCES_SHA)
    validate_sources(report, sources["originals"], base)
    manifests = [read_frozen(results / name, sha) for name, sha in (
        ("qfo_corrected_orthomcl_perl_runtime_20260918.json", RUNTIME_SHA),
        ("qfo_corrected_orthomcl_system_helpers_20260918.json", HELPERS_SHA))]
    for manifest in manifests:
        verify(manifest)
    checked = unique_records([report_record, *report["checked_records"], *evidence["checked_records"],
        report["staging"], report["native_log"], report["native_timing"], *report["staged_inputs"].values(),
        record(source_path), *sources["checked_records"], *sources["originals"],
        *report["native_sources"]["configured_sources"],
        *[record(p) for p in sorted(Path(__file__).parent.glob("*.py"))],
        *[record(Path(__file__).with_name(n)) for n in ("run_orthomcl_perl_script.pl", "validate_orthomcl_bpo_indexes.pl")]])
    for item in checked:
        check(item)
    native = validate_outputs(report, base)
    output.mkdir(parents=True, exist_ok=False)
    result = dict(status="recovered_native_validation_running", **identity, source=record(__file__),
        scheduler=scheduler, accounting=accounting, execution=report_record, checked_records=checked,
        runtime_before=runtime, query_coverage=report["query_coverage"], conversion_authorized=False,
        accuracy_admitted=False, publication_ready=False)
    save_status(output / "report.json", result)
    try:
        validation = audit(native / "all_orthomcl.out", native / "tmp/all_ortho.mcl", native / "tmp/all_ortho.idx",
                           base / "inputs/all.gg", output / "groups.json")
        content = validation["content"]
        if (content["input_proteins"] != 984137 or content["input_species"] != 78
                or report["content"] != dict(groups=content["final_groups"], grouped_proteins=content["grouped_proteins"],
                                              ungrouped_proteins=content["ungrouped_input_proteins"])):
            raise ValueError("Independent recovered final-group counts differ")
        helpers = Path(__file__).resolve().parent
        argv = [str(TOOL / "venv_orthomcl/bin/perl"), str(helpers / "run_orthomcl_perl_script.pl"),
                str(helpers / "validate_orthomcl_bpo_indexes.pl"),
                *[str(base / "inputs" / NAMES[k]) for k in ("bpo", "offsets", "ranges")]]
        native_step(argv, output, environment(), "indexes")
        if json.loads((output / "indexes.stdout").read_text()) != evidence["index_validation"]:
            raise ValueError("Independent recovered staged index check differs")
        for item in checked:
            check(item)
        validate_staging(staging, report, evidence["inputs"], evidence["index_validation"], base / "inputs", executor)
        validate_sources(report, sources["originals"], base)
        validate_outputs(report, base)
        for manifest in manifests:
            verify(manifest)
        result.update(status="recovered_orthomcl_native_outputs_admitted", conversion_authorized=True,
            content=content, runtime_after=verify_runtime(root), native_groups=report["native_groups"],
            pair_semantics="cross_species_final_group_cliques", index_validation_command=argv,
            outputs=[record(output / n) for n in ("groups.json", "indexes.stdout", "indexes.stderr")],
            limitations=["Final-group structural admission, not biological accuracy or publication readiness.",
                         "Failed-query diagnostics remain retained; no missing search hits are repaired.",
                         "Shared-host accounting is not matched end-to-end timing."])
    except BaseException as error:
        result.update(status="recovered_native_validation_failed", error_type=type(error).__name__, error=str(error))
        raise
    finally:
        save_status(output / "report.json", result)
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for flag in ("root", "executor", "bpo-executor", "output"):
        parser.add_argument("--" + flag, type=Path, required=True)
    parser.add_argument("--job", type=int, required=True)
    parser.add_argument("--commit", required=True)
    parser.add_argument("--bpo-commit", required=True)
    args = parser.parse_args()
    result = admit(args.root.resolve(), args.job, args.executor.resolve(), args.commit,
                   args.bpo_executor.resolve(), args.bpo_commit, args.output.absolute())
    print(result["status"])
