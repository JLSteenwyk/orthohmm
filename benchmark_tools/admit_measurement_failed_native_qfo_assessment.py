"""Independently validate recovered native QfO scores without admitting failed timing."""

import argparse
import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.admit_native_factorial_qfo_assessment import fas_sample
from benchmark_tools.admit_qfo_corrected_comparator_assessment import accounting, validate_completion
from benchmark_tools.admit_qfo_recovered_assessment import validate_trace
from benchmark_tools.prepare_ob_candidate_neighborhood import check, record
from benchmark_tools.probe_dgx_step_separation import save
from benchmark_tools.run_measurement_failed_native_qfo_assessment import conversion_binding, execution_spec
from benchmark_tools.run_native_factorial_cost import read
from benchmark_tools.validate_native_factorial_outputs import require
from benchmark_tools.validate_qfo_native_assessment import validate_directory


def admit(root, pairs_ref, conversion_job, assessment_job, destination):
    destination = Path(destination)
    require(destination.is_absolute() and destination.resolve() == destination
        and destination.parent == root / "benchmarks/results/measurement_failed_native_qfo_admission_v1",
        "Require separate direct recovered admission destination")
    if destination.exists() or destination.is_symlink():
        raise FileExistsError(destination)
    destination.mkdir(parents=True, exist_ok=False)
    report = dict(schema="measurement_failed_native_qfo_admission_v1", status="validating",
        source=record(__file__), pairs_manifest=pairs_ref, accuracy_admitted=False, publication_ready=False,
        next_identity_authorized=False, automatic_retry=False, original_native_scheduler_success=False,
        scientific_timings_admitted=False, eligible_for_timing_comparison=False, resources=None,
        native_inference_reexecuted=False)
    try:
        text, scheduler = accounting(assessment_job)
        stage, manifest, binding = conversion_binding(root, pairs_ref, conversion_job)
        spec = execution_spec(root, pairs_ref, stage, manifest, binding)
        directory, results = Path(spec["cwd"]), Path(spec["results"])
        require(not destination.is_relative_to(directory) and not directory.is_relative_to(destination),
            "Recovered admission overlaps execution directory")
        execution_ref, preflight_ref = record(directory / "results.json"), record(directory / "preflight.json")
        execution, preflight = read(execution_ref), read(preflight_ref)
        validate_completion(execution, preflight, scheduler)
        # Fresh scheduler queries may have different timestamps, not different outcomes.
        for key, value in spec.items():
            if key not in {"status", "native_scheduler"}:
                require(execution.get(key) == value, "Changed recovered execution provenance: " + key)
        require(execution["log"]["path"] == str(directory / "scoring.log"), "Recovered assessment log path differs")
        require(type(execution.get("started_monotonic_ns")) is int and type(execution.get("finished_monotonic_ns")) is int
            and execution["finished_monotonic_ns"] >= execution["started_monotonic_ns"], "Invalid recovered assessment interval")
        observed = [record(p) for p in sorted(results.rglob("*")) if p.is_file()]
        require(observed == execution["outputs"], "Recovered assessment output inventory changed")
        traces = list((results / "stats").glob("trace_*.txt"))
        require(len(traces) == 1, "Require exactly one fresh recovered scoring trace")
        tasks = validate_trace(traces[0].read_text())
        assessment, paths = validate_directory(results, stage["participant"], Path(manifest["pipeline"]) / "reference_data")
        inventoried = {r["path"] for r in observed}
        require(all(str(p.resolve()) in inventoried for p in paths if p.is_relative_to(results)),
            "Recovered endpoint absent from output inventory")
        sample = fas_sample(results, stage, assessment)
        validators = [record(Path(__file__).with_name(name)) for name in (
            "admit_native_factorial_qfo_assessment.py", "admit_qfo_recovered_assessment.py", "audit_qfo_fas_samples.py",
            "validate_qfo_native_assessment.py", "qfo_summarize_scores.py", "run_measurement_failed_native_qfo_assessment.py")]
        checked = [execution_ref, preflight_ref, execution["source"], execution["log"], *binding["verified_records"],
            *observed, *[record(p) for p in paths], *validators, report["source"]]
        for ref in checked:
            check(ref)
        report.update(status="measurement_failed_native_qfo_assessment_admitted", native_index=stage["native_index"],
            native_job_id=stage["native_job_id"], cell=stage["cell"], participant=stage["participant"],
            scheduler=scheduler, accounting=text, conversion=stage, conversion_scheduler=binding["conversion_scheduler"],
            conversion_accounting=binding["conversion_accounting"], native_scheduler=binding["native_scheduler"],
            execution_report=execution_ref, preflight=preflight_ref, environment_manifest=binding["environment_manifest"],
            scientific_recovery=stage["scientific_recovery"], assessment=assessment, native_trace=record(traces[0]),
            native_tasks=tasks, fas_protocol=binding["fas_protocol"], fas_sample=sample, checked_records=checked,
            accuracy_admitted=True, limitations=[
                "Scientific accuracy assessment of recovered outputs; original allocation/timing remains failed and excluded.",
                "Development-exposed QfO, not independent biological validation or new species-tree correctness evidence.",
                "Only VGNC, SwissTrees and TreeFam-A summaries are F1; GO/EC similarity and FAS are not F1.",
                "The six-endpoint mean is a project-defined secondary summary, not official QfO F1.",
                "Native SEMs are not independent-family or paired-method confidence intervals.",
                "FAS eligibility, unseeded sampling and omitted missing scores limit representativeness; no imputation.",
                "Manifest/output checks are not new hermetic or continuous-runtime closure certification.",
                "Scoring is separate from failed inference timing; shared-host effects are unknown and tool-dependent."])
    except BaseException as error:
        report.update(status="measurement_failed_native_qfo_admission_failed_retained",
            error_type=type(error).__name__, error=str(error))
        raise
    finally:
        save(destination / "results.json", report)
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("root", "pairs", "output-directory"):
        parser.add_argument("--" + name, type=Path, required=True)
    for name in ("pairs-sha256", "conversion-job", "assessment-job"):
        parser.add_argument("--" + name, required=True)
    args = parser.parse_args()
    pairs_ref = record(args.pairs)
    require(pairs_ref["sha256"] == args.pairs_sha256, "Recovered pair-manifest checksum differs")
    result = admit(args.root.resolve(), pairs_ref, args.conversion_job, args.assessment_job, args.output_directory.absolute())
    print(json.dumps(dict(status=result["status"], cell=result["cell"]), sort_keys=True))
