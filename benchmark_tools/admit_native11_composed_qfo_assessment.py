"""Independently admit unchanged QfO endpoints for composed native11 provenance."""

import argparse
import json
from pathlib import Path

from benchmark_tools.admit_qfo_corrected_comparator_assessment import accounting, validate_completion
from benchmark_tools.admit_qfo_recovered_assessment import validate_trace
from benchmark_tools.admit_native_factorial_qfo_assessment import fas_sample
from benchmark_tools.prepare_ob_candidate_neighborhood import record, check
from benchmark_tools.probe_dgx_step_separation import save
from benchmark_tools.run_native11_composed_qfo_assessment import conversion_binding, execution_spec
from benchmark_tools.run_native_factorial_cost import ROOT, read
from benchmark_tools.validate_native_factorial_outputs import require
from benchmark_tools.validate_qfo_native_assessment import validate_directory


SCHEMA = "native11_composed_qfo_admission_v1"
DESTINATION = ROOT / "benchmarks/results/native11_composed_qfo_admission_v1"


def admit(root, pairs_ref, conversion_job, assessment_job, destination, expected_source_sha256):
    source = record(__file__)
    require(source["sha256"] == expected_source_sha256, "Prospective admission worker changed")
    destination = Path(destination)
    require(root == ROOT and destination == DESTINATION and destination.is_absolute()
        and destination.resolve() == destination and not destination.exists(),
        "Require separate fresh composed admission destination")
    destination.mkdir(parents=True, exist_ok=False)
    report = dict(schema=SCHEMA, status="validating", source=source, pairs_manifest=pairs_ref,
        accuracy_admitted=False, publication_ready=False, next_identity_authorized=False,
        automatic_retry=False, original_review_translated=False)
    try:
        raw, scheduler = accounting(assessment_job)
        stage, manifest, binding = conversion_binding(root, pairs_ref, conversion_job)
        spec = execution_spec(root, pairs_ref, stage, manifest, binding)
        directory, results = Path(spec["cwd"]), Path(spec["results"])
        require(not directory.is_relative_to(destination) and not destination.is_relative_to(directory),
            "Admission overlaps execution")
        execution_ref, preflight_ref = record(directory / "results.json"), record(directory / "preflight.json")
        execution, preflight = read(execution_ref), read(preflight_ref)
        validate_completion(execution, preflight, scheduler)
        for key, value in spec.items():
            if key not in {"status", "native_scheduler"}:
                require(execution.get(key) == value, "Changed composed QfO execution provenance: " + key)
        require(execution["log"]["path"] == str(directory / "scoring.log"), "Assessment log path differs")
        require(type(execution.get("started_monotonic_ns")) is int
            and type(execution.get("finished_monotonic_ns")) is int
            and execution["finished_monotonic_ns"] >= execution["started_monotonic_ns"], "Invalid assessment interval")
        observed = [record(p) for p in sorted(results.rglob("*")) if p.is_file()]
        require(observed == execution["outputs"], "Assessment output inventory changed")
        traces = list((results / "stats").glob("trace_*.txt"))
        require(len(traces) == 1, "Require exactly one fresh scoring trace")
        tasks = validate_trace(traces[0].read_text())
        assessment, paths = validate_directory(results, stage["participant"], Path(manifest["pipeline"]) / "reference_data")
        inventoried = {r["path"] for r in observed}
        require(all(str(p.resolve()) in inventoried for p in paths if p.is_relative_to(results)),
            "Endpoint absent from output inventory")
        sample = fas_sample(results, stage, assessment)
        validators = [record(ROOT / "benchmark_tools" / name) for name in (
            "admit_qfo_recovered_assessment.py", "audit_qfo_fas_samples.py", "validate_qfo_native_assessment.py",
            "qfo_summarize_scores.py", "admit_native_factorial_qfo_assessment.py", "run_native11_composed_qfo_assessment.py")]
        checked = [execution_ref, preflight_ref, execution["source"], execution["log"],
            *binding["verified_records"], *observed, *[record(p) for p in paths], *validators, source]
        for ref in checked:
            check(ref)
        report.update(status="native11_composed_qfo_assessment_admitted", native_index=11, native_job_id=23985,
            cell=stage["cell"], participant=stage["participant"], amendment=stage["amendment"],
            scheduler=scheduler, accounting=raw, conversion=stage,
            conversion_scheduler=binding["conversion_scheduler"], conversion_accounting=binding["conversion_accounting"],
            native_scheduler=binding["native_scheduler"], composed_binding=binding["composed_binding"],
            execution_report=execution_ref, preflight=preflight_ref, environment_manifest=binding["environment_manifest"],
            assessment=assessment, native_trace=record(traces[0]), native_tasks=tasks, fas_protocol=binding["fas_protocol"],
            fas_sample=sample, checked_records=checked, accuracy_admitted=True,
            limitations=["Six unchanged endpoints on development-exposed QfO, not independent test validation.",
                "VGNC/SwissTrees/TreeFam-A F1 differs from GO/EC similarity and FAS; do not relabel the latter F1.",
                "Six-metric mean is a project-defined secondary summary, not official QfO F1.",
                "Native SEMs and FAS sampling do not establish paired independent-family uncertainty or representativeness.",
                "Explicit new composed lineage; original failures retained and ordinary review unchanged.",
                "Scoring/conversion separate from shared-host inference; no isolated efficiency claim.",
                "No missing-score imputation, inference repeat, history authorization or publication readiness."])
    except Exception as error:
        report.update(status="native11_composed_qfo_admission_failed_retained",
            error_type=type(error).__name__, error=str(error))
        raise
    finally:
        save(destination / "results.json", report)
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pairs", type=Path, required=True)
    parser.add_argument("--pairs-sha256", required=True)
    parser.add_argument("--conversion-job", required=True)
    parser.add_argument("--assessment-job", required=True)
    parser.add_argument("--source-sha256", required=True)
    args = parser.parse_args()
    ref = record(args.pairs)
    require(ref["sha256"] == args.pairs_sha256, "Conversion checksum differs")
    result = admit(ROOT, ref, args.conversion_job, args.assessment_job, DESTINATION, args.source_sha256)
    print(json.dumps(dict(status=result["status"], accuracy_admitted=result["accuracy_admitted"]), sort_keys=True))
