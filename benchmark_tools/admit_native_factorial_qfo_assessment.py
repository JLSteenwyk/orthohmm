"""Independently admit six native QfO endpoints with retained failure receipts."""

import argparse
import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.admit_qfo_corrected_comparator_assessment import accounting, validate_completion
from benchmark_tools.admit_qfo_recovered_assessment import validate_trace
from benchmark_tools.audit_qfo_fas_samples import read_sample, summarize
from benchmark_tools.prepare_ob_candidate_neighborhood import check, record
from benchmark_tools.probe_dgx_step_separation import save
from benchmark_tools.run_native_factorial_qfo_assessment import conversion_binding, execution_spec
from benchmark_tools.run_native_factorial_cost import read
from benchmark_tools.validate_native_factorial_outputs import require
from benchmark_tools.validate_qfo_native_assessment import validate_directory


def fas_sample(results, stage, assessment):
    paths = list((results / "results/FAS").glob("*raw.txt.gz"))
    require(len(paths) == 1, "Require exactly one native raw FAS sample")
    raw_ref = record(paths[0])
    sample = read_sample(paths[0])
    missing = set(sample)
    with Path(stage["filtered_pairs"]["path"]).open() as stream:
        for line in stream:
            missing.discard(tuple(line.rstrip("\n").split("\t")))
    require(not missing, "Raw FAS sample contains pairs absent from submitted predictions")
    summary = summarize(sample, assessment["endpoints"]["FAS"]["native_participant"])
    require(summary["reported_eligible_pairs"] <= stage["retained_pairs"], "FAS eligible count exceeds submissions")
    check(raw_ref)
    return dict(raw=raw_ref, sample_membership_verified=True, **summary)


def admit(root, pairs_ref, conversion_job, assessment_job, destination):
    destination = Path(destination)
    require(destination.is_absolute() and destination.resolve() == destination
        and destination.parent == root / "benchmarks/results/full_native_qfo_admission_v1",
        "Require separate direct admission destination")
    if destination.exists() or destination.is_symlink():
        raise FileExistsError(destination)
    destination.mkdir(parents=True, exist_ok=False)
    report = dict(schema="full_native_factorial_qfo_admission_v1", status="validating",
        source=record(__file__), pairs_manifest=pairs_ref, accuracy_admitted=False,
        publication_ready=False, next_identity_authorized=False, automatic_retry=False)
    try:
        text, scheduler = accounting(assessment_job)
        stage, manifest, binding = conversion_binding(root, pairs_ref, conversion_job)
        spec = execution_spec(root, pairs_ref, stage, manifest, binding)
        directory, results = Path(spec["cwd"]), Path(spec["results"])
        require(not destination.is_relative_to(directory) and not directory.is_relative_to(destination),
            "Admission overlaps native execution directory")
        execution_ref, preflight_ref = record(directory / "results.json"), record(directory / "preflight.json")
        execution, preflight = read(execution_ref), read(preflight_ref)
        validate_completion(execution, preflight, scheduler)
        # A fresh native scheduler observation may have a different query timestamp.
        for key, value in spec.items():
            if key not in {"status", "native_scheduler"}:
                require(execution.get(key) == value, "Changed native QfO execution provenance: " + key)
        require(execution["log"]["path"] == str(directory / "scoring.log"), "Assessment log path differs")
        require(type(execution.get("started_monotonic_ns")) is int and type(execution.get("finished_monotonic_ns")) is int
            and execution["finished_monotonic_ns"] >= execution["started_monotonic_ns"], "Invalid assessment interval")
        observed = [record(p) for p in sorted(results.rglob("*")) if p.is_file()]
        require(observed == execution["outputs"], "Native assessment output inventory changed")
        traces = list((results / "stats").glob("trace_*.txt"))
        require(len(traces) == 1, "Require exactly one fresh native scoring trace")
        tasks = validate_trace(traces[0].read_text())
        assessment, paths = validate_directory(results, stage["participant"], Path(manifest["pipeline"]) / "reference_data")
        inventoried = {r["path"] for r in observed}
        require(all(str(p.resolve()) in inventoried for p in paths if p.is_relative_to(results)),
            "Native endpoint absent from admitted output inventory")
        sample = fas_sample(results, stage, assessment)
        validators = [record(Path(__file__).with_name(n)) for n in (
            "admit_qfo_recovered_assessment.py", "audit_qfo_fas_samples.py", "validate_qfo_native_assessment.py",
            "qfo_summarize_scores.py", "run_native_factorial_qfo_assessment.py")]
        checked = [execution_ref, preflight_ref, execution["source"], execution["log"],
            *binding["verified_records"], *observed, *[record(p) for p in paths], *validators, report["source"]]
        for item in checked:
            check(item)
        report.update(status="full_native_factorial_qfo_assessment_admitted", native_index=stage["native_index"],
            native_job_id=stage["native_job_id"], cell=stage["cell"], participant=stage["participant"],
            scheduler=scheduler, accounting=text, conversion=stage, conversion_scheduler=binding["conversion_scheduler"],
            conversion_accounting=binding["conversion_accounting"], native_scheduler=binding["native_scheduler"],
            execution_report=execution_ref, preflight=preflight_ref, environment_manifest=binding["environment_manifest"],
            assessment=assessment, native_trace=record(traces[0]), native_tasks=tasks,
            fas_protocol=binding["fas_protocol"], fas_sample=sample, checked_records=checked, accuracy_admitted=True,
            limitations=["Six native endpoints on development-exposed QfO; not independent biological validation.",
                "Only VGNC, SwissTrees and TreeFam-A summaries are F1; GO/EC similarity and FAS are not F1.",
                "The six-endpoint mean is a project-defined secondary summary, not official QfO F1.",
                "Native SEMs do not establish independent-family or paired-method confidence intervals.",
                "FAS eligible population is not all predictions; sample membership/arithmetic does not prove annotation completeness or sample representativeness.",
                "Frozen manifest files are rechecked; this is not a new hermetic host-runtime closure certification.",
                "Empty predictions and failed assessments are retained; missing scores are not imputed.",
                "Scoring interval is separate from inference/conversion and is not an isolated-performance estimate."])
    except BaseException as error:
        report.update(status="full_native_factorial_qfo_admission_failed_retained", error_type=type(error).__name__, error=str(error))
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
    require(pairs_ref["sha256"] == args.pairs_sha256, "Pair preparation checksum differs")
    result = admit(args.root.resolve(), pairs_ref, args.conversion_job, args.assessment_job, args.output_directory.absolute())
    print(json.dumps(dict(status=result["status"], cell=result["cell"]), sort_keys=True))
