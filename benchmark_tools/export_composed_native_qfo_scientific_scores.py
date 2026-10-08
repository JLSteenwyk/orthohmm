"""Extend frozen native scores only with explicit final-identity admission."""

import argparse
from copy import deepcopy
import csv
import json
from pathlib import Path

from benchmark_tools import admit_native12_composed_qfo_assessment as admission
from benchmark_tools import export_allocated_native_qfo_scientific_scores as prior
from benchmark_tools import export_native_qfo_scientific_scores as arithmetic
from benchmark_tools import native12_composed_review_binding as binding
from benchmark_tools import run_native12_composed_qfo_assessment as assessment
from benchmark_tools.prepare_ob_candidate_neighborhood import record, check
from benchmark_tools.probe_dgx_step_separation import save
from benchmark_tools.run_native_factorial_cost import ROOT, read


original = arithmetic.original
SCHEMA = "composed_native_qfo_scientific_reporting_v1"
BASELINE = ROOT / "benchmark_tools/results/native_qfo_scientific_scores_20261007_v3/report.json"
BASELINE_SHA = "7aab1cbd31fb650df42e6e80a14ee0167b35a810f08202768c89c514755211a2"
FAILURE = ROOT / "benchmark_tools/results/native11_qfo_scoring_failure_addendum_20261008_v1/report.json"
FAILURE_SHA = "4f1a5c536af8ce992d46866374b3ab44ee3158077316d69ba2e9f5fb34cee203"


def extract_final(report, run, plan_ref, evidence):
    original.require(report.get("schema") == admission.SCHEMA
        and report.get("status") == "native12_composed_qfo_assessment_admitted"
        and report.get("accuracy_admitted") is True and report.get("native_index") == 12
        and type(report["native_index"]) is int and report.get("native_job_id") == 24036
        and report.get("cell") == run["cell"] == "p1_c1_r1" and run["index"] == 12
        and run.get("dataset") == "qfo_corrected"
        and all(report.get(key) is False for key in ("publication_ready", "next_identity_authorized",
            "automatic_retry", "original_review_translated", "native_inference_reexecuted")),
        "Require independently admitted explicit final native score")
    source = record(admission.__file__)
    original.require(report.get("source") == source, "Final admission source differs")
    evidence.append(source)
    stage = original.bound(report["pairs_manifest"], evidence)
    original.require(report["conversion"] == stage and stage["plan"] == plan_ref
        and stage["input_fastas"] == run["inputs"] and stage["participant"] == report["participant"],
        "Final conversion/plan/inputs differ")
    assessment.validate_stage(stage, report["conversion_scheduler"])
    for key, path in (("source", assessment.converter.__file__), ("binding_source", binding.__file__),
            ("conversion_kernel_source", ROOT / "benchmark_tools/prepare_native_factorial_qfo_pairs.py")):
        ref = record(path)
        original.require(stage[key] == ref, "Final conversion helper differs: " + key)
        evidence.append(ref)
    request = original.bound(stage["request"], evidence)
    original.require(stage["request"]["path"] == str(binding.executor.REQUEST)
        and stage["request"]["sha256"] == binding.REQUEST_SHA
        and request.get("schema") == binding.executor.REQUEST_SCHEMA
        and request.get("plan") == plan_ref and request.get("index") == 12
        and request.get("job_id") == 24036 and request.get("cell") == run["cell"]
        and request.get("amendment") == stage["amendment"]
        and request.get("source") == record(binding.executor.__file__)
        and request.get("original_review_translated") is False,
        "Final request is not the exact new identity/source")
    review = original.bound(stage["terminal_review"], evidence)
    provenance = stage["composed_binding"]
    kind = binding.admit_review(review, stage["request"], request, run,
        provenance["review_producer_scheduler"], provenance["review_producer_job_id"])
    reviewer_source = record(binding.reviewer.__file__)
    original.require(reviewer_source["sha256"] == binding.REVIEWER_SHA
        and review["source"] == reviewer_source and reviewer_source in request["new_sources"]
        and kind == stage["conversion_kind"] == "native", "Final reviewer or native pair semantics differ")
    evidence.extend([reviewer_source, request["source"]])
    outputs = original.bound(review["reviews"]["outputs_or_failure"], evidence)
    ready_ref = record(Path(run["output_root"]) / "measurement/ready.json")
    original.require(outputs.get("schema") == "native12_composed_output_review_v1"
        and outputs.get("source") == reviewer_source and outputs.get("native_outputs_validated") is True
        and outputs.get("accuracy_evaluated") is False and outputs.get("index") == 12
        and outputs.get("job_id") == 24036 and outputs.get("cell") == run["cell"]
        and outputs.get("request") == stage["request"] and outputs.get("plan") == plan_ref
        and outputs.get("amendment") == stage["amendment"]
        and outputs.get("semantic_validator_source") == record(ROOT / "benchmark_tools/validate_native_factorial_outputs.py")
        and outputs.get("allocated_ready") == stage["allocated_ready"] == ready_ref
        and outputs.get("gene_ownership_sha256") == stage["gene_ownership_sha256"]
        and stage["native_input"] in outputs["checked_files"]
        and stage["native_input"]["path"] == str(Path(run["output_root"]) /
            "native/orthohmm_phylogeny/orthohmm_pairwise_orthologs.tsv")
        and type(outputs.get("phylogeny", {}).get("native_pair_rows")) is int
        and outputs["phylogeny"]["native_pair_rows"] == stage["retained_pairs"],
        "Final independently reviewed native prediction/count differs")
    cpu_ids = review["native_cpu_ids"]
    original.require(type(cpu_ids) is list and len(cpu_ids) == 32
        and all(type(cpu) is int and cpu >= 0 for cpu in cpu_ids)
        and cpu_ids == sorted(set(cpu_ids)) and outputs.get("native_cpu_ids") == stage["native_cpu_ids"] == cpu_ids,
        "Final selected native CPU identities differ")
    count = stage["retained_pairs"]
    coverage = stage["pair_coverage"]
    total, covered = coverage["input_accessions"], coverage["accessions_in_any_pair"]
    original.require(type(run["genes"]) is int and total == run["genes"] and covered <= 2 * count
        and (count == 0 or covered >= 2), "Final coverage denominator differs")
    execution = original.bound(report["execution_report"], evidence)
    preflight = original.bound(report["preflight"], evidence)
    original.scheduler(report["scheduler"], execution["job_id"], 8)
    original.require(report["scheduler"].get("ReqMem") in {"128G", "128Gn", "131072M", "131072Mn"},
        "Final scoring memory allocation differs")
    env_ref = stage["environment_manifest"]
    manifest = original.bound(env_ref, evidence)
    original.require(env_ref["sha256"] == assessment.converter.ENV_SHA
        and report["environment_manifest"] == env_ref, "Final frozen QfO environment differs")
    expected = assessment.execution_spec(ROOT, report["pairs_manifest"], stage, manifest, {})
    original.require(all(execution.get(key) == value for key, value in expected.items() if key != "status"),
        "Final exact scoring spec differs")
    original.require(execution.get("status") == "process_succeeded_pending_independent_admission"
        and type(execution.get("exit_code")) is int and execution["exit_code"] == 0
        and preflight.get("status") == "running"
        and all(execution.get(key) == value for key, value in preflight.items() if key != "status")
        and report.get("assessment_resource_limits") == execution["assessment_resource_limits"]
        and report.get("resource_scope") == execution["resource_scope"]
        and report.get("composed_binding") == stage["composed_binding"] == execution["composed_binding"]
        and report.get("fas_protocol") == execution["fas_protocol"]
        and report.get("fas_sample", {}).get("sample_membership_verified") is True,
        "Final scoring preflight/lineage/FAS evidence differs")
    checked = report["checked_records"]
    original.require(all(ref in checked for ref in (report["execution_report"], stage["pairs"],
        stage["filtered_pairs"], stage["native_input"], report["preflight"], execution["log"])),
        "Final reporting input was not checked by independent admission")
    evidence.extend([ready_ref, stage["pairs"], stage["filtered_pairs"], stage["native_input"], execution["source"],
        execution["log"], outputs["semantic_validator_source"], *request["new_sources"]])
    scores, details, mean = arithmetic.endpoints(report)
    resources = review["resources"]
    original.require(set(resources) == set(original.SCOPES), "Final native resource scopes differ")
    for key, value in resources.items():
        original.finite(value, key)
    original.require(type(resources["peak_memory_bytes"]) is int, "Final native peak is not integer bytes")
    return dict(index=12, cell=run["cell"], status="supplied_composed_native_admission", native_job_id=24036,
        conversion_job_id=stage["job_id"], assessment_job_id=execution["job_id"], scores=scores,
        endpoint_details=details, secondary_mean=mean, submitted_pairs=count, input_accessions=total,
        relation_accessions=covered, relation_coverage=covered / total, prediction_semantics=stage["semantics"],
        participant=stage["participant"], fas_protocol=report["fas_protocol"], fas_sample=report["fas_sample"],
        measurement_status="successful_composed_native_terminal_review", accuracy_admitted=True,
        amendment=stage["amendment"], native_cpu_ids=cpu_ids, resources=resources, resource_scopes=original.SCOPES,
        assessment_resource_limits=execution["assessment_resource_limits"], scientific_timings_admitted=False,
        resource_observation_scope="terminal-reviewed shared-host native command; not scoring or isolated algorithm cost")


def collect(final_admission_ref=None):
    baseline_ref, failure_ref = record(BASELINE), record(FAILURE)
    original.require(baseline_ref["sha256"] == BASELINE_SHA and failure_ref["sha256"] == FAILURE_SHA,
        "Frozen baseline or scoring-failure addendum changed")
    baseline, failure = read(baseline_ref), read(failure_ref)
    original.require(baseline.get("schema") == "allocated_native_qfo_scientific_reporting_snapshot_v1"
        and baseline.get("source") == record(prior.__file__) and baseline.get("supplied_admissions") == 4
        and [row["index"] for row in baseline["rows"]] == list(range(6, 13))
        and [row["index"] for row in baseline["rows"] if row["accuracy_admitted"] is True] == [6, 7, 8, 10],
        "Require exact retained four-cell native snapshot")
    original.require(failure.get("schema") == "composed_native_qfo_failure_reporting_v1"
        and failure.get("new_scientific_admission") is False and failure.get("frozen_manuscript_replaced") is False,
        "Require separate retained failure-reporting evidence")
    failed = failure["row"]
    original.require(failed["native_index"] == 11 and failed["cell"] == "p1_c1_r0"
        and failed["scoring_status"] == "OUT_OF_MEMORY" and failed["accuracy_admitted"] is False
        and failed["admitted_endpoint_score_count"] == 0 and failed["secondary_six_metric_mean"] is None,
        "Failure row supplies no admitted score")
    evidence = [baseline_ref, failure_ref, baseline["source"], failure["source"],
        *baseline["outputs"], *failure["outputs"], *failure["references"].values()]
    plan = original.bound(baseline["plan"], evidence)
    rows = deepcopy(baseline["rows"])
    rows[5].update(status="retained_composed_scoring_failure", native_job_id=failed["native_job_id"],
        conversion_job_id=failed["conversion_job_id"], assessment_job_id=failed["assessment_job_id"],
        measurement_status=failed["inference_status"], resources=failed["resources"], resource_scopes=failed["resource_scopes"],
        submitted_pairs=failed["submitted_pairs"], input_accessions=failed["input_accessions"],
        relation_accessions=failed["relation_accessions"], relation_coverage=failed["relation_coverage"],
        prediction_semantics=failed["prediction_semantics"], scoring_status=failed["scoring_status"],
        failure_addendum=failure_ref, scores={endpoint: None for endpoint in original.ENDPOINTS},
        secondary_mean=None, accuracy_admitted=False)
    if final_admission_ref is not None:
        original.require(final_admission_ref["path"] == str(admission.DESTINATION / "results.json"),
            "Final admission namespace differs")
        report = original.bound(final_admission_ref, evidence)
        rows[6] = dict(extract_final(report, plan["runs"][12], baseline["plan"], evidence), admission=final_admission_ref)
    source = record(__file__)
    helpers = [record(module.__file__) for module in (arithmetic, binding, assessment, admission)]
    for ref in [*evidence, source, *helpers]:
        check(ref)
    return dict(schema=SCHEMA, source=source, baseline_snapshot=baseline_ref, failure_addendum=failure_ref,
        plan=baseline["plan"], rows=rows, supplied_admissions=sum(row["accuracy_admitted"] for row in rows),
        final_admission=final_admission_ref, historical_four_cell_rows_preserved=True,
        new_scoring_or_admission=False, publication_ready=False, evidence=evidence, helpers=helpers,
        timing_disclosure=baseline["timing_disclosure"], limitations=[
            "Retains the frozen four admitted-cell rows exactly; failed scoring adds no score or mean.",
            "Final row requires explicit new independent native-pair admission; no original-schema translation.",
            "Direct report/source/arithmetic checks, not raw scientific readmission or new score computation.",
            "Final scoring uses a prospective128GiB allocation; inference, conversion and scoring resource scopes remain separate.",
            "No old uncertainty is attached by aggregate agreement; actual final Swiss-family counts must be audited.",
            "QfO is development-exposed; six-metric mean is secondary and GO/EC/FAS are not F1."])


def export(destination, final_admission_ref=None):
    destination = Path(destination).absolute()
    original.require(destination.is_relative_to(ROOT) and destination.resolve() == destination
        and not destination.exists() and not destination.is_symlink(), "Require one fresh successor export directory")
    report = collect(final_admission_ref)
    headers = ["Index", "Cell", "Accuracy Status", "VGNC F1", "SwissTrees F1", "TreeFam-A F1",
        "GO Similarity", "EC Similarity", "FAS", "Secondary Mean", "Relation Coverage", "Prediction Semantics"]
    values = [[row["index"], row["cell"], row["status"], *[row["scores"][endpoint] for endpoint in original.ENDPOINTS],
        row["secondary_mean"], row["relation_coverage"], row["prediction_semantics"]] for row in report["rows"]]
    destination.mkdir(parents=True, exist_ok=False)
    with (destination / "scores.tsv").open("x") as stream:
        writer = csv.writer(stream, delimiter="\t", lineterminator="\n")
        writer.writerow(headers)
        writer.writerows(["Unavailable" if value is None else value for value in row] for row in values)
    lines = ["# Native QfO Scientific Scores", "", "| " + " | ".join(headers) + " |",
        "| " + " | ".join(["---"] * len(headers)) + " |"]
    for row in values:
        lines.append("| " + " | ".join("Unavailable" if value is None else f"{value:.6f}"
            if isinstance(value, float) else str(value) for value in row) + " |")
    lines.extend(["", report["timing_disclosure"], "", *["- " + item for item in report["limitations"]], ""])
    with (destination / "scores.md").open("x") as stream:
        stream.write("\n".join(lines))
    report["outputs"] = [record(destination / name) for name in ("scores.tsv", "scores.md")]
    save(destination / "report.json", report)
    return record(destination / "report.json")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--final-admission", nargs=2, metavar=("PATH", "SHA256"))
    args = parser.parse_args()
    ref = None
    if args.final_admission:
        ref = record(args.final_admission[0])
        original.require(ref["sha256"] == args.final_admission[1], "Final admission digest differs")
    print(json.dumps(export(args.output, ref), sort_keys=True))
