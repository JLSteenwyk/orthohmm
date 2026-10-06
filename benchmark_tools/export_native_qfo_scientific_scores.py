"""Combine native accuracy admissions while distinguishing failed timing recovery."""

import argparse
import csv
import json
import math
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools import export_native_qfo_factorial_scores as original
from benchmark_tools.prepare_measurement_failed_native_qfo_pairs import admit_recovery
from benchmark_tools.prepare_ob_candidate_neighborhood import check, record
from benchmark_tools.review_native_factorial_attempt import scheduler_fields
from benchmark_tools.run_measurement_failed_native_qfo_assessment import validate_stage


def endpoints(report):
    assessment, participant = report["assessment"], report["participant"]
    original.require(assessment["participant"] == participant
        and set(assessment["endpoints"]) == set(original.ENDPOINTS), "Incomplete or mixed recovered endpoints")
    scores, details = {}, {}
    for endpoint in original.ENDPOINTS:
        item = assessment["endpoints"][endpoint]
        native = item["native_participant"]
        original.require(native["participant_id"] == participant
            and (item["axes"]["x_axis"], item["axes"]["y_axis"]) == original.AXES[endpoint],
            "Recovered endpoint participant or axes differ")
        x = original.finite(native["metric_x"], endpoint + " x",
            high=1 if endpoint in original.F1_ENDPOINTS else math.inf)
        y = original.finite(native["metric_y"], endpoint + " y", high=1)
        errors = {k: original.finite(native[k], endpoint + " " + k) for k in ("stderr_x", "stderr_y")}
        if endpoint in original.F1_ENDPOINTS:
            value, statistic = original.harmonic_mean(x, y), "F1"
            detail = dict(precision=y, recall=x)
            expected = "harmonic mean of native TPR and PPV"
        else:
            original.require(x == int(x), "Noninteger recovered assessed relation count")
            value = y
            statistic = expected = original.AXES[endpoint][1]
            detail = dict(assessed_relations=int(x))
        original.require(item["score_semantics"] == expected
            and math.isclose(original.finite(item["score"], endpoint + " score", high=1), value,
                rel_tol=0, abs_tol=1e-12), "Recovered score semantics or arithmetic differ")
        scores[endpoint] = value
        details[endpoint] = dict(statistic=statistic, native_standard_errors=errors, **detail)
    mean = sum(scores.values()) / 6
    original.require(math.isclose(original.finite(assessment["secondary_six_metric_mean"], "secondary mean", high=1),
        mean, rel_tol=0, abs_tol=1e-12), "Recovered secondary mean differs")
    return scores, details, mean


def extract_recovered(report, run, plan_ref, evidence):
    original.require(report.get("schema") == "measurement_failed_native_qfo_admission_v1"
        and report.get("status") == "measurement_failed_native_qfo_assessment_admitted"
        and report.get("accuracy_admitted") is True
        and all(report.get(k) is False for k in ("publication_ready", "next_identity_authorized", "automatic_retry",
            "native_inference_reexecuted", "original_native_scheduler_success", "scientific_timings_admitted",
            "eligible_for_timing_comparison"))
        and "resources" in report and report["resources"] is None,
        "Require independently admitted recovered accuracy with failed timing retained")
    index, cell, job = run["index"], run["cell"], report["native_job_id"]
    original.require(type(job) is int and job > 0 and type(report.get("native_index")) is int
        and report["native_index"] == index and report["cell"] == cell, "Recovered admission identity differs")
    stage = original.bound(report["pairs_manifest"], evidence)
    original.require(report["conversion"] == stage and stage["plan"] == plan_ref
        and stage["native_index"] == index and stage["native_job_id"] == job and stage["cell"] == cell
        and stage["participant"] == report["participant"], "Recovered embedded conversion differs")
    validate_stage(stage, report["conversion_scheduler"])
    total, covered = stage["pair_coverage"]["input_accessions"], stage["pair_coverage"]["accessions_in_any_pair"]
    count = stage["retained_pairs"]
    original.require(type(run["genes"]) is int and total == run["genes"] and covered <= 2 * count
        and (count == 0 or covered >= 2)
        and original.finite(stage["pair_coverage"]["fraction_inputs_in_any_pair"], "relation coverage", high=1)
            == covered / total, "Recovered relation coverage differs")
    request = original.bound(stage["request"], evidence)
    recovery = original.bound(stage["scientific_recovery"], evidence)
    original.require(request["plan"] == plan_ref and request["index"] == index and request["job_id"] == job
        and report["scientific_recovery"] == stage["scientific_recovery"], "Recovered request/disposition binding differs")
    admit_recovery(recovery, stage["request"], request, dict(run, repeat=run.get("repeat", 0)))
    execution = original.bound(report["execution_report"], evidence)
    preflight = original.bound(report["preflight"], evidence)
    original.scheduler(report["scheduler"], execution["job_id"], 8)
    for terminal in (report["native_scheduler"], execution["native_scheduler"]):
        fields = scheduler_fields(terminal)
        original.require(fields.get("JobState", fields.get("State")) == "FAILED"
            and fields.get("ExitCode") == "1:0", "Recovered report relabels failed native scheduler")
    original.require(execution.get("schema") == "measurement_failed_native_qfo_execution_v1"
        and execution.get("status") == "process_succeeded_pending_independent_admission"
        and type(execution.get("exit_code")) is int and execution["exit_code"] == 0
        and execution.get("native_index") == index and execution.get("native_job_id") == job
        and execution.get("cell") == cell and execution.get("stage") == stage
        and execution.get("pairs_manifest") == report["pairs_manifest"]
        and all(execution.get(k) is False for k in ("accuracy_admitted", "original_native_scheduler_success",
            "scientific_timings_admitted", "eligible_for_timing_comparison", "native_inference_reexecuted",
            "automatic_retry", "next_identity_authorized", "publication_ready"))
        and "resources" in execution and execution["resources"] is None
        and preflight.get("status") == "running"
        and all(execution.get(k) == v for k, v in preflight.items() if k != "status"),
        "Recovered execution/preflight binding differs")
    source = record(Path(__file__).with_name("admit_measurement_failed_native_qfo_assessment.py"))
    original.require(report["source"] == source and report["fas_sample"]["sample_membership_verified"] is True,
        "Recovered admission source or FAS validation differs")
    evidence.append(source)
    scores, details, mean = endpoints(report)
    return dict(index=index, cell=cell, status="supplied_recovered_scientific_admission", native_job_id=job,
        conversion_job_id=stage["job_id"], assessment_job_id=execution["job_id"], scores=scores,
        endpoint_details=details, secondary_mean=mean, submitted_pairs=count, input_accessions=total,
        relation_accessions=covered, relation_coverage=covered / total, prediction_semantics=stage["semantics"],
        participant=stage["participant"], fas_protocol=report["fas_protocol"], fas_sample=report["fas_sample"],
        measurement_status="failed_timing_scientific_outputs_recovered", resources=None,
        timing_eligible=False, timing_admitted=False, accuracy_admitted=True,
        scientific_recovery=stage["scientific_recovery"])


def collect(plan_path, plan_sha, admissions, recovered_admissions):
    snapshot = original.collect(plan_path, plan_sha, admissions)
    evidence, rows = snapshot["evidence"], snapshot["rows"]
    helpers = [record(Path(__file__).with_name(name)) for name in (
        "export_native_qfo_factorial_scores.py", "prepare_measurement_failed_native_qfo_pairs.py",
        "run_measurement_failed_native_qfo_assessment.py", "review_native_factorial_attempt.py")]
    evidence.extend(helpers)
    plan = original.bound(snapshot["plan"], evidence)
    for row in rows:
        row["measurement_status"] = ("successful_native_terminal_review" if row["status"] == "supplied_native_admission"
            else "no_supplied_accuracy_admission")
        row["accuracy_admitted"] = row["status"] == "supplied_native_admission"
    seen = {r["index"] for r in rows if r["accuracy_admitted"]}
    for path, digest in recovered_admissions:
        report, ref = original.load(path, digest, evidence)
        index = report.get("native_index")
        original.require(type(index) is int and 6 <= index < 13 and index not in seen,
            "Invalid or duplicate native scientific admission index")
        seen.add(index)
        rows[index - 6] = dict(extract_recovered(report, plan["runs"][index], snapshot["plan"], evidence), admission=ref)
    for ref in evidence:
        check(ref)
    snapshot.update(schema="native_qfo_scientific_reporting_snapshot_v1", supplied_admissions=len(seen),
        supplied_recovered_admissions=len(recovered_admissions), recovered_inference_resources_admitted=False,
        recovered_reporting_helpers=helpers)
    snapshot["limitations"] = [
        "Includes supplied successful-measurement and explicit recovered-science admissions, not live or cached scores.",
        "Recovered accuracy does not admit failed native timing; resource values remain null and eligibility false.",
        *snapshot["limitations"][1:],
        "These are direct report/source/arithmetic checks, not a repeated transitive raw scientific or score admission."]
    return snapshot


def export(plan_path, plan_sha, admissions, recovered_admissions, output):
    output = Path(output).absolute()
    original.require(not output.exists() and not output.is_symlink(), "Output already exists")
    result = collect(plan_path, plan_sha, admissions, recovered_admissions)
    headers = ["Index", "Cell", "Accuracy status", "Measurement status", "VGNC F1", "SwissTrees F1",
        "TreeFam-A F1", "GO similarity", "EC similarity", "FAS", "Secondary mean", "Submitted pairs",
        "Input accessions", "Relation accessions", "Relation coverage", "Prediction semantics"]
    values = [[r["index"], r["cell"], r["status"], r["measurement_status"],
        *[r["scores"][e] for e in original.ENDPOINTS], r["secondary_mean"], r["submitted_pairs"],
        r["input_accessions"], r["relation_accessions"], r["relation_coverage"], r["prediction_semantics"]]
        for r in result["rows"]]
    output.mkdir(parents=True, exist_ok=False)
    with (output / "scores.tsv").open("x") as stream:
        writer = csv.writer(stream, delimiter="\t", lineterminator="\n")
        writer.writerow(headers)
        writer.writerows(values)
    lines = ["# Native QfO Scientific Results", "", "| " + " | ".join(headers) + " |",
        "| " + " | ".join(["---"] * len(headers)) + " |"]
    for row in values:
        lines.append("| " + " | ".join("Unavailable" if v is None else f"{v:.6f}" if isinstance(v, float)
            else str(v) for v in row) + " |")
    lines += ["", result["timing_disclosure"], "", *["- " + line for line in result["limitations"]]]
    (output / "scores.md").write_text("\n".join(lines) + "\n")
    result.update(source=record(__file__), outputs=[record(output / name) for name in ("scores.tsv", "scores.md")])
    (output / "report.json").write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--plan", required=True)
    parser.add_argument("--plan-sha256", required=True)
    parser.add_argument("--admission", nargs=2, action="append", default=[], metavar=("PATH", "SHA256"))
    parser.add_argument("--recovered-admission", nargs=2, action="append", default=[], metavar=("PATH", "SHA256"))
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = export(args.plan, args.plan_sha256, args.admission, args.recovered_admission, args.output)
    print(json.dumps(dict(supplied_admissions=result["supplied_admissions"], rows=len(result["rows"]))))
