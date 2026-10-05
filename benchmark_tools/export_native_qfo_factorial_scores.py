"""Report supplied full-native QfO admissions, never cached scores or live outputs."""

import argparse
import csv
import json
import math
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.export_native_factorial_progress import DISCLOSURE, SCOPES, finite, load, require
from benchmark_tools.prepare_ob_candidate_neighborhood import check, record
from benchmark_tools.qfo_summarize_scores import harmonic_mean
from benchmark_tools.validate_qfo_native_assessment import AXES


CELLS = ("p0_c0_r0", "p0_c0_r1", "p0_c1_r0", "p0_c1_r1",
         "p1_c0_r1", "p1_c1_r0", "p1_c1_r1")
ENDPOINTS = ("VGNC", "SwissTrees", "TreeFam-A", "GO", "EC", "FAS")
F1_ENDPOINTS = frozenset(ENDPOINTS[:3])


def bound(ref, evidence):
    value, observed = load(ref["path"], ref["sha256"], evidence)
    require(observed == ref, "Report reference metadata differs")
    return value


def scheduler(value, job, cpus):
    require(value.get("JobIDRaw") == str(job) and value.get("State") == "COMPLETED"
            and value.get("ExitCode") == "0:0" and value.get("NodeList") == "bizon"
            and value.get("AllocCPUS") == str(cpus), "Unsuccessful or mismatched scheduler")


def extract(report, run, plan_ref, evidence):
    require(report.get("schema") == "full_native_factorial_qfo_admission_v1"
            and report.get("status") == "full_native_factorial_qfo_assessment_admitted"
            and report.get("accuracy_admitted") is True
            and all(report.get(k) is False for k in
                    ("publication_ready", "next_identity_authorized", "automatic_retry")),
            "Require supplied successful full-native QfO admission")
    index, cell = run["index"], run["cell"]
    job = report["native_job_id"]
    require(type(job) is int and job > 0 and report["native_index"] == index
            and report["cell"] == cell, "Admission identity differs")
    stage = bound(report["pairs_manifest"], evidence)
    require(report["conversion"] == stage, "Embedded conversion differs from manifest")
    kind = "native" if cell.endswith("r1") else "group"
    semantics = "native phylogenetically inferred pairs" if kind == "native" else "cross-species group-derived clique pairs"
    require(stage["schema"] == "full_native_factorial_qfo_conversion_v1"
            and stage["status"] == "full_native_factorial_qfo_pairs_prepared_unscored"
            and stage["plan"] == plan_ref and stage["native_index"] == index
            and stage["native_job_id"] == job and stage["cell"] == cell
            and stage["conversion_kind"] == kind and stage["semantics"] == semantics
            and stage["participant"] == report["participant"] == "ohmm_qfo_full_native_" + cell
            and all(stage.get(k) is False for k in ("accuracy_evaluated", "publication_ready",
                    "native_inference_reexecuted", "next_identity_authorized", "automatic_retry")),
            "Conversion native identity or semantics differ")
    scheduler(report["conversion_scheduler"], stage["job_id"], 2)
    counts = [stage[k] for k in ("total_pairs", "retained_pairs", "expected_pairs", "removed_mapping_pairs")]
    require(all(type(v) is int for v in counts) and 0 <= counts[0] == counts[1] == counts[2]
            and counts[3] == 0 and stage["empty_predictions"] is (counts[0] == 0), "Pair counts differ")
    coverage = stage["pair_coverage"]
    total, covered = coverage["input_accessions"], coverage["accessions_in_any_pair"]
    require(type(total) is int and type(covered) is int and 0 <= covered <= total
            and type(run["genes"]) is int and total == run["genes"] and total > 0
            and type(coverage["pair_rows"]) is int and coverage["pair_rows"] == counts[0]
            and covered <= 2 * counts[0] and (counts[0] == 0 or covered >= 2)
            and finite(coverage["fraction_inputs_in_any_pair"], "relation coverage", high=1) == covered / total,
            "Relation coverage differs")

    request = bound(stage["request"], evidence)
    review = bound(stage["terminal_review"], evidence)
    require(request["plan"] == plan_ref and request["index"] == index and request["job_id"] == job
            and review["request"] == stage["request"] and review["plan"] == plan_ref
            and review["schema"] == "native_factorial_terminal_review_v1"
            and review["index"] == index and review["cell"] == cell
            and review["dataset"] == "qfo_corrected" and review["job_id"] == job
            and review["status"] == "native_success" and review["scheduler_state"] == "COMPLETED"
            and review["scheduler_exit_code"] == "0:0"
            and review["execution_scope"] == "shared_host_matched_resources"
            and review["resource_scopes"] == SCOPES and review["uncontended_timing"] is False
            and all(review.get(k) is True for k in ("terminal_reviewed", "native_outputs_validated",
                    "primary_resources_replayed", "shared_host_resources_reviewed")),
            "Native terminal review/request binding differs")
    execution = bound(report["execution_report"], evidence)
    preflight = bound(report["preflight"], evidence)
    scheduler(report["scheduler"], execution["job_id"], 8)
    require(execution["schema"] == "full_native_factorial_qfo_execution_v1"
            and execution["status"] == "process_succeeded_pending_independent_admission"
            and execution["exit_code"] == 0 and execution["accuracy_admitted"] is False
            and execution["native_index"] == index and execution["native_job_id"] == job
            and execution["cell"] == cell and execution["stage"] == stage
            and execution["pairs_manifest"] == report["pairs_manifest"]
            and preflight["status"] == "running"
            and all(execution.get(k) == v for k, v in preflight.items() if k != "status"),
            "Assessment execution/preflight binding differs")
    source = record(Path(__file__).with_name("admit_native_factorial_qfo_assessment.py"))
    require(report["source"] == source, "Native admission source differs")
    evidence.append(source)

    assessment = report["assessment"]
    participant = stage["participant"]
    require(assessment["participant"] == participant
            and set(assessment["endpoints"]) == set(ENDPOINTS), "Incomplete or mixed endpoints")
    scores, details = {}, {}
    for endpoint in ENDPOINTS:
        result = assessment["endpoints"][endpoint]
        native = result["native_participant"]
        require(native["participant_id"] == participant, "Mixed endpoint participants")
        require((result["axes"]["x_axis"], result["axes"]["y_axis"]) == AXES[endpoint], "Wrong metric axes")
        x = finite(native["metric_x"], endpoint + " x", high=1 if endpoint in F1_ENDPOINTS else math.inf)
        y = finite(native["metric_y"], endpoint + " y", high=1)
        errors = {k: finite(native[k], endpoint + " " + k) for k in ("stderr_x", "stderr_y")}
        if endpoint in F1_ENDPOINTS:
            value, statistic = harmonic_mean(x, y), "F1"
            detail = dict(precision=y, recall=x)
            expected_semantics = "harmonic mean of native TPR and PPV"
        else:
            require(x == int(x), "Noninteger challenge-assessed relation count")
            value, statistic = y, AXES[endpoint][1]
            detail = dict(assessed_relations=int(x))
            expected_semantics = statistic
        require(result["score_semantics"] == expected_semantics, "Wrong endpoint score semantics")
        require(math.isclose(finite(result["score"], endpoint + " score", high=1), value,
                             rel_tol=0, abs_tol=1e-12), "Endpoint score arithmetic differs")
        scores[endpoint] = value
        details[endpoint] = dict(statistic=statistic, native_standard_errors=errors, **detail)
    mean = sum(scores.values()) / 6
    require(math.isclose(finite(assessment["secondary_six_metric_mean"], "secondary mean", high=1),
                         mean, rel_tol=0, abs_tol=1e-12), "Secondary mean arithmetic differs")
    require(report["fas_sample"]["sample_membership_verified"] is True, "Unverified native FAS sample")
    return dict(index=index, cell=cell, status="supplied_native_admission", native_job_id=job,
                conversion_job_id=stage["job_id"], assessment_job_id=execution["job_id"],
                scores=scores, endpoint_details=details, secondary_mean=mean,
                submitted_pairs=counts[0], input_accessions=total, relation_accessions=covered,
                relation_coverage=covered / total, prediction_semantics=semantics,
                participant=participant, fas_protocol=report["fas_protocol"], fas_sample=report["fas_sample"])


def collect(plan_path, plan_sha, admissions):
    helpers = [record(Path(__file__).with_name(name)) for name in (
        "export_native_factorial_progress.py", "prepare_ob_candidate_neighborhood.py",
        "qfo_summarize_scores.py", "validate_qfo_native_assessment.py")]
    evidence = list(helpers)
    plan, plan_ref = load(plan_path, plan_sha, evidence)
    runs = plan["runs"]
    require(len(runs) == 13 and all(type(r["index"]) is int for r in runs)
            and [r["index"] for r in runs] == list(range(13))
            and [(r["dataset"], r["cell"]) for r in runs[6:]] == [("qfo_corrected", c) for c in CELLS],
            "Wrong seven-identity native QfO plan")
    rows = [dict(index=i + 6, cell=cell, status="no_supplied_native_admission",
                 scores={e: None for e in ENDPOINTS}, secondary_mean=None, submitted_pairs=None,
                 input_accessions=None, relation_accessions=None, relation_coverage=None,
                 native_job_id=None, conversion_job_id=None, assessment_job_id=None,
                 prediction_semantics="native phylogenetically inferred pairs" if cell.endswith("r1")
                 else "cross-species group-derived clique pairs") for i, cell in enumerate(CELLS)]
    seen = set()
    for path, digest in admissions:
        admission, ref = load(path, digest, evidence)
        index = admission.get("native_index")
        require(type(index) is int and 6 <= index < 13 and index not in seen,
                "Invalid or duplicate native admission index")
        seen.add(index)
        rows[index - 6] = dict(extract(admission, runs[index], plan_ref, evidence), admission=ref)
    for ref in evidence:
        check(ref)
    return dict(schema="native_qfo_reporting_snapshot_v1", plan=plan_ref, evidence=evidence, helpers=helpers, rows=rows,
                supplied_admissions=len(seen), new_scoring_or_admission=False, publication_ready=False,
                timing_disclosure=DISCLOSURE, limitations=[
                    "Only supplied full-native admissions; missing values are neither zero nor live-job statuses.",
                    "The seven fresh identities exclude reused P1C0R0; historical/cached scores are not substituted.",
                    "Initial HMM search remains on; R-off group cliques and R-on resolved pairs have different semantics.",
                    "VGNC/SwissTrees/TreeFam-A are F1; GO/EC similarity and FAS are not F1.",
                    "The six-metric mean is a project-defined secondary summary, not official QfO F1.",
                    "Relation coverage uses all inference inputs, not reference-family coverage or accuracy.",
                    "Direct report checks and endpoint arithmetic do not repeat transitive raw admission or biological validation.",
                    "FAS sample/population limits and native SEM semantics are retained; no paired uncertainty or superiority claim.",
                    "No inference, scoring, admission, job launch, timing correction or isolated-performance ranking."])


def export(plan_path, plan_sha, admissions, output):
    output = Path(output).absolute()
    require(not output.exists() and not output.is_symlink(), "Output already exists")
    result = collect(plan_path, plan_sha, admissions)
    headers = ["Index", "Cell", "Status", "VGNC F1", "SwissTrees F1", "TreeFam-A F1",
               "GO similarity", "EC similarity", "FAS", "Secondary mean", "Submitted pairs",
               "Input accessions", "Relation accessions", "Relation coverage", "Prediction semantics"]
    values = [[r["index"], r["cell"], r["status"], *[r["scores"][e] for e in ENDPOINTS],
               r["secondary_mean"], r["submitted_pairs"], r["input_accessions"], r["relation_accessions"],
               r["relation_coverage"], r["prediction_semantics"]] for r in result["rows"]]
    output.mkdir(parents=True, exist_ok=False)
    with (output / "scores.tsv").open("x") as stream:
        writer = csv.writer(stream, delimiter="\t", lineterminator="\n")
        writer.writerow(headers)
        writer.writerows(values)
    lines = ["# Full-Native QfO Factorial", "", "| " + " | ".join(headers) + " |",
             "| " + " | ".join(["---"] * len(headers)) + " |"]
    for row in values:
        lines.append("| " + " | ".join("Unavailable" if v is None else f"{v:.6f}" if isinstance(v, float)
                                       else str(v) for v in row) + " |")
    lines += ["", result["timing_disclosure"], "", *["- " + text for text in result["limitations"]]]
    (output / "scores.md").write_text("\n".join(lines) + "\n")
    result.update(source=record(__file__), outputs=[record(output / name) for name in ("scores.tsv", "scores.md")])
    (output / "report.json").write_text(json.dumps(result, sort_keys=True, indent=2) + "\n")
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--plan", required=True)
    parser.add_argument("--plan-sha256", required=True)
    parser.add_argument("--admission", nargs=2, action="append", default=[], metavar=("PATH", "SHA256"))
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    report = export(args.plan, args.plan_sha256, args.admission, args.output)
    print(json.dumps(dict(supplied_admissions=report["supplied_admissions"], rows=len(report["rows"]))))
