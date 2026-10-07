"""Report both native execution routes without replaying raw scientific admission."""

import argparse
import csv
import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools import export_native_qfo_scientific_scores as historical
from benchmark_tools import native_factorial_allocated_execution as contract
from benchmark_tools.prepare_allocated_native_factorial_qfo_pairs import admit_conversion
from benchmark_tools.run_allocated_native_factorial_qfo_assessment import validate_stage, execution_spec, ENV_SHA
from benchmark_tools.prepare_ob_candidate_neighborhood import check, record
from benchmark_tools.review_native_factorial_attempt import scheduler_fields
from benchmark_tools.verify_lineage_native_provenance import same

original = historical.original


def native_terminal(value, job, request_ref, root):
    if value.get("source") == "live_controller":
        observed = value["observation"]
        original.require(observed.get("command") == ["scontrol", "show", "job", str(job), "--oneliner"]
            and type(observed.get("returncode")) is int and observed["returncode"] == 0,
            "Retained native controller query differs")
        parsed = contract.validate_controller(observed["stdout"], job, "terminal", command=str(contract.SCRIPT),
            cwd=str(root), time_limit="1-02:00:00", allocation_mode="shared")
        original.require(parsed["fields"].get("Comment") == request_ref["sha256"]
            and parsed["fields"].get("JobName") == "orthohmm_allocated_factorial",
            "Retained native request comment or job name differs")
    else:
        original.require(value.get("source") == "fresh_accounting_after_controller_expiry",
            "Unknown native terminal observation route")
        controller = value["controller_observation"]
        original.require(controller.get("command") == ["scontrol", "show", "job", str(job), "--oneliner"]
            and type(controller.get("returncode")) is int and controller["returncode"] != 0
            and "Invalid job id specified" in controller.get("stderr", ""),
            "Accounting fallback lacks retained specific controller expiry")
        observed = value["observation"]
        command = ["sacct", "-X", "-j", str(job), "-n", "-P",
            "--format=JobIDRaw,State,ExitCode,AllocCPUS,ReqMem,NodeList,Partition,TimelimitRaw,Submit,Start,End,JobName"]
        original.require(observed.get("command") == command, "Retained native accounting query differs")
        parsed = contract.terminal_accounting(observed["stdout"], job)
    original.require(same(parsed, value["verified"]), "Retained native scheduler parser differs")
    fields = scheduler_fields(value)
    original.require(fields.get("JobState", fields.get("State")) == "COMPLETED" and fields.get("ExitCode") == "0:0",
        "Allocated native scheduler is not successful")


def extract(report, run, plan_ref, evidence):
    original.require(report.get("schema") == "allocated_native_factorial_qfo_admission_v1"
        and report.get("status") == "allocated_native_factorial_qfo_assessment_admitted"
        and report.get("accuracy_admitted") is True
        and all(report.get(k) is False for k in ("publication_ready", "next_identity_authorized", "automatic_retry")),
        "Require independently admitted allocated native accuracy")
    index, job = report.get("native_index"), report.get("native_job_id")
    original.require(type(index) is int and index in (10, 11, 12) and index == run["index"]
        and type(job) is int and job > 0 and report.get("cell") == run["cell"], "Allocated admission identity differs")
    amendment_ref = report["amendment"]
    execution_amendment = original.bound(amendment_ref, evidence)
    plan = contract.validate_amendment(execution_amendment)
    original.require(execution_amendment["historical_plan"] == plan_ref and plan["runs"][index] == run,
        "Allocated amendment is outside the reported frozen plan")
    root = Path(execution_amendment["root"])
    evidence.extend(execution_amendment["new_sources"])
    stage = original.bound(report["pairs_manifest"], evidence)
    original.require(report["conversion"] == stage and stage["amendment"] == amendment_ref
        and stage["plan"] == plan_ref and stage["native_index"] == index
        and stage["native_job_id"] == job and stage["cell"] == run["cell"]
        and stage["participant"] == report["participant"] and stage["input_fastas"] == run["inputs"],
        "Allocated conversion differs from admission/inputs")
    validate_stage(stage, report["conversion_scheduler"])
    request = original.bound(stage["request"], evidence)
    contract.validate_request(request, amendment_ref, execution_amendment, job)
    original.require(request["index"] == index, "Allocated request index differs")
    review = original.bound(stage["terminal_review"], evidence)
    admit_conversion(review, stage["request"], request, run)
    expected_review_source = record(root / "benchmark_tools/review_allocated_native_factorial_attempt.py")
    original.require(review["source"] == expected_review_source
        and review["common_reviewer_source"] == record(root / "benchmark_tools/review_native_factorial_attempt.py"),
        "Allocated terminal reviewer source differs")
    evidence.extend([review["source"], review["common_reviewer_source"]])
    outputs = original.bound(review["reviews"]["outputs_or_failure"], evidence)
    cpu_ids = review["native_cpu_ids"]
    original.require(type(cpu_ids) is list and len(cpu_ids) == 32
        and all(type(cpu) is int and cpu >= 0 for cpu in cpu_ids)
        and len(set(cpu_ids)) == 32 and cpu_ids == sorted(cpu_ids)
        and outputs.get("native_cpu_ids") == stage.get("native_cpu_ids") == cpu_ids,
        "Allocated reporting CPU selection differs")
    ready = record(Path(run["output_root"]) / "measurement/ready.json")
    original.require(outputs.get("schema") == "allocated_native_factorial_output_review_v1"
        and outputs.get("native_outputs_validated") is True and outputs.get("request") == stage["request"]
        and outputs.get("plan") == plan_ref and outputs.get("amendment") == amendment_ref
        and outputs.get("index") == index and outputs.get("job_id") == job and outputs.get("cell") == run["cell"]
        and outputs.get("execution_scope") == contract.SCOPE
        and outputs.get("source") == record(root / "benchmark_tools/validate_allocated_native_factorial_outputs.py")
        and outputs.get("semantic_validator_source") == record(root / "benchmark_tools/validate_native_factorial_outputs.py")
        and outputs.get("allocated_ready") == stage.get("allocated_ready") == ready
        and outputs.get("gene_ownership_sha256") == stage["gene_ownership_sha256"]
        and stage["native_input"] in outputs["checked_files"], "Allocated scientific output/placement binding differs")
    evidence.extend([ready, outputs["source"], outputs["semantic_validator_source"]])
    for name, expected in (("source", "prepare_allocated_native_factorial_qfo_pairs.py"),
            ("conversion_kernel_source", "prepare_native_factorial_qfo_pairs.py")):
        canonical = record(root / "benchmark_tools" / expected)
        original.require(stage[name] == canonical, "Allocated conversion source differs: " + name)
        evidence.append(canonical)
    suffix = "native/orthohmm_phylogeny/orthohmm_pairwise_orthologs.tsv" if stage["conversion_kind"] == "native" else \
        "native/orthohmm_working_res/orthohmm_edges_clustered.txt"
    original.require(stage["native_input"]["path"] == str(Path(run["output_root"]) / suffix),
        "Allocated converted prediction path differs")
    count = stage["retained_pairs"]
    total, covered = stage["pair_coverage"]["input_accessions"], stage["pair_coverage"]["accessions_in_any_pair"]
    original.require(type(run["genes"]) is int and total == run["genes"]
        and covered <= 2*count and (count == 0 or covered >= 2), "Allocated relation coverage differs")
    if stage["conversion_kind"] == "native":
        original.require(type(outputs["phylogeny"]["native_pair_rows"]) is int
            and outputs["phylogeny"]["native_pair_rows"] == count, "Allocated native pair count differs")
    execution = original.bound(report["execution_report"], evidence)
    preflight = original.bound(report["preflight"], evidence)
    env_ref = record(root / "benchmark_tools/results/qfo_assessment_environment_20260917.json")
    original.require(env_ref["sha256"] == ENV_SHA
        and report.get("environment_manifest") == stage.get("environment_manifest") == env_ref,
        "Allocated QfO reference environment differs")
    manifest = original.bound(env_ref, evidence)
    expected = execution_spec(root, report["pairs_manifest"], stage, manifest, {})
    original.require(all(execution.get(key) == value for key, value in expected.items() if key != "status"),
        "Allocated scoring command, namespace or environment differs")
    original.scheduler(report["scheduler"], execution["job_id"], 8)
    original.require(execution.get("schema") == "allocated_native_factorial_qfo_execution_v1"
        and execution.get("status") == "process_succeeded_pending_independent_admission"
        and type(execution.get("exit_code")) is int and execution["exit_code"] == 0
        and execution.get("native_index") == index and execution.get("native_job_id") == job
        and execution.get("cell") == run["cell"] and execution.get("stage") == stage
        and execution.get("pairs_manifest") == report["pairs_manifest"] and execution.get("amendment") == amendment_ref
        and execution.get("source") == record(root / "benchmark_tools/run_allocated_native_factorial_qfo_assessment.py")
        and all(execution.get(k) is False for k in ("accuracy_admitted", "publication_ready", "automatic_retry",
            "native_inference_reexecuted", "next_identity_authorized"))
        and preflight.get("status") == "running"
        and all(execution.get(k) == v for k, v in preflight.items() if k != "status"),
        "Allocated scoring execution/preflight binding differs")
    for terminal in (report["native_scheduler"], execution["native_scheduler"]):
        native_terminal(terminal, job, stage["request"], root)
    source = record(root / "benchmark_tools/admit_allocated_native_factorial_qfo_assessment.py")
    original.require(report["source"] == source and report["fas_sample"]["sample_membership_verified"] is True,
        "Allocated admission source or FAS sample validation differs")
    evidence.extend([source, execution["source"]])
    scores, details, mean = historical.endpoints(report)
    resources = review["resources"]
    original.require(isinstance(resources, dict) and set(resources) == set(original.SCOPES),
        "Missing or invalid terminal native resource observations")
    for key, value in resources.items():
        original.finite(value, key)
    original.require(type(resources["peak_memory_bytes"]) is int, "Noninteger native-step peak memory")
    foreign = original.finite(review["whole_run_maximum_foreign_average_cores"], "foreign CPU")
    return dict(index=index, cell=run["cell"], status="supplied_allocated_native_admission", native_job_id=job,
        conversion_job_id=stage["job_id"], assessment_job_id=execution["job_id"], scores=scores,
        endpoint_details=details, secondary_mean=mean, submitted_pairs=count, input_accessions=total,
        relation_accessions=covered, relation_coverage=covered/total, prediction_semantics=stage["semantics"],
        participant=stage["participant"], fas_protocol=report["fas_protocol"], fas_sample=report["fas_sample"],
        measurement_status="successful_allocated_native_terminal_review", accuracy_admitted=True,
        amendment=amendment_ref, native_cpu_ids=cpu_ids, resources=resources, resource_scopes=original.SCOPES,
        whole_run_maximum_foreign_average_cores=foreign, scientific_timings_admitted=False,
        resource_observation_scope="terminal-reviewed shared-host native command; not isolated algorithm cost")


def collect(plan_path, plan_sha, admissions, recovered_admissions, allocated_admissions):
    result = historical.collect(plan_path, plan_sha, admissions, recovered_admissions)
    evidence = result["evidence"]
    plan = original.bound(result["plan"], evidence)
    seen = {row["index"] for row in result["rows"] if row["accuracy_admitted"]}
    for path, digest in allocated_admissions:
        report, ref = original.load(path, digest, evidence)
        index = report.get("native_index")
        original.require(type(index) is int and index in (10,11,12) and index not in seen,
            "Invalid or duplicate allocated scientific admission index")
        row = extract(report, plan["runs"][index], result["plan"], evidence)
        result["rows"][index-6] = dict(row, admission=ref)
        seen.add(index)
    helpers = [record(Path(__file__).with_name(name)) for name in (
        "export_native_qfo_scientific_scores.py", "native_factorial_allocated_execution.py",
        "prepare_allocated_native_factorial_qfo_pairs.py", "run_allocated_native_factorial_qfo_assessment.py")]
    evidence.extend([*helpers, record(__file__)])
    for ref in evidence:
        check(ref)
    result.update(schema="allocated_native_qfo_scientific_reporting_snapshot_v1", supplied_admissions=len(seen),
        supplied_allocated_admissions=len(allocated_admissions), allocated_reporting_helpers=helpers)
    result["limitations"] += [
        "Allocated-route rows require separate successful terminal and independent QfO admissions; startup/live results are excluded.",
        "This rechecks direct metadata, source bindings and arithmetic, not transitive raw files or a new scientific admission.",
        "Historical rows are retained exactly; failed-timing recovered science remains timing-ineligible with null resources.",
        "Reported resource observations preserve native wrapper/launcher scopes and unknown tool-dependent contention."]
    return result


def export(plan_path, plan_sha, admissions, recovered_admissions, allocated_admissions, output):
    output = Path(output).absolute()
    original.require(not output.exists() and not output.is_symlink(), "Output already exists")
    result = collect(plan_path, plan_sha, admissions, recovered_admissions, allocated_admissions)
    headers = ["Index", "Cell", "Accuracy status", "Measurement status", "VGNC F1", "SwissTrees F1",
        "TreeFam-A F1", "GO similarity", "EC similarity", "FAS", "Secondary mean", "Submitted pairs",
        "Input accessions", "Relation accessions", "Relation coverage", "Prediction semantics"]
    values = [[row["index"], row["cell"], row["status"], row["measurement_status"],
        *[row["scores"][e] for e in original.ENDPOINTS], row["secondary_mean"], row["submitted_pairs"],
        row["input_accessions"], row["relation_accessions"], row["relation_coverage"], row["prediction_semantics"]]
        for row in result["rows"]]
    output.mkdir(parents=True, exist_ok=False)
    with (output / "scores.tsv").open("x") as stream:
        writer = csv.writer(stream, delimiter="\t", lineterminator="\n")
        writer.writerow(headers)
        writer.writerows(values)
    lines = ["# Native QfO Scientific Results", "", "| " + " | ".join(headers) + " |",
        "| " + " | ".join(["---"]*len(headers)) + " |"]
    for row in values:
        lines.append("| " + " | ".join("Unavailable" if value is None else f"{value:.6f}"
            if isinstance(value, float) else str(value) for value in row) + " |")
    lines += ["", result["timing_disclosure"], "", *["- " + item for item in result["limitations"]]]
    with (output / "scores.md").open("x") as stream:
        stream.write("\n".join(lines) + "\n")
    result.update(source=record(__file__), outputs=[record(output / name) for name in ("scores.tsv", "scores.md")])
    with (output / "report.json").open("x") as stream:
        json.dump(result, stream, sort_keys=True, indent=2, allow_nan=False)
        stream.write("\n")
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--plan", required=True)
    parser.add_argument("--plan-sha256", required=True)
    for name in ("admission", "recovered-admission", "allocated-admission"):
        parser.add_argument("--" + name, nargs=2, action="append", default=[], metavar=("PATH", "SHA256"))
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = export(args.plan, args.plan_sha256, args.admission, args.recovered_admission, args.allocated_admission, args.output)
    print(json.dumps(dict(supplied_admissions=result["supplied_admissions"],
        supplied_allocated_admissions=result["supplied_allocated_admissions"], rows=len(result["rows"]))))
