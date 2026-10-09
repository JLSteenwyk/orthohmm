"""Integrate reviewed terminal failures without changing four admitted QfO cells."""

import argparse
from copy import deepcopy
import csv
import json
import math
import os
from pathlib import Path
import subprocess

from benchmark_tools.derive_threadripper_resources import SCOPES
from benchmark_tools.prepare_ob_candidate_neighborhood import check, record
from benchmark_tools.probe_dgx_step_separation import save
from benchmark_tools.review_native_factorial_attempt import DISCLOSURE, scheduler_fields
from benchmark_tools.run_native_factorial_cost import ROOT, read
from benchmark_tools.validate_native_factorial_outputs import require


PINS = {
    "parent": ("benchmark_tools/results/PUBLICATION_MAIN_TEXT_20261007_v4.md",
        "0b7012e3dd46bd3b028857bdb3cbe3a9e9048936b90f3df9e4833cb342d5705a"),
    "snapshot": ("benchmark_tools/results/native_qfo_scientific_scores_20261007_v3/report.json",
        "7aab1cbd31fb650df42e6e80a14ee0167b35a810f08202768c89c514755211a2"),
    "native11": ("benchmark_tools/results/native11_qfo_scoring_failure_addendum_20261008_v1/report.json",
        "4f1a5c536af8ce992d46866374b3ab44ee3158077316d69ba2e9f5fb34cee203"),
}
REVIEW_PATH = "benchmarks/work/native12_composed_terminal_review_20261008_v1/review.json"
REVIEWER_SHA = "08b170ba1ee5aabd5665dd12fc359fc7f5ffe3d0834b0afb2e73289fa61ed649"
REQUEST_SHA = "70bde16142da526de6f198aaf255a224b0f4e10c2bd79f98f0b1fd7d70269f89"
ENDPOINTS = ("VGNC", "SwissTrees", "TreeFam-A", "GO", "EC", "FAS")
CELLS = ("p0_c0_r0", "p0_c0_r1", "p0_c1_r0", "p0_c1_r1", "p1_c0_r1", "p1_c1_r0", "p1_c1_r1")
ADMITTED = {6, 7, 8, 10}
ANCHOR = "### Fixed-Stratum Errors And Estimated Sequence Divergence\n"


def finite_resources(resources):
    require(set(resources) == set(SCOPES)
        and all(type(value) in (int, float) and math.isfinite(value) and value > 0
            for value in resources.values())
        and type(resources["peak_memory_bytes"]) is int, "Invalid reviewed resource observations")


def validate(snapshot, native11, review, components):
    require(snapshot.get("schema") == "allocated_native_qfo_scientific_reporting_snapshot_v1"
        and snapshot.get("supplied_admissions") == 4 and snapshot.get("publication_ready") is False,
        "Require frozen four-cell snapshot")
    rows = deepcopy(snapshot["rows"])
    require(len(rows) == 7 and [row["index"] for row in rows] == list(range(6, 13))
        and [row["cell"] for row in rows] == list(CELLS), "Changed cell inventory")
    for row in rows:
        require(set(row["scores"]) == set(ENDPOINTS)
            and row["accuracy_admitted"] is (row["index"] in ADMITTED), "Changed admitted cells")
        if row["accuracy_admitted"]:
            require(all(type(value) in (int, float) and math.isfinite(value) and 0 <= value <= 1
                for value in row["scores"].values())
                and type(row["secondary_mean"]) in (int, float)
                and math.isclose(sum(row["scores"].values()) / 6, row["secondary_mean"], abs_tol=1e-12, rel_tol=0),
                "Changed endpoint scores or secondary mean")
        else:
            require(all(value is None for value in row["scores"].values())
                and row["secondary_mean"] is None, "Missing scores cannot be imputed")
    require(rows[1].get("resources") is None and rows[1].get("timing_eligible") is False,
        "Recovered accuracy cannot repair failed timing")
    require(native11.get("schema") == "composed_native_qfo_failure_reporting_v1"
        and all(native11.get(key) is False for key in ("publication_ready", "new_scientific_admission",
            "automatic_retry", "native_inference_reexecuted", "frozen_manuscript_replaced")),
        "Require retained native11 failure report")
    failed = native11["row"]
    require(failed.get("native_index") == 11 and failed.get("cell") == CELLS[5]
        and failed.get("native_job_id") == 23985 and failed.get("assessment_job_id") == 24038
        and failed.get("conversion_job_id") == 24035
        and failed.get("inference_status") == "successful_composed_terminal_review"
        and failed.get("conversion_status") == "successful_unscored_group_cliques"
        and failed.get("scoring_status") == "OUT_OF_MEMORY" and failed.get("failed_endpoint") == "FAS"
        and failed.get("failed_endpoint_exit_code") == 137 and failed.get("completed_endpoint_task_count") == 5
        and failed.get("admitted_endpoint_score_count") == 0 and failed.get("accuracy_admitted") is False
        and failed.get("secondary_six_metric_mean") is None, "Changed native11 scoring failure")
    require(failed.get("submitted_pairs") == 11755521 and failed.get("input_accessions") == 984137
        and failed.get("relation_accessions") == 585180
        and failed.get("relation_coverage") == 585180 / 984137
        and failed.get("prediction_semantics") == "cross-species group-derived clique pairs"
        and failed.get("resource_scopes") == SCOPES, "Changed native11 coverage or resource scope")
    finite_resources(failed["resources"])
    require(review.get("schema") == "native12_composed_terminal_review_v1"
        and review.get("status") == "native_failure_retained" and review.get("index") == 12
        and review.get("job_id") == 24036 and review.get("cell") == CELLS[6]
        and review.get("scheduler_state") == "FAILED" and review.get("scheduler_exit_code") == "1:0"
        and review.get("resource_scopes") == SCOPES
        and all(review.get(key) is True for key in ("terminal_reviewed", "primary_resources_replayed",
            "shared_host_resources_reviewed"))
        and all(review.get(key) is False for key in ("native_outputs_validated", "accuracy_evaluated",
            "automatic_retry", "publication_ready", "uncontended_timing", "scientific_timings_admitted")),
        "Require independently reviewed native12 failure, not success or incomplete review")
    fields = scheduler_fields(components["scheduler"])
    require(fields.get("JobState", fields.get("State")) == "FAILED" and fields.get("ExitCode") == "1:0",
        "Native terminal accounting differs")
    outcome, replay, resources, environment = (components[key] for key in (
        "outputs_or_failure", "resource_replay", "resources", "environment"))
    require(outcome == dict(status="native_failure_retained", native_outcome="exited_nonzero",
            native_exit_code=-11, native_outputs_validated=False, accuracy_evaluated=False, automatic_retry=False)
        and replay.get("schema") == "native12_composed_resource_replay_summary_v1"
        and replay.get("native_outcome") == resources.get("native_outcome") == "exited_nonzero"
        and type(replay.get("native_exit_code")) is int
        and replay["native_exit_code"] == resources.get("native_exit_code") == -11
        and replay.get("full_replay_executed") is True
        and replay.get("measured_matches_retained_wrapper") is True,
        "Require actual native SIGSEGV, not observer wrapper success")
    finite_resources(review["resources"])
    require(resources.get("primary") == review["resources"] and resources.get("primary_scopes") == SCOPES
        and resources.get("shared_host_observation") is True and resources.get("uncontended_timing") is False
        and environment.get("status") == "shared_environment_replayed"
        and environment.get("sampled_environment_evidence_valid") is True
        and all(environment.get(key) is False for key in ("uncontended_timing",
            "background_cpu_used_for_eligibility", "pressure_thresholds_used_for_eligibility")),
        "Changed resource scope or contention policy")
    done, native = components["native_done"], components["native_receipt"]
    require(type(done.get("exit_code")) is int and done["exit_code"] == -11
        and done.get("timed_out") is False
        and native.get("schema") == "allocated_native_factorial_execution_v1"
        and native.get("status") == "native_factorial_running" and native.get("index") == 12
        and native.get("cell") == CELLS[6]
        and components["missing_final_outputs"] == {"metrics.json": True, "native/orthohmm_orthogroups.txt": True},
        "Native receipt or final-output absence does not support failure description")
    rows[5].update(status="scoring_OUT_OF_MEMORY_no_admission", native_job_id=23985,
        conversion_job_id=24035, assessment_job_id=24038,
        submitted_pairs=failed["submitted_pairs"], input_accessions=failed["input_accessions"],
        relation_accessions=failed["relation_accessions"], relation_coverage=failed["relation_coverage"],
        prediction_semantics=failed["prediction_semantics"], inference_status=failed["inference_status"],
        conversion_status=failed["conversion_status"], scoring_status="OUT_OF_MEMORY",
        resources=failed["resources"], resource_scopes=failed["resource_scopes"])
    rows[6].update(status="native_SIGSEGV_no_admission", native_job_id=24036,
        inference_status="FAILED_SIGSEGV", conversion_status="not_run_failed_inference",
        scoring_status="not_run_failed_inference", resources=review["resources"], resource_scopes=SCOPES,
        resource_observation_scope="failed native command; not successful inference timing",
        submitted_pairs=None, relation_accessions=None, relation_coverage=None)
    return rows


def reviewer_accounting():
    command = ["sacct", "-X", "-j", "24080", "-n", "-P",
        "--format=JobIDRaw,State,ExitCode,AllocCPUS,ReqMem"]
    result = subprocess.run(command, capture_output=True, text=True, timeout=5, check=True)
    rows = [line.split("|") for line in result.stdout.splitlines() if line.strip()]
    require(rows == [["24080", "COMPLETED", "0:0", "2", "128G"]], "Independent reviewer has not completed successfully")
    return dict(command=command, stdout=result.stdout, stderr=result.stderr, returncode=result.returncode)


def inputs(review_path, review_sha):
    refs, docs = {}, {}
    for key, (relative, digest) in PINS.items():
        ref = record(ROOT / relative)
        require(ref["sha256"] == digest, "Changed frozen reporting input: " + key)
        refs[key] = ref
        docs[key] = Path(ref["path"]).read_text() if key == "parent" else read(ref)
    require(Path(review_path).absolute() == ROOT / REVIEW_PATH, "Wrong terminal review path")
    refs["review"] = record(review_path)
    require(refs["review"]["sha256"] == review_sha, "Terminal review digest differs")
    docs["review"] = review = read(refs["review"])
    require(review.get("source") == record(ROOT / "benchmark_tools/review_native12_composed_attempt.py")
        and review["source"]["sha256"] == REVIEWER_SHA
        and review.get("request", {}).get("sha256") == REQUEST_SHA, "Changed terminal review source/request")
    refs["reviewer_source"], refs["request"] = review["source"], review["request"]
    components = {}
    for key, ref in {**review["reviews"], "scheduler": review["scheduler"], "resource_replay": review["resource_replay"]}.items():
        refs["review_" + key] = ref
        components[key] = read(ref)
    native_root = ROOT / "benchmarks/results/native_factorial_cost_v2_20261004/run_12"
    for key, relative in (("native_done", "measurement/done.json"), ("native_receipt", "native_execution.json")):
        refs[key] = record(native_root / relative)
        components[key] = read(refs[key])
    components["missing_final_outputs"] = {name: not (native_root / name).exists()
        and not (native_root / name).is_symlink() for name in ("metrics.json", "native/orthohmm_orthogroups.txt")}
    require(docs["native11"]["references"]["manuscript"] == refs["parent"], "Mixed manuscript lineage")
    for ref in refs.values():
        check(ref)
    return docs, refs, components


def section(rows, refs, directory):
    failed = rows[6]
    resources = failed["resources"]
    lines = ["### Terminal Native QfO Outcomes", "",
        "The later terminal reviews do not add any admitted scores: four of seven",
        "fresh cells retain admitted accuracy and three retain missing scores.",
        "The preceding four-cell endpoint table, figure, intervals and localization",
        "remain unchanged. This update is failure reporting, not a complete factorial,",
        "independent confirmation or superiority over OrthoFinder.", "",
        "| Cell | Accuracy status | Six-metric mean |", "| --- | --- | ---: |"]
    for row in rows:
        mean = "Unavailable" if row["secondary_mean"] is None else f"{row['secondary_mean']:.8f}"
        lines.append(f"| {row['cell']} | {row['status']} | {mean} |")
    lines.extend(["", "Native11 (P1/C1/R0) completed inference and group-clique conversion,",
        "retaining 11,755,521 pairs and relation coverage 585,180/984,137.",
        "Its assessment ended OUT_OF_MEMORY at 8 CPU slots/32 GiB; FAS exited137.",
        "Five other endpoint tasks completed but supply no admitted score or mean.",
        "Coverage is not accuracy; its exact scoring peak memory remains unknown.", "",
        "Native12 (P1/C1/R1; job24036) ended FAILED1:0. The observer native-step",
        "wrapper completed0:0, but the actual native command exited-11 (SIGSEGV).",
        "The final native receipt remained running and final metrics/orthogroup",
        "outputs were absent. Retained downstream phylogeny files are partial",
        "evidence, not successful inference or admitted pair predictions. Conversion",
        "and scoring were not run. No stack trace was captured; the root cause and",
        "crash location are unknown. The last buffered progress message does not",
        "identify the crashing operation. No automatic inference retry was performed.", "",
        f"Reviewed failed-command observations: {resources['wall_seconds']:.6f} wall seconds,",
        f"{resources['cpu_seconds']:.6f} native-task CPU seconds and {resources['peak_memory_bytes']} bytes",
        "native-step lifetime peak memory, with the documented wrapper/launcher scopes.",
        "These describe a failed attempt, not successful inference or scoring cost.", "", DISCLOSURE, "",
        "Missing outcomes remain null, not zero. The three supported SwissTrees",
        "contrasts and all 42 planned adjusted endpoints retain their original scope;",
        "no new bootstrap draws, endpoints, functional similarity claims or default",
        "selection follow from this update. Prepared success-only five-cell reporting",
        "was not executed. This revision remains not submission-ready.", "",
        "Evidence:", *[f"- [{key}]({Path(os.path.relpath(ref['path'], directory)).as_posix()}) (SHA256 `{ref['sha256']}`)"
            for key, ref in refs.items() if key in ("snapshot", "native11", "review")], ""])
    return "\n".join(lines)


def export(review_path, review_sha, destination):
    destination = Path(destination).absolute()
    require(destination.is_relative_to(ROOT) and destination.resolve() == destination
        and not destination.exists() and not destination.is_symlink(), "Require fresh direct reporting destination")
    docs, refs, components = inputs(review_path, review_sha)
    rows = validate(docs["snapshot"], docs["native11"], docs["review"], components)
    accounting = reviewer_accounting()
    require(docs["parent"].count(ANCHOR) == 1 and docs["parent"].count("## Abstract\n") == 1,
        "Changed manuscript insertion anchor")
    # Keep inherited relative links beside the unchanged parent manuscript.
    manuscript_path = Path(refs["parent"]["path"]).parent / (destination.name + "_manuscript.md")
    require(not manuscript_path.exists() and not manuscript_path.is_symlink(), "Require fresh manuscript path")
    source = record(__file__)
    helpers = [record(Path(__file__).with_name(name)) for name in (
        "prepare_ob_candidate_neighborhood.py", "probe_dgx_step_separation.py",
        "review_native_factorial_attempt.py", "derive_threadripper_resources.py",
        "run_native_factorial_cost.py", "validate_native_factorial_outputs.py")]
    destination.mkdir(parents=True, exist_ok=False)
    with (destination / "scores.tsv").open("x") as stream:
        writer = csv.writer(stream, delimiter="\t", lineterminator="\n")
        writer.writerow(["Cell", "Endpoint", "Statistic", "Value", "Accuracy Admitted"])
        for row in rows:
            for endpoint in ENDPOINTS:
                statistic = "F1" if endpoint in ENDPOINTS[:3] else "Sample mean" if endpoint == "FAS" else "Similarity"
                value = row["scores"][endpoint]
                writer.writerow([row["cell"], endpoint, statistic, "Unavailable" if value is None else value, row["accuracy_admitted"]])
    with (destination / "status.tsv").open("x") as stream:
        writer = csv.writer(stream, delimiter="\t", lineterminator="\n")
        writer.writerow(["Cell", "Status", "Accuracy Admitted", "Secondary Mean", "Inference", "Conversion", "Scoring"])
        for row in rows:
            writer.writerow([row["cell"], row["status"], row["accuracy_admitted"],
                "Unavailable" if row["secondary_mean"] is None else row["secondary_mean"],
                *[row.get(key, "See frozen snapshot") for key in ("inference_status", "conversion_status", "scoring_status")]])
    text = section(rows, refs, manuscript_path.parent)
    (destination / "addendum.md").write_text("# Native QfO Terminal-Failure Addendum\n\n" + section(rows, refs, destination))
    body = docs["parent"][docs["parent"].index("## Abstract\n"):]
    body = body.replace(ANCHOR, text + "\n" + ANCHOR, 1)
    manuscript = ("# OrthoHMM: HMM-Centered Group Inference With Phylogenetic Refinement\n\n"
        "Terminal-failure evidence revision of frozen v4. Not submission-ready.\n"
        "Four admitted native QfO cells are unchanged; later inference/scoring\n"
        "failures supply no score. This manuscript source has not yet received a\n"
        "new render, manual page/citation review or archival inclusion.\n"
        f"Preserves [frozen v4]({Path(refs['parent']['path']).name}); "
        f"[terminal-failure reporting]({Path(os.path.relpath(destination / 'report.json', manuscript_path.parent)).as_posix()}).\n\n") + body
    manuscript_path.write_text(manuscript)
    for ref in [source, *helpers, *refs.values()]:
        check(ref)
    report = dict(schema="native_qfo_terminal_failure_reporting_v1", source=source, helpers=helpers,
        references=refs, independent_reviewer_accounting=accounting, rows=rows,
        missing_final_output_observation=components["missing_final_outputs"],
        admitted_cells=4, missing_cells=3, admitted_endpoints=24, missing_endpoints=18,
        outputs=[*[record(destination / name) for name in ("scores.tsv", "status.tsv", "addendum.md")],
            record(manuscript_path)],
        new_scientific_admission=False, automatic_retry=False, native_inference_reexecuted=False,
        new_bootstrap_draws=0, frozen_manuscript_replaced=False, publication_ready=False,
        limitations=["Direct retained review metadata/reporting, not a new raw scientific admission or raw-resource replay.",
            "Native12 SIGSEGV root cause and crash location unknown; partial phylogeny outputs not scored.",
            "Native11 scoring OOM supplies no partial scores or six-metric mean.",
            "Four-cell uncertainty and figures unchanged; incomplete factorial, development-exposed evidence.",
            "New manuscript render, manual review and archive integration remain unverified."])
    save(destination / "report.json", report)
    return record(destination / "report.json")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--review", type=Path, required=True)
    parser.add_argument("--review-sha256", required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    print(json.dumps(export(args.review, args.review_sha256, args.output), sort_keys=True))
