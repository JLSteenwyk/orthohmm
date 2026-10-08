"""Report the retained scoring failure separately from frozen manuscript v4."""

import argparse
import csv
import json
import math
from pathlib import Path

from benchmark_tools.derive_threadripper_resources import SCOPES
from benchmark_tools.prepare_ob_candidate_neighborhood import record, check
from benchmark_tools.probe_dgx_step_separation import save
from benchmark_tools.review_native_factorial_attempt import DISCLOSURE
from benchmark_tools.run_native_factorial_cost import ROOT, read
from benchmark_tools.validate_native_factorial_outputs import require


INPUTS = {
    "failure": ("benchmark_tools/results/native11_composed_assessment_failure_24038_20261008_v1.json",
        "b40f8075fd0168e3e7e1b895bcdba2576e2b903cc58e3484a7f3ce8b53b93e6b"),
    "review": ("benchmarks/work/native11_composed_terminal_review_20261008_v1/review.json",
        "a87d65439a3f65dd80b199602231e8f95d10f26d2dbbee54c2f0d0eaa50abf2e"),
    "conversion": ("benchmarks/results/native11_composed_qfo_pairs_v1/results.json",
        "77d33a77eb5eb1c85736ba6e7b303aa269907c8a930d74647d3b8bddff36fac8"),
    "manuscript": ("benchmark_tools/results/PUBLICATION_MAIN_TEXT_20261007_v4.md",
        "0b7012e3dd46bd3b028857bdb3cbe3a9e9048936b90f3df9e4833cb342d5705a"),
}


def validate_reports(failure, review, conversion):
    require(failure.get("schema") == "native11_composed_assessment_failure_readback_v1"
        and failure.get("status") == "failed_scoring_retained_not_admitted"
        and failure.get("job_id") == 24038 and failure.get("native_job_id") == 23985
        and failure.get("native_index") == 11 and failure.get("cell") == "p1_c1_r0"
        and failure.get("output_inventory_reproduced") is True
        and all(failure.get(key) is False for key in ("accuracy_admitted", "admission_submitted",
            "automatic_retry", "native_inference_reexecuted", "publication_ready")),
        "Require retained actual scoring failure, not accuracy admission")
    rows = failure["accounting"]["rows"]
    parent = [row for row in rows if row.get("JobIDRaw") == "24038"]
    require(len(parent) == 1 and parent[0].get("State") == "OUT_OF_MEMORY"
        and parent[0].get("ExitCode") == "0:125" and parent[0].get("ReqMem") == "32G"
        and parent[0].get("AllocCPUS") == "8", "Failure accounting differs")
    failed = failure["failed_tasks"]
    tasks = failure["trace_tasks"]
    require(len(failed) == 1 and failed[0].get("name") == "fas_benchmark (1)"
        and failed[0].get("status") == "FAILED" and failed[0].get("exit") == "137"
        and [row for row in tasks if row["status"] != "COMPLETED"] == failed
        and failure.get("execution_status") == "failed" and failure.get("execution_exit_code") == 1,
        "Failed FAS task or scoring outcome differs")
    metrics = [row for row in tasks if row["name"].startswith(("vgnc_benchmark", "reference_genetrees_benchmark",
        "ec_benchmark", "go_benchmark"))]
    require(len(metrics) == 5 and all(row["status"] == "COMPLETED" and row["exit"] == "0" for row in metrics),
        "Require exact retained five completed endpoint tasks, not five admitted scores")
    require(review.get("schema") == "native11_composed_terminal_review_v1"
        and review.get("status") == "native_success" and review.get("index") == 11
        and review.get("job_id") == 23985 and review.get("cell") == "p1_c1_r0"
        and all(review.get(key) is True for key in ("terminal_reviewed", "composed_full_review_complete",
            "native_outputs_validated", "primary_resources_replayed", "shared_host_resources_reviewed"))
        and review.get("resource_scopes") == SCOPES
        and review.get("accuracy_evaluated") is False and review.get("uncontended_timing") is False,
        "Require successful reviewed inference separate from failed scoring")
    resources = review["resources"]
    require(set(resources) == set(SCOPES) and all(type(value) in (int, float)
        and math.isfinite(value) and value > 0 for value in resources.values())
        and type(resources["peak_memory_bytes"]) is int, "Invalid scoped inference resource observations")
    require(conversion.get("schema") == "native11_composed_qfo_conversion_v1"
        and conversion.get("status") == "native11_composed_qfo_pairs_prepared_unscored"
        and conversion.get("native_index") == 11 and conversion.get("native_job_id") == 23985
        and conversion.get("job_id") == "24035" and conversion.get("cell") == "p1_c1_r0"
        and conversion.get("conversion_kind") == "group"
        and conversion.get("semantics") == "cross-species group-derived clique pairs"
        and conversion.get("accuracy_evaluated") is False and conversion.get("native_inference_reexecuted") is False
        and conversion.get("removed_mapping_pairs") == 0
        and type(conversion.get("retained_pairs")) is int
        and conversion["retained_pairs"] == conversion.get("total_pairs") == conversion.get("expected_pairs"),
        "Require successful unscored group conversion")
    coverage = conversion["pair_coverage"]
    require(type(coverage.get("input_accessions")) is int and coverage["input_accessions"] == 984137
        and type(coverage.get("accessions_in_any_pair")) is int
        and 0 <= coverage["accessions_in_any_pair"] <= coverage["input_accessions"]
        and coverage.get("pair_rows") == conversion["retained_pairs"]
        and coverage.get("fraction_inputs_in_any_pair") == coverage["accessions_in_any_pair"] / coverage["input_accessions"],
        "Invalid full-input relation coverage")
    return dict(native_index=11, cell="p1_c1_r0", native_job_id=23985, conversion_job_id=24035,
        assessment_job_id=24038, inference_status="successful_composed_terminal_review",
        conversion_status="successful_unscored_group_cliques", scoring_status="OUT_OF_MEMORY",
        failed_endpoint="FAS", failed_endpoint_exit_code=137, completed_endpoint_task_count=5,
        admitted_endpoint_score_count=0, secondary_six_metric_mean=None, accuracy_admitted=False,
        submitted_pairs=conversion["retained_pairs"], input_accessions=coverage["input_accessions"],
        relation_accessions=coverage["accessions_in_any_pair"], relation_coverage=coverage["fraction_inputs_in_any_pair"],
        prediction_semantics=conversion["semantics"], resources=resources, resource_scopes=SCOPES,
        scoring_cpu_slots=8, scoring_memory_limit_bytes=32 * 1024 ** 3,
        scoring_peak_memory_bytes=None, scoring_scheduler_elapsed=parent[0]["Elapsed"])


def inputs():
    refs, values = {}, {}
    for name, (relative, digest) in INPUTS.items():
        ref = record(ROOT / relative)
        require(ref["sha256"] == digest, "Frozen addendum input changed: " + name)
        check(ref)
        refs[name] = ref
        if name != "manuscript":
            values[name] = read(ref)
    require(values["conversion"]["terminal_review"] == refs["review"], "Conversion/full review reference differs")
    require(values["failure"]["references"]["pairs_manifest"] == refs["conversion"], "Failure/conversion reference differs")
    return refs, values


def export(destination):
    destination = Path(destination).absolute()
    require(destination.is_relative_to(ROOT) and destination.resolve() == destination
        and not destination.exists() and not destination.is_symlink(), "Require new direct addendum destination")
    source = record(__file__)
    refs, values = inputs()
    row = validate_reports(values["failure"], values["review"], values["conversion"])
    destination.mkdir(parents=True, exist_ok=False)
    headers = ["Cell", "Inference", "Conversion", "Scoring", "Admitted Score", "Six-Metric Mean",
        "Submitted Pairs", "Relation Coverage", "Scoring CPU Slots", "Scoring Memory Limit GiB"]
    data = [row["cell"], row["inference_status"], row["conversion_status"], row["scoring_status"],
        "Unavailable", "Unavailable", row["submitted_pairs"], row["relation_coverage"],
        row["scoring_cpu_slots"], row["scoring_memory_limit_bytes"] / 1024 ** 3]
    with (destination / "status.tsv").open("x") as stream:
        writer = csv.writer(stream, delimiter="\t", lineterminator="\n")
        writer.writerow(headers)
        writer.writerow(data)
    r = row["resources"]
    text = ["# Native QfO Scoring Failure Addendum", "",
        "This new evidence update supplements the frozen manuscript v4; it does not replace its bytes or four admitted-cell scores.", "",
        "| " + " | ".join(headers) + " |", "| " + " | ".join(["---"] * len(headers)) + " |",
        "| " + " | ".join(str(value) for value in data) + " |", "",
        "Inference job 23985 completed and passed its explicit full composed review. Conversion job 24035 succeeded. "
        "Assessment job 24038 ended OUT_OF_MEMORY under 8 CPUs / 32 GiB. Its FAS task exited 137; Slurm reported an OOM kill. "
        "Five other endpoint tasks completed, but the six-endpoint assessment did not complete. No score or secondary mean is admitted.", "",
        f"Conversion retained {row['submitted_pairs']:,} cross-species group-clique pairs, with no mapping loss. "
        f"Relation coverage is {row['relation_accessions']:,}/{row['input_accessions']:,} "
        f"({100 * row['relation_coverage']:.4f}%). This is coverage, not accuracy, and not native phylogenetic-pair semantics.", "",
        f"Reviewed inference observations: {r['wall_seconds']:.6f} wall seconds, {r['cpu_seconds']:.6f} native-task CPU seconds, "
        f"{r['peak_memory_bytes']} bytes native-step lifetime peak memory, including the documented wrapper/launcher scopes. "
        "These are not scoring measurements. The scoring scheduler elapsed time is " + row["scoring_scheduler_elapsed"] +
        "; its exact peak RSS and the failing allocation are unavailable.", "", DISCLOSURE, "",
        "The failed scoring attempt is retained without retry, missing-score imputation or partial six-metric mean. "
        "A 128-GiB limit is prospectively specified for the genuinely unrun native12 scoring workflow only; "
        "it does not change inference, endpoint definitions or FAS sampling and does not guarantee success.", "",
        "Evidence:", *[f"- {name}: `{ref['path']}` (SHA256 `{ref['sha256']}`)" for name, ref in refs.items()], ""]
    with (destination / "addendum.md").open("x") as stream:
        stream.write("\n".join(text))
    for ref in [source, *refs.values()]:
        check(ref)
    report = dict(schema="composed_native_qfo_failure_reporting_v1", source=source, references=refs, row=row,
        outputs=[record(destination / name) for name in ("status.tsv", "addendum.md")],
        publication_ready=False, new_scientific_admission=False, native_inference_reexecuted=False,
        automatic_retry=False, frozen_manuscript_replaced=False,
        limitations=["Direct reporting of retained diagnostic/review/conversion evidence, not a new raw scientific admission.",
            "Five completed endpoint tasks do not supply a successful six-endpoint assessment.",
            "Coverage is not accuracy; failed scoring contributes no zero or six-metric mean.",
            "Historical four-cell scores and uncertainty remain unchanged; final native outcome is still pending."])
    save(destination / "report.json", report)
    return record(destination / "report.json")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    print(json.dumps(export(args.output), sort_keys=True))
