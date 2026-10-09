"""Consolidate terminal native observations without new timing or scientific admission."""

import argparse
import copy
import csv
import io
import json
from pathlib import Path

from benchmark_tools.bundle_publication_package import identity, pin, save
from benchmark_tools.export_native_factorial_progress import SCOPES, DISCLOSURE, finite, load, require


SOURCES = {
    "baseline": ("benchmark_tools/results/native_factorial_progress_20261005_v7/report.json",
                 "7db3f188739f67fb254b6fb8b299153799a0cd6c18e20d6428bf234863c61648"),
    "qfo": ("benchmark_tools/results/native_qfo_terminal_failures_20261009_v1/report.json",
            "863ebc334dbdb5a41d1a9e9bba56be82896f0ff26c98366047928c7ef802966f"),
    7: ("benchmarks/work/native_factorial_measurement_failure_review_22437/review.json",
        "650b69e0965e0a667359efaa43c1d814108f613547cbe66cf12964bd1e64e459"),
    8: ("benchmarks/work/native_factorial_terminal_review_22444_fault_reported_v1/review.json",
        "1ddc14f79f0f501e4189382cdfe613a883310a6fceacbd94b37f342bb33827ef"),
    9: ("benchmarks/work/native_factorial_launch_failure_review_23891/review.json",
        "1f6167cf997f37dcff5f795c9479475bd065aa7359a68a73994f7f53abf4ed9a"),
    10: ("benchmarks/work/native10_review_finalization_20261007_v1/review/review.json",
         "7eab212d9788deb031232e4641089ff96503c0e5727e1ec369067f93ffb518b5"),
    11: ("benchmarks/work/native11_composed_terminal_review_20261008_v1/review.json",
         "a87d65439a3f65dd80b199602231e8f95d10f26d2dbbee54c2f0d0eaa50abf2e"),
    12: ("benchmarks/work/native12_composed_terminal_review_20261008_v1/review.json",
         "a55a6da5ffb8dfbb362883937c915b161796a9ddd7275f4ac48a1c4d049e5c4a"),
}
REVIEW_SCHEMAS = {
    7: "native_factorial_measurement_failure_review_v1",
    8: "native_factorial_terminal_review_v1",
    9: "native_factorial_launch_failure_review_v1",
    10: "allocated_native_factorial_terminal_review_v1",
    11: "native11_composed_terminal_review_v1",
    12: "native12_composed_terminal_review_v1",
}
FIELDS = ["index", "dataset", "cell", "job_id", "native_outcome", "accuracy_admitted",
          "accuracy_status", "resource_observation", "wall_seconds", "cpu_seconds",
          "peak_memory_bytes", "new_timing_admission"]


def collect(baseline, qfo, reviews):
    require(baseline["schema"] == "native_factorial_reporting_snapshot_v1"
            and baseline["resource_scopes"] == SCOPES
            and baseline["publication_ready"] is False, "Wrong baseline scope")
    require(qfo["schema"] == "native_qfo_terminal_failure_reporting_v1"
            and qfo["publication_ready"] is False
            and qfo["new_scientific_admission"] is False, "Wrong QfO snapshot scope")
    original = baseline["rows"]
    require([row["index"] for row in original] == list(range(13)), "Wrong thirteen-row inventory")
    newer = {row["index"]: row for row in qfo["rows"]}
    require(len(newer) == len(qfo["rows"]) == 7 and set(newer) == set(range(6, 13))
            and set(reviews) == set(range(7, 13)), "Missing or duplicate terminal identities")
    rows = []
    for old in original:
        index = old["index"]
        require(old["dataset"] == ("orthobench" if index < 6 else "qfo_corrected"), "Dataset differs")
        row = dict(index=index, dataset=old["dataset"], cell=old["cell"], job_id=old["job_id"],
                   native_outcome=old["outcome"], accuracy_admitted=old["f1"] is not None,
                   accuracy_status="retained_orthobench_score" if index < 6 else "not_joined",
                   resource_observation="reviewed_native_command", new_timing_admission=False)
        resources = {key: old[key] for key in SCOPES}
        if index == 0:
            row["resource_observation"] = "failed_wrapper_observation_not_clean_success"
        if index >= 6:
            score = newer[index]
            require(score["cell"] == old["cell"] and type(score["accuracy_admitted"]) is bool,
                    "Accuracy identity differs")
            require(score["accuracy_admitted"] == all(value is not None for value in score["scores"].values())
                    and len(score["scores"]) == 6, "Missing scores cannot be admitted")
            row.update(accuracy_admitted=score["accuracy_admitted"], accuracy_status=score["status"])
            if index == 6:
                require(score["native_job_id"] == old["job_id"], "Retained native job differs")
            else:
                review = reviews[index]
                require(review["schema"] == REVIEW_SCHEMAS[index]
                        and review["index"] == index and review["cell"] == old["cell"]
                        and review["dataset"] == old["dataset"]
                        and review["terminal_reviewed"] is True
                        and review["execution_scope"] == "shared_host_matched_resources"
                        and review["uncontended_timing"] is False
                        and review["scientific_timings_admitted"] is False, "Wrong terminal review scope")
                require(score["native_job_id"] in (None, review["job_id"]), "Score/review job differs")
                resources = review["resources"]
                row.update(job_id=review["job_id"], native_outcome=review["status"])
                if resources is None:
                    require(index in (7, 9) and review["scheduler_state"] == "FAILED"
                            and review["primary_resources_replayed"] is False
                            and review.get("resource_scopes") is None, "Unexpected missing measurements")
                    require(review["status"] == (
                        "native_scientific_outputs_recovered_measurement_failure_retained" if index == 7
                        else "pre_native_cpu_binding_failure_reviewed_retained"), "Wrong missing-resource outcome")
                    row["resource_observation"] = (
                        "measurement_failure_no_valid_resources" if index == 7 else "pre_native_failure_no_resources")
                    require(score.get("resources") is None, "Failed measurement has supplied resources")
                else:
                    require(review.get("resource_scopes") == SCOPES
                            and review["primary_resources_replayed"] is True
                            and review["shared_host_resources_reviewed"] is True, "Unreviewed resource scope")
                    if index >= 10:
                        require(score["resources"] == resources and score["resource_scopes"] == SCOPES,
                                "Snapshot/review resources differ")
                    if index == 12:
                        require(review["status"] == "native_failure_retained"
                                and review["scheduler_state"] == "FAILED"
                                and review["native_outputs_validated"] is False
                                and score["accuracy_admitted"] is False, "Failed inference cannot supply accuracy")
                        row["resource_observation"] = "failed_native_command_not_successful_inference"
                    else:
                        require(review["status"] == "native_success"
                                and review["scheduler_state"] == "COMPLETED"
                                and review["native_outputs_validated"] is True, "Inconsistent native success")
        if resources is None:
            resources = {key: None for key in SCOPES}
        require(set(resources) == set(SCOPES), "Wrong resource endpoints")
        require(all(value is None for value in resources.values())
                or all(value is not None for value in resources.values()), "Partial resource vector")
        for key, value in resources.items():
            if value is not None:
                finite(value, key)
                require(key != "peak_memory_bytes" or type(value) is int, "Noninteger peak bytes")
        row.update(resources)
        rows.append(row)
    return rows


def run(root, output):
    root, output = Path(root).resolve(), Path(output).absolute()
    if output.exists() or output.is_symlink():
        raise FileExistsError(output)
    evidence, values = [], {}
    for name, (path, sha) in SOURCES.items():
        values[name], _ = load(root / path, sha, evidence)
    helper = Path(load.__code__.co_filename).resolve()
    require(identity(helper) == pin(values["baseline"]["source"]), "Frozen prior reporter helper changed")
    rows = collect(values["baseline"], values["qfo"], {key: values[key] for key in REVIEW_SCHEMAS})
    table = io.StringIO()
    writer = csv.DictWriter(table, fieldnames=FIELDS, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows({key: "NA" if row[key] is None else row[key] for key in FIELDS} for row in rows)
    markdown = ["# Terminal Native Ablation Resource Observations", "", DISCLOSURE, "",
        "Native inference only; excludes preparation, conversion and scoring. CPU includes the task wrapper;",
        "peak is native-step lifetime memory including its launcher, not pure algorithm RSS.", "",
        "| Dataset | Cell | Job | Native Outcome | Accuracy Admitted | Resource Status | Wall (s) | CPU (s) | Peak (GiB) |",
        "| --- | --- | --- | --- | --- | --- | ---: | ---: | ---: |"]
    for row in rows:
        numbers = ["NA" if row[key] is None else f"{row[key] / (2**30 if key == 'peak_memory_bytes' else 1):.4f}"
                   for key in SCOPES]
        markdown.append("| " + " | ".join(str(row[key]) for key in (
            "dataset", "cell", "job_id", "native_outcome", "accuracy_admitted", "resource_observation"))
            + " | " + " | ".join(numbers) + " |")
    markdown.extend(["", "No new timing/accuracy admission, failed-attempt retry or isolated-speed ranking.",
        "Index 7 has recovered accuracy but no valid resources; index 9 failed before native inference.",
        "Index 11 has native observations but scoring OOM; index 12 has failed-command resources only.",
        "Original cached factorial full costs and separate historical configuration associations remain unchanged.", ""])
    output.mkdir(parents=True)
    (output / "resources.tsv").write_text(table.getvalue())
    (output / "resources.md").write_text("\n".join(markdown))
    result = dict(schema="native_factorial_terminal_resource_snapshot_v1", rows=rows, evidence=evidence,
        source=dict(path=str(Path(__file__).resolve()), **identity(Path(__file__))),
        helpers=[dict(path=str(helper), **identity(helper)),
                 dict(path=str(Path(save.__code__.co_filename).resolve()),
                      **identity(Path(save.__code__.co_filename)))],
        baseline_rows_preserved=copy.deepcopy(values["baseline"]["rows"]), resource_scopes=SCOPES,
        timing_disclosure=DISCLOSURE, publication_ready=False, native_inference_repeated=False,
        scoring_repeated=False, new_scientific_admission=False, new_timing_admission=False,
        limitations=["Direct retained-report projection, not fresh raw accounting or scientific admission.",
                     "Failed/missing observations are retained; no comparable isolated efficiency is established.",
                     "Separate historical configuration costs are not original cached factorial full costs."])
    save(output / "report.json", result)
    return dict(rows=len(rows), measured_rows=sum(row["wall_seconds"] is not None for row in rows),
                report=identity(output / "report.json"), publication_ready=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    print(json.dumps(run(args.root, args.output)))
