"""Report supplied native factorial reviews without executing or admitting runs."""

import argparse
import csv
import hashlib
import io
import json
import math
from pathlib import Path


SCOPES = {
    "wall_seconds": "native_command_monotonic_interval",
    "cpu_seconds": "native_task_subtree_cpu_stat_bracket_including_wrapper",
    "peak_memory_bytes": "native_step_lifetime_memory_peak_including_launcher",
}
DISCLOSURE = (
    "Timing measurements were collected on a shared Threadripper while other "
    "analyses were running. Competition for CPU, memory bandwidth and I/O may "
    "have affected elapsed times, with an unknown and potentially tool-dependent "
    "impact. These are observed shared-host timings, not estimates of isolated performance."
)
FIELDS = ["index", "dataset", "cell", "job_id", "outcome", "f1", "precision",
          "recall", "reference_gene_coverage", "partition_equal", "wall_seconds",
          "cpu_seconds", "peak_memory_bytes", "maximum_foreign_average_cores"]


def require(condition, message):
    if not condition:
        raise ValueError(message)


def load(path, sha256, evidence):
    path = Path(path).absolute()
    require(path.is_file() and not path.is_symlink(), "Nonregular report: " + str(path))
    data = path.read_bytes()
    ref = {"path": str(path), "bytes": len(data),
           "sha256": hashlib.sha256(data).hexdigest()}
    require(ref["sha256"] == sha256, "Report checksum mismatch: " + str(path))
    evidence.append(ref)
    return json.loads(data), ref


def same_bytes(left, right):
    return all(k in left and k in right and left[k] == right[k]
               for k in ("bytes", "sha256"))


def finite(value, name, low=0, high=math.inf):
    require(type(value) in (int, float) and math.isfinite(value) and low <= value <= high,
            "Invalid numeric field: " + name)
    return value


def score_row(review, review_ref, score, score_ref, plan, plan_ref, evidence):
    require(same_bytes(score["terminal_review"], review_ref), "Score/review binding differs")
    for key in ("index", "cell", "job_id"):
        require(score[key] == review[key], "Score identity differs: " + key)
    require(score["native_outputs_validated"] is True
            and score["accuracy_evaluated"] is True, "Unvalidated score")
    success = review["status"] == "native_success"
    if success:
        require(score["schema"] == "native_factorial_orthobench_score_v1"
                and score["status"] == "terminal_native_orthobench_scored"
                and score["dataset"] == "OrthoBench"
                and same_bytes(score["plan"], plan_ref), "Wrong successful score schema/plan")
        require(score["resource_observation"] == review["resources"]
                and score["resource_scopes"] == SCOPES, "Score resource binding differs")
        data, partition = score, score["canonical_partition_comparison"]
        count = score["score_percent"]["refogs"]
        covered, total = score["covered_reference_genes"], score["reference_genes"]
        require(score["reference_families"] == count, "Reference family count differs")
    else:
        adoption = plan.get("history_adoption", {})
        require(score["schema"] == "native_factorial_receipt_failure_recovery_v1"
                and score["status"] == "scientific_outputs_recovered_from_failed_wrapper"
                and score["scheduler_success"] is False
                and score["native_command_success"] is False
                and score["timing_success_established"] is False
                and same_bytes(adoption["recovery"], score_ref)
                and same_bytes(adoption["prior_review"], review_ref), "Unbound failed recovery")
        require(score["resources"] == review["resources"]
                and score["resource_scopes"] == SCOPES, "Recovery resource binding differs")
        ref = score["score"]
        data, observed = load(ref["path"], ref["sha256"], evidence)
        require(observed == ref and data["schema"] == "recovered_native_orthobench_score_v1"
                and data["cell"] == review["cell"] and data["dataset"] == "OrthoBench"
                and data["timing_success_established"] is False, "Wrong recovered score")
        partition = data["original_partition_comparison"]
        count = data["score"]["refogs"]
        # The recovery has no explicit reference-coverage numerator; do not infer one.
        covered, total = None, None
    require(count == 70 and data["development_exposed"] is True
            and data["independent_validation"] is False, "Wrong accuracy scope")
    metrics = data["score_fraction"]
    values = {key: finite(metrics[key], key, high=1)
              for key in ("f_score", "precision", "recall")}
    require(type(partition["partition_equal"]) is bool, "Invalid partition comparison")
    coverage = None
    if covered is not None:
        require(type(covered) is int and type(total) is int and 0 <= covered <= total
                and total > 0, "Invalid reference coverage")
        coverage = covered / total
    return {"f1": values["f_score"], "precision": values["precision"],
            "recall": values["recall"], "reference_gene_coverage": coverage,
            "partition_equal": partition["partition_equal"]}


def collect(plan_path, plan_sha256, attempts):
    evidence = []
    plan, plan_ref = load(plan_path, plan_sha256, evidence)
    runs = plan["runs"]
    require(len(runs) == 13 and [r["index"] for r in runs] == list(range(13)),
            "Wrong thirteen-identity plan")
    require(len({(r["dataset"], r["cell"], r["repeat"]) for r in runs}) == 13,
            "Duplicate planned identity")
    rows = [{**{k: None for k in FIELDS},
             **{k: run[k] for k in ("index", "dataset", "cell")},
             "outcome": "no_supplied_terminal_review"} for run in runs]
    seen = set()
    for review_path, review_sha, score_path, score_sha in attempts:
        review, review_ref = load(review_path, review_sha, evidence)
        index = review["index"]
        require(type(index) is int and index in range(13) and index not in seen,
                "Invalid or duplicate attempt index")
        seen.add(index)
        run = runs[index]
        require(review["schema"] == "native_factorial_terminal_review_v1"
                and all(review[k] == run[k] for k in ("index", "dataset", "cell", "repeat"))
                and review["terminal_reviewed"] is True
                and review["primary_resources_replayed"] is True
                and review["shared_host_resources_reviewed"] is True
                and review["execution_scope"] == "shared_host_matched_resources"
                and review["uncontended_timing"] is False
                and review["resource_scopes"] == SCOPES, "Invalid terminal review scope")
        require(same_bytes(review["plan"], plan_ref)
                or (same_bytes(review_ref, plan.get("history_adoption", {}).get("prior_review", {}))
                    and same_bytes(review["plan"], plan["history_adoption"]["prior_plan"])),
                "Terminal review is outside plan/history adoption")
        success = review["status"] == "native_success"
        require((success and review["scheduler_state"] == "COMPLETED"
                 and review["scheduler_exit_code"] == "0:0"
                 and review["native_outputs_validated"] is True)
                or (review["status"] == "native_failure_retained"
                    and review["scheduler_state"] in ("FAILED", "TIMEOUT", "CANCELLED", "OUT_OF_MEMORY")),
                "Inconsistent terminal outcome")
        resources = review["resources"]
        require(not success or resources is not None, "Successful native run has no resources")
        if resources is not None:
            require(set(resources) == set(SCOPES), "Wrong resource endpoints")
            for key, value in resources.items():
                finite(value, key)
            require(type(resources["peak_memory_bytes"]) is int, "Noninteger memory bytes")
        row = rows[index]
        foreign = review["whole_run_maximum_foreign_average_cores"]
        if foreign is not None:
            finite(foreign, "foreign CPU")
        require(resources is None or foreign is not None, "Measured run lacks contention observation")
        row.update({"job_id": review["job_id"], "outcome": review["status"],
                    "maximum_foreign_average_cores": foreign})
        if resources:
            row.update(resources)
        require((score_path == "-") == (score_sha == "-"), "Incomplete optional score arguments")
        if score_path != "-":
            require(run["dataset"] == "orthobench", "QfO scores require native endpoint reporting")
            score, score_ref = load(score_path, score_sha, evidence)
            row.update(score_row(review, review_ref, score, score_ref, plan, plan_ref, evidence))
            if not success:
                row["outcome"] = "failed_wrapper_science_recovered"
    for ref in evidence:
        data = Path(ref["path"]).read_bytes()
        require(len(data) == ref["bytes"] and hashlib.sha256(data).hexdigest() == ref["sha256"],
                "Report changed during export")
    return {"schema": "native_factorial_reporting_snapshot_v1", "plan": plan_ref,
            "evidence": evidence, "rows": rows, "resource_scopes": SCOPES,
            "timing_disclosure": DISCLOSURE, "publication_ready": False,
            "new_scoring_or_admission": False,
            "limitations": [
                "Only explicitly supplied terminal reviews are summarized; blanks are not live-job status or zeros.",
                "Recovered failed-wrapper resources are failed-attempt observations, not clean-success timings.",
                "Missing accuracy is not inherited from cached predictions; QfO metrics require separate native admission.",
                "This checks direct report identities, not transitive raw artifacts or new independent accuracy.",
                "No averages, causal component overhead, isolated ranking, retries or job release.",
                "The three reused configurations and historical cached stage costs remain separately reported."]}


def markdown(report):
    lines = ["# Full-Native Factorial Progress", "", report["timing_disclosure"], "",
             "P toggles downstream profiles, C candidate expansion, R reconciliation. Initial HMM search remains on.",
             "OrthoBench F1 is reference-group co-membership, not native resolved-pair accuracy.", "",
             "| Index | Dataset | Cell | Outcome | F1 (%) | Wall (s) | CPU (s) | Peak (GiB) |",
             "| ---: | --- | --- | --- | ---: | ---: | ---: | ---: |"]
    def display(value, multiplier=1):
        return "Unavailable" if value is None else f"{value * multiplier:.4f}"
    for row in report["rows"]:
        lines.append(f"| {row['index']} | {row['dataset']} | {row['cell']} | {row['outcome']} | "
                     f"{display(row['f1'], 100)} | {display(row['wall_seconds'])} | "
                     f"{display(row['cpu_seconds'])} | {display(row['peak_memory_bytes'], 1 / 2**30)} |")
    lines += ["", "CPU includes the wrapper; peak includes the native-step launcher, not pure algorithm RSS.",
              "Preparation, conversion and scoring are outside the native interval.", ""]
    lines += ["- " + item for item in report["limitations"]]
    return "\n".join(lines) + "\n"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--plan", required=True)
    parser.add_argument("--plan-sha256", required=True)
    parser.add_argument("--attempt", nargs=4, action="append", default=[],
                        metavar=("REVIEW", "REVIEW_SHA256", "SCORE_OR_RECOVERY", "SCORE_SHA256"))
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    output = Path(args.output).absolute()
    require(not output.exists() and not output.is_symlink(), "Output already exists")
    report = collect(args.plan, args.plan_sha256, args.attempt)
    source = Path(__file__).resolve()
    report["source"] = {"path": str(source), "bytes": source.stat().st_size,
                        "sha256": hashlib.sha256(source.read_bytes()).hexdigest()}
    table = io.StringIO()
    writer = csv.DictWriter(table, fieldnames=FIELDS, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(report["rows"])
    output.mkdir(parents=True, exist_ok=False)
    (output / "report.json").write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    (output / "rows.tsv").write_text(table.getvalue())
    (output / "report.md").write_text(markdown(report))
    print(json.dumps({"output": str(output), "supplied_reviews": len(args.attempt),
                      "identities": len(report["rows"])}))


if __name__ == "__main__":
    main()
