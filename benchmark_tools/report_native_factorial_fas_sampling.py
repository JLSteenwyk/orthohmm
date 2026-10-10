"""Report conditional design ranges for four retained native QfO FAS cells."""

from itertools import combinations
import json
import math
from pathlib import Path
import re
import traceback

import scipy

from benchmark_tools.audit_qfo_fas_samples import read_sample
from benchmark_tools import native_fas_two_sample_interval as kernel
from benchmark_tools.native_fas_sampling_interval import difference_interval, require
from benchmark_tools.prepare_ob_candidate_neighborhood import record, check

CELLS = ((6, "p0_c0_r0"), (7, "p0_c0_r1"), (8, "p0_c1_r0"), (10, "p1_c0_r1"))
LOGS = {
    6: ("qfo_benchmark/w/nq06/60/f411bd4592afe3d639882d25ba58a8/.command.log", "f568e601c2c54ebffc7b45e2045dfafe15935bbbd10f7b0fdb8f461f930e087d"),
    7: ("qfo_benchmark/w/mr07/0e/84a604de7a9c2efdd1194ed7b06dd3/.command.log", "f2206e988f2118cd8f106350e5a5cf0d2232a6dfc437c2bc4494a9b08bd8ea08"),
    8: ("qfo_benchmark/w/nq08/e3/7104306d6074c7d9753190515364e8/.command.log", "ac036f9cc32f622a0f672c58c0cf27c8f3e38917b5485a9071c133b58a765353"),
    10: ("qfo_benchmark/w/aq10/50/ae65cc225b9a9933d7bd4f132f0851/.command.log", "3b3af37c68c34eba2b4e00d2f1bbc7a81ff9d23551bb684ada8a3741b672cde5"),
}
PINS = {
    "snapshot": ("benchmark_tools/results/native_qfo_scientific_scores_20261007_v3/report.json", "7aab1cbd31fb650df42e6e80a14ee0167b35a810f08202768c89c514755211a2"),
    "protocol": ("benchmark_tools/results/NATIVE_FACTORIAL_FAS_SAMPLING_PROTOCOL_20261010.md", "ef29124c338194469d6f39fff6e03a8150cbae5820fedfd6f1b73280ea0c6be8"),
    "native_scorer": ("qfo_benchmark/benchmark-webservice/fas_benchmark.py", "1045c57f4d0f4787bec3d1f0690799df63c68dcc67ccfd5a925c338e3d33661d"),
    "count_kernel": ("benchmark_tools/native_fas_sampling_interval.py", "803ff0f455213e6a55b83e30c586ab1c57566830f8f877c916b78e84bf3658f9"),
    "two_sample_kernel": ("benchmark_tools/native_fas_two_sample_interval.py", "7acd70146adcd60672c76a733896820f3fe25b69c8408b54ed7db881555d3d21"),
    "sample_reader": ("benchmark_tools/audit_qfo_fas_samples.py", "e21f6865062685651fed6e16e30dc9935c53de08aca5cd96baf6dce7331deac2"),
    "context": ("benchmark_tools/results/native_fas_context_probe_20261010_v2.json", "845b4baea528f7633830e5a26ec58267da95165d8c4f79f2a0112c6063826f33"),
}
COMPONENT_ERROR = .05 / 12


def one(pattern, text):
    matches = re.findall(pattern, text)
    require(len(matches) == 1, "Missing or ambiguous native log field")
    return matches[0]


def sample_summary(text, pairs, row):
    require("participant='" + row["participant"] + "'" in text and "limited_species=False" in text,
            "Native participant/species scope differs")
    require("Computing fas.runMultiTaxa failed:" not in text, "Batch failure cannot admit fixed pair outcomes")
    P, M, no_features = map(int, one(r"(\d+) pairs precomputed, (\d+) missing \(will compute\); (\d+) no feature annotations", text))
    c, requested_k = map(int, one(r"we will compute (\d+) new pairs and sample (\d+) precomputed pairs", text))
    old_mean, _, old_n = one(r"FAS score\[precomputed\]: ([0-9.e+-]+) \+- ([0-9.e+-]+) \[N=(\d+)\]", text)
    new_mean, _, new_n = one(r"FAS score\[missing\]: ([0-9.e+-]+) \+- ([0-9.e+-]+) \[N=(\d+)\]", text)
    k, r = int(old_n), int(new_n)
    require(P > 0 and M > 0 and c == min(M, 9000), "Native eligible/count support differs")
    fraction = P / (P + M)
    require(requested_k == round(c * fraction / (1 - fraction)) and k == min(P, requested_k)
            and k > 0 and 0 < r <= c, "Native float rounding or returned counts differ")
    require(P + M == row["fas_sample"]["reported_eligible_pairs"] and len(pairs) == k + r
            and len(pairs) == row["fas_sample"]["sample_pairs"], "Admitted native sample counts differ")
    values = list(pairs.values())
    a, b = math.fsum(values[:k]) / k, math.fsum(values[k:]) / r
    require(abs(a - float(old_mean)) <= 5.00001e-7 and abs(b - float(new_mean)) <= 5.00001e-7,
            "Precomputed-first raw strata do not reproduce rounded native means")
    z = math.fsum(values) / len(values)
    log_z, eligible, size, size2 = one(r"FAS_mean: ([0-9.e+-]+) \+- [0-9.e+-]+; nr_orthologs: (\d+); sample_size: (\d+) vs (\d+)", text)
    require(int(eligible) == P + M and int(size) == int(size2) == k + r
            and abs(z - float(log_z)) <= 1e-12 and abs(z - row["scores"]["FAS"]) <= 1e-12,
            "Native raw/aggregate/log mean differs")
    return dict(P=P, M=M, k=k, c=c, r=r, omitted_new_numeric_scores=c-r,
                no_feature_annotation_pairs=no_features, precomputed_sample_mean=a,
                returned_sample_mean=b, observed_native_mean=row["scores"]["FAS"])


def read_cell(row, log_ref, scorer_ref, evidence):
    admission_ref = row["admission"]
    check(admission_ref)
    admission = json.loads(Path(admission_ref["path"]).read_text())
    require(admission["accuracy_admitted"] is True and admission["cell"] == row["cell"]
            and admission["native_index"] == row["index"] and admission["native_job_id"] == row["native_job_id"]
            and admission["participant"] == row["participant"]
            and admission["fas_protocol"]["source"] == scorer_ref, "Changed native scientific admission")
    execution_ref = admission["execution_report"]
    require(execution_ref in admission["checked_records"], "Execution not bound by admission")
    check(execution_ref)
    execution = json.loads(Path(execution_ref["path"]).read_text())
    require(execution["exit_code"] == 0 and execution["status"] == "process_succeeded_pending_independent_admission",
            "Native assessment did not succeed")
    raw = row["fas_sample"]["raw"]
    require(raw == admission["fas_sample"]["raw"] and raw in admission["checked_records"]
            and raw in execution["outputs"], "Raw FAS sample not bound to native execution/admission")
    tasks = [t for t in admission["native_tasks"] if t["name"].startswith("fas_benchmark")]
    require(len(tasks) == 1 and tasks[0]["status"] == "COMPLETED" and tasks[0]["exit"] == "0",
            "FAS native task is not uniquely successful")
    task = tasks[0]
    parent = Path(log_ref["path"]).parent
    require(parent.parent.parent == Path(execution["work"]) and
            (parent.parent.name + "/" + parent.name).startswith(task["hash"]), "Native task/log path differs")
    if row["index"] == 7:
        require(row["resources"] is None and row["timing_admitted"] is False and row["timing_eligible"] is False,
                "Recovered science must not repair failed timing")
    check(raw)
    check(log_ref)
    summary = sample_summary(Path(log_ref["path"]).read_text(), read_sample(Path(raw["path"])), row)
    evidence.extend([admission_ref, execution_ref, raw, log_ref])
    return dict(cell=row["cell"], index=row["index"], prediction_semantics=row["prediction_semantics"],
                raw=raw, native_log=log_ref, native_log_historically_hash_bound=log_ref in admission["checked_records"],
                failed_timing_preserved=row["index"] == 7, **summary)


def panel(rows):
    require([(r["index"], r["cell"]) for r in rows] == list(CELLS), "Changed four-cell inventory")
    methods = []
    for row in rows:
        result = dict(row)
        try:
            interval = kernel.method_interval(row["M"], row["c"], row["P"], row["k"],
                row["precomputed_sample_mean"], row["r"], row["returned_sample_mean"], COMPONENT_ERROR)
            result.update(interval=interval, status="conditional_design_range_computed", error=None)
        except Exception as exc:
            result.update(interval=None, status="conditional_numerical_failure",
                          error=dict(type=type(exc).__name__, message=str(exc), traceback=traceback.format_exc()))
        methods.append(result)
    contrasts = []
    for left, right in combinations(methods, 2):
        bounds = (difference_interval(left["interval"], right["interval"])
                  if left["interval"] and right["interval"] else None)
        contrasts.append(dict(left=left["cell"], right=right["cell"],
            observed_difference=left["observed_native_mean"] - right["observed_native_mean"],
            conditional_expected_difference_bounds=bounds,
            zero_included=bounds[0] <= 0 <= bounds[1] if bounds else None))
    complete = all(r["interval"] is not None for r in methods)
    return dict(schema="native_factorial_conditional_fas_sampling_v1", methods=methods, contrasts=contrasts,
                status="conditional_design_ranges_complete" if complete else "conditional_design_ranges_partial",
                alpha=.05, component_error=COMPONENT_ERROR, components=12, joint_error_bound=.05,
                conditional_design_ranges_computed=complete, historical_scores_rerun=False,
                observed_scores_changed=False, biological_generalization_intervals=False,
                unconditional_historical_interval_admission=False, other_endpoint_uncertainty_admitted=False,
                failed_factorial_cells_repaired=False, publication_ready=False,
                target="expected_native_post_attrition_ratio_under_fixed_design",
                limitations=[
                    "Conditional uniform sampling and fixed pair outcomes are assumptions, not historical or universal certification.",
                    "Both population means are unknown; observed Z is descriptive, not an exact expected-ratio point.",
                    "No biological pair/family IID, cross-cell independence or missing-at-random assumption.",
                    "Raw stratum order follows pinned native writer and rounded logs, not a new full lookup classification.",
                    "Finite native context checks do not prove all historical batch/worker outcomes invariant.",
                    "These are four development-exposed ablations, not the older selected-default eight-tool comparison.",
                    "No biological/other-endpoint intervals, repaired failures/timing, overall superiority or publication readiness."])


def build(repo):
    require(scipy.__version__ == "1.15.3", "Numerical runtime differs from tested SciPy")
    references = {"driver": record(__file__)}
    for name, (relative, sha) in PINS.items():
        ref = record(repo / relative)
        require(ref["sha256"] == sha, "Changed fixed input/source: " + name)
        references[name] = ref
    snapshot = json.loads(Path(references["snapshot"]["path"]).read_text())
    require(snapshot["schema"] == "allocated_native_qfo_scientific_reporting_snapshot_v1"
            and snapshot["publication_ready"] is False and snapshot["new_scoring_or_admission"] is False,
            "Changed incomplete native snapshot scope")
    check(snapshot["source"])
    context = json.loads(Path(references["context"]["path"]).read_text())
    require(context["status"] == "controlled_native_context_invariance_passed"
            and context["assessment"]["all_tested_contexts_invariant"] is True
            and context["assessment"]["pair_evaluations"] == 40,
            "Require retained finite native context check")
    by_index = {r["index"]: r for r in snapshot["rows"]}
    require(len(by_index) == len(snapshot["rows"]), "Duplicate native cell")
    evidence, rows = [snapshot["source"]], []
    for index, cell in CELLS:
        row = by_index[index]
        require(row["cell"] == cell and row["accuracy_admitted"] is True, "Cell not admitted")
        log_path, sha = LOGS[index]
        log_ref = record(repo / log_path)
        require(log_ref["sha256"] == sha, "Changed fixed native log")
        rows.append(read_cell(row, log_ref, references["native_scorer"], evidence))
    result = panel(rows)
    for ref in list(references.values()) + evidence:
        check(ref)
    result.update(references=references, evidence=evidence, scipy_version=scipy.__version__)
    return result


def render(result):
    lines = ["# Four Native Cells: Conditional FAS Sampling Ranges", "",
             "Simultaneous conditional 95% expected-ratio ranges; observed Z is separate.",
             "Both population means are unknown. These are not biological/generalization intervals.", "",
             "| Cell | Observed Z | Conditional Theta Range | k | c | r |",
             "| --- | ---: | --- | ---: | ---: | ---: |"]
    for row in result["methods"]:
        bounds = row["interval"]["expected_native_mean_bounds"] if row["interval"] else None
        value = "[%.9f, %.9f]" % tuple(bounds) if bounds else "Unavailable"
        lines.append("| %s | %.9f | %s | %d | %d | %d |" %
                     (row["cell"], row["observed_native_mean"], value, row["k"], row["c"], row["r"]))
    lines += ["", "| Left - Right | Observed Difference | Conditional Difference Range | Zero Included |",
              "| --- | ---: | --- | --- |"]
    for row in result["contrasts"]:
        bounds = row["conditional_expected_difference_bounds"]
        value = "[%.9f, %.9f]" % tuple(bounds) if bounds else "Unavailable"
        lines.append("| %s - %s | %+.9f | %s | %s |" %
                     (row["left"], row["right"], row["observed_difference"], value, row["zero_included"]))
    return "\n".join(lines + ["", *result["limitations"], ""])


def run(repo, output):
    require(not output.exists() and not output.is_symlink(), "Output already exists")
    output.mkdir(parents=True)
    try:
        result = build(repo)
        with (output / "report.json").open("x") as stream:
            json.dump(result, stream, indent=2, sort_keys=True, allow_nan=False)
            stream.write("\n")
        (output / "summary.md").write_text(render(result))
        return result
    except Exception as exc:
        failure = dict(status="conditional_reporting_failed", type=type(exc).__name__, message=str(exc),
                       traceback=traceback.format_exc(), historical_scores_rerun=False, publication_ready=False)
        with (output / "failure.json").open("x") as stream:
            json.dump(failure, stream, indent=2)
        raise
