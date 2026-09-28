"""Bound omitted-score effects on intended FAS samples, not population uncertainty."""

import argparse
import json
import math
from pathlib import Path
import re
import shlex

from benchmark_tools.audit_qfo_fas_samples import read_sample
from benchmark_tools.prepare_ob_candidate_neighborhood import check, record


def one(pattern, text):
    matches = re.findall(pattern, text)
    if len(matches) != 1:
        raise ValueError("Missing or ambiguous FAS log record")
    return matches[0]


def scoring_task(command, raw_name):
    tokens = shlex.split(command, comments=True)
    if not tokens or tokens[0] != "fas_benchmark.py":
        return False
    if tokens.count("--participant") != 1:
        raise ValueError("Ambiguous FAS participant")
    index = tokens.index("--participant") + 1
    if index >= len(tokens) or tokens[index].startswith("--"):
        raise ValueError("Missing FAS participant")
    participant = tokens[index].replace(" ", "-").replace("_", "-")
    return raw_name == f"FAS_{participant}_raw.txt.gz"


def analyze(text, values, eligible):
    pre, missing, unannotated = map(int, one(
        r"(\d+) pairs precomputed, (\d+) missing \(will compute\); (\d+) no feature annotations", text))
    requested_new, requested_pre = map(int, one(
        r"we will compute (\d+) new pairs and sample (\d+) precomputed pairs", text))
    if missing <= 0 or pre + missing != eligible or requested_new != min(missing, 9000):
        raise ValueError("Inconsistent eligible population or request")
    fraction = pre / (pre + missing)
    expected_pre = min(pre, round(requested_new * fraction / (1 - fraction)))
    if requested_pre != expected_pre:
        raise ValueError("Inconsistent requested precomputed stratum")
    strata = {}
    offset = 0
    for name, requested in (("precomputed", requested_pre), ("missing", requested_new)):
        mean, count = one(r"FAS score\[" + name + r"\]: ([0-9.]+) \+- [0-9.]+ \[N=(\d+)\]", text)
        count = int(count)
        if count <= 0 or count > requested or (name == "precomputed" and count != requested):
            raise ValueError("Invalid realized stratum count")
        part = values[offset:offset + count]
        if len(part) != count or any(not math.isfinite(v) or not 0 <= v <= 1 for v in part):
            raise ValueError("Invalid stratum scores")
        observed = math.fsum(part) / count
        if not math.isclose(observed, float(mean), rel_tol=0, abs_tol=5.00001e-7):
            raise ValueError("Ordered raw stratum disagrees with rounded log mean")
        strata[name] = dict(population=pre if name == "precomputed" else missing,
                            requested=requested, saved=count, mean=observed)
        offset += count
    if offset != len(values):
        raise ValueError("Raw sample count differs from logged strata")
    omitted = requested_new - strata["missing"]["saved"]
    intended = requested_pre + requested_new
    score_sum = math.fsum(values)
    return dict(strata=strata, unannotated_pairs_logged=unannotated,
        omitted_requested_new_scores=omitted, intended_sample_size=intended,
        saved_mean=score_sum / len(values),
        intended_sample_mean_bounds=[score_sum / intended, (score_sum + omitted) / intended])


def run(repo):
    source = repo / "benchmark_tools/results/qfo_corrected_fas_sample_audit_20260928.json"
    audit_ref = record(source)
    audit = json.loads(source.read_text())
    if audit["status"] != "corrected_fas_sample_arithmetic_verified" or len(audit["methods"]) != 8:
        raise ValueError("Require corrected eight-method sample audit")
    checked = [audit_ref, record(__file__), record(Path(__file__).with_name("audit_qfo_fas_samples.py")),
               record(repo / "qfo_benchmark/benchmark-webservice/fas_benchmark.py")]
    for ref in audit["checked_records"]:
        check(ref)
    rows = []
    work = repo / "qfo_benchmark/w"
    for method in audit["methods"]:
        raw = method["raw"]
        check(raw)
        candidates = list(work.glob("*/*/*/results/FAS/" + Path(raw["path"]).name))
        matches = [p for p in candidates
                   if scoring_task((p.parents[2] / ".command.sh").read_text(), p.name)
                   and record(p)["sha256"] == raw["sha256"]]
        if len(matches) != 1:
            raise ValueError("Ambiguous matching native task output")
        task = matches[0].parents[2]
        log, command = record(task / ".command.log"), record(task / ".command.sh")
        checked.extend([raw, record(matches[0]), log, command])
        values = list(read_sample(Path(raw["path"])).values())
        result = analyze(Path(log["path"]).read_text(), values, method["reported_eligible_pairs"])
        if len(values) != method["sample_pairs"] or not math.isclose(
                result["saved_mean"], method["mean"], rel_tol=0, abs_tol=1e-12):
            raise ValueError("Sample differs from corrected audit")
        rows.append(dict(method=method["method"], log=log, command=command, raw=raw, **result))
    for ref in [*checked, *audit["checked_records"]]:
        check(ref)
    return dict(status="fas_requested_sample_attrition_bounded", methods=rows,
        checked_records=checked, uncertainty_admitted=False, benchmark_scores_changed=False,
        limitations=[
            "Bounds assume each omitted requested score lies in [0,1]; they condition on intended sampled pairs.",
            "Not confidence intervals, population-mean bounds, sampling-error estimates or family-aware comparisons.",
            "Native logs are newly hashed and linked by a byte-identical historically pinned task output; logs have no independent historical hash pin.",
            "Strata are attributed by native output order and checked against rounded log means, not independently reclassified using precomputed lookup.",
            "Eligible stratum counts are logged values, not independently recounted database predictions.",
            "Missing score identities, causes, RNG state and representativeness remain unresolved."])


def render_table(result):
    lines = ["# FAS Requested-Sample Attrition", "",
        "Generated by `audit_fas_sample_attrition.py`. Bounds condition on the intended",
        "sample and assume omitted scores lie in [0,1]. They are not confidence intervals.", "",
        "| Method | Requested new | Saved new | Omitted new | Saved mean | Intended-sample bound |",
        "| --- | ---: | ---: | ---: | ---: | ---: |"]
    for row in result["methods"]:
        new = row["strata"]["missing"]
        lo, hi = row["intended_sample_mean_bounds"]
        lines.append(f"| {row['method']} | {new['requested']:,} | {new['saved']:,} | "
                     f"{row['omitted_requested_new_scores']:,} | {row['saved_mean']:.6f} | {lo:.6f} - {hi:.6f} |")
    return "\n".join(lines) + "\n"


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", type=Path, default=Path("."))
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--table", type=Path, required=True)
    args = parser.parse_args()
    if any(p.exists() or p.is_symlink() for p in (args.output, args.table)):
        raise FileExistsError("Outputs must be fresh")
    if args.output.resolve() == args.table.resolve():
        raise ValueError("Outputs must be distinct")
    result = run(args.repo.resolve())
    with args.output.open("x") as handle:
        json.dump(result, handle, indent=2, sort_keys=True, allow_nan=False)
        handle.write("\n")
    with args.table.open("x") as handle:
        handle.write(render_table(result))
