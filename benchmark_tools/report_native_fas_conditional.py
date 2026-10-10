"""Report retained native FAS design uncertainty, not biological error bars."""

import argparse
from fractions import Fraction
from itertools import combinations
import json
import math
from pathlib import Path
import traceback

import scipy

from benchmark_tools import audit_fas_stratum_weights as audit
from benchmark_tools import native_fas_sampling_interval as kernel
from benchmark_tools.probe_native_fas_context import fingerprint, require


PROTOCOL_SHA = "c22a7b602614b37327a3d64d50b9fcb36053ca9d7f643ebe92390c057f99c347"
PINS = {
    "attrition": ("benchmark_tools/results/qfo_fas_sample_attrition_20260928.json", "1f6ddfd5aadf9f5be8147164b23904d2e5db1ffee69a6fc18b654a54f3f25d36"),
    "population": ("benchmark_tools/results/qfo_fas_population_completed_22383/report.json", "05f1ab3f7f0ffed8ddcffe2a6bfd4ce26ef46b974dbb946cb60883857803eeaf"),
    "validation": ("benchmark_tools/results/native_fas_sampling_validation_20261010_v1.json", "236ffb82887d816904d08df0e5b350a358badd8e9f5f6fc088004cd0d021ac57"),
    "context": ("benchmark_tools/results/native_fas_context_probe_20261010_v2.json", "845b4baea528f7633830e5a26ec58267da95165d8c4f79f2a0112c6063826f33"),
    "kernel": ("benchmark_tools/native_fas_sampling_interval.py", "803ff0f455213e6a55b83e30c586ab1c57566830f8f877c916b78e84bf3658f9"),
    "audit": ("benchmark_tools/audit_fas_stratum_weights.py", "8712b7e36fc756b61c161847493f7c65b7b376d5ef675c6c803f52ab21d5fbc3"),
    "provenance_helper": ("benchmark_tools/probe_native_fas_context.py", "d722944a7cadeee3e25433fc472cdf724ba58309aaac29952338e5c7575bae06"),
}
ALPHA = .05
COMPONENT_ERROR = ALPHA / 16


def prerequisites(validation, context):
    require(validation["status"] == "prospective_validation_complete" and len(validation["cells"]) == 15,
            "Require completed prospective validation")
    require(validation["kernel"]["sha256"] == PINS["kernel"][1] and
            Fraction(validation["component_error_fraction"]) == Fraction(1, 320), "Validation binding/allocation differs")
    for row in validation["cells"]:
        require(Fraction(row["count_coverage_fraction"]) >= 1-Fraction(1,320) and
                Fraction(row["mean_coverage_fraction"]) >= 1-Fraction(1,320) and
                Fraction(row["target_coverage_fraction"]) >= 1-Fraction(1,160) and
                row["numerical_max_error"] <= 1e-12, "Validation coverage/numerics failed")
    joint = validation["joint"]
    require(joint["method_count"] == 8 and joint["differences"] == 28 and
            Fraction(joint["arbitrary_dependence_union_failure_bound_fraction"]) <= Fraction(1,20) and
            Fraction(joint["coupled_difference_coverage_fraction"]) >= Fraction(19,20), "Joint validation failed")
    require(validation["context_control"]["status"] == "fixed_return_assumption_rejected" and
            validation["context_control"]["coverage_claimed"] is False, "Invalid control admission")
    require(context["status"] == "controlled_native_context_invariance_passed" and
            context["assessment"]["complete"] is True and context["assessment"]["pair_evaluations"] == 40 and
            context["assessment"]["all_tested_contexts_invariant"] is True and
            context["assessment"]["variations"] == [], "Require completed finite context check")
    for report in (validation, context):
        require(all(report[key] is False for key in ("historical_scores_rerun", "native_sampling_law_admitted",
                                                     "native_intervals_admitted", "publication_ready")),
                "Prerequisite scope changed")


def panel(attrition, population):
    audit.panel(attrition, population)
    require(type(population["reuse"]["historical_parser_hash_identity_established"]) is bool,
            "Missing historical parser scope")
    methods = []
    for sample, full in zip(attrition["methods"], population["methods"]):
        new, old = sample["strata"]["missing"], sample["strata"]["precomputed"]
        muP = full["precomputed_mean"]
        require(type(muP) in (int,float) and math.isfinite(muP) and 0 <= muP <= 1 and
                math.isclose(full["precomputed_score_sum"]/full["precomputed"], muP, rel_tol=0, abs_tol=1e-12),
                "Invalid full precomputed population mean")
        require(type(full["database_historically_hash_bound"]) is bool, "Missing historical database binding scope")
        row = dict(method=sample["method"], observed_native_mean=sample["saved_mean"],
                   population=dict(P=full["precomputed"],M=full["missing"]),
                   sample=dict(k=old["requested"],c=new["requested"],r=new["saved"],
                               omitted_numeric_returns=sample["omitted_requested_new_scores"]),
                   precomputed_population_mean=muP, returned_sample_mean=new["mean"],
                   database_historically_hash_bound=full["database_historically_hash_bound"])
        methods.append(row)
    for row in methods:
        try:
            interval = kernel.method_interval(row["population"]["M"],row["sample"]["c"],row["sample"]["k"],
                                              row["precomputed_population_mean"],row["sample"]["r"],
                                              row["returned_sample_mean"],COMPONENT_ERROR)
            lower, upper = interval["expected_native_mean_bounds"]
            require(math.isfinite(lower) and math.isfinite(upper) and 0 <= lower <= upper <= 1,
                    "Invalid conditional expected-ratio range")
            row.update(status="conditional_design_range_computed",interval=interval,error=None)
        except Exception as exc:
            row.update(status="conditional_numerical_failure",interval=None,
                       error=dict(type=type(exc).__name__,message=str(exc),traceback=traceback.format_exc()))
        print(row["method"],row["status"],flush=True)
    contrasts = []
    for left,right in combinations(methods,2):
        available = left["interval"] is not None and right["interval"] is not None
        bounds = kernel.difference_interval(left["interval"],right["interval"]) if available else None
        contrasts.append(dict(left=left["method"],right=right["method"],
                              observed_difference=left["observed_native_mean"]-right["observed_native_mean"],
                              conditional_expected_difference_bounds=bounds,
                              status="conditional_design_range_computed" if available else "unavailable_dependent_method",
                              zero_included=bounds[0] <= 0 <= bounds[1] if available else None))
    complete = all(row["interval"] is not None for row in methods)
    return dict(schema="retained_conditional_native_fas_v1",
                status="conditional_expected_native_ranges_complete" if complete else "conditional_expected_native_ranges_partial",
                methods=methods,contrasts=contrasts,alpha=ALPHA,component_error=COMPONENT_ERROR,
                allocation="Two components per each of eight methods; all 28 differences use the same joint rectangle.",
                all_eight_ranges_computed=complete,conditional_design_ranges_computed=any(row["interval"] is not None for row in methods),
                observed_native_scores_changed=False,historical_scores_rerun=False,
                historical_parser_hash_identity_established=population["reuse"]["historical_parser_hash_identity_established"],
                unconditional_historical_interval_admission=False,biological_generalization_intervals=False,
                other_endpoint_uncertainty_admitted=False,publication_ready=False,
                target="expected_native_post_attrition_ratio_under_fixed_design",
                limitations=[
                    "Conditional uniform fixed-population sampling and fixed pair outcomes are assumptions, not historical certification.",
                    "Observed fixed means are separate descriptive scores, not point values or centers of the expected-ratio ranges.",
                    "Known precomputed population means are conditioned on the retained audit; historical parser/database hash gaps persist.",
                    "Finite context checks do not establish every historical annotation, file race, batch failure or worker choice.",
                    "No biological pair/family or cross-method independence, missing-at-random claim, new score, sample or rescore.",
                    "No interval for another endpoint, overall tool superiority, independent biological confirmation or publication readiness."])


def run(repo, protocol):
    references = {"driver":fingerprint(__file__),"protocol":fingerprint(protocol)}
    require(references["protocol"]["sha256"] == PROTOCOL_SHA,"Changed prospective reporting protocol")
    inputs = {}
    for name,(relative,sha) in PINS.items():
        ref = fingerprint(repo/relative)
        require(ref["sha256"] == sha,"Changed frozen input/source: " + name)
        references[name] = ref
        if name in ("attrition","population","validation","context"):
            inputs[name] = json.loads(Path(ref["path"]).read_text())
    prerequisites(inputs["validation"],inputs["context"])
    require(scipy.__version__ == inputs["validation"]["scipy_version"],"Validated numerical runtime differs")
    result = panel(inputs["attrition"],inputs["population"])
    require(all(fingerprint(ref["path"]) == ref for ref in references.values()),"Inputs changed during reporting")
    result.update(references=references,scipy_version=scipy.__version__)
    return result


def bounds_text(bounds):
    return "Unavailable" if bounds is None else "[%.9f, %.9f]" % tuple(bounds)


def render(result):
    lines = ["# Conditional Native FAS Design Ranges", "",
             "Observed Z is preserved separately from the simultaneous conditional 95% range",
             "for expected repeated native ratio theta. Unknown G and mu_G preclude an exact",
             "theta point value. These are not biological error bars on Z or unconditional",
             "historical certificates. Table decimals are rounded; JSON retains full values.", "",
             "| Method | Observed Z | Conditional Theta Range | G Bounds |",
             "| --- | ---: | --- | --- |"]
    for row in result["methods"]:
        interval = row["interval"]
        g = str(interval["success_count_bounds"]) if interval else "Unavailable"
        lines.append("| %s | %.9f | %s | %s |" % (row["method"],row["observed_native_mean"],
                     bounds_text(interval["expected_native_mean_bounds"] if interval else None),g))
    lines += ["", "All 28 contrasts use left-minus-right in the frozen method order.", "",
              "| Left | Right | Observed Z Difference | Conditional Expected Difference |",
              "| --- | --- | ---: | --- |"]
    for row in result["contrasts"]:
        lines.append("| %s | %s | %.9f | %s |" % (row["left"],row["right"],row["observed_difference"],
                     bounds_text(row["conditional_expected_difference_bounds"])))
    lines += ["", "## Scope", ""] + ["- " + item for item in result["limitations"]]
    if not result["all_eight_ranges_computed"]:
        lines += ["", "This is a partial numerical report; unavailable methods and dependent contrasts are not filled."]
    return "\n".join(lines)+"\n"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo",type=Path,default=Path.cwd())
    parser.add_argument("--protocol",type=Path,required=True)
    parser.add_argument("--output",type=Path,required=True)
    args = parser.parse_args()
    if args.output.exists() or args.output.is_symlink():
        raise FileExistsError(args.output)
    args.output.mkdir()
    try:
        result = run(args.repo,args.protocol)
    except Exception as exc:
        failure = dict(status="conditional_report_input_or_execution_failure",type=type(exc).__name__,
                       message=str(exc),traceback=traceback.format_exc(),publication_ready=False)
        (args.output/"failure.json").write_text(json.dumps(failure,indent=2)+"\n")
        print(json.dumps(failure))
        return 1
    with (args.output/"report.json").open("x") as stream:
        json.dump(result,stream,indent=2,sort_keys=True,allow_nan=False)
        stream.write("\n")
    (args.output/"summary.md").write_text(render(result))
    print(json.dumps(dict(status=result["status"],methods=len(result["methods"]),contrasts=len(result["contrasts"]))))
    return 0 if result["all_eight_ranges_computed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
