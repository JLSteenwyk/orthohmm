"""Decompose FAS sample-mixture weights without replacing native endpoints."""

import argparse
from fractions import Fraction
import json
import math
from pathlib import Path

from benchmark_tools.prepare_ob_candidate_neighborhood import check, record

INPUTS = (
    ("benchmark_tools/results/qfo_fas_sample_attrition_20260928.json",
     "1f6ddfd5aadf9f5be8147164b23904d2e5db1ffee69a6fc18b654a54f3f25d36"),
    ("benchmark_tools/results/qfo_fas_population_completed_22383/report.json",
     "05f1ab3f7f0ffed8ddcffe2a6bfd4ce26ef46b974dbb946cb60883857803eeaf"),
)
METHODS = (
    "orthohmm_high_sensitivity", "orthohmm_phylogeny_satellite_v2",
    "orthofinder_3_1_5_full", "orthofinder_3_1_5_sequence_only",
    "sonicparanoid_2_0_9", "proteinortho_6_3_6", "fastoma_0_3_5", "orthomcl_1_4",
)


def decompose(sample, population):
    if sample["method"] != population["method"]:
        raise ValueError("Method identities differ")
    pre, new = (sample["strata"][key] for key in ("precomputed", "missing"))
    for key, stratum in (("precomputed", pre), ("missing", new)):
        counts = [stratum[field] for field in ("population", "requested", "saved")]
        if any(type(v) is not int for v in counts) or not 0 < counts[2] <= counts[1] <= counts[0]:
            raise ValueError("Require two observed, nonempty sample strata")
        if stratum["population"] != population[key]:
            raise ValueError("Sample and population recount differ")
        value = stratum["mean"]
        if type(value) not in (int, float) or not math.isfinite(value) or not 0 <= value <= 1:
            raise ValueError("Invalid saved stratum mean")
    p, m, k, c, r = (pre["population"], new["population"], pre["requested"],
                      new["requested"], new["saved"])
    if (pre["saved"] != k or c != min(m, 9000) or population["eligible_pairs"] != p + m
            or population["native_counts_match"] is not True
            or population["saved_lookup_strata_and_values_match"] is not True):
        raise ValueError("Retained native sample design or population checks differ")
    fraction = p / (p + m)
    if k != min(p, round(c * fraction / (1 - fraction))):
        raise ValueError("Requested precomputed count differs from native float arithmetic")
    if (type(sample["omitted_requested_new_scores"]) is not int
            or sample["omitted_requested_new_scores"] != c - r
            or sample["intended_sample_size"] != k + c):
        raise ValueError("Sample attrition totals differ")

    # Exact rational weights avoid cancellation in the small rounding component.
    eligible = Fraction(p, p + m)
    intended = Fraction(k, k + c)
    saved = Fraction(k, k + r)
    rounding, attrition = intended - eligible, saved - intended
    gap = pre["mean"] - new["mean"]
    reconstructed = float(saved) * pre["mean"] + float(1 - saved) * new["mean"]
    if not math.isclose(reconstructed, sample["saved_mean"], rel_tol=0, abs_tol=1e-12):
        raise ValueError("Saved endpoint does not match retained stratum arithmetic")
    reweighted = float(eligible) * pre["mean"] + float(1 - eligible) * new["mean"]
    shift = float(saved - eligible) * gap
    if not math.isclose(reconstructed - reweighted, shift, rel_tol=0, abs_tol=1e-12):
        raise ValueError("Weight decomposition does not reproduce numerical difference")
    return {
        "method": sample["method"], "eligible_pairs": p + m,
        "requested_precomputed": k, "requested_new": c, "saved_new": r,
        "eligible_precomputed_weight": float(eligible),
        "intended_precomputed_weight": float(intended),
        "saved_precomputed_weight": float(saved),
        "rounding_weight_shift": float(rounding),
        "attrition_weight_shift": float(attrition),
        "rounding_half_unit_condition_met": abs(Fraction(k) - Fraction(c * p, m)) <= Fraction(1, 2),
        "precomputed_minus_new_saved_mean": gap,
        "native_saved_mean": sample["saved_mean"],
        "reweighted_saved_strata_diagnostic": reweighted,
        "native_minus_diagnostic": shift,
        "rounding_component": float(rounding) * gap,
        "attrition_component": float(attrition) * gap,
        "weight_only_absolute_difference_bound": float(abs(saved - eligible)),
    }


def panel(attrition, population):
    if (attrition["status"] != "fas_requested_sample_attrition_bounded"
            or population["status"] != "retained_fas_eligible_populations_completed_with_reuse"
            or any(report["benchmark_scores_changed"] is not False
                   or report["uncertainty_admitted"] is not False for report in (attrition, population))
            or population["publication_ready"] is not False):
        raise ValueError("Require unchanged, unadmitted retained diagnostic reports")
    if any([row["method"] for row in report["methods"]] != list(METHODS)
           for report in (attrition, population)):
        raise ValueError("Require all eight fixed methods in their original order")
    return [decompose(a, p) for a, p in zip(attrition["methods"], population["methods"])]


def run(repo):
    reports, records = [], [record(__file__), record(Path(__file__).with_name("prepare_ob_candidate_neighborhood.py"))]
    for relative, sha in INPUTS:
        ref = record(repo / relative)
        if ref["sha256"] != sha:
            raise ValueError("Retained audit report identity changed")
        reports.append(json.loads(Path(ref["path"]).read_text()))
        records.append(ref)
    methods = panel(*reports)
    for ref in records:
        check(ref)
    return dict(status="retained_fas_stratum_weight_decomposition_verified", methods=methods,
        checked_records=records, benchmark_scores_changed=False, uncertainty_admitted=False,
        publication_ready=False, limitations=[
            "Aggregate readback only: no new lookup scan, database join, sample reconstruction or FAS rescore.",
            "Existing sample/population audit limitations and historical provenance gaps remain unchanged.",
            "Reweighted saved strata are a numerical diagnostic, not a corrected benchmark or population estimate.",
            "Weight effects hold the observed stratum means fixed; omitted scores can also change the new-stratum mean.",
            "Uniform sampling and ignorable score omissions are not established; RNG state and omitted pair IDs are absent.",
            "No family-aware paired confidence interval, significance claim or accuracy ranking is admitted."])


def render(result):
    lines = ["# Retained FAS Stratum Weights", "",
        "All eight native FAS scores remain unchanged. This is aggregate arithmetic,",
        "not a corrected benchmark, a population estimate or a confidence interval.", "",
        "Weights are precomputed-stratum percentages; differences are in FAS score units.", "",
        "| Method | Eligible % | Intended % | Saved % | Rounding component | Attrition component | Total weight effect |",
        "| --- | ---: | ---: | ---: | ---: | ---: | ---: |"]
    for row in result["methods"]:
        weights = [100 * row[key] for key in ("eligible_precomputed_weight",
                   "intended_precomputed_weight", "saved_precomputed_weight")]
        shifts = [row[key] for key in ("rounding_component", "attrition_component", "native_minus_diagnostic")]
        lines.append("| " + row["method"] + " | " + " | ".join(
            [*(f"{value:.6f}" for value in weights), *(f"{value:+.6e}" for value in shifts)]) + " |")
    return "\n".join(lines) + "\n"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", type=Path, default=Path("."))
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--table", type=Path, required=True)
    args = parser.parse_args()
    if args.output.resolve() == args.table.resolve():
        raise ValueError("Outputs must be distinct")
    if any(path.exists() or path.is_symlink() for path in (args.output, args.table)):
        raise FileExistsError("Outputs must be fresh")
    result = run(args.repo.resolve())
    with args.output.open("x") as handle:
        json.dump(result, handle, indent=2, sort_keys=True, allow_nan=False)
        handle.write("\n")
    with args.table.open("x") as handle:
        handle.write(render(result))


if __name__ == "__main__":
    main()
