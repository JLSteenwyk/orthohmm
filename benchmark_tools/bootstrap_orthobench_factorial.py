"""Prespecified paired RefOG statistics for the eight-cell publication factorial.

Consumes already-validated family sufficient statistics, not native tool
outputs. This module does not certify execution, conversion or reference gates.
"""

from itertools import combinations, product

import numpy as np

from benchmark_tools.bootstrap_orthobench import METRICS, statistics, weighted_records


FACTORS = ("profile_expansion", "candidate_expansion", "reconciliation")
CELLS = tuple(f"p{p}_c{c}_r{r}" for p, c, r in product((0, 1), repeat=3))


def label(values):
    return f"p{values[0]}_c{values[1]}_r{values[2]}"


def conditional_contrasts():
    result = []
    for index, factor in enumerate(FACTORS):
        other = [i for i in range(3) if i != index]
        for fixed in product((0, 1), repeat=2):
            off = [0, 0, 0]
            for i, value in zip(other, fixed):
                off[i] = value
            on = off.copy()
            on[index] = 1
            result.append({"factor": factor, "off": label(off), "on": label(on),
                           "fixed": {FACTORS[i]: off[i] for i in other}})
    return result


def descriptive_interactions(observed):
    results = []
    for first, second in combinations(range(3), 2):
        third = next(i for i in range(3) if i not in (first, second))
        for fixed in (0, 1):
            coefficients = {}
            for a, b in product((0, 1), repeat=2):
                values = [0, 0, 0]
                values[first], values[second], values[third] = a, b, fixed
                coefficients[label(values)] = 1 if a == b else -1
            missing = [cell for cell in coefficients if cell not in observed]
            difference = None if missing else sum(observed[cell] * coefficient
                                                 for cell, coefficient in coefficients.items())
            results.append({"factors": [FACTORS[first], FACTORS[second]], "fixed": {FACTORS[third]: fixed},
                            "coefficients": coefficients, "status": "unavailable" if missing else "complete",
                            "missing_cells": missing,
                            "difference_of_differences_percentage_points": None if missing else dict(zip(METRICS, difference.tolist()))})
    return results


def factorial_bootstrap(cells, replicates=20000, seed=20260918, alpha=0.05):
    if set(cells) != set(CELLS):
        raise ValueError("Require exactly the eight prespecified factorial cells")
    if not isinstance(replicates, int) or isinstance(replicates, bool) or replicates < 100 or not 0 < alpha < 1:
        raise ValueError("Need at least 100 integer replicates and 0 < alpha < 1")
    families, sizes, weights = None, None, {}
    failed = {}
    for cell in CELLS:
        data = cells[cell]
        if data.get("status") == "failed":
            if not data.get("reason") or "refog_records" in data or "score" in data:
                raise ValueError("Failed cells need a reason and must not contain scores")
            failed[cell] = data["reason"]
            continue
        if data.get("status") != "complete":
            raise ValueError("Cannot analyze a running, pending or unknown cell")
        names, current_sizes, counts = weighted_records(data["refog_records"])
        if families is None:
            families, sizes = names, current_sizes
        elif names != families or not np.array_equal(sizes, current_sizes):
            raise ValueError("Completed cells must share identical RefOGs and sizes")
        weights[cell] = counts
    draws, observed = {}, {}
    if weights:
        rng = np.random.Generator(np.random.PCG64(seed))
        multiplicities = rng.multinomial(len(families), np.full(len(families), 1 / len(families)), size=replicates)
        draws = {cell: statistics(multiplicities @ counts) for cell, counts in weights.items()}
        observed = {cell: statistics(counts.sum(axis=0)) for cell, counts in weights.items()}
    comparisons = []
    for contrast in conditional_contrasts():
        on, off = contrast["on"], contrast["off"]
        missing = [cell for cell in (on, off) if cell not in weights]
        result = {**contrast, "status": "unavailable" if missing else "complete", "missing_cells": missing}
        if missing:
            result["metrics"] = None
        else:
            differences = draws[on] - draws[off]
            result["metrics"] = {}
            for index, metric in enumerate(METRICS):
                result["metrics"][metric] = {
                    "difference_percentage_points": float(observed[on][index] - observed[off][index]),
                    "paired_percentile_ci": np.quantile(differences[:, index], [alpha / 2, 1 - alpha / 2]).tolist(),
                    "bonferroni_percentile_ci": np.quantile(differences[:, index], [alpha / 72, 1 - alpha / 72]).tolist()}
            individual = statistics(weights[on])[:, 0] - statistics(weights[off])[:, 0]
            result.update(family_f1_wins=int(np.sum(individual > 1e-10)),
                          family_f1_ties=int(np.sum(np.abs(individual) <= 1e-10)),
                          family_f1_losses=int(np.sum(individual < -1e-10)))
        comparisons.append(result)
    return {"schema_version": 1, "families": families or [], "replicates": replicates, "seed": seed,
            "draws_generated": replicates if weights else 0,
            "alpha": alpha, "rng": "numpy.PCG64 multinomial, shared across all completed cells",
            "numpy_version": np.__version__, "planned_contrasts": 12, "multiplicity_endpoints": 36,
            "point_estimates_percent": {cell: dict(zip(METRICS, values.tolist())) for cell, values in observed.items()},
            "failed_cells": failed, "comparisons": comparisons,
            "descriptive_interactions": descriptive_interactions(observed),
            "limitations": [
                "Development-exposed evidence; intervals are not selection-adjusted independent confirmation.",
                "Resampling units are RefOGs, not gene pairs; shared history and merged predictions can violate exchangeability.",
                "All 36 planned contrast/metric endpoints remain in the Bonferroni adjustment even when cells fail.",
                "Family wins and interaction differences are descriptive, without additional inferential claims.",
                "Weighted F1 is recomputed within every draw, not averaged across family F1 values.",
                "Percentile intervals are approximate; failed cells have no imputed score.",
                "Profile-off retains the HMM-based initial search; these contrasts do not measure the total HMM contribution.",
                "Reconciliation contrasts also change the output from candidate co-membership to final root HOGs.",
                "Execution, native-output conversion, reference coverage and provenance require separate validation."]}


def render_report(result):
    lines = ["# OrthoBench Factorial Uncertainty", "",
             "Development-exposed component analysis, not independent confirmation.", "",
             f"{len(result['families'])} RefOGs; {result['draws_generated']:,} shared paired draws; seed {result['seed']}.", "",
             "| Cell | F1 (%) | Precision (%) | Recall (%) |", "| --- | ---: | ---: | ---: |"]
    for cell in CELLS:
        values = result["point_estimates_percent"].get(cell)
        scores = " | ".join(f"{values[key]:.6f}" for key in METRICS) if values else "NA | NA | NA"
        lines.append(f"| {cell} | {scores} |")
    lines += ["", "Conditional effects are factor on minus factor off, in percentage points.", "",
              "| Factor | On - Off | Metric | Difference | Nominal CI | Adjusted CI |",
              "| --- | --- | --- | ---: | --- | --- |"]
    for comparison in result["comparisons"]:
        for metric in METRICS:
            prefix = f"| {comparison['factor']} | {comparison['on']} - {comparison['off']} | {metric} |"
            if comparison["metrics"] is None:
                lines.append(prefix + " NA | NA | NA |")
            else:
                value = comparison["metrics"][metric]
                intervals = [", ".join(f"{v:.3f}" for v in value[k]) for k in ("paired_percentile_ci", "bonferroni_percentile_ci")]
                lines.append(prefix + f" {value['difference_percentage_points']:.3f} | [{intervals[0]}] | [{intervals[1]}] |")
    lines += ["", "Bonferroni adjustment retains all 36 prespecified contrast/metric endpoints.", "",
              "## Descriptive Family Counts", "",
              "| On - Off | Family F1 Wins | Ties | Losses |", "| --- | ---: | ---: | ---: |"]
    for row in result["comparisons"]:
        counts = " | ".join(str(row[key]) for key in ("family_f1_wins", "family_f1_ties", "family_f1_losses")) if row["metrics"] is not None else "NA | NA | NA"
        lines.append(f"| {row['on']} - {row['off']} | {counts} |")
    lines += ["",
              "## Descriptive Interactions", "",
              "Difference of conditional effects; no interaction confidence intervals or significance claims.", "",
              "| Factors | Fixed Setting | F1 | Precision | Recall |", "| --- | --- | ---: | ---: | ---: |"]
    for row in result["descriptive_interactions"]:
        values = row["difference_of_differences_percentage_points"]
        scores = " | ".join(f"{values[k]:.3f}" for k in METRICS) if values else "NA | NA | NA"
        fixed = ", ".join(f"{key}={value}" for key, value in row["fixed"].items())
        lines.append(f"| {' x '.join(row['factors'])} | {fixed} | {scores} |")
    lines += ["", "## Failed Cells", ""]
    lines += [f"- {cell}: {reason}" for cell, reason in result["failed_cells"].items()] or ["None."]
    lines += ["", "## Limitations", "", *["- " + note for note in result["limitations"]]]
    return "\n".join(lines) + "\n"
