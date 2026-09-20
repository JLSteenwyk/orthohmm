"""Apply prespecified complete-pair engineering budgets, not timing admission."""

import math
from statistics import median

from benchmark_tools.prepare_frontier_overhead_panel import METHODS
from benchmark_tools.verify_root_context_native_provenance import same


def within_budget(value, limit):
    # Ratios at an exact decimal boundary can round a few ULPs above it.
    return value <= limit or math.isclose(value, limit, rel_tol=0., abs_tol=1e-12)


def summarize(plan, rows, panel_issues):
    tasks = plan["runs"]
    if len(tasks) != 18 or len(rows) != 18:
        raise ValueError("Require all 18 prescribed outcomes, including failed and unrun tasks")
    identity = ("index", "block", "pair", "arm", "method")
    for task, row in zip(tasks, rows):
        if any(not same(task[key], row[key]) for key in identity):
            raise ValueError("Outcome order or task identity differs")
        if row["scientific_timings_admitted"] is not False:
            raise ValueError("Engineering comparison cannot admit scientific timing")
    budget = plan["engineering_budget"]
    if not same(budget, dict(every_pair_max=.1, per_method_median_max=.05, required_pairs_per_method=3)):
        raise ValueError("Engineering budgets differ from protocol")
    pairs = []
    for pair in range(9):
        selected = rows[2*pair:2*pair+2]
        arms = {row["arm"]: row for row in selected}
        if (set(arms) != {"lineage", "root_context"} or any(row["pair"] != pair for row in selected)
                or selected[0]["method"] != selected[1]["method"]):
            raise ValueError("Invalid paired design")
        result = dict(pair=pair, block=selected[0]["block"], method=selected[0]["method"],
            order=[row["arm"] for row in selected], indices=[row["index"] for row in selected],
            statuses=[row["status"] for row in selected], signed_ratio=None, pair_budget_passed=None)
        if any(row["status"] != "validated" for row in selected):
            result.update(status="invalid_pair", reason="At least one task failed, was unrun, mismatched or invalid")
        elif any(row.get("output_equivalent") is not True for row in selected):
            result.update(status="invalid_pair", reason="Prior same-method output equivalence not established")
        elif not same(selected[0]["work_identity"], selected[1]["work_identity"]):
            result.update(status="invalid_pair", reason="Paired native output identities differ")
        else:
            walls = [arms[arm]["native_wall_s"] for arm in ("lineage", "root_context")]
            if any(type(v) not in (int, float) or not math.isfinite(v) or v <= 0 for v in walls):
                result.update(status="invalid_pair", reason="Nonpositive or invalid native wall time")
            else:
                ratio = walls[1]/walls[0]-1
                if not math.isfinite(ratio):
                    result.update(status="invalid_pair", reason="Nonfinite elapsed-time ratio")
                else:
                    result.update(status="valid_pair", lineage_native_wall_s=walls[0],
                        root_context_native_wall_s=walls[1], signed_ratio=ratio,
                        pair_budget_passed=within_budget(ratio, budget["every_pair_max"]))
        pairs.append(result)
    methods = []
    for method in METHODS:
        selected = [pair for pair in pairs if pair["method"] == method]
        if len(selected) != 3:
            raise ValueError("Require three pairs for every method")
        valid = [pair for pair in selected if pair["status"] == "valid_pair"]
        value = median(pair["signed_ratio"] for pair in valid) if len(valid) == 3 else None
        methods.append(dict(method=method, valid_pairs=len(valid), pair_indices=[p["pair"] for p in selected],
            median_signed_ratio=value, median_budget_passed=None if value is None else within_budget(value, budget["per_method_median_max"])))
    complete = not panel_issues and all(pair["status"] == "valid_pair" for pair in pairs)
    passed = (all(pair["pair_budget_passed"] for pair in pairs)
              and all(method["median_budget_passed"] for method in methods)) if complete else None
    return dict(status="root_context_incremental_elapsed_comparison", pairs=pairs, methods=methods,
        budget_comparison_absolute_tolerance=1e-12,
        panel_issues=panel_issues, all_pairs_and_panel_valid=complete, engineering_budget_passed=passed,
        scientific_timings_admitted=False, publication_ready=False,
        limitations=["Consumes independently audited outcomes; does not replace source, raw, scheduler or output validation.",
            "Signed elapsed-time changes retain negative values; no correction of scientific timings.",
            "Complete three-pair medians only; failed/unrun/mismatched observations are not discarded.",
            "Engineering budgets are not confidence bounds, causal overhead estimates or scientific timing eligibility."])
