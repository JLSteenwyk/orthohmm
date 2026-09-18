"""Pure paired arithmetic; callers must independently audit native evidence."""

import csv
import io
import math
import re
import statistics

METHODS = ("orthohmm_high_sensitivity", "orthohmm_satellite_v2", "orthofinder_full")
STATUSES = {"validated", "failed", "missing_evidence", "invalid_evidence"}


def terminal_scheduler_rows(text, array_id):
    """Check complete terminal task inventory before inspecting native outcomes."""
    if type(array_id) is not int or array_id <= 0:
        raise ValueError("Invalid array identity")
    pattern = re.compile(rf"{array_id}_(\d+)")
    rows = {}
    terminal = {"COMPLETED", "FAILED", "TIMEOUT", "OUT_OF_MEMORY", "CANCELLED",
                "NODE_FAIL", "PREEMPTED", "BOOT_FAIL", "DEADLINE", "REVOKED"}
    for fields in csv.reader(io.StringIO(text), delimiter="|"):
        if not fields:
            continue
        match = pattern.fullmatch(fields[0])
        if not match:
            continue
        index = int(match.group(1))
        if len(fields) != 7 or index in rows or not 0 <= index < 18:
            raise ValueError("Malformed or duplicate scheduler task")
        if fields[1].split(" ", 1)[0] not in terminal:
            raise ValueError("Panel has nonterminal tasks; do not inspect outcomes")
        if not re.fullmatch(r"\d+:\d+", fields[2]):
            raise ValueError("Malformed scheduler exit status")
        rows[index] = dict(job_id=fields[0], state=fields[1], exit_code=fields[2],
                           elapsed=fields[3], allocated_cpus=fields[4], requested_memory=fields[5], node=fields[6])
    if set(rows) != set(range(18)):
        raise ValueError("Incomplete terminal scheduler inventory")
    return rows


def finite_positive(value):
    return type(value) in (int, float) and math.isfinite(value) and value > 0


def summarize(plan, observations):
    tasks = plan["runs"]
    if ([task["index"] for task in tasks] != list(range(18))
            or plan["engineering_budget"] != {"per_method_median_max": .05, "every_pair_max": .10}):
        raise ValueError("Require complete frozen panel and unchanged budgets")
    by_index = {}
    for row in observations:
        index = row["index"]
        if type(index) is not int or not 0 <= index < 18 or index in by_index:
            raise ValueError("Duplicate or invalid task index")
        if row["status"] not in STATUSES:
            raise ValueError("Unknown audit status")
        task = tasks[index]
        if any(row[key] != task[key] for key in ("method", "mode", "pair")):
            raise ValueError("Task assignment differs from plan")
        if row["status"] == "validated":
            if not finite_positive(row["native_wall_s"]):
                raise ValueError("Invalid successful native duration")
            if not isinstance(row["work_identity"], dict) or not row["work_identity"]:
                raise ValueError("Require canonical output identity from independent audit")
            if type(row["whole_command_screen_passed"]) is not bool:
                raise ValueError("Require explicit whole-command screening result")
            flags = row["flagged_intervals"]
            if row["mode"] == "boundary":
                if flags is not None:
                    raise ValueError("Boundary-only interval screening is unavailable")
            elif (not isinstance(flags, list) or any(type(i) is not int or i < 0 for i in flags)
                  or len(set(flags)) != len(flags)):
                raise ValueError("Invalid periodic interval flags")
        by_index[index] = row
    pairs = []
    for method in METHODS:
        for pair in range(3):
            assignment = [task for task in tasks if task["method"] == method and task["pair"] == pair]
            if len(assignment) != 2 or {task["mode"] for task in assignment} != {"boundary", "periodic"}:
                raise ValueError("Incomplete paired assignment")
            arms = {task["mode"]: by_index.get(task["index"]) for task in assignment}
            reasons = []
            for mode, row in arms.items():
                if row is None:
                    reasons.append(mode + ":missing_audit")
                elif row["status"] != "validated":
                    reasons.append(mode + ":" + row["status"])
            ratio, duration_gate, work_equal = None, None, None
            if not reasons:
                boundary, periodic = arms["boundary"], arms["periodic"]
                work_equal = boundary["work_identity"] == periodic["work_identity"]
                duration_gate = min(boundary["native_wall_s"], periodic["native_wall_s"]) >= 60
                # Retain the descriptive ratio even if equivalent-work/duration gates fail.
                ratio = (periodic["native_wall_s"] - boundary["native_wall_s"]) / boundary["native_wall_s"]
                if not math.isfinite(ratio):
                    raise ValueError("Nonfinite paired ratio")
            pairs.append(dict(method=method, pair=pair, arms=arms, unavailable_reasons=reasons,
                              wall_ratio_minus_one=ratio, equivalent_work=work_equal,
                              duration_gate_met=duration_gate,
                              numerical_pair_budget_met=None if ratio is None else ratio <= .10))
    summaries = []
    for method in METHODS:
        selected = [pair for pair in pairs if pair["method"] == method]
        values = [pair["wall_ratio_minus_one"] for pair in selected if pair["wall_ratio_minus_one"] is not None]
        complete = len(values) == 3
        summaries.append(dict(method=method, available_pairs=len(values), expected_pairs=3,
            median=statistics.median(values) if complete else None,
            minimum=min(values) if complete else None, maximum=max(values) if complete else None,
            numerical_budget_met=(statistics.median(values) <= .05 and max(values) <= .10) if complete else None))
    complete = all(row["available_pairs"] == 3 for row in summaries)
    return dict(status="paired_overhead_arithmetic_only", pairs=pairs, methods=summaries,
        complete_numeric_panel=complete,
        numerical_budget_met=all(row["numerical_budget_met"] for row in summaries) if complete else None,
        equivalent_work_and_duration_met=(all(pair["equivalent_work"] and pair["duration_gate_met"] for pair in pairs)
                                         if complete else None),
        raw_evidence_audited=False, environmental_validity_established=False,
        scientific_timings_admitted=False, publication_ready=False,
        limitations=["Pure arithmetic over caller-supplied audited outcomes; not a native-output or provenance audit.",
            "Missing/failed pairs are retained; incomplete panels have no complete-panel budget result.",
            "Original flags remain in arm records; numerical budget results do not override them.",
            "Boundary arms have no interval coverage; environmental validity is not established here.",
            "Negative ratios reflect variability, not proof of negative overhead; no timing correction."])
