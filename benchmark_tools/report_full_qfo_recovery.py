"""Independently admit recovered fresh-search QfO cells without replacing failed attempts."""

import json
from pathlib import Path
import subprocess

from benchmark_tools import recover_full_qfo_search as recovery
from benchmark_tools.prepare_ob_candidate_neighborhood import check
from benchmark_tools.report_hillclimb import endpoint_vector, readback_scores, snapshot


def audit_recovered_cell(repo, row, plan, producer):
    if (row.get("status") != "independent_content_readback_verified"
            or row.get("inputs_sources_runtime_reverified") is not True):
        raise ValueError("Missing recovery postflight validation")
    if Path(plan["repo"]).resolve() != repo.resolve():
        raise ValueError("Recovery belongs to another repository")
    check(plan["parent_plan"])
    parent = json.loads(Path(plan["parent_plan"]["path"]).read_text())
    recovery.verify_inputs(plan, parent, [producer["execution_plan"]])
    receipt = producer["resolved_parent_inputs"]
    if Path(receipt["path"]).resolve() != Path(plan["output"]) / "resolved_parent_inputs.json":
        raise ValueError("Wrong terminal-parent resolution receipt")
    check(receipt)
    resolution = json.loads(Path(receipt["path"]).read_text())
    if resolution.get("status") != "terminal_parent_native_outputs_resolved":
        raise ValueError("Recovery parent inputs were not resolved")
    if Path(resolution["parent"]["path"]).resolve() != Path(plan["parent_results"]):
        raise ValueError("Wrong original failed-attempt report")
    checked = [receipt, resolution["snapshot"], resolution["parent"], *resolution["checked"]]
    for item in checked:
        check(item)
    original = json.loads(Path(resolution["snapshot"]["path"]).read_text())
    accounting = subprocess.check_output(["sacct", "-j", plan["parent_job"], "--parsable2",
        "--format=JobIDRaw,State,ExitCode,Elapsed,NodeList,AllocCPUS"], text=True)
    scheduler = recovery.terminal_scheduler(accounting, plan["parent_job"])
    if scheduler != resolution["scheduler"]:
        raise ValueError("Terminal-parent scheduler receipt differs from current accounting")
    actions = recovery.recovery_actions(original, plan["parent_plan"], plan["parent_job"], scheduler)
    expected_actions = [{"label": label, "action": action} for label, action in zip(recovery.LABELS, actions)]
    if resolution["actions"] != expected_actions or tuple(c["label"] for c in plan["cases"]) != recovery.LABELS:
        raise ValueError("Recovery action or case inventory differs")
    index = recovery.LABELS.index(row["label"])
    action, case = actions[index], plan["cases"][index]
    original_case = recovery.parent_cases(parent)[index]
    if case != recovery.planned_case(repo, Path(plan["output"]), original_case, index) or row["action"] != action:
        raise ValueError("Recovered case differs from the prespecified native/scoring handoff")
    expected_native = case["fresh_case"] if action == "run_missing" else case["original_case"]
    expected_folder = case["original_case"]["output"] if action == "reuse_scored" else case["scoring_output"]
    expected_id = case["original_case"]["scoring_id"] if action == "reuse_scored" else case["scoring_id"]
    if (row["native_case"] != expected_native or row["scoring_output"] != expected_folder
            or row["scoring_id"] != expected_id):
        raise ValueError("Native inference or scoring namespace changed")
    native = recovery.verified_native(expected_native)
    if native["products"] != row["native"]["products"]:
        raise ValueError("Native recovery product identities changed")
    scores, details, files = readback_scores(repo, "QfO", row, native, Path(expected_folder),
                                           "ohmm_hillclimb_r1_" + expected_id)
    checked.extend([*native["checked"], *files])
    recovery.verify_inputs(plan, parent, [producer["execution_plan"], *checked])
    return {"status": "independent_content_readback_verified", "scores": endpoint_vector("QfO", scores),
        "details": details, "coverage": {k: native[k] for k in ("genes", "species", "groups", "pairs")},
        "products": native["products"], "checked": checked, "native_summary": native["summary"],
        "native_source_revision": recovery.REVISION, "recovery_worker_revision": plan["source_revision"],
        "parent_scheduler": scheduler, "original_attempt": resolution["parent"], "action": action,
        "resources": {"native_metrics": native["metrics"], **{k: row[k] for k in
            ("fresh_native_wall_s", "scoring_recovery_wall_s", "total_recovery_wall_s") if k in row}},
        "scope": "Fresh native inference independently re-read; scoring recovery wall time is not full-pipeline runtime. Shared-host timing remains confounded."}


def collect_recovered_cells(repo, root, snapshots):
    path = root / "results.json"
    if not path.exists():
        return [{"method_key": "full_search:" + label, "dataset": "QfO", "scores": {},
                 "status": "recovery_report_missing_execution_state_not_inferred", "attempt_role": "scoring_recovery"}
                for label in recovery.LABELS], []
    producer, capture = snapshot(path, snapshots, "full_search_QfO_recovery_producer.json")
    captures = [capture]
    plan, error = None, None
    try:
        check(producer["execution_plan"])
        plan, plan_capture = snapshot(producer["execution_plan"]["path"], snapshots, "full_search_QfO_recovery_plan.json")
        captures.append(plan_capture)
        if plan["source_revision"] != producer["source_revision"] or Path(plan["output"]).resolve() != root.resolve():
            raise ValueError("Recovery producer/execution plan identity differs")
    except Exception as problem:
        error, plan = str(problem), None
    cells, observed, control_valid = [], set(), False
    for row in producer["arms"]:
        label = row["label"]
        if label not in recovery.LABELS or label in observed:
            raise ValueError("Duplicate or unknown recovery arm")
        observed.add(label)
        cell = {"method_key": "full_search:" + label, "dataset": "QfO", "scores": {},
            "status": row["status"], "producer_status": row["status"], "source_snapshot": capture["snapshot"],
            "claimed_scores": row.get("scores", {}), "producer_job_id": producer.get("job_id"),
            "attempt_role": "scoring_recovery", **{k: row[k] for k in ("error", "error_type") if k in row}}
        if row["status"] == "independent_content_readback_verified":
            try:
                if plan is None:
                    raise ValueError(error)
                if label != "control100" and not control_valid:
                    raise ValueError("Independent recovered matched-control admission failed or is missing")
                cell.update(audit_recovered_cell(repo, row, plan, producer))
                if label == "control100":
                    control_valid = True
            except Exception as problem:
                cell.update(status="readback_failed", error_type=type(problem).__name__, error=str(problem))
        cells.append(cell)
    terminal = producer.get("status", "").startswith(("failed", "finished"))
    cells.extend({"method_key": "full_search:" + label, "dataset": "QfO", "scores": {},
        "status": "not_run_terminal_recovery_gap" if terminal else "not_yet_reported_execution_state_not_inferred",
        "producer_status": producer.get("status"), "attempt_role": "scoring_recovery",
        "source_snapshot": capture["snapshot"], **{k: producer[k] for k in ("error", "error_type") if k in producer}}
        for label in recovery.LABELS if label not in observed)
    return cells, captures
