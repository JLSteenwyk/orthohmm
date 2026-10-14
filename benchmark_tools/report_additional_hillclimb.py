"""Incrementally admit the separate QfO continuation and shared-hit trials."""

import argparse
from copy import deepcopy
import csv
import hashlib
import json
from pathlib import Path
import subprocess
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.metric_frontier import ENDPOINTS
from benchmark_tools.prepare_ob_candidate_neighborhood import check, record
from benchmark_tools.recover_full_qfo_search import verified_native
from benchmark_tools.report_current_hillclimb import canonical_cell, upsert_cell
from benchmark_tools.report_hillclimb import endpoint_vector, readback_scores, screen, snapshot
from benchmark_tools.report_hillclimb import replay_native
from benchmark_tools.run_full_search_hillclimb import runtime, validate_native, verify_input_inventory
from benchmark_tools.run_membership_hillclimb import DATASETS, write_json
from benchmark_tools.run_shared_hit_pipeline import selected_control


def check_plan(repo, root, report):
    check(report["execution_plan"])
    plan = json.loads(Path(report["execution_plan"]["path"]).read_text())
    if Path(plan["repo"]).resolve() != repo or Path(plan["output"]).resolve() != root:
        raise ValueError("Producer plan belongs to another repository or output")
    pins = [report["execution_plan"], *plan["inputs"], *plan["sources"], *plan["runtime"]["files"]]
    for pin in pins:
        check(pin)
    if runtime() != plan["runtime"]:
        raise ValueError("Scientific runtime differs from the pinned producer")
    return plan, pins


def fresh_cell(repo, dataset, row, native, folder, participant, revision):
    scores, details, pins = readback_scores(repo, dataset, row, native, folder, participant)
    checked = [*native["checked"], *pins]
    for pin in checked:
        check(pin)
    return {"dataset": dataset, "status": "independent_content_readback_verified",
            "scores": endpoint_vector(dataset, scores), "details": details,
            "products": native["products"], "checked": checked,
            "coverage": {k: native[k] for k in ("genes", "species", "groups", "pairs")},
            "source_revision": revision, "resources": row["execution"],
            "execution_kind": "fresh_full_pipeline_with_separate_scoring",
            "scope": "Fresh variant inference. Shared-host cost is not comparable to cached control replay cost."}


def continuation_cell(repo, plan, report):
    case = plan["case"]["fresh_case"]
    if (case["dataset"] != "QfO" or case["label"] != "cap400"
            or case["search_max_candidates_per_query"] != 400 or report["execution"]["returncode"] != 0):
        raise ValueError("Require the successful original cap400 continuation")
    verify_input_inventory(case)
    native = verified_native(case)
    cell = fresh_cell(repo, "QfO", report, native, Path(plan["case"]["scoring_output"]),
                      "ohmm_hillclimb_r1_" + plan["case"]["scoring_id"], plan["scientific_revision"])
    for pin in report["control"]["checked"]:
        check(pin)
    cell["execution_kind"] = "fresh_original_arrival_cap400_continuation"
    return cell


def shared_native(case, plan):
    if case["search_min_total_hits"] != 1 or case["search_max_candidates_per_query"] != 100:
        raise ValueError("Wrong shared-hit intervention")
    verify_input_inventory(case)
    native = validate_native(case, 16)
    metadata = native["metrics"]["metadata"]
    if metadata.get("search_result_order") != "canonical" or metadata.get("search_min_total_hits") != 1:
        raise ValueError("Actual native shared-hit mechanism differs from the plan")
    harness = native["metrics"]["harness"]
    source = Path(plan["source_root"])
    if harness.get("exit_code") != 0:
        raise ValueError("Shared-hit native harness exit differs")
    # Git in a plain frozen copy discovers the mutable ancestor checkout.
    if (source / ".git").exists() and harness.get("git_commit") != plan["source_revision"]:
        raise ValueError("Shared-hit frozen worktree revision differs")
    expected_sources = [record(p) for p in sorted((source / "orthohmm").rglob("*.py"))]
    expected_sources.append(record(source / "benchmark_tools/benchmark_production.py"))
    actual_sources = {r["path"]: (r["bytes"], r["sha256"]) for r in harness["source_manifest"]}
    expected = {str(Path(r["path"]).relative_to(source)): (r["bytes"], r["sha256"]) for r in expected_sources}
    expected_inputs = {Path(r["path"]).name: (r["bytes"], r["sha256"]) for r in case["inputs"]}
    if (len(actual_sources) != len(harness["source_manifest"]) or actual_sources != expected
            or len(harness["input_manifest"]) != len(case["inputs"])
            or {r["path"]: (r["bytes"], r["sha256"]) for r in harness["input_manifest"]} != expected_inputs):
        raise ValueError("Shared-hit native source/input manifests differ")
    return native


def shared_cell(repo, root, plan, report, row):
    dataset = row["dataset"]
    cases = [case for case in plan["cases"] if case["dataset"] == dataset and case["label"] == "shared_hit1"]
    if len(cases) != 1 or report["source_revision"] != plan["source_revision"]:
        raise ValueError("Ambiguous shared-hit case or changed source revision")
    case = cases[0]
    if row["execution"]["returncode"] != 0 or row["control"]["initial_min_total_hits"] != 4:
        raise ValueError("Wrong shared-hit intervention or control threshold")
    native = shared_native(case, plan)
    harness = native["metrics"]["harness"]
    control_path = report["control_report"]
    check(control_path)
    control_report = json.loads(Path(control_path["path"]).read_text())
    control_root = Path(plan["control_root"])
    control_pin = control_report["execution_plan"]
    if (Path(control_path["path"]).resolve() != control_root / "results.json"
            or Path(control_pin["path"]).resolve() != control_root / "execution_plan.json"
            or control_pin not in plan["inputs"]):
        raise ValueError("Canonical control is not bound to the frozen trial plan")
    check(control_pin)
    control_plan = json.loads(Path(control_pin["path"]).read_text())
    if not (control_plan["source_revision"] == control_report["source_revision"] == row["control"]["source_revision"]):
        raise ValueError("Canonical control source revision differs")
    control = selected_control(control_report, dataset)
    control_cell = canonical_cell(repo, control_root, control, control_report)
    if (endpoint_vector(dataset, row["control"]["scores"]) != control_cell["scores"]
            or (native["genes"], native["species"]) !=
            (control_cell["coverage"]["genes"], control_cell["coverage"]["species"])):
        raise ValueError("Independent canonical control scores or universe differ")
    folder = Path(case["output"])
    native.update(prediction=str(folder / "output/orthohmm_working_res/orthohmm_edges_clustered.txt"),
                  phylogeny=str(folder / "output/orthohmm_phylogeny"),
                  checked=[*native["products"].values(), record(folder / "metrics.json"),
                           record(folder / "output/orthohmm_phylogeny/provenance_manifest.json"),
                           record(folder / "output/orthohmm_phylogeny/reconciliation_summary.json"),
                           *[record(p) for p in sorted((folder / "output/orthohmm_working_res/high_sensitivity_checkpoint").iterdir())]])
    nonce = hashlib.sha256(str(root).encode()).hexdigest()[:16]
    cell = fresh_cell(repo, dataset, row, native, folder, "ohmm_hillclimb_r1_sharedhit_" + nonce,
                      plan["source_revision"])
    cell["matched_control"] = {"report": control_path, "cell": control_cell}
    cell["harness_git_metadata"] = {
        "git_commit": harness.get("git_commit"), "git_dirty": harness.get("git_dirty"),
        "scope": "Contextual Git metadata; frozen source identity is verified by source manifests and execution-plan pins."}
    cell["execution_kind"] = "fresh_canonical_initial_shared_hit1_cap100"
    return cell


def canonical_recovery_cell(repo, plan, report, root):
    from benchmark_tools.recover_canonical_qfo_scoring import REVISION, failed_row, verify_parent
    if (report.get("status") != "canonical_qfo_scoring_recovered_independently_verified"
            or report.get("original_attempt_preserved") is not True
            or report.get("original_status") != "failed" or plan["scientific_revision"] != REVISION):
        raise ValueError("Require independently verified scoring-only canonical recovery")
    check(plan["parent_snapshot"])
    check(plan["protocol"])
    parent = json.loads(Path(plan["parent_snapshot"]["path"]).read_text())
    if failed_row(parent) != plan["original_failed_row"]:
        raise ValueError("Archived original failed row differs")
    native, pins = verify_parent(repo, Path(plan["parent_root"]), parent, plan["original_failed_row"])
    if native != plan["native"] or native != report["native"]:
        raise ValueError("Recovered native predictions differ from the original completed arm")
    scores, details, score_pins = readback_scores(repo, "QfO", report, native, root,
        "ohmm_hillclimb_r1_" + plan["scoring_id"])
    pins.extend([plan["parent_snapshot"], plan["protocol"], *score_pins, *report["checked"]])
    for pin in pins:
        check(pin)
    original_plan = parent["execution_plan"]
    compact = [original_plan, plan["parent_snapshot"], plan["protocol"],
               report["execution_plan"], record(root / "results.json"),
               *native["checked"], *score_pins]
    return {"dataset": "QfO", "status": "independent_content_readback_verified",
            "scores": endpoint_vector("QfO", scores), "details": details,
            "products": native["products"], "checked": list({p["path"]: p for p in compact}.values()),
            "verification_record_count": len({p["path"] for p in pins}),
            "verification_scope": "All original-plan and producer verification records rehashed during admission; complete inventories retained in the pinned execution plan and immutable producer report, not duplicated in the score table.",
            "coverage": {k: native[k] for k in ("genes", "species", "groups", "pairs")},
            "source_revision": REVISION, "original_failed_attempt": plan["parent_snapshot"],
            "resources": {"original_native_execution": plan["original_failed_row"]["native_execution"],
                          "scoring_recovery_wall_s": report["scoring_recovery_wall_s"]},
            "execution_kind": "scoring_only_recovery_of_completed_canonical_downstream",
            "scope": "Original failed scoring attempt retained. Cached downstream inference plus separate scoring recovery is not fresh full-pipeline timing; FAS is unseeded."}


def recovery_capture(path, captures):
    identity = record(path)
    report = json.loads(path.read_text())
    check(identity)
    job = str(report.get("job_id", ""))
    if report.get("status") == "canonical_qfo_scoring_recovered_independently_verified" and job.isdigit():
        accounting = subprocess.check_output(["sacct", "-j", job, "-X", "-n", "-P",
                                             "-o", "JobID,State%30,ExitCode"], text=True)
        if accounting.strip() != job + "|COMPLETED|0:0":
            raise ValueError("Compact recovery reference requires terminal successful producer")
        reference = {"producer": identity, "scheduler_accounting": accounting,
                     "scope": "External immutable terminal producer; full payload and verification inventory remain at this hash-pinned source path."}
        write_json(captures / "canonical_recovery_reference.json", reference)
        return report, {"snapshot": identity, "reference": record(captures / "canonical_recovery_reference.json")}
    return snapshot(path, captures, "canonical_recovery.json")


def singleton_policy_cell(repo, root, plan, report, row):
    from benchmark_tools.run_canonical_downstream import commands
    from benchmark_tools.run_membership_hillclimb import compare_golden
    policy = row.get("label")
    if (row.get("dataset") != "OrthoBench" or policy not in ("first_hit", "defer_ambiguous")
            or row.get("status") != "native_and_score_independently_verified"
            or report["source_revision"] != plan["source_revision"] or plan["cpu"] != 16):
        raise ValueError("Require verified singleton-policy OrthoBench arm and frozen source")
    folder = root / policy
    preflight_path = folder / "preflight.json"
    preflight = json.loads(preflight_path.read_text())
    graph, command = commands(Path(plan["source_root"]), plan, plan["case"], folder,
                              Path(plan["checkpoint"]), plan["checkpoint_sha"])
    graph.extend(["--singleton-tie-policy", policy])
    if preflight != {"graph": graph, "native": command}:
        raise ValueError("Singleton-policy commands differ from the frozen pilot plan")
    for key, expected in (("graph_execution", graph), ("native_execution", command)):
        if row[key]["returncode"] != 0 or row[key]["command"] != expected:
            raise ValueError("Singleton-policy actual execution differs from planned command")
    graph_path = folder / "graph_metrics.json"
    graph_report = json.loads(graph_path.read_text())
    parameters = graph_report["parameters"]
    if (parameters.get("singleton_tie_policy") != policy or parameters.get("cpm_resolution") != .1
            or parameters.get("leiden_seed") != 4 or parameters.get("profile_iterations") != 1):
        raise ValueError("Actual singleton graph mechanism or unchanged settings differ")
    native = replay_native(row, {"command": command})
    if (native["genes"], native["species"]) != (251378, 12):
        raise ValueError("Singleton-policy complete gene universe differs")
    scores, details, pins = readback_scores(repo, "OrthoBench", row, native, folder, None)
    if policy == "first_hit":
        compare_golden(Path(native["phylogeny"]), plan["golden"])
        if scores != plan["control"]["scores"]:
            raise ValueError("Independent singleton-policy control score differs")
        for key in ("candidate_partition", "membership_constraints"):
            a, b = row["candidate"][key], plan["control"]["candidate"][key]
            check(a)
            check(b)
            if (a["bytes"], a["sha256"]) != (b["bytes"], b["sha256"]):
                raise ValueError("Singleton-policy candidate control differs")
    checked = [record(preflight_path), record(graph_path), *native["checked"], *pins]
    for pin in checked:
        check(pin)
    cell = {"dataset": "OrthoBench", "status": "independent_content_readback_verified",
            "scores": endpoint_vector("OrthoBench", scores), "details": details,
            "products": native["products"], "checked": checked,
            "coverage": {k: native[k] for k in ("genes", "species", "groups", "pairs")},
            "source_revision": report["source_revision"], "singleton_tie_policy": policy,
            "resources": {k: row[k] for k in ("graph_execution", "native_execution")},
            "scope": "Matched canonical accepted-hit downstream replay with selective tree reuse, not fresh HMM search or isolated timing."}
    if policy == "defer_ambiguous":
        controls = [r for r in report["arms"] if r.get("label") == "first_hit"]
        if len(controls) != 1:
            raise ValueError("Ambiguous matched singleton-policy control")
        cell["matched_control"] = singleton_policy_cell(repo, root, plan, report, controls[0])
    return cell


def build(repo, parent_path, output, continuation_root=None, shared_root=None, canonical_recovery_root=None, singleton_root=None):
    if output.exists() or all(root is None for root in (continuation_root, shared_root, canonical_recovery_root, singleton_root)):
        raise ValueError("Require a fresh output and at least one producer namespace")
    parent = json.loads(parent_path.read_text())
    check(parent["score_table"])
    methods = deepcopy(parent["methods"])
    output.mkdir()
    captures = output / "snapshots"
    captures.mkdir()
    observations = []
    for kind, root in (("continuation", continuation_root), ("shared_hit", shared_root),
                       ("canonical_recovery", canonical_recovery_root), ("singleton_policy", singleton_root)):
        if root is None:
            continue
        path = root / "results.json"
        if not path.exists():
            observations.append({"kind": kind, "path": str(path), "status": "report_missing_execution_state_not_inferred"})
            continue
        report, capture = (recovery_capture(path, captures) if kind == "canonical_recovery"
                           else snapshot(path, captures, kind + ".json"))
        plan, pins = check_plan(repo, root, report)
        observations.append({"kind": kind, "producer_status": report["status"], "snapshot": capture["snapshot"]})
        if kind == "singleton_policy":
            labels = [row.get("label") for row in report.get("arms", [])]
            if len(set(labels)) != len(labels) or any(label not in ("first_hit", "defer_ambiguous") for label in labels):
                raise ValueError("Wrong or duplicate singleton-policy arm inventory")
            rows = [row for row in report["arms"] if row["label"] == "defer_ambiguous"]
            cell = {"dataset": "OrthoBench", "status": "not_yet_reported_execution_state_not_inferred", "scores": {},
                    "source_snapshot": capture["snapshot"], "producer_job_id": report["job_id"]}
            if rows:
                row = rows[0]
                cell.update(status=row["status"], producer_status=row["status"])
                cell.update({k: row[k] for k in ("error", "error_type") if k in row})
                if row["status"] == "native_and_score_independently_verified":
                    cell.update(singleton_policy_cell(repo, root, plan, report, row))
            upsert_cell(methods, "singleton_tie:defer_ambiguous", "OrthoHMM ambiguous singleton deferral (experimental)", cell)
            for pin in pins:
                check(pin)
            if runtime() != plan["runtime"]:
                raise ValueError("Scientific runtime changed during singleton admission")
            continue
        if kind == "canonical_recovery":
            from benchmark_tools.recover_canonical_qfo_scoring import failed_row
            check(plan["parent_snapshot"])
            archived = json.loads(Path(plan["parent_snapshot"]["path"]).read_text())
            if failed_row(archived) != plan["original_failed_row"]:
                raise ValueError("Recovery archive does not preserve the original failed row")
            key = "canonical_collection:canonical"
            failed = {"dataset": "QfO", "status": "failed", "scores": {},
                      "producer_status": "failed", "source_snapshot": plan["parent_snapshot"],
                      "producer_job_id": "24240",
                      "error": plan["original_failed_row"]["error"],
                      "error_type": plan["original_failed_row"]["error_type"]}
            upsert_cell(methods, key, "OrthoHMM canonical collection (experimental)", failed)
            cell = {"dataset": "QfO", "status": report["status"], "scores": {},
                    "producer_status": report["status"], "source_snapshot": capture["snapshot"],
                    "producer_job_id": report["job_id"]}
            cell.update({k: report[k] for k in ("error", "error_type") if k in report})
            if report["status"] == "canonical_qfo_scoring_recovered_independently_verified":
                cell.update(canonical_recovery_cell(repo, plan, report, root))
            upsert_cell(methods, key, "OrthoHMM canonical collection (experimental)", cell)
            for pin in pins:
                check(pin)
            if runtime() != plan["runtime"]:
                raise ValueError("Scientific runtime changed during recovery admission")
            continue
        rows = [("QfO", report)] if kind == "continuation" else [(row["dataset"], row) for row in report["arms"]]
        if any(dataset not in DATASETS for dataset, _ in rows) or len({dataset for dataset, _ in rows}) != len(rows):
            raise ValueError("Wrong or duplicate producer dataset inventory")
        for dataset, row in rows:
            key = "full_search:cap400" if kind == "continuation" else "shared_hit:shared_hit1"
            cell = {"dataset": dataset, "status": row["status"], "producer_status": row["status"],
                    "scores": {}, "source_snapshot": capture["snapshot"], "producer_job_id": report["job_id"]}
            cell.update({k: row[k] for k in ("error", "error_type") if k in row})
            success = ("original_cap400_native_and_all_six_scores_verified" if kind == "continuation"
                       else "independent_content_and_score_readback_verified")
            if row["status"] == success:
                cell.update(continuation_cell(repo, plan, row) if kind == "continuation"
                            else shared_cell(repo, root, plan, report, row))
                if kind == "shared_hit":
                    control = {**cell["matched_control"]["cell"], "dataset": dataset}
                    upsert_cell(methods, "canonical_collection:canonical", "OrthoHMM canonical collection (experimental)", control)
            upsert_cell(methods, key, key, cell)
        for pin in pins:
            check(pin)
        if runtime() != plan["runtime"]:
            raise ValueError("Scientific runtime changed during admission")
    allowed = {"full_search:cap400", "shared_hit:shared_hit1", "canonical_collection:canonical", "singleton_tie:defer_ambiguous"}
    for previous in parent["methods"]:
        if previous["key"] not in allowed and previous != next(m for m in methods if m["key"] == previous["key"]):
            raise ValueError("Unrelated method changed")
    with (output / "scores.tsv").open("x", newline="") as stream:
        writer = csv.writer(stream, delimiter="\t", lineterminator="\n")
        writer.writerow(["method", *ENDPOINTS, "complete_eight_endpoints"])
        for method in methods:
            writer.writerow([method["label"], *(method["scores"].get(k, "NA") for k in ENDPOINTS),
                             set(method["scores"]) == set(ENDPOINTS)])
    result = deepcopy(parent)
    # Earlier reporting-code identities remain in the pinned parent, not this run.
    result.pop("reporter_revision", None)
    result["reporter_files"] = [record(Path(__file__).parent / name) for name in (
        "report_additional_hillclimb.py", "report_current_hillclimb.py", "report_hillclimb.py",
        "recover_full_qfo_search.py", "run_full_search_hillclimb.py", "run_membership_hillclimb.py",
        "run_shared_hit_pipeline.py", "recover_canonical_qfo_scoring.py", "run_canonical_downstream.py",
        "run_singleton_tie_orthobench.py")]
    result.update(status="incremental_separate_trial_readback_no_default_promotion", parent_readback=record(parent_path),
                  reporter_source=record(__file__), methods=methods, score_table=record(output / "scores.tsv"),
                  additional_trial_observations=observations, defaults_promoted=False, goal_completed=False)
    vectors = {m["key"]: m["scores"] for m in methods
               if m["key"] in {"shared_hit:shared_hit1", "canonical_collection:canonical"}}
    result["round_screening"]["shared_hit"] = screen(vectors, "canonical_collection:canonical")
    singleton = {m["key"]: m["scores"] for m in methods
                 if m["key"] in {"singleton_tie:defer_ambiguous", "canonical_collection:canonical"}}
    result["round_screening"]["singleton_tie"] = screen(singleton, "canonical_collection:canonical")
    canonical = {m["key"].split(":", 1)[1]: m["scores"] for m in methods
                 if m["key"].startswith("canonical_collection:")}
    result["round_screening"]["canonical_collection"] = screen(canonical, "arrival")
    full = {m["key"].split(":", 1)[1]: m["scores"] for m in methods if m["key"].startswith("full_search:")}
    result["round_screening"]["full_search"] = screen(full, "control100")
    result["limitations"].extend(["Newly completed trial cells are independently read back; unrelated earlier admissions are inherited, not rescored.",
        "Pending, failed, absent and unverified cells do not supply scores. Cell admission is not a claim that the producer job is terminal.",
        "Fresh shared-hit inference and cached canonical controls have different timing scopes. FAS remains unseeded."])
    write_json(output / "summary.json", result)
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("repo", "parent", "output"):
        parser.add_argument("--" + name, type=Path, required=True)
    for name in ("continuation-root", "shared-root", "canonical-recovery-root", "singleton-root"):
        parser.add_argument("--" + name, type=Path)
    args = parser.parse_args()
    result = build(args.repo.resolve(), args.parent.resolve(), args.output.resolve(),
                   args.continuation_root.resolve() if args.continuation_root else None,
                   args.shared_root.resolve() if args.shared_root else None,
                   args.canonical_recovery_root.resolve() if args.canonical_recovery_root else None,
                   args.singleton_root.resolve() if args.singleton_root else None)
    print(result["status"], len(result["methods"]))
