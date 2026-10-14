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


def shared_qfo_cell(repo, root, plan, report):
    from benchmark_tools.continue_shared_qfo import REVISION, original_skip, relocated_case, terminal_dependencies
    if (report.get("status") != "original_shared_qfo_native_and_all_six_scores_verified"
            or plan["scientific_revision"] != REVISION or plan["cpu"] != 16
            or plan["dependency_jobs"] != ["24242", "24246"]):
        raise ValueError("Require successful original-science QfO sensitivity continuation")
    terminal_dependencies(report["dependency_accounting"])
    check(plan["parent_plan"])
    parent = json.loads(Path(plan["parent_plan"]["path"]).read_text())
    if (Path(plan["parent_plan"]["path"]).resolve() != Path(plan["parent_root"]) / "execution_plan.json"
            or plan["parent_plan"] not in plan["inputs"]
            or parent["source_root"] != plan["scientific_source"]
            or relocated_case(parent, root) != plan["case"]):
        raise ValueError("Continuation changes the original scientific case or source")
    check(report["original_skip"])
    if Path(report["original_skip"]["path"]).resolve() != root / "original_skip_snapshot.json":
        raise ValueError("Original skipped attempt belongs to another continuation")
    archived = json.loads(Path(report["original_skip"]["path"]).read_text())
    skipped = original_skip(archived, parent, plan["parent_plan"])
    identifier = "sh_" + hashlib.sha256(str(root).encode()).hexdigest()[:12]
    if (plan["scoring_id"] != identifier or report["execution"]["returncode"] != 0
            or report["execution"]["command"] != plan["case"]["command"]):
        raise ValueError("Continuation actual command or scoring namespace differs")
    recovery_root = Path(plan["recovery_root"])
    control_pin = report["control"]["report"]
    check(control_pin)
    if Path(control_pin["path"]).resolve() != recovery_root / "results.json":
        raise ValueError("Recovered control belongs to another producer")
    recovery = json.loads(Path(control_pin["path"]).read_text())
    recovery_plan, recovery_pins = check_plan(repo, recovery_root, recovery)
    if recovery.get("job_id") != "24246" or recovery["execution_plan"] not in plan["inputs"]:
        raise ValueError("Recovered control differs from the planned dependency")
    control = canonical_recovery_cell(repo, recovery_plan, recovery, recovery_root)
    for key in ("scores", "products", "coverage", "source_revision"):
        if report["control"]["cell"][key] != control[key]:
            raise ValueError("Independent recovered control content differs")
    native = shared_native(plan["case"], parent)
    if (native["genes"], native["species"]) != (984137, 78) or (
            native["genes"], native["species"]) != (control["coverage"]["genes"], control["coverage"]["species"]):
        raise ValueError("QfO sensitivity continuation gene universe differs")
    folder = Path(plan["case"]["output"])
    phylo = folder / "output/orthohmm_phylogeny"
    native.update(phylogeny=str(phylo), prediction=str(folder / "output/orthohmm_working_res/orthohmm_edges_clustered.txt"),
                  checked=[*native["products"].values(), record(folder / "metrics.json"),
                           record(phylo / "provenance_manifest.json"), record(phylo / "reconciliation_summary.json"),
                           *[record(p) for p in sorted((folder / "output/orthohmm_working_res/high_sensitivity_checkpoint").iterdir())]])
    cell = fresh_cell(repo, "QfO", report, native, root / "scoring",
                      "ohmm_hillclimb_r1_" + identifier, REVISION)
    checked = [plan["parent_plan"], report["original_skip"], control_pin, *recovery_pins,
               *control["checked"], *report["checked"]]
    for pin in checked:
        check(pin)
    compact = [plan["parent_plan"], report["original_skip"], control_pin,
               recovery["execution_plan"], *control["checked"]]
    cell["checked"].extend(list({pin["path"]: pin for pin in compact}.values()))
    cell["verification_record_count"] = len({pin["path"] for pin in checked})
    cell["verification_scope"] = "Full producer and recovery inventories rehashed; retained in pinned producer snapshots rather than duplicated in the score table."
    cell["matched_control"] = {"report": control_pin, "cell": control}
    cell["original_unrun_attempt"] = {**skipped, "scores": {}, "source_snapshot": report["original_skip"],
                                    "producer_job_id": archived["job_id"]}
    cell["execution_kind"] = "fresh_original_canonical_shared_hit1_cap100_continuation"
    return cell


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


def singleton_followup_plan(repo, root, report):
    from benchmark_tools.run_singleton_tie_followup import REVISION, scoring_id, selected_case
    check(report["execution_plan"])
    plan = json.loads(Path(report["execution_plan"]["path"]).read_text())
    check(plan["prepared_inputs"])
    prepared = json.loads(Path(plan["prepared_inputs"]["path"]).read_text())
    case = selected_case(prepared, plan["dataset"])
    if (Path(plan["repo"]).resolve() != repo or Path(plan["output"]).resolve() != root
            or plan["case"] != case or plan["scientific_source"] != prepared["scientific_source"]
            or plan["scientific_revision"] != REVISION or report["source_revision"] != REVISION
            or plan["driver_revision"] != report["driver_revision"] or plan["cpu"] != 16
            or plan["runtime"] != prepared["runtime"] or plan["tool_path"] != prepared["tool_path"]
            or plan["scoring_id"] != scoring_id(root)):
        raise ValueError("Singleton follow-up changes the fixed scientific plan")
    captures = {item["original"]["path"]: item for item in plan["frozen_reporting_evidence"]}
    if len(captures) != len(plan["frozen_reporting_evidence"]):
        raise ValueError("Duplicate frozen reporting evidence")
    inventory = {pin["path"]: (pin["bytes"], pin["sha256"]) for pin in plan["inputs"]}
    if len(inventory) != len(plan["inputs"]):
        raise ValueError("Duplicate singleton follow-up verification records")
    for required in [prepared["preparer_source"], *prepared["inputs"]]:
        captured = captures.get(required["path"])
        pin = required
        if captured:
            original, frozen = captured["original"], captured["frozen"]
            path = Path(required["path"])
            if (original != required or path.parent != repo / "benchmark_tools" or path.suffix != ".py"
                    or Path(frozen["path"]) != Path(plan["driver_root"]) / "evidence/benchmark_tools" / path.name
                    or (original["bytes"], original["sha256"]) != (frozen["bytes"], frozen["sha256"])):
                raise ValueError("Invalid frozen reporting evidence substitution")
            pin = frozen
        if inventory.get(pin["path"]) != (pin["bytes"], pin["sha256"]):
            raise ValueError("Singleton follow-up verification inventory is incomplete")
    pins = [report["execution_plan"], plan["prepared_inputs"], *plan["inputs"], *plan["driver_sources"],
            *plan["runtime"]["files"]]
    for pin in pins:
        check(pin)
    if runtime() != plan["runtime"]:
        raise ValueError("Scientific runtime differs from the singleton follow-up")
    return plan, pins


def singleton_followup_cell(repo, root, plan, report, row):
    from benchmark_tools.run_canonical_downstream import commands
    from benchmark_tools.run_membership_hillclimb import compare_golden
    from benchmark_tools.run_singleton_tie_followup import UNIVERSES
    dataset, policy = row.get("dataset"), row.get("label")
    if (dataset != plan["dataset"] or dataset not in UNIVERSES or policy not in ("first_hit", "defer_ambiguous")
            or row.get("status") != "native_and_score_independently_verified"):
        raise ValueError("Require verified singleton follow-up arm")
    case, folder = plan["case"], root / policy
    graph, command = commands(Path(plan["scientific_source"]), plan, case["case"], folder,
                              Path(case["checkpoint"]), case["checkpoint_sha"])
    graph.extend(["--singleton-tie-policy", policy])
    preflight_path = folder / "preflight.json"
    if json.loads(preflight_path.read_text()) != {"graph": graph, "native": command}:
        raise ValueError("Singleton follow-up preflight commands differ")
    for key, expected in (("graph_execution", graph), ("native_execution", command)):
        if row[key]["returncode"] != 0 or row[key]["command"] != expected:
            raise ValueError("Singleton follow-up actual execution differs")
    graph_path = folder / "graph_metrics.json"
    parameters = json.loads(graph_path.read_text())["parameters"]
    expected = {"singleton_tie_policy": policy, "cpm_resolution": .1, "leiden_seed": 4,
                "profile_iterations": 1, "profile_min_species": 1, "matrix": "BLOSUM62"}
    if any(parameters.get(key) != value for key, value in expected.items()):
        raise ValueError("Actual singleton follow-up graph settings differ")
    native = replay_native(row, {"command": command})
    if (native["genes"], native["species"]) != UNIVERSES[dataset]:
        raise ValueError("Singleton follow-up gene universe differs")
    scoring_folder, participant, extra = folder, None, []
    if dataset == "QfO":
        if policy == "first_hit":
            scoring_folder = repo / "benchmarks/results/canonical_qfo_scoring_recovery_r1_v1_20261014"
            pin = row["control_scoring_reuse"]
            check(pin)
            recovery_path = scoring_folder / "execution_plan.json"
            recovery_pin = record(recovery_path)
            if (Path(pin["path"]).resolve() != scoring_folder / "results.json"
                    or pin not in plan["inputs"] or recovery_pin not in plan["inputs"]):
                raise ValueError("Reused QfO assessment is not bound to the follow-up")
            recovery = json.loads(recovery_path.read_text())
            participant = "ohmm_hillclimb_r1_" + recovery["scoring_id"]
            extra = [pin, recovery_pin]
        else:
            participant = "ohmm_hillclimb_r1_" + plan["scoring_id"]
    scores, details, pins = readback_scores(repo, dataset, row, native, scoring_folder, participant)
    vector = endpoint_vector(dataset, scores)
    if policy == "first_hit":
        compare_golden(Path(native["phylogeny"]), case["golden"])
        if vector != case["control"]["scores"]:
            raise ValueError("Independent singleton follow-up control scores differ")
        for key in ("candidate_partition", "membership_constraints"):
            a, b = row["candidate"][key], case["candidate"][key]
            check(a)
            check(b)
            if (a["bytes"], a["sha256"]) != (b["bytes"], b["sha256"]):
                raise ValueError("Singleton follow-up candidate control differs")
    checked = [record(preflight_path), record(graph_path), *native["checked"], *pins, *extra]
    for pin in checked:
        check(pin)
    cell = {"dataset": dataset, "status": "independent_content_readback_verified", "scores": vector,
            "details": details, "products": native["products"], "checked": checked,
            "coverage": {key: native[key] for key in ("genes", "species", "groups", "pairs")},
            "source_revision": plan["scientific_revision"], "driver_revision": plan["driver_revision"],
            "singleton_tie_policy": policy,
            "resources": {key: row[key] for key in ("graph_execution", "native_execution")},
            "scope": "Matched cached canonical accepted-hit downstream contrast. QfO control assessment reused after native equivalence; variant FAS unseeded. Not fresh search or isolated timing."}
    if policy == "defer_ambiguous":
        controls = [arm for arm in report["arms"] if arm.get("label") == "first_hit"]
        if len(controls) != 1:
            raise ValueError("Ambiguous singleton follow-up matched control")
        cell["matched_control"] = singleton_followup_cell(repo, root, plan, report, controls[0])
    return cell


def build(repo, parent_path, output, continuation_root=None, shared_root=None, canonical_recovery_root=None,
          singleton_root=None, shared_qfo_root=None, singleton_qfo_root=None, singleton_tk_root=None):
    if output.exists() or all(root is None for root in (
            continuation_root, shared_root, canonical_recovery_root, singleton_root, shared_qfo_root,
            singleton_qfo_root, singleton_tk_root)):
        raise ValueError("Require a fresh output and at least one producer namespace")
    parent = json.loads(parent_path.read_text())
    check(parent["score_table"])
    methods = deepcopy(parent["methods"])
    output.mkdir()
    captures = output / "snapshots"
    captures.mkdir()
    observations = []
    for kind, root in (("continuation", continuation_root), ("shared_hit", shared_root),
                       ("canonical_recovery", canonical_recovery_root), ("singleton_policy", singleton_root),
                       ("shared_qfo", shared_qfo_root), ("singleton_qfo", singleton_qfo_root),
                       ("singleton_tk", singleton_tk_root)):
        if root is None:
            continue
        path = root / "results.json"
        if not path.exists():
            observations.append({"kind": kind, "path": str(path), "status": "report_missing_execution_state_not_inferred"})
            continue
        report, capture = (recovery_capture(path, captures) if kind == "canonical_recovery"
                           else snapshot(path, captures, kind + ".json"))
        followup = kind in ("singleton_qfo", "singleton_tk")
        plan, pins = (singleton_followup_plan(repo, root, report) if followup else check_plan(repo, root, report))
        observations.append({"kind": kind, "producer_status": report["status"], "snapshot": capture["snapshot"]})
        if kind == "singleton_policy" or followup:
            dataset = "QfO" if kind == "singleton_qfo" else "ThreeKingdoms" if kind == "singleton_tk" else "OrthoBench"
            labels = [row.get("label") for row in report.get("arms", [])]
            if (len(set(labels)) != len(labels) or any(label not in ("first_hit", "defer_ambiguous") for label in labels)
                    or any(row.get("dataset") != dataset for row in report.get("arms", []))
                    or (followup and plan["dataset"] != dataset)):
                raise ValueError("Wrong or duplicate singleton-policy arm inventory")
            rows = [row for row in report["arms"] if row["label"] == "defer_ambiguous"]
            cell = {"dataset": dataset, "status": "not_yet_reported_execution_state_not_inferred", "scores": {},
                    "source_snapshot": capture["snapshot"], "producer_job_id": report["job_id"]}
            if rows:
                row = rows[0]
                cell.update(status=row["status"], producer_status=row["status"])
                cell.update({k: row[k] for k in ("error", "error_type") if k in row})
                if row["status"] == "native_and_score_independently_verified":
                    cell.update(singleton_followup_cell(repo, root, plan, report, row) if followup
                                else singleton_policy_cell(repo, root, plan, report, row))
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
        rows = [("QfO", report)] if kind in ("continuation", "shared_qfo") else [(row["dataset"], row) for row in report["arms"]]
        if any(dataset not in DATASETS for dataset, _ in rows) or len({dataset for dataset, _ in rows}) != len(rows):
            raise ValueError("Wrong or duplicate producer dataset inventory")
        for dataset, row in rows:
            key = "full_search:cap400" if kind == "continuation" else "shared_hit:shared_hit1"
            cell = {"dataset": dataset, "status": row["status"], "producer_status": row["status"],
                    "scores": {}, "source_snapshot": capture["snapshot"], "producer_job_id": report["job_id"]}
            cell.update({k: row[k] for k in ("error", "error_type") if k in row})
            success = ("original_cap400_native_and_all_six_scores_verified" if kind == "continuation"
                       else "original_shared_qfo_native_and_all_six_scores_verified" if kind == "shared_qfo"
                       else "independent_content_and_score_readback_verified")
            if row["status"] == success:
                cell.update(continuation_cell(repo, plan, row) if kind == "continuation"
                            else shared_qfo_cell(repo, root, plan, report) if kind == "shared_qfo"
                            else shared_cell(repo, root, plan, report, row))
                if kind in ("shared_hit", "shared_qfo"):
                    control = {**cell["matched_control"]["cell"], "dataset": dataset}
                    upsert_cell(methods, "canonical_collection:canonical", "OrthoHMM canonical collection (experimental)", control)
                if kind == "shared_qfo":
                    upsert_cell(methods, key, key, cell["original_unrun_attempt"])
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
        "run_singleton_tie_orthobench.py", "continue_shared_qfo.py", "run_singleton_tie_followup.py")]
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
    for name in ("continuation-root", "shared-root", "canonical-recovery-root", "singleton-root", "shared-qfo-root",
                 "singleton-qfo-root", "singleton-three-kingdoms-root"):
        parser.add_argument("--" + name, type=Path)
    args = parser.parse_args()
    result = build(args.repo.resolve(), args.parent.resolve(), args.output.resolve(),
                   args.continuation_root.resolve() if args.continuation_root else None,
                   args.shared_root.resolve() if args.shared_root else None,
                   args.canonical_recovery_root.resolve() if args.canonical_recovery_root else None,
                   args.singleton_root.resolve() if args.singleton_root else None,
                   args.shared_qfo_root.resolve() if args.shared_qfo_root else None,
                   args.singleton_qfo_root.resolve() if args.singleton_qfo_root else None,
                   args.singleton_three_kingdoms_root.resolve() if args.singleton_three_kingdoms_root else None)
    print(result["status"], len(result["methods"]))
