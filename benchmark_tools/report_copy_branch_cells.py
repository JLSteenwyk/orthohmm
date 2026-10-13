"""Independently read copy-branch ablations, keeping structural reuse explicit."""

from collections import Counter
import json
from pathlib import Path
import subprocess

from benchmark_tools import run_copy_branch_ablation as runner
from benchmark_tools.prepare_ob_candidate_neighborhood import check
from benchmark_tools.report_deferred_cells import audit_deferred_cell
from benchmark_tools.report_hillclimb import EXPECTED_GENES, endpoint_vector, readback_scores, replay_native, snapshot
from benchmark_tools.run_membership_hillclimb import argument
from benchmark_tools.run_three_kingdoms_membership_supplement import candidate_control_matches

LABELS = (*runner.PARENT_LABELS, *(label for label, _ in runner.ARMS))


def check_candidate(candidate):
    for key in ("seed_partition", "candidate_partition", "membership_constraints"):
        check(candidate[key])


def trace_counts(row, dataset):
    check(row["copy_split_trace"])
    events = json.loads(Path(row["copy_split_trace"]["path"]).read_text())["events"]
    counts = Counter()
    scope = {"control": None, "deferred500": "both", "large_only500": "large", "medium_only500": "medium"}[row["label"]]
    for event in events:
        members = event["members"]
        if (not isinstance(members, list) or members != sorted(set(members))
                or any(type(g) is not int or not 0 <= g < EXPECTED_GENES[dataset] for g in members)
                or event["size"] != len(members)
                or type(event["maximum_species_copies"]) is not int
                or not 1 <= event["maximum_species_copies"] <= len(members)
                or event["dataset_species_count"] != (78 if dataset == "QfO" else 12)
                or any(type(event[key]) is not bool for key in ("large_copy_eligible", "medium_copy_eligible", "deferred"))
                or not (event["large_copy_eligible"] or event["medium_copy_eligible"])):
            raise ValueError("Invalid post-merge copy-split event")
        selected = scope == "both" or scope == "large" and event["large_copy_eligible"] or scope == "medium" and event["medium_copy_eligible"]
        if event["deferred"] != bool(selected and event["size"] <= 500):
            raise ValueError("Executed deferral differs from the prespecified branch/ceiling")
        if dataset != "QfO" and event["large_copy_eligible"]:
            raise ValueError("Large-copy branch cannot activate on a 12-species panel")
        for branch in ("large", "medium"):
            if event[branch + "_copy_eligible"]:
                counts[branch + "_eligible_families"] += 1
                if event["deferred"]:
                    counts[branch + "_deferred_families"] += 1
                    counts[branch + "_deferred_genes"] += event["size"]
    if dict(counts) != row["copy_split_counts"]:
        raise ValueError("Claimed executed copy-split counts disagree with the trace")
    return dict(counts)


def admit_parent(repo, plan, producer):
    if Path(plan["repo"]).resolve() != repo.resolve():
        raise ValueError("Branch experiment belongs to another repository")
    runner.verify(plan, [producer["execution_plan"]])
    receipt = producer["parent_admission"]
    if Path(receipt["path"]).resolve() != Path(plan["output"]) / "resolved_parent_inputs.json":
        raise ValueError("Wrong branch parent admission receipt")
    check(receipt)
    resolution = json.loads(Path(receipt["path"]).read_text())
    if resolution["status"] != "terminal_matched_parent_independently_admitted":
        raise ValueError("Parent was not independently admitted")
    checked = [receipt, *resolution["checked"]]
    for capture in resolution["captures"]:
        checked.append(capture["snapshot"])
    for item in checked:
        check(item)
    accounting = subprocess.check_output(["sacct", "-j", plan["parent_job"], "--parsable2",
        "--format=JobIDRaw,State,ExitCode,Elapsed,NodeList,AllocCPUS"], text=True)
    scheduler = runner.terminal_job(accounting, plan["parent_job"])
    if scheduler != resolution["scheduler"]:
        raise ValueError("Terminal parent scheduler identity changed")
    parent = json.loads(Path(plan["parent_plan"]["path"]).read_text())
    top = json.loads((Path(plan["parent_root"]) / "results.json").read_text())
    child = json.loads(Path(plan["parent_results"]).read_text())
    if (top.get("job_id") != plan["parent_job"] or top.get("execution_plan") != plan["parent_plan"]
            or top.get("source_revision") != parent["source_revision"]
            or not top.get("status", "").startswith(("finished", "failed"))
            or resolution["parent_rows"] != child["arms"]
            or tuple(r["label"] for r in child["arms"]) != runner.PARENT_LABELS
            or plan["parent_cases"] != runner.cases_for(parent, plan["dataset"])):
        raise ValueError("Wrong terminal matched-parent attempt or arm inventory")
    if child["arms"][0]["status"] != "independent_content_readback_verified":
        raise ValueError("Missing verified matched control")
    if plan["dataset"] != "QfO":
        if child["status"] != "finished_all_arms_verified":
            raise ValueError("Narrow-panel reuse requires both verified parent arms")
    elif child["status"] == "finished_all_arms_verified":
        if scheduler["State"] != "COMPLETED":
            raise ValueError("Successful QfO parent disagrees with scheduler")
    elif (child["status"] != "finished_with_explicit_gaps" or scheduler["State"] != "FAILED"
            or child["arms"][1]["status"] != "failed" or not child["arms"][1].get("error")):
        raise ValueError("Require a specific terminal QfO variant failure, not an unexplained gap")
    equivalence = runner.replay_equivalence(Path(parent["source_root"]), Path(plan["source_root"]))
    if equivalence != plan["replay_equivalence"]:
        raise ValueError("Downstream reuse equivalence differs from the frozen plan")
    admitted = {}
    for row in child["arms"]:
        if row["status"] == "independent_content_readback_verified":
            admitted[row["label"]] = audit_deferred_cell(repo, "deferred_inferred", plan["dataset"], row, parent, top)
            checked.extend(admitted[row["label"]]["checked"])
    return {"resolution": resolution, "parent": parent, "producer": top, "rows": child["arms"],
            "admitted": admitted, "scheduler": scheduler, "checked": checked}


def audit_structural_controls(plan, producer):
    rows = producer["structural_controls"]
    if set(rows) != set(runner.PARENT_LABELS):
        raise ValueError("Missing exact structural controls")
    for index, label in enumerate(runner.PARENT_LABELS):
        row = rows[label]
        if row["label"] != label:
            raise ValueError("Structural control label differs")
        check(row["seed"])
        check_candidate(row["candidate"])
        trace_counts(row, plan["dataset"])
        candidate = row["candidate"]
        candidate_control_matches(Path(candidate["candidate_partition"]["path"]),
            Path(candidate["membership_constraints"]["path"]),
            plan["parent_candidates"][index]["candidate"], plan["parent_candidates"][index]["constraints"])
    if "golden_refinement" in plan["data"] and rows["control"]["seed"]["sha256"] != plan["data"]["golden_refinement"]["sha256"]:
        raise ValueError("Default seed differs from the admitted refinement")
    return rows


def audit_branch_cell(repo, row, plan, producer, parent):
    if row.get("inputs_reverified_unchanged") is not True:
        raise ValueError("Missing branch postflight verification")
    index = tuple(label for label, _ in runner.ARMS).index(row["label"])
    arm = plan["arms"][index]
    if (row["profile"] != runner.ARMS[index][1]
            or (arm["label"], arm["profile"]) != runner.ARMS[index]):
        raise ValueError("Branch profile differs from the prespecified arm")
    controls = audit_structural_controls(plan, producer)
    check(row["seed"])
    check_candidate(row["candidate"])
    counts = trace_counts(row, plan["dataset"])
    if row["candidate"]["expansion"]["parameters"] != controls["control"]["candidate"]["expansion"]["parameters"]:
        raise ValueError("Branch satellite expansion parameters changed")
    if plan["dataset"] == "QfO":
        if row["execution_kind"] != "new_cached_graph_branch_inference_own_tree" or row["native_execution"]["returncode"] != 0:
            raise ValueError("New QfO branch inference did not succeed")
        case, native_row, scoring_row = arm, row, row
        folder, participant = Path(arm["output"]), "ohmm_hillclimb_r1_" + arm["scoring_id"]
        fixed = json.loads(Path(parent["parent"]["parent_plan"]["path"]).read_text())
        checkpoint = Path(fixed["cases"][1]["output"]) / "output"
        candidate_path = folder / "candidate/orthohmm_working_res/orthohmm_edges_clustered.txt"
        expected = runner.native_command(plan["parent_cases"][0]["command"], Path(plan["source_root"]), folder,
            candidate_path, candidate_path.with_name("phylogeny_candidate_merges.json"), checkpoint)
        if case["command"] != expected:
            raise ValueError("QfO branch command changed the matched own-tree inference protocol")
        for flag, key in (("--candidate-clusters", "candidate_partition"), ("--membership-constraints", "membership_constraints")):
            if Path(argument(case["command"], flag)).resolve() != Path(row["candidate"][key]["path"]).resolve():
                raise ValueError("Branch candidate identity differs from inference command")
    else:
        if (row["execution_kind"] != "exact_candidate_structural_reuse_no_native_rerun"
                or row["reused_parent_label"] != runner.PARENT_LABELS[index]
                or any(key in row for key in ("native_execution", "native_and_scoring_wall_s"))):
            raise ValueError("Narrow-panel structural reuse mislabeled as new inference")
        reference = controls[runner.PARENT_LABELS[index]]["candidate"]
        candidate_control_matches(Path(row["candidate"]["candidate_partition"]["path"]),
            Path(row["candidate"]["membership_constraints"]["path"]),
            reference["candidate_partition"], reference["membership_constraints"])
        case, native_row = plan["parent_cases"][index], parent["rows"][index]
        scoring_row = native_row
        folder, participant = Path(case["output"]), "unused"
    if row["native_case"] != case or Path(row["scoring_output"]).resolve() != folder.resolve():
        raise ValueError("Wrong branch native or scoring handoff")
    native = replay_native(native_row, case)
    if native["genes"] != EXPECTED_GENES[plan["dataset"]] or native["species"] != (78 if plan["dataset"] == "QfO" else 12):
        raise ValueError("Branch inference lost the benchmark gene/species universe")
    runner.require_markers(native["summary"])
    if row["native_products"] != native["products"] or row["native_summary"] != native["summary"]:
        raise ValueError("Claimed branch native products or summary differ")
    scores, details, files = readback_scores(repo, plan["dataset"], scoring_row, native, folder, participant)
    if scores != row["scores"]:
        raise ValueError("Independent branch scores differ from producer claims")
    checked = [*parent["checked"], *native["checked"], *files, row["seed"], row["copy_split_trace"],
               *(row["candidate"][key] for key in ("seed_partition", "candidate_partition", "membership_constraints"))]
    runner.verify(plan, checked)
    return {"status": "independent_content_readback_verified", "scores": endpoint_vector(plan["dataset"], scores),
        "details": details, "products": native["products"], "checked": checked,
        "coverage": {key: native[key] for key in ("genes", "species", "groups", "pairs")},
        "native_summary": native["summary"], "execution_kind": row["execution_kind"],
        "copy_split_counts": counts, "parent_scheduler": parent["scheduler"],
        "resources": {key: row[key] for key in ("candidate_wall_s", "native_execution") if key in row},
        "scope": "Structural reuse is not a new inference timing. New QfO inference uses cached HMM graph inputs and own inferred trees; shared-host timing is confounded."}


def collect_branch_cells(repo, dataset, root, snapshots):
    path = root / "results.json"
    def missing(label, status, **extra):
        return {"method_key": "copy_branch:" + label, "dataset": dataset, "scores": {}, "status": status, **extra}
    if not path.exists():
        return [missing(label, "report_missing_execution_state_not_inferred") for label in LABELS], []
    producer, capture = snapshot(path, snapshots, f"copy_branch_{dataset}_producer.json")
    captures, cells = [capture], []
    plan, parent, error = None, None, None
    try:
        check(producer["execution_plan"])
        plan, plan_capture = snapshot(producer["execution_plan"]["path"], snapshots, f"copy_branch_{dataset}_plan.json")
        captures.append(plan_capture)
        if (plan["source_revision"] != producer["source_revision"] or plan["dataset"] != dataset
                or producer["dataset"] != dataset or Path(plan["output"]).resolve() != root.resolve()
                or tuple((arm["label"], arm["profile"]) for arm in plan["arms"]) != runner.ARMS):
            raise ValueError("Branch producer/execution plan identity differs")
        if "parent_admission" in producer:
            parent = admit_parent(repo, plan, producer)
    except Exception as problem:
        error = str(problem)
    for label in runner.PARENT_LABELS:
        cell = missing(label, "parent_not_admitted_execution_state_not_inferred", source_snapshot=capture["snapshot"])
        if error:
            cell.update(status="readback_failed", error=error)
        elif parent is not None:
            if label in parent["admitted"]:
                cell.update(parent["admitted"][label], execution_kind="independently_reused_matched_parent",
                            parent_scheduler=parent["scheduler"])
            else:
                old = parent["rows"][runner.PARENT_LABELS.index(label)]
                cell.update(status=old["status"], **{key: old[key] for key in ("error", "error_type") if key in old})
        cells.append(cell)
    observed = set()
    for row in producer["arms"]:
        label = row["label"]
        if label not in LABELS[2:] or label in observed:
            raise ValueError("Duplicate or unknown branch arm")
        observed.add(label)
        cell = missing(label, row["status"], producer_status=row["status"], producer_job_id=producer.get("job_id"),
                       source_snapshot=capture["snapshot"], claimed_scores=row.get("scores", {}),
                       **{key: row[key] for key in ("error", "error_type") if key in row})
        if row["status"] == "independent_content_readback_verified":
            try:
                if error or parent is None:
                    raise ValueError(error or "No independently admitted matched parent")
                cell.update(audit_branch_cell(repo, row, plan, producer, parent))
            except Exception as problem:
                cell.update(status="readback_failed", error_type=type(problem).__name__, error=str(problem))
        cells.append(cell)
    terminal = producer.get("status", "").startswith(("failed", "finished"))
    cells.extend(missing(label, "not_run_terminal_producer_gap" if terminal else "not_yet_reported_execution_state_not_inferred")
                 for label in LABELS[2:] if label not in observed)
    return cells, captures
