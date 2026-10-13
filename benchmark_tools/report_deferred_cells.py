"""Independent readback of fixed-tree and own-inferred copy-deferral cells."""

from __future__ import annotations

import json
from pathlib import Path

from benchmark_tools.prepare_ob_candidate_neighborhood import check
from benchmark_tools.report_hillclimb import (
    EXPECTED_GENES, endpoint_vector, readback_scores, replay_native, snapshot,
)
from benchmark_tools.run_membership_hillclimb import argument, compare_golden

LABELS = ("control", "deferred500")


def planned_case(plan, dataset, label):
    if dataset == "QfO":
        cases = plan["cases"]
        settings = plan
    else:
        datasets = plan["cases"]
        if tuple(c["dataset"] for c in datasets) != ("OrthoBench", "ThreeKingdoms"):
            raise ValueError("Unexpected inferred supplement dataset inventory")
        settings = next(c for c in datasets if c["dataset"] == dataset)
        cases = settings["arms"]
    if tuple(c["label"] for c in cases) != LABELS:
        raise ValueError("Unexpected copy-deferral arm inventory")
    return next(c for c in cases if c["label"] == label), settings


def audit_deferred_cell(repo, kind, dataset, row, plan, producer):
    from benchmark_tools.deferred_qfo_inferred import require_inferred_summary
    from benchmark_tools.run_held_tree_hillclimb import require_control_scores, require_reuse, tree_contrast

    if (row["status"] != "independent_content_readback_verified"
            or row.get("inputs_reverified_unchanged") is not True):
        raise ValueError("Missing copy-deferral producer postflight validation")
    if kind not in ("deferred_fixed", "deferred_inferred"):
        raise ValueError("Unknown copy-deferral comparison kind")
    case, settings = planned_case(plan, dataset, row["label"])
    argv = case["command"]
    fixed = kind == "deferred_fixed"
    if (fixed and dataset != "QfO"
            or argument(argv, "--species-tree-mode") != ("supplied" if fixed else "infer")
            or ("--species-tree" in argv) != fixed
            or argument(argv, "--root-rule") != "species_overlap"
            or argument(argv, "--pair-rule") != "positive_paralogy"
            or argument(argv, "--species-tree-rooting") != "min_variance"
            or argument(argv, "--satellite-membership-rule") != "high_confidence_pair"
            or "--unconstrained-membership" in argv):
        raise ValueError("Fixed/own-inferred tree or reconciliation protocol differs")
    if dataset == "QfO" and plan.get("evaluation_mode", "fixed_species_tree") != (
            "fixed_species_tree" if fixed else "inferred_parent"):
        raise ValueError("Fixed/own-inferred evaluation mode differs")
    if fixed and Path(argument(argv, "--species-tree")).resolve() != Path(plan["tree"]["path"]).resolve():
        raise ValueError("Fixed species-tree input differs from plan")
    checked = [producer["execution_plan"], *plan["inputs"], *plan["sources"], *plan["runtime"]["files"]]
    if not fixed and dataset == "QfO":
        receipt = producer["parent_admission"]
        check(receipt)
        resolution = json.loads(Path(receipt["path"]).read_text())
        if (resolution["status"] != "terminal_parent_native_inputs_resolved_and_verified"
                or Path(receipt["path"]).resolve() != Path(plan["output"]) / "resolved_parent_inputs.json"):
            raise ValueError("Own-inferred parent inputs were not resolved")
        from benchmark_tools.deferred_qfo_inferred import require_terminal_parent
        check(resolution["snapshot"])
        require_terminal_parent(json.loads(Path(resolution["snapshot"]["path"]).read_text()), plan)
        checked.extend([receipt, *resolution["checked"]])
    candidate = row["candidate"]
    for flag, key in (("--candidate-clusters", "candidate_partition"),
                      ("--membership-constraints", "membership_constraints")):
        identity = candidate[key]
        if Path(argument(argv, flag)).resolve() != Path(identity["path"]).resolve():
            raise ValueError("Generated candidate/constraint identity differs from command")
        checked.append(identity)
    checked.append(candidate["seed_partition"])
    if "seed" in row:
        checked.append(row["seed"])
    for item in checked:
        check(item)
    native = replay_native(row, case)
    if native["genes"] != EXPECTED_GENES[dataset] or native["species"] != (78 if dataset == "QfO" else 12):
        raise ValueError("Wrong benchmark gene/species universe")
    if fixed:
        contrast = tree_contrast(plan["tree"]["path"], Path(native["phylogeny"]) / "species_tree.rooted.nwk")
        if contrast["rooted_clade_distance"]:
            raise ValueError("Fixed species-tree handoff changed")
        if row["label"] == "control":
            require_reuse(native["summary"])
    else:
        if dataset == "QfO":
            require_inferred_summary(native["summary"])
        else:
            markers = native["summary"].get("species_tree_families")
            if isinstance(markers, bool) or not isinstance(markers, int) or markers <= 0:
                raise ValueError("Own-inferred species tree has no native marker families")
    if row["label"] == "control":
        compare_golden(Path(native["phylogeny"]), settings["golden"])
        if dataset != "QfO" or fixed:
            from benchmark_tools.run_three_kingdoms_membership_supplement import candidate_control_matches
            candidate_control_matches(Path(candidate["candidate_partition"]["path"]),
                Path(candidate["membership_constraints"]["path"]),
                settings["baseline_candidate"], settings["baseline_constraints"])
    participant = "ohmm_hillclimb_r1_" + case["scoring_id"] if dataset == "QfO" else "unused"
    scores, details, files = readback_scores(repo, dataset, row, native, Path(case["output"]), participant)
    if row["label"] == "control":
        if dataset == "QfO":
            require_control_scores(scores, settings["expected_control_scores"])
        elif endpoint_vector(dataset, scores) != settings["expected_control_scores"]:
            raise ValueError("Deterministic native control scores differ")
    checked.extend([*native["checked"], *files])
    for item in checked:
        check(item)
    return {"status": "independent_content_readback_verified", "scores": endpoint_vector(dataset, scores),
        "details": details, "coverage": {k: native[k] for k in ("genes", "species", "groups", "pairs")},
        "products": native["products"], "checked": checked, "native_summary": native["summary"],
        "tree_mode": "supplied_fixed_diagnostic" if fixed else "candidate_specific_inferred",
        "resources": {k: row[k] for k in ("native_execution", "native_and_scoring_wall_s", "candidate_wall_s") if k in row},
        "scope": "Independent content readback of cached-stage accuracy; not fresh full-pipeline timing or terminal scheduler admission."}


def collect_deferred_cells(repo, kind, dataset, root, snapshots):
    cells, captures = [], []
    top = root / "results.json"
    if not top.exists():
        return [{"method_key": kind + ":" + label, "dataset": dataset,
                 "status": "report_missing_execution_state_not_inferred", "scores": {}} for label in LABELS], captures
    producer, capture = snapshot(top, snapshots, f"{kind}_{dataset}_producer.json")
    captures.append(capture)
    plan, error = None, None
    try:
        check(producer["execution_plan"])
        plan, plan_capture = snapshot(producer["execution_plan"]["path"], snapshots, f"{kind}_{dataset}_plan.json")
        captures.append(plan_capture)
        if (plan["source_revision"] != producer["source_revision"]
                or Path(plan["output"]).resolve() != root.resolve()):
            raise ValueError("Copy-deferral producer plan identity differs")
    except Exception as problem:
        error = str(problem)
        plan = None
    report = producer
    if dataset != "QfO":
        child = root / dataset / "results.json"
        if child.exists():
            report, capture = snapshot(child, snapshots, f"{kind}_{dataset}_cells.json")
            captures.append(capture)
            if report["dataset"] != dataset:
                raise ValueError("Copy-deferral child dataset identity differs")
        else:
            report = {"arms": [], "status": "dataset_report_missing"}
    observed = set()
    control_valid = False
    for row in report["arms"]:
        label = row["label"]
        if label not in LABELS or label in observed:
            raise ValueError("Duplicate or unknown copy-deferral arm")
        observed.add(label)
        cell = {"method_key": kind + ":" + label, "dataset": dataset, "scores": {},
            "status": row["status"], "producer_status": row["status"],
            "producer_job_id": producer.get("job_id"), "producer_scope": producer.get("scope"),
            "source_snapshot": capture["snapshot"], "claimed_scores": row.get("scores", {})}
        if row["status"] == "independent_content_readback_verified":
            try:
                if plan is None:
                    raise ValueError(error)
                if label != "control" and not control_valid:
                    raise ValueError("Independent matched-control admission failed or is missing")
                cell.update(audit_deferred_cell(repo, kind, dataset, row, plan, producer))
                if label == "control":
                    control_valid = True
            except Exception as problem:
                cell.update(status="readback_failed", error_type=type(problem).__name__, error=str(problem))
        cells.append(cell)
    cells.extend({"method_key": kind + ":" + label, "dataset": dataset, "scores": {},
                  "status": "not_yet_reported_execution_state_not_inferred"} for label in LABELS if label not in observed)
    return cells, captures
