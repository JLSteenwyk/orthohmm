"""Resolve completed fixed-tree candidates/raw trees for native own-tree inference."""

from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys

from orthohmm.files import fetch_fasta_files
from benchmark_tools.prepare_ob_candidate_neighborhood import check, record
from benchmark_tools.report_hillclimb import readback_scores, replay_native
from benchmark_tools.run_held_tree_hillclimb import require_control_scores, require_reuse, source_records
from benchmark_tools.run_membership_hillclimb import PRODUCTS, argument, write_json
from benchmark_tools.run_search_recovery import runtime_records
from benchmark_tools.verify_ygob_validation import require_completed_job

LABELS = ("control", "deferred500")
PROTOCOL = "hillclimb_deferred_qfo_inferred_20261013.json"


def inferred_command(original, source, folder, cpu):
    command = list(original)
    if (argument(command, "--species-tree-mode") != "supplied" or command.count("--species-tree") != 1
            or argument(command, "--satellite-membership-rule") != "high_confidence_pair"
            or "--unconstrained-membership" in command):
        raise ValueError("Require the unchanged fixed-tree direct-policy parent")
    argument(command, "--candidate-clusters")
    argument(command, "--membership-constraints")
    checkpoint = argument(command, "--output-directory")
    index = command.index("--species-tree")
    del command[index:index + 2]
    command[:2] = [sys.executable, str(source / "benchmark_tools/replay_phylogeny.py")]
    for flag, value in (("--species-tree-mode", "infer"), ("--cpu", cpu),
                        ("--output-directory", folder / "output"), ("--json", folder / "metrics.json"),
                        ("--checkpoint-source", checkpoint)):
        argument(command, flag)
        command[command.index(flag) + 1] = str(value)
    return command


def prepare(repo, source, output, parent_plan, parent_results, cpu, parent_job):
    if (output.exists() or isinstance(cpu, bool) or not 1 <= cpu <= 16
            or not str(parent_job).isdigit() or int(parent_job) <= 0):
        raise ValueError("Require fresh output, 1-16 CPUs and an explicit numeric parent job")
    subprocess.run(["git", "-C", str(source), "diff", "--exit-code", "HEAD", "--",
                    "orthohmm", "benchmark_tools"], check=True)
    parent = json.loads(parent_plan.read_text())
    if tuple(c["label"] for c in parent["cases"]) != LABELS:
        raise ValueError("Unexpected parent arm inventory")
    if parent_results != Path(parent["output"]) / "results.json":
        raise ValueError("Parent report does not belong to the supplied execution plan")
    inputs = [record(parent_plan), record(source / "benchmark_tools" / PROTOCOL), parent["tree"],
              *parent["golden"].values(), parent["baseline_candidate"], parent["baseline_constraints"]]
    fasta = Path(parent["fasta_directory"])
    inputs.extend(record(fasta / name) for name in sorted(fetch_fasta_files(str(fasta))))
    cases = []
    nonce = hashlib.sha256(str(output).encode()).hexdigest()[:16]
    for i, case in enumerate(parent["cases"]):
        for flag in ("--candidate-clusters", "--membership-constraints"):
            inputs.append(record(Path(argument(case["command"], flag))))
        inputs.append(record(Path(case["output"]) / "refined_seed.txt"))
        for flag in ("--aligner", "--tree-builder"):
            tool = argument(case["command"], flag)
            inputs.append(record(shutil.which(tool, path=parent["tool_path"]) or tool))
        folder = output / case["label"]
        cases.append({"label": case["label"], "output": str(folder),
                      "command": inferred_command(case["command"], source, folder, cpu),
                      "scoring_id": f"deferred_own_{nonce}_{i}"})
    inputs = list({item["path"]: item for item in inputs}.values())
    for item in inputs:
        check(item)
    result = {"status": "prepared_waiting_for_terminal_parent", "evaluation_mode": "inferred_parent",
        "repo": str(repo), "source_root": str(source), "output": str(output), "cpu": cpu,
        "parent_plan": record(parent_plan), "parent_results": str(parent_results), "parent_job": str(parent_job),
        "inputs": inputs, "sources": [*source_records(source), record(source / "benchmark_tools/run_deferred_copy_split.slurm")],
        "runtime": runtime_records(), "tool_path": parent["tool_path"], "cases": cases,
        "source_revision": subprocess.check_output(["git", "-C", str(source), "rev-parse", "HEAD"], text=True).strip(),
        "tree": parent["tree"], "golden": {**parent["golden"], "species_tree.rooted.nwk": parent["tree"]},
        "expected_control_scores": parent["expected_control_scores"],
        "scope": "Development-exposed cached HMM-derived candidates/raw gene trees with candidate-specific native inferred species trees; not fresh pipeline timing.",
        "deferred_inputs": "Terminal parent report, independently verified parent native products/assessments, and matching completed native alignment/raw-tree checkpoints. Resolve and freeze before inference."}
    output.mkdir(parents=True)
    write_json(output / "execution_plan.json", result)
    return result


def require_terminal_parent(report, plan):
    if (report.get("status") != "finished_all_arms_verified" or report.get("job_id") != plan["parent_job"]
            or report.get("execution_plan") != plan["parent_plan"]
            or tuple(r.get("label") for r in report.get("arms", [])) != LABELS
            or any(r.get("status") != "independent_content_readback_verified"
                   or r.get("inputs_reverified_unchanged") is not True for r in report["arms"])):
        raise ValueError("Require successful terminal parent with both independently verified native arms")


def require_inferred_summary(summary):
    require_reuse(summary)
    families = summary.get("species_tree_families")
    if isinstance(families, bool) or not isinstance(families, int) or families <= 0:
        raise ValueError("Own inferred species tree has no valid native marker-family count")


def resolve(plan, output):
    parent_path = Path(plan["parent_results"])
    parent_identity = record(parent_path)
    parent = json.loads(parent_path.read_text())
    require_terminal_parent(parent, plan)
    accounting = subprocess.check_output(["sacct", "-j", plan["parent_job"], "--parsable2",
        "--format=JobIDRaw,State,ExitCode,Elapsed,NodeList,AllocCPUS"], text=True)
    scheduler = require_completed_job(accounting, plan["parent_job"])
    if scheduler["NodeList"] != "bizon" or scheduler["AllocCPUS"] != "16":
        raise ValueError("Wrong parent host/resources")
    check(plan["parent_plan"])
    previous = json.loads(Path(plan["parent_plan"]["path"]).read_text())
    snapshot = output / "parent_results_snapshot.json"
    write_json(snapshot, parent)
    checked = [parent_identity, record(snapshot), plan["parent_plan"]]
    readbacks, arms = [], []
    for row, old_case, new_case in zip(parent["arms"], previous["cases"], plan["cases"]):
        if row["label"] != old_case["label"] or row["label"] != new_case["label"]:
            raise ValueError("Parent/current arm inventories differ")
        for item in (row["seed"], row["candidate"]["seed_partition"], row["candidate"]["candidate_partition"],
                     row["candidate"]["membership_constraints"]):
            check(item)
            checked.append(item)
        for flag in ("--candidate-clusters", "--membership-constraints"):
            if argument(old_case["command"], flag) != argument(new_case["command"], flag):
                raise ValueError("Own-tree inference changed parent candidates or constraints")
        if argument(new_case["command"], "--checkpoint-source") != argument(old_case["command"], "--output-directory"):
            raise ValueError("Wrong matching raw-tree checkpoint source")
        native = replay_native(row, old_case)
        if native["genes"] != 984137 or native["species"] != 78:
            raise ValueError("Wrong corrected QfO parent universe")
        scores, details, files = readback_scores(Path(plan["repo"]), "QfO", row, native, Path(old_case["output"]),
                                                "ohmm_hillclimb_r1_" + old_case["scoring_id"])
        if row["label"] == "control":
            require_control_scores(scores, plan["expected_control_scores"])
        phylo = Path(native["phylogeny"])
        checked.extend([*native["checked"], *files])
        for name in ("alignments", "checkpoints", "gene_trees"):
            paths = [p for p in sorted((phylo / name).rglob("*")) if p.is_file()]
            if not paths:
                raise ValueError("Missing completed matching raw-tree cache")
            checked.extend(record(path) for path in paths)
        readbacks.append({"label": row["label"], "scores": scores, "details": details,
                          "summary": native["summary"], "status": "independent_parent_readback_verified"})
        arms.append({"label": row["label"], "status": "candidates_verified_unscored", "scores": {},
                     "seed": row["seed"], "candidate": row["candidate"]})
    checked = list({item["path"]: item for item in checked}.values())
    for item in checked:
        check(item)
    resolution = {"status": "terminal_parent_native_inputs_resolved_and_verified", "parent": parent_identity,
                  "snapshot": record(snapshot), "scheduler": scheduler, "readbacks": readbacks,
                  "checked": checked, "scope": "Input admission, not new own-inferred scores or default promotion."}
    path = output / "resolved_parent_inputs.json"
    write_json(path, resolution)
    return {"arms": arms, "checked": [record(path), *checked], "receipt": record(path), "readbacks": readbacks}
