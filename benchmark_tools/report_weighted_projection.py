"""Separately frozen six-arm readback with the registered constituent gates."""

import argparse
import hashlib
import json
from pathlib import Path
import re
import subprocess
import sys

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from benchmark_tools import report_profile_evidence_pilot as base
from benchmark_tools import report_projection_recovery as recovery
from benchmark_tools import weighted_projection_experiment as experiment
from benchmark_tools.prepare_ob_candidate_neighborhood import check, record
from benchmark_tools.run_full_search_hillclimb import runtime
from benchmark_tools.run_membership_hillclimb import write_json
from benchmark_tools.run_profile_weighting_pilot import advancement
from benchmark_tools.run_three_kingdoms_membership_supplement import sources


LEDGER = base.LEDGER
ENTRY = "weighted_projection_pilot"


def launch(repo, revision):
    if not isinstance(revision, str) or re.fullmatch(r"[0-9a-f]{40}", revision) is None:
        raise ValueError("Require an immutable full weighted projection launch revision")
    if subprocess.check_output(["git", "cat-file", "-t", revision], cwd=repo, text=True).strip() != "commit":
        raise ValueError("Weighted projection launch is not a commit")
    blob = subprocess.check_output(["git", "show", revision + ":" + LEDGER], cwd=repo)
    entry = json.loads(blob)[ENTRY]
    job = entry["inference"]
    reader = entry["independent_reader"]
    if (tuple(entry["arms"]) != experiment.LABELS or job.get("submitted") is not True
            or not str(job.get("job_id", "")).isdigit()
            or (job.get("planned_cpus"), job.get("planned_memory"), job.get("planned_wall_limit"))
            != (16, "128G", "1-00:00:00")):
        raise ValueError("Changed weighted projection scientific panel/resources")
    for revision_key, path_key, item in (("runner_commit", "scientific_source", entry),
                                         ("runner_commit", "scientific_reader_source", reader)):
        path = Path(item[path_key])
        if (re.fullmatch(r"[0-9a-f]{40}", item.get(revision_key, "")) is None
                or path.is_absolute() or not (repo / path).resolve().is_relative_to(repo.resolve())):
            raise ValueError("Require bound frozen scientific and separate reader sources")
    pin = {k: entry["execution_plan"][k] for k in ("path", "bytes", "sha256")}
    pin["path"] = str((repo / pin["path"]).resolve())
    return entry, pin, {"revision": revision, "path": LEDGER, "entry_key": ENTRY,
                        "git_blob_sha256": hashlib.sha256(blob).hexdigest()}


def verify_panel(cells, rows):
    if tuple(c["label"] for c in cells) != experiment.LABELS or tuple(r["label"] for r in rows) != experiment.LABELS:
        raise ValueError("Require every ordered arm of the registered six-arm panel")
    controls, gates = {}, {}
    for cell, row in zip(cells, rows):
        label = cell["label"]
        arm = experiment.configuration(label)
        activation = cell["activation"]
        if (cell.get("status") != "independent_content_readback_verified" or cell.get("dataset") != "OrthoBench"
                or cell.get("arm_parameters") != arm or row.get("arm_parameters") != arm
                or activation.get("arm") != label or activation.get("role") != arm["role"]
                or activation.get("mode") != arm["profile_sequence_weighting"]
                or activation.get("projection_policy") != arm["profile_projection_policy"]
                or activation.get("profiles_built", 0) <= 0):
            raise ValueError("Incomplete or mislabeled weighted projection activation")
        if arm["role"] != "candidate":
            if ("advancement" in row or activation.get("joint_profiles") != 0
                    or experiment.point(cell) != experiment.specification()["control_points"][label]
                    or row.get("graph_control_matches") is not True):
                raise ValueError("Constituent reference changed or became an advancing candidate")
            controls[label] = experiment.point(cell)
        else:
            if activation.get("joint_profiles") != activation["profiles_built"]:
                raise ValueError("Joint candidate has no complete actual construction evidence")
            gate = experiment.matched_gate(label, experiment.point(cell), controls, advancement)
            if row.get("advancement") != gate:
                raise ValueError("Candidate gates differ from the original three-control rule")
            cell["advancement"] = gate
            gates[label] = gate
    return gates


def uncertainty(cells):
    spec = experiment.specification()["uncertainty"]
    if (tuple(c["label"] for c in cells) != experiment.LABELS
            or spec["replicates"] != 100000 or spec["seed"] != 20261019
            or spec["alpha"] != .05 or spec["metrics"] != list(base.METRICS)):
        raise ValueError("Changed weighted projection uncertainty plan or incomplete panel")
    records = {}
    for cell in cells:
        details = cell["details"]["family_details"]
        names, _, counts = base.weighted_records(details["refog_records"])
        if len(names) != 70 or not np.allclose(base.statistics(counts.sum(axis=0)),
                [details[k] for k in base.METRICS], atol=1e-8, rtol=0):
            raise ValueError("Weighted family statistics do not reproduce official scores")
        records[cell["label"]] = details["refog_records"]
    expected = [[label, control] for label in experiment.LABELS[4:]
                for control in experiment.configuration(label)["matched_controls"]]
    if spec["contrasts"] != expected:
        raise ValueError("Changed six prespecified paired contrasts")
    contrasts = []
    for candidate, control in expected:
        result = base.paired_bootstrap({k: records[k] for k in (control, candidate)}, control,
            replicates=100000, seed=20261019, alpha=.05, multiplicity_endpoints=18)
        contrasts.append({"candidate": candidate, "control": control, "result": result})
    return {"contrasts": contrasts, "multiplicity_endpoints": 18,
            "scope": "Prespecified within-phase paired70-family intervals, in percentage points. Development-exposed, not campaign-wide correction, independent confirmation, equivalence or a point-gate amendment."}


def reader_inventory(driver):
    revision = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=driver, text=True).strip()
    pins = sources(driver) + [record(driver / "benchmark_tools" / p) for p in (
        experiment.PROTOCOL, base.PROTOCOL, "report_weighted_projection.slurm")]
    for pin in pins:
        relative = Path(pin["path"]).relative_to(driver)
        if relative.suffix != ".so" and hashlib.sha256(subprocess.check_output(
                ["git", "show", revision + ":" + str(relative)], cwd=driver)).hexdigest() != pin["sha256"]:
            raise ValueError("Weighted projection reader differs from committed source")
    return revision, pins


def run(repo, revision, output):
    if output.exists():
        raise FileExistsError("Preserve every weighted projection readback attempt")
    entry, plan_pin, binding = launch(repo, revision)
    check(plan_pin)
    plan = json.loads(Path(plan_pin["path"]).read_text())
    root, source = Path(plan["output"]), Path(plan["source_root"])
    if (Path(plan_pin["path"]).parent != root or Path(plan["repo"]) != repo
            or source != (repo / entry["scientific_source"]).resolve()
            or plan["source_revision"] != entry["runner_commit"]
            or plan.get("arm_parameters") != experiment.specification()["arms"]):
        raise ValueError("Weighted projection plan differs from immutable launch")
    producer = record(root / "results.json")
    report = json.loads(Path(producer["path"]).read_text())
    if report["execution_plan"] != plan_pin:
        raise ValueError("Producer uses a different weighted projection plan")
    base.require_completed(report, plan, experiment.LABELS, "weighted_projection", entry["inference"]["job_id"])
    accounting = base.terminal(entry["inference"]["job_id"])
    if subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=source, text=True).strip() != plan["source_revision"]:
        raise ValueError("Frozen weighted projection scientific source moved")
    if sources(source) + [record(source / "benchmark_tools/run_profile_weighting_pilot.slurm")] != plan["sources"]:
        raise ValueError("Weighted projection scientific source inventory changed")
    controls, constituent_pins = experiment.control_references(repo, source)
    if controls != plan.get("constituent_controls") or runtime() != plan["runtime"]:
        raise ValueError("Original constituent references or runtime changed")
    checked = [plan_pin, producer, *plan["inputs"], *plan["sources"], *plan["runtime"]["files"], *constituent_pins]
    for pin in checked: check(pin)
    names = recovery.ownership(plan, report)
    driver = Path(__file__).resolve().parents[1]
    reader_revision, reader_files = reader_inventory(driver)
    if (reader_revision != entry["independent_reader"]["runner_commit"]
            or driver != (repo / entry["independent_reader"]["scientific_reader_source"]).resolve()
            or driver == source):
        raise ValueError("Actual separately frozen reader differs from launch binding")
    context, context_pins = base.retained_context(repo)
    checked += [*reader_files, *context_pins]
    output.mkdir(parents=True)
    result = {"status": "running", "experiment": "weighted_projection", "launch_ledger": binding,
        "execution_plan": plan_pin, "producer": producer, "scheduler": accounting,
        "source_revision": plan["source_revision"], "reporter_revision": reader_revision,
        "reporter_files": reader_files, "protocol": record(driver / "benchmark_tools" / experiment.PROTOCOL),
        "context": context, "cells": [], "accuracy_result_admitted": False,
        "defaults_promoted": False, "goal_completed": False,
        "scope": "Complete-proteome controlled development pilot with shared-helper separate source/native/official-score readback, not independent biological confirmation. Seven endpoint slots per configuration remain unrun until gated extensions; cached shared-host phase costs are not fresh isolated tool runtimes."}
    write_json(output / "summary.json", result)
    try:
        for row in report["arms"]:
            cell = base.admit_arm(repo, plan, row, "weighted_projection", names)
            result["cells"].append(cell)
            checked += cell["checked"]
            write_json(output / "summary.json", result)
        result["matched_gates"] = verify_panel(result["cells"], report["arms"])
        result["uncertainty"] = uncertainty(result["cells"])
        result["advancing_labels"] = [label for label, gate in result["matched_gates"].items() if gate["eligible"]]
        result["required_next_actions"] = (
            "Run unchanged matched full QfO/Three Kingdoms closure for every admitted passing combination."
            if result["advancing_labels"] else
            "Stop both fixed combinations without retuning/full-dataset extension; preserve the bounded negative interaction result.")
        base.write_table(output / "scores.tsv", result["cells"], context)
        result["score_table"] = record(output / "scores.tsv")
        checked.append(result["score_table"])
        for pin in checked: check(pin)
        if runtime() != plan["runtime"]:
            raise ValueError("Runtime changed during weighted projection readback")
        result["checked"] = list({pin["path"]: pin for pin in checked}.values())
        result.update(status="weighted_projection_independently_read_back", accuracy_result_admitted=True)
    except Exception as error:
        result.update(status="failed", error_type=type(error).__name__, error=str(error))
        raise
    finally:
        write_json(output / "summary.json", result)
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", required=True, type=Path)
    parser.add_argument("--launch-revision", required=True)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    print(run(args.repo.resolve(), args.launch_revision, args.output.resolve())["status"])
