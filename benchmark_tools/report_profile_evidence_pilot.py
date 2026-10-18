"""Post-run native/score readback and paired-family uncertainty for profile pilots."""

from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
from pathlib import Path
import subprocess
import sys

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from benchmark_tools.audit_historical_profile_ablation import read_partition
from benchmark_tools.bootstrap_orthobench import METRICS, paired_bootstrap, statistics, weighted_records
from benchmark_tools.metric_frontier import ENDPOINTS
from benchmark_tools.prepare_ob_candidate_neighborhood import check, record
from benchmark_tools.replay_phylogeny import load_membership_constraints
from benchmark_tools.report_hillclimb import endpoint_vector, readback_scores, replay_native
from benchmark_tools.run_full_search_hillclimb import runtime
from benchmark_tools.run_membership_hillclimb import compare_golden, validate_partition, write_json
from benchmark_tools.run_profile_weighting_pilot import advancement, arm_commands, experiment_settings, normalized, validate_activation
from benchmark_tools.run_three_kingdoms_membership_supplement import candidate_control_matches, sources
from benchmark_tools.summarize_matched_resources import verbose_time
from orthohmm.accuracy import load_accuracy_checkpoint


LEDGER = "benchmark_tools/results/hmm_evidence_campaign_status_20261018.json"
LAUNCH_REVISION = "79e69e6c23abb9dee732dd613bf98e5d88e4deea"
ADMISSION = "benchmark_tools/results/hillclimb_comparison_readback_v27_20261016/summary.json"
PROTOCOL = "hmm_profile_evidence_readback_20261018.json"
RETAINED_KEYS = ("orthohmm_high_sensitivity", "orthohmm_phylogeny_satellite_v2", "orthofinder_3_1_5_full",
                 "orthofinder_3_1_5_sequence_only", "sonicparanoid_2_0_9", "proteinortho_6_3_6",
                 "fastoma_0_3_5", "orthomcl_1_4", "canonical_collection:canonical")


def terminal(job_id):
    command = ["sacct", "-j", str(job_id), "-X", "-P", "-o", "JobIDRaw,State,ExitCode,AllocCPUS,ReqMem,NodeList,Elapsed"]
    rows = list(csv.DictReader(io.StringIO(subprocess.check_output(command, text=True)), delimiter="|"))
    expected = {"JobIDRaw": str(job_id), "State": "COMPLETED", "ExitCode": "0:0",
                "AllocCPUS": "16", "ReqMem": "128G", "NodeList": "bizon"}
    if len(rows) != 1 or any(rows[0].get(key) != value for key, value in expected.items()):
        raise ValueError("Profile pilot is not uniquely terminal successful on planned resources")
    return {"command": command, "accounting": rows[0]}


def launch(repo, experiment):
    payload = subprocess.check_output(["git", "show", LAUNCH_REVISION + ":" + LEDGER], cwd=repo)
    ledger = json.loads(payload)
    key = "weighting_orthobench_pilot" if experiment == "weighting" else "alignment_guard_pilot"
    entry = ledger[key]
    plan_pin = dict(entry["execution_plan"])
    plan_pin["path"] = str((repo / plan_pin["path"]).resolve())
    return entry, plan_pin, {"revision": LAUNCH_REVISION, "path": LEDGER,
                            "git_blob_sha256": hashlib.sha256(payload).hexdigest()}


def require_completed(report, plan, arms, experiment, job_id):
    if (report.get("status") != "finished_all_arms_verified" or report.get("job_id") != str(job_id)
            or report.get("experiment", "weighting") != experiment or plan.get("experiment", "weighting") != experiment
            or tuple(plan.get("arms", [])) != arms or tuple(row.get("label") for row in report.get("arms", [])) != arms
            or report.get("source_revision") != plan["source_revision"] or plan.get("cpu") != 16
            or any(row.get("status") != "independent_content_readback_verified" or row.get("dataset") != "OrthoBench"
                   or "error" in row for row in report["arms"])
            or report.get("defaults_promoted") is not False or report.get("goal_completed") is not False):
        raise ValueError("Require the complete unchanged profile pilot, not partial scores")


def execution_receipt(execution, command, folder, stage):
    expected = {"log": record(folder / (stage + ".log")), "time": record(folder / (stage + ".time.txt"))}
    resources = verbose_time(Path(expected["time"]["path"]).read_text())
    if (execution.get("returncode") != 0 or execution.get("command") != command
            or any(execution.get(key) != value for key, value in expected.items())
            or execution.get("resources") != resources or resources["exit_status"] != 0):
        raise ValueError("Profile execution command, resource record or raw log differs")
    return list(expected.values())


def admit_arm(repo, plan, row, experiment, names):
    source, root = Path(plan["source_root"]), Path(plan["output"])
    folder, mode = root / row["label"], row["label"]
    graph_command, native_command = arm_commands(source, plan, folder, mode)
    preflight = folder / "preflight.json"
    if json.loads(preflight.read_text()) != {"graph": graph_command, "native": native_command}:
        raise ValueError("Arm preflight differs from fixed contrast")
    checked = [record(preflight)]
    for stage, command in (("graph", graph_command), ("native", native_command)):
        checked.extend(execution_receipt(row[stage + "_execution"], command, folder, stage))
    graph_path = folder / "graph_metrics.json"
    if row["graph_report"] != record(graph_path):
        raise ValueError("Graph report differs from producer pin")
    graph = json.loads(graph_path.read_text())
    activated = validate_activation(graph, mode, experiment)
    if activated != row["activation"]:
        raise ValueError("Reported profile activation differs from evidence")
    checked.append(record(graph_path))
    for stage in graph["stages"]:
        check(stage["output"])
        read_partition(Path(stage["output"]["path"]), names)
        checked.append(stage["output"])
    candidate, constraints = (row["candidate"][key] for key in ("candidate_partition", "membership_constraints"))
    expected = folder / "candidate/orthohmm_working_res"
    if (Path(candidate["path"]) != expected / "orthohmm_edges_clustered.txt"
            or Path(constraints["path"]) != expected / "phylogeny_candidate_merges.json"):
        raise ValueError("Candidate artifacts are outside the planned arm")
    for pin in (candidate, constraints, row["candidate"]["seed_partition"]):
        check(pin)
        checked.append(pin)
    read_partition(Path(candidate["path"]), names)
    load_membership_constraints(Path(constraints["path"]), Path(candidate["path"]))
    prediction = folder / "native/orthohmm_working_res/orthohmm_edges_clustered.txt"
    if validate_partition(Path(candidate["path"]), prediction) != row["coverage"]:
        raise ValueError("Reported complete native coverage differs")
    native = replay_native(row, {"command": native_command})
    if ((native["genes"], native["species"]) != (251378, 12) or native["summary"] != row["native_summary"]):
        raise ValueError("Native universe or summary differs")
    scores, details, pins = readback_scores(repo, "OrthoBench", row, native, folder, "profile_evidence_readback")
    if details != row["details"] or endpoint_vector("OrthoBench", scores) != row["metric_vector"]:
        raise ValueError("Per-family/aggregate score readback differs")
    arms, _ = experiment_settings(experiment)
    if mode == arms[0]:
        control = plan["case"]
        if (candidate_control_matches(Path(candidate["path"]), Path(constraints["path"]), control["candidate"], control["constraints"])
                != row["candidate_control_matches"]
                or compare_golden(Path(native["phylogeny"]), control["golden"]) != row["golden_matches"]
                or normalized(row) != control["expected"]):
            raise ValueError("Unchanged off control differs from canonical control")
    checked += [*native["checked"], *pins]
    for pin in checked:
        check(pin)
    return {"dataset": "OrthoBench", "label": mode, "status": "independent_content_readback_verified",
        "scores": endpoint_vector("OrthoBench", scores), "details": details,
        "coverage": {key: native[key] for key in ("genes", "species", "groups", "pairs")},
        "products": native["products"], "activation": activated, "source_revision": plan["source_revision"],
        "candidate": row["candidate"], "checked": list({pin["path"]: pin for pin in checked}.values()),
        "resources": {key: row[key] for key in ("graph_execution", "candidate_resources", "native_execution", "total_incremental_wall_s")}}


def uncertainty(cells):
    protocol = json.loads(Path(__file__).with_name(PROTOCOL).read_text())["uncertainty"]
    if any(protocol[key] != value for key, value in {
        "replicates": 100000, "seed": 20261018, "alpha": .05, "families": 70,
        "metrics": list(METRICS), "multiplicity": "three_metrics_per_noncontrol_arm"}.items()):
        raise ValueError("Changed post hoc uncertainty plan")
    records = {}
    for cell in cells:
        score = cell["details"]["family_details"]
        names, _, counts = weighted_records(score["refog_records"])
        if len(names) != 70 or not np.allclose(statistics(counts.sum(axis=0)), [score[key] for key in METRICS], atol=1e-8, rtol=0):
            raise ValueError("Family sufficient statistics do not reproduce full OrthoBench scores")
        records[cell["label"]] = score["refog_records"]
    options = {"replicates": 100000, "seed": 20261018, "alpha": .05, "multiplicity_endpoints": 3 * (len(cells) - 1)}
    result = paired_bootstrap(records, cells[0]["label"], **options)
    if result != paired_bootstrap(records, cells[0]["label"], **options):
        raise ValueError("Paired-family uncertainty does not reproduce")
    result["scope"] = "Post hoc descriptive uncertainty after observed development point results. Does not change the original point gates or provide independent confirmation."
    return result


def retained_context(repo):
    path = repo / ADMISSION
    committed = subprocess.check_output(["git", "show", LAUNCH_REVISION + ":" + ADMISSION], cwd=repo)
    if hashlib.sha256(committed).hexdigest() != record(path)["sha256"]:
        raise ValueError("Retained comparison snapshot differs from committed admission")
    data = json.loads(path.read_text())
    check(data["score_table"])
    rows = []
    for key in RETAINED_KEYS:
        matches = [row for row in data["methods"] if row["key"] == key]
        if len(matches) != 1 or set(matches[0]["scores"]) != set(ENDPOINTS):
            raise ValueError("Require complete unique retained comparator vectors")
        method = matches[0]
        for endpoint, score in method["scores"].items():
            if type(score) not in (int, float) or not np.isfinite(score) or not 0 <= score <= 1:
                raise ValueError("Invalid retained comparator score: " + endpoint)
        populations = [{name: cell[name] for name in ("dataset", "status", "retained_details", "details", "coverage", "semantics", "limitations")
                        if name in cell} for cell in method["cells"]
                       if cell["status"] in ("retained_registered_result", "independent_content_readback_verified")]
        if len(populations) != 3 or {cell["dataset"] for cell in populations} != {"OrthoBench", "QfO", "ThreeKingdoms"}:
            raise ValueError("Retained comparator population details are incomplete or ambiguous")
        rows.append({"key": key, "label": method["label"], "scores": method["scores"],
                     "population_details": populations, "status": "retained_context_not_reaudited_here"})
    return rows, [record(path), data["score_table"]]


def write_table(path, cells, context):
    with path.open("x") as stream:
        writer = csv.writer(stream, delimiter="\t", lineterminator="\n")
        writer.writerow(["method", "evidence", *ENDPOINTS])
        for row in context:
            writer.writerow([row["label"], row["status"], *[row["scores"][key] for key in ENDPOINTS]])
        for cell in cells:
            writer.writerow(["OrthoHMM profile " + cell["label"] + " (new pilot)", cell["status"],
                             *[cell["scores"].get(key, "NA") for key in ENDPOINTS]])


def run(repo, experiment, output):
    if output.exists():
        raise FileExistsError("Preserve every readback attempt")
    arms, _ = experiment_settings(experiment)
    entry, plan_pin, ledger_pin = launch(repo, experiment)
    check(plan_pin)
    plan = json.loads(Path(plan_pin["path"]).read_text())
    root = Path(plan["output"])
    producer = record(root / "results.json")
    report = json.loads(Path(producer["path"]).read_text())
    check(producer)
    if (report["execution_plan"] != plan_pin or Path(plan_pin["path"]).parent != root
            or Path(plan["repo"]) != repo or plan["source_revision"] != entry["runner_commit"]):
        raise ValueError("Pilot differs from immutable launch ledger")
    require_completed(report, plan, arms, experiment, entry["inference"]["job_id"])
    accounting = terminal(entry["inference"]["job_id"])
    source = Path(plan["source_root"])
    if subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=source, text=True).strip() != plan["source_revision"]:
        raise ValueError("Frozen scientific checkout revision differs from launch")
    expected_sources = sources(source) + [record(source / "benchmark_tools/run_profile_weighting_pilot.slurm")]
    if expected_sources != plan["sources"] or runtime() != plan["runtime"]:
        raise ValueError("Frozen scientific source inventory or runtime changed")
    checked = [plan_pin, producer, *plan["inputs"], *plan["sources"], *plan["runtime"]["files"]]
    for pin in checked:
        check(pin)
    names, _, _, _, _ = load_accuracy_checkpoint(Path(plan["case"]["checkpoint"]), verify=True)
    if len(names) != 251378 or len(set(names)) != len(names):
        raise ValueError("Incorrect complete OrthoBench gene universe")
    driver = Path(__file__).resolve().parents[1]
    reporter_revision = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=driver, text=True).strip()
    reporter_files = sources(driver) + [record(Path(__file__).with_name(PROTOCOL))]
    for pin in reporter_files:
        relative = Path(pin["path"]).relative_to(driver)
        if relative.suffix != ".so":
            committed = subprocess.check_output(["git", "show", reporter_revision + ":" + str(relative)], cwd=driver)
            if hashlib.sha256(committed).hexdigest() != pin["sha256"]:
                raise ValueError("Readback source differs from committed bytes: " + str(relative))
    context, context_pins = retained_context(repo)
    checked += [*reporter_files, *context_pins]
    output.mkdir(parents=True)
    result = {"status": "running", "experiment": experiment, "producer": producer, "execution_plan": plan_pin,
        "source_revision": plan["source_revision"], "launch_ledger": ledger_pin, "scheduler": accounting,
        "reporter_revision": reporter_revision, "protocol": record(Path(__file__).with_name(PROTOCOL)),
        "reporter_files": reporter_files, "cells": [], "context": context, "accuracy_result_admitted": False,
        "defaults_promoted": False, "goal_completed": False,
        "scope": "Separate post-run content/score readback with shared validated native/scoring helpers, not an independent scoring implementation, independent biological panel or full-vector improvement. Timings are shared-host cached downstream replay, not isolated fresh inference. Functional EC/GO/FAS means are not F1 and evaluated populations can differ."}
    write_json(output / "summary.json", result)
    try:
        for row in report["arms"]:
            cell = admit_arm(repo, plan, row, experiment, set(names))
            if row["label"] != arms[0]:
                control = {"f1": result["cells"][0]["scores"][ENDPOINTS[0]],
                           "precision": result["cells"][0]["details"]["precision"], "recall": result["cells"][0]["details"]["recall"]}
                actual = {"f1": cell["scores"][ENDPOINTS[0]], "precision": cell["details"]["precision"], "recall": cell["details"]["recall"]}
                cell["advancement"] = advancement(control, actual)
                if cell["advancement"] != row["advancement"]:
                    raise ValueError("Independent original advancement gate differs")
            result["cells"].append(cell)
            checked += cell["checked"]
            write_json(output / "summary.json", result)
        result["uncertainty"] = uncertainty(result["cells"])
        write_table(output / "scores.tsv", result["cells"], context)
        result["score_table"] = record(output / "scores.tsv")
        for pin in checked:
            check(pin)
        if runtime() != plan["runtime"]:
            raise ValueError("Runtime changed during independent readback")
        result["checked"] = list({pin["path"]: pin for pin in checked}.values())
        result["advancing_labels"] = [cell["label"] for cell in result["cells"][1:] if cell["advancement"]["eligible"]]
        result.update(status="profile_pilot_independently_read_back", accuracy_result_admitted=True)
    except Exception as error:
        result.update(status="failed", error_type=type(error).__name__, error=str(error))
        raise
    finally:
        write_json(output / "summary.json", result)
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", type=Path, required=True)
    parser.add_argument("--experiment", choices=("weighting", "alignment"), required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    print(run(args.repo.resolve(), args.experiment, args.output.resolve())["status"])
