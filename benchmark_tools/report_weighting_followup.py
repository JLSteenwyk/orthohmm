"""Read complete frozen weighting follow-ups and assemble eight-endpoint vectors."""

import argparse
import csv
import hashlib
import json
from pathlib import Path
import subprocess
import sys
from types import SimpleNamespace

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from Bio import SeqIO
from benchmark_tools import run_profile_weighting_followup as followup
from benchmark_tools.audit_historical_profile_ablation import read_partition
from benchmark_tools.metric_frontier import ENDPOINTS, metric_frontier
from benchmark_tools.prepare_ob_candidate_neighborhood import check, record
from benchmark_tools.replay_phylogeny import load_membership_constraints
from benchmark_tools.report_hillclimb import endpoint_vector, readback_scores, replay_native
from benchmark_tools.report_profile_evidence_pilot import execution_receipt, retained_context, terminal
from benchmark_tools.run_canonical_downstream import commands
from benchmark_tools.run_full_search_hillclimb import runtime
from benchmark_tools.run_membership_hillclimb import compare_golden, validate_partition, write_json
from benchmark_tools.run_profile_weighting_pilot import advancement, validate_activation
from benchmark_tools.run_three_kingdoms_membership_supplement import candidate_control_matches, sources, validate_ownership
from orthohmm.accuracy import load_accuracy_checkpoint


LAUNCH = "93a4a96d7782889703cf8325871d6f281bf0f0ee"
LEDGER = "benchmark_tools/results/hmm_evidence_campaign_status_20261018.json"
PROTOCOL = "hmm_profile_weighting_postrun_20261018.json"
DATASETS = ("QfO", "ThreeKingdoms")


def launch(repo):
    blob = subprocess.check_output(["git", "show", LAUNCH + ":" + LEDGER], cwd=repo)
    entry = json.loads(blob)["profile_weighting_followup_execution"]
    if (tuple(row["dataset"] for row in entry["datasets"]) != DATASETS
            or entry["scientific_revision"] != followup.REVISION):
        raise ValueError("Wrong immutable weighting follow-up launch")
    return entry["datasets"], {"revision": LAUNCH, "path": LEDGER, "git_blob_sha256": hashlib.sha256(blob).hexdigest()}


def require_completed(plan, report, entry):
    followup.validate_plan(plan)
    if (plan["driver_revision"] != "6e684fb2d050563d03578d636e3de9325ecec73c"
            or report.get("status") != "all_weighting_followup_arms_scored"
            or report.get("job_id") != entry["inference"]["job_id"]
            or report.get("dataset") != entry["dataset"] or plan["dataset"] != entry["dataset"]
            or report.get("source_revision") != followup.REVISION or report.get("driver_revision") != plan["driver_revision"]
            or [row.get("label") for row in report.get("arms", [])] != list(followup.ARMS)
            or any(row.get("status") != "native_and_score_independently_verified"
                   or row.get("dataset") != plan["dataset"] or "error" in row for row in report["arms"])
            or report.get("defaults_promoted") is not False or report.get("goal_completed") is not False):
        raise ValueError("Require complete successful ordered weighting follow-up, not partial scores")


def ownership(plan, report):
    names, species, q, t, scores = load_accuracy_checkpoint(Path(plan["case"]["checkpoint"]), verify=True)
    del q, t, scores
    owners = {}
    for pin in plan["case"]["case"]["fastas"]:
        check(pin)
        for sequence in SeqIO.parse(pin["path"], "fasta"):
            if sequence.id in owners or not sequence.seq:
                raise ValueError("Duplicate or empty follow-up input sequence")
            owners[sequence.id] = Path(pin["path"]).name
    actual = json.loads(json.dumps(validate_ownership(names, species, owners)))
    if actual != report["ownership"] or (actual["genes"], actual["species"]) != followup.UNIVERSES[plan["dataset"]]:
        raise ValueError("Complete follow-up checkpoint/FASTA ownership differs")
    return names


def admit_arm(repo, plan, row, names):
    dataset, mode = plan["dataset"], row["label"]
    if mode not in followup.ARMS: raise ValueError("Unknown weighting follow-up arm")
    case, source = plan["case"], Path(plan["scientific_source"])
    folder = Path(plan["output"]) / mode
    graph_command, native_command = commands(source, plan, case["case"], folder,
        Path(case["checkpoint"]), case["checkpoint_sha"])
    if mode != "none": graph_command += ["--profile-sequence-weighting", mode]
    preflight = folder / "preflight.json"
    if json.loads(preflight.read_text()) != {"graph": graph_command, "native": native_command}:
        raise ValueError("Follow-up preflight differs from fixed scientific contrast")
    checked = [record(preflight), *row["checked"]]
    for stage, command in (("graph", graph_command), ("native", native_command)):
        checked += execution_receipt(row[stage + "_execution"], command, folder, stage)
    graph_path = folder / "graph_metrics.json"
    if record(graph_path) != row["graph_report"]: raise ValueError("Graph receipt differs")
    graph = json.loads(graph_path.read_text())
    activation = validate_activation(graph, mode)
    if activation != row["activation"]: raise ValueError("Weighted profile activation differs")
    checked.append(record(graph_path))
    for stage in graph["stages"]:
        check(stage["output"])
        read_partition(Path(stage["output"]["path"]), set(names))
        checked.append(stage["output"])
    working = folder / "candidate/orthohmm_working_res"
    for key, filename in (("candidate_partition", "orthohmm_edges_clustered.txt"),
                          ("membership_constraints", "phylogeny_candidate_merges.json")):
        if row["candidate"][key] != record(working / filename):
            raise ValueError("Follow-up candidate or constraints differ from planned namespace")
        checked.append(row["candidate"][key])
    seed = row["candidate"].get("seed_partition")
    if seed is None or seed != record(folder / "graph/orthogroups_profiles_refined.txt"):
        raise ValueError("Follow-up candidate seed differs from final graph partition")
    checked.append(seed)
    candidate = working / "orthohmm_edges_clustered.txt"
    constraints = working / "phylogeny_candidate_merges.json"
    read_partition(candidate, set(names))
    load_membership_constraints(constraints, candidate)
    prediction = folder / "native/orthohmm_working_res/orthohmm_edges_clustered.txt"
    if validate_partition(candidate, prediction) != row["coverage"]:
        raise ValueError("Follow-up complete partition coverage differs")
    native = replay_native(row, {"command": native_command})
    if (native["genes"], native["species"]) != followup.UNIVERSES[dataset] or native["summary"] != row["native_summary"]:
        raise ValueError("Follow-up native universe or summary differs")
    scoring_folder, participant = folder, None
    reused = None
    if mode == "none":
        candidate_control_matches(candidate, constraints, case["candidate"]["candidate_partition"], case["candidate"]["membership_constraints"])
        if compare_golden(Path(native["phylogeny"]), case["golden"]) != row["golden_matches"]:
            raise ValueError("Canonical native control differs")
        if dataset == "QfO":
            scoring_folder = repo / "benchmarks/results/canonical_qfo_scoring_recovery_r1_v1_20261014"
            reused = record(scoring_folder / "results.json")
            recovery_plan = record(scoring_folder / "execution_plan.json")
            recovered = json.loads(Path(reused["path"]).read_text())
            if (row.get("control_scoring_reuse") != reused or reused not in plan["inputs"]
                    or recovered["execution_plan"] != recovery_plan or recovery_plan not in plan["inputs"]
                    or recovered["status"] != "canonical_qfo_scoring_recovered_independently_verified"
                    or any(row.get(key) != recovered.get(key) for key in ("scores", "conversion", "metric_files", "references"))):
                raise ValueError("QfO control assessment reuse is not exactly bound to canonical recovery")
            participant = "ohmm_hillclimb_r1_" + json.loads(Path(recovery_plan["path"]).read_text())["scoring_id"]
            checked += [reused, recovery_plan]
    elif dataset == "QfO":
        identifier = followup.scoring_id(Path(plan["output"]), mode)
        if plan["scoring_ids"][mode] != identifier: raise ValueError("Changed QfO scoring namespace")
        participant = "ohmm_hillclimb_r1_" + identifier
        if "control_scoring_reuse" in row: raise ValueError("Weighted QfO arm must use a fresh assessment")
    scores, details, pins = readback_scores(repo, dataset, row, native, scoring_folder, participant)
    vector = endpoint_vector(dataset, scores)
    if scores != row["scores"] or details != row["details"] or vector != row["metric_vector"]:
        raise ValueError("Follow-up scores, populations or endpoint vector do not reproduce")
    if mode == "none" and vector != case["control"]["scores"]:
        raise ValueError("Follow-up control score vector differs from original admitted canonical scores")
    checked += [*native["checked"], *pins]
    for pin in checked: check(pin)
    return {"label": mode, "dataset": dataset, "status": "independent_content_readback_verified",
        "source_revision": plan["scientific_revision"], "scores": vector, "details": details,
        "coverage": {key: native[key] for key in ("genes", "species", "groups", "pairs")},
        "products": native["products"], "activation": activation, "candidate": row["candidate"],
        "resources": {key: row[key] for key in ("graph_execution", "candidate_resources", "native_execution", "total_incremental_wall_s")},
        "control_assessment_reuse": reused, "checked": list({p["path"]: p for p in checked}.values())}


def dataset_deltas(cells):
    if [cell["label"] for cell in cells] != list(followup.ARMS) or len({cell["dataset"] for cell in cells}) != 1:
        raise ValueError("Require all ordered weighting arms for one dataset")
    control = cells[0]["scores"]
    if any(set(cell["scores"]) != set(control) for cell in cells):
        raise ValueError("Follow-up endpoint sets differ")
    return {cell["label"]: {key: value - control[key] for key, value in cell["scores"].items()} for cell in cells[1:]}


def assemble(pilot, cells):
    methods = []
    for original in pilot["cells"]:
        label = original["label"]
        selected = [original, *[c for c in cells if c["label"] == label]]
        expected = {"OrthoBench", *DATASETS} if label in followup.ARMS else {"OrthoBench"}
        if len(selected) != len(expected) or {cell["dataset"] for cell in selected} != expected:
            raise ValueError("Weighting vectors require exactly every declared dataset, without selected winners")
        vector = {key: value for cell in selected for key, value in cell["scores"].items()}
        metric_frontier({label: vector})
        methods.append({"key": "profile_weighting:" + label, "label": "OrthoHMM profile weighting " + label,
            "status": "complete_development_vector" if label in followup.ARMS else "stopped_after_orthobench_point_gate",
            "scores": vector, "cells": selected, "missing_endpoints": [key for key in ENDPOINTS if key not in vector]})
    control = methods[0]["scores"]
    for method in methods[1:]:
        if method["missing_endpoints"]: continue
        deltas = {key: method["scores"][key] - control[key] for key in ENDPOINTS}
        primary_gain = any(deltas[key] > 0 for key in ENDPOINTS[:3])
        losses = [key for key in ENDPOINTS if method["scores"][key] < control[key] - .005]
        method["full_vector_screen"] = {"eligible": primary_gain and not losses, "deltas": deltas,
            "primary_f1_improvement": primary_gain, "declines_exceeding_0_005": losses,
            "scope": "Original development point screen, not significance, independent generalizability or default promotion. Functional populations require separate interpretation."}
    return methods


def write_table(path, methods, context):
    with path.open("x") as stream:
        writer = csv.writer(stream, delimiter="\t", lineterminator="\n")
        writer.writerow(["method", "evidence", *ENDPOINTS])
        for method in [*context, *methods]:
            writer.writerow([method["label"], method["status"], *[method["scores"].get(key, "NA") for key in ENDPOINTS]])


def run(repo, output):
    if output.exists(): raise FileExistsError("Preserve every weighting readback namespace")
    entries, binding = launch(repo)
    api = SimpleNamespace(check=check, record=record, advancement=advancement)
    _, pilot, _, checked, _ = followup.pilot_inputs(api, repo)
    context, context_pins = retained_context(repo)
    checked += context_pins
    reader = Path(__file__).resolve().parents[1]
    protocol = json.loads((reader / "benchmark_tools" / PROTOCOL).read_text())
    if (protocol["launch_revision"] != LAUNCH or protocol["scientific_revision"] != followup.REVISION
            or protocol["arms"] != list(followup.ARMS) or protocol["datasets"] != list(DATASETS)
            or protocol["screen"]["maximum_endpoint_point_decline"] != .005
            or protocol["screen"]["require_complete_eight_endpoints"] is not True
            or protocol["screen"]["require_at_least_one_primary_f1_improvement"] is not True):
        raise ValueError("Changed fixed weighting post-run protocol")
    revision = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=reader, text=True).strip()
    reader_files = sources(reader) + [record(reader / "benchmark_tools" / PROTOCOL)]
    for pin in reader_files:
        relative = Path(pin["path"]).relative_to(reader)
        if relative.suffix != ".so" and hashlib.sha256(subprocess.check_output(
                ["git", "show", revision + ":" + str(relative)], cwd=reader)).hexdigest() != pin["sha256"]:
            raise ValueError("Weighting reader differs from committed source")
    checked += reader_files
    output.mkdir(parents=True)
    result = {"status": "running", "launch_ledger": binding, "reporter_revision": revision, "reporter_files": reader_files,
        "protocol": record(reader / "benchmark_tools" / PROTOCOL), "datasets": [], "cells": [], "context": context,
        "orthobench_uncertainty": pilot["uncertainty"], "accuracy_result_admitted": False,
        "defaults_promoted": False, "goal_completed": False,
        "scope": "Shared-algorithm post-run readback of development-exposed panels, not independent biological confirmation. Functional means are not orthology F1 and assessed populations may differ. FAS is unseeded; shared-host cached timings and max-process RSS are not isolated fresh-pipeline timing or aggregate memory. Unavailable defensible QfO matched-family uncertainty remains explicit, not a reason to fabricate pair-independent intervals."}
    write_json(output / "summary.json", result)
    try:
        for entry in entries:
            plan_pin = {key: entry["execution_plan"][key] for key in ("path", "bytes", "sha256")}
            plan_pin["path"] = str((repo / plan_pin["path"]).resolve())
            check(plan_pin)
            plan = json.loads(Path(plan_pin["path"]).read_text())
            root = Path(plan["output"])
            producer = record(root / "results.json")
            report = json.loads(Path(producer["path"]).read_text())
            require_completed(plan, report, entry)
            if report["execution_plan"] != plan_pin or Path(plan_pin["path"]).parent != root or Path(plan["repo"]) != repo:
                raise ValueError("Follow-up producer is not bound to its launched plan")
            accounting = terminal(entry["inference"]["job_id"])
            for key, field, inventory in (("scientific_source", "scientific_revision", plan["scientific_sources"]),
                                          ("driver_root", "driver_revision", plan["driver_sources"])):
                source = Path(plan[key])
                if subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=source, text=True).strip() != plan[field]:
                    raise ValueError("Frozen weighting scientific or driver checkout moved")
                expected = (sources(source) + [record(source / "benchmark_tools/run_profile_weighting_pilot.slurm")]
                    if key == "scientific_source" else [record(source / relative) for relative in (
                        "benchmark_tools/run_profile_weighting_followup.py", "benchmark_tools/run_profile_weighting_followup.slurm", followup.PROTOCOL)])
                if expected != inventory: raise ValueError("Frozen weighting source inventory differs")
            pins = [plan_pin, producer, *plan["inputs"], *plan["scientific_sources"], *plan["driver_sources"], *plan["runtime"]["files"]]
            for pin in pins: check(pin)
            if runtime() != plan["runtime"]: raise ValueError("Weighting runtime differs")
            names = ownership(plan, report)
            cells = [admit_arm(repo, plan, row, names) for row in report["arms"]]
            if dataset_deltas(cells) != report["deltas"]: raise ValueError("Original per-dataset endpoint deltas differ")
            result["datasets"].append({"dataset": plan["dataset"], "execution_plan": plan_pin, "producer": producer, "scheduler": accounting})
            result["cells"] += cells
            checked += [*pins, *[p for cell in cells for p in cell["checked"]]]
            write_json(output / "summary.json", result)
        result["methods"] = assemble(pilot, result["cells"])
        alignment_path = repo / "benchmark_tools/results/aligned_profile_guard_admission_20261018.json"
        blob = subprocess.check_output(["git", "show", LAUNCH + ":" + str(alignment_path.relative_to(repo))], cwd=repo)
        if hashlib.sha256(blob).hexdigest() != record(alignment_path)["sha256"]: raise ValueError("Stopped alignment admission changed")
        admission = json.loads(blob)
        alignment_pin = dict(admission["readback"])
        alignment_pin["path"] = str((repo / alignment_pin["path"]).resolve())
        check(alignment_pin)
        alignment = json.loads(Path(alignment_pin["path"]).read_text())
        if alignment["advancing_labels"] or alignment["accuracy_result_admitted"] is not True:
            raise ValueError("Alignment branch is not an admitted stopped branch")
        guard, = [cell for cell in alignment["cells"] if cell["label"] == "strict_span_guard"]
        result["methods"].append({"key": "alignment:strict_span_guard", "label": "OrthoHMM alignment guard",
            "status": "stopped_after_orthobench_point_gate", "scores": guard["scores"], "cells": [guard],
            "missing_endpoints": [key for key in ENDPOINTS if key not in guard["scores"]]})
        checked += [record(alignment_path), alignment_pin]
        write_table(output / "scores.tsv", result["methods"], context)
        result["score_table"] = record(output / "scores.tsv")
        result["frontier"] = metric_frontier({m["key"]: m["scores"] for m in [*context, *result["methods"]]})
        for pin in checked: check(pin)
        result["checked"] = list({p["path"]: p for p in checked}.values())
        result.update(status="weighting_full_vectors_independently_read_back", accuracy_result_admitted=True)
    except Exception as error:
        result.update(status="failed", error_type=type(error).__name__, error=str(error))
        raise
    finally:
        write_json(output / "summary.json", result)
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    print(run(args.repo.resolve(), args.output.resolve())["status"])
