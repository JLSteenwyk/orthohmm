"""Independently admit fixed projection transfers and retain all eight endpoints."""

import argparse
from copy import deepcopy
import csv
from datetime import datetime
import hashlib
import io
import json
from pathlib import Path
import re
import subprocess
import sys
from types import SimpleNamespace

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from benchmark_tools import plan_projection_followup as planner
from benchmark_tools import report_profile_evidence_pilot as pilot_reader
from benchmark_tools import report_projection_recovery as recovery
from benchmark_tools import report_weighting_followup as weighting
from benchmark_tools import run_projection_followup as followup
from benchmark_tools.audit_historical_profile_ablation import read_partition
from benchmark_tools.audit_qfo_fas_samples import read_sample, summarize
from benchmark_tools.compare_qfo_scored_pairs import read_scores
from benchmark_tools.metric_frontier import ENDPOINTS, metric_frontier
from benchmark_tools.native_pair_turnover import native_pair_turnover
from benchmark_tools.prepare_ob_candidate_neighborhood import check, record
from benchmark_tools.replay_phylogeny import load_membership_constraints
from benchmark_tools.report_hillclimb import endpoint_vector, readback_scores, replay_native
from benchmark_tools.report_hmm_evidence_campaign import unique_pins
from benchmark_tools.run_canonical_downstream import commands
from benchmark_tools.run_full_search_hillclimb import runtime
from benchmark_tools.run_membership_hillclimb import compare_golden, validate_partition, write_json
from benchmark_tools.run_profile_weighting_pilot import recovery_configuration, validate_activation
from benchmark_tools.run_three_kingdoms_membership_supplement import sources


PROTOCOL = "hmm_projection_followup_readback_20261019.json"
DRIVER = "80974fba96897fccc8ec2eea36ba08bb6830bd85"
DATASETS = ("QfO", "ThreeKingdoms")
LEDGER = "benchmark_tools/results/hmm_evidence_campaign_status_20261018.json"
DRIVER_FILES = ("benchmark_tools/run_projection_followup.py", "benchmark_tools/run_projection_followup.slurm",
                "benchmark_tools/run_profile_weighting_followup.py", followup.PROTOCOL)


def configuration(campaign_name="projection"):
    if campaign_name == "projection":
        return {"driver": DRIVER, "protocol": PROTOCOL, "entry": "projection_followup_execution"}
    if campaign_name == "stockholm_annotation":
        return {"driver": "abe393ded90143880378136bf60691581ed24061",
            "protocol": "hmm_stockholm_projection_followup_readback_20261019.json",
            "entry": "stockholm_projection_followup_execution"}
    raise ValueError("Unknown explicit projection readback campaign")


def relative_path(repo, value):
    path = Path(value)
    if path.is_absolute() or not (repo / path).resolve().is_relative_to(repo.resolve()):
        raise ValueError("Launch path must remain relative inside repository")
    return (repo / path).resolve()


def launch(repo, revision, *, campaign_name="projection"):
    if not isinstance(revision, str) or re.fullmatch(r"[0-9a-f]{40}", revision) is None:
        raise ValueError("Require explicit full immutable transfer launch revision")
    if subprocess.check_output(["git", "cat-file", "-t", revision], cwd=repo, text=True).strip() != "commit":
        raise ValueError("Transfer launch must be a commit")
    binding, transfer = configuration(campaign_name), followup.configuration(campaign_name)
    blob = subprocess.check_output(["git", "show", revision + ":" + LEDGER], cwd=repo)
    entry = json.loads(blob)[binding["entry"]]
    if (entry.get("runner_commit") != binding["driver"] or entry.get("scientific_revision") != transfer["science"]
            or entry.get("campaign", "projection") != campaign_name
            or tuple(row["dataset"] for row in entry.get("datasets", [])) != DATASETS):
        raise ValueError("Changed original transfer driver/science or dataset inventory")
    relative_path(repo, entry["source"])
    reader = entry["independent_reader"]
    if re.fullmatch(r"[0-9a-f]{40}", reader.get("runner_commit", "")) is None:
        raise ValueError("Require separately frozen transfer reader revision")
    relative_path(repo, reader["scientific_reader_source"])
    if relative_path(repo, entry["prepared_inputs"]["path"]) != repo / transfer["prepared"]:
        raise ValueError("Changed original transfer preparation namespace")
    jobs = []
    for row in entry["datasets"]:
        inference = row["inference"]
        if (inference.get("submitted") is not True or not str(inference.get("job_id", "")).isdigit()
                or inference.get("planned_cpus") != 16 or inference.get("planned_memory") != "128G"
                or inference.get("planned_wall_limit") != "1-00:00:00"):
            raise ValueError("Require original submitted transfer jobs and planned resources")
        relative_path(repo, row["execution_plan"]["path"])
        jobs.append(str(inference["job_id"]))
    if len(set(jobs)) != len(jobs): raise ValueError("Transfer datasets cannot share a scientific job")
    receipt = {"revision": revision, "path": LEDGER, "git_blob_sha256": hashlib.sha256(blob).hexdigest()}
    if campaign_name != "projection": receipt["entry_key"] = binding["entry"]
    return entry, receipt


def require_completed(plan, report, entry, chosen):
    followup.validate_plan(plan)
    campaign_name = plan.get("campaign", "projection")
    binding, transfer = configuration(campaign_name), followup.configuration(campaign_name)
    if (plan["driver_revision"] != binding["driver"] or plan["selection"] != chosen
            or report.get("status") != "all_projection_followup_arms_scored"
            or report.get("job_id") != entry["inference"]["job_id"]
            or report.get("source_revision") != transfer["science"] or report.get("driver_revision") != binding["driver"]
            or report.get("campaign", "projection") != campaign_name
            or report.get("dataset") != entry["dataset"] or plan["dataset"] != entry["dataset"]
            or report.get("selection") != chosen or [row["label"] for row in report.get("arms", [])] != chosen["execution_arms"]
            or any(row.get("status") != followup.SUCCESS or row.get("dataset") != plan["dataset"]
                   or row.get("role") != chosen["roles"][row["label"]] or "error" in row for row in report["arms"])
            or report.get("defaults_promoted") is not False or report.get("goal_completed") is not False):
        raise ValueError("Require complete unchanged successful transfer, not partial scores or relabeled controls")


def sequential_science(entries):
    identifiers = [str(entry["inference"]["job_id"]) for entry in entries]
    command = ["sacct", "-j", ",".join(identifiers), "-X", "-P", "-o", "JobIDRaw,State,ExitCode,Start,End"]
    rows = list(csv.DictReader(io.StringIO(subprocess.check_output(command, text=True)), delimiter="|"))
    by_id = {row["JobIDRaw"]: row for row in rows}
    if len(rows) != len(identifiers) or set(by_id) != set(identifiers):
        raise ValueError("Require unique accounting for every ordered transfer job")
    ordered = [by_id[identifier] for identifier in identifiers]
    prior_end = None
    for row in ordered:
        if (row["State"], row["ExitCode"]) != ("COMPLETED", "0:0"):
            raise ValueError("Sequential transfer accounting is not terminal successful")
        start, end = datetime.fromisoformat(row["Start"]), datetime.fromisoformat(row["End"])
        if end < start or prior_end is not None and start < prior_end:
            raise ValueError("QfO/Three Kingdoms inference did not run sequentially in original order")
        prior_end = end
    return {"command": command, "accounting": ordered,
        "scope": "Original QfO then Three Kingdoms scientific jobs do not overlap. Does not imply a quiet host or isolated full-pipeline timing."}


def arm_commands(plan, folder, label):
    case = plan["case"]
    graph, native = commands(Path(plan["scientific_source"]), plan, case["case"], folder, Path(case["checkpoint"]), case["checkpoint_sha"])
    policy, mode = recovery_configuration(label)
    if policy != "legacy": graph += ["--profile-projection-policy", policy]
    if mode != "off": graph += ["--selective-profile-mode", mode]
    return graph, native


def scoring_context(repo, plan, row, projection):
    folder = Path(plan["output"]) / row["label"]
    if ("control_scoring_reuse" in row and (row["label"], plan["dataset"]) != ("legacy_off", "QfO")
            or "measurement_scoring_reuse" in row and row["label"] != "projection_hmmer_score_only"):
        raise ValueError("Assessment reuse is attached to the wrong control")
    if row["label"] == "legacy_off" and plan["dataset"] == "QfO":
        folder = repo / "benchmarks/results/canonical_qfo_scoring_recovery_r1_v1_20261014"
        pin, plan_pin = record(folder / "results.json"), record(folder / "execution_plan.json")
        original = json.loads(Path(pin["path"]).read_text())
        if (row.get("control_scoring_reuse") != pin or pin not in plan["inputs"] or plan_pin not in plan["inputs"]
                or original.get("status") != "canonical_qfo_scoring_recovered_independently_verified"
                or original["execution_plan"] != plan_pin
                or any(row.get(key) != original.get(key) for key in ("scores", "conversion", "metric_files", "references"))):
            raise ValueError("Canonical QfO assessment reuse is not exactly bound")
        identifier = json.loads(Path(plan_pin["path"]).read_text())["scoring_id"]
        return folder, "ohmm_hillclimb_r1_" + identifier, [pin, plan_pin]
    if row["label"] == "projection_hmmer_score_only":
        if projection is None or projection["label"] != "projection_off":
            raise ValueError("Measurement readback requires admitted repaired control")
        folder = Path(plan["output"]) / "projection_off"
        receipt = record(folder / "assessment_receipt.json")
        if (row.get("measurement_scoring_reuse") != receipt or projection["assessment_receipt"] != receipt
                or any(row.get(key) != projection["assessment"].get(key) for key in ("scores", "conversion", "metric_files", "references"))):
            raise ValueError("Measurement assessment must exactly reuse repaired control")
        participant = "ohmm_hillclimb_r1_" + plan["scoring_ids"]["projection_off"] if plan["dataset"] == "QfO" else None
        return folder, participant, [receipt]
    if "control_scoring_reuse" in row or "measurement_scoring_reuse" in row:
        raise ValueError("Only exact controls may reuse assessments")
    participant = "ohmm_hillclimb_r1_" + plan["scoring_ids"][row["label"]] if plan["dataset"] == "QfO" else None
    return folder, participant, []


def functional_populations(folder, scores, details, admitted):
    result, checked = {}, []
    for metric in ("EC", "GO", "FAS"):
        directory = folder / "qfo_scoring/results" / metric
        paths = list(directory.glob("*raw.txt.gz"))
        if len(paths) != 1: raise ValueError("Require unique saved numeric " + metric + " population")
        raw = record(paths[0])
        aggregate = record(directory / (metric + ".json"))
        if aggregate not in admitted: raise ValueError("Functional aggregate absent from independent score admission")
        binding = "Saved raw population newly pinned/read back against the admitted aggregate; legacy aggregate admission did not itself pin raw numeric samples."
        if metric == "FAS":
            inline = json.loads(Path(aggregate["path"]).read_text())["datalink"]["inline_data"]
            participants = inline["challenge_participants"]
            if len(participants) != 1: raise ValueError("Ambiguous FAS participant")
            sample = summarize(read_sample(paths[0]), participants[0])
            if sample["reported_eligible_pairs"] != details[metric]["assessed_relations"] or abs(sample["mean"] - scores[metric]) > 1e-12:
                raise ValueError("FAS numeric sample differs from admitted score/population")
            result[metric] = {"raw": raw, "aggregate": aggregate, "binding": binding, **sample,
                "scope": "Saved unseeded numeric sample; eligible pool is not sample size. Native IID SEM is arithmetic only, not valid biological uncertainty with repeated proteins/families."}
        else:
            values = read_scores(paths[0], metric)
            mean = sum(values.values()) / (len(values) * 1000000)
            if len(values) != details[metric]["assessed_relations"] or abs(mean - scores[metric]) > 5.051e-7:
                raise ValueError("Rounded functional pair count/mean differs from admitted " + metric)
            result[metric] = {"raw": raw, "aggregate": aggregate, "binding": binding, "numeric_scored_pairs": len(values), "rounded_raw_mean": mean,
                "native_mean": scores[metric], "mean_tolerance": 5.051e-7,
                "scope": "Annotation-eligible scored pairs, not all predictions; six-decimal serialized scores."}
        checked += [aggregate, raw]
    for pin in checked: check(pin)
    return result, checked


def admit_arm(repo, plan, row, names, projection=None):
    label, dataset = row["label"], plan["dataset"]
    if label not in plan["arms"] or row.get("role") != plan["selection"]["roles"][label]:
        raise ValueError("Unexpected transfer arm or changed role")
    folder = Path(plan["output"]) / label
    graph_command, native_command = arm_commands(plan, folder, label)
    preflight = folder / "preflight.json"
    if json.loads(preflight.read_text()) != {"graph": graph_command, "native": native_command}:
        raise ValueError("Projection transfer preflight differs from fixed contrast")
    checked = [record(preflight), *row["checked"]]
    for stage, command in (("graph", graph_command), ("native", native_command)):
        checked += pilot_reader.execution_receipt(row[stage + "_execution"], command, folder, stage)
    graph_path = folder / "graph_metrics.json"
    if record(graph_path) != row["graph_report"]: raise ValueError("Transfer graph receipt differs")
    graph = json.loads(graph_path.read_text())
    activation = validate_activation(graph, label, "projection_recovery", folder=folder, gene_names=names, hmmer_tools=plan["hmmer_tools"])
    if activation != row["activation"]: raise ValueError("Transfer profile/backend activation differs")
    stages = [{"label": stage["label"], **{key: stage["output"][key] for key in ("bytes", "sha256")}} for stage in graph["stages"]]
    if [[s["label"], s["bytes"], s["sha256"]] for s in stages] != row["graph_stage_signature"]:
        raise ValueError("Reported graph stage signature differs")
    checked += [record(graph_path), *activation.get("checked", [])]
    for stage in graph["stages"]:
        check(stage["output"])
        read_partition(Path(stage["output"]["path"]), set(names))
        checked.append(stage["output"])
    working = folder / "candidate/orthohmm_working_res"
    for key, filename in (("candidate_partition", "orthohmm_edges_clustered.txt"), ("membership_constraints", "phylogeny_candidate_merges.json")):
        if row["candidate"][key] != record(working / filename): raise ValueError("Transfer candidate/constraints differ from arm namespace")
        checked.append(row["candidate"][key])
    if row["candidate"].get("seed_partition") != record(folder / "graph/orthogroups_profiles_refined.txt"):
        raise ValueError("Transfer candidate seed differs from final graph")
    checked.append(row["candidate"]["seed_partition"])
    candidate, constraints = (Path(row["candidate"][key]["path"]) for key in ("candidate_partition", "membership_constraints"))
    read_partition(candidate, set(names))
    load_membership_constraints(constraints, candidate)
    prediction = folder / "native/orthohmm_working_res/orthohmm_edges_clustered.txt"
    if validate_partition(candidate, prediction) != row["coverage"]: raise ValueError("Transfer partition coverage differs")
    native = replay_native(row, {"command": native_command})
    if (native["genes"], native["species"]) != followup.UNIVERSES[dataset] or native["summary"] != row["native_summary"]:
        raise ValueError("Transfer native universe/summary differs")
    tree = record(Path(native["phylogeny"]) / "species_tree.rooted.nwk")
    if tree != row["rooted_species_tree"]: raise ValueError("Transfer rooted tree receipt differs")
    reference = plan["case"] if label == "legacy_off" else projection if label == "projection_hmmer_score_only" else None
    if label == "projection_hmmer_score_only" and reference is None:
        raise ValueError("Measurement requires independently admitted repaired control")
    if reference is not None:
        followup.candidate_equivalent(row, reference)
        golden = reference["golden"] if label == "legacy_off" else {**reference["products"], "species_tree.rooted.nwk": reference["rooted_species_tree"]}
        if compare_golden(Path(native["phylogeny"]), golden) != row["golden_matches"]:
            raise ValueError("Exact transfer control native products differ")
    scoring_folder, participant, reuse_pins = scoring_context(repo, plan, row, projection)
    actual, details, score_pins = readback_scores(repo, dataset, row, native, scoring_folder, participant)
    vector = endpoint_vector(dataset, actual)
    if actual != row["scores"] or details != row["details"] or vector != row["metric_vector"]:
        raise ValueError("Transfer scores/populations/vector do not reproduce")
    if label == "legacy_off" and vector != plan["case"]["control"]["scores"]:
        raise ValueError("Canonical transfer vector differs from original admitted control")
    receipt = record(folder / "assessment_receipt.json")
    expected = {key: row[key] for key in ("scores", "details", "metric_vector")}
    expected["checked"] = [pin for pin in row["checked"] if pin != receipt]
    if row["checked"].count(receipt) != 1 or json.loads(Path(receipt["path"]).read_text()) != expected:
        raise ValueError("Saved assessment receipt differs from actual producer evidence")
    turnover = native_pair_turnover(Path(plan["case"]["case"]["baseline_output"]) / "orthohmm_phylogeny/orthohmm_pairwise_orthologs_confidence.tsv",
        Path(native["phylogeny"]) / "orthohmm_pairwise_orthologs_confidence.tsv")
    if {k: v for k, v in turnover.items() if k != "wall_s"} != {k: v for k, v in row["native_pair_turnover"].items() if k != "wall_s"}:
        raise ValueError("Transfer native pair/confidence turnover differs")
    populations, population_pins = functional_populations(scoring_folder, actual, details, score_pins) if dataset == "QfO" else ({}, [])
    checked += [receipt, tree, *native["checked"], *score_pins, *reuse_pins, *turnover["files"], *population_pins]
    checked = unique_pins(checked)
    for pin in checked: check(pin)
    return {"label": label, "dataset": dataset, "role": row["role"], "status": "independent_content_readback_verified",
        "source_revision": plan["scientific_revision"], "scores": vector, "details": details,
        "coverage": {key: native[key] for key in ("genes", "species", "groups", "pairs")},
        "products": native["products"], "rooted_species_tree": tree, "candidate": row["candidate"], "graph_stages": stages,
        "activation": activation, "native_pair_turnover": turnover, "functional_populations": populations,
        "assessment": {key: deepcopy(row[key]) for key in ("scores", "conversion", "metric_files", "references") if key in row},
        "assessment_receipt": receipt, "assessment_reuse": row.get("measurement_scoring_reuse", row.get("control_scoring_reuse")),
        "resources": {key: row[key] for key in ("graph_execution", "candidate_resources", "native_execution", "total_incremental_wall_s")},
        "checked": checked}


def verify_dataset(cells, chosen, *, campaign_name="projection"):
    scientific_revision = followup.configuration(campaign_name)["science"]
    if ([cell["label"] for cell in cells] != chosen["execution_arms"] or len({cell["dataset"] for cell in cells}) != 1
            or any(cell.get("status") != "independent_content_readback_verified" or cell.get("source_revision") != scientific_revision
                   or cell.get("role") != chosen["roles"][cell["label"]] or set(cell["scores"]) != set(cells[0]["scores"])
                   or cell["activation"].get("arm") != cell["label"]
                   or cell["activation"].get("projection_policy") != recovery_configuration(cell["label"])[0] for cell in cells)):
        raise ValueError("Require complete ordered source-bound transfer cells with original roles/settings")
    active = [cell for cell in cells if cell["label"] not in ("legacy_off", "projection_off")]
    if active and len({cell["activation"]["candidate_digest"] for cell in active}) != 1:
        raise ValueError("Active transfer candidate nominations differ")
    by_label = {cell["label"]: cell for cell in cells}
    if "projection_hmmer_score_only" in by_label:
        measurement, projection = by_label["projection_hmmer_score_only"], by_label["projection_off"]
        if (measurement["graph_stages"] != projection["graph_stages"]
                or recovery.signature(measurement["products"]) != recovery.signature(projection["products"])
                or recovery.signature({"tree": measurement["rooted_species_tree"]}) != recovery.signature({"tree": projection["rooted_species_tree"]})
                or recovery.signature({k: measurement["candidate"][k] for k in ("candidate_partition", "membership_constraints")})
                != recovery.signature({k: projection["candidate"][k] for k in ("candidate_partition", "membership_constraints")})
                or any(measurement[key] != projection[key] for key in ("scores", "details", "coverage", "functional_populations", "assessment"))
                or measurement["assessment_reuse"] != projection["assessment_receipt"]):
            raise ValueError("Measurement transfer differs from repaired control or reused assessment")
    return {label: {key: cell["scores"][key] - cells[0]["scores"][key] for key in cells[0]["scores"]}
            for label, cell in by_label.items() if label != "legacy_off"}


def assemble(pilot, cells, chosen):
    if ([cell["label"] for cell in pilot["cells"]] != list(followup.ARMS)
            or [(cell["dataset"], cell["label"]) for cell in cells] != [(dataset, label) for dataset in DATASETS for label in chosen["execution_arms"]]):
        raise ValueError("Require all original pilot arms and every declared transfer, without selected winners")
    methods = []
    for original in pilot["cells"]:
        label = original["label"]
        selected = [original, *[cell for cell in cells if cell["label"] == label]]
        vector = {key: value for cell in selected for key, value in cell["scores"].items()}
        metric_frontier({label: vector})
        role = chosen["roles"].get(label, "measurement_control" if label == "projection_hmmer_score_only" else "stopped_candidate")
        methods.append({"key": "projection_recovery:" + label, "label": "OrthoHMM projection " + label, "role": role,
            "status": "complete_development_vector" if len(selected) == 3 else "not_extended_after_original_pilot",
            "scores": vector, "cells": selected, "missing_endpoints": [key for key in ENDPOINTS if key not in vector]})
    control = methods[0]["scores"]
    for method in methods[1:]:
        if method["missing_endpoints"]: continue
        deltas = {key: method["scores"][key] - control[key] for key in ENDPOINTS}
        primary_gain = any(deltas[key] > 0 for key in ENDPOINTS[:3])
        losses = [key for key in ENDPOINTS if method["scores"][key] < control[key] - .005]
        original_pass = method["role"] == "advancing_candidate" and chosen["original_gates"][method["cells"][0]["label"]]["eligible"]
        method["full_vector_screen"] = {"eligible": original_pass and primary_gain and not losses,
            "original_pilot_eligible": bool(original_pass), "point_screen_pass": primary_gain and not losses,
            "deltas": deltas, "primary_f1_improvement": primary_gain, "declines_exceeding_0_005": losses,
            "scope": "Original development screen only. Applicability/control-only roles cannot become passing candidates. Conditional functional means, unseeded FAS and development exposure prevent automatic default promotion or superiority claims."}
    by_label = {method["cells"][0]["label"]: method for method in methods}
    contrasts = {"versus_projection": {}, "hmmer_versus_builtin": None}
    for label in ("projection_builtin_reselect", "projection_hmmer_reselect"):
        left, right = by_label[label]["scores"], by_label["projection_off"]["scores"]
        contrasts["versus_projection"][label] = {key: left[key] - right[key] for key in ENDPOINTS if key in left and key in right}
    left, right = by_label["projection_hmmer_reselect"]["scores"], by_label["projection_builtin_reselect"]["scores"]
    contrasts["hmmer_versus_builtin"] = {key: left[key] - right[key] for key in ENDPOINTS if key in left and key in right}
    contrasts["scope"] = "Matched construction/decision/backend point contrasts, not independent confirmation. Legacy-relative gains alone do not establish HMMER benefit. Missing transfer effects stay absent, not zero."
    return methods, contrasts


def reporter_inventory(driver, *, campaign_name="projection"):
    revision = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=driver, text=True).strip()
    files = sources(driver) + [record(driver / "benchmark_tools" / configuration(campaign_name)["protocol"]),
                              record(driver / followup.configuration(campaign_name)["protocol"])]
    if campaign_name != "projection":
        files += [record(driver / path) for path in ("benchmark_tools/" + PROTOCOL, followup.PROTOCOL,
                                                    "benchmark_tools/hmm_stockholm_annotation_recovery_20261019.json")]
    for pin in files:
        relative = Path(pin["path"]).relative_to(driver)
        if relative.suffix != ".so" and hashlib.sha256(subprocess.check_output(["git", "show", revision + ":" + str(relative)], cwd=driver)).hexdigest() != pin["sha256"]:
            raise ValueError("Transfer reader differs from committed bytes")
    return revision, files


def write_table(path, methods, context):
    with path.open("x") as stream:
        writer = csv.writer(stream, delimiter="\t", lineterminator="\n")
        writer.writerow(["method", "role", "evidence", *ENDPOINTS])
        for method in [*context, *methods]:
            writer.writerow([method["label"], method.get("role", "contextual_comparator"), method["status"], *[method["scores"].get(key, "NA") for key in ENDPOINTS]])


def run(repo, revision, output, *, campaign_name="projection"):
    if output.exists(): raise FileExistsError("Preserve every projection transfer readback")
    config, transfer = configuration(campaign_name), followup.configuration(campaign_name)
    entry, binding = (launch(repo, revision) if campaign_name == "projection"
        else launch(repo, revision, campaign_name=campaign_name))
    driver = Path(__file__).resolve().parents[1]
    reader_revision, reader_files = (reporter_inventory(driver) if campaign_name == "projection"
        else reporter_inventory(driver, campaign_name=campaign_name))
    if (reader_revision != entry["independent_reader"]["runner_commit"]
            or driver != relative_path(repo, entry["independent_reader"]["scientific_reader_source"])):
        raise ValueError("Actual transfer reader differs from immutable launch")
    spec = json.loads((driver / "benchmark_tools" / config["protocol"]).read_text())
    if (spec["scientific_revision"] != transfer["science"] or spec["driver_revision"] != config["driver"]
            or spec["preparer_revision"] != transfer["preparer"] or spec["datasets"] != list(DATASETS)
            or spec["screen"] != {"require_complete_eight_endpoints": True, "primary_f1_endpoints": list(ENDPOINTS[:3]),
                "require_at_least_one_primary_f1_improvement": True, "maximum_endpoint_point_decline": .005,
                "preserve_original_pilot_eligibility": True, "interpretation": spec["screen"]["interpretation"]}):
        raise ValueError("Changed prespecified transfer readback protocol")
    if campaign_name != "projection" and (spec.get("campaign") != campaign_name or spec.get("parent") != PROTOCOL):
        raise ValueError("Changed explicit Stockholm readback contract or original parent")
    prepared_pin = {key: entry["prepared_inputs"][key] for key in ("path", "bytes", "sha256")}
    prepared_pin["path"] = str(relative_path(repo, prepared_pin["path"]))
    check(prepared_pin)
    prepared = json.loads(Path(prepared_pin["path"]).read_text())
    if prepared.get("campaign", "projection") != campaign_name:
        raise ValueError("Transfer readback and preparation belong to different explicit campaigns")
    audit = followup.authenticate_prepared(repo, prepared)
    pilot, _, chosen, checked, provenance = (planner.authenticate(repo) if campaign_name == "projection"
        else planner.authenticate(repo, campaign_name=campaign_name))
    if audit["selection"] != chosen or chosen != entry["selection"]: raise ValueError("Transfer launch changes original eligibility/control roles")
    context, context_pins = pilot_reader.retained_context(repo)
    checked += [prepared_pin, *prepared["inputs"], *reader_files, *context_pins]
    output.mkdir(parents=True)
    result = {"status": "running", "launch_ledger": binding, "reporter_revision": reader_revision,
        "reporter_files": reader_files, "protocol": record(driver / "benchmark_tools" / config["protocol"]),
        "prepared_inputs": prepared_pin, "pilot_provenance": provenance, "selection": chosen, "datasets": [], "cells": [],
        "context": context, "orthobench_uncertainty": pilot["uncertainty"], "accuracy_result_admitted": False,
        "defaults_promoted": False, "goal_completed": False,
        "scope": "Independent content/score readback using shared tested helpers, not independent biological confirmation. Development-exposed panels; TreeFam supporting, Three Kingdoms supplementary. Conditional GO/EC/FAS populations and unseeded FAS limit claims. No fabricated QfO family or IID pair intervals. Cached shared-host incremental timings are not isolated fresh-pipeline costs or aggregate memory. Other campaign mechanisms and original failures remain retained in their separate source-bound reports; this transfer report alone does not close the goal."}
    if campaign_name != "projection": result["campaign"] = campaign_name
    write_json(output / "summary.json", result)
    try:
        for case_entry in entry["datasets"]:
            plan_pin = {key: case_entry["execution_plan"][key] for key in ("path", "bytes", "sha256")}
            plan_pin["path"] = str(relative_path(repo, plan_pin["path"]))
            check(plan_pin)
            plan = json.loads(Path(plan_pin["path"]).read_text())
            root = Path(plan["output"])
            producer = record(root / "results.json")
            report = json.loads(Path(producer["path"]).read_text())
            require_completed(plan, report, case_entry, chosen)
            if (report["execution_plan"] != plan_pin or Path(plan_pin["path"]).parent != root or Path(plan["repo"]) != repo
                    or Path(plan["driver_root"]) != relative_path(repo, entry["source"]) or plan["prepared_inputs"] != prepared_pin
                    or any(plan[key] != value for key, value in audit.items())):
                raise ValueError("Transfer plan/producer differs from immutable launch and original preparation")
            accounting = pilot_reader.terminal(case_entry["inference"]["job_id"])
            for key, field in (("scientific_source", "scientific_revision"), ("driver_root", "driver_revision")):
                source = Path(plan[key])
                if subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=source, text=True).strip() != plan[field]:
                    raise ValueError("Frozen transfer source moved")
            source = Path(plan["scientific_source"])
            if (sources(source) + [record(source / "benchmark_tools/run_profile_weighting_pilot.slurm")] != plan["scientific_sources"]
                    or [record(Path(plan["driver_root"]) / path) for path in (
                        DRIVER_FILES if campaign_name == "projection" else followup.driver_paths(campaign_name))] != plan["driver_sources"]
                    or runtime() != plan["runtime"]):
                raise ValueError("Frozen transfer source/runtime inventory differs")
            expected_case, _ = planner.selected_case(SimpleNamespace(record=record, check=check,
                verified_cell=planner.verified_cell, PRODUCTS=planner.PRODUCTS), repo, plan["dataset"])
            if expected_case != plan["case"]: raise ValueError("Transfer case differs from original admitted canonical case")
            pins = [plan_pin, producer, *plan["inputs"], *plan["driver_sources"]]
            for pin in unique_pins(pins): check(pin)
            names = weighting.ownership(plan, report)
            cells = []
            for row in report["arms"]:
                projection = next((cell for cell in cells if cell["label"] == "projection_off"), None)
                cells.append(admit_arm(repo, plan, row, names, projection))
            deltas = verify_dataset(cells, chosen, campaign_name=campaign_name)
            result["cells"] += cells
            result["datasets"].append({"dataset": plan["dataset"], "execution_plan": plan_pin, "producer": producer,
                "scheduler": accounting, "deltas": deltas})
            checked += [*pins, *[p for cell in cells for p in cell["checked"]]]
            write_json(output / "summary.json", result)
        result["sequential_science"] = sequential_science(entry["datasets"])
        result["methods"], result["matched_contrasts"] = assemble(pilot, result["cells"], chosen)
        write_table(output / "scores.tsv", result["methods"], context)
        result["score_table"] = record(output / "scores.tsv")
        checked.append(result["score_table"])
        result["checked"] = unique_pins(checked)
        for pin in result["checked"]: check(pin)
        if runtime() != prepared["runtime"]: raise ValueError("Runtime changed during transfer readback")
        result["frontier"] = metric_frontier({method["key"]: method["scores"] for method in [*context, *result["methods"]]})
        result["eligible_labels"] = [method["cells"][0]["label"] for method in result["methods"] if method.get("full_vector_screen", {}).get("eligible")]
        result.update(status="projection_full_vectors_independently_read_back", accuracy_result_admitted=True)
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
    parser.add_argument("--campaign", choices=("projection", "stockholm_annotation"), default="projection")
    args = parser.parse_args()
    print(run(args.repo.resolve(), args.launch_revision, args.output.resolve(), campaign_name=args.campaign)["status"])
