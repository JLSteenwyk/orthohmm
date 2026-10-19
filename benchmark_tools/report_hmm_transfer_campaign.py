"""Join admitted individual HMM experiments and complete projection transfers."""

import argparse
import csv
import hashlib
import io
import json
from pathlib import Path
import subprocess
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from benchmark_tools import plan_projection_followup as planner
from benchmark_tools import report_hmm_evidence_campaign as base
from benchmark_tools import report_hmm_recovery_campaign as recovery
from benchmark_tools import report_projection_followup as transfer
from benchmark_tools.prepare_ob_candidate_neighborhood import check, record
from benchmark_tools.report_completed_hmm_resources import FIELDS, cell_rows
from benchmark_tools.run_membership_hillclimb import write_json
from benchmark_tools.run_three_kingdoms_membership_supplement import sources


PROTOCOL = "hmm_transfer_campaign_report_20261019.json"
READBACK_LAUNCH = "4db3320036c388af184c39246d49b37860602dec"
SCIENCE = "8f0cd5605842e470e64f40c31bc608c37945b8b3"
READER = "840746675e9dc1e953c0adee7dba158182922b19"


def admit_transfer(repo):
    blob = subprocess.check_output(["git", "show", READBACK_LAUNCH + ":" + base.LEDGER], cwd=repo)
    entry = json.loads(blob)["stockholm_projection_followup_postrun_readback"]
    if (entry.get("runner_commit") != READER
            or entry.get("scientific_reader_source") != "benchmarks/work/stockholm_projection_followup_readback_v1_20261019"
            or entry["job"].get("submitted") is not True
            or (entry["job"].get("cpus"), entry["job"].get("memory_allocation"), entry["job"].get("wall_limit")) != (1, "16G", "02:00:00")):
        raise ValueError("Changed separately frozen transfer reader submission")
    report, pins, provenance = base.read_admission(repo, entry, "projection_full_vectors_independently_read_back")
    science, binding = transfer.launch(repo, entry["launch_binding"]["revision"], campaign_name="stockholm_annotation")
    if (report.get("campaign") != "stockholm_annotation" or report["launch_ledger"] != binding
            or binding != entry["launch_binding"]
            or entry["job"].get("dependency") != "afterany:" + science["datasets"][-1]["inference"]["job_id"]
            or science["independent_reader"]["runner_commit"] != entry["runner_commit"]
            or science["independent_reader"]["scientific_reader_source"] != entry["scientific_reader_source"]):
        raise ValueError("Transfer report/reader uses a different immutable scientific launch")
    return report, pins, {"reader": provenance, "science": binding,
        "submission": {"revision": READBACK_LAUNCH, "path": base.LEDGER, "git_blob_sha256": hashlib.sha256(blob).hexdigest()}}


def merge_transfer(result, pilot, producer, following):
    if (following.get("status") != "projection_full_vectors_independently_read_back"
            or following.get("accuracy_result_admitted") is not True
            or following.get("defaults_promoted") is not False or following.get("goal_completed") is not False
            or following["context"] != result["context"] or following["orthobench_uncertainty"] != pilot["uncertainty"]):
        raise ValueError("Require complete transfer admission with the same context and pilot uncertainty")
    selected = planner.selection(pilot, producer, scientific_revision=SCIENCE, reader_revision=SCIENCE)
    if following["selection"] != selected:
        raise ValueError("Transfer roles or selection differ from original pilot gates")
    methods, contrasts = transfer.assemble(pilot, following["cells"], selected)
    if any(m["missing_endpoints"] for m in methods if m["cells"][0]["label"] in selected["execution_arms"]):
        raise ValueError("Every selected transfer configuration requires all eight endpoints")
    if (methods != following["methods"] or contrasts != following["matched_contrasts"]
            or following["eligible_labels"] != [m["cells"][0]["label"] for m in methods if m.get("full_vector_screen", {}).get("eligible")]):
        raise ValueError("Transfer vectors, contrasts or original screens do not reproduce")
    replacement = {m["key"]: m for m in methods}
    original = {m["key"]: m for m in result["methods"] if m["key"].startswith("projection_recovery:")}
    if set(replacement) != set(original) or any(replacement[k]["cells"][0] != original[k]["cells"][0] for k in original):
        raise ValueError("Transfer does not retain every exact original pilot cell")
    result["methods"] = [{**m, **replacement[m["key"]], "pilot_role": m["role"]} if m["key"] in replacement else m
                         for m in result["methods"]]
    by_key = {m["key"]: m for m in result["methods"]}
    if by_key["projection_recovery:legacy_off"]["scores"] != by_key["profile_weighting:none"]["scores"]:
        raise ValueError("Canonical control differs between individually admitted mechanisms")
    result["required_next_actions"] = [a for a in result["required_next_actions"] if a["method"] not in replacement]
    for method in methods:
        if method["role"] != "advancing_candidate": continue
        screen = method["full_vector_screen"]
        result["required_next_actions"].append({"method": method["key"], "action": (
            "Original full-vector development screen passes. Consider only justified combinations with ablations after individual effects; freeze candidate before independent confirmation or document unavailable data. No automatic default promotion."
            if screen["eligible"] else "Original full-vector screen fails. Stop this fixed branch without retuning; retain complete negative evidence.")})
    qfo = {c["label"]: c for c in following["cells"] if c["dataset"] == "QfO"}
    control = qfo["legacy_off"]
    for label, cell in qfo.items():
        if label == "legacy_off": continue
        for metric in ("EC", "GO", "FAS"):
            before, after = control["details"][metric], cell["details"][metric]
            old_pop, new_pop = control["functional_populations"][metric], cell["functional_populations"][metric]
            count = "sample_pairs" if metric == "FAS" else "numeric_scored_pairs"
            result["functional_population_comparisons"].append({"method": "projection_recovery:" + label,
                "role": cell["role"], "endpoint": "QfO: " + metric,
                "control_assessed_relations": before["assessed_relations"], "candidate_assessed_relations": after["assessed_relations"],
                "control_numeric_scored_pairs": old_pop[count], "candidate_numeric_scored_pairs": new_pop[count],
                "control_mean": before["similarity"], "candidate_mean": after["similarity"],
                "scope": "Inherited independently read numeric populations; FAS eligible pools are not sample sizes. Conditional unseeded means with repeated proteins, not matched-population orthology accuracy or IID biological uncertainty."})
    result["matched_transfer_contrasts"] = contrasts
    result["frontier"] = base.metric_frontier({m["key"]: m["scores"] for m in [*result["context"], *result["methods"]]})
    result["recommendation"] = "Keep production defaults unchanged. Retain original stopped branches and source-bound development candidates; complete justified-combination assessment, post-freeze confirmation or documented data limits and final recommendation/remote audit. Full-vector point screens alone do not establish generalizability or superiority."
    result["scope"] += " Complete projection transfers now replace pilot-only rows without reviving stopped candidates. This compilation does not complete the overall goal."
    return result


def projection_failures(repo, spec):
    pins = [{**p, "path": str(repo / p["path"])} for p in spec["projection_failure_records"]]
    for pin in pins: check(pin)
    failed = json.loads(Path(pins[0]["path"]).read_text())
    if failed.get("job_id") != "24514" or failed.get("status") != "finished_with_explicit_gaps":
        raise ValueError("Original projection failure history changed")
    jobs = []
    for job in ("24514", "24515", "24516", "24517"):
        command = ["sacct", "-j", job, "-X", "-P", "-o", "JobIDRaw,State,ExitCode"]
        rows = list(csv.DictReader(io.StringIO(subprocess.check_output(command, text=True)), delimiter="|"))
        if rows != [{"JobIDRaw": job, "State": "FAILED", "ExitCode": "1:0"}]:
            raise ValueError("Original projection failure accounting changed")
        jobs.append({"command": command, "accounting": rows[0]})
    return {"jobs": jobs, "files": pins, "scope": "Failed execution, never negative accuracy or inactivity."}, pins


def resource_table(output, methods):
    rows, pins = [], []
    for method in methods:
        for cell in method["cells"]:
            values, checked = cell_rows(method["key"], cell, cell["checked"], local=method["key"].startswith("local:"))
            rows += values
            pins += checked
    with (output / "resources.tsv").open("x") as stream:
        writer = csv.DictWriter(stream, fieldnames=FIELDS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows({k: "NA" if row[k] is None else row[k] for k in FIELDS} for row in rows)
    return rows, pins


def reporter_inventory(driver):
    revision = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=driver, text=True).strip()
    pins = sources(driver) + [record(driver / "benchmark_tools" / name) for name in (PROTOCOL, recovery.PROTOCOL)]
    for pin in pins:
        relative = Path(pin["path"]).relative_to(driver)
        if relative.suffix != ".so" and hashlib.sha256(subprocess.check_output(
                ["git", "show", revision + ":" + str(relative)], cwd=driver)).hexdigest() != pin["sha256"]:
            raise ValueError("Combined campaign reporter differs from committed source")
    return revision, pins


def run(repo, output):
    if output.exists(): raise FileExistsError("Preserve every combined campaign readback")
    driver = Path(__file__).resolve().parents[1]
    revision, reporter_files = reporter_inventory(driver)
    spec = json.loads((driver / "benchmark_tools" / PROTOCOL).read_text())
    if (spec["historical_launch"] != base.LAUNCH or spec["recovery_reader_submission"] != planner.STOCKHOLM_LAUNCH
            or spec["transfer_reader_submission"] != READBACK_LAUNCH or spec["scientific_revision"] != SCIENCE
            or spec["transfer_reader_revision"] != READER):
        raise ValueError("Changed prespecified complete campaign bindings")
    following, checked, transfer_provenance = admit_transfer(repo)
    old, old_binding = base.launch(repo)
    new, science, new_binding = planner.launch(repo, "stockholm_annotation")
    entries = {"weighting": old["weighting"], "local": old["local"], "recovery": new["recovery"]}
    reports, provenance = {}, {"transfer": transfer_provenance}
    checked += reporter_files
    for key, status in recovery.REQUIRED.items():
        reports[key], pins, provenance[key] = base.read_admission(repo, entries[key], status)
        checked += pins
    pilot = reports["recovery"]
    check(pilot["producer"])
    producer = json.loads(Path(pilot["producer"]["path"]).read_text())
    if (pilot["launch_ledger"] != entries["recovery"]["launch_binding"]
            or pilot["source_revision"] != science["runner_commit"] or producer["job_id"] != science["inference"]["job_id"]):
        raise ValueError("Complete recovery admission changed its immutable scientific launch")
    result = recovery.collate(reports["weighting"], reports["local"], pilot, producer)
    result = merge_transfer(result, pilot, producer, following)
    fas, fas_pins = recovery.retained_fas_samples(repo, reports["weighting"])
    history, history_pins = recovery.failure_history(repo, spec)
    prior, prior_pins = projection_failures(repo, spec)
    result["methods"].append({"key": "projection_original:failed_panel", "label": "OrthoHMM prior projection recovery (failed execution)",
        "role": "failed_execution", "status": "failed_incomplete_panel_not_accuracy_evidence", "scores": {}, "cells": [],
        "missing_endpoints": list(base.ENDPOINTS)})
    checked += [pilot["producer"], *fas_pins, *history_pins, *prior_pins]
    for pin in base.unique_pins(checked): check(pin)
    output.mkdir(parents=True)
    base.write_table(output / "scores.tsv", result["methods"], result["context"])
    costs, cost_pins = resource_table(output, result["methods"])
    result.update(status="admitted_individual_hmm_methods_and_projection_transfers_collated",
        inputs=provenance, launch_ledger={"historical": old_binding, "recovery": new_binding, "transfer": transfer_provenance},
        retained_fas_samples=fas, original_failure_history=history, projection_failure_history=prior,
        reporter_revision=revision, reporter_files=reporter_files, protocol=record(driver / "benchmark_tools" / PROTOCOL),
        score_table=record(output / "scores.tsv"), resource_table=record(output / "resources.tsv"), resource_rows=costs,
        resource_scope="Cached shared-host phase costs; cache reuse differs. No end-to-end speed ratios, summed RSS or simultaneous aggregate memory. Local earlier nomination/joint-tree costs excluded; producer parent/total metadata not independently retimed.",
        accuracy_result_admitted=True, defaults_promoted=False, goal_completed=False)
    result["frontier"] = base.metric_frontier({m["key"]: m["scores"] for m in [*result["context"], *result["methods"]]})
    result["checked"] = base.unique_pins([*checked, *cost_pins, result["score_table"], result["resource_table"]])
    for pin in result["checked"]: check(pin)
    write_json(output / "summary.json", result)
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    print(run(args.repo.resolve(), args.output.resolve())["status"])
