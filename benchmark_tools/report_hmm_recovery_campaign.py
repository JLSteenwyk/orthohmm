"""Collate admitted five-arm recovery with original successful HMM experiments."""

import argparse
import csv
import hashlib
import io
import json
from pathlib import Path
import subprocess
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from benchmark_tools import report_hmm_evidence_campaign as base
from benchmark_tools.prepare_ob_candidate_neighborhood import check, record
from benchmark_tools.run_membership_hillclimb import write_json
from benchmark_tools.run_three_kingdoms_membership_supplement import sources


LAUNCH = "a884d9c0d2d7b3762e0ef2960bfb411bc847837a"
PROTOCOL = "hmm_recovery_campaign_report_20261019.json"
REQUIRED = {"weighting": "weighting_full_vectors_independently_read_back",
            "local": "local_joint_pilot_independently_read_back",
            "recovery": "projection_recovery_independently_read_back"}
FAS_EVIDENCE = "benchmark_tools/results/completed_weighting_local_readback_20261019.json"


def launch(repo):
    old, old_binding = base.launch(repo)
    blob = subprocess.check_output(["git", "show", LAUNCH + ":" + base.LEDGER], cwd=repo)
    ledger = json.loads(blob)
    entry = ledger["projection_recovery_postrun_readback"]
    science = ledger["projection_recovery_pilot"]
    if (entry["job"].get("submitted") is not True or science["inference"].get("submitted") is not True
            or entry["runner_commit"] != science["independent_reader"]["runner_commit"]
            or entry["scientific_reader_source"] != science["independent_reader"]["scientific_reader_source"]
            or entry["job"]["dependency"] != "afterany:" + science["inference"]["job_id"]):
        raise ValueError("Recovery reader/science launch binding differs")
    return {"weighting": old["weighting"], "local": old["local"], "recovery": entry}, science, {
        "historical": old_binding, "recovery": {"revision": LAUNCH, "path": base.LEDGER,
                                                 "git_blob_sha256": hashlib.sha256(blob).hexdigest()}}


def failure_history(repo, spec):
    pins = [{**pin, "path": str((repo / pin["path"]).resolve())} for pin in spec["old_failure_records"]]
    for pin in pins: check(pin)
    producer = json.loads(Path(pins[0]["path"]).read_text())
    if producer.get("job_id") != "24493" or producer.get("status") != "finished_with_explicit_gaps":
        raise ValueError("Original failed HMMER history differs")
    jobs = []
    for job in ("24493", "24498", "24501"):
        command = ["sacct", "-j", job, "-X", "-P", "-o", "JobIDRaw,State,ExitCode"]
        rows = list(csv.DictReader(io.StringIO(subprocess.check_output(command, text=True)), delimiter="|"))
        if rows != [{"JobIDRaw": job, "State": "FAILED", "ExitCode": "1:0"}]:
            raise ValueError("Original failure accounting differs")
        jobs.append({"command": command, "accounting": rows[0]})
    return {"status": "original_failures_preserved_not_accuracy_evidence", "jobs": jobs, "files": pins}, pins


def collate(weighting, local, recovery, producer):
    result = base.collate(weighting, local, recovery, producer, recovery=True)
    for population in result["functional_population_comparisons"]:
        if population["endpoint"] == "QfO: FAS":
            population["population_count_semantics"] = (
                "Native NR_ORTHOLOGS eligible pool, not the numeric scored-sample count. "
                "FAS means use saved unseeded samples with repeated proteins; no matched-population or pair-independent inference."
            )
    result["methods"].append({"key": "selective_original:failed_panel",
        "label": "OrthoHMM selective HMMER original pilot (unadmitted)", "role": "failed_execution",
        "status": "failed_incomplete_panel_not_accuracy_evidence", "scores": {}, "cells": [],
        "missing_endpoints": list(base.ENDPOINTS),
        "interpretation": "Original source-integrity failure, not inactivity or a negative accuracy result. Replaced only by a separately controlled recovery, not overwritten."})
    result["frontier"] = base.metric_frontier({m["key"]: m["scores"] for m in [*result["context"], *result["methods"]]})
    result["scope"] += " Original failed HMMER pilot remains explicitly unadmitted; report does not redirect the original failed collator."
    return result


def retained_fas_samples(repo, weighting):
    blob = subprocess.check_output(["git", "show", LAUNCH + ":" + FAS_EVIDENCE], cwd=repo)
    pin = record(repo / FAS_EVIDENCE)
    if hashlib.sha256(blob).hexdigest() != pin["sha256"]:
        raise ValueError("Previously audited FAS population evidence changed")
    populations = json.loads(blob)["weighting"]["functional_populations"]
    qfo = {cell["label"]: cell for cell in weighting["cells"] if cell["dataset"] == "QfO"}
    if [row["mode"] for row in populations] != ["none", "henikoff", "henikoff_cluster"] or set(qfo) != {row["mode"] for row in populations}:
        raise ValueError("FAS population inventory differs from complete weighting panel")
    for row in populations:
        if (row["FAS_reported_eligible_pool"] != qfo[row["mode"]]["details"]["FAS"]["assessed_relations"]
                or type(row["FAS_numeric_sample_pairs"]) is not int
                or not 0 < row["FAS_numeric_sample_pairs"] <= row["FAS_reported_eligible_pool"]):
            raise ValueError("FAS eligible population or numeric sample count differs")
    return {"source": pin, "populations": populations,
            "scope": "Previously audited saved numeric samples, not new sampling or independent confidence evidence. Eligible pools differ from scored pairs; FAS is unseeded and repeated proteins/families prevent simple IID pair inference."}, [pin]


def reporter_inventory(driver):
    revision = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=driver, text=True).strip()
    pins = sources(driver) + [record(driver / "benchmark_tools" / PROTOCOL)]
    for pin in pins:
        relative = Path(pin["path"]).relative_to(driver)
        if relative.suffix != ".so" and hashlib.sha256(subprocess.check_output(
                ["git", "show", revision + ":" + str(relative)], cwd=driver)).hexdigest() != pin["sha256"]:
            raise ValueError("Recovery campaign reporter differs from committed source")
    return revision, pins


def run(repo, output):
    if output.exists(): raise FileExistsError("Preserve every recovery campaign report namespace")
    driver = Path(__file__).resolve().parents[1]
    revision, reporter_files = reporter_inventory(driver)
    spec = json.loads((driver / "benchmark_tools" / PROTOCOL).read_text())
    if spec["historical_launch"] != base.LAUNCH or spec["recovery_launch"] != LAUNCH:
        raise ValueError("Changed recovery campaign launch protocol")
    entries, science, binding = launch(repo)
    reports, checked, provenance = {}, list(reporter_files), {}
    for name, status in REQUIRED.items():
        report, pins, receipt = base.read_admission(repo, entries[name], status)
        reports[name], provenance[name] = report, receipt
        checked += pins
    recovery = reports["recovery"]
    if (recovery["launch_ledger"] != entries["recovery"]["launch_binding"]
            or recovery["source_revision"] != science["runner_commit"]):
        raise ValueError("Admitted recovery reader uses a different original science launch")
    checked.append(recovery["producer"])
    for pin in base.unique_pins(checked): check(pin)
    producer = json.loads(Path(recovery["producer"]["path"]).read_text())
    if producer.get("job_id") != science["inference"]["job_id"]:
        raise ValueError("Admitted recovery producer is not the original scientific job")
    result = collate(reports["weighting"], reports["local"], recovery, producer)
    fas, fas_pins = retained_fas_samples(repo, reports["weighting"])
    result["retained_fas_samples"] = fas
    checked += fas_pins
    history, history_pins = failure_history(repo, spec)
    checked += history_pins
    output.mkdir(parents=True)
    base.write_table(output / "scores.tsv", result["methods"], result["context"])
    result.update(status="admitted_recovery_campaign_collated_followups_explicit", launch_ledger=binding,
        original_failure_history=history, inputs=provenance, reporter_revision=revision,
        reporter_files=reporter_files, protocol=record(driver / "benchmark_tools" / PROTOCOL),
        score_table=record(output / "scores.tsv"), defaults_promoted=False, goal_completed=False)
    checked.append(result["score_table"])
    result["checked"] = base.unique_pins(checked)
    for pin in result["checked"]: check(pin)
    write_json(output / "summary.json", result)
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    print(run(args.repo.resolve(), args.output.resolve())["status"])
