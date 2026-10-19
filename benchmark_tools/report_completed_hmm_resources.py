"""Read completed HMM experiment costs without presenting a timing leaderboard."""

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import subprocess
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from benchmark_tools.prepare_ob_candidate_neighborhood import check, record
from benchmark_tools.report_hmm_evidence_campaign import unique_pins
from benchmark_tools.report_profile_evidence_pilot import execution_receipt
from benchmark_tools.run_membership_hillclimb import write_json
from benchmark_tools.run_profile_weighting_followup import successful_job


ADMISSION = "benchmark_tools/results/completed_weighting_local_readback_20261019.json"
ADMISSION_REVISION = "5aa0b86ad246e0e9165fefdf344fe1a3f4aa1727"
STATUS = {"weighting": "weighting_full_vectors_independently_read_back",
          "local": "local_joint_pilot_independently_read_back"}
FIELDS = ("method", "dataset", "phase", "elapsed_seconds", "cpu_seconds", "cpu_scope",
          "max_process_rss_gib", "parent_lifetime_max_rss_gib", "evidence", "scope")


def nonnegative(value):
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value) or value < 0:
        raise ValueError("Require finite nonnegative recorded resource values")
    return value


def cell_rows(method, cell, admitted, *, local=False):
    if cell.get("status") != "independent_content_readback_verified":
        raise ValueError("Cannot report costs as completed from an unadmitted cell")
    resources = cell["resources"]
    expected = {"native_execution"} if local else {
        "graph_execution", "native_execution", "candidate_resources", "total_incremental_wall_s"}
    if set(resources) != expected:
        raise ValueError("Unknown or incomplete cost scope; do not fill absent phases with zero")
    rows, checked = [], []
    for stage in (("native",) if local else ("graph", "native")):
        execution = resources[stage + "_execution"]
        command = execution["command"]
        if command.count("--cpu") != 1 or command[command.index("--cpu") + 1] != "16":
            raise ValueError("Recorded experiment CPU allocation differs")
        folder = Path(execution["time"]["path"]).parent
        pins = execution_receipt(execution, command, folder, stage)
        if any(pin not in admitted and pin not in cell.get("checked", []) for pin in pins):
            raise ValueError("Raw resource/log receipt absent from prior independent admission")
        checked += pins
        values = execution["resources"]
        scope = ("Cached final native replay only; excludes nomination, joint alignment/tree preparation and initial search."
                 if local else "Cached graph/profile replay; excludes initial search." if stage == "graph"
                 else "Native phylogeny replay with compatible retained caches; not fresh full-pipeline inference.")
        rows.append(dict(method=method, dataset=cell["dataset"], phase=stage,
            elapsed_seconds=nonnegative(values["elapsed_seconds"]),
            cpu_seconds=nonnegative(values["user_seconds"]) + nonnegative(values["system_seconds"]),
            cpu_scope="native_child_and_waited_descendants",
            max_process_rss_gib=nonnegative(values["max_process_rss_kib"]) / 1024 ** 2,
            parent_lifetime_max_rss_gib=None, evidence="raw_GNU_time_reparsed_and_log_rehashed", scope=scope))
    if not local:
        candidate = resources["candidate_resources"]
        rows.append(dict(method=method, dataset=cell["dataset"], phase="candidate_parent",
            elapsed_seconds=nonnegative(candidate["wall_s"]), cpu_seconds=nonnegative(candidate["cpu_s"]),
            cpu_scope="parent_process_only", max_process_rss_gib=None,
            parent_lifetime_max_rss_gib=nonnegative(candidate["peak_parent_process_rss_GiB"]),
            evidence="retained_producer_metadata_not_independently_retimed",
            scope="In-process candidate preparation; RSS is parent lifetime maximum, not this phase's incremental or aggregate memory."))
        rows.append(dict(method=method, dataset=cell["dataset"], phase="arm_total",
            elapsed_seconds=nonnegative(resources["total_incremental_wall_s"]), cpu_seconds=None,
            cpu_scope="unmeasured", max_process_rss_gib=None, parent_lifetime_max_rss_gib=None,
            evidence="retained_producer_metadata_not_independently_retimed",
            scope="Shared-host incremental arm including verification and assessment/reuse; not end-to-end tool runtime. Do not sum with component phases."))
    return rows, checked


def run(repo, output):
    if output.exists(): raise FileExistsError("Preserve every resource readback")
    source = Path(__file__).resolve()
    relative = source.relative_to(repo)
    revision = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=repo, text=True).strip()
    if hashlib.sha256(subprocess.check_output(["git", "show", revision + ":" + str(relative)], cwd=repo)).hexdigest() != record(source)["sha256"]:
        raise ValueError("Resource reader must match committed source")
    blob = subprocess.check_output(["git", "show", ADMISSION_REVISION + ":" + ADMISSION], cwd=repo)
    proof_pin = record(repo / ADMISSION)
    if hashlib.sha256(blob).hexdigest() != proof_pin["sha256"]:
        raise ValueError("Original completed-experiment admission record changed")
    proof = json.loads(blob)["integration_readback"]
    rows, checked, inputs = [], [proof_pin, record(source)], {}
    for kind, status in STATUS.items():
        original = proof[kind]
        pin = {**original["summary"], "path": str(repo / original["summary"]["path"])}
        check(pin)
        report = json.loads(Path(pin["path"]).read_text())
        if (report.get("status") != status or report.get("accuracy_result_admitted") is not True
                or report.get("reporter_revision") != original["reader_revision"]
                or report.get("defaults_promoted") is not False or report.get("goal_completed") is not False):
            raise ValueError("Require the original complete independently admitted report")
        accounting = successful_job(original["job_id"], 1, "16G")
        checked.append(pin)
        inputs[kind] = {"report": pin, "reader_revision": report["reporter_revision"], "scheduler": accounting}
        methods = report["methods"] if kind == "weighting" else [
            {"key": "local:" + cell["label"], "cells": [cell]} for cell in report["cells"]]
        for method in methods:
            for cell in method["cells"]:
                extra_rows, pins = cell_rows(method["key"], cell, report["checked"], local=kind == "local")
                rows += extra_rows
                checked += pins
        del report
    identities = [(row["method"], row["dataset"], row["phase"]) for row in rows]
    if len(set(identities)) != len(identities): raise ValueError("Duplicate resource population")
    checked = unique_pins(checked)
    for pin in checked: check(pin)
    output.mkdir(parents=True)
    with (output / "resources.tsv").open("x") as stream:
        writer = csv.DictWriter(stream, fieldnames=FIELDS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows({key: "NA" if row[key] is None else row[key] for key in FIELDS} for row in rows)
    result = {"status": "completed_hmm_incremental_resource_receipts_read_back",
        "admission_revision": ADMISSION_REVISION, "admission": proof_pin,
        "source_revision": revision, "source": record(source), "inputs": inputs, "rows": rows, "checked": checked,
        "resource_table": record(output / "resources.tsv"), "new_biological_inference": False,
        "defaults_promoted": False, "goal_completed": False,
        "scope": "Prior source-bound accuracy admission retained, not a new accuracy/native-output re-audit. Raw saved GNU-time/log receipts freshly verified; candidate/arm totals are retained producer metadata, not independent retiming. Shared-host cached phase costs, not isolated fresh-pipeline costs, tool rankings or aggregate simultaneous memory. No sums of RSS or arm/component times; local nomination/joint-tree preparation excluded. Unadmitted/in-progress HMMER and contextual tools are not silently assigned these costs."}
    write_json(output / "summary.json", result)
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    print(run(args.repo.resolve(), args.output.resolve())["status"])
