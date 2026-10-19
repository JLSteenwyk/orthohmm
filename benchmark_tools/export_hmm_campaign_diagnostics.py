"""Export admitted campaign endpoint details without new scoring or sampling."""

import argparse
import csv
import hashlib
import json
from pathlib import Path
import subprocess
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from benchmark_tools import export_hillclimb_diagnostics as diagnostics
from benchmark_tools import report_hmm_evidence_campaign as base
from benchmark_tools import report_hmm_transfer_campaign as campaign
from benchmark_tools.prepare_ob_candidate_neighborhood import check, record
from benchmark_tools.run_membership_hillclimb import write_json


LAUNCH = "6ed21295ab918a4e3de5f65aeef3943545941b13"
COMPILER = "d21d076cf0ce817a7c6807b4c9512532c9b36dcd"
RETAINED = "benchmark_tools/results/current_benchmark_scores_20260926_v2/manifest.json"
STATUS = "admitted_individual_hmm_methods_and_projection_transfers_collated"
COLUMNS = (*diagnostics.COLUMNS, "method_role", "method_status", "evidence_scope", "reference_groups",
           "numeric_scored_pairs", "numeric_population_source_document", "numeric_population_source_location")


def admit(repo):
    blob = subprocess.check_output(["git", "show", LAUNCH + ":" + base.LEDGER], cwd=repo)
    entry = json.loads(blob)["complete_individual_effects_compilation"]
    if (entry["runner_commit"] != COMPILER
            or entry["scientific_reader_source"] != "benchmarks/work/hmm_transfer_campaign_report_v1_20261019"
            or entry["transfer_reader_submission"] != campaign.READBACK_LAUNCH
            or any(entry["job"].get(k) != v for k, v in {"submitted": True, "job_id": "24527",
                "dependency": "afterany:24525", "cpus": 1, "memory_allocation": "16G", "wall_limit": "02:00:00"}.items())):
        raise ValueError("Changed immutable complete campaign report submission")
    report, pins, provenance = base.read_admission(repo, entry, STATUS)
    retained_blob = subprocess.check_output(["git", "show", LAUNCH + ":" + RETAINED], cwd=repo)
    retained_pin = record(repo / RETAINED)
    if hashlib.sha256(retained_blob).hexdigest() != retained_pin["sha256"]:
        raise ValueError("Retained comparator diagnostics differ from committed evidence")
    pins += [retained_pin, report["resource_table"]]
    provenance.update(launch={"revision": LAUNCH, "path": base.LEDGER,
        "git_blob_sha256": hashlib.sha256(blob).hexdigest()}, retained=retained_pin)
    return report, json.loads(retained_blob), pins, provenance


def build(report, source, retained, retained_source):
    contexts = []
    for method in report["context"]:
        cells = []
        for cell in method["population_details"]:
            scores = {k: v for k, v in method["scores"].items() if k.startswith(cell["dataset"] + ": ")}
            cells.append({**cell, "scores": cell.get("scores", scores)})
        contexts.append({**method, "cells": cells})
    methods = [*contexts, *report["methods"]]
    if len({m["key"] for m in methods}) != len(methods):
        raise ValueError("Duplicate campaign diagnostic configuration")
    current, history = diagnostics.build({"methods": methods}, source, retained, retained_source)
    by_key = {m["key"]: (i, m) for i, m in enumerate(methods)}
    populations = report["retained_fas_samples"]["populations"]
    modes = [p["mode"] for p in populations]
    if modes != ["none", "henikoff", "henikoff_cluster"]:
        raise ValueError("Missing or reordered retained weighting FAS populations")
    for row in [*current, *history]:
        i, method = by_key[row["method"]]
        contextual = i < len(contexts)
        target = ["context", i] if contextual else ["methods", i - len(contexts)]
        for key in ("score_source_location", "detail_source_location", "coverage_source_location"):
            location = row.get(key)
            if not location: continue
            document = row["detail_source_document"] if key == "detail_source_location" else row["score_source_document"]
            if document == str(source) and location[:2] == ["methods", i]:
                tail = location[2:]
                if contextual and tail[:1] == ["cells"]:
                    tail = tail[2:] if key == "score_source_location" else ["population_details", *tail[1:]]
                row[key] = [*target, *tail]
        row.update(method_role=method.get("role", "contextual_comparator" if contextual else "not_recorded"),
            method_status=method["status"],
            evidence_scope="retained_context_not_reaudited" if contextual else "inherited_admitted_experiment_not_new_scoring")
        if row.get("score") is None: continue
        # Missing diagnostics point to their existing parent, not an invented field.
        document = retained if row["detail_source_document"] == str(retained_source) else report
        anchored = []
        for key in row["detail_source_location"]:
            try:
                document = document[key]
            except (KeyError, IndexError):
                break
            anchored.append(key)
        row["detail_source_location"] = anchored
        cell_index = row["cell_index"]
        cell = method["cells"][cell_index]
        details = cell.get("details", {})
        row["reference_groups"] = diagnostics.count(details.get("family_details", {}).get("refogs")
            if row["dataset"] == "OrthoBench" else details.get("reference_orthogroups")
            if row["dataset"] == "ThreeKingdoms" else None)
        if row["dataset"] != "QfO" or row["statistic"] != "conditional_similarity": continue
        endpoint = row["endpoint"]
        numeric = cell.get("functional_populations", {}).get(endpoint)
        if numeric is not None:
            field = "sample_pairs" if endpoint == "FAS" else "numeric_scored_pairs"
            count = numeric[field]
            location = [*target, "population_details" if contextual else "cells", cell_index,
                        "functional_populations", endpoint, field]
        elif not contextual and method["key"].startswith("profile_weighting:") and endpoint == "FAS":
            mode = cell["label"]
            index = modes.index(mode)
            population = populations[index]
            if population["FAS_reported_eligible_pool"] != row["assessed_relations"]:
                raise ValueError("Weighting FAS eligible pool differs from admitted diagnostics")
            count = population["FAS_numeric_sample_pairs"]
            location = ["retained_fas_samples", "populations", index, "FAS_numeric_sample_pairs"]
        else:
            continue
        row["numeric_scored_pairs"] = diagnostics.count(count)
        if row["assessed_relations"] is None or row["numeric_scored_pairs"] > row["assessed_relations"]:
            raise ValueError("Numeric scored population exceeds the eligible assessed pool")
        row.update(numeric_population_source_document=str(source), numeric_population_source_location=location)
    return current, history


def run(repo, output):
    if output.exists(): raise FileExistsError("Preserve every endpoint diagnostic export")
    revision, source_files = campaign.reporter_inventory(Path(__file__).resolve().parents[1])
    report, retained, checked, provenance = admit(repo)
    checked += source_files
    source = provenance["reader"]["path"]
    current, history = build(report, source, retained, provenance["retained"]["path"])
    checked = base.unique_pins(checked)
    for pin in checked: check(pin)
    output.mkdir(parents=True)
    tables = []
    for name, rows in (("endpoint_diagnostics.tsv", current), ("historical_cell_diagnostics.tsv", history)):
        path = output / name
        with path.open("x", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=COLUMNS, delimiter="\t", lineterminator="\n")
            writer.writeheader()
            writer.writerows({k: "NA" if row.get(k) is None else json.dumps(row[k]) if isinstance(row[k], list)
                else row[k] for k in COLUMNS} for row in rows)
        tables.append(record(path))
    result = dict(status="admitted_hmm_campaign_endpoint_diagnostics_exported", inputs=provenance,
        reporter_revision=revision, reporter_files=source_files, current=current, history=history, tables=tables,
        checked=checked, new_scoring_or_sampling=False, defaults_promoted=False, goal_completed=False,
        limitations=["Inherited admitted results, not new raw scoring, independent biological confirmation or a fresh historical provenance audit.",
            "Scores and precision/recall are 0-to-1. OrthoBench weighted group F1, QfO pair F1 and Three Kingdoms reference co-membership including within-species pairs are distinct statistics.",
            "Native universe coverage, benchmark-reported assignment coverage, reference coverage, submitted/retained pairs and assessed relations remain distinct; missing fields stay NA.",
            "Functional similarities have no orthology precision/recall. FAS eligible pools are not numeric sample sizes; saved samples are unseeded with repeated proteins/families. Unknown numeric counts stay NA.",
            "Configuration/control/history rows are not distinct tools. Earlier failed executions and stopped unrun endpoints stay NA, not zero; original fixed gate decisions are unchanged.",
            "Retained comparator limitations and source semantics remain in the source report. No heterogeneous mean, universal ranking, significance claim, default promotion or goal closure."])
    for pin in [*checked, *tables]: check(pin)
    write_json(output / "summary.json", result)
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    print(run(args.repo.resolve(), args.output.resolve())["status"])
