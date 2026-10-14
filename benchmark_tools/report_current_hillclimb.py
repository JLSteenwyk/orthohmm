"""Extend admitted comparison tables with current full-search/canonical cells."""

import argparse
from copy import deepcopy
import csv
import hashlib
import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.metric_frontier import ENDPOINTS
from benchmark_tools.prepare_ob_candidate_neighborhood import check, record
from benchmark_tools.report_hillclimb import audit_cell, endpoint_vector, readback_scores, replay_native, screen, snapshot
from benchmark_tools.run_membership_hillclimb import DATASETS, PRODUCTS, write_json


def upsert_cell(methods, key, label, cell):
    method = next((method for method in methods if method["key"] == key), None)
    if method is None:
        method = {"key": key, "label": label, "scores": {}, "cells": []}
        methods.append(method)
    old = [row for row in method["cells"] if row["dataset"] == cell["dataset"]]
    verified = [row for row in old if row["status"] == "independent_content_readback_verified"]
    if verified and cell["status"] != "independent_content_readback_verified":
        if len(verified) != 1:
            raise ValueError("Ambiguous existing verified cell")
        method["cells"].append({**cell, "scope": "Later execution state does not erase earlier verified results."})
        return
    if any(row["scores"] != cell["scores"] for row in verified):
        raise ValueError("Conflicting scores for an already verified configuration")
    for previous in verified:
        if "products" not in previous:
            continue
        for name in PRODUCTS:
            old = previous["products"].get(name)
            new = cell.get("products", {}).get(name)
            if old is None or new is None or (old["bytes"], old["sha256"]) != (new["bytes"], new["sha256"]):
                raise ValueError("Conflicting or missing native products despite equal scores: " + name)
    method["cells"] = [row for row in method["cells"] if row["dataset"] != cell["dataset"]]
    method["cells"].append(cell)
    if cell["status"] == "independent_content_readback_verified":
        method["scores"].update(cell["scores"])
    method["complete_eight_endpoints"] = set(method["scores"]) == set(ENDPOINTS)


def canonical_cell(repo, root, row, report):
    dataset, label = row["dataset"], row["label"]
    folder = root / dataset / label
    preflight = json.loads((folder / "preflight.json").read_text())
    native = replay_native(row, {"command": preflight["native"]})
    nonce = hashlib.sha256(str(root).encode()).hexdigest()[:16]
    scores, details, files = readback_scores(repo, dataset, row, native, folder,
                                            f"ohmm_hillclimb_r1_canon_{nonce}_{label}")
    checked = [record(folder / "preflight.json"), *native["checked"], *files]
    for pin in checked: check(pin)
    return {"status": "independent_content_readback_verified", "scores": endpoint_vector(dataset, scores),
            "details": details, "products": native["products"], "checked": checked,
            "coverage": {name: native[name] for name in ("genes", "species", "groups", "pairs")},
            "resources": {name: row[name] for name in ("graph_execution", "native_execution", "candidate_resources",
                          "total_incremental_wall_s") if name in row},
            "scope": "Complete downstream accepted-hit replay with selective tree-cache reuse; not fresh full-pipeline timing.",
            "source_revision": report["source_revision"]}


def build(repo, parent_path, full_root, canonical_root, output):
    if output.exists():
        raise FileExistsError(output)
    parent = json.loads(parent_path.read_text())
    check(parent["score_table"])
    methods = deepcopy(parent["methods"])
    output.mkdir(parents=True)
    captures = output / "snapshots"
    captures.mkdir()
    for dataset in DATASETS:
        path = full_root / dataset / "results.json"
        if not path.exists(): continue
        report, capture = snapshot(path, captures, f"full_search_{dataset}.json")
        check(report["execution_plan"])
        plan = json.loads(Path(report["execution_plan"]["path"]).read_text())
        for row in report["arms"]:
            cell = {"dataset": dataset, "status": row["status"], "producer_status": row["status"],
                    "scores": {}, "source_snapshot": capture["snapshot"], "producer_job_id": report["job_id"]}
            cell.update({name: row[name] for name in ("error", "error_type") if name in row})
            if row["status"] == "scored_development_result":
                cell.update(audit_cell(repo, "full_search", dataset, row, plan))
            upsert_cell(methods, "full_search:" + row["label"], "full_search:" + row["label"], cell)
    path = canonical_root / "results.json"
    report, capture = snapshot(path, captures, "canonical_collection.json")
    check(report["execution_plan"])
    plan = json.loads(Path(report["execution_plan"]["path"]).read_text())
    for pin in [*plan["inputs"], *plan["sources"], *plan["runtime"]["files"]]: check(pin)
    for label in ("arrival", "canonical"):
        for dataset in DATASETS:
            rows = [row for row in report["arms"] if row["dataset"] == dataset and row["label"] == label]
            if len(rows) > 1: raise ValueError("Duplicate canonical producer cell")
            cell = {"dataset": dataset, "scores": {}, "status": "not_yet_reported_execution_state_not_inferred",
                    "source_snapshot": capture["snapshot"], "producer_job_id": report["job_id"]}
            if rows:
                row = rows[0]
                cell.update(status=row["status"], producer_status=row["status"])
                cell.update({name: row[name] for name in ("error", "error_type") if name in row})
                if row["status"] == "independent_content_readback_verified":
                    cell.update(canonical_cell(repo, canonical_root, row, report))
            upsert_cell(methods, "canonical_collection:" + label,
                        "OrthoHMM " + ("canonical collection (experimental)" if label == "canonical" else "arrival replay control"), cell)
    for old in parent["methods"]:
        if not old["key"].startswith("full_search:"):
            new = next(method for method in methods if method["key"] == old["key"])
            if old != new: raise ValueError("Unrelated method or scores changed")
    with (output / "scores.tsv").open("x", newline="") as stream:
        writer = csv.writer(stream, delimiter="\t", lineterminator="\n")
        writer.writerow(["method", *ENDPOINTS, "complete_eight_endpoints"])
        for method in methods:
            writer.writerow([method["label"], *(method["scores"].get(name, "NA") for name in ENDPOINTS),
                             set(method["scores"]) == set(ENDPOINTS)])
    vectors = {method["key"].split(":", 1)[1]: method["scores"] for method in methods
               if method["key"].startswith("canonical_collection:")}
    full_vectors = {method["key"].split(":", 1)[1]: method["scores"] for method in methods
                    if method["key"].startswith("full_search:")}
    result = {"status": "verified_current_partial_vectors_no_default_promotion", "parent_readback": record(parent_path),
              "reporter_source": record(__file__), "methods": methods, "score_table": record(output / "scores.tsv"),
              "round_screening": {**parent["round_screening"], "full_search": screen(full_vectors, "control100"),
                                  "canonical_collection": screen(vectors, "arrival")},
              "limitations": [*parent["limitations"],
                  "Canonical OrthoBench F1 equals its arrival control, but native output file hashes differ; numerical equality is not assumed to mean identical predictions.",
                  "Original fresh QfO control crashed during native phylogeny; its failed status and unrun cap400 remain visible. Cached native recovery is separate and unscored until independently verified.",
                  "Missing endpoint values are NA, not zero, inherited scores or inferred success. Fresh and cached replay costs must not be pooled."],
              "defaults_promoted": False, "goal_completed": False}
    write_json(output / "summary.json", result)
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("repo", "parent", "full-root", "canonical-root", "output"):
        parser.add_argument("--" + name, type=Path, required=True)
    args = parser.parse_args()
    result = build(args.repo.resolve(), args.parent.resolve(), args.full_root.resolve(), args.canonical_root.resolve(), args.output.resolve())
    print(result["status"], len(result["methods"]))
