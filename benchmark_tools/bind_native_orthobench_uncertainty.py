"""Reuse retained intervals only for identical native RefOG sufficient statistics."""

import argparse
import hashlib
import json
from pathlib import Path

from benchmark_tools.export_native_factorial_progress import collect, load, require


FACTORIAL_SHA = "6a0d588b5cb47c60fc6bc8bae8aa0c83e5f2aadb11de970919d8c6527c387141"
SCORE_SCHEMAS = {"native_factorial_orthobench_score_v1",
                 "native_factorial_receipt_failure_recovery_v1"}


def source_record(path):
    source = Path(path).resolve()
    data = source.read_bytes()
    return {"path": str(source), "bytes": len(data),
            "sha256": hashlib.sha256(data).hexdigest()}


def named(records):
    names = [row["refog"] for row in records]
    require(len(names) == 70 and len(set(names)) == 70, "Require 70 unique RefOGs")
    return {row["refog"]: row for row in records}


def project(scores, factorial):
    require(scores and set(scores) <= set(factorial["scores"]), "Unknown or absent native cells")
    require(factorial["planned_contrasts"] == 12
            and factorial["multiplicity_endpoints"] == 36, "Changed multiplicity scope")
    for cell, score in scores.items():
        require(named(score["refog_records"]) == named(factorial["scores"][cell]["refog_records"]),
                "Native family records differ: " + cell)
        require(all(score[k] == factorial["scores"][cell][k]
                    for k in ("f_score", "precision", "recall")), "Native point estimates differ")
    contrasts = []
    for row in factorial["comparisons"]:
        missing = [cell for cell in (row["on"], row["off"]) if cell not in scores]
        if missing:
            contrasts.append({**{k: row[k] for k in ("factor", "on", "off", "fixed")},
                              "status": "native_records_unavailable", "missing_cells": missing,
                              "metrics": None})
        else:
            require(row["status"] == "complete" and row["metrics"] is not None,
                    "Retained contrast is unavailable")
            contrasts.append({**row, "status": "native_records_matched", "missing_cells": []})
    require(len(contrasts) == 12, "Wrong contrast count")
    return contrasts


def bind(snapshot_path, snapshot_sha, factorial_path):
    evidence = []
    snapshot, snapshot_ref = load(snapshot_path, snapshot_sha, evidence)
    factorial, factorial_ref = load(factorial_path, FACTORIAL_SHA, evidence)
    exporter = source_record(collect.__code__.co_filename)
    require(exporter == snapshot["source"], "Snapshot exporter source differs")
    documents = []
    for ref in snapshot["evidence"]:
        data, observed = load(ref["path"], ref["sha256"], evidence)
        require(observed == ref, "Snapshot evidence identity differs")
        documents.append((data, ref))
    attempts, native, bindings = [], {}, {}
    for review, ref in documents:
        if review.get("schema") != "native_factorial_terminal_review_v1":
            continue
        matches = [(data, item) for data, item in documents
                   if data.get("schema") in SCORE_SCHEMAS and data["index"] == review["index"]]
        require(len(matches) <= 1, "Duplicate score binding")
        score, score_ref = matches[0] if matches else (None, None)
        attempts.append([ref["path"], ref["sha256"],
                         score_ref["path"] if score_ref else "-",
                         score_ref["sha256"] if score_ref else "-"])
        if score is None:
            continue
        cell = review["cell"]
        require(cell not in native, "Duplicate native cell")
        if score["schema"] == "native_factorial_receipt_failure_recovery_v1":
            recovered, recovered_ref = load(score["score"]["path"], score["score"]["sha256"], evidence)
            require(recovered_ref == score["score"], "Recovery score identity differs")
            native[cell] = recovered["score"]
        else:
            native[cell] = score["score_percent"]
        bindings[cell] = {"index": review["index"], "job_id": review["job_id"],
                          "terminal_review": ref, "score_or_recovery": score_ref,
                          "scheduler_state": review["scheduler_state"],
                          "scheduler_exit_code": review["scheduler_exit_code"],
                          "whole_job_succeeded": review["scheduler_state"] == "COMPLETED"}
    replay = collect(snapshot["plan"]["path"], snapshot["plan"]["sha256"], attempts)
    require(all(snapshot[k] == value for k, value in replay.items()), "Snapshot replay differs")
    result = {"schema": "native_orthobench_retained_uncertainty_binding_v1",
              "snapshot": snapshot_ref, "factorial": factorial_ref, "bound_cells": bindings,
              "contrasts": project(native, factorial), "families": factorial["families"],
              "replicates_reused": factorial["replicates"], "seed_reused": factorial["seed"],
              "alpha": factorial["alpha"], "planned_contrasts": 12, "multiplicity_endpoints": 36,
              "new_bootstrap_draws": 0, "independent_confirmation": False,
              "new_accuracy_or_resource_admission": False, "publication_ready": False,
              "limitations": [
                  "All family records, including metadata, must match; aggregate agreement alone is insufficient.",
                  "Intervals are reused from retained paired RefOG draws, not new independent evidence.",
                  "Development exposure, family exchangeability and approximate percentile-coverage limits remain.",
                  "All 36 planned endpoints remain in the multiplicity adjustment; unavailable native contrasts are not imputed.",
                  "Reconciliation changes candidate co-membership to root-HOG co-membership, not resolved-pair accuracy.",
                  "Profile-off retains initial HMM search; these effects do not measure its total contribution.",
                  "Failed-wrapper scientific recovery is retained, not converted to successful timing.",
                  "Direct report replay is not a new transitive raw review, runtime closure or speed comparison."]}
    for ref in evidence:
        _, observed = load(ref["path"], ref["sha256"], [])
        require(observed == ref, "Evidence changed during binding")
    result["evidence"] = evidence
    result["exporter_source"] = exporter
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--snapshot", required=True)
    parser.add_argument("--snapshot-sha256", required=True)
    parser.add_argument("--factorial", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    result = bind(args.snapshot, args.snapshot_sha256, args.factorial)
    result["source"] = source_record(__file__)
    with Path(args.output).open("x") as handle:
        json.dump(result, handle, indent=2, sort_keys=True)
        handle.write("\n")
    print(json.dumps({"bound_cells": len(result["bound_cells"]),
                      "matched_contrasts": sum(row["status"] == "native_records_matched"
                                               for row in result["contrasts"])}))


if __name__ == "__main__":
    main()
