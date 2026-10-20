"""Consolidate admitted individual effects and a stopped interaction panel."""

import argparse
from copy import deepcopy
import csv
import hashlib
import json
from pathlib import Path
import subprocess
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from benchmark_tools import export_hmm_campaign_diagnostics as diagnostics
from benchmark_tools import plan_weighted_projection_followup as planner
from benchmark_tools import report_hmm_evidence_campaign as base
from benchmark_tools import report_hmm_transfer_campaign as campaign
from benchmark_tools.prepare_ob_candidate_neighborhood import check, record
from benchmark_tools.run_membership_hillclimb import write_json
from benchmark_tools.weighted_projection_experiment import configuration, point, specification


FREEZE = "benchmark_tools/hmm_campaign_candidate_freeze_20261020.json"
FREEZE_REVISION = "9a1f10ca8a403eb2c8811152167c0b255e67b0ef"
REVIEW = "benchmark_tools/results/hmm_confirmation_review_20261020.json"
CANDIDATES = ("profile_weighting:henikoff", "profile_weighting:henikoff_cluster",
              "projection_recovery:projection_off")
INDIVIDUAL_KEYS = (
    "profile_weighting:none", *CANDIDATES[:2], "profile_weighting:species_henikoff_cluster",
    "alignment:strict_span_guard", "local:off", "local:local_joint_v1",
    "projection_recovery:legacy_off", CANDIDATES[-1],
    "projection_recovery:projection_builtin_reselect", "projection_recovery:projection_hmmer_score_only",
    "projection_recovery:projection_hmmer_reselect", "selective_original:failed_panel",
    "projection_original:failed_panel")


def validate_review(freeze, review):
    if (freeze.get("status") != "frozen_experimental_portfolio_before_postfreeze_confirmation_review"
            or tuple(c["key"] for c in freeze["retained_candidates"]) != CANDIDATES
            or review.get("status") != "bounded_postfreeze_confirmation_review_complete_no_independent_panel_admitted"
            or review.get("freeze_revision") != FREEZE_REVISION
            or review.get("logical_probes") != 6
            or len(review["requests"]) != 6
            or [r["url"] for r in review["requests"]] != freeze["postfreeze_confirmation_review"]["head_urls"]
            or any(r.get("method") != "HEAD" or r.get("body_bytes_read") != 0 for r in review["requests"])
            or review.get("request_policy") != freeze["postfreeze_confirmation_review"]["request_bounds"]
            or review.get("limitations_nonblocking") is not True):
        raise ValueError("Require the frozen portfolio and bounded label-free post-freeze review")
    policies = [(c["profile_sequence_weighting"], c["profile_projection_policy"], c["profile_alignment_policy"])
                for c in freeze["retained_candidates"]]
    if policies != [("henikoff", "legacy", "off"), ("henikoff_cluster", "legacy", "off"),
                    ("none", "residue_preserving", "off")]:
        raise ValueError("Frozen standalone policies changed")
    for document in (freeze, review):
        for key in ("defaults_promoted", "goal_completed", "independent_confirmation_established"):
            if document.get(key) is not False:
                raise ValueError("Confirmation limits cannot become promotion or completion")
    if (review.get("new_accuracy_evaluation") is not False
            or review.get("new_reference_assignments_read") is not False
            or review.get("new_reference_assignments_downloaded") is not False):
        raise ValueError("Metadata review cannot contain new label access or accuracy evaluation")


def assemble(individual, pilot, producer, freeze, review):
    validate_review(freeze, review)
    if (individual.get("status") != diagnostics.STATUS or individual.get("accuracy_result_admitted") is not True
            or individual.get("defaults_promoted") is not False or individual.get("goal_completed") is not False
            or tuple(m["key"] for m in individual["methods"]) != INDIVIDUAL_KEYS
            or individual["context"] != pilot["context"]):
        raise ValueError("Require the complete individual-effects admission and unchanged comparator context")
    selected = planner.selection(pilot, producer)
    if selected["execution_arms"]:
        raise ValueError("Passing interactions require matched full-vector follow-up, not plateau closure")
    if (pilot["uncertainty"]["multiplicity_endpoints"] != 18
            or [[c["candidate"], c["control"]] for c in pilot["uncertainty"]["contrasts"]]
            != specification()["uncertainty"]["contrasts"]):
        raise ValueError("The original six interaction uncertainty contrasts must remain intact")
    methods = deepcopy(individual["methods"])
    by_key = {m["key"]: m for m in methods}
    if tuple(m["key"] for m in methods if m.get("full_vector_screen", {}).get("eligible")) != CANDIDATES:
        raise ValueError("Frozen candidates differ from the original full-vector screens")
    if any(m["status"] in ("full_vector_extension_required", "applicability_fallback_required") for m in methods):
        raise ValueError("Unfinished individual mechanisms cannot be closed by an interaction report")
    canonical, = [c for c in individual["context"] if c["key"] == freeze["baseline"]["context_alias"]]
    if by_key[freeze["baseline"]["key"]]["scores"] != canonical["scores"]:
        raise ValueError("Full-vector canonical control changed")
    for key in (freeze["baseline"]["key"], *CANDIDATES):
        if set(by_key[key]["scores"]) != set(base.ENDPOINTS) or by_key[key]["missing_endpoints"]:
            raise ValueError("Every frozen standalone configuration requires all eight endpoints")
    for cell in pilot["cells"]:
        arm = configuration(cell["label"])
        if arm["role"] != "candidate":
            original = by_key[arm["reference_key"]]["cells"][0]
            if point(cell) != point(original):
                raise ValueError("Interaction control differs from its admitted standalone constituent")
        methods.append({"key": "weighted_projection:" + cell["label"],
            "label": "OrthoHMM weighting/projection " + cell["label"], "role": arm["role"],
            "status": "stopped_after_stronger_constituent_point_gate" if arm["role"] == "candidate"
                      else "interaction_control_reproduced",
            "scores": deepcopy(cell["scores"]), "cells": [deepcopy(cell)],
            "missing_endpoints": [k for k in base.ENDPOINTS if k not in cell["scores"]],
            "matched_point_gate": deepcopy(pilot["matched_gates"].get(cell["label"]))})
    if len({m["key"] for m in [*individual["context"], *methods]}) != len(individual["context"]) + len(methods):
        raise ValueError("Duplicate contextual or experimental configuration")
    result = {**individual, "methods": methods, "context": deepcopy(individual["context"]),
        "candidate_freeze": deepcopy(freeze), "confirmation_review": deepcopy(review),
        "interaction": {"selection": selected, "uncertainty": deepcopy(pilot["uncertainty"]),
                        "matched_gates": deepcopy(pilot["matched_gates"])},
        "recommendation": freeze["final_recommendation_scope"],
        "required_next_actions": [{"method": "campaign", "action":
            "Finish final documentation, requirement-by-requirement audit and verified remote delivery. No further inference or retuning is required by the admitted stopping rules."}],
        "scope": "Completed bounded development investigation, not independent biological confirmation or universal superiority. Original individual effects, failures and stopping gates are preserved; both fixed interactions fail the stronger-constituent gate. Functional means remain conditional, TreeFam supporting and Three Kingdoms supplementary. No mixed-endpoint average or default promotion."}
    result.pop("checked", None)
    result["frontier"] = base.metric_frontier({m["key"]: m["scores"] for m in [*result["context"], *methods]})
    return result


def verify_diagnostic_pointers(report, retained, rows, source, retained_source):
    documents = {str(source): report, str(retained_source): retained}
    def resolve(document, location):
        value = documents[document]
        for key in location:
            value = value[key]
        return value
    for row in rows:
        if row.get("score") is None:
            continue
        if resolve(row["score_source_document"], row["score_source_location"]) != row["score"]:
            raise ValueError("Diagnostic score pointer differs from the admitted value")
        if not isinstance(resolve(row["detail_source_document"], row["detail_source_location"]), dict):
            raise ValueError("Diagnostic detail pointer has no admitted parent")
        if row.get("coverage_source_location") and not isinstance(resolve(
                row["score_source_document"], row["coverage_source_location"]), dict):
            raise ValueError("Diagnostic coverage pointer has no admitted parent")
        if row.get("numeric_population_source_location") and resolve(
                row["numeric_population_source_document"], row["numeric_population_source_location"]) != row["numeric_scored_pairs"]:
            raise ValueError("Numeric population pointer differs from the admitted count")


def run(repo, output):
    if output.exists():
        raise FileExistsError("Preserve every final campaign compilation attempt")
    driver = Path(__file__).resolve().parents[1]
    revision = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=driver, text=True).strip()
    if driver == repo:
        raise ValueError("Require a frozen committed report checkout")
    files = base.sources(driver) + [record(driver / p) for p in (
        FREEZE, REVIEW, "benchmark_tools/report_hmm_campaign_plateau.slurm")]
    for pin in files:
        relative = Path(pin["path"]).relative_to(driver)
        if relative.suffix != ".so" and hashlib.sha256(subprocess.check_output(
                ["git", "show", revision + ":" + str(relative)], cwd=driver)).hexdigest() != pin["sha256"]:
            raise ValueError("Final report source differs from committed bytes")
    freeze_pin, review_pin = record(driver / FREEZE), record(driver / REVIEW)
    freeze, review = json.loads((driver / FREEZE).read_text()), json.loads((driver / REVIEW).read_text())
    original_freeze = subprocess.check_output(["git", "show", FREEZE_REVISION + ":" + FREEZE], cwd=repo)
    if ((driver / FREEZE).read_bytes() != original_freeze
            or any(review["freeze"][k] != freeze_pin[k] for k in ("bytes", "sha256"))):
        raise ValueError("Final candidate freeze changed after confirmation review")
    individual, retained, checked, individual_provenance = diagnostics.admit(repo)
    pilot, plan, selected, pilot_pins, pilot_provenance = planner.authenticate(repo)
    producer = json.loads(Path(pilot["producer"]["path"]).read_text())
    result = assemble(individual, pilot, producer, freeze, review)
    if result["interaction"]["selection"] != selected:
        raise ValueError("Authenticated interaction selection changed during collation")
    checked = base.unique_pins([*checked, *pilot_pins, *files, *review["previous_evidence"]])
    for pin in checked:
        check(pin)
    output.mkdir(parents=True)
    result.update(status="running", reporter_revision=revision, reporter_files=files,
        inputs={"individual": individual_provenance, "interaction": pilot_provenance},
        frozen_candidate_record=freeze_pin, confirmation_record=review_pin,
        accuracy_result_admitted=False, defaults_promoted=False, goal_completed=False)
    try:
        base.write_table(output / "scores.tsv", result["methods"], result["context"])
        costs, cost_pins = campaign.resource_table(output, result["methods"])
        current, history = diagnostics.build(result, str(output / "summary.json"), retained,
                                              individual_provenance["retained"]["path"])
        verify_diagnostic_pointers(result, retained, [*current, *history], output / "summary.json",
                                   individual_provenance["retained"]["path"])
        for name, rows in (("endpoint_diagnostics.tsv", current), ("historical_cell_diagnostics.tsv", history)):
            with (output / name).open("x", newline="") as stream:
                writer = csv.DictWriter(stream, fieldnames=diagnostics.COLUMNS, delimiter="\t", lineterminator="\n")
                writer.writeheader()
                writer.writerows({k: "NA" if row.get(k) is None else json.dumps(row[k]) if isinstance(row[k], list)
                                  else row[k] for k in diagnostics.COLUMNS} for row in rows)
        tables = [record(output / name) for name in (
            "scores.tsv", "resources.tsv", "endpoint_diagnostics.tsv", "historical_cell_diagnostics.tsv")]
        result.update(tables=tables, score_table=tables[0], resource_table=tables[1],
            resource_rows=costs, endpoint_diagnostics=current, historical_diagnostics=history,
            configurations=len(result["context"]) + len(result["methods"]), endpoint_slots=len(current),
            scored_endpoint_observations=sum(r.get("score") is not None for r in current),
            resource_scope="Shared-host cached incremental phase costs; compatible cache reuse differs. Candidate/arm metadata not independently retimed; local nomination/joint-tree preparation excluded. No full-tool timing leaderboard, summed RSS or summation of arm totals with components.")
        result["checked"] = base.unique_pins([*checked, *cost_pins, *tables])
        for pin in result["checked"]:
            check(pin)
        result.update(status="complete_admitted_hmm_campaign_plateau_report", accuracy_result_admitted=True)
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
