"""Collate separately admitted HMM experiments without changing their gates."""

import argparse
from copy import deepcopy
import hashlib
import json
from pathlib import Path
import subprocess
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from benchmark_tools.metric_frontier import ENDPOINTS, metric_frontier
from benchmark_tools.prepare_ob_candidate_neighborhood import check, record
from benchmark_tools.report_profile_evidence_pilot import verify_selective_panel
from benchmark_tools.report_weighting_followup import assemble, write_table
from benchmark_tools.run_membership_hillclimb import write_json
from benchmark_tools.run_profile_weighting_followup import successful_job
from benchmark_tools.run_profile_weighting_pilot import advancement
from benchmark_tools.run_three_kingdoms_membership_supplement import sources


LAUNCH = "23abc652b56990063552572ab98549aeb4a2269e"
LEDGER = "benchmark_tools/results/hmm_evidence_campaign_status_20261018.json"
READERS = {
    "weighting": ("profile_weighting_postrun_readback", "weighting_full_vectors_independently_read_back"),
    "local": ("local_joint_postrun_readback", "local_joint_pilot_independently_read_back"),
    "hmmer": ("selective_hmmer_postrun_readback", "profile_pilot_independently_read_back"),
}


def launch(repo):
    blob = subprocess.check_output(["git", "show", LAUNCH + ":" + LEDGER], cwd=repo)
    ledger = json.loads(blob)
    return {key: ledger[name] for key, (name, _) in READERS.items()}, {
        "revision": LAUNCH, "path": LEDGER, "git_blob_sha256": hashlib.sha256(blob).hexdigest()}


def unique_pins(pins):
    result = {}
    for pin in pins:
        path = str(Path(pin["path"]).resolve())
        normalized = {"path": path, "bytes": pin["bytes"], "sha256": pin["sha256"]}
        if path in result and result[path] != normalized:
            raise ValueError("Conflicting original content records for the same path")
        result[path] = normalized
    return list(result.values())


def read_admission(repo, entry, expected_status):
    accounting = successful_job(entry["job"]["job_id"], 1, "16G")
    source = repo / entry["scientific_reader_source"]
    if subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=source, text=True).strip() != entry["runner_commit"]:
        raise ValueError("Campaign input reader source moved")
    path = repo / entry["job"]["summary"]
    pin = record(path)
    report = json.loads(path.read_text())
    if (report.get("status") != expected_status or report.get("accuracy_result_admitted") is not True
            or report.get("reporter_revision") != entry["runner_commit"]
            or report.get("defaults_promoted") is not False or report.get("goal_completed") is not False
            or not report.get("checked") or not report.get("reporter_files")):
        raise ValueError("Require a complete separately admitted source report, not partial or failed science")
    expected_files = sources(source)
    reported = {p["path"]: p for p in report["reporter_files"]}
    if any(reported.get(p["path"]) != p for p in expected_files):
        raise ValueError("Input report source inventory differs from frozen reader")
    pins = [pin, report["score_table"], *report["checked"], *report["reporter_files"]]
    return report, pins, {"reader": pin, "scheduler": accounting, "revision": entry["runner_commit"]}


def point(cell):
    return {"f1": cell["scores"][ENDPOINTS[0]], "precision": cell["details"]["precision"],
            "recall": cell["details"]["recall"]}


def pilot_method(prefix, cell, role, control):
    label = cell["label"]
    gate = advancement(point(control), point(cell))
    applicable = cell.get("activation", {}).get("completely_measured_queries", 1) > 0
    if role == "candidate":
        if not applicable:
            gate.update(eligible=False, reason="No completely measured queries; applicability fallback remains required")
        if gate != cell["advancement"]:
            raise ValueError("Campaign point/applicability gate differs from admitted pilot")
        status = ("applicability_fallback_required" if not applicable else
                  "full_vector_extension_required" if gate["eligible"] else "stopped_after_orthobench_point_gate")
    else:
        status = role
    return {"key": prefix + ":" + label, "label": "OrthoHMM " + prefix + " " + label,
        "role": role, "status": status, "scores": cell["scores"], "cells": [cell],
        "point_gate": gate if role == "candidate" else None,
        "missing_endpoints": [key for key in ENDPOINTS if key not in cell["scores"]]}


def collate(weighting, local, hmmer, hmmer_producer, *, recovery=False):
    context = weighting["context"]
    if local["context"] != context or hmmer["context"] != context:
        raise ValueError("The three experiment reports use different retained comparator snapshots")
    canonical, = [m for m in context if m["key"] == "canonical_collection:canonical"]
    methods = deepcopy(weighting["methods"])
    expected = ["profile_weighting:" + label for label in ("none", "henikoff", "henikoff_cluster", "species_henikoff_cluster")]
    if [m["key"] for m in methods] != [*expected, "alignment:strict_span_guard"]:
        raise ValueError("Require all weighting and stopped alignment rows, not selected winners")
    pilot = {"cells": [m["cells"][0] for m in methods[:4]]}
    if assemble(pilot, weighting["cells"]) != methods[:4]:
        raise ValueError("Complete weighting vectors or original full-vector screens differ")
    if methods[0]["scores"] != canonical["scores"]:
        raise ValueError("Weighting none control differs from full canonical comparator vector")
    for method in methods[1:3]:
        method["campaign_decision"] = ("development_screen_pass_confirmation_required" if method["full_vector_screen"]["eligible"]
                                       else "fixed_full_vector_branch_stopped_without_retuning")
    control = methods[0]["cells"][0]
    if [c["label"] for c in local["cells"]] != ["off", "local_joint_v1"]:
        raise ValueError("Require both ordered local pilot cells")
    if recovery:
        from benchmark_tools.report_projection_recovery import ARMS, verify_panel
        labels, verify = list(ARMS), verify_panel
    else:
        labels, verify = ["off", "builtin_reselect", "hmmer_score_only", "hmmer_reselect"], verify_selective_panel
    if [c["label"] for c in hmmer["cells"]] != labels:
        raise ValueError("Require the complete ordered HMMER panel")
    if any(point(report["cells"][0]) != point(control) for report in (local, hmmer)):
        raise ValueError("Local/HMMER off precision, recall or F1 differs from the weighting canonical control")
    if hmmer["backend_contrast"] != verify(hmmer["cells"], hmmer_producer["arms"]):
        raise ValueError("HMMER-specific contrast differs from matched builtin decision control")
    hmmer_prefix = "projection_recovery" if recovery else "selective"
    for prefix, report in (("local", local), (hmmer_prefix, hmmer)):
        for cell in report["cells"]:
            if cell["status"] != "independent_content_readback_verified" or cell["dataset"] != "OrthoBench":
                raise ValueError("Unadmitted pilot cell in campaign report")
            role = ("unchanged_control" if cell["label"] == report["cells"][0]["label"] else
                    "measurement_control" if cell["label"] in ("hmmer_score_only", "projection_hmmer_score_only") else "candidate")
            method = pilot_method(prefix, cell, role, report["cells"][0])
            if recovery and prefix == hmmer_prefix:
                method["construction_role"] = ("legacy_control" if cell["label"] == "legacy_off" else
                                               "projection_only" if cell["label"] == "projection_off" else "matched_repaired_profile")
                if "matched_projection_advancement" in cell:
                    method["matched_projection_gate"] = cell["matched_projection_advancement"]
                method["interpretation"] = "Original advancement uses legacy canonical; repair-only gain is not backend-specific benefit."
            methods.append(method)
    vectors = {m["key"]: m["scores"] for m in [*context, *methods]}
    if len(vectors) != len(context) + len(methods): raise ValueError("Duplicate campaign method key")
    frontier = metric_frontier(vectors)
    actions = []
    for method in methods:
        if method["status"] == "full_vector_extension_required":
            actions.append({"method": method["key"], "action": "Run unchanged corrected QfO and Three Kingdoms, then complete tested production opt-in interface if promising; no default promotion from an OrthoBench pilot."})
        elif method["status"] == "applicability_fallback_required":
            actions.append({"method": method["key"], "action": "Probe corrected QfO then Three Kingdoms in the original fixed order without retuning. Tool failures are not inactivity; unmeasured queries are not negative accuracy evidence."})
        elif method.get("full_vector_screen", {}).get("eligible"):
            actions.append({"method": method["key"], "action": "Candidate passes original development full-vector point screen. Freeze before seeking independent confirmation; inspect conditional functional populations, uncertainty and costs. This does not establish significance or justify default promotion."})
    population_changes = []
    qfo = [c for c in weighting["cells"] if c["dataset"] == "QfO"]
    none, = [c for c in qfo if c["label"] == "none"]
    for cell in qfo:
        if cell["label"] == "none": continue
        for name in ("EC", "GO", "FAS"):
            base, actual = none["details"][name], cell["details"][name]
            population_changes.append({"method": "profile_weighting:" + cell["label"], "endpoint": "QfO: " + name,
                "control_assessed_relations": base["assessed_relations"], "candidate_assessed_relations": actual["assessed_relations"],
                "control_mean": base["similarity"], "candidate_mean": actual["similarity"],
                "scope": "Conditional similarities over assessed relation populations, not orthology F1 or matched-population causal effects."})
    return {"methods": methods, "context": context, "frontier": frontier, "required_next_actions": actions,
        "functional_population_comparisons": population_changes,
        "uncertainty": {"weighting": weighting["orthobench_uncertainty"], "local": local["uncertainty"],
                        "projection_recovery" if recovery else "selective": hmmer["uncertainty"]},
        "backend_contrast": hmmer["backend_contrast"],
        "recommendation": "Keep production defaults unchanged. Preserve fixed negative branches and experimental candidates. Execute listed follow-ups/confirmation or document unavailable independent data before a final recommendation. Do not combine mechanisms before their individual full-vector effects are understood.",
        "scope": "Development-exposed point comparisons and inherited comparator context; no mixed-endpoint average, universal superiority, independent confirmation or timing leaderboard. TreeFam supports primary OrthoBench/QfO evidence; Three Kingdoms remains supplementary. FAS is unseeded. Shared-host cached incremental time and max-process RSS are not isolated fresh pipeline costs or simultaneous aggregate memory."}


def run(repo, output):
    if output.exists(): raise FileExistsError("Preserve every campaign report namespace")
    entries, binding = launch(repo)
    reports, pins, provenance = {}, [], {}
    for key, entry in entries.items():
        report, checked, receipt = read_admission(repo, entry, READERS[key][1])
        reports[key], provenance[key] = report, receipt
        pins += checked
    hmmer_source = reports["hmmer"]["producer"]
    pins.append(hmmer_source)
    for pin in unique_pins(pins): check(pin)
    producer = json.loads(Path(hmmer_source["path"]).read_text())
    result = collate(reports["weighting"], reports["local"], reports["hmmer"], producer)
    reader = Path(__file__).resolve().parents[1]
    revision = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=reader, text=True).strip()
    source_pins = sources(reader)
    for pin in source_pins:
        relative = Path(pin["path"]).relative_to(reader)
        if relative.suffix != ".so" and hashlib.sha256(subprocess.check_output(
                ["git", "show", revision + ":" + str(relative)], cwd=reader)).hexdigest() != pin["sha256"]:
            raise ValueError("Campaign reporter differs from committed source")
    pins += source_pins
    output.mkdir(parents=True)
    write_table(output / "scores.tsv", result["methods"], result["context"])
    result.update(status="admitted_campaign_results_collated_followup_decisions_explicit", launch_ledger=binding,
        inputs=provenance, reporter_revision=revision, reporter_files=source_pins, score_table=record(output / "scores.tsv"),
        defaults_promoted=False, goal_completed=False)
    result["checked"] = unique_pins(pins)
    for pin in result["checked"]: check(pin)
    write_json(output / "summary.json", result)
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    print(run(args.repo.resolve(), args.output.resolve())["status"])
