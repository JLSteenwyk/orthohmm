"""Authenticate complete recovery evidence and prepare bounded transfer inputs."""

import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import sys
from types import SimpleNamespace

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from benchmark_tools import report_hmm_evidence_campaign as campaign
from benchmark_tools import report_hmm_recovery_campaign as recovery_campaign
from benchmark_tools import report_profile_evidence_pilot as pilot_reader
from benchmark_tools import report_projection_recovery as recovery_reader
from benchmark_tools.plan_singleton_tie_followup import verified_cell
from benchmark_tools.prepare_ob_candidate_neighborhood import check, record
from benchmark_tools.run_full_search_hillclimb import runtime
from benchmark_tools.run_membership_hillclimb import PRODUCTS, write_json
from benchmark_tools.run_profile_weighting_followup import selected_case, successful_job
from benchmark_tools.run_three_kingdoms_membership_supplement import sources


PROTOCOL = "hmm_projection_followup_20261019.json"
SCIENCE = "a6a56ac50ef7d89a2b70428b28914f8bcdce876d"
READER = "84ac3933a759e8cf30251c173b74b848b6512961"
DATASETS = ("QfO", "ThreeKingdoms")
CANDIDATES = ("projection_off", "projection_builtin_reselect", "projection_hmmer_reselect")


def selection(report, producer):
    cells = report["cells"]
    if (report.get("status") != "projection_recovery_independently_read_back"
            or report.get("accuracy_result_admitted") is not True or report.get("source_revision") != SCIENCE
            or report.get("reporter_revision") != READER
            or report.get("defaults_promoted") is not False or report.get("goal_completed") is not False):
        raise ValueError("Require complete independently admitted original recovery, not partial scores")
    if any(cell.get("source_revision") != SCIENCE
           or (cell["coverage"]["genes"], cell["coverage"]["species"]) != (251378, 12) for cell in cells):
        raise ValueError("Recovery source or complete pilot universe differs")
    contrast = recovery_reader.verify_panel(cells, producer["arms"])
    if contrast != report["backend_contrast"]:
        raise ValueError("Recovery backend contrast differs from complete original panel")
    by_label = {cell["label"]: cell for cell in cells}
    advancing = [label for label in CANDIDATES if by_label[label]["advancement"]["eligible"]]
    if advancing != report["advancing_labels"]:
        raise ValueError("Declared candidate selection differs from original gates")
    fallback = [label for label in CANDIDATES[1:]
                if by_label[label]["activation"]["completely_measured_queries"] == 0]
    selected = set(advancing) | set(fallback)
    needed = set()
    if selected:
        needed.update(("legacy_off", "projection_off"))
        needed.update(selected)
        if "projection_hmmer_reselect" in selected:
            needed.update(("projection_builtin_reselect", "projection_hmmer_score_only"))
    arms = [label for label in recovery_reader.ARMS if label in needed]
    roles = {}
    for label in arms:
        roles[label] = ("advancing_candidate" if label in advancing else "applicability_followup" if label in fallback else
                        "canonical_control" if label == "legacy_off" else "matched_construction_control" if label == "projection_off" else
                        "measurement_control" if label == "projection_hmmer_score_only" else "matched_decision_control")
    dependencies = {"legacy_off": [], "projection_off": ["legacy_off"],
        "projection_builtin_reselect": ["legacy_off", "projection_off"],
        "projection_hmmer_score_only": ["projection_off"],
        "projection_hmmer_reselect": ["legacy_off", "projection_off", "projection_builtin_reselect", "projection_hmmer_score_only"]}
    return {"advancing_candidates": advancing, "applicability_candidates": fallback,
        "stopped_candidates": [label for label in CANDIDATES if label not in selected],
        "execution_arms": arms, "roles": roles, "control_dependencies": {label: dependencies[label] for label in arms},
        "backend_contrast": contrast,
        "original_gates": {label: by_label[label]["advancement"] for label in CANDIDATES},
        "matched_projection_gates": {label: by_label[label]["matched_projection_advancement"] for label in CANDIDATES[1:]},
        "scope": "Original fixed pilot gates and declared inactivity only; required controls do not become advancing candidates."}


def authenticate(repo):
    entries, science, bindings = recovery_campaign.launch(repo)
    entry = entries["recovery"]
    if science["runner_commit"] != SCIENCE or entry["runner_commit"] != READER:
        raise ValueError("Wrong original recovery science or reader revision")
    report, pins, provenance = campaign.read_admission(repo, entry, "projection_recovery_independently_read_back")
    if report["launch_ledger"] != entry["launch_binding"]:
        raise ValueError("Recovery readback uses a different immutable launch")
    for pin in (report["producer"], report["execution_plan"]): check(pin)
    producer = json.loads(Path(report["producer"]["path"]).read_text())
    pilot = json.loads(Path(report["execution_plan"]["path"]).read_text())
    pilot_reader.require_completed(producer, pilot, recovery_reader.ARMS, "projection_recovery", science["inference"]["job_id"])
    accounting = successful_job(science["inference"]["job_id"], 16, "128G")
    if producer["execution_plan"] != report["execution_plan"] or pilot["source_revision"] != SCIENCE:
        raise ValueError("Recovery producer/plan source binding differs")
    source = (repo / science["scientific_source"]).resolve()
    if (Path(pilot["source_root"]) != source
            or subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=source, text=True).strip() != SCIENCE
            or sources(source) + [record(source / "benchmark_tools/run_profile_weighting_pilot.slurm")] != pilot["sources"]
            or runtime() != pilot["runtime"]):
        raise ValueError("Original recovery scientific source/runtime changed")
    chosen = selection(report, producer)
    pins += [report["producer"], report["execution_plan"], *pilot["sources"], *pilot["runtime"]["files"]]
    for pin in campaign.unique_pins(pins): check(pin)
    return report, pilot, chosen, campaign.unique_pins(pins), {
        "launches": bindings, "reader": provenance, "science_accounting": accounting}


def prepare(repo, driver, output):
    if output.exists(): raise FileExistsError("Preserve every projection follow-up preparation")
    if driver == repo or Path(__file__).resolve() != driver / "benchmark_tools/plan_projection_followup.py":
        raise ValueError("Require committed frozen preparation source")
    report, pilot, chosen, pins, provenance = authenticate(repo)
    protocol_path = driver / "benchmark_tools" / PROTOCOL
    protocol = json.loads(protocol_path.read_text())
    if (protocol["pilot_launch"] != recovery_campaign.LAUNCH or protocol["scientific_revision"] != SCIENCE
            or protocol["reader_revision"] != READER or protocol["datasets"] != list(DATASETS)
            or protocol["candidate_order"] != list(CANDIDATES)
            or protocol["resources"]["cpus_per_run"] != 16 or protocol["resources"]["memory_gib"] != 128
            or protocol["resources"]["concurrent_new_scientific_runs"] != 1):
        raise ValueError("Changed prespecified projection transfer protocol")
    revision = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=driver, text=True).strip()
    driver_pins = sources(driver) + [record(protocol_path)]
    for pin in driver_pins:
        relative = Path(pin["path"]).relative_to(driver)
        if relative.suffix != ".so" and hashlib.sha256(subprocess.check_output(
                ["git", "show", revision + ":" + str(relative)], cwd=driver)).hexdigest() != pin["sha256"]:
            raise ValueError("Projection follow-up planner differs from committed bytes")
    cases = []
    api = SimpleNamespace(record=record, check=check, verified_cell=verified_cell, PRODUCTS=PRODUCTS)
    if chosen["execution_arms"]:
        for dataset in DATASETS:
            case, case_pins = selected_case(api, repo, dataset)
            cases.append(case)
            pins += case_pins
            checkpoint = Path(case["checkpoint"])
            pins += [record(path) for path in sorted(checkpoint.iterdir()) if path.is_file()]
            phylo = Path(case["case"]["baseline_output"]) / "orthohmm_phylogeny"
            for name in ("alignments", "gene_trees", "checkpoints", "species_tree_inference"):
                pins += [record(path) for path in sorted((phylo / name).rglob("*")) if path.is_file()]
    pins += driver_pins
    pins = campaign.unique_pins(pins)
    for pin in pins: check(pin)
    result = {"status": "projection_transfer_inputs_prepared_execution_driver_required" if cases else "no_projection_followup_required",
        "scientific_revision": SCIENCE, "scientific_source": pilot["source_root"], "runtime": pilot["runtime"],
        "tool_path": pilot["tool_path"], "pilot_readback": provenance["reader"]["reader"],
        "pilot_plan": report["execution_plan"], "prerequisites": provenance, "selection": chosen,
        "protocol": record(protocol_path), "preparer_revision": revision, "preparer_files": driver_pins,
        "repo": str(repo), "output": str(output), "cpu": 16, "cases": cases, "inputs": pins,
        "defaults_promoted": False, "goal_completed": False,
        "scope": "Authenticated input/role preparation only; no scientific inference or scoring submitted. Complete/test/freeze execution driver and independent full-vector reader before launch. Original gates, fixed scientific settings and conditional population limitations remain unchanged."}
    output.mkdir(parents=True)
    write_json(output / "prepared_inputs.json", result)
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("repo", "driver", "output"): parser.add_argument("--" + name, required=True, type=Path)
    args = parser.parse_args()
    print(prepare(args.repo.resolve(), args.driver.resolve(), args.output.resolve())["status"])
