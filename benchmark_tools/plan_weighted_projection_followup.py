"""Admit the fixed factorial pilot and select only registered matched transfers."""

import argparse
from copy import deepcopy
import hashlib
import json
from pathlib import Path
import subprocess
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from benchmark_tools import report_hmm_evidence_campaign as campaign
from benchmark_tools import report_profile_evidence_pilot as base
from benchmark_tools import report_weighted_projection as reader
from benchmark_tools import weighted_projection_experiment as experiment
from benchmark_tools.prepare_ob_candidate_neighborhood import check, record
from benchmark_tools.run_full_search_hillclimb import runtime
from benchmark_tools.run_membership_hillclimb import write_json
from benchmark_tools.run_profile_weighting_followup import successful_job
from benchmark_tools.run_three_kingdoms_membership_supplement import sources


SCIENCE = "86e39314ba462698b726f569a5d98a144c8d0086"
SUBMISSION = "ddabca65c610855d86fa21402296e8920a9ad35b"
DATASETS = ("QfO", "ThreeKingdoms")


def launch(repo):
    blob = subprocess.check_output(["git", "show", SUBMISSION + ":" + reader.LEDGER], cwd=repo)
    entry = json.loads(blob)["weighted_projection_postrun_readback"]
    science, plan, binding = reader.launch(repo, entry["launch_binding"]["revision"])
    job = entry["job"]
    if (entry["launch_binding"] != binding or science["runner_commit"] != SCIENCE
            or entry["runner_commit"] != SCIENCE
            or entry["scientific_reader_source"] != science["independent_reader"]["scientific_reader_source"]
            or entry["scientific_reader_source"] == science["scientific_source"]
            or science["independent_reader"]["runner_commit"] != SCIENCE
            or science["inference"]["job_id"] != "24538"
            or (job.get("submitted"), job.get("job_id"), job.get("dependency")) != (True, "24539", "afterany:24538")
            or (job.get("cpus"), job.get("memory_allocation"), job.get("wall_limit")) != (1, "16G", "02:00:00")):
        raise ValueError("Changed submitted weighted projection scientific/reader binding")
    return entry, science, plan, {"reader_submission": {
        "revision": SUBMISSION, "path": reader.LEDGER, "git_blob_sha256": hashlib.sha256(blob).hexdigest()},
        "science": binding}


def selection(report, producer):
    if (report.get("status") != "weighted_projection_independently_read_back"
            or report.get("accuracy_result_admitted") is not True
            or report.get("source_revision") != SCIENCE or report.get("reporter_revision") != SCIENCE
            or report.get("defaults_promoted") is not False or report.get("goal_completed") is not False):
        raise ValueError("Require complete separately admitted joint pilot, not partial scores")
    cells = deepcopy(report["cells"])
    if any(cell.get("source_revision") != SCIENCE
           or (cell["coverage"]["genes"], cell["coverage"]["species"]) != (251378, 12) for cell in cells):
        raise ValueError("Wrong joint pilot source or complete gene/species universe")
    gates = reader.verify_panel(cells, producer["arms"])
    passing = [label for label in experiment.LABELS[4:] if gates[label]["eligible"]]
    if report.get("matched_gates") != gates or report.get("advancing_labels") != passing:
        raise ValueError("Joint transfer selection changes original matched gates or hides a passing arm")
    needed = set(passing)
    for label in passing:
        needed.update(experiment.configuration(label)["matched_controls"])
    arms = [label for label in experiment.LABELS if label in needed]
    controls = [label for label in arms if experiment.configuration(label)["role"] != "candidate"]
    return {"advancing_candidates": passing,
        "stopped_candidates": [label for label in experiment.LABELS[4:] if label not in passing],
        "execution_arms": arms,
        "roles": {label: experiment.configuration(label)["role"] for label in arms},
        "arm_parameters": [experiment.configuration(label) for label in arms],
        "control_dependencies": {label: list(controls) if label in passing else [] for label in arms},
        "original_gates": gates,
        "scope": "Only unchanged separately admitted point-gate passes advance; all passing combinations and their matched constituent closure are retained. No fallback, retuning, apparent-winner selection or promotion of controls."}


def authenticate(repo):
    entry, science, plan_pin, binding = launch(repo)
    report, pins, provenance = campaign.read_admission(repo, entry, "weighted_projection_independently_read_back")
    if report.get("launch_ledger") != binding["science"] or report.get("execution_plan") != plan_pin:
        raise ValueError("Joint admission uses a different immutable scientific plan or launch")
    for pin in (report["producer"], plan_pin): check(pin)
    producer = json.loads(Path(report["producer"]["path"]).read_text())
    plan = json.loads(Path(plan_pin["path"]).read_text())
    base.require_completed(producer, plan, experiment.LABELS, "weighted_projection", science["inference"]["job_id"])
    accounting = successful_job(science["inference"]["job_id"], 16, "128G")
    source = (repo / science["scientific_source"]).resolve()
    if (producer.get("execution_plan") != plan_pin or plan.get("source_revision") != SCIENCE
            or Path(plan["source_root"]) != source
            or subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=source, text=True).strip() != SCIENCE
            or sources(source) + [record(source / "benchmark_tools/run_profile_weighting_pilot.slurm")] != plan["sources"]
            or runtime() != plan["runtime"]):
        raise ValueError("Original scientific pilot, frozen source inventory or runtime changed")
    chosen = selection(report, producer)
    pins = campaign.unique_pins([*pins, report["producer"], plan_pin, *plan["sources"], *plan["runtime"]["files"]])
    for pin in pins: check(pin)
    return report, plan, chosen, pins, {"launches": binding, "reader": provenance, "science_accounting": accounting}


def prepare(repo, driver, output):
    if output.exists():
        raise FileExistsError("Preserve every weighted projection transfer selection")
    if driver == repo or Path(__file__).resolve() != driver / "benchmark_tools/plan_weighted_projection_followup.py":
        raise ValueError("Require frozen committed selection planner")
    report, pilot, chosen, pins, provenance = authenticate(repo)
    protocol = driver / "benchmark_tools" / experiment.PROTOCOL
    if protocol.read_bytes() != subprocess.check_output(
            ["git", "show", experiment.PROTOCOL_REVISION + ":benchmark_tools/" + experiment.PROTOCOL], cwd=repo):
        raise ValueError("Changed original factorial/transfer protocol")
    revision = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=driver, text=True).strip()
    planner_files = sources(driver) + [record(protocol), record(driver / "benchmark_tools/plan_weighted_projection_followup.slurm")]
    for pin in planner_files:
        relative = Path(pin["path"]).relative_to(driver)
        if relative.suffix != ".so" and hashlib.sha256(subprocess.check_output(
                ["git", "show", revision + ":" + str(relative)], cwd=driver)).hexdigest() != pin["sha256"]:
            raise ValueError("Selection planner differs from committed source")
    pins = campaign.unique_pins([*pins, *planner_files])
    for pin in pins: check(pin)
    result = {"status": "weighted_projection_followup_selection_prepared_execution_driver_required" if chosen["execution_arms"]
              else "no_weighted_projection_followup_required",
        "selection": chosen, "datasets": list(DATASETS) if chosen["execution_arms"] else [],
        "scientific_revision": SCIENCE, "scientific_source": pilot["source_root"],
        "runtime": pilot["runtime"], "tool_path": pilot["tool_path"],
        "pilot_readback": provenance["reader"]["reader"], "pilot_plan": report["execution_plan"],
        "prerequisites": provenance, "protocol": record(protocol), "planner_revision": revision,
        "planner_files": planner_files, "inputs": pins, "repo": str(repo), "output": str(output),
        "defaults_promoted": False, "goal_completed": False,
        "scope": "Authenticated selection only, no new scientific inference/scoring, case preparation or follow-up submission. Any passing combinations still require tested/frozen execution and separate full-vector readback with all matched controls and unchanged settings. If neither passes, retain the complete informative negative and finish final recommendation/confirmation limits/remote audit."}
    output.mkdir(parents=True)
    write_json(output / "selection.json", result)
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("repo", "driver", "output"): parser.add_argument("--" + name, required=True, type=Path)
    args = parser.parse_args()
    print(prepare(args.repo.resolve(), args.driver.resolve(), args.output.resolve())["status"])
