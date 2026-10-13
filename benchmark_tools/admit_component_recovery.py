"""Independently admit split-process component recovery without rerunning inference."""

import argparse
import csv
import json
from pathlib import Path
import subprocess
import sys


def require_terminal_job(text, job, cpus):
    rows = list(csv.DictReader(text.splitlines(), delimiter="|"))
    if (len(rows) != 1 or rows[0]["JobID"] != str(job) or rows[0]["State"] != "COMPLETED"
            or rows[0]["ExitCode"] != "0:0" or rows[0]["AllocCPUS"] != str(cpus)):
        raise ValueError("Recovery phase has no matching successful terminal job")
    return rows[0]


def check_phase_links(original_pin, producer_pin, probe_pin, candidates, scoring_plan, scored):
    if (candidates["original_execution_plan"] != original_pin or candidates["original_producer"] != producer_pin
            or scoring_plan["original_execution_plan"] != original_pin
            or scoring_plan["original_failed_producer"] != producer_pin or scoring_plan["native_probe"] != probe_pin):
        raise ValueError("Recovery phases refer to different original artifacts")
    if (candidates["status"] != "isolated_component_candidates_exactly_reconstructed"
            or candidates.get("inputs_reverified_unchanged") is not True
            or scored["status"] != "scored_recovery_pending_independent_admission"):
        raise ValueError("Recovery candidate or scoring phase is incomplete")
    if any(phase[key] for phase in (candidates, scored) for key in ("native_inference_rerun", "defaults_promoted", "goal_completed")):
        raise ValueError("Recovery phase exceeded its nonduplicate scientific scope")


def run(original_path, candidate_path, scoring_path, output):
    if output.exists():
        raise ValueError("Require a fresh independent admission namespace")
    original = json.loads(original_path.read_text())
    source = Path(original["source_root"])
    sys.path.insert(0, str(source))
    from benchmark_tools import run_large_copy_components as runner
    from benchmark_tools.prepare_ob_candidate_neighborhood import check, record
    from benchmark_tools.report_hillclimb import endpoint_vector, readback_scores, replay_native
    from benchmark_tools.run_membership_hillclimb import write_json
    if Path(runner.__file__).resolve() != source / "benchmark_tools/run_large_copy_components.py":
        raise ValueError("Admission loaded the wrong original scientific helpers")
    candidate = json.loads(candidate_path.read_text())
    scoring = json.loads(scoring_path.read_text())
    scoring_plan_path = Path(scoring["execution_plan"]["path"])
    scoring_plan = json.loads(scoring_plan_path.read_text())
    probe_path = Path(scoring_plan["native_probe"]["path"])
    probe = json.loads(probe_path.read_text())
    producer_path = Path(original["output"]) / "results.json"
    producer = json.loads(producer_path.read_text())
    check_phase_links(record(original_path), record(producer_path), record(probe_path), candidate, scoring_plan, scoring)
    if (probe["status"] != "clean_original_native_readback_verified_scoring_still_pending"
            or any(phase["source_revision"] != original["source_revision"] for phase in (candidate, scoring, probe))):
        raise ValueError("Recovery native verification or source identity differs")
    schedulers = {}
    for phase, job, cpus in (("native_probe", scoring_plan["native_probe_scheduler"]["JobID"], 1),
                             ("candidate", candidate["job_id"], 16), ("scoring", scoring["job_id"], 16)):
        text = subprocess.check_output(["sacct", "-j", str(job), "-X",
            "--format=JobID,State,ExitCode,AllocCPUS", "--parsable2"], text=True)
        schedulers[phase] = require_terminal_job(text, job, cpus)
    checked = [record(original_path), record(producer_path), record(probe_path), record(candidate_path),
               record(scoring_path), record(scoring_plan_path), record(Path(__file__)),
               *probe["checked"], *probe["native"]["checked"], *candidate["checked"]]
    if candidate["component_counts"] != producer["component_counts"] or candidate["ownership"] != producer["ownership"]:
        raise ValueError("Independently reconstructed component counts or ownership differ")
    for label in ("control", "components500"):
        old, new = producer["structural_controls"][label], candidate["structural_controls"][label]
        for key in ("seed", "copy_split_trace"):
            if old[key]["sha256"] != new[key]["sha256"]:
                raise ValueError("Independent reconstructed seed or trace differs")
            checked.extend((old[key], new[key]))
        for key in ("seed_partition", "candidate_partition", "membership_constraints"):
            if old["candidate"][key]["sha256"] != new["candidate"][key]["sha256"]:
                raise ValueError("Independent candidate or membership constraints differ")
            checked.extend((old["candidate"][key], new["candidate"][key]))
    runner.verify(original, checked)
    output.mkdir(parents=True)
    report = {"status": "independent_split_recovery_admission_running", "scores": {}, "schedulers": schedulers,
              "native_inference_rerun": False, "defaults_promoted": False, "goal_completed": False}
    write_json(output / "report.json", report)
    try:
        row = {"outputs": list(probe["native"]["products"].values()), "coverage": probe["coverage"]}
        native = replay_native(row, original["case"])
        runner.require_markers(native["summary"])
        if native["genes"] != 984137 or native["species"] != 78:
            raise ValueError("Independent native replay lost full dataset coverage")
        scores, details, files = readback_scores(Path(original["repo"]), "QfO", scoring["scoring"], native,
            scoring_path.parent, "ohmm_hillclimb_r1_" + scoring_plan["scoring_id"])
        vector = endpoint_vector("QfO", scores)
        coverage = {k: native[k] for k in ("genes", "species", "groups", "pairs")}
        if (vector != scoring["scores"] or details != scoring["details"] or native["products"] != scoring["native_products"]
                or native["summary"] != scoring["native_summary"] or coverage != scoring["coverage"]):
            raise ValueError("Independently reconstructed native products or scores differ")
        checked.extend([*native["checked"], *files])
        runner.verify(original, checked)
        report.update(status="terminal_component_candidates_and_raw_scores_independently_verified", dataset="QfO",
            method_key="large_copy_components:components500", scores=vector, details=details, coverage=coverage,
            products=native["products"], native_summary=native["summary"], component_counts=candidate["component_counts"],
            execution_kind="recovered_scoring_on_retained_candidate_specific_own_tree", execution_plan=record(original_path),
            original_failed_producer=record(producer_path), native_probe=record(probe_path),
            independent_candidates=record(candidate_path), scoring_producer=record(scoring_path),
            checked=checked, reader=record(Path(__file__)),
            scope="Independent clean native replay, separate exact candidate reconstruction and raw six-endpoint scoring admission; no repeated native inference or fresh full-pipeline timing.")
    except Exception as error:
        report.update(status="failed_independent_component_recovery_admission", error_type=type(error).__name__, error=str(error))
    write_json(output / "report.json", report)
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("original-plan", "candidate-report", "scoring-report", "output"):
        parser.add_argument("--" + name, type=Path, required=True)
    args = parser.parse_args()
    result = run(args.original_plan.resolve(), args.candidate_report.resolve(), args.scoring_report.resolve(), args.output.resolve())
    print(result["status"])
    raise SystemExit(0 if result["status"] == "terminal_component_candidates_and_raw_scores_independently_verified" else 1)
