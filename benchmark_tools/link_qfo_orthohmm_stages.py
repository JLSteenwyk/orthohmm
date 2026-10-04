"""Bind selected QfO predictions to retained, explicitly incremental stage evidence."""

import argparse
from copy import deepcopy
import csv
import json
import math
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from benchmark_tools.assemble_benchmark_provenance import Reader, bind_conversion
from benchmark_tools.audit_ob_orthofinder_provenance import log_command
from benchmark_tools.collect_factorial_resources import metrics_values, require
from benchmark_tools.prepare_ob_candidate_neighborhood import check, record
from benchmark_tools.summarize_matched_resources import verbose_time


SOURCES = {
    "register": ("all_benchmark_metadata_integrated_20261004_v3/register.json",
                 "4c6f01afb28ef36c25991155c2484a6e2db273556529506350c2f7912083ba80"),
    "selected": ("qfo_corrected_comparison_20260926_v7/manifest.json",
                 "042aa221d01554114a1b0b413ca3e2ff56dad5f8c06362c69a1bf49f2de642fc"),
    "costs": ("factorial_retained_resources_20261004/resources.json",
              "4c7801197a9bc89927964f952d03325187341a0ce2fe2cc5f6c4aa7700e581a6"),
}
CELLS = {"orthohmm_high_sensitivity": "p1_c0_r0",
         "orthohmm_phylogeny_satellite_v2": "p1_c1_r1"}
ADDED_FIELD = "qfo_stage_provenance"


def command(argv):
    require(isinstance(argv, list) and bool(argv)
            and all(isinstance(x, str) and x for x in argv), "Invalid recorded stage command")
    return argv


def stage(scope, evidence, wall, user=None, system=None, memory=None, unit=None, memory_scope=None):
    for value in (wall, user, system, memory):
        require(value is None or (type(value) in (int, float) and math.isfinite(value) and value >= 0),
                "Invalid stage measurement")
    require(wall is not None and wall > 0, "Missing stage wall time")
    return {"scope": scope, "evidence": evidence, "measurement": {
        "elapsed_seconds": wall, "user_seconds": user, "system_seconds": system},
        "memory": {"value": memory, "unit": unit, "scope": memory_scope},
        "full_inference": False, "independent_repeat": False}


def link(row, selected, conversion, candidates, prepared, replay_admission,
         replay_report, replay, replay_plan, timing, costs, native=None, metrics=None):
    cell = CELLS[row["key"]]
    arm_name = cell[:-3]
    require(row["dataset"] == "QfO" and selected["status"] == "admitted"
            and selected["factorial_cell"] == cell, "Wrong selected method/cell")
    require(row["conversion"] == selected["conversion"]
            and row["score_admission"] == selected["admission"]
            and row["scores"] == selected["scores"]
            and row["secondary_mean"] == selected["secondary_mean"], "Selected score binding differs")
    require(conversion["cell"] == cell and row["output_records"] == [bind_conversion(selected, conversion)],
            "Selected pair output differs")
    require(candidates["status"] == "corrected_qfo_candidates_admitted"
            and candidates["prepared_manifest"] == conversion["prepared"], "Wrong candidate preparation binding")
    require(prepared["status"] == "corrected_qfo_four_candidate_arms_prepared_unscored"
            and candidates["candidate_arms"][arm_name] == prepared["candidate_arms"][arm_name]["content_audit"],
            "Candidate arm differs from admission")
    arm = prepared["candidate_arms"][arm_name]
    require(replay_admission["status"] == "corrected_checked_replay_admitted"
            and replay_report["status"] == "corrected_checked_replay_complete_pending_admission"
            and replay_report["exit_code"] == 0, "Cached replay not admitted/completed")
    require(replay_admission["plan"] == prepared["plan"] == replay_report["plan"],
            "Replay plan chain differs")
    require(replay_admission["coverage"] == replay_report["coverage"], "Replay coverage chain differs")
    seed = [entry for entry in replay_admission["coverage"] if entry["label"] == "strict_profiles_refined"]
    require(len(seed) == 1 and seed[0]["output"] == arm["seed_partition"], "Arm uses different profile seed")
    require(replay["stages"][-1]["output"] == seed[0]["output"]
            and replay["command"] == replay_plan["native_command"]
            and replay["input"]["manifest"] == replay_plan["checkpoint_manifest"],
            "Replay command/checkpoint/output differs")
    require(prepared["input_fastas"] == replay_plan["input_fastas"] == conversion["input_fastas"],
            "Corrected FASTA manifest differs")
    require(row["input_records"] == conversion["input_fastas"], "Register input manifest differs")
    require(timing["argv"] == replay_report["worker_command"] and timing["measurement"]["exit_status"] == 0,
            "Cached wrapper timing command differs")
    require(costs["dataset"] == "Corrected QfO" and costs["cell"] == cell
            and costs["candidate_preparation_wall_s"] == arm["incremental_preparation_seconds"]
            and costs["full_pipeline_wall_s"] is None, "Stage cost table binding differs")
    require(row["commands"] == [] and row["resources"] == [], "Unexpected prior QfO stage evidence")
    require(row["native_admission"] == conversion.get("native_admission", conversion["candidate_admission"]),
            "Register admission chain differs")
    intervals = [stage("Shared checked cached replay worker; excludes initial search, parent validation and downstream stages",
                       timing["evidence"], timing["measurement"]["elapsed_seconds"],
                       timing["measurement"]["user_seconds"], timing["measurement"]["system_seconds"],
                       timing["measurement"]["max_process_rss_kib"], "KiB",
                       "GNU-time maximum process RSS; not aggregate workflow memory"),
                 stage("Candidate arm preparation from cached seed; shared with the arm's other R cell",
                       conversion["prepared"], arm["incremental_preparation_seconds"])]
    commands = [{"scope": "Checked cached replay worker", "argv": command(replay_report["worker_command"]),
                 "evidence": replay_admission["source_report"]},
                {"scope": "Replay's recorded internal entrypoint, not a separate timed process",
                 "argv": command(replay["command"]), "evidence": replay_report["replay"]}]
    reconciliation = None
    if cell.endswith("r0"):
        require(conversion["status"] == "corrected_factorial_group_pairs_prepared_unscored"
                and conversion["candidate_partition"] == arm["candidate_partition"],
                "Group-pair conversion uses a different candidate partition")
        require(costs["reconciliation_measurement_status"] == "not_applicable", "Unexpected R-off cost")
        commands.append({"scope": "Recorded group-to-pair conversion; no retained native interval",
                         "argv": command(conversion["command"]), "evidence": selected["conversion"]})
    else:
        require(native is not None and metrics is not None, "Missing native reconciliation evidence")
        require(conversion["status"] == "corrected_factorial_native_pairs_prepared_unscored"
                and native["status"] == "corrected_qfo_native_pair_output_verified"
                and native["cell"] == cell and native["native_pairs"] == conversion["native_input"]
                and native["candidate_admission"] == conversion["candidate_admission"],
                "Native pair/admission binding differs")
        group = native["native_group_integrity"]
        require(group["status"] == "native_group_output_verified" and group["cell"] == cell,
                "Wrong native group admission")
        values = metrics_values(metrics)
        require(metrics["input"]["candidate_clusters"] == arm["candidate_partition"]
                and metrics["input"]["membership_constraints"] == arm["membership_constraints"],
                "Reconciliation candidate/constraint binding differs")
        require(metrics["outputs"]["manifest"] == group["native_manifest"], "Native metrics output differs")
        require(all(costs[k] == values[v] for k, v in (
            ("reconciliation_wall_s", "wall_s"), ("reconciliation_user_cpu_s", "user_cpu_s"),
            ("reconciliation_system_cpu_s", "system_cpu_s"),
            ("reconciliation_peak_sampled_tree_rss_bytes", "peak_process_tree_rss_bytes"))),
            "Reconciliation cost differs from metrics")
        intervals.append(stage("Cached-candidate native reconciliation; excludes initial search, replay, preparation, conversion and scoring",
                               group["native_metrics"], values["wall_s"], values["user_cpu_s"],
                               values["system_cpu_s"], values["peak_process_tree_rss_bytes"], "bytes",
                               metrics["rss_measurement"]))
        commands.append({"scope": "Cached-candidate native reconciliation", "argv": command(metrics["command"]),
                         "evidence": group["native_metrics"]})
        reconciliation = {"metrics": group["native_metrics"], "native_pairs": native["native_pairs"],
                          "candidate_partition": arm["candidate_partition"],
                          "membership_constraints": arm["membership_constraints"],
                          "scheduler": group["integrity"]["scheduler"], "parameters": metrics["parameters"]}
    return {"cell": cell, "candidate_arm": arm_name, "selected_prediction": row["output_records"][0],
            "candidate_partition": arm["candidate_partition"], "profile_seed": arm["seed_partition"],
            "candidate_admission": conversion["candidate_admission"], "prepared": conversion["prepared"],
            "replay_admission": prepared["admission"], "replay_report": replay_admission["source_report"],
            "replay_metrics": replay_report["replay"], "replay_plan": prepared["plan"],
            "replay_internal_wall_s": replay["wall_s"], "replay_internal_timings": replay["timings"],
            "replay_scheduler": replay_admission["scheduler"], "reconciliation": reconciliation,
            "commands": commands, "resource_intervals": intervals,
            "full_pipeline_wall_s": None, "full_pipeline_cpu_s": None, "full_pipeline_peak_memory_bytes": None,
            "observations_are_not_summed": True, "historical_consumption_reaudited": False}


def collect(repo):
    repo = Path(repo).resolve()
    reader, docs = Reader(), {}
    for name, (path, digest) in SOURCES.items():
        pin = record(repo / "benchmark_tools/results" / path)
        require(pin["sha256"] == digest, "Fixed input changed: " + name)
        docs[name] = reader.read(pin)
    baseline = docs["register"]
    require(baseline["method_dataset_cells"] == 24 and len(baseline["rows"]) == 24,
            "Incomplete integrated register")
    rows = deepcopy(baseline["rows"])
    selected = {r["key"]: r for r in docs["selected"]["methods"]}
    new_intervals = []
    for row in rows:
        require(ADDED_FIELD not in row, "Already linked metadata")
        if row["dataset"] != "QfO" or row["key"] not in CELLS:
            continue
        selection = selected[row["key"]]
        conversion = reader.read(selection["conversion"])
        score = reader.read(selection["admission"])
        require(score["status"] == "corrected_factorial_assessment_admitted"
                and score["accuracy_admitted"] is True and score["conversion"] == conversion,
                "Score admission does not bind conversion")
        candidates = reader.read(conversion["candidate_admission"])
        prepared = reader.read(candidates["prepared_manifest"])
        replay_admission = reader.read(prepared["admission"])
        replay_report = reader.read(replay_admission["source_report"])
        replay = reader.read(replay_report["replay"])
        replay_plan = reader.read(replay_admission["plan"])
        timing_pin = record(Path(replay_admission["source_report"]["path"]).with_name("time.txt"))
        require(timing_pin in replay_admission["checked_records"], "Timing log not bound by replay admission")
        text = reader.text(timing_pin)
        timing = {"evidence": timing_pin, "argv": log_command(text, "Command being timed: "),
                  "measurement": verbose_time(text)}
        costs = [r for r in docs["costs"]["rows"] if r["dataset"] == "Corrected QfO"
                 and r["cell"] == CELLS[row["key"]]]
        require(len(costs) == 1, "Missing/duplicate selected stage-cost row")
        native = metrics = None
        if "native_admission" in conversion:
            native = reader.read(conversion["native_admission"])
            metrics = reader.read(native["native_group_integrity"]["native_metrics"])
        linked = link(row, selection, conversion, candidates, prepared, replay_admission,
                      replay_report, replay, replay_plan, timing, costs[0], native, metrics)
        row[ADDED_FIELD] = linked
        for interval in linked["resource_intervals"]:
            values, memory = interval["measurement"], interval["memory"]
            new_intervals.append({"dataset": "QfO", "method": row["key"],
                "declared_version": row["declared_version"], "origin": "retained_qfo_stage_linkage",
                "scope": interval["scope"], "wall_seconds": values["elapsed_seconds"],
                "cpu_user_seconds": values["user_seconds"], "cpu_system_seconds": values["system_seconds"],
                "memory_value": memory["value"], "memory_unit": memory["unit"], "memory_scope": memory["scope"],
                "wall_measurement_status": "recorded", "independent_repeat": False})
    require(sum(ADDED_FIELD in r for r in rows) == 2 and len(new_intervals) == 5, "Incomplete stage linkage")
    require([{k: v for k, v in row.items() if k != ADDED_FIELD} for row in rows] == baseline["rows"],
            "Earlier register fields changed")
    for pin in reader.checked.values():
        require(Path(pin["path"]).resolve().is_relative_to(repo), "Direct evidence outside repository")
        check(pin)
    intervals = deepcopy(baseline["resource_intervals"]) + new_intervals
    return {"schema": "qfo_orthohmm_stage_metadata_v1", "status": "selected_qfo_stage_provenance_linked_partial",
            "source": record(__file__), "checked_records": list(reader.checked.values()), "rows": rows,
            "method_dataset_cells": 24, "qfo_rows_supplemented": 2, "resource_intervals": intervals,
            "resource_table_entries": len(intervals),
            "recorded_wall_observations": sum(r["wall_seconds"] is not None for r in intervals),
            "previous_register_fields_unchanged": True, "native_inference_or_scoring_repeated": False,
            "complete_transitive_provenance": False, "controlled_comparative_resources": False,
            "publication_ready": False, "limitations": baseline["limitations"] + [
                "Two selected QfO OrthoHMM chains are linked through exact score/conversion/candidate/replay/admission metadata. Raw predictions, FASTAs, numeric hits and source/runtime inventories are inherited from prior admissions, not re-audited here.",
                "Five stage associations add a shared replay worker in both rows, two shared candidate-arm intervals and one reconciliation interval. These are not five new runs or independent repetitions. The replay wall/CPU is a checked worker wrapper; internal replay timings are a distinct scope.",
                "Initial HMM search is absent from cached replay costs. No stage sum, memory sum or full-pipeline cost is inferred. Full costs for the exact selected cached executions remain unavailable; upstream native equivalence is not used as a cost substitute.",
                "Candidate preparation CPU and memory, and native conversion costs are unestablished. Recorded argv is not a new executable/dependency or historical consumption attestation. All shared-host distortion remains unknown and potentially method dependent.",
            ]}


def render(report):
    lines = ["# Selected QfO OrthoHMM Stage Provenance", "",
             "All 24 earlier rows and their scores remain unchanged. New records are incremental, not full inference.",
             "The same checked replay observation appears in both rows; it is not two independent runs.", "",
             "| Method | Cell | Cached Replay Worker (s) | Arm Preparation (s) | Reconciliation (s) | Full Pipeline |",
             "| --- | --- | ---: | ---: | ---: | --- |"]
    for row in report["rows"]:
        if ADDED_FIELD not in row:
            continue
        linked = row[ADDED_FIELD]
        values = [i["measurement"]["elapsed_seconds"] for i in linked["resource_intervals"]]
        lines.append("| %s | %s | %.3f | %.3f | %s | NA |" % (
            row["label"], linked["cell"], values[0], values[1],
            "NA (R off)" if len(values) == 2 else "%.3f" % values[2]))
    lines += ["", "## Limits", "", *["- " + s for s in report["limitations"][-4:]]]
    return "\n".join(lines) + "\n"


def write(repo, output):
    output = Path(output)
    if output.exists() or output.is_symlink():
        raise FileExistsError(output)
    report = collect(repo)
    output.mkdir(parents=True)
    (output / "register.json").write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    (output / "stages.md").write_text(render(report))
    with (output / "resources.tsv").open("x", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(report["resource_intervals"][0]),
                                delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows({k: "NA" if v is None else v for k, v in r.items()}
                        for r in report["resource_intervals"])
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    result = write(args.repo, args.output)
    print(json.dumps({"status": result["status"], "direct_records": len(result["checked_records"]),
                      "resource_entries": result["resource_table_entries"]}))
