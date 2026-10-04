"""Integrate retained historical supplements without replacing score-run evidence."""

import argparse
from copy import deepcopy
import csv
import itertools
import json
import math
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from benchmark_tools.assemble_benchmark_provenance import KEYS, memory_info, same_score
from benchmark_tools.prepare_ob_candidate_neighborhood import record


DATASETS = ("OrthoBench", "QfO", "ThreeKingdoms")
SOURCES = {
    "register": ("all_benchmark_provenance_20261004_v2/register.json",
                 "48be4b216f9f3776354246894c8c2a38f790fd4f075528e88f56b230a0a3c47c"),
    "supplement": ("three_kingdoms_run_metadata_20261004.json",
                   "bb6d910021583da61ba4892690c5bf50796fcf81cafcfaeb0d09742e41b8a093"),
    "scores": ("current_benchmark_scores_20260926_v2/manifest.json",
               "8d56dae74721c4530ab2b86c5008b2d0ebac07878e574fd85b1ad84e4b908f07"),
}
NEW_FIELDS = ("historical_metadata_supplement", "supplemental_commands", "supplemental_resources")
TIME_SCOPES = {
    "time.log": "historical native-command GNU-time interval",
    "conversion_time.log": "checkpoint conversion only; not separate sequence-only inference",
    "bpo_converter.time.log": "BLAST-to-BPO conversion only",
    "downstream_time.log": "recovered mode4 downstream only; excludes BLAST/BPO",
}


def require(condition, message):
    if not condition:
        raise ValueError(message)


def argv(command):
    require(isinstance(command, list) and bool(command)
            and all(isinstance(x, str) and x for x in command), "Invalid retained command argv")
    return command


def measurements(values):
    require("elapsed_seconds" in values, "Missing interval")
    for key in ("elapsed_seconds", "user_seconds", "system_seconds"):
        require(type(values.get(key)) in (int, float) and math.isfinite(values[key]) and values[key] >= 0,
                "Invalid retained resource value")
    return values


def supplement_details(row, baseline):
    require(isinstance(baseline["output_records"], list) and len(baseline["output_records"]) == 1
            and row["selected_prediction"] == baseline["output_records"][0],
            "Supplement prediction differs from score row")
    require(row["declared_version"] == baseline["declared_version"], "Conflicting version declaration")
    commands, resources = [], []
    metrics = row["metrics"]
    if metrics is not None:
        require(any((p["bytes"], p["sha256"]) == (row["selected_prediction"]["bytes"],
                                                  row["selected_prediction"]["sha256"])
                    for p in metrics["matched_output_manifest"]), "Metrics output binding differs")
        evidence = row["evidence"]["metrics.json"]
        commands.extend([
            {"scope": "Recorded metrics entrypoint", "argv": argv(metrics["entrypoint_argv"]), "evidence": evidence},
            {"scope": "Harness-recorded launch argv; manifest collected after execution",
             "argv": argv(metrics["harness_argv"]), "evidence": evidence},
        ])
        resources.append({"scope": "Historical pipeline metrics; not newly attested native execution/accounting",
                          "measurement": measurements(metrics["measurement"]), "memory": metrics["memory"],
                          "cpu_scope": metrics["cpu_scope"], "evidence": evidence,
                          "independent_repeat": False})
    for timing in row["timing_records"]:
        name = Path(timing["evidence"]["path"]).name
        require(name in TIME_SCOPES and timing["evidence"] == row["evidence"][name],
                "Unknown or conflicting timing-record binding")
        values = measurements(timing["measurement"])
        require(values["exit_status"] == 0, "Retained native interval failed")
        scope = TIME_SCOPES[name]
        if row["key"] == "fastoma_0_3_5":
            scope += "; supplied-tree workflow; CPU/RSS exclude aggregate Docker tasks"
        commands.append({"scope": scope, "argv": argv(timing["argv"]), "evidence": timing["evidence"]})
        resources.append({"scope": scope, "measurement": values, "evidence": timing["evidence"],
                          "memory_scope": timing["memory_scope"], "independent_repeat": False})
    return commands, resources


def integrate(register, supplement, scores):
    require(register["status"] == "all_benchmark_provenance_register_partial"
            and register["complete_transitive_provenance"] is False
            and register["scores_recomputed"] is False, "Unexpected baseline register scope")
    require(supplement["status"] == "historical_three_kingdoms_metadata_supplement"
            and supplement["historical_consumption_proven"] is False
            and supplement["native_runs_repeated"] is False
            and supplement["scores_recomputed"] is False, "Unexpected supplement scope")
    rows = deepcopy(register["rows"])
    identities = [(row["dataset"], row["key"]) for row in rows]
    require(len(rows) == 24 and set(identities) == set(itertools.product(DATASETS, KEYS)),
            "Require exactly 24 unique method/dataset cells")
    require(all(not any(field in row for field in NEW_FIELDS) for row in rows), "Already integrated register")
    current = {row["key"]: row for row in scores["rows"]}
    require(len(scores["rows"]) == 8 and set(current) == set(KEYS), "Incomplete selected score inventory")
    supplements = {row["key"]: row for row in supplement["rows"]}
    require(len(supplement["rows"]) == 7 and set(supplements) == set(KEYS) - {"sonicparanoid_2_0_9"},
            "Require seven historical rows; exclude replaced Sonic run")
    for row in rows:
        selected = current[row["key"]]["scores"]
        if row["dataset"] == "QfO":
            for metric, value in row["scores"].items():
                same_score(value, selected[metric])
            same_score(row["secondary_mean"], selected["QfO_secondary_mean"])
        else:
            require(len(row["scores"]) == 1, "Unexpected score endpoint inventory")
            same_score(next(iter(row["scores"].values())), selected[row["dataset"]])
        if row["dataset"] != "ThreeKingdoms" or row["key"] not in supplements:
            continue
        raw = supplements[row["key"]]
        commands, resources = supplement_details(raw, row)
        row.update(historical_metadata_supplement=deepcopy(raw), supplemental_commands=commands,
                   supplemental_resources=resources)
    # Strip only additions: every earlier score, command, resource and limitation must survive.
    require([{k: v for k, v in row.items() if k not in NEW_FIELDS} for row in rows] == register["rows"],
            "Historical register fields changed")
    return rows


def resource_rows(rows):
    table = []
    for row in rows:
        previous = row["resources"] or [{"scope": "Unestablished; no full-inference interval",
                                         "measurement": {}}]
        for origin, records in (("previous_register", previous),
                                ("historical_directory_supplement", row.get("supplemental_resources", []))):
            for interval in records:
                values, memory = interval["measurement"], memory_info(interval)
                table.append({"dataset": row["dataset"], "method": row["key"],
                              "declared_version": row["declared_version"], "origin": origin,
                              "scope": interval["scope"], "wall_seconds": values.get("elapsed_seconds"),
                              "cpu_user_seconds": values.get("user_seconds"),
                              "cpu_system_seconds": values.get("system_seconds"),
                              "memory_value": memory["value"], "memory_unit": memory["unit"],
                              "memory_scope": memory["scope"],
                              "wall_measurement_status": "unavailable" if values.get("elapsed_seconds") is None else "recorded",
                              "independent_repeat": False})
    return table


def render(report):
    lines = ["# Integrated All-Tool Benchmark Metadata", "",
             "Eight methods across three datasets. Scores and all previous register fields are unchanged.",
             "Historical command/resource supplements are retained associations, not new execution attestations.", "",
             "| Dataset | Method | Declared Version | Output Semantics | Command Evidence | Supplemental Intervals |",
             "| --- | --- | --- | --- | --- | ---: |"]
    for row in report["rows"]:
        command = "previous register" if row["commands"] else "unavailable"
        if row.get("supplemental_commands"):
            command += "; retained metadata supplement"
            if command.startswith("unavailable;"):
                command = "retained metadata supplement"
        lines.append("| %s | %s | %s | %s | %s | %d |" % (row["dataset"], row["label"],
                     row["declared_version"], row["output_semantics"], command,
                     len(row.get("supplemental_resources", []))))
    lines += ["", "## Limits", "", *["- " + text for text in report["limitations"]]]
    return "\n".join(lines) + "\n"


def collect(repo):
    repo = Path(repo).resolve()
    refs, documents = [], {}
    for name, (path, digest) in SOURCES.items():
        ref = record(repo / "benchmark_tools/results" / path)
        require(ref["sha256"] == digest, "Retained metadata source changed: " + name)
        refs.append(ref)
        documents[name] = json.loads(Path(ref["path"]).read_text())
    rows = integrate(documents["register"], documents["supplement"], documents["scores"])
    for ref in refs:
        require(record(ref["path"]) == ref, "Source changed during integration")
    intervals = resource_rows(rows)
    return {"schema": "integrated_benchmark_metadata_v1", "status": "retained_all_tool_metadata_integrated_partial",
            "inputs": refs, "source": record(__file__), "rows": rows, "resource_intervals": intervals,
            "resource_table_entries": len(intervals),
            "recorded_wall_observations": sum(row["wall_seconds"] is not None for row in intervals),
            "method_dataset_cells": 24, "historical_rows_supplemented": 7,
            "previous_register_fields_unchanged": True, "native_inference_or_scoring_repeated": False,
            "historical_input_consumption_proven": False, "complete_transitive_provenance": False,
            "controlled_comparative_resources": False, "publication_ready": False,
            "limitations": [
                "Three fixed metadata/score receipts are directly rehashed; their prior raw prediction/input/runtime audits are reused, not rerun. Transitive artifact identities remain inherited.",
                "All 24 original rows, including exact scores, versions, semantics, commands, input notes, gaps and resources, are preserved. Version text is not universal runtime attestation.",
                "Seven historical Three Kingdoms supplements add recorded argv, metrics, source/harness chronology and recovery events. Directory/output association is not immutable historical execution or input-consumption proof.",
                "The contemporary matched Sonic row is unchanged; do not attach its replaced historical run's timing or commands.",
                "Thirty-seven table entries preserve 29 earlier resource entries, including unavailable placeholders, and add eight supplemental intervals. Not 37 independent runs: wrappers, native intervals, cached replay, checkpoint conversion and recovered stages are not summed, replaced or treated as independent repeats.",
                "GNU-time process RSS, sampled tree RSS, driver/container scopes and historical allocations remain distinct. Unknown is not zero. Shared-host contention is unknown and potentially method dependent; no isolated speed ranking.",
                "The sequence-only OrthoFinder output remains a same-full-run MCL checkpoint, not a separately timed sequence-only inference. FastOMA supplied-tree dependence remains disclosed.",
                "The completed replacement timing panel and original factorial resource observations remain separate. No new cost is imputed to unmatched configurations or historical score runs.",
                "QfO GO/EC/FAS are similarities and the six-metric mean is secondary. OrthoBench group F1 and Three Kingdoms BUSCO co-membership are not interchangeable endpoints or universal generalization evidence.",
            ]}


def write(repo, output):
    output = Path(output)
    if output.exists() or output.is_symlink():
        raise FileExistsError(output)
    report = collect(repo)
    output.mkdir(parents=True)
    (output / "register.json").write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    (output / "register.md").write_text(render(report))
    with (output / "resources.tsv").open("x", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(report["resource_intervals"][0]),
                                delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows({k: "NA" if v is None else v for k, v in row.items()}
                        for row in report["resource_intervals"])
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    report = write(args.repo, args.output)
    print(json.dumps({"status": report["status"], "cells": len(report["rows"]),
                      "resource_intervals": len(report["resource_intervals"])}))


if __name__ == "__main__":
    main()
