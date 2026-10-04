"""Replay retained reporting summaries without opening historical evidence paths."""

import argparse
from collections import Counter
import csv
import hashlib
import json
import math
from pathlib import Path
from statistics import median


INPUTS = {
    "inventory": ("inventory.json", "d5a9b2b778164fd86ed069d379a96eaaedf7be541c021d10c0534acaba154761"),
    "trace": ("diagnostic.json", "11d576cbdf671dc777bdd3c7a0926d80e15922e095000785c82f9f3e332fffbb"),
    "native": ("linkage.json", "896063b253c06f8691068231d241f8cc6eb2c7f588baee4efa1a6e72a42d4cac"),
    "qfo": ("qfo_register.json", "a1cf40e172614f7c9a5e7906cb4121dce17ac7f6c8ff516977a421de05ac867a"),
    "factorial": ("factorial_resources.json", "4c7801197a9bc89927964f952d03325187341a0ce2fe2cc5f6c4aa7700e581a6"),
    "all_tools": ("all_tool_register.json", "4c6f01afb28ef36c25991155c2484a6e2db273556529506350c2f7912083ba80"),
}


def require(ok, message):
    if not ok:
        raise ValueError(message)


def associations(inventory):
    return [{"dataset": e["dataset"], "family": family, "evidence_file": e["file"],
             "json_pointer": e["pointer"], "reported_partition": e["reported_partition"],
             "count_payload_sha256": e["count_payload_sha256"], "evidence_kind": e["evidence_kind"],
             "independent_native_experiment": "unestablished"}
            for e in inventory["evaluations"] for family in e["families"]]


def exposure(inventory):
    require(inventory["publication_ready"] is False and inventory["family_disjoint_validation_established"] is False,
            "Unsupported inventory scope")
    files = inventory["files"]
    require(len(files) == inventory["json_files_scanned"] and len({f["path"] for f in files}) == len(files),
            "Duplicate/missing inventory file")
    require(sum(f["bytes"] for f in files if f["origin"] == "git_snapshot") == inventory["snapshot_json_bytes"]
            and sum(f["bytes"] for f in files if f["origin"] == "frozen_local_only") == inventory["local_json_bytes"],
            "Inventory byte summary differs")
    blocks = inventory["evaluations"]
    require(len({(e["file"], e["pointer"], e["dataset"]) for e in blocks}) == len(blocks), "Duplicate score block")
    summary = {}
    names = {}
    for dataset in ("OrthoBench", "QfO_SwissTrees"):
        rows = [r for r in inventory["family_rows"] if r["dataset"] == dataset]
        names[dataset] = {r["family"] for r in rows}
        require(len(rows) == len(names[dataset]), "Duplicate canonical family")
        subset = [e for e in blocks if e["dataset"] == dataset]
        for e in subset:
            require(len(set(e["families"])) == len(e["families"]) and set(e["families"]) <= names[dataset],
                    "Duplicate/unknown scored family")
        summary[dataset] = {"reference_families": len(rows), "scored_blocks": len(subset),
            "scored_files": len({e["file"] for e in subset}),
            "family_block_associations": sum(len(e["families"]) for e in subset),
            "distinct_count_vectors": len({e["count_payload_sha256"] for e in subset}),
            "families_with_scored_evidence": len({f for e in subset for f in e["families"]})}
        for row in rows:
            hits = [e for e in subset if row["family"] in e["families"]]
            require(row["scored_evidence_blocks"] == len(hits)
                    and row["scored_evidence_files"] == len({e["file"] for e in hits})
                    and row["reported_partition_block_counts"] == dict(Counter(e["reported_partition"] or "not_declared" for e in hits)),
                    "Family exposure projection differs")
            require(row["independent_native_experiments"] is None and row["causal_tuning_influence_established"] is False,
                    "Unsupported experiment/tuning attribution")
    require(summary == inventory["summary"] and len(blocks) == sum(v["scored_blocks"] for v in summary.values()),
            "Inventory summary differs")
    return {"files": len(files), "families": sum(len(n) for n in names.values()), "summary": summary,
            "associations": len(associations(inventory)),
            "unresolved_kinds": dict(Counter(r["kind"] for r in inventory["unresolved_family_containers"]))}


def candidate_trace(trace):
    fixture = trace["fixture"]
    require(fixture["genes"] == 19 and len(fixture["cases"]) == 3, "Unexpected fixture")
    cases = []
    for row in fixture["cases"]:
        partition = row["partition"]
        require(set(g for group in partition for g in group) == set(range(19))
                and sum(len(group) for group in partition) == 19, "Invalid fixture partition")
        anchor = [set(group) for group in partition if set(range(9, 19)) <= set(group)]
        require(len(anchor) == 1, "Lost anchor identity")
        unattached = sorted(set(range(9)) - anchor[0])
        require(unattached == row["unattached_satellites"] and len(row["selections"]) == row["merges"],
                "Satellite identity or merge count differs")
        cases.append({"case": row["case"], "unattached_ids": unattached,
                      "unattached_count": len(unattached), "merges": row["merges"]})
    return {"fixture_cases": cases, "retained_trace_common_merges": [p["comparison"]["common_semantic_merges"] for p in trace["points"]],
            "retained_support_maximum_delta": max(p["comparison"]["common_feature_differences"]["support"]["maximum_absolute_delta"] for p in trace["points"]),
            "raw_trace_replayed": False, "fixture_engine_rerun": False}


def native_costs(native):
    summaries = []
    seen = set()
    for row in native["summaries"]:
        cell = row["cell"]
        require(cell not in seen, "Duplicate native cell")
        seen.add(cell)
        points = [p for p in native["points"] if p["cell"] == cell]
        require(len(points) == row["repeats"] == 3 and {p["repeat"] for p in points} == {0, 1, 2}
                and all(p["comparative_timing_eligible"] is True for p in points), "Incomplete native repeats")
        resources = {}
        for metric in ("wall_seconds", "cpu_seconds", "peak_memory_bytes"):
            values = [p["resources"][metric] for p in points]
            require(all(type(v) in (int, float) and math.isfinite(v) and v >= 0 for v in values), "Invalid resource value")
            resources[metric] = {"median": median(values), "minimum": min(values), "maximum": max(values)}
        exact = sum(p["partition"]["partition_equal"] for p in points)
        require(resources == row["resources"] and exact == row["exact_partition_repeats"], "Native summary differs")
        summaries.append({"cell": cell, "resources": resources, "exact_partition_repeats": exact})
    require(sum(r["repeats"] for r in native["summaries"]) == len(native["points"]), "Unaccounted native point")
    return {"summaries": summaries, "new_measurements": False, "original_cached_costs_established": False}


def register(qfo, all_tools, factorial):
    old = {(r["dataset"], r["key"]): r for r in all_tools["rows"]}
    require(len(old) == len(all_tools["rows"]) == 24 and len(qfo["rows"]) == 24, "Duplicate/incomplete method register")
    require({(r["dataset"], r["key"]) for r in qfo["rows"]} == set(old), "Method register identity set differs")
    selected = []
    intervals = []
    for row in qfo["rows"]:
        key = row["dataset"], row["key"]
        require({k: v for k, v in row.items() if k != "qfo_stage_provenance"} == old[key], "Prior register changed")
        if "qfo_stage_provenance" not in row:
            continue
        stage = row["qfo_stage_provenance"]
        require(all(stage[k] is None for k in ("full_pipeline_wall_s", "full_pipeline_cpu_s", "full_pipeline_peak_memory_bytes"))
                and stage["observations_are_not_summed"] is True, "Unsupported complete cost")
        selected.append({"cell": stage["cell"], "full_pipeline_wall_s": None})
        for i, interval in enumerate(stage["resource_intervals"]):
            # Preparation shares a manifest across arms; arm identity is part of its observation key.
            arm = stage["candidate_arm"] if i == 1 else None
            intervals.append({"evidence_sha256": interval["evidence"]["sha256"], "scope": interval["scope"],
                              "candidate_arm": arm, "measurement": interval["measurement"], "memory": interval["memory"]})
    keys = {(i["evidence_sha256"], i["scope"], i["candidate_arm"]) for i in intervals}
    require(len(selected) == 2 and len(intervals) == 5 and len(keys) == 4, "Incorrect stage association scope")
    require(len(factorial["rows"]) == factorial["cells"] == 16 and factorial["full_pipeline_cost_cells_available"] == 0
            and all(row["full_pipeline_wall_s"] is None and row["full_pipeline_cpu_s"] is None
                    and row["full_pipeline_peak_memory_bytes"] is None for row in factorial["rows"]), "Original missing costs changed")
    score_positions = sum(len(r["scores"]) for r in old.values())
    means = {r["key"]: sum(r["scores"].values()) / len(r["scores"]) for r in old.values() if r["dataset"] == "QfO"}
    for row in old.values():
        if row["dataset"] == "QfO":
            require(math.isclose(means[row["key"]], row["secondary_mean"], rel_tol=0, abs_tol=1e-12), "QfO secondary mean differs")
    return {"method_dataset_rows": len(old), "metric_positions": score_positions,
            "secondary_mean_positions": len(means), "secondary_means": means,
            "selected_qfo_cells": selected, "stage_associations": len(intervals), "distinct_stage_observations": len(keys),
            "stage_observations": intervals, "original_full_cost_cells_unavailable": len(factorial["rows"])}


def summarize(inventory, trace, native, qfo, factorial, all_tools):
    return {"exposure": exposure(inventory), "candidate_trace": candidate_trace(trace),
            "native_costs": native_costs(native), "register": register(qfo, all_tools, factorial)}


def run(directory, output):
    if output.exists() or output.is_symlink():
        raise FileExistsError(output)
    reports = {}
    for role, (name, expected) in INPUTS.items():
        data = (directory / name).read_bytes()
        require(hashlib.sha256(data).hexdigest() == expected, "Changed retained report: " + name)
        reports[role] = json.loads(data)
    result = {"schema": "publication_addenda_reporting_replay_v1", "summary": summarize(**reports),
              "input_sha256": {name: expected for name, expected in INPUTS.values()},
              "source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
              "raw_or_historical_paths_accessed": False, "native_or_scoring_repeated": False,
              "publication_ready": False, "scope": "Retained reporting metadata/arithmetic, not raw admissions or independent scientific validation"}
    rows = associations(reports["inventory"])
    require(bool(rows), "No family associations")
    output.mkdir(parents=True)
    (output / "summary.json").write_text(json.dumps(result, indent=2, sort_keys=True, allow_nan=False) + "\n")
    with (output / "family_evidence.tsv").open("x", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows({k: "NA" if v is None else v for k, v in row.items()} for row in rows)
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--inputs", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    result = run(args.inputs.resolve(), args.output.absolute())
    print(json.dumps({"families": result["summary"]["exposure"]["families"],
                      "associations": result["summary"]["exposure"]["associations"],
                      "register_rows": result["summary"]["register"]["method_dataset_rows"]}))
