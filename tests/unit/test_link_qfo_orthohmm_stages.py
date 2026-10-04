"""Require exact selected-output chains without promoting incremental costs."""

from copy import deepcopy
import csv
import json
from pathlib import Path

import pytest

from benchmark_tools import link_qfo_orthohmm_stages as linkage


def pin(name):
    return {"path": "/fixture/" + name, "bytes": 12, "sha256": name * 4}


def documents(phy=False):
    key = "orthohmm_phylogeny_satellite_v2" if phy else "orthohmm_high_sensitivity"
    cell = linkage.CELLS[key]
    seed, candidate, pair = pin("seed"), pin("candidate"), pin("pair")
    admission = pin("native") if phy else pin("candidate_admission")
    selected = {"status": "admitted", "factorial_cell": cell, "conversion": pin("conversion"),
                "admission": pin("score"), "scores": {"GO": .4}, "secondary_mean": .6,
                "prediction_semantics": "native pairs" if phy else "group-derived pairs",
                "submitted_pairs": 123, "retained_pairs": 123, "removed_mapping_pairs": 0,
                "participant": "fixture"}
    row = {"dataset": "QfO", "key": key, "conversion": selected["conversion"],
           "score_admission": selected["admission"], "scores": selected["scores"],
           "secondary_mean": .6, "output_records": [pair], "native_admission": admission,
           "input_records": [pin("fasta")], "commands": [], "resources": []}
    conversion = {"cell": cell, "semantics": selected["prediction_semantics"],
                  "total_pairs": 123, "retained_pairs": 123, "removed_mapping_pairs": 0,
                  "participant": "fixture", "filtered_pairs": pair, "prepared": pin("prepared"),
                  "candidate_admission": pin("candidate_admission"), "input_fastas": row["input_records"],
                  "candidate_partition": candidate, "command": ["python", "conversion.py"],
                  "status": "corrected_factorial_native_pairs_prepared_unscored" if phy else
                            "corrected_factorial_group_pairs_prepared_unscored"}
    audit = {"status": "candidate_arm_content_verified"}
    candidates = {"status": "corrected_qfo_candidates_admitted", "prepared_manifest": pin("prepared"),
                  "candidate_arms": {cell[:-3]: audit}}
    arm = {"content_audit": audit, "seed_partition": seed, "candidate_partition": candidate,
           "incremental_preparation_seconds": 3, "membership_constraints": pin("constraints")}
    prepared = {"status": "corrected_qfo_four_candidate_arms_prepared_unscored",
                "candidate_arms": {cell[:-3]: arm}, "plan": pin("plan"), "admission": pin("replay_admission"),
                "input_fastas": row["input_records"]}
    coverage = [{"label": "strict_profiles_refined", "output": seed}]
    replay_admission = {"status": "corrected_checked_replay_admitted", "plan": pin("plan"),
                        "source_report": pin("replay_report"), "coverage": coverage, "scheduler": {"JobIDRaw": "1"}}
    report = {"status": "corrected_checked_replay_complete_pending_admission", "exit_code": 0,
              "plan": pin("plan"), "coverage": coverage, "worker_command": ["python", "worker.py"],
              "replay": pin("replay")}
    replay = {"stages": coverage, "command": ["python", "replay.py"], "input": {"manifest": pin("checkpoint")},
              "wall_s": 8, "timings": {"profile_s": 5}}
    plan = {"native_command": replay["command"], "checkpoint_manifest": pin("checkpoint"),
            "input_fastas": row["input_records"]}
    timing = {"argv": report["worker_command"], "evidence": pin("time"), "measurement": {
        "exit_status": 0, "elapsed_seconds": 10, "user_seconds": 20, "system_seconds": 1,
        "max_process_rss_kib": 500}}
    costs = {"dataset": "Corrected QfO", "cell": cell, "candidate_preparation_wall_s": 3,
             "full_pipeline_wall_s": None, "reconciliation_measurement_status":
                 "recorded_incremental" if phy else "not_applicable"}
    native = metrics = None
    if phy:
        conversion.update(native_admission=admission, native_input=pin("native_pairs"))
        native = {"status": "corrected_qfo_native_pair_output_verified", "cell": cell,
                  "native_pairs": conversion["native_input"], "candidate_admission": conversion["candidate_admission"],
                  "native_group_integrity": {"status": "native_group_output_verified", "cell": cell,
                      "native_manifest": pin("manifest"), "native_metrics": pin("metrics"),
                      "integrity": {"scheduler": {"JobIDRaw": "2"}}}}
        metrics = {"status": "complete", "metadata": {"cpu_budget": 32},
                   "rss_measurement": "sampled_sum_of_linux_proc_tree_rss",
                   "wall_s": 12, "user_cpu_s": 30, "system_cpu_s": 2,
                   "peak_process_tree_rss_bytes": 1000, "command": ["python", "reconcile.py"],
                   "input": {"candidate_clusters": candidate, "membership_constraints": pin("constraints")},
                   "outputs": {"manifest": pin("manifest")}, "parameters": {"species_tree_mode": "infer"}}
        costs.update(reconciliation_wall_s=12, reconciliation_user_cpu_s=30,
                     reconciliation_system_cpu_s=2, reconciliation_peak_sampled_tree_rss_bytes=1000)
    return [row, selected, conversion, candidates, prepared, replay_admission,
            report, replay, plan, timing, costs, native, metrics]


@pytest.mark.parametrize("phy", [False, True])
def test_link_is_additive_and_does_not_sum_promote_or_impute(phy):
    data = documents(phy)
    before = deepcopy(data)
    result = linkage.link(*data)
    assert data == before
    assert len(result["resource_intervals"]) == (3 if phy else 2)
    assert result["full_pipeline_wall_s"] is None and result["full_pipeline_cpu_s"] is None
    assert result["full_pipeline_peak_memory_bytes"] is None
    assert result["observations_are_not_summed"] is True
    assert all(r["full_inference"] is False and r["independent_repeat"] is False
               for r in result["resource_intervals"])
    assert result["resource_intervals"][0]["memory"]["unit"] == "KiB"
    assert result["resource_intervals"][1]["measurement"]["user_seconds"] is None
    assert result["resource_intervals"][1]["memory"]["value"] is None
    assert result["replay_internal_wall_s"] == 8 != result["resource_intervals"][0]["measurement"]["elapsed_seconds"]
    if phy:
        assert result["resource_intervals"][2]["memory"]["unit"] == "bytes"
    else:
        assert result["reconciliation"] is None


@pytest.mark.parametrize("damage", ["cell", "score", "secondary", "conversion", "pair", "prepared",
    "arm", "replay_status", "replay_failure", "plan", "seed", "internal_command", "checkpoint",
    "inputs", "register_inputs", "timing_command", "timing_failure", "cost", "full_cost", "already",
    "group_candidate", "R_off", "native_pair", "native_cell", "constraints", "metrics_manifest", "RSS", "cpu_budget"])
def test_wrong_source_output_and_cost_bindings_rejected(damage):
    phy = damage in {"native_pair", "native_cell", "constraints", "metrics_manifest", "RSS", "cpu_budget"}
    data = documents(phy)
    row, selected, conversion, candidates, prepared, admission, report, replay, plan, timing, costs, native, metrics = data
    if damage == "cell": selected["factorial_cell"] = "p0_c0_r0"
    elif damage == "score": row["scores"] = {"GO": .41}
    elif damage == "secondary": row["secondary_mean"] = .61
    elif damage == "conversion": row["conversion"] = pin("wrong")
    elif damage == "pair": row["output_records"] = [pin("wrong")]
    elif damage == "prepared": candidates["prepared_manifest"] = pin("wrong")
    elif damage == "arm": prepared["candidate_arms"][selected["factorial_cell"][:-3]]["content_audit"] = {}
    elif damage == "replay_status": admission["status"] = "running"
    elif damage == "replay_failure": report["exit_code"] = 1
    elif damage == "plan": report["plan"] = pin("wrong")
    elif damage == "seed": prepared["candidate_arms"][selected["factorial_cell"][:-3]]["seed_partition"] = pin("wrong")
    elif damage == "internal_command": replay["command"] = ["wrong"]
    elif damage == "checkpoint": plan["checkpoint_manifest"] = pin("wrong")
    elif damage == "inputs": plan["input_fastas"] = [pin("wrong")]
    elif damage == "register_inputs": row["input_records"] = []
    elif damage == "timing_command": timing["argv"] = ["wrong"]
    elif damage == "timing_failure": timing["measurement"]["exit_status"] = 1
    elif damage == "cost": costs["candidate_preparation_wall_s"] = 4
    elif damage == "full_cost": costs["full_pipeline_wall_s"] = 15
    elif damage == "already": row["resources"] = [{"scope": "new"}]
    elif damage == "group_candidate": conversion["candidate_partition"] = pin("wrong")
    elif damage == "R_off": costs["reconciliation_measurement_status"] = "recorded"
    elif damage == "native_pair": native["native_pairs"] = pin("wrong")
    elif damage == "native_cell": native["cell"] = "p0_c0_r1"
    elif damage == "constraints": metrics["input"]["membership_constraints"] = pin("wrong")
    elif damage == "metrics_manifest": metrics["outputs"]["manifest"] = pin("wrong")
    elif damage == "RSS": costs["reconciliation_peak_sampled_tree_rss_bytes"] = 1
    elif damage == "cpu_budget": metrics["metadata"]["cpu_budget"] = 8
    with pytest.raises(ValueError):
        linkage.link(*data)


@pytest.mark.parametrize("value", [True, -1, float("nan"), float("inf"), "5", None])
def test_bad_wall_measurement(value):
    with pytest.raises(ValueError):
        linkage.stage("scope", pin("stage"), value)


def test_existing_output_is_never_overwritten(tmp_path):
    with pytest.raises(FileExistsError):
        linkage.write(tmp_path, tmp_path)
    target = tmp_path / "broken-link"
    target.symlink_to(tmp_path / "missing")
    with pytest.raises(FileExistsError):
        linkage.write(tmp_path, target)


def test_selected_actual_readback():
    repo = Path(__file__).resolve().parents[2]
    out = repo / "benchmark_tools/results/qfo_orthohmm_stage_metadata_20261004"
    if not (out / "register.json").exists():
        pytest.skip("Source is prepared before actual selected export")
    result = json.loads((out / "register.json").read_text())
    earlier = json.loads((repo / "benchmark_tools/results/all_benchmark_metadata_integrated_20261004_v3/register.json").read_text())
    assert [{k: v for k, v in row.items() if k != linkage.ADDED_FIELD} for row in result["rows"]] == earlier["rows"]
    assert result["resource_intervals"][:37] == earlier["resource_intervals"]
    assert result["resource_table_entries"] == 42 and result["recorded_wall_observations"] == 39
    stages = [r[linkage.ADDED_FIELD] for r in result["rows"] if linkage.ADDED_FIELD in r]
    assert len(stages) == 2
    assert stages[0]["resource_intervals"][0] == stages[1]["resource_intervals"][0]
    assert stages[0]["resource_intervals"][0]["measurement"]["elapsed_seconds"] == 2987.81
    assert stages[0]["resource_intervals"][1]["measurement"]["elapsed_seconds"] == 3.035759819991654
    assert stages[1]["resource_intervals"][1]["measurement"]["elapsed_seconds"] == 77.64100861502811
    assert stages[1]["resource_intervals"][2]["measurement"]["elapsed_seconds"] == 6822.365357
    with (out / "resources.tsv").open(newline="") as stream:
        table = list(csv.DictReader(stream, delimiter="\t"))
    assert len(table) == 42
    for rendered, actual in zip(table, result["resource_intervals"]):
        assert rendered == {k: "NA" if v is None else str(v) for k, v in actual.items()}
