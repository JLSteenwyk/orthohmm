from copy import deepcopy
import csv
import json
from pathlib import Path

import numpy as np
import pytest

from benchmark_tools import report_profile_evidence_pilot as readback
from benchmark_tools import run_profile_weighting_pilot as pilot
from benchmark_tools.bootstrap_orthobench import METRICS, statistics, weighted_records
from benchmark_tools.metric_frontier import ENDPOINTS
from benchmark_tools.prepare_ob_candidate_neighborhood import record


def completed(experiment="weighting"):
    arms, _ = pilot.experiment_settings(experiment)
    plan = {"arms": list(arms), "cpu": 16, "source_revision": "frozen", "experiment": experiment}
    report = {"status": "finished_all_arms_verified", "job_id": "24477", "source_revision": "frozen",
              "defaults_promoted": False, "goal_completed": False,
              "arms": [{"label": mode, "dataset": "OrthoBench", "status": "independent_content_readback_verified"} for mode in arms]}
    if experiment != "weighting":
        report["experiment"] = experiment
    return report, plan, arms


@pytest.mark.parametrize("experiment", ["weighting", "alignment", "selective_hmmer"])
def test_full_ordered_terminal_producer_required_and_legacy_weighting_header_supported(experiment):
    report, plan, arms = completed(experiment)
    readback.require_completed(report, plan, arms, experiment, "24477")


@pytest.mark.parametrize("kind", ["partial", "missing_arm", "reordered", "failed_arm", "wrong_dataset", "source", "job", "budget", "promoted", "error"])
def test_partial_or_changed_pilot_cannot_enter_admission(kind):
    report, plan, arms = completed()
    if kind == "partial":
        report["status"] = "running"
    elif kind == "missing_arm":
        report["arms"].pop()
    elif kind == "reordered":
        report["arms"].reverse()
    elif kind == "failed_arm":
        report["arms"][1]["status"] = "failed"
    elif kind == "wrong_dataset":
        report["arms"][1]["dataset"] = "QfO"
    elif kind == "source":
        report["source_revision"] = "other"
    elif kind == "job":
        report["job_id"] = "other"
    elif kind == "budget":
        plan["cpu"] = 32
    elif kind == "promoted":
        report["defaults_promoted"] = True
    else:
        report["arms"][1]["error"] = "retained failure"
    with pytest.raises(ValueError, match="complete unchanged"):
        readback.require_completed(report, plan, arms, "weighting", "24477")


@pytest.mark.parametrize("accounting", ["24477|COMPLETED|0:0|16|128G|bizon|01:33:31", "24477|RUNNING|0:0|16|128G|bizon|01:33:31",
    "24477|COMPLETED|1:0|16|128G|bizon|01:33:31", "24477|COMPLETED|0:0|32|128G|bizon|01:33:31"])
def test_scheduler_terminal_resources_checked_independently(monkeypatch, accounting):
    text = "JobIDRaw|State|ExitCode|AllocCPUS|ReqMem|NodeList|Elapsed\n" + accounting + "\n"
    monkeypatch.setattr(readback.subprocess, "check_output", lambda *args, **kwargs: text)
    if "COMPLETED|0:0|16|" in accounting:
        assert readback.terminal("24477")["accounting"]["Elapsed"] == "01:33:31"
    else:
        with pytest.raises(ValueError, match="terminal"):
            readback.terminal("24477")


def timing(folder, stage, command):
    (folder / (stage + ".log")).write_text("complete\n")
    (folder / (stage + ".time.txt")).write_text(
        "Elapsed (wall clock) time (h:mm:ss or m:ss): 0:01.50\nUser time (seconds): 1.0\n"
        "System time (seconds): 0.1\nMaximum resident set size (kbytes): 2000\nExit status: 0\n")
    return {"returncode": 0, "command": command, "log": record(folder / (stage + ".log")),
            "time": record(folder / (stage + ".time.txt")),
            "resources": readback.verbose_time((folder / (stage + ".time.txt")).read_text())}


@pytest.mark.parametrize("kind", [None, "command", "resources", "log", "returncode"])
def test_resource_receipts_come_from_exact_raw_logs_not_allocations(tmp_path, kind):
    command = ["python", "inference"]
    execution = timing(tmp_path, "graph", command)
    if kind == "command":
        execution["command"] = ["other"]
    elif kind == "resources":
        execution["resources"]["max_process_rss_kib"] = 128 * 1024 ** 2
    elif kind == "log":
        (tmp_path / "graph.log").write_text("changed")
    elif kind == "returncode":
        execution["returncode"] = 1
    if kind:
        with pytest.raises(ValueError, match="resource record or raw log"):
            readback.execution_receipt(execution, command, tmp_path, "graph")
    else:
        assert len(readback.execution_receipt(execution, command, tmp_path, "graph")) == 2


def arm_fixture(tmp_path, monkeypatch, mode="none", *, experiment="weighting", activation=None):
    source, root = tmp_path / "source", tmp_path / "producer"
    folder = root / mode
    working = folder / "candidate/orthohmm_working_res"
    working.mkdir(parents=True)
    seed = folder / "seed.txt"
    seed.write_text("a b\n")
    candidate, constraints = working / "orthohmm_edges_clustered.txt", working / "phylogeny_candidate_merges.json"
    candidate.write_text("a b\n")
    constraints.write_text("[]\n")
    prediction = folder / "native/orthohmm_working_res/orthohmm_edges_clustered.txt"
    prediction.parent.mkdir(parents=True)
    prediction.write_text("a b\n")
    phylo = folder / "native/orthohmm_phylogeny"
    phylo.mkdir()
    (phylo / "product.tsv").write_text("fixture\n")
    case = {"checkpoint": "/hits", "checkpoint_manifest": {"sha256": "pin"}, "fasta_directory": "/proteomes",
            "baseline_output": "/old/native", "aligner": "mafft", "tree_builder": "FastTree",
            "candidate": record(candidate), "constraints": record(constraints),
            "golden": {"product.tsv": record(phylo / "product.tsv")}, "expected": {"f1": .7, "precision": .8, "recall": .6}}
    plan = {"source_root": str(source), "output": str(root), "cpu": 16, "source_revision": "frozen", "case": case}
    if experiment != "weighting": plan.update(experiment=experiment, hmmer_tools={})
    joint = experiment == "weighted_projection"
    if joint:
        from benchmark_tools.weighted_projection_experiment import configuration
        arm = configuration(mode)
        plan["constituent_controls"] = {mode: case} if arm["role"] != "candidate" else {}
    weighting = arm["profile_sequence_weighting"] if joint else mode if experiment == "weighting" else "none"
    graph_command, native_command = pilot.arm_commands(source, plan, folder, mode)
    (folder / "preflight.json").write_text(json.dumps({"graph": graph_command, "native": native_command}))
    graph = {"parameters": {"profile_sequence_weighting": weighting, "profile_iterations": 1, "profile_expansion": True,
        "matrix": "BLOSUM62", "cpm_resolution": .1, "leiden_seed": 4, "jackknife_profile_thresholds": False,
        "jackknife_single_copy_profiles": False, "profile_min_species": 1, "singleton_tie_policy": "first_hit"},
        "counts": {"profiles_built": 1, "strict_profile_edges": 1}, "profile_iterations": [{"profiles_built": 1,
            "profile_effective_counts": [] if weighting == "none" else [{"cluster_id": 0, "effective_sequence_count": 3.}]}],
        "stages": [{"label": "fixture", "output": record(seed)}]}
    if joint:
        graph["parameters"]["profile_projection_policy"] = arm["profile_projection_policy"]
        if arm["role"] == "candidate":
            graph["profile_iterations"][0]["profile_construction_policies"] = [
                {**r, "sequence_weighting": weighting, "projection_policy": arm["profile_projection_policy"]}
                for r in graph["profile_iterations"][0]["profile_effective_counts"]]
    if experiment == "selective_hmmer": graph["parameters"]["selective_profile_mode"] = mode
    (folder / "graph_metrics.json").write_text(json.dumps(graph))
    details = {"precision": .8, "recall": .6, "family_details": {"fixture": True}}
    row = {"label": mode, "graph_execution": timing(folder, "graph", graph_command),
        "native_execution": timing(folder, "native", native_command), "graph_report": record(folder / "graph_metrics.json"),
        "activation": pilot.validate_activation(graph, mode, experiment) if activation is None else activation, "scores": {"OrthoBench F1": .7}, "details": details,
        "metric_vector": {ENDPOINTS[0]: .7}, "native_summary": {"fixture": True},
        "coverage": readback.validate_partition(candidate, prediction),
        "candidate": {"candidate_partition": record(candidate), "membership_constraints": record(constraints), "seed_partition": record(seed)},
        "candidate_control_matches": {"candidate_partition": True, "membership_constraints": True},
        "golden_matches": {"product.tsv": True}, "candidate_resources": {}, "total_incremental_wall_s": 2.}
    native = {"genes": 251378, "species": 12, "groups": 1, "pairs": 1, "summary": row["native_summary"],
              "phylogeny": str(phylo), "products": {}, "checked": []}
    if joint:
        (phylo / "species_tree.rooted.nwk").write_text("(a,b);\n")
        row.update(rooted_species_tree=record(phylo / "species_tree.rooted.nwk"), arm_parameters=arm,
                   graph_control_matches=True)
        case["coverage"] = {k: native[k] for k in ("genes", "species", "groups", "pairs")}
        case["graph_stages"] = [{"label": "fixture", **{k: record(seed)[k] for k in ("bytes", "sha256")}}]
        if mode != "legacy_none":
            row["native_pair_turnover"] = {"files": [], "added_pairs": 0, "lost_pairs": 0}
            monkeypatch.setattr(readback, "native_pair_turnover", lambda *a: row["native_pair_turnover"])
    monkeypatch.setattr(readback, "replay_native", lambda *args: native)
    monkeypatch.setattr(readback, "readback_scores", lambda *args: (row["scores"], details, []))
    return plan, row, native, details, folder


@pytest.mark.parametrize("mode", ["none", "henikoff", "henikoff_cluster", "species_henikoff_cluster"])
def test_arm_admission_revalidates_planned_commands_activation_and_outputs(tmp_path, monkeypatch, mode):
    plan, row, _, _, _ = arm_fixture(tmp_path, monkeypatch, mode)
    result = readback.admit_arm(tmp_path, plan, row, "weighting", {"a", "b"})
    assert result["status"] == "independent_content_readback_verified"
    assert result["source_revision"] == "frozen" and result["scores"] == {ENDPOINTS[0]: .7}
    assert result["activation"]["mode"] == mode and result["checked"]


@pytest.mark.parametrize("kind", ["preflight", "activation", "graph_pin", "candidate_path", "summary", "universe", "coverage", "details", "metric_vector"])
def test_claimed_arm_fields_cannot_override_raw_readback(tmp_path, monkeypatch, kind):
    plan, row, native, details, folder = arm_fixture(tmp_path, monkeypatch)
    if kind == "preflight":
        (folder / "preflight.json").write_text("{}")
    elif kind == "activation":
        row["activation"]["profiles_built"] = 0
    elif kind == "graph_pin":
        (folder / "graph_metrics.json").write_text("{}")
    elif kind == "candidate_path":
        row["candidate"]["candidate_partition"]["path"] = str(tmp_path / "other")
    elif kind == "summary":
        row["native_summary"] = {}
    elif kind == "universe":
        native["genes"] = 251377
    elif kind == "coverage":
        row["coverage"]["prediction_genes"] = 1
    elif kind == "details":
        row["details"] = {**details, "precision": .9}
    else:
        row["metric_vector"] = {ENDPOINTS[0]: .8}
    with pytest.raises(ValueError):
        readback.admit_arm(tmp_path, plan, row, "weighting", {"a", "b"})


def family_cells():
    records = [{"refog": f"RefOG{index:03d}", "genes": 3, "true_positive": 2,
                "false_positive": 1, "false_negative": 1} for index in range(70)]
    cells = []
    for index, label in enumerate(pilot.ARMS):
        data = deepcopy(records)
        if index:
            data[0].update(true_positive=3, false_negative=0)
        _, _, counts = weighted_records(data)
        aggregate = dict(zip(METRICS, statistics(counts.sum(axis=0)).tolist()))
        cells.append({"label": label, "details": {"family_details": {**aggregate, "refog_records": data}}})
    return cells


def test_real_paired_bootstrap_includes_all_variants_and_does_not_change_point_gates():
    result = readback.uncertainty(family_cells())
    assert result["replicates"] == 100000 and result["seed"] == 20261018
    assert set(result["comparisons"]) == set(pilot.ARMS[1:])
    assert "9 planned endpoints" in result["multiplicity"]
    assert "Post hoc" in result["scope"] and len(result["families"]) == 70


def test_inconsistent_family_counts_do_not_enter_uncertainty():
    cells = family_cells()
    cells[0]["details"]["family_details"]["f_score"] = 99
    with pytest.raises(ValueError, match="sufficient statistics"):
        readback.uncertainty(cells)


def test_table_never_invents_missing_full_vector_endpoints(tmp_path):
    path = tmp_path / "scores.tsv"
    context = [{"label": "OrthoFinder full", "status": "retained", "scores": dict.fromkeys(ENDPOINTS, .8)}]
    cells = [{"label": "henikoff", "status": "independent_content_readback_verified", "scores": {ENDPOINTS[0]: .74}}]
    readback.write_table(path, cells, context)
    with path.open() as stream:
        rows = list(csv.DictReader(stream, delimiter="\t"))
    assert rows[0][ENDPOINTS[-1]] == "0.8"
    assert rows[1][ENDPOINTS[0]] == "0.74" and all(rows[1][key] == "NA" for key in ENDPOINTS[1:])
    with pytest.raises(FileExistsError):
        readback.write_table(path, cells, context)


def test_existing_readback_is_never_overwritten(tmp_path):
    with pytest.raises(FileExistsError, match="Preserve"):
        readback.run(tmp_path, "weighting", tmp_path)


@pytest.mark.parametrize("kind", [None, "changed_snapshot", "invalid_score", "missing_endpoint", "duplicate_method", "missing_population"])
def test_retained_context_is_immutable_complete_and_keeps_assessment_populations(tmp_path, monkeypatch, kind):
    score_table = tmp_path / "old_scores.tsv"
    score_table.write_text("retained fixture\n")
    methods = [{"key": key, "label": key, "scores": dict.fromkeys(ENDPOINTS, .8), "cells": [
        {"dataset": dataset, "status": "retained_registered_result", "retained_details": {"assessed_relations": 42}}
        for dataset in ("OrthoBench", "QfO", "ThreeKingdoms")]} for key in readback.RETAINED_KEYS]
    if kind == "invalid_score":
        methods[0]["scores"][ENDPOINTS[0]] = True
    elif kind == "missing_endpoint":
        methods[0]["scores"].pop(ENDPOINTS[0])
    elif kind == "duplicate_method":
        methods.append(deepcopy(methods[0]))
    elif kind == "missing_population":
        methods[0]["cells"].pop()
    payload = json.dumps({"methods": methods, "score_table": record(score_table)}).encode()
    path = tmp_path / readback.ADMISSION
    path.parent.mkdir(parents=True)
    path.write_bytes(payload + (b"\n" if kind == "changed_snapshot" else b""))
    monkeypatch.setattr(readback.subprocess, "check_output", lambda *args, **kwargs: payload)
    if kind:
        with pytest.raises(ValueError):
            readback.retained_context(tmp_path)
    else:
        rows, pins = readback.retained_context(tmp_path)
        assert len(rows) == 9 and len(pins) == 2
        assert all(row["status"] == "retained_context_not_reaudited_here" for row in rows)
        assert rows[0]["population_details"][1]["retained_details"]["assessed_relations"] == 42
