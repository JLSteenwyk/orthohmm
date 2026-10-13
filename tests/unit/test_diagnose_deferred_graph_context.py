from copy import deepcopy
import csv
import json
from pathlib import Path

import pytest

from benchmark_tools import diagnose_deferred_graph_context as diagnostic
from benchmark_tools.prepare_ob_candidate_neighborhood import record


def test_exec_memory_counter_is_not_an_inherited_rusage_peak():
    assert diagnostic.exec_peak_rss_kib("VmRSS:\t12 kB\nVmHWM:\t42 kB\n") == 42


@pytest.mark.parametrize("status", ["VmRSS: 12 kB", "VmHWM: -1 kB", "VmHWM: 1 MB", "VmHWM: NaN kB", "VmHWM: 1 kB extra"])
def test_missing_or_invalid_exec_memory_accounting_is_explicit(status):
    with pytest.raises(ValueError):
        diagnostic.exec_peak_rss_kib(status)


def row(left, right, category="TP", stage="predicted", cross="True"):
    return {"left": left, "right": right, "reference_label": "family", "native_gene_left": "sp|" + left,
        "native_gene_right": "sp|" + right, "reference_is_cross_species": cross, "category": category,
        "stage": stage, "initial_hit_directions": "2", "forward_normalized_score": "1.5",
        "reverse_normalized_score": "1.6", "same_root_hog": "False" if stage == "different_candidate_families" else "True"}


def cohort():
    control = [row("A", "B", "FN", "different_candidate_families"), row("C", "D"),
        row("E", "F", "FN", "reference_within_species", "False"),
        row("G", "H", "FN", "reference_within_species", "False")]
    variant = deepcopy(control)
    variant[0].update(category="TP", stage="predicted", same_root_hog="True")
    variant[1].update(category="FN", stage="within_candidate_not_predicted")
    baseline = deepcopy(control)
    for graph in baseline:
        for stage in diagnostic.GRAPHS:
            graph[stage + "_edge"], graph[stage + "_edge_weight"] = "True", "0.7"
        for stage in diagnostic.PARTITIONS:
            graph["same_" + stage] = "True"
        graph["same_root"] = graph["same_root_hog"]
        graph["same_candidate"] = "False" if graph["stage"] == "different_candidate_families" else "True"
    baseline[0]["same_multipass_refined"] = "False"
    return control, variant, baseline


def test_complete_cohort_strata_retain_baseline_context_and_both_views():
    result = diagnostic.diagnose(*cohort())
    assert result["paired_transitions"]["rescued_asserted_pairs"] == 1
    assert result["paired_transitions"]["lost_asserted_pairs"] == 1
    rescued = result["strata"]["transition:FN->TP"]
    assert rescued["pairs"] == 1 and rescued["initial_hit_counts"] == {"2": 1}
    assert rescued["baseline_edges"]["initial/present_True"] == 1
    assert rescued["baseline_partitions"]["candidate/same_False"] == 1
    assert rescued["baseline_partition_boundaries"]["multipass:True/multipass_refined:False"] == 1
    assert result["strata"]["variant:FN/reference_within_species"]["pairs"] == 2
    assert sum(v["pairs"] for k, v in result["strata"].items() if k.startswith("transition:")) == 4
    assert sum(v["pairs"] for k, v in result["strata"].items() if k.startswith("variant:")) == 4


@pytest.mark.parametrize("mutation", ["length", "order", "baseline_category", "reference_label", "alias", "hit", "graph_boolean", "nan", "negative", "absent_weight", "partition_boolean", "root", "candidate", "variant_stage", "hit_count"])
def test_changed_cohorts_malformed_context_and_variant_stage_are_rejected(mutation):
    control, variant, baseline = cohort()
    if mutation == "length":
        baseline.pop()
    elif mutation == "order":
        baseline.reverse()
    elif mutation in ("baseline_category", "reference_label", "alias", "hit"):
        key = {"baseline_category": "category", "reference_label": "reference_label", "alias": "native_gene_left", "hit": "forward_normalized_score"}[mutation]
        baseline[0][key] = "changed"
    elif mutation == "graph_boolean":
        baseline[0]["initial_edge"] = "yes"
    elif mutation in ("nan", "negative"):
        baseline[0]["initial_edge_weight"] = "NaN" if mutation == "nan" else "-1"
    elif mutation == "absent_weight":
        baseline[0]["initial_edge"] = "False"
    elif mutation == "partition_boolean":
        baseline[0]["same_initial"] = "yes"
    elif mutation == "root":
        baseline[0]["same_root"] = "True"
    elif mutation == "candidate":
        baseline[0]["same_candidate"] = "True"
    elif mutation == "variant_stage":
        variant[0]["stage"] = "different_candidate_families"
    else:
        for rows in (control, variant, baseline):
            rows[0]["initial_hit_directions"] = "3"
    with pytest.raises(ValueError):
        diagnostic.diagnose(control, variant, baseline)


def write_json(path, data):
    path.write_text(json.dumps(data))
    return record(path)


def write_table(path, rows):
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    return record(path)


def setup_driver(tmp_path, monkeypatch):
    control, variant, baseline = cohort()
    control_record = write_table(tmp_path / "control.tsv", control)
    variant_record = write_table(tmp_path / "variant.tsv", variant)
    baseline_record = write_table(tmp_path / "baseline.tsv", baseline)
    graph_plan = write_json(tmp_path / "graph_plan.json", {"cohort": control_record})
    pair_plan = write_json(tmp_path / "pair_plan.json", {"original_control_table": control_record,
        "cases": [{"label": label} for label in ("control", "deferred500")]})
    arms = [{"label": label, "status": "complete_reference_trace_verified", "table": table,
        "raw_categories": {"counts": {"TP": 1, "FN": 3, "FP": 0}}}
        for label, table in (("control", control_record), ("deferred500", variant_record))]
    graph = write_json(tmp_path / "graph_report.json", {"status": "admitted_native_graph_partition_stages_traced",
        "checked": [graph_plan], "table": baseline_record, "reference_pairs": 4})
    pair = write_json(tmp_path / "pair_report.json", {"status": "complete_paired_reference_localization_verified",
        "execution_plan": pair_plan, "arms": arms, "paired_transitions": diagnostic.paired_transitions(control, variant)})
    graph_receipt = tmp_path / "graph_receipt.json"
    pair_receipt = tmp_path / "pair_receipt.json"
    write_json(graph_receipt, {"status": "complete_native_graph_stage_table_readback_verified", "checked": [graph]})
    write_json(pair_receipt, {"status": "independently_recounted_paired_tables_and_raw_vgnc_verified", "checked": [pair],
        "arms": [{"label": a["label"], "raw_categories": a["raw_categories"]["counts"]} for a in arms]})
    monkeypatch.setattr(diagnostic, "EXPECTED_PAIRS", 4)
    monkeypatch.setattr(diagnostic.subprocess, "check_output", lambda *a, **k: "fixture\n")
    monkeypatch.setattr(diagnostic, "exec_peak_rss_kib", lambda status: 1024)
    return Path(diagnostic.__file__).resolve().parent.parent, graph_receipt, pair_receipt


def test_driver_rechecks_inputs_writes_new_receipt_and_never_overwrites(tmp_path, monkeypatch):
    source, graph, pair = setup_driver(tmp_path, monkeypatch)
    before = {p: p.read_bytes() for p in tmp_path.iterdir()}
    output = tmp_path / "result.json"
    result = diagnostic.run(source, graph, pair, output)
    assert result["status"] == "complete_paired_baseline_graph_context_verified"
    assert result["defaults_promoted"] is False and result["goal_completed"] is False
    assert result["paired_transitions"]["asserted_pairs"] == 4
    assert all(p.read_bytes() == data for p, data in before.items())
    with pytest.raises(FileExistsError):
        diagnostic.run(source, graph, pair, output)


@pytest.mark.parametrize("mutation", ["graph_receipt", "paired_receipt", "producer", "order", "native_counts", "control_hash", "file_changed", "full_cohort", "source_root", "memory", "wall"])
def test_driver_rejects_unadmitted_changed_or_incomplete_inputs_before_output(tmp_path, monkeypatch, mutation):
    source, graph, pair = setup_driver(tmp_path, monkeypatch)
    receipt = json.loads(pair.read_text())
    if mutation == "graph_receipt":
        payload = json.loads(graph.read_text())
        payload["status"] = "unverified"
        write_json(graph, payload)
    elif mutation == "paired_receipt":
        receipt["status"] = "unverified"
    elif mutation in ("producer", "order", "native_counts", "control_hash"):
        path = Path(receipt["checked"][0]["path"])
        producer = json.loads(path.read_text())
        if mutation == "producer":
            producer["status"] = "running"
        elif mutation == "order":
            producer["arms"].reverse()
        elif mutation == "native_counts":
            producer["arms"][0]["raw_categories"]["counts"]["TP"] += 1
        else:
            producer["arms"][0]["table"] = producer["arms"][1]["table"]
        receipt["checked"][0] = write_json(path, producer)
    elif mutation == "file_changed":
        (tmp_path / "variant.tsv").write_text("changed")
    elif mutation == "full_cohort":
        monkeypatch.setattr(diagnostic, "EXPECTED_PAIRS", 5)
    elif mutation == "memory":
        monkeypatch.setattr(diagnostic, "exec_peak_rss_kib", lambda status: 2 * 1024 ** 2 + 1)
    elif mutation == "wall":
        clock = iter((0, 601))
        monkeypatch.setattr(diagnostic.time, "monotonic", lambda: next(clock))
    else:
        source = tmp_path
    write_json(pair, receipt)
    with pytest.raises(ValueError):
        diagnostic.run(source, graph, pair, tmp_path / "result.json")
    assert not (tmp_path / "result.json").exists()
