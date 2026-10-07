"""Candidate-only VGNC diagnostics preserve the original baseline audit."""

import gzip
import json
from pathlib import Path

import pytest

from benchmark_tools import export_native_qfo_candidate_vgnc as exporter
from benchmark_tools import readback_native_qfo_candidate_vgnc as reader
from tests.unit.test_native_qfo_vgnc_blocks import inputs as historical_inputs, run as run_historical, dump, compressed


@pytest.fixture
def inputs(historical_inputs, tmp_path):
    previous = run_historical(historical_inputs)
    baseline_ref = exporter.old.record(historical_inputs["output"] / "report.json")
    readback = reader.independent.review(baseline_ref["path"], baseline_ref["sha256"])
    readback_ref = dump(tmp_path / "baseline_readback.json", readback)
    snapshot = json.loads(Path(historical_inputs["snapshot"]["path"]).read_text())
    baseline = snapshot["rows"][0]
    admission = json.loads(Path(baseline["admission"]["path"]).read_text())
    root = tmp_path / "candidate"
    (root / "other").mkdir(parents=True)
    vg = root / "results/VGNC"
    vg.mkdir(parents=True)
    participant = "native_p0_c1_r0"
    database = root / "other" / (participant + ".db")
    database.write_bytes(Path(previous["methods"][0]["database"]["path"]).read_bytes())
    database_ref = exporter.old.record(database)
    with gzip.open(previous["methods"][0]["raw"]["path"], "rt") as stream:
        lines = stream.read().replace("g4\tg5\tFN", "g4\tg5\tTP")
    lines += "g1\tg5\tFP\tB\tC\ts1\ts2\n"
    raw = compressed(vg / ("VGNC_" + participant.replace("_", "-") + "_raw.txt.gz"), lines)
    counts = dict(TP=4, FP=4, FN=0)
    ratios = exporter.old.ratios(counts)
    native = dict(participant_id=participant, metric_x=ratios["recall"], metric_y=ratios["precision"], stderr_x=0, stderr_y=0)
    metric = dump(vg / "VGNC.json", dict(datalink=dict(inline_data=dict(
        visualization=dict(x_axis="TPR", y_axis="PPV"), challenge_participants=[native]))))
    execution = dump(root / "execution.json", dict(exit_code=0, native_index=8, cell="p0_c1_r0",
        status="process_succeeded_pending_independent_admission", outputs=[database_ref, metric, raw]))
    admission.update(native_index=8, cell="p0_c1_r0", participant=participant, execution_report=execution,
        checked_records=[execution, database_ref, metric, raw], assessment=dict(endpoints=dict(VGNC=dict(native_participant=native))))
    admission_ref = dump(root / "admission.json", admission)
    candidate = dict(baseline, index=8, cell="p0_c1_r0", participant=participant, admission=admission_ref,
        endpoint_details=dict(VGNC=dict(precision=ratios["precision"], recall=ratios["recall"])), scores=dict(VGNC=ratios["f1"]))
    snapshot["rows"] = snapshot["rows"][:2] + [candidate] + [dict(index=i, accuracy_admitted=False) for i in range(9, 13)]
    snapshot["evidence"].append(admission_ref)
    snapshot_ref = dump(tmp_path / "candidate_snapshot.json", snapshot)
    return dict(snapshot=snapshot_ref, baseline=baseline_ref, readback=readback_ref, output=tmp_path / "candidate_out")


def run(inputs):
    args = [value for key in ("snapshot", "baseline", "readback") for value in (inputs[key]["path"], inputs[key]["sha256"])]
    return exporter.export(*args, inputs["output"])


def test_complete_candidate_export_and_independent_reader(inputs):
    report = run(inputs)
    candidate = report["candidate"]
    assert candidate["counts"] == dict(TP=4, FP=4, FN=0)
    assert candidate["validation"]["tp_fp_overlap"] == 1
    assert candidate["alias_rows"] == dict(alias=1, identical=1)
    assert candidate["within_block_false_positives"] == 2
    assert candidate["cross_block_false_positives"] == 2
    assert dict(baseline="FN", candidate="TP", pairs=1) in report["transition_counts"]
    assert dict(baseline="not_scored", candidate="FP", pairs=1) in report["transition_counts"]
    assert dict(baseline="TP+FP", candidate="TP+FP", pairs=1) in report["transition_counts"]
    assert report["baseline_raw_audit_repeated"] is False
    ref = exporter.old.record(inputs["output"] / "report.json")
    result = reader.review(ref["path"], ref["sha256"])
    assert result["candidate_raw_rows"] == 8
    assert result["cached_baseline_union_rows"] == 7
    assert result["transition_pairs"] == 7
    assert result["primary_exporter_imported"] is False
    assert result["uncertainty_admitted"] is False
    assert result["failed_r1_timing_remains_ineligible"] is True
    with pytest.raises(ValueError, match="fresh direct"):
        run(inputs)


@pytest.mark.parametrize("fault", ("snapshot_hash", "baseline_hash", "readback_hash", "scope", "cohort",
    "duplicate_identity", "timing", "baseline_admission", "baseline_metrics", "readback_link", "readback_scope",
    "candidate_semantics", "admission", "execution", "inventory", "raw_identity", "database_identity", "aggregate_identity"))
def test_export_refuses_changed_bound_evidence(inputs, fault):
    if fault.endswith("_hash"):
        inputs[fault.removesuffix("_hash")]["sha256"] = "0" * 64
    else:
        snapshot = json.loads(Path(inputs["snapshot"]["path"]).read_text())
        row = snapshot["rows"][2]
        if fault in ("scope", "cohort", "duplicate_identity", "timing", "baseline_admission", "baseline_metrics", "candidate_semantics"):
            if fault == "scope": snapshot["publication_ready"] = True
            elif fault == "cohort": snapshot["rows"][3]["accuracy_admitted"] = True
            elif fault == "duplicate_identity": snapshot["rows"][3]["index"] = 8
            elif fault == "timing": snapshot["rows"][1]["timing_eligible"] = True
            elif fault == "baseline_admission": snapshot["rows"][0]["admission"] = row["admission"]
            elif fault == "baseline_metrics": snapshot["rows"][0]["scores"]["VGNC"] += .1
            else: row["prediction_semantics"] = "tree-derived pairs"
            inputs["snapshot"] = dump(Path(inputs["snapshot"]["path"]), snapshot)
        elif fault in ("readback_link", "readback_scope"):
            value = json.loads(Path(inputs["readback"]["path"]).read_text())
            if fault == "readback_link": value["report"] = inputs["snapshot"]
            else: value["uncertainty_admitted"] = True
            inputs["readback"] = dump(Path(inputs["readback"]["path"]), value)
        else:
            admission = json.loads(Path(row["admission"]["path"]).read_text())
            if fault.endswith("_identity"):
                suffix = {"raw_identity": "raw.txt.gz", "database_identity": ".db", "aggregate_identity": "VGNC.json"}[fault]
                ref = exporter.old.unique_record(admission["checked_records"], suffix)
                with Path(ref["path"]).open("ab") as stream: stream.write(b"changed")
            else:
                if fault == "admission": admission["native_index"] = 9
                elif fault == "inventory": admission["checked_records"] = []
                else:
                    value = json.loads(Path(admission["execution_report"]["path"]).read_text())
                    value["exit_code"] = 1
                    old_ref = admission["execution_report"]
                    admission["execution_report"] = dump(Path(old_ref["path"]), value)
                    admission["checked_records"] = [admission["execution_report"] if r == old_ref else r for r in admission["checked_records"]]
                row["admission"] = dump(Path(row["admission"]["path"]), admission)
                snapshot["evidence"] = [r for r in snapshot["evidence"] if r["path"] != row["admission"]["path"]] + [row["admission"]]
                inputs["snapshot"] = dump(Path(inputs["snapshot"]["path"]), snapshot)
    with pytest.raises((ValueError, KeyError)):
        run(inputs)
    assert not inputs["output"].exists()


@pytest.mark.parametrize("fault", ("scope", "source", "counts", "metrics", "transition_count", "union",
    "block_summary", "alias", "mapping", "input_binding", "baseline", "table", "raw", "report_hash"))
def test_reader_refuses_corruption(inputs, fault):
    report = run(inputs)
    path = inputs["output"] / "report.json"
    if fault == "scope": report["uncertainty_admitted"] = True
    elif fault == "source": report["source"]["sha256"] = "0" * 64
    elif fault == "counts": report["candidate"]["counts"]["FP"] += 1
    elif fault == "metrics": report["differences"]["f1"] += .1
    elif fault == "transition_count": report["transition_counts"][0]["pairs"] += 1
    elif fault == "union": report["union_scored_pairs"] += 1
    elif fault == "block_summary": report["candidate"]["within_block_false_positives"] += 1
    elif fault == "alias": report["candidate"]["alias_rows"]["alias"] += 1
    elif fault == "mapping": report["candidate"]["selected_reference_mapping_sha256"] = "0" * 64
    elif fault == "input_binding": report["checked_records"] = []
    elif fault == "baseline": report["baseline"]["counts"]["TP"] += 1
    elif fault in ("table", "raw"):
        ref = report["candidate"]["table" if fault == "table" else "raw"]
        with Path(ref["path"]).open("ab") as stream: stream.write(b"changed")
    updated = dump(path, report)
    if fault == "report_hash": updated["sha256"] = "0" * 64
    with pytest.raises((ValueError, KeyError)):
        reader.review(updated["path"], updated["sha256"])


@pytest.mark.parametrize("value", ("", "FP+TP", "TP+TP", "TN", "TP+FN+TN", "not_scored+TP"))
def test_noncanonical_states_are_rejected(value):
    with pytest.raises(ValueError): exporter.categories(value)
    with pytest.raises(ValueError): reader.state(value)


def test_preserves_historical_kernels_and_frozen_route():
    assert exporter.old.__name__.endswith("export_native_qfo_vgnc_blocks")
    assert reader.independent.__name__.endswith("readback_native_qfo_vgnc_blocks")
    assert "export_native_qfo_candidate_vgnc" not in reader.__dict__
    assert "bootstrap" not in reader.__dict__


@pytest.mark.parametrize("fault", ("duplicate", "annotation", "unasserted", "ineligible", "missing_truth"))
def test_rebound_raw_still_requires_native_semantics(inputs, fault):
    snapshot = json.loads(Path(inputs["snapshot"]["path"]).read_text())
    row = snapshot["rows"][2]
    admission = json.loads(Path(row["admission"]["path"]).read_text())
    ref = exporter.old.unique_record(admission["checked_records"], "raw.txt.gz")
    with gzip.open(ref["path"], "rt") as stream:
        lines = stream.read().splitlines()
    if fault == "duplicate": lines.append(lines[0])
    elif fault == "annotation": lines[0] = lines[0].replace("\tB\tA\t", "\tA\tA\t")
    elif fault == "unasserted": lines[0] = "g1\tg6\tTP\tB\tA\ts1\ts1"
    elif fault == "ineligible": lines[-1] = "g1\tg3\tFP\tB\tB\ts1\ts2"
    else: lines = [line for line in lines if not line.startswith("g4\tg5\t")]
    raw_ref = compressed(Path(ref["path"]), "\n".join(lines) + "\n")
    execution = json.loads(Path(admission["execution_report"]["path"]).read_text())
    execution["outputs"] = [raw_ref if r == ref else r for r in execution["outputs"]]
    old_execution = admission["execution_report"]
    admission["execution_report"] = dump(Path(old_execution["path"]), execution)
    admission["checked_records"] = [raw_ref if r == ref else admission["execution_report"] if r == old_execution else r
                                    for r in admission["checked_records"]]
    old_admission = row["admission"]
    row["admission"] = dump(Path(old_admission["path"]), admission)
    snapshot["evidence"] = [row["admission"] if r == old_admission else r for r in snapshot["evidence"]]
    inputs["snapshot"] = dump(Path(inputs["snapshot"]["path"]), snapshot)
    with pytest.raises(ValueError): run(inputs)
    assert not inputs["output"].exists()


def test_database_mutation_during_mapping_is_retained_as_failure(inputs, monkeypatch):
    original = exporter.old.mapped_reference
    def mutate(database, *args):
        result = original(database, *args)
        with database.open("ab") as stream: stream.write(b"changed during mapping")
        return result
    monkeypatch.setattr(exporter.old, "mapped_reference", mutate)
    with pytest.raises(ValueError, match="identity changed"):
        run(inputs)
    assert not inputs["output"].exists()
