"""Native VGNC decomposition fixtures; no scientific admission or rescore."""

import gzip
import json
from pathlib import Path
import sqlite3

import pytest

from benchmark_tools import export_native_qfo_vgnc_blocks as exporter
from benchmark_tools import readback_native_qfo_vgnc_blocks as reader


def compressed(path, lines):
    with gzip.open(path, "wt") as stream:
        stream.write(lines)
    return exporter.record(path)


def dump(path, value):
    path.write_text(json.dumps(value))
    return exporter.record(path)


@pytest.fixture
def inputs(tmp_path):
    reference = compressed(tmp_path / "vgnc-orthologs.txt.gz", "1\t2\tA\n2\t6\tA\n1\t3\tB\n4\t5\tC\n")
    truth, labels = exporter.reference_data(Path(reference["path"]))
    mapping, summary = exporter.reference_blocks(truth)
    database = tmp_path / "db.sqlite"
    with sqlite3.connect(database) as connection:
        connection.execute("CREATE TABLE proteomes (prot_nr INTEGER, uniprot_id TEXT, species TEXT)")
        connection.executemany("INSERT INTO proteomes VALUES (?, ?, ?)",
            [(1, "old", "s1"), (1, "g1", "s1"), (1, "g1", "s1"),
             (2, "g2", "s2"), (3, "g3", "s2"), (4, "g4", "s1"), (5, "g5", "s2"), (6, "g6", "s1")])
    _, _, digest, _ = exporter.mapped_reference(database, truth, labels)
    reference_table = tmp_path / "reference_blocks.tsv"
    exporter.write_table(reference_table, ["block", "proteins", "asserted_pairs", "species"],
        [["A", 4, 3, "s1,s2"], ["C", 2, 1, "s1,s2"]])
    historical = dict(status="corrected_vgnc_native_rows_mapped_to_reference_blocks", uncertainty_admitted=False,
        checked_records=[reference], reference=summary, reference_pairs=4,
        reference_table=exporter.record(reference_table), reference_mapping_sha256=digest)
    historical_ref = dump(tmp_path / "historical.json", historical)
    base_rows = [
        "g1\tg2\tTP\tB\tA\ts1\ts2", "g2\tg6\tTP\tA\tA\ts2\ts1",
        "g1\tg3\tTP\tB\tB\ts1\ts2", "g4\tg5\tFN\tC\tC\ts1\ts2",
        "g1\tg2\tFP\tB\tA\ts1\ts2", "g3\tg6\tFP\tB\tA\ts2\ts1", "g5\tg6\tFP\tC\tA\ts2\ts1"]
    other_rows = [base_rows[0], base_rows[1].replace("TP", "FN"), base_rows[2].replace("TP", "FN"),
        base_rows[3], base_rows[4], "g1\tg5\tFP\tB\tC\ts1\ts2"]
    snapshot_rows = []
    for index, cell, lines in zip((6, 7), exporter.CELLS, (base_rows, other_rows)):
        participant = "native_" + cell
        root = tmp_path / cell
        (root / "other").mkdir(parents=True)
        vg = root / "results/VGNC"
        vg.mkdir(parents=True)
        db = root / "other" / (participant + ".db")
        db.write_bytes(database.read_bytes())
        db_ref = exporter.record(db)
        raw = compressed(vg / ("VGNC_" + participant.replace("_", "-") + "_raw.txt.gz"), "\n".join(lines) + "\n")
        counts = dict(TP=3, FP=3, FN=1) if index == 6 else dict(TP=1, FP=2, FN=3)
        ratios = exporter.ratios(counts)
        native = dict(participant_id=participant, metric_x=ratios["recall"], metric_y=ratios["precision"], stderr_x=0, stderr_y=0)
        metric_ref = dump(vg / "VGNC.json", dict(datalink=dict(inline_data=dict(
            visualization=dict(x_axis="TPR", y_axis="PPV"), challenge_participants=[native]))))
        execution_ref = dump(root / "execution.json", dict(exit_code=0, native_index=index, cell=cell,
            status="process_succeeded_pending_independent_admission", outputs=[db_ref, metric_ref, raw]))
        admission = dict(schema="full_native_factorial_qfo_admission_v1" if index == 6 else "measurement_failed_native_qfo_admission_v1",
            status="full_native_factorial_qfo_assessment_admitted" if index == 6 else "measurement_failed_native_qfo_assessment_admitted",
            accuracy_admitted=True, publication_ready=False, native_index=index, cell=cell, participant=participant,
            execution_report=execution_ref, checked_records=[execution_ref, db_ref, metric_ref, raw, metric_ref],
            assessment=dict(endpoints=dict(VGNC=dict(native_participant=native))), resources=None,
            scientific_timings_admitted=False, eligible_for_timing_comparison=False, original_native_scheduler_success=False)
        admission_ref = dump(root / "admission.json", admission)
        snapshot_rows.append(dict(index=index, cell=cell, participant=participant, accuracy_admitted=True,
            admission=admission_ref, prediction_semantics="clique" if index == 6 else "native phylogenetic pairs",
            measurement_status="successful_native_terminal_review" if index == 6 else "failed_timing_scientific_outputs_recovered",
            endpoint_details=dict(VGNC=dict(precision=ratios["precision"], recall=ratios["recall"])),
            scores=dict(VGNC=ratios["f1"]), timing_eligible=False, timing_admitted=False))
    plan_ref = dump(tmp_path / "plan.json", dict(frozen=True))
    snapshot = dict(schema="native_qfo_scientific_reporting_snapshot_v1", new_scoring_or_admission=False,
        publication_ready=False, recovered_inference_resources_admitted=False,
        source=exporter.record(Path(exporter.__file__).with_name("export_native_qfo_scientific_scores.py")),
        plan=plan_ref, rows=snapshot_rows + [dict(accuracy_admitted=False, index=8, cell="unfinished")],
        evidence=[r["admission"] for r in snapshot_rows])
    snapshot_ref = dump(tmp_path / "snapshot.json", snapshot)
    return dict(snapshot=snapshot_ref, historical=historical_ref, output=tmp_path / "out")


def run(inputs):
    return exporter.export(inputs["snapshot"]["path"], inputs["snapshot"]["sha256"],
        inputs["historical"]["path"], inputs["historical"]["sha256"], inputs["output"])


def test_complete_export_and_independent_readback(inputs):
    result = run(inputs)
    assert result["reference"] == dict(family_labels=3, reference_proteins=6, shared_proteins=1,
        reference_blocks=2, merged_label_groups=[["A", "B"]])
    assert result["union_scored_pairs"] == 7
    left, right = result["methods"]
    assert left["counts"] == dict(TP=3, FP=3, FN=1)
    assert left["validation"]["tp_fp_overlap"] == 1
    assert left["within_block_false_positives"] == 2 and left["cross_block_false_positives"] == 1
    assert right["within_block_false_positives"] == 1 and right["cross_block_false_positives"] == 1
    assert left["alias_rows"] == right["alias_rows"] == dict(identical=1, alias=1)
    assert dict(r0="TP+FP", r1="TP+FP", pairs=1) in result["transition_counts"]
    assert dict(r0="not_scored", r1="FP", pairs=1) in result["transition_counts"]
    ref = exporter.record(inputs["output"] / "report.json")
    checked = reader.review(ref["path"], ref["sha256"])
    assert checked["raw_rows"] == 13 and checked["transition_pairs"] == 7
    assert checked["reference_rows"] == 2 and checked["uncertainty_admitted"] is False
    with pytest.raises(ValueError, match="Output already exists"):
        run(inputs)


@pytest.mark.parametrize("fault", ("snapshot_hash", "historical_hash", "snapshot_scope", "snapshot_source", "plan",
    "reference", "reference_summary", "historical_scope", "mapping", "reference_table", "missing_cell", "extra_cell",
    "snapshot_admission", "accuracy", "index", "schema", "status", "timing", "execution_inventory", "execution",
    "database", "raw", "aggregate", "score", "semantics_axes", "native_participant"))
def test_export_refuses_changed_evidence(inputs, fault):
    snapshot = json.loads(Path(inputs["snapshot"]["path"]).read_text())
    historical = json.loads(Path(inputs["historical"]["path"]).read_text())
    row = snapshot["rows"][0]
    admission = json.loads(Path(row["admission"]["path"]).read_text())
    if fault.endswith("_hash"):
        inputs[fault.removesuffix("_hash")]["sha256"] = "0" * 64
    elif fault == "snapshot_scope":
        snapshot["new_scoring_or_admission"] = True
    elif fault in ("snapshot_source", "plan"):
        snapshot["source" if fault == "snapshot_source" else "plan"]["sha256"] = "0" * 64
    elif fault in ("reference", "reference_summary", "historical_scope", "mapping", "reference_table"):
        if fault == "reference":
            historical["checked_records"][0]["sha256"] = "0" * 64
        elif fault == "reference_summary":
            historical["reference"]["reference_blocks"] += 1
        elif fault == "historical_scope":
            historical["uncertainty_admitted"] = True
        elif fault == "mapping":
            historical["reference_mapping_sha256"] = "0" * 64
        else:
            Path(historical["reference_table"]["path"]).write_text("block\tproteins\tasserted_pairs\tspecies\n")
            historical["reference_table"] = exporter.record(historical["reference_table"]["path"])
        inputs["historical"] = dump(Path(inputs["historical"]["path"]), historical)
    elif fault in ("missing_cell", "extra_cell", "snapshot_admission", "score"):
        if fault == "missing_cell":
            snapshot["rows"].pop(1)
        elif fault == "extra_cell":
            snapshot["rows"][2]["accuracy_admitted"] = True
        elif fault == "snapshot_admission":
            snapshot["evidence"] = []
        else:
            row["scores"]["VGNC"] += .01
    elif fault in ("database", "raw", "aggregate"):
        suffix = ".db" if fault == "database" else "raw.txt.gz" if fault == "raw" else "/VGNC.json"
        ref = exporter.unique_record(admission["checked_records"], suffix)
        with Path(ref["path"]).open("ab") as stream:
            stream.write(b"changed")
    elif fault in ("semantics_axes", "native_participant"):
        ref = exporter.unique_record(admission["checked_records"], "/VGNC.json")
        native = json.loads(Path(ref["path"]).read_text())
        inline = native["datalink"]["inline_data"]
        if fault == "semantics_axes":
            inline["visualization"]["x_axis"] = "wrong"
        else:
            inline["challenge_participants"][0]["participant_id"] = "wrong"
        new_ref = dump(Path(ref["path"]), native)
        admission["checked_records"] = [new_ref if r == ref else r for r in admission["checked_records"]]
        execution = json.loads(Path(admission["execution_report"]["path"]).read_text())
        execution["outputs"] = [new_ref if r == ref else r for r in execution["outputs"]]
        old_execution = admission["execution_report"]
        admission["execution_report"] = dump(Path(old_execution["path"]), execution)
        admission["checked_records"] = [admission["execution_report"] if r == old_execution else r for r in admission["checked_records"]]
    else:
        if fault == "timing":
            row = snapshot["rows"][1]
            admission = json.loads(Path(row["admission"]["path"]).read_text())
            row["timing_eligible"] = True
        elif fault == "accuracy":
            admission["accuracy_admitted"] = False
        elif fault == "index":
            admission["native_index"] = 8
        elif fault in ("schema", "status"):
            admission[fault] = "wrong"
        elif fault == "execution_inventory":
            admission["checked_records"] = []
        elif fault == "execution":
            admission["execution_report"]["sha256"] = "0" * 64
    if fault in ("accuracy", "index", "schema", "status", "timing", "execution_inventory", "execution", "semantics_axes", "native_participant"):
        old = row["admission"]
        row["admission"] = dump(Path(old["path"]), admission)
        snapshot["evidence"] = [row["admission"] if r == old else r for r in snapshot["evidence"]]
    if fault != "snapshot_hash":
        inputs["snapshot"] = dump(Path(inputs["snapshot"]["path"]), snapshot)
    with pytest.raises((ValueError, KeyError)):
        run(inputs)


@pytest.mark.parametrize("fault", ("scope", "source", "cohort", "counts", "metrics", "validation", "within_fp", "cross_fp",
    "nonzero", "alias", "mapping", "database_scope", "prediction_scope", "transition_counts", "union", "difference",
    "sparse_table", "reference_table", "transition_table", "input_identity", "report_digest"))
def test_independent_readback_refuses_corruption(inputs, fault):
    report = run(inputs)
    method = report["methods"][0]
    if fault == "scope":
        report["uncertainty_admitted"] = True
    elif fault == "source":
        report["source"]["sha256"] = "0" * 64
    elif fault == "cohort":
        report["methods"].reverse()
    elif fault in ("counts", "metrics", "validation"):
        target = method[fault]
        if fault == "validation":
            target["counts"]["TP"] += 1
        else:
            target["TP" if fault == "counts" else "f1"] += 1
    elif fault in ("within_fp", "cross_fp", "nonzero"):
        method[{"within_fp": "within_block_false_positives", "cross_fp": "cross_block_false_positives", "nonzero": "nonzero_cells"}[fault]] += 1
    elif fault == "alias":
        method["alias_rows"]["alias"] += 1
    elif fault == "mapping":
        method["selected_reference_mapping_sha256"] = "0" * 64
    elif fault in ("database_scope", "prediction_scope"):
        method["full_database_hash_checked" if fault == "database_scope" else "prediction_edges_requeried"] = fault != "database_scope"
    elif fault == "transition_counts":
        report["transition_counts"][0]["pairs"] += 1
    elif fault == "union":
        report["union_scored_pairs"] += 1
    elif fault == "difference":
        report["differences"]["f1"] += .01
    elif fault in ("sparse_table", "reference_table", "transition_table"):
        ref = method["table"] if fault == "sparse_table" else report[fault]
        path = Path(ref["path"])
        text = path.read_text()
        path.write_text(text + text.splitlines()[-1] + "\n")
        new_ref = exporter.record(path)
        if fault == "sparse_table":
            method["table"] = new_ref
        else:
            report[fault] = new_ref
    elif fault == "input_identity":
        with Path(method["database"]["path"]).open("ab") as stream:
            stream.write(b"changed")
    ref = dump(inputs["output"] / "report.json", report)
    with pytest.raises(ValueError):
        reader.review(ref["path"], "0" * 64 if fault == "report_digest" else ref["sha256"])


@pytest.mark.parametrize("counts", (dict(TP=True, FP=1, FN=1), dict(TP=-1, FP=1, FN=1),
    dict(TP=1., FP=1, FN=1), dict(TP=0, FP=0, FN=1), dict(TP=1, FP=1)))
def test_invalid_counts(counts):
    with pytest.raises(ValueError):
        exporter.ratios(counts)


def test_raw_duplicate_and_cross_block_truth_refuse(tmp_path):
    path = tmp_path / "raw.gz"
    compressed(path, "a\tb\tTP\tA\tA\ts1\ts2\na\tb\tTP\tA\tA\ts1\ts2\n")
    with pytest.raises(ValueError, match="Duplicate"):
        exporter.aggregate(path, dict(A="A", B="B"))
    compressed(path, "a\tb\tTP\tA\tB\ts1\ts2\n")
    with pytest.raises(ValueError, match="Truth crosses"):
        exporter.aggregate(path, dict(A="A", B="B"))


def test_independent_components_keep_transitive_shared_labels(tmp_path):
    reference = compressed(tmp_path / "ref.gz", "1\t2\tZ\n2\t3\tY\n3\t4\tX\n5\t6\tW\n")
    truth, _ = exporter.reference_data(Path(reference["path"]))
    mapping, summary = exporter.reference_blocks(truth)
    _, _, independent, independent_summary = reader.reconstruct(reference["path"])
    assert mapping == independent == dict(W="W", X="X", Y="X", Z="X")
    assert summary == independent_summary and summary["shared_proteins"] == 2


def test_readback_has_no_exporter_or_helper_imports():
    import ast

    tree = ast.parse(Path(reader.__file__).read_text())
    imports = [n for n in ast.walk(tree) if isinstance(n, (ast.Import, ast.ImportFrom))]
    names = [n.module if isinstance(n, ast.ImportFrom) else a.name for n in imports
        for a in ([None] if isinstance(n, ast.ImportFrom) else n.names)]
    assert all(not name.startswith(("benchmark_tools", "igraph", "numpy", "scipy")) for name in names)


@pytest.mark.parametrize("text", (
    "a\tb\tTP\tA\tA\ts1\ts2\na\tb\tTP\tA\tA\ts1\ts2\n",
    "a\tb\tTP\tA\tA\ts1\ts2\na\tb\tFN\tA\tA\ts1\ts2\n",
    "a\tb\tFP\tA\tA\ts1\ts2\n",
    "a\tb\tTP\twrong\tA\ts1\ts2\n",
    "a\tc\tTP\tA\tB\ts1\ts2\n",
    "a\ta\tTP\tA\tA\ts1\ts1\n",
    "a\tb\tTN\tA\tA\ts1\ts2\n",
    "a\tb\tTP\n",
    "",
))
def test_independent_raw_refusals(tmp_path, text):
    ref = compressed(tmp_path / "raw.gz", text)
    with pytest.raises(ValueError):
        reader.raw_counts(ref["path"], {("a", "b")},
            dict(a=("A", "s1"), b=("A", "s2"), c=("B", "s2")), dict(A="A", B="B"))


@pytest.mark.parametrize("fault", ("missing", "nonbijective", "species"))
def test_independent_database_refusals(tmp_path, fault):
    path = tmp_path / "db.sqlite"
    rows = [(1, "a", "s1"), (2, "b", "s2")]
    if fault == "missing":
        rows.pop()
    elif fault == "nonbijective":
        rows[1] = (2, "a", "s2")
    else:
        rows.append((1, "other", "s2"))
    with sqlite3.connect(path) as connection:
        connection.execute("CREATE TABLE proteomes (prot_nr INT, uniprot_id TEXT, species TEXT)")
        connection.executemany("INSERT INTO proteomes VALUES (?, ?, ?)", rows)
    with pytest.raises(ValueError):
        reader.mapped_database(path, {1: "A", 2: "A"})


def test_database_mutation_during_primary_mapping_is_detected(inputs, monkeypatch):
    original = exporter.mapped_reference

    def change_after_read(database, *args):
        value = original(database, *args)
        with Path(database).open("ab") as stream:
            stream.write(b"changed-after-read")
        return value

    monkeypatch.setattr(exporter, "mapped_reference", change_after_read)
    with pytest.raises(ValueError, match="identity changed"):
        run(inputs)
