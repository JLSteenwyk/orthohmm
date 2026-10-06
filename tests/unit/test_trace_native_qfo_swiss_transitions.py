"""Pair identity and truth safeguards for the native SwissTrees diagnostic."""

from collections import Counter
import copy
import csv
import gzip
import json
from pathlib import Path

import numpy as np
import pytest

from benchmark_tools import trace_native_qfo_swiss_transitions as trace
from benchmark_tools.prepare_ob_candidate_neighborhood import record


@pytest.fixture
def example():
    families = [f"family{i:02d}" for i in range(18)]
    before, after = {}, {}
    for i, f in enumerate(families):
        for j, label in enumerate(("TP", "FP", "FN")):
            key = (f, f + str(j * 2), f + str(j * 2 + 1))
            before[key] = label
            after[key] = ({"TP": "FN", "FP": "TN", "FN": "FN"}[label] if i < 9 else label)
    return families, before, after


def audited(labels, families):
    rows, values = [], []
    for f in families:
        counts = Counter({label: 0 for label in trace.LABELS})
        genes = set()
        for (family, a, b), label in labels.items():
            if family == f:
                counts[label] += 1
                genes.update((a, b))
        scores = trace.statistics(counts)
        rows.append(dict(family=f, counts_without_prior=dict(counts), represented_genes=sorted(genes),
                         statistics_with_prior=scores))
        values.append([scores["PPV"], scores["TPR"]])
    return dict(families=rows, aggregate=dict(zip(trace.METRICS,
        trace.aggregate(np.asarray(values).mean(axis=0)).tolist())))


def write_raw(path, labels, header=trace.HEADER):
    with gzip.open(path, "wt", newline="") as handle:
        handle.write(header + "\n")
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerows((*key, label) for key, label in labels.items())


def test_exact_transitions_and_marginals(example):
    families, before, after = example
    result, changes = trace.compare_labels(before, after, families, audited(before, families), audited(after, families))
    assert result["reference_relations"] == 54 and result["changed_relations"] == len(changes) == 18
    assert result["removed_true_positives"] == result["removed_false_positives"] == 9
    assert result["added_true_positives"] == result["added_false_positives"] == 0
    assert result["after_predictions_subset_on_reference"] is True
    assert len(result["transitions"]) == 16 and sum(result["transitions"].values()) == 54
    assert changes == sorted(changes)
    assert [row["family"] for row in result["families"]] == families


def test_added_predictions_are_not_hidden(example):
    families, before, after = example
    before, after = after, before
    result, _ = trace.compare_labels(before, after, families, audited(before, families), audited(after, families))
    assert result["added_true_positives"] == result["added_false_positives"] == 9
    assert result["after_predictions_subset_on_reference"] is False


def test_equal_counts_do_not_prove_equal_decisions(example):
    families, before, _ = example
    after = dict(before)
    keys = list(after)
    after[keys[0]], after[keys[2]] = after[keys[2]], after[keys[0]]
    assert audited(before, families) == audited(after, families)
    result, changes = trace.compare_labels(before, after, families, audited(before, families), audited(after, families))
    assert result["changed_relations"] == len(changes) == 2
    assert result["removed_true_positives"] == result["added_true_positives"] == 1


@pytest.mark.parametrize("label", ("FP", "TN"))
def test_changed_reference_truth_refuses(example, label):
    families, before, after = example
    after[next(iter(after))] = label
    with pytest.raises(ValueError, match="reference truth"):
        trace.compare_labels(before, after, families, audited(before, families), audited(after, families))


def test_missing_relation_refuses(example):
    families, before, after = example
    after.pop(next(iter(after)))
    with pytest.raises(ValueError, match="reference relation universe"):
        trace.compare_labels(before, after, families, audited(before, families), audited(after, families))


@pytest.mark.parametrize("field", ("counts_without_prior", "represented_genes", "statistics_with_prior"))
def test_changed_family_audit_refuses(example, field):
    families, before, _ = example
    report = audited(before, families)
    report["families"][0][field] = {}
    with pytest.raises(ValueError):
        trace.verified_counts(before, families, report)


def test_changed_macro_refuses(example):
    families, before, _ = example
    report = audited(before, families)
    report["aggregate"]["F1"] += .01
    with pytest.raises(ValueError, match="macro statistic"):
        trace.verified_counts(before, families, report)


def test_reader_exact_and_orientation_invariant(tmp_path, example):
    families, before, _ = example
    path = tmp_path / "raw.gz"
    write_raw(path, {(f, b, a): label for (f, a, b), label in before.items()})
    assert trace.read_labels(path, families) == before


@pytest.mark.parametrize("fault", ("header", "duplicate", "self", "label", "family", "empty_id", "width", "missing_family"))
def test_bad_raw_refuses(tmp_path, example, fault):
    families, before, _ = example
    path = tmp_path / "raw.gz"
    write_raw(path, before, "wrong" if fault == "header" else trace.HEADER)
    first = next(iter(before))
    with gzip.open(path, "at") as handle:
        row = (*first, "TP")
        if fault == "duplicate":
            row = (first[0], first[2], first[1], "TP")
        elif fault == "self":
            row = (first[0], first[1], first[1], "TP")
        elif fault == "label":
            row = (*first, "XX")
        elif fault == "family":
            row = ("unknown", "a", "b", "TP")
        elif fault == "empty_id":
            row = (first[0], "", "b", "TP")
        elif fault == "width":
            row = ("a", "b")
        handle.write("\t".join(row) + "\n")
    if fault == "missing_family":
        write_raw(path, {key: label for key, label in before.items() if key[0] != families[0]})
    with pytest.raises(ValueError):
        trace.read_labels(path, families)


@pytest.mark.parametrize("families", (["x"] * 18, ["x"], []))
def test_bad_inventory_refuses(tmp_path, families):
    with pytest.raises(ValueError, match="18 distinct"):
        trace.read_labels(tmp_path / "absent.gz", families)


@pytest.mark.parametrize("fault", (None, "audit_source", "audit_scope", "audit_identity", "raw_not_admitted",
                                 "count_mismatch", "truth_changed", "timing_relabel", "duplicate_cell",
                                 "reference_changed", "snapshot_source"))
def test_integration_with_explicit_snapshot_replay_stub(tmp_path, monkeypatch, example, fault):
    families, before, after = example
    def dump(name, value):
        path = tmp_path / name
        path.write_text(json.dumps(value))
        return record(path)

    rows, audit_refs = [], []
    reference = dump("reference.json", {"reference": True})
    for index, (cell, labels, (schema, status, source)) in enumerate(zip(trace.CELLS, (before, after), trace.AUDITS)):
        if fault == "truth_changed" and index == 1:
            labels = dict(labels)
            labels[next(iter(labels))] = "FP"
        raw_path = tmp_path / str(index) / "SwissTrees" / "raw.txt.gz"
        raw_path.parent.mkdir(parents=True)
        write_raw(raw_path, labels)
        raw = record(raw_path)
        execution = dump(f"execution{index}.json", dict(outputs=[raw]))
        admission = dump(f"admission{index}.json", dict(execution_report=execution,
            checked_records=[execution] if fault == "raw_not_admitted" and index == 1 else [raw, execution]))
        count = dict(audited(labels, families), cell=cell, index=index + 6, native_job_id=index + 22435,
                     admission=admission, raw_file=raw, native_endpoint_f1=.7)
        audit = dict(schema=schema, status=status, source=record(source.__file__),
            families=families, reference=reference, cells=[count], reference_relation_count=54,
            checked_inputs=[raw], helpers=[], new_accuracy_or_resource_admission=False,
            historical_intervals_attached=False, independent_confirmation=False, publication_ready=False,
            new_bootstrap_draws=0)
        if index == 1:
            if fault == "audit_source":
                audit["source"] = reference
            elif fault == "audit_scope":
                audit["new_accuracy_or_resource_admission"] = True
            elif fault == "audit_identity":
                count["native_job_id"] += 1
            elif fault == "count_mismatch":
                count["families"][0]["counts_without_prior"]["TP"] += 1
            elif fault == "reference_changed":
                audit["reference"] = dump("other_reference.json", {"reference": False})
        audit_ref = dump(f"audit{index}.json", audit)
        audit_refs.append((audit_ref["path"], audit_ref["sha256"]))
        rows.append(dict(cell=cell, admission=admission, index=index + 6, native_job_id=index + 22435,
            accuracy_admitted=True, scores=dict(SwissTrees=.7),
            **(dict(timing_eligible=False, timing_admitted=False, resources=None) if index else {}),
            status="supplied_native_admission" if index == 0 else "supplied_recovered_scientific_admission"))
    snapshot = dict(schema="native_qfo_scientific_reporting_snapshot_v1", source=record(trace.reporter.__file__),
        rows=rows, publication_ready=False, plan=dict(path="stub", sha256="stub"))
    if fault == "timing_relabel":
        rows[1]["timing_eligible"] = True
    elif fault == "duplicate_cell":
        rows.append(copy.deepcopy(rows[0]))
    elif fault == "snapshot_source":
        snapshot["source"] = reference
    snapshot_ref = dump("snapshot.json", snapshot)
    monkeypatch.setattr(trace.reporter, "collect", lambda *args: copy.deepcopy(snapshot))
    output, changes = tmp_path / "result.json", tmp_path / "changes.tsv"
    if fault is not None:
        with pytest.raises(ValueError):
            trace.run(snapshot_ref["path"], snapshot_ref["sha256"], audit_refs, output, changes)
        assert not output.exists() and not changes.exists()
        return
    result = trace.run(snapshot_ref["path"], snapshot_ref["sha256"], audit_refs, output, changes)
    assert json.loads(output.read_text()) == result
    assert result["comparison"]["removed_true_positives"] == 9
    assert result["changed_relations_ledger"] == record(changes)
    with changes.open() as handle:
        assert len(list(csv.DictReader(handle, delimiter="\t"))) == 18
    for flag in ("new_scoring_or_admission", "uncertainty_admitted", "scientific_timings_admitted",
                 "independent_confirmation", "publication_ready"):
        assert result[flag] is False
    assert result["cells"][1]["timing_eligible"] is False
    assert not {"timing_eligible", "timing_admitted", "resources"}.intersection(result["cells"][0])
    with pytest.raises(ValueError, match="fresh output"):
        trace.run(snapshot_ref["path"], snapshot_ref["sha256"], audit_refs, output, changes)
