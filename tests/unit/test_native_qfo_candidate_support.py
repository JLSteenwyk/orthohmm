"""Support associations are descriptive event units, never pair-IID confidence."""

import copy
import csv
import json
import math
from pathlib import Path
import subprocess
import sys

import pytest

from benchmark_tools import export_native_qfo_candidate_support as primary
from benchmark_tools import readback_native_qfo_candidate_support as independent


def event(source, target, iteration=0):
    result = dict(source_genes=source, target_genes=target, source_cluster=1, target_cluster=2, iteration=iteration)
    result.update({name: 1 for name in primary.FEATURES})
    result.update(source_size=len(source), target_size=len(target), margin=math.inf, support=0.1)
    return result


def row(a, b, state, index, iteration=0, path="direct_cross_endpoint"):
    return [a, b, a, b, "FN" if state == "TP" else "not_scored", state,
            "old_" + a, "old_" + b, "new", str(iteration), str(index) if index is not None else "", path]


@pytest.fixture
def specimen():
    trace = [event(["a", "b"], ["c", "d"]), event(["e"], ["f"], 1),
             event(["g"], ["h"]), event(["i"], ["j"])]
    trace[1]["margin"] = 2.0
    trace[2]["support"] = 0.2
    trace[3]["support"] = 0.3
    rows = [row("a", "c", "TP", 0), row("b", "c", "FP", 0), row("a", "d", "FP", 0),
            row("e", "f", "TP", 1, 1), row("g", "h", "FP", 2),
            row("e", "h", "FP", None, 1, "transitive_union")]
    return trace, rows


def test_independent_algorithms_agree_with_mixed_unlabeled_and_transitive_cases(specimen):
    annotated, events, summary = primary.analyze(*specimen)
    other_pairs, other_events, other_summary = independent.reconstruct(*specimen)
    assert annotated == other_pairs and events == other_events
    independent.equivalent(summary, other_summary)
    assert summary["totals"] == dict(accepted_events=4, changed_pairs=6, implicated_events=3,
        direct_tp_pairs=2, direct_fp_pairs=3, transitive_tp_pairs=0, transitive_fp_pairs=1)
    assert [r["events"] for r in summary["cohort_summaries"][:4]] == [1, 1, 1, 1]
    assert annotated[-1] == specimen[1][-1] + ["transitive_no_direct_event"] + [""] * len(primary.FEATURES)
    assert events[0][4:7] == ["mixed", "1", "2"]


def test_each_event_once_not_weighted_by_dependent_pairs(specimen):
    summary = primary.analyze(*specimen)[2]
    assert summary["cohort_summaries"][2]["features"]["support"]["finite"] == 1
    assert summary["cohort_summaries"][2]["features"]["support"]["mean"] == 0.1
    assert summary["cohort_summaries"][2]["features"]["margin"]["positive_infinity"] == 1
    assert summary["cohort_summaries"][2]["features"]["margin"]["mean"] is None
    empty = summary["cohort_summaries"][-1]
    assert empty["events"] == 0
    assert all(s["missing"] == 0 and s["minimum"] is None and s["mean"] is None for s in empty["features"].values())


@pytest.mark.parametrize("field,value", [
    ("support", -1), ("support", math.inf), ("support", math.nan), ("support", True), ("support", "1"),
    ("margin", -math.inf), ("margin", math.nan), ("margin", None), ("forward_hits", 1.5),
    ("source_seed_families", False), ("iteration", True), ("iteration", 2),
    ("source_size", 3), ("source_size", 2.0), ("source_cluster", -1), ("source_cluster", True),
    ("source_genes", []), ("source_genes", ["a", "a"]), ("source_genes", ["a", " c"]),
    ("source_genes", ["a", "c"]), ("source_genes", ["a", "b\tc"]),
])
def test_both_reject_invalid_event_features_and_memberships(specimen, field, value):
    trace, rows = specimen
    trace[0][field] = value
    with pytest.raises(ValueError):
        primary.analyze(trace, rows)
    with pytest.raises(ValueError):
        independent.reconstruct(trace, rows)


def test_both_reject_duplicate_semantic_events_even_if_order_differs(specimen):
    trace, rows = specimen
    repeated = copy.deepcopy(trace[0])
    repeated["source_genes"].reverse()
    repeated["target_genes"].reverse()
    trace.append(repeated)
    with pytest.raises(ValueError, match="Duplicate"):
        primary.analyze(trace, rows)
    with pytest.raises(ValueError, match="Repeated"):
        independent.reconstruct(trace, rows)


@pytest.mark.parametrize("column,value", [(2, "z"), (3, "a"), (4, "FP"), (5, "TN"), (7, "old_a"),
    (8, ""), (9, "1"), (9, "2"), (10, "9"), (10, "00"), (10, "-1"), (10, ""), (11, "invented")])
def test_both_reject_invalid_direct_pair_references(specimen, column, value):
    trace, rows = specimen
    rows[0][column] = value
    with pytest.raises(ValueError):
        primary.analyze(trace, rows)
    with pytest.raises(ValueError):
        independent.reconstruct(trace, rows)


def test_transitive_cannot_be_assigned_direct_event(specimen):
    trace, rows = specimen
    rows[-1][10] = "1"
    with pytest.raises(ValueError, match="Transitive"):
        primary.analyze(trace, rows)
    with pytest.raises(ValueError, match="Transitive"):
        independent.reconstruct(trace, rows)


@pytest.mark.parametrize("collapse", [False, True])
def test_both_reject_duplicate_or_collapsed_native_pairs(specimen, collapse):
    trace, rows = specimen
    new = rows[0].copy()
    if collapse:
        new[:2] = ["alias_a", "alias_c"]
    rows.append(new)
    with pytest.raises(ValueError):
        primary.analyze(trace, rows)
    with pytest.raises(ValueError):
        independent.reconstruct(trace, rows)


def test_finite_even_medians_and_rational_means():
    events = [event([str(i)], ["t" + str(i)]) for i in range(4)]
    for item, value in zip(events, (0.1, 0.2, 0.3, 0.4)):
        item["support"] = value
    reported = primary.summaries(events)
    expected = independent.rational_features(events)
    independent.equivalent(reported, expected)
    assert expected["support"]["median"] == expected["support"]["mean"] == 0.25


def test_rational_comparison_requires_exact_integer_counts_and_identities():
    independent.equivalent({"count": 1, "mean": 0.30000000000000004}, {"count": 1, "mean": 0.3})
    for actual in ({"count": True, "mean": 0.3}, {"count": 1.0, "mean": 0.3}, {"count": 1, "mean": 0.31}):
        with pytest.raises(ValueError):
            independent.equivalent(actual, {"count": 1, "mean": 0.3})


@pytest.mark.parametrize("module", [primary, independent])
def test_json_duplicate_keys_rejected(tmp_path, module):
    path = tmp_path / "ambiguous.json"
    path.write_text('{"a":1,"a":2}')
    loader = primary.read_json if module is primary else independent.json_data
    with pytest.raises(ValueError, match="Duplicate"):
        loader(path)


def evidence_fixture(tmp_path, specimen):
    trace_path, ledger_path = tmp_path / "trace.json", tmp_path / "ledger.tsv"
    trace_path.write_text(json.dumps(specimen[0]))
    primary.write_table(ledger_path, primary.HEADER, specimen[1])
    trace_ref, ledger_ref = primary.identity(trace_path), primary.identity(ledger_path)
    shared = dict(genes=10, baseline_groups=10, candidate_groups=6, accepted_merges=4, changed_pairs=6,
        bridges=[], localized_summary=primary.analyze(*specimen)[2]["localized_summary"],
        whole_candidate_partition_reconstructed=True, complete_pair_localization=True,
        original_scored_identifiers_preserved=True, failed_r1_timing_remains_ineligible=True,
        accuracy_rescored=False, native_inference_reexecuted=False, new_scoring_or_admission=False,
        uncertainty_admitted=False, publication_ready=False)
    prior = dict(shared, schema="native_qfo_candidate_alias_group_trace_v1",
        status="original_protein_identity_join_and_all_changed_pairs_localized",
        native_trace=trace_ref, ledger=ledger_ref, checked_records=[trace_ref])
    report_path = tmp_path / "prior.json"
    report_path.write_text(json.dumps(prior))
    readback = dict(shared, schema="native_qfo_candidate_alias_group_readback_v1",
        status="original_protein_identity_join_and_all_paths_independently_verified", report=primary.identity(report_path))
    readback_path = tmp_path / "readback.json"
    readback_path.write_text(json.dumps(readback))
    kwargs = dict(report_sha=primary.identity(report_path)["sha256"], readback_sha=primary.identity(readback_path)["sha256"],
                  trace_sha=trace_ref["sha256"], ledger_sha=ledger_ref["sha256"])
    return report_path, readback_path, kwargs


def test_bound_inputs_and_current_byte_checks(tmp_path, specimen):
    a, b, kwargs = evidence_fixture(tmp_path, specimen)
    prior, refs, trace, rows = primary.inputs(a, b, **kwargs)
    assert trace == specimen[0] and rows == specimen[1]
    (tmp_path / "ledger.tsv").write_text("changed")
    with pytest.raises(ValueError, match="Changed current"):
        primary.inputs(a, b, **kwargs)
    with pytest.raises(ValueError, match="changed"):
        independent.unchanged(refs[3])


def test_selected_defaults_reject_unanchored_fixture(tmp_path, specimen):
    a, b, kwargs = evidence_fixture(tmp_path, specimen)
    with pytest.raises(ValueError, match="anchor"):
        primary.inputs(a, b)


def test_export_preserves_full_tables_and_refuses_overwrite(tmp_path, specimen, monkeypatch):
    a, b, kwargs = evidence_fixture(tmp_path, specimen)
    loader = primary.inputs
    monkeypatch.setattr(primary, "inputs", lambda a, b: loader(a, b, **kwargs))
    destination = tmp_path / "result"
    result = primary.export(a, b, destination)
    assert result["totals"]["implicated_events"] == 3
    assert independent.table_data(result["annotated_pairs"]["path"], primary.PAIR_HEADER) == primary.analyze(*specimen)[0]
    for flag in ("grouping_replayed", "aliases_rederived", "accuracy_rescored", "defaults_changed", "uncertainty_admitted", "publication_ready"):
        assert result[flag] is False
    original = (destination / "report.json").read_bytes()
    with pytest.raises(FileExistsError):
        primary.export(a, b, destination)
    assert (destination / "report.json").read_bytes() == original


def test_failure_retained_and_never_retried(tmp_path, specimen):
    a, b, _ = evidence_fixture(tmp_path, specimen)
    destination = tmp_path / "failed"
    with pytest.raises(ValueError):
        primary.export(a, b, destination)
    assert json.loads((destination / "failure.json").read_text())["automatic_retry"] is False
    with pytest.raises(FileExistsError):
        primary.export(a, b, destination)
    assert not (destination / "report.json").exists()


def test_partial_tables_rejected(tmp_path):
    path = tmp_path / "partial.tsv"
    primary.write_table(path, primary.PAIR_HEADER, [["a"]])
    with pytest.raises(ValueError, match="Partial"):
        independent.table_data(path, primary.PAIR_HEADER)


def test_independent_cli_works_without_site_or_project_imports():
    reader = Path(independent.__file__)
    run = subprocess.run([sys.executable, "-I", "-S", "-B", str(reader), "--help"], capture_output=True, text=True)
    assert run.returncode == 0 and "--report-sha256" in run.stdout
    assert "benchmark_tools" not in reader.read_text().split("if __name__")[0].replace(
        '"benchmark_tools/results/NATIVE_QFO_CANDIDATE_SUPPORT_PROTOCOL_20261006.md"', "")
