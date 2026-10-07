"""Invented grouping/event handoffs plus unchanged topology-kernel coverage."""

from copy import deepcopy
import csv
import hashlib
from pathlib import Path
import sys

import pytest

from benchmark_tools import trace_allocated_native_qfo_profile as trace
from benchmark_tools.prepare_ob_candidate_neighborhood import record
from tests.unit.test_export_native_qfo_factorial_scores import write
from tests.unit.test_trace_native_qfo_swiss_reconciliation import tree


def selection():
    return {name: dict(gene=name, source_family="Family0000000", root_hog="RootHOG0000000")
        for name in "abcd"}


def test_observed_positive_and_negative_pair_rules(tree):
    rows = selection()
    groups = {"Family0000000": set("abcd")}
    nodes = {"Family0000000": tree}
    positive = trace.describe("a", "b", rows, groups, nodes, {("a", "b")})
    negative = trace.describe("a", "c", rows, groups, nodes, set())
    assert positive["predicted"] and positive["lca"]["pair_event"] == "speciation"
    assert not negative["predicted"] and negative["lca"]["pair_event"] == "duplication"
    assert positive["candidate_members_a_sha256"] == negative["candidate_members_a_sha256"]


def test_different_candidate_exclusion_does_not_claim_lca():
    rows = selection()
    rows["b"]["source_family"] = "Family0000001"
    point = trace.describe("a", "b", rows, {"Family0000000": {"a"}, "Family0000001": {"b"}}, {}, set())
    assert not point["same_source_family"] and point["lca"] is None and not point["predicted"]


@pytest.mark.parametrize("a,b,predictions", [("a", "b", set()), ("a", "c", {("a", "c")})])
def test_prediction_annotation_disagreement_rejected(tree, a, b, predictions):
    with pytest.raises(ValueError, match="Native pair disagrees"):
        trace.describe(a, b, selection(), {"Family0000000": set("abcd")}, {"Family0000000": tree}, predictions)


def test_membership_hash_ignores_order_not_membership():
    assert trace.member_hash(["a", "b"]) == trace.member_hash(["b", "a"])
    assert trace.member_hash(["a", "b"]) != trace.member_hash(["a", "c"])


def fixture(tmp_path, monkeypatch, tree):
    cells = []
    for index in (7, 10):
        directory = tmp_path / str(index)
        directory.mkdir()
        roots = directory / "orthohmm_root_hogs.tsv"
        rows = [("a,b", 0), ("c,d", 1)] if index == 7 else [("a,d", 0), ("b,c", 1)]
        roots.write_text("root_hog\tsource_family\tgenes\n" + "".join(
            f"RootHOG{i:07d}\tFamily0000000\t{genes}\n" for genes, i in rows))
        manifest_ref = write(directory / "provenance_manifest.json", dict(membership_reconciliation=None,
            pair_orthology_rule="positive_paralogy", input_cluster_sha256=hashlib.sha256(b"a b c d\n").hexdigest(),
            species_tree_sha256="invented_species_tree"))
        nodes = deepcopy(tree)
        if index == 10:
            nodes["b"]["parent_node_id"], nodes["d"]["parent_node_id"] = "G1", "G0"
            nodes["G0"]["genes"], nodes["G1"]["genes"] = set("ad"), set("bc")
        node_path = directory / "orthohmm_reconciliation_nodes.tsv"
        with node_path.open("w", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=trace.topology.NODE_HEADER, delimiter="\t", lineterminator="\n")
            writer.writeheader()
            for node in nodes.values():
                writer.writerow(dict(node, genes=",".join(sorted(node["genes"])), species=",".join(sorted(node["species"]))))
        pair_path = directory / "pairs.tsv"
        pair_path.write_text("gene_a\tspecies_a\tgene_b\tspecies_b\n" + (
            "a\tA\tb\tB\nc\tA\td\tB\n" if index == 7 else "a\tA\td\tB\nb\tB\tc\tA\n"))
        cell = trace.CELLS[0 if index == 7 else 1]
        reviewed_ref = write(directory / "outputs.json", dict(native_outputs_validated=True, cell=cell,
            input_genes=4, checked_files=[record(roots), manifest_ref, record(pair_path)],
            phylogeny=dict(root_hogs=len(rows), native_pair_rows=2), source=record(roots)))
        review = dict(source=record(roots))
        review.update(outputs=reviewed_ref) if index == 7 else review.update(reviews=dict(outputs_or_failure=reviewed_ref))
        review_ref = write(directory / "review.json", review)
        stage = dict(native_input=record(pair_path), conversion_kind="native")
        stage["scientific_recovery" if index == 7 else "terminal_review"] = review_ref
        admission = dict(accuracy_admitted=True, native_index=index, cell=cell, native_job_id=index,
            conversion=stage, resources=None, eligible_for_timing_comparison=False)
        admission_ref = write(directory / "admission.json", admission)
        audit_ref = write(directory / "counts.json", dict(cells=[dict(cell=cell, index=index,
            native_job_id=index, admission=admission_ref, raw_file=record(roots))]))
        cells.append((cell, audit_ref))
    binding_ref = write(tmp_path / "binding.json", dict(schema="allocated_native_qfo_retained_swiss_uncertainty_binding_v1",
        bound_cells={cell: dict(count_audit=ref, status="native_records_matched") for cell, ref in cells},
        families=[f"family{i}" for i in range(18)]))
    readback = dict(schema="allocated_native_qfo_profile_swiss_rational_readback_v1",
        source=record(Path(trace.__file__).with_name("readback_allocated_native_qfo_profile_swiss.py")),
        cells=list(trace.CELLS), profile_pair_labels_matched=10765, families_checked=18,
        prior_matched_contrasts_checked=2, prior_matched_contrasts_unchanged=True,
        new_accuracy_or_resource_admission=False, scientific_timings_admitted=False,
        independent_confirmation=False, publication_ready=False, checked_inputs=[], binding=binding_ref)
    # These small interface fixtures do not pretend to recreate a production
    # raw reference; those primitives and the readback have separate tests.
    monkeypatch.setattr(trace.transitions, "read_labels", lambda *args: {})
    monkeypatch.setattr(trace.transitions, "compare_labels", lambda *args:
        (dict(changed_relations=1), [("family0", "a", "b", "TP", "FN")]))
    return write(tmp_path / "readback.json", readback), tmp_path / "trace.json"


def test_complete_streaming_handoff_identifies_event_change_with_same_members(tmp_path, monkeypatch, tree):
    ref, output = fixture(tmp_path, monkeypatch, tree)
    result = trace.run(ref["path"], ref["sha256"], output)
    assert result["changed_pairs_traced"] == 1 and result["summary"] == {"positive_paralogy_exclusion": 1}
    row = result["changed_pairs"][0]
    assert row["before"]["lca"]["pair_event"] == "speciation"
    assert row["after"]["lca"]["pair_event"] == "duplication"
    assert row["candidate_members_a_unchanged"] and row["candidate_members_b_unchanged"]
    assert result["species_tree_bytes_identical"] and not result["node_annotations_previously_inventoried"]
    assert not result["new_scoring_or_admission"] and not result["publication_ready"]


@pytest.mark.parametrize("key,value", [("schema", "cached"), ("source", {}), ("cells", []),
    ("profile_pair_labels_matched", 1), ("families_checked", 1), ("prior_matched_contrasts_checked", 0),
    ("prior_matched_contrasts_unchanged", False), ("publication_ready", True), ("scientific_timings_admitted", True)])
def test_readback_scope_cannot_authorize_trace(tmp_path, monkeypatch, tree, key, value):
    import json
    ref, output = fixture(tmp_path, monkeypatch, tree)
    value_dict = json.loads(Path(ref["path"]).read_text())
    value_dict[key] = value
    changed = write(Path(ref["path"]), value_dict)
    with pytest.raises(ValueError, match="read-back profile contrast"):
        trace.run(changed["path"], changed["sha256"], output)


def test_existing_output_not_overwritten(tmp_path, monkeypatch):
    output = tmp_path / "existing"
    output.write_text("retain me")
    monkeypatch.setattr(sys, "argv", ["trace", "--readback", "absent", "--readback-sha256", "unused",
        "--output", str(output)])
    with pytest.raises(ValueError, match="fresh"):
        trace.main()
    assert output.read_text() == "retain me"
