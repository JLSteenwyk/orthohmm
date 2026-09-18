import csv
import json

import pytest

from benchmark_tools.audit_reconciled_pair_events import audit, reconstruct, record


def rows():
    return [
        dict(node_id="a", parent_node_id="r", event="leaf", pair_event="leaf", genes="a", species="A"),
        dict(node_id="b", parent_node_id="r", event="leaf", pair_event="leaf", genes="b", species="B"),
        dict(node_id="r", parent_node_id="", event="speciation", pair_event="speciation", genes="a,b",
             species="A,B", mapping_conflict="false", species_overlap_count="0")]


def test_speciation_and_membership_filter():
    data = rows()
    groups = {"a": "g1", "b": "g2"}
    assert reconstruct(data, groups, False)[0] == {("a", "b")}
    pairs, summary = reconstruct(data, groups, True)
    assert pairs == set() and summary["group_filtered_pairs"] == 1


def test_mapping_conflict_is_uncertain_not_pair_duplication():
    data = rows()
    data[-1].update(event="duplication", pair_event="uncertain", mapping_conflict="true")
    assert reconstruct(data, {"a": "g", "b": "g"}, False)[0] == {("a", "b")}


def test_species_overlap_blocks_cross_child_pairs():
    data = rows()
    data[1]["species"] = "A"
    data[-1].update(event="duplication", pair_event="duplication", species="A", species_overlap_count="1")
    assert not reconstruct(data, {"a": "g", "b": "g"}, False)[0]


def test_duplication_blocks_cross_species_pairs_between_subtrees():
    data = rows()
    data[0]["parent_node_id"] = data[1]["parent_node_id"] = "left"
    data[-1].update(node_id="left", parent_node_id="root")
    data.extend([
        dict(node_id="c", parent_node_id="right", event="leaf", pair_event="leaf", genes="c", species="A"),
        dict(node_id="d", parent_node_id="right", event="leaf", pair_event="leaf", genes="d", species="C"),
        dict(node_id="right", parent_node_id="root", event="speciation", pair_event="speciation", genes="c,d",
             species="A,C", mapping_conflict="false", species_overlap_count="0"),
        dict(node_id="root", parent_node_id="", event="duplication", pair_event="duplication", genes="a,b,c,d",
             species="A,B,C", mapping_conflict="false", species_overlap_count="1")])
    assert reconstruct(data, dict.fromkeys("abcd", "g"), False)[0] == {("a", "b"), ("c", "d")}


@pytest.mark.parametrize("field,value", [("genes", "a,a"), ("genes", "a"), ("species", "A"),
    ("event", "duplication"), ("pair_event", "uncertain"), ("mapping_conflict", "unknown"),
    ("species_overlap_count", "1"), ("parent_node_id", "missing")])
def test_invalid_node_fields(field, value):
    data = rows()
    data[-1][field] = value
    with pytest.raises(ValueError):
        reconstruct(data, {"a": "g", "b": "g"}, False)


def test_incomplete_groups():
    with pytest.raises(ValueError):
        reconstruct(rows(), {"a": "g"}, False)


def test_duplicate_node():
    data = rows()
    with pytest.raises(ValueError):
        reconstruct(data + [data[0]], {"a": "g", "b": "g"}, False)


def test_disconnected_cycle():
    data = rows()
    data.extend([dict(node_id="x", parent_node_id="y"), dict(node_id="y", parent_node_id="x")])
    with pytest.raises(ValueError, match="Disconnected"):
        reconstruct(data, {"a": "g", "b": "g"}, False)


def write_table(path, fields, entries):
    with path.open("w") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(entries)


@pytest.mark.parametrize("case", ["unexpanded", "expanded", "missing_pair", "table_changed", "wrong_admission", "duplicate_pair"])
def test_file_audit(tmp_path, case):
    cell = "p0_c1_r1" if case == "expanded" else "p0_c0_r1"
    nodes = tmp_path / "orthohmm_reconciliation_nodes.tsv"
    groups = tmp_path / "orthohmm_root_hogs.tsv"
    pairs = tmp_path / "orthohmm_pairwise_orthologs.tsv"
    data = rows()
    for row in data:
        row["source_family"] = "Family1"
    fields = sorted(set().union(*(set(row) for row in data)))
    write_table(nodes, fields, data)
    group_rows = ([dict(root_hog="g", source_family="Family1", genes="a,b")]
                  if case != "expanded" else
                  [dict(root_hog=g, source_family="Family1", genes=g) for g in "ab"])
    write_table(groups, ["root_hog", "source_family", "genes"], group_rows)
    pair = dict(gene_a="a", species_a="A", gene_b="b", species_b="B")
    pair_rows = [] if case in {"expanded", "missing_pair"} else [pair]
    if case == "duplicate_pair":
        pair_rows.append(pair)
    write_table(pairs, list(pair), pair_rows)
    artifacts = [{"absolute_path": str(path), **{k: record(path)[k] for k in ("bytes", "sha256")}}
                 for path in (nodes, groups)]
    artifacts.append({"absolute_path": str(tmp_path / "gene_trees/Family1.reconciled.nwk")})
    execution = tmp_path / "execution.json"
    execution.write_text(json.dumps({"methods": {cell: {"outputs": artifacts}}}))
    native_manifest = tmp_path / "manifest.json"
    native_manifest.write_text("{}")
    admission = tmp_path / "admission.json"
    admission.write_text(json.dumps({"status": "qfo_native_pair_output_verified", "cell": cell,
        "native_pairs": record(pairs), "native_group_integrity": {
            "native_manifest": record(native_manifest), "integrity": {"execution_status": record(execution)}}}))
    expected_sha = record(admission)["sha256"]
    if case == "table_changed":
        with groups.open("a") as stream:
            stream.write("\n")
    if case == "wrong_admission":
        expected_sha = "wrong"
    if case in {"unexpanded", "expanded"}:
        result = audit(admission, expected_sha, 1)
        assert result["selected_families"] == ["Family1"]
        assert result["families"]["Family1"]["retained_pairs"] == (case == "unexpanded")
        assert result["accuracy_evaluated"] is False
    else:
        with pytest.raises(ValueError):
            audit(admission, expected_sha, 1)
