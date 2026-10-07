"""Invented saved-tree grouping transitions; no selected production execution."""

from copy import deepcopy
import hashlib
import json
from pathlib import Path
import sys

import pytest

from benchmark_tools import readback_allocated_native_qfo_profile_newick as reader
from tests.unit.test_export_native_qfo_factorial_scores import write


def trees(directory, family, genes, newick):
    checkpoint = dict(schema_version=2, status="complete", family_id=family,
                      species_tree_sha256=directory.name, genes=list(genes))
    for suffix, key in (("raw", "raw_tree_sha256"), ("rooted", "rooted_tree_sha256"),
                        ("reconciled", "annotated_tree_sha256")):
        path = directory / "gene_trees" / (family + "." + suffix + ".nwk")
        path.write_text(newick)
        checkpoint[key] = reader.record(path)["sha256"]
    write(directory / "checkpoints" / (family + ".json"), checkpoint)


def fixture(tmp_path, category="candidate_separation"):
    observed, evidence = [], []
    points = []
    for index, cell in enumerate(reader.CELLS):
        directory = tmp_path / cell
        (directory / "checkpoints").mkdir(parents=True)
        (directory / "gene_trees").mkdir()
        separated = index == 1 and category == "candidate_separation"
        groups = {"Family0000000": "ac", "Family0000001": "bd"} if separated else {"Family0000000": "abcd"}
        text = "root_hog\tsource_family\tgenes\n"
        for ordinal, (family, genes) in enumerate(groups.items()):
            text += f"RootHOG{ordinal:07d}\t{family}\t{','.join(genes)}\n"
            code = "D" if index == 1 and not separated else "S"
            newick = f"(a,b,c,d)'G00000|{code}@S0';\n" if not separated else f"({','.join(genes)})'G00000|S@S0';\n"
            trees(directory, family, genes, newick)
        root = directory / "orthohmm_root_hogs.tsv"
        root.write_text(text)
        manifest = write(directory / "provenance_manifest.json", dict(species_tree_sha256=cell,
            membership_reconciliation=None, pair_orthology_rule="positive_paralogy"))
        node = directory / "nodes.tsv"
        node.write_text("not consumed by the independent reader\n")
        root_ref = reader.record(root)
        evidence.extend((root_ref, manifest, reader.record(node)))
        observed.append(dict(cell=cell, root_hogs=root_ref, manifest=manifest,
            observed_node_annotations=reader.record(node), species_tree_sha256=cell,
            selected_source_families=len(groups)))
        point = dict(gene_a="a", gene_b="b", same_source_family=not separated,
            predicted=index == 0, lca=None if separated else dict(node_id="G00000", pair_event="S"))
        if point["lca"]:
            point["lca"]["pair_event"] = "duplication" if index == 1 else "speciation"
        for side, gene in (("a", "a"), ("b", "b")):
            family = next(f for f, genes in groups.items() if gene in genes)
            ordinal = list(groups).index(family)
            point["source_family_" + side] = family
            point["root_hog_" + side] = f"RootHOG{ordinal:07d}"
            point["candidate_size_" + side] = len(groups[family])
            point["candidate_members_" + side + "_sha256"] = hashlib.sha256(
                ("\n".join(sorted(groups[family])) + "\n").encode()).hexdigest()
        points.append(point)
    pair = dict(family="invented", protein_a="a", protein_b="b", before_label="TP", after_label="FN",
        before=points[0], after=points[1], localization=category)
    for side in ("a", "b"):
        pair["candidate_members_" + side + "_unchanged"] = points[0]["candidate_members_" + side + "_sha256"] == points[1]["candidate_members_" + side + "_sha256"]
    dummy = tmp_path / "earlier_readback.json"
    dummy.write_text("{}\n")
    report = dict(schema="allocated_native_qfo_profile_pair_localization_v1",
        source=reader.record(Path(reader.__file__).with_name("trace_allocated_native_qfo_profile.py")),
        readback=reader.record(dummy), helpers=[reader.record(Path(reader.__file__).with_name(name)) for name in (
            "trace_native_qfo_swiss_reconciliation.py", "trace_native_qfo_swiss_transitions.py")],
        evidence=evidence, states=observed, changed_pairs=[pair], changed_pairs_traced=1,
        comparison=dict(families=[dict(family="invented", transitions={"TP->FN": 1})]),
        summary={category: 1}, species_tree_bytes_identical=False, node_annotations_previously_inventoried=False)
    report.update({flag: False for flag in reader.SCOPE_FLAGS})
    return report


def run(tmp_path, report):
    ref = write(tmp_path / "report.json", report)
    return reader.verify(ref["path"], ref["sha256"])


@pytest.mark.parametrize("category,lcas", [("candidate_separation", 1), ("positive_paralogy_exclusion", 2)])
def test_separation_and_same_family_exclusion(tmp_path, category, lcas):
    result = run(tmp_path, fixture(tmp_path, category))
    assert result["changed_pairs_checked"] == 1 and result["tree_leaves_checked"] == 8
    assert result["source_families_checked"] == (3 if category == "candidate_separation" else 2)
    assert result["distinct_lcas_checked"] == lcas and result["summary"] == {category: 1}
    assert result["species_tree_bytes_identical"] is False
    assert result["original_node_annotation_admission_established"] is False
    assert all(result[key] is False for key in reader.SCOPE_FLAGS)


@pytest.mark.parametrize("fault", ["source", "scope", "inventoried", "cells", "helper", "count", "duplicate",
    "transition", "family_count", "manifest", "member_hash", "member_size", "root", "gene", "family",
    "same_family", "invented_lca", "prediction", "lca", "event", "equality", "label", "category", "species", "summary"])
def test_changed_trace_refused(tmp_path, fault):
    report = fixture(tmp_path)
    pair = report["changed_pairs"][0]
    if fault == "source":
        report["source"] = reader.record(__file__)
    elif fault == "scope":
        report["scientific_timings_admitted"] = True
    elif fault == "inventoried":
        report["node_annotations_previously_inventoried"] = True
    elif fault == "cells":
        report["states"].reverse()
    elif fault == "helper":
        report["helpers"].reverse()
    elif fault == "count":
        report["changed_pairs_traced"] = 2
    elif fault == "duplicate":
        report["changed_pairs"].append(deepcopy(pair))
        report["changed_pairs_traced"] = 2
    elif fault == "transition":
        report["comparison"]["families"][0]["transitions"]["TP->FN"] = 2
    elif fault == "family_count":
        report["states"][1]["selected_source_families"] = 1
    elif fault == "manifest":
        report["states"][1]["species_tree_sha256"] = "different"
    elif fault == "member_hash":
        pair["after"]["candidate_members_a_sha256"] = "0" * 64
    elif fault == "member_size":
        pair["after"]["candidate_size_a"] = 1
    elif fault == "root":
        pair["after"]["root_hog_a"] = "RootHOG0000001"
    elif fault == "gene":
        pair["after"]["gene_a"] = "x"
    elif fault == "family":
        pair["after"]["source_family_a"] = "Family0000009"
    elif fault == "same_family":
        pair["after"]["same_source_family"] = True
    elif fault == "invented_lca":
        pair["after"]["lca"] = dict(node_id="G00000", pair_event="duplication")
    elif fault == "prediction":
        pair["after"]["predicted"] = True
    elif fault == "lca":
        pair["before"]["lca"]["node_id"] = "G99999"
    elif fault == "event":
        pair["before"]["lca"]["pair_event"] = "duplication"
    elif fault == "equality":
        pair["candidate_members_a_unchanged"] = True
    elif fault == "label":
        pair["after_label"] = "TN"
        report["comparison"]["families"][0]["transitions"] = {"TP->TN": 1}
    elif fault == "category":
        pair["localization"] = "positive_paralogy_exclusion"
    elif fault == "species":
        report["species_tree_bytes_identical"] = True
    elif fault == "summary":
        report["summary"] = {"candidate_separation": 2}
    with pytest.raises(ValueError):
        run(tmp_path, report)


@pytest.mark.parametrize("fault", ["checkpoint_species", "tree_bytes", "tree_members", "root_duplicate", "root_missing"])
def test_saved_evidence_refused(tmp_path, fault):
    report = fixture(tmp_path)
    state = report["states"][0]
    directory = Path(state["root_hogs"]["path"]).parent
    checkpoint = directory / "checkpoints/Family0000000.json"
    if fault == "checkpoint_species":
        value = json.loads(checkpoint.read_text())
        value["species_tree_sha256"] = "different"
        write(checkpoint, value)
    elif fault == "tree_bytes":
        (directory / "gene_trees/Family0000000.raw.nwk").write_text("(x,y);\n")
    elif fault == "tree_members":
        trees(directory, "Family0000000", "abcx", "(a,b,c,x)'G00000|S@S0';\n")
    else:
        path = Path(state["root_hogs"]["path"])
        path.write_text("root_hog\tsource_family\tgenes\nRootHOG0000000\tFamily0000000\t"
                       + ("a,a,b,c,d\n" if fault == "root_duplicate" else "a,b,c\n"))
        changed = reader.record(path)
        report["evidence"][report["evidence"].index(state["root_hogs"])] = changed
        state["root_hogs"] = changed
    with pytest.raises(ValueError):
        run(tmp_path, report)


def test_changed_digest_and_existing_output_refused(tmp_path, monkeypatch):
    report = fixture(tmp_path)
    ref = write(tmp_path / "report.json", report)
    with pytest.raises(ValueError, match="Changed profile localization report"):
        reader.verify(ref["path"], "0" * 64)
    path = tmp_path / "retain.json"
    path.write_text("retain\n")
    monkeypatch.setattr(sys, "argv", ["reader", "--report", ref["path"], "--report-sha256", ref["sha256"], "--output", str(path)])
    with pytest.raises(ValueError, match="already exists"):
        reader.main()
    assert path.read_text() == "retain\n"
