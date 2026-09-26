from collections import Counter
import json

import pytest

from benchmark_tools.audit_phylogeny_hierarchy import (
    HIERARCHY_COLUMNS, audit, expected_hierarchy, needs_tree,
)
from benchmark_tools.derive_phylogeny_events import NODE_COLUMNS
from benchmark_tools.prepare_ob_candidate_neighborhood import record


@pytest.mark.parametrize("species,expected", [(["A"], False), (["A", "B"], False),
    (["A", "A", "A"], False), (["A", "B", "C"], False), (["A", "A", "B"], True)])
def test_reconciliation_selection(species, expected):
    mapping = {str(i): taxon for i, taxon in enumerate(species)}
    assert needs_tree(mapping.keys(), mapping) is expected


def test_expected_internal_and_bypass_rows():
    families = {"F": {"a", "b", "c"}, "Z": {"d"}}
    mapping = {"a": "A", "b": "A", "c": "B", "d": "C"}
    node_rows = [dict(source_family="F", event="leaf"),
                 dict(source_family="F", event="duplication", node_id="G0",
                      parent_node_id="G1", species_tree_node="S0", genes="a,b"),
                 dict(source_family="F", event="speciation", node_id="G1",
                      parent_node_id="", species_tree_node="S0", genes="a,b,c")]
    totals = Counter()
    result = list(expected_hierarchy(families, mapping, node_rows, totals))
    assert [row["hog_id"] for row in result] == ["F.G0", "F.G1", "Z.root"]
    assert result[0]["parent_hog_id"] == "F.G1"
    assert result[2]["event"] == "unambiguous"
    assert totals == dict(reconciled_families=1, bypassed_families=1, internal_node_rows=2, bypass_rows=1)


@pytest.mark.parametrize("families,mapping,records", [
    ({"F": {"a", "b", "c"}}, {"a": "A", "b": "A", "c": "B"}, []),
    ({"F": {"a"}}, {"a": "A"}, [dict(source_family="F", event="leaf")]),
])
def test_wrong_selection(families, mapping, records):
    with pytest.raises(ValueError):
        list(expected_hierarchy(families, mapping, records, Counter()))


@pytest.fixture
def bypass_fixture(tmp_path):
    directory = tmp_path / "phylo"
    directory.mkdir()
    fasta = tmp_path / "A.fa"
    fasta.write_text(">a\nACD\n>b\nACE\n")
    roots = directory / "orthohmm_root_hogs.tsv"
    roots.write_text("root_hog\tsource_family\tgenes\nH\tF\ta,b\n")
    nodes = directory / "orthohmm_reconciliation_nodes.tsv"
    nodes.write_text("\t".join(NODE_COLUMNS) + "\n")
    manifest = directory / "provenance_manifest.json"
    totals = dict(reconciled_families=0, bypassed_families=1)
    manifest.write_text(json.dumps(dict(results=totals, input_proteomes=[
        dict(filename="A.fa", taxon="A", sha256=record(fasta)["sha256"])])))
    hierarchy = directory / "orthohmm_hierarchical_orthogroups.tsv"
    hierarchy.write_text("\t".join(HIERARCHY_COLUMNS) + "\nF.root\t\t\tF\tunambiguous\ta,b\n")
    report = tmp_path / "events.json"
    report.write_text(json.dumps(dict(status="frozen_phylogeny_event_pair_semantics_verified",
        totals=totals, checked_records=[record(path) for path in (fasta, roots, nodes, manifest)])))
    return directory, report


def test_hierarchy_readback(bypass_fixture):
    result = audit(*bypass_fixture)
    assert result["totals"]["hierarchical_groups"] == 1
    assert not result["scientific_scores_admitted"]


@pytest.mark.parametrize("mutation", ["event", "parent", "members", "missing", "extra", "schema"])
def test_hierarchy_mutations(bypass_fixture, mutation):
    directory, report = bypass_fixture
    path = directory / "orthohmm_hierarchical_orthogroups.tsv"
    text = path.read_text()
    if mutation == "event":
        text = text.replace("unambiguous", "speciation")
    elif mutation == "parent":
        text = text.replace("F.root\t", "F.root\tF.parent", 1)
    elif mutation == "members":
        text = text.replace("a,b", "a")
    elif mutation == "missing":
        text = text.splitlines()[0] + "\n"
    elif mutation == "extra":
        text += text.splitlines()[1] + "\n"
    else:
        text = text.replace("hog_id", "bad_id", 1)
    path.write_text(text)
    with pytest.raises(ValueError):
        audit(directory, report)
