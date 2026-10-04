"""Independent original-history/rule consistency checks, including tampering."""

from copy import deepcopy
import json
import xml.etree.ElementTree as ET

from Bio import Phylo
import dendropy
from io import StringIO
import pytest

from benchmark_tools import readback_simulation_oracle_residuals as readback
from benchmark_tools import trace_simulation_oracle_residuals as primary
from orthohmm import phylogeny


def dtree(text):
    return dendropy.Tree.get(data=text, schema="newick", preserve_underscores=True, rooting="force-rooted")


@pytest.mark.parametrize("relative", ["prepared/a/truth.json", "another/root/truth.json"])
def test_retained_absolute_path_takes_precedence_over_manifest_relative_path(tmp_path, relative):
    actual = tmp_path / "truth.json"
    assert readback.record_path({"path": relative, "absolute_path": str(actual)}) == actual


def test_records_without_absolute_path_resolve_recorded_path(tmp_path):
    assert readback.record_path({"path": str(tmp_path / "truth.json")}) == tmp_path / "truth.json"


def test_empty_optional_absolute_path_falls_back_to_recorded_path(tmp_path):
    actual = tmp_path / "truth.json"
    assert readback.record_path({"path": str(actual), "absolute_path": ""}) == actual


def fixture(tmp_path, bypass=False, constraint=False):
    last = "a_2" if constraint else "d_1"
    graph = {"Root_1": ("D", ("Left_1", "Right_1")),
             "Left_1": ("S", ("a_1", "b_1")), "Right_1": ("S", ("c_1", last)),
             **{g: ("F", ()) for g in ("a_1", "b_1", "c_1", last)}}
    codes = {"S": "speciation", "D": "duplication", "F": "P"}
    root = ET.Element("recGeneTree")
    phy = ET.SubElement(root, "phylogeny")
    def clade(parent, name):
        element = ET.SubElement(parent, "clade")
        ET.SubElement(element, "name").text = name
        annotation = ET.SubElement(element, "eventsRec")
        event, children = graph[name]
        ET.SubElement(annotation, codes[event])
        for child in children:
            clade(element, child)
    clade(phy, "Root_1")
    xml_path = tmp_path / "original.xml"
    ET.ElementTree(root).write(xml_path)
    genes = {"F1__" + g for g in ("a_1", "b_1", "c_1", last)}
    owners = {g: g.removeprefix("F1__").rsplit("_", 1)[0] for g in genes}
    truth = {("F1__a_1", "F1__b_1"), tuple(sorted(("F1__c_1", "F1__" + last)))}
    species_text = "(a,(b,c));" if constraint else "((a,b),(c,d));"
    event = [(5, {"source_genes": ["F1__a_1"], "target_genes": ["F1__a_2"]})] if constraint else []
    result = primary.trace_candidate(phylogeny, dtree(species_text),
        dtree(f"((F1__a_1,F1__b_1),(F1__c_1,F1__{last}));"), graph, "1", genes, owners,
        owners, truth, event, constraint, bypass, "Family0000000")
    row = json.loads(json.dumps({"cell": "baseline_20261101", "family": "Family0000000",
        "ancestor": "1", "genes": sorted(genes), "status": "unambiguous_bypass" if bypass else "oracle_eligible",
        "native_constraint_policy_active": constraint, **result}))
    expected = {"status": row["status"], "ancestral_families": ["1"], "arms": {"generating_root": row["counts"]}}
    return row, expected, readback.xml_lineages(xml_path), owners, Phylo.read(StringIO(species_text), "newick"), event, genes, truth


def check(data):
    row, expected, (paths, descendants, events), owners, species, constraints, genes, truth = data
    return readback.check_candidate(row, expected, paths, descendants, events, owners,
        species, constraints, bool(constraints), genes, truth)


@pytest.mark.parametrize("bypass", [False, True])
def test_independent_xml_and_rule_readback_reproduces_complete_pair_cohort(tmp_path, bypass):
    result = check(fixture(tmp_path, bypass=bypass))
    assert result["pair_rows_verified"] == 6
    assert result["counts"] == {"tp": 2, "fp": 4, "fn": 0}


def test_independent_native_constraint_replay_identifies_true_pair_loss(tmp_path):
    result = check(fixture(tmp_path, constraint=True))
    assert result["counts"] == {"tp": 1, "fp": 0, "fn": 1}
    assert result["error_classes"] == {"unsupported_satellite_constraint": 1}


@pytest.mark.parametrize("change", ["pair_missing", "pair_duplicate", "flag", "history", "clade", "parent",
    "call", "overlap", "root_groups", "pair_node", "class", "aggregate", "eligibility", "support"])
def test_independent_readback_rejects_changed_pair_or_tree_evidence(tmp_path, change):
    data = fixture(tmp_path)
    row = data[0]
    if change == "pair_missing":
        row["pairs"].pop()
    elif change == "pair_duplicate":
        row["pairs"].append(deepcopy(row["pairs"][0]))
    elif change == "flag":
        row["pairs"][0]["raw_predicted"] = False
    elif change == "history":
        row["pairs"][0]["history_node"] = "other"
    elif change == "clade":
        row["reconciliation_nodes"][0]["genes"] = ["unknown"]
    elif change == "parent":
        row["reconciliation_nodes"][0]["parent_node_id"] = None
    elif change == "call":
        next(n for n in row["reconciliation_nodes"] if n["event"] != "leaf")["pair_event"] = "duplication"
    elif change == "overlap":
        row["pairs"][0]["candidate_species_overlap"] = ["a"]
    elif change == "root_groups":
        row["root_groups"] = [[g] for g in row["genes"]]
    elif change == "pair_node":
        row["pairs"][0]["pair_node"]["node_id"] = "other"
    elif change == "class":
        row["pairs"][0]["error_class"] = "other"
    elif change == "aggregate":
        row["counts"]["tp"] += 1
    elif change == "eligibility":
        row["status"] = "mixed_ancestry_ineligible"
    else:
        next(n for n in row["reconciliation_nodes"] if n["event"] != "leaf")["branch_support"] = 0.99
    with pytest.raises(ValueError):
        check(data)


def test_independent_readback_checks_original_constraint_support_not_reported_boolean(tmp_path):
    data = fixture(tmp_path, constraint=True)
    data[0]["constraint_evidence"][0]["supported"] = True
    with pytest.raises(ValueError, match="constraint evidence"):
        check(data)


@pytest.mark.parametrize("text", ["<root/>", "<recGeneTree><phylogeny/></recGeneTree>",
    "<recGeneTree><phylogeny><clade><name>Root_1</name><eventsRec><duplication/></eventsRec></clade></phylogeny></recGeneTree>"])
def test_original_xml_rejects_missing_root_or_invalid_event_arity(tmp_path, text):
    path = tmp_path / "bad.xml"
    path.write_text(text)
    with pytest.raises(ValueError):
        readback.xml_lineages(path)
