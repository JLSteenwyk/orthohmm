import json

import pytest

from benchmark_tools.audit_phylogeny_events import audit
from benchmark_tools.audit_phylogeny_structure import audit as structural_audit
from benchmark_tools.derive_phylogeny_events import NODE_COLUMNS
from benchmark_tools.prepare_ob_candidate_neighborhood import record
from tests.unit.test_audit_phylogeny_structure import tree_fixture  # noqa: F401


@pytest.fixture
def semantic_fixture(tree_fixture):
    directory, inputs = tree_fixture
    path = directory / "gene_trees/F.reconciled.nwk"
    path.write_text("[&R] (a:1,b:1)G00000|S@S0000;\n")
    cp = directory / "checkpoints/F.json"
    checkpoint = json.loads(cp.read_text())
    checkpoint["annotated_tree_sha256"] = record(path)["sha256"]
    cp.write_text(json.dumps(checkpoint))
    summary_path = directory / "reconciliation_summary.json"
    summary = json.loads(summary_path.read_text())
    summary.update(duplications=0, speciations=1, uncertain_events=0)
    summary_path.write_text(json.dumps(summary))
    path = directory / "provenance_manifest.json"
    manifest = json.loads(path.read_text())
    manifest.update(results=summary, root_duplication_rule="species_overlap",
                    pair_orthology_rule="positive_paralogy", membership_reconciliation=None)
    path.write_text(json.dumps(manifest))
    node_rows = [
        ["F", "a", "G00000", "leaf", "S1", "S1", "a", "leaf", "not_applicable", "0", "false", ""],
        ["F", "b", "G00000", "leaf", "S2", "S2", "b", "leaf", "not_applicable", "0", "false", ""],
        ["F", "G00000", "", "speciation", "S0000", "S1,S2", "a,b", "speciation", "high", "0", "false", ""],
    ]
    (directory / "orthohmm_reconciliation_nodes.tsv").write_text(
        "\t".join(NODE_COLUMNS) + "\n" + "\n".join("\t".join(row) for row in node_rows) + "\n")
    return directory, inputs


def evidence(directory, inputs):
    report = directory.parent / "structure.json"
    report.write_text(json.dumps(structural_audit(directory, inputs)))
    return report


def test_full_event_readback(semantic_fixture):
    directory, inputs = semantic_fixture
    result = audit(directory, evidence(directory, inputs))
    assert result["totals"]["nodes"] == 3
    assert result["totals"]["ortholog_pairs"] == 1
    assert result["totals"]["speciations"] == 1


@pytest.mark.parametrize("mutation", ["node", "pair_missing", "confidence", "groups", "annotation", "rule"])
def test_semantic_mutations_pass_structure_but_fail_oracle(semantic_fixture, mutation):
    directory, inputs = semantic_fixture
    summary_path = directory / "reconciliation_summary.json"
    manifest_path = directory / "provenance_manifest.json"
    summary = json.loads(summary_path.read_text())
    manifest = json.loads(manifest_path.read_text())
    if mutation == "node":
        path = directory / "orthohmm_reconciliation_nodes.tsv"
        path.write_text(path.read_text().replace("speciation", "duplication"))
    elif mutation == "pair_missing":
        for name in ("orthohmm_pairwise_orthologs.tsv", "orthohmm_pairwise_orthologs_confidence.tsv"):
            path = directory / name
            path.write_text(path.read_text().splitlines()[0] + "\n")
        summary["ortholog_pairs"] = 0
    elif mutation == "confidence":
        path = directory / "orthohmm_pairwise_orthologs_confidence.tsv"
        path.write_text(path.read_text().replace("high", "medium"))
    elif mutation == "groups":
        path = directory / "orthohmm_root_hogs.tsv"
        path.write_text("root_hog\tsource_family\tgenes\nH1\tF\ta\nH2\tF\tb\n")
        summary["root_hogs"] = 2
    elif mutation == "annotation":
        path = directory / "gene_trees/F.reconciled.nwk"
        path.write_text(path.read_text().replace("|S@", "|D@"))
        cp = directory / "checkpoints/F.json"
        data = json.loads(cp.read_text())
        data["annotated_tree_sha256"] = record(path)["sha256"]
        cp.write_text(json.dumps(data))
    else:
        manifest["pair_orthology_rule"] = "other"
    summary_path.write_text(json.dumps(summary))
    manifest["results"] = summary
    manifest_path.write_text(json.dumps(manifest))
    report = evidence(directory, inputs)
    with pytest.raises(ValueError):
        audit(directory, report)
