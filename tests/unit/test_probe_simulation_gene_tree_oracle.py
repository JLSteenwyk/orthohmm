"""Scientific-scope and input guards for the fixed-candidate oracle control."""

import json
from types import SimpleNamespace

import dendropy
import pytest

from benchmark_tools import probe_simulation_gene_tree_oracle as probe
from orthohmm import phylogeny


def tree(newick):
    return dendropy.Tree.get(data=newick, schema="newick", preserve_underscores=True,
                             rooting="force-rooted")


def test_fixed_panel_contains_every_condition_and_seed():
    rows = [{"condition": c, "seed": s, "method": probe.METHOD, "variant": "generating"}
            for c in probe.CONDITIONS for s in range(20261101, 20261111)]
    assert len(probe.select_cells({"records": rows})) == 70
    assert len(probe.select_cells({"records": rows + [{**rows[0], "variant": "nni1"}]})) == 70


@pytest.mark.parametrize("change", ["missing", "duplicate", "seed", "condition"])
def test_incomplete_or_changed_fixed_panel_rejects(change):
    rows = [{"condition": c, "seed": s, "method": probe.METHOD, "variant": "generating"}
            for c in probe.CONDITIONS for s in range(20261101, 20261111)]
    if change == "missing":
        rows.pop()
    elif change == "duplicate":
        rows[-1] = rows[0]
    else:
        rows[0][change] = 1 if change == "seed" else "unfrozen"
    with pytest.raises(ValueError, match="70-cell"):
        probe.select_cells({"records": rows})


def test_partition_preserves_native_family_numbering(tmp_path):
    path = tmp_path / "groups"
    path.write_text("b a\nc\n")
    assert probe.partition(path, {"a", "b", "c"}) == {
        "Family0000000": {"a", "b"}, "Family0000001": {"c"}}


@pytest.mark.parametrize("text", ["a a\nb\n", "a\na b\n", "a\n\nb\n", "a\n", "a d\nb\n"])
def test_partition_rejects_invalid_complete_membership(tmp_path, text):
    path = tmp_path / "groups"
    path.write_text(text)
    with pytest.raises(ValueError):
        probe.partition(path, {"a", "b"})


def test_induction_preserves_gene_labels_root_and_distances(tmp_path):
    path = tmp_path / "generating.nwk"
    path.write_text("(((a_1:1,b_1:1)X:1,c_1:2)Y:1,d_1:3)Root:0;")
    parent = {"F7__" + g for g in ("a_1", "b_1", "c_1", "d_1")}
    genes = parent - {"F7__d_1"}
    induced = probe.induced_tree(path, "7", parent, genes)
    assert probe.tree_labels(induced) == genes
    assert all(n.label is None for n in induced.preorder_node_iter() if not n.is_leaf())
    distances = induced.phylogenetic_distance_matrix()
    a, b = [induced.taxon_namespace.get_taxon(label="F7__" + g) for g in ("a_1", "c_1")]
    assert distances(a, b) == 4
    expected = tree("((F7__a_1,F7__b_1),F7__c_1);")
    assert probe.distances(induced, expected)["rooted_clade_distance"] == 0


@pytest.mark.parametrize("parent,genes", [({"F7__a_1"}, {"F7__a_1"}),
    ({"F7__a_1", "F7__b_1"}, {"F7__a_1", "F7__c_1"})])
def test_induction_rejects_parent_or_candidate_mismatch(tmp_path, parent, genes):
    path = tmp_path / "generating.nwk"
    path.write_text("(a_1,b_1);")
    with pytest.raises(ValueError, match="parent family"):
        probe.induced_tree(path, "7", parent, genes)


def test_root_disagreement_is_not_unrooted_topology_disagreement():
    first = tree("((a,b),(c,d));")
    other_root = tree("(a,(b,(c,d)));")
    assert probe.distances(first, other_root) == {
        "rooted_clade_distance": 2, "unrooted_split_distance": 0}
    assert probe.distances(first, tree("((a,c),(b,d));"))["unrooted_split_distance"] == 2


def test_distance_comparison_rejects_different_gene_universes():
    with pytest.raises(ValueError, match="identical"):
        probe.distances(tree("(a,b);"), tree("(a,c);"))


def test_no_constraint_policy_does_not_filter_root_groups():
    value = SimpleNamespace(ortholog_pairs=(("a", "b"),), root_groups=(("a",), ("b",)))
    assert probe.constrained_pairs(value, [], {"a", "b"}, False) == {("a", "b")}
    with pytest.raises(ValueError, match="without native"):
        probe.constrained_pairs(value, [(0, {})], {"a", "b"}, False)


def test_active_constraint_policy_filters_root_groups_even_without_local_event():
    value = SimpleNamespace(ortholog_pairs=(("a", "b"),), root_groups=(("a",), ("b",)),
                            ortholog_pair_confidence=(("a", "b", "high"),))
    assert probe.constrained_pairs(value, [], {"a", "b"}, True) == set()


def test_unsupported_satellite_is_detached():
    value = SimpleNamespace(ortholog_pairs=(("a", "b"),), root_groups=(("a", "b"),),
                            ortholog_pair_confidence=(("a", "b", "medium"),))
    events = [(4, {"source_genes": ["a"], "target_genes": ["b"]})]
    assert probe.constrained_pairs(value, events, {"a", "b"}, True) == set()
    value.ortholog_pair_confidence = (("a", "b", "high"),)
    assert probe.constrained_pairs(value, events, {"a", "b"}, True) == {("a", "b")}


def test_input_checks_do_not_accept_equal_size_changed_content(tmp_path):
    path = tmp_path / "source"
    path.write_text("a")
    item = probe.record(path)
    path.write_text("b")
    with pytest.raises(ValueError, match="Input changed"):
        probe.checked(item, {})


def test_frozen_source_is_checked_before_import(tmp_path):
    path = tmp_path / "phylogeny.py"
    path.write_text("raise AssertionError('must not import')")
    with pytest.raises(ValueError, match="source changed"):
        probe.frozen_module(path, {})


def test_unavailable_cell_is_preserved_without_reading_artifacts():
    row = {"status": "failed", "label": "baseline_20261101", "condition": "baseline", "seed": 20261101}
    result = probe.cell(row, {}, None, {})
    assert result["status"] == "native_unavailable" and result["native_record"] == row


@pytest.mark.parametrize("mixed", [False, True])
def test_complete_cell_reproduces_baseline_and_retains_mixed_candidates(tmp_path, mixed):
    cell_dir = tmp_path / "cell"
    native = cell_dir / probe.METHOD
    inputs = {}
    def write(path, value):
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(value) if not isinstance(value, str) else value)
        return probe.record(path)
    genes = {"F1__a_1", "F1__a_2", "F2__b_1" if mixed else "F1__b_1"}
    owners = {g: g.split("__")[1].split("_")[0] for g in genes}
    fasta = [write(tmp_path / "input" / f"{s}.fasta", "".join(f">{g}\nAAA\n" for g in sorted(genes) if owners[g] == s)) for s in ("a", "b")]
    species_record = write(native / "orthohmm_phylogeny/species_tree.rooted.nwk", "(a,b);")
    g_b = "F2__b_1" if mixed else "F1__b_1"
    rooted = f"(F1__a_1,(F1__a_2,{g_b}));"
    raw_record = write(native / "orthohmm_phylogeny/gene_trees/Family0000000.raw.nwk", rooted)
    rooted_record = write(native / "orthohmm_phylogeny/gene_trees/Family0000000.rooted.nwk", rooted)
    value = phylogeny.reconcile_gene_tree(tree(rooted), tree("(a,b);"), owners,
        family_id="Family0000000", root_duplication_rule="species_overlap", pair_orthology_rule="positive_paralogy")
    native_pairs = "gene_a\tspecies_a\tgene_b\tspecies_b\n" + "".join(
        f"{a}\t{owners[a]}\t{b}\t{owners[b]}\n" for a, b in value.ortholog_pairs)
    prediction = write(native / "orthohmm_phylogeny/orthohmm_pairwise_orthologs.tsv", native_pairs)
    candidates = write(native / "orthohmm_working_res/phylogeny_candidate_superfamilies.txt", " ".join(sorted(genes)) + "\n")
    checkpoint = write(native / "orthohmm_phylogeny/checkpoints/Family0000000.json", {
        "status": "complete", "family_id": "Family0000000", "genes": sorted(genes),
        "raw_tree_sha256": raw_record["sha256"], "rooted_tree_sha256": rooted_record["sha256"],
        "species_tree_sha256": species_record["sha256"]})
    manifest = write(native / "orthohmm_phylogeny/provenance_manifest.json", {
        "root_duplication_rule": "species_overlap", "pair_orthology_rule": "positive_paralogy",
        "membership_reconciliation": None, "species_tree_sha256": species_record["sha256"]})
    ancestral = {"1": sorted(g for g in genes if g.startswith("F1__"))}
    if mixed:
        ancestral["2"] = [g_b]
    true_tree = write(tmp_path / "parent/G/Gene_trees/1_prunedtree.nwk", "(a_1,(a_2,b_1));")
    truth = write(tmp_path / "truth.json", {"families": ancestral, "extant_genes": len(genes),
        "ortholog_pairs": list(value.ortholog_pairs), "inputs": [{**true_tree, "path": "G/Gene_trees/1_prunedtree.nwk"}]})
    original = write(cell_dir / "results.json", {})
    preflight = write(cell_dir / "preflight.json", {})
    outputs = [species_record, raw_record, rooted_record, candidates, checkpoint, manifest]
    execution = write(cell_dir / "execution.json", {"methods": {probe.METHOD: {
        "status": "process_succeeded", "exit_code": 0, "outputs": [{**r, "absolute_path": r["path"]} for r in outputs]}}})
    row = {"status": "admitted", "label": "baseline_20261101", "condition": "baseline", "seed": 20261101,
        "native_report": original, "preflight": preflight, "execution": execution, "truth": truth,
        "prediction_files": [prediction], "tree": species_record}
    dataset = {"condition": "baseline", "seed": 20261101, "parent": "baseline_20261101",
        "input_evidence": {"truth": truth, "inputs": fasta}, "generating_tree": {"path": str(tmp_path / "parent/T/ExtantTree.nwk")}}
    result = probe.cell(row, {"datasets": [dataset]}, phylogeny, inputs)
    expected = "mixed_ancestry_ineligible" if mixed else "oracle_eligible"
    assert result["candidate_status_counts"] == {expected: 1}
    assert result["arms"]["inferred"]["predicted_pairs"] == len(value.ortholog_pairs)
    if mixed:
        assert len({json.dumps(value, sort_keys=True) for value in result["arms"].values()}) == 1
        assert result["changes"]["generating_root"] == {"added": 0, "removed": 0}
