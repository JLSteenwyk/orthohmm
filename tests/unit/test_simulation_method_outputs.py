import pytest

from benchmark_tools.simulation_method_outputs import load_predictions, orthofinder_pairs, orthohmm_pairs


def of_tables(root, empty=False):
    for a, b, ga, gb in (("A", "B", "a", "b"), ("B", "A", "b", "a")):
        path = root / f"Orthologues/Orthologues_{a}/{a}__v__{b}.tsv"
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(f"Orthogroup\t{a}\t{b}\n" + ("" if empty else f"OG0\t{ga}\t{gb}\n"))


def test_complete_native_of_and_empty_tables(tmp_path):
    owners = {"a": "A", "b": "B"}
    of_tables(tmp_path)
    assert orthofinder_pairs(tmp_path, owners, ["A", "B"]) == [("a", "b")]
    assert load_predictions("orthofinder_full", tmp_path, owners, ["A", "B"])[0] == [("a", "b")]
    of_tables(tmp_path, empty=True)
    assert orthofinder_pairs(tmp_path, owners, ["A", "B"]) == []


@pytest.mark.parametrize("fault", ["missing", "orientation", "unknown", "species", "empty_endpoint"])
def test_of_partial_or_inconsistent_output_is_not_an_empty_success(tmp_path, fault):
    of_tables(tmp_path)
    reverse = tmp_path / "Orthologues/Orthologues_B/B__v__A.tsv"
    if fault == "missing":
        reverse.unlink()
    elif fault == "orientation":
        reverse.write_text("Orthogroup\tB\tA\n")
    elif fault == "unknown":
        reverse.write_text("Orthogroup\tB\tA\nOG0\tunknown\ta\n")
    elif fault == "species":
        reverse.write_text("Orthogroup\tB\tA\nOG0\ta\tb\n")
    else:
        reverse.write_text("Orthogroup\tB\tA\nOG0\t\ta\n")
    with pytest.raises(ValueError):
        orthofinder_pairs(tmp_path, {"a": "A", "b": "B"}, ["A", "B"])


def test_orthohmm_native_ids_and_headers(tmp_path):
    path = tmp_path / "pairs.tsv"
    path.write_text("gene_a\tspecies_a\tgene_b\tspecies_b\na\tA\tb\tB\n")
    assert list(orthohmm_pairs(path, {"a": "A", "b": "B"})) == [("a", "b")]
    path.write_text("gene_a\tspecies_a\tgene_b\tspecies_b\na\tB\tb\tA\n")
    with pytest.raises(ValueError):
        list(orthohmm_pairs(path, {"a": "A", "b": "B"}))


def test_high_sensitivity_cliques_and_no_same_species_pairs(tmp_path):
    (tmp_path / "orthohmm_orthogroups.txt").write_text("OG0: a aa b\n")
    pairs, sources = load_predictions("orthohmm_high_sensitivity", tmp_path, {"a": "A", "aa": "A", "b": "B"}, ["A", "B"])
    assert set(pairs) == {("a", "b"), ("aa", "b")}
    assert len(sources) == 1


def test_sequence_checkpoint_uses_complete_id_mapping(tmp_path):
    (tmp_path / "SequenceIDs.txt").write_text("0_0: a\n1_0: b\n")
    (tmp_path / "clusters_OrthoFinder_I1.5.txt_id_pairs.txt").write_text("begin\n0 0_0 1_0 $\n)\n")
    assert load_predictions("orthofinder_sequence_only", tmp_path, {"a": "A", "b": "B"}, ["A", "B"])[0] == [("a", "b")]
    (tmp_path / "second").mkdir()
    (tmp_path / "second/SequenceIDs.txt").write_text("0_0: a\n")
    with pytest.raises(ValueError, match="exactly one"):
        load_predictions("orthofinder_sequence_only", tmp_path, {"a": "A", "b": "B"}, ["A", "B"])
