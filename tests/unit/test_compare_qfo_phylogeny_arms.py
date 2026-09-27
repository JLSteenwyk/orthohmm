import pytest

from benchmark_tools.compare_qfo_phylogeny_arms import species_tree_signature, pair_family_contrast


def test_species_topology_ignores_serialization_but_detects_rooting(tmp_path):
    a, b, c = [tmp_path / name for name in ("a", "b", "c")]
    a.write_text("[&R] ((A:1,B:2):3,(C:4,D:5):6);\n")
    b.write_text("[&R] ((D:9,C:8):7,(B:6,A:5):4);\n")
    c.write_text("[&R] (A,(B,(C,D)));\n")
    assert species_tree_signature(a) == species_tree_signature(b)
    assert species_tree_signature(a) != species_tree_signature(c)


@pytest.mark.parametrize("text", ["[&U] ((A,B),(C,D));", "[&R] (A,B);[&R] (C,D);"])
def test_invalid_species_tree(tmp_path, text):
    path = tmp_path / "tree"
    path.write_text(text)
    with pytest.raises(ValueError):
        species_tree_signature(path)


def test_pair_family_attribution():
    a = [(("a", "b"), "high"), (("a", "c"), "low")]
    b = [(("a", "b"), "low"), (("b", "c"), "high")]
    result = pair_family_contrast(a, b, dict.fromkeys("abc", "F1"), dict.fromkeys("abc", "F2"))
    assert result["shared"] == 1
    assert result["left_only"] == result["right_only"] == result["annotations_changed"] == 1
    assert result["affected_families"] == [
        dict(side="left", family="F1", unique_pairs=1, annotation_changes=1),
        dict(side="right", family="F2", unique_pairs=1, annotation_changes=1)]


def test_shared_pair_relabel_is_not_prediction_difference():
    pairs = [(("a", "b"), "high")]
    result = pair_family_contrast(pairs, pairs, dict.fromkeys("ab", "F1"), dict.fromkeys("ab", "F2"))
    assert result["confidence_tables_equal"]
    assert result["affected_families"] == []


def test_cross_family_pairs_rejected():
    with pytest.raises(ValueError):
        pair_family_contrast([(("a", "b"), "high")], [], dict(a="F1", b="F2"), {})
