import csv
import json

import pytest

from benchmark_tools.compare_native_membership import compare
from benchmark_tools import compare_native_membership as diagnostic


def native(folder, groups, pairs):
    folder.mkdir()
    with (folder / "orthohmm_root_hogs.tsv").open("w", newline="") as stream:
        writer = csv.writer(stream, delimiter="\t")
        writer.writerow(["root_hog", "source_family", "genes"])
        for name, members in groups.items():
            writer.writerow([name, "arbitrary", ",".join(members)])
    with (folder / "orthohmm_pairwise_orthologs.tsv").open("w", newline="") as stream:
        writer = csv.writer(stream, delimiter="\t")
        writer.writerow(["gene_a", "species_a", "gene_b", "species_b"])
        writer.writerows(pairs)
    (folder / "reconciliation_summary.json").write_text(json.dumps({"root_hogs": len(groups), "ortholog_pairs": len(pairs)}))
    return folder


def test_serialization_group_ids_and_pair_orientation_do_not_change_semantics(tmp_path):
    a = native(tmp_path / "a", {"root1": ["a", "b"], "root2": ["c"]},
               [("a", "species1", "b", "species2")])
    b = native(tmp_path / "b", {"different2": ["c"], "different1": ["b", "a"]},
               [("b", "species2", "a", "species1")])
    result = compare(a, b, 3)
    assert result["root_membership"]["semantic_multisets_equal"]
    assert result["pair_membership"]["semantic_multisets_equal"]
    assert result["pair_membership"]["affected_genes"] == 0
    assert result["accuracy_evaluated"] is False


def test_replacement_changes_memberships_despite_equal_record_counts(tmp_path):
    a = native(tmp_path / "a", {"one": ["a", "b"], "two": ["c"]}, [("a", "s1", "b", "s2")])
    b = native(tmp_path / "b", {"one": ["a", "c"], "two": ["b"]}, [("a", "s1", "c", "s3")])
    result = compare(a, b, 3)
    assert result["root_membership"]["removed_records"] == result["root_membership"]["added_records"] == 2
    assert result["pair_membership"]["removed_records"] == result["pair_membership"]["added_records"] == 1
    assert result["pair_membership"]["affected_genes"] == 3


def test_duplicate_pairs_are_not_silently_deduplicated(tmp_path):
    pair = ("a", "s1", "b", "s2")
    a = native(tmp_path / "a", {"one": ["a", "b"]}, [pair, pair])
    b = native(tmp_path / "b", {"one": ["a", "b"]}, [pair])
    result = compare(a, b, 2)
    assert result["pair_membership"]["baseline_duplicate_rows"] == 1
    assert result["pair_membership"]["removed_records"] == 1
    assert result["pair_membership"]["common_records"] == 1


@pytest.mark.parametrize("pair", [("a", "s1", "b", "s1"), ("a", "s1", "a", "s2"),
                                  ("a", "s1", "unknown", "s2")])
def test_invalid_pairs_rejected(tmp_path, pair):
    a = native(tmp_path / "a", {"one": ["a", "b"]}, [pair])
    b = native(tmp_path / "b", {"one": ["a", "b"]}, [])
    with pytest.raises(ValueError, match="Self-pair, same-species"):
        compare(a, b)


def test_species_ownership_changes_rejected(tmp_path):
    a = native(tmp_path / "a", {"one": ["a", "b"]}, [("a", "s1", "b", "s2")])
    b = native(tmp_path / "b", {"one": ["a", "b"]}, [("a", "s3", "b", "s2")])
    with pytest.raises(ValueError, match="ownership differs"):
        compare(a, b)


def test_gene_universe_changes_and_wrong_expected_size_rejected(tmp_path):
    a = native(tmp_path / "a", {"one": ["a", "b"]}, [])
    b = native(tmp_path / "b", {"one": ["a", "b"]}, [])
    with pytest.raises(ValueError, match="gene universes"):
        compare(a, b, 3)
    native_b = tmp_path / "b/orthohmm_root_hogs.tsv"
    native_b.write_text("root_hog\tsource_family\tgenes\none\tarbitrary\ta,c\n")
    with pytest.raises(ValueError, match="gene universes"):
        compare(a, b)


def test_duplicate_root_membership_rejected(tmp_path):
    a = native(tmp_path / "a", {"one": ["a", "b"], "two": ["b"]}, [])
    b = native(tmp_path / "b", {"one": ["a", "b"]}, [])
    with pytest.raises(ValueError, match="Duplicate gene membership"):
        compare(a, b)


def test_reconciliation_count_mismatch_rejected(tmp_path):
    a = native(tmp_path / "a", {"one": ["a", "b"]}, [])
    b = native(tmp_path / "b", {"one": ["a", "b"]}, [])
    (b / "reconciliation_summary.json").write_text(json.dumps({"root_hogs": 1, "ortholog_pairs": 10}))
    with pytest.raises(ValueError, match="reconciliation counts"):
        compare(a, b)


def test_input_change_during_readback_rejected(tmp_path, monkeypatch):
    a = native(tmp_path / "a", {"one": ["a", "b"]}, [("a", "s1", "b", "s2")])
    b = native(tmp_path / "b", {"one": ["a", "b"]}, [])
    original = diagnostic.pair_multiset

    def change_after_read(path, index, owners):
        result = original(path, index, owners)
        if path.parent == a:
            path.write_text("changed after parsing\n")
        return result

    monkeypatch.setattr(diagnostic, "pair_multiset", change_after_read)
    with pytest.raises(ValueError):
        compare(a, b)


def test_unexpected_pair_columns_rejected(tmp_path):
    a = native(tmp_path / "a", {"one": ["a", "b"]}, [])
    b = native(tmp_path / "b", {"one": ["a", "b"]}, [])
    (a / "orthohmm_pairwise_orthologs.tsv").write_text("wrong\theader\n")
    with pytest.raises(ValueError, match="pair header"):
        compare(a, b)


def test_partition_contrast_ignores_group_and_member_row_order(tmp_path):
    a, b = tmp_path / "a", tmp_path / "b"
    a.write_text("a b\nc\n")
    b.write_text("c\nb a\n")
    result = diagnostic.partition_contrast(a, b, {"a", "b", "c"})
    assert result["semantic_multisets_equal"]
    assert result["genes_in_changed_groups"] == 0
    b.write_text("a c\nb\n")
    result = diagnostic.partition_contrast(a, b, {"a", "b", "c"})
    assert not result["semantic_multisets_equal"]
    assert result["added_records"] == result["removed_records"] == 2
    assert result["genes_in_changed_groups"] == 3


@pytest.mark.parametrize("text", ["a b\n", "a b\nb c\n", "a b\nx\n", "a a b\nc\n"])
def test_partition_contrast_checks_full_unique_gene_universe(tmp_path, text):
    a, b = tmp_path / "a", tmp_path / "b"
    a.write_text("a b\nc\n")
    b.write_text(text)
    with pytest.raises(ValueError):
        diagnostic.partition_contrast(a, b, {"a", "b", "c"})
