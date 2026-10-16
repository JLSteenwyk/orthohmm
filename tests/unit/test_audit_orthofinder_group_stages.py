import csv
import pytest

from benchmark_tools.audit_orthofinder_group_stages import audit, padded_hogs


def fixture(tmp_path):
    working = tmp_path / "WorkingDirectory"
    working.mkdir()
    (tmp_path / "Orthogroups").mkdir()
    (working / "SequenceIDs.txt").write_text("0_0: a:pep\n0_1: b\n1_0: c\n1_1: d\n")
    (working / "SpeciesIDs.txt").write_text("0: first.fa\n1: second.fa\n")
    (working / "clusters_OrthoFinder_I1.2.txt_id_pairs.txt").write_text(
        "(mclmatrix\nbegin\n0 0_0 0_1 1_0 1_1 $\n)\n")
    (tmp_path / "Orthogroups/Orthogroups.txt").write_text("OG0: a_pep b c\nOG1: d\n")
    with (working / "N0.ids.tsv").open("w", newline="") as stream:
        writer = csv.writer(stream, delimiter="\t")
        writer.writerow(["HOG", "OG", "Gene Tree Parent Clade", "first", "second"])
        writer.writerow(["N0.HOG0", "OG0", "n0", "0_0", "1_0"])
    return tmp_path


def test_official_groups_can_differ_from_later_hogs_and_from_mcl(tmp_path):
    root = fixture(tmp_path)
    result = audit(root)
    assert (result["genes"], result["species"], result["official_groups"], result["mcl_groups"]) == (4, 2, 2, 1)
    assert result["latest_root_hog_padded_groups"] == 3
    assert result["latest_root_hog_singleton_padding"] == 2
    assert result["identifiers_restored"] == 1
    assert result["official_differs_from_mcl"]
    assert not result["official_equals_latest_root_hog_partition"]
    assert result["latest_root_hogs_refine_official_partition"]
    assert result["genes_in_changed_groups"] == 3


def test_equal_partitions_do_not_require_inventing_a_difference(tmp_path):
    root = fixture(tmp_path)
    (root / "Orthogroups/Orthogroups.txt").write_text("OG0: a_pep c\nOG1: b\nOG2: d\n")
    assert audit(root)["official_equals_latest_root_hog_partition"]


def test_matching_mcl_does_not_prove_absence_of_phylogenetic_processing(tmp_path):
    root = fixture(tmp_path)
    (root / "WorkingDirectory/clusters_OrthoFinder_I1.2.txt_id_pairs.txt").write_text(
        "(mclmatrix\nbegin\n0 0_0 0_1 1_0 $\n1 1_1 $\n)\n")
    assert not audit(root)["official_differs_from_mcl"]


@pytest.mark.parametrize("row", [
    ["N0.HOG1", "OG0", "n0", "0_0", "0_1"],
    ["N0.HOG1", "OG0", "n0", "0_0, 0_0", "1_0"],
    ["N0.HOG1", "OG0", "n0", "missing", "1_0"],
    ["N0.HOG1", "OG0", "n0", "", ""],
    ["N0.HOG1", "OG0", "n0", "0_0"],
])
def test_invalid_hog_membership_or_width_rejected(tmp_path, row):
    root = fixture(tmp_path)
    with (root / "WorkingDirectory/N0.ids.tsv").open("a", newline="") as stream:
        csv.writer(stream, delimiter="\t").writerow(row)
    with pytest.raises(ValueError):
        audit(root)


def test_duplicate_hog_identifier_is_rejected(tmp_path):
    root = fixture(tmp_path)
    with (root / "WorkingDirectory/N0.ids.tsv").open("a", newline="") as stream:
        csv.writer(stream, delimiter="\t").writerow(["N0.HOG0", "OG0", "n0", "0_1", "1_1"])
    with pytest.raises(ValueError, match="duplicate root-HOG"):
        audit(root)


def test_cross_parent_later_hog_is_not_accepted_as_a_refinement(tmp_path):
    root = fixture(tmp_path)
    (root / "Orthogroups/Orthogroups.txt").write_text("OG0: a_pep b\nOG1: c d\n")
    with pytest.raises(ValueError, match="crosses official"):
        audit(root)


@pytest.mark.parametrize("text,match", [
    ("OG0: a_pep b c\nOG1: b d\n", "exact gene universe"),
    ("OG0: a_pep b c\n", "exact gene universe"),
    ("OG0: a_pep b c\nOG1: unknown\n", "Unknown official"),
])
def test_official_universe_must_be_exact(tmp_path, text, match):
    root = fixture(tmp_path)
    (root / "Orthogroups/Orthogroups.txt").write_text(text)
    with pytest.raises(ValueError, match=match):
        audit(root)


def test_ambiguous_sanitization_is_rejected(tmp_path):
    root = fixture(tmp_path)
    (root / "WorkingDirectory/SequenceIDs.txt").write_text("0_0: a:pep\n0_1: a_pep\n1_0: c\n1_1: d\n")
    with pytest.raises(ValueError, match="Ambiguous identifier"):
        audit(root)


def test_duplicate_original_identifiers_are_rejected(tmp_path):
    root = fixture(tmp_path)
    (root / "WorkingDirectory/SequenceIDs.txt").write_text("0_0: a\n0_1: a\n1_0: c\n1_1: d\n")
    with pytest.raises(ValueError, match="original identifiers"):
        audit(root)


def test_wrong_species_header_is_rejected(tmp_path):
    root = fixture(tmp_path)
    path = root / "WorkingDirectory/N0.ids.tsv"
    path.write_text("HOG\tOG\tGene Tree Parent Clade\twrong\tsecond\n")
    with pytest.raises(ValueError, match="columns"):
        padded_hogs(path, {"0_0": "a", "1_0": "b"}, {0: "first.fa", 1: "second.fa"})
