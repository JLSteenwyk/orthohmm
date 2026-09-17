import pytest

from benchmark_tools.read_wgd_native_groups import read_orthohmm, read_species_table

OWNERS = {"a": "Scerevisiae", "b": "Smikatae", "c": "Scerevisiae"}
COLUMNS = {"Scerevisiae": "Scerevisiae", "Smikatae": "Smikatae"}
OF_HEADER = "HOG\tOG\tGene Tree Parent Clade\tScerevisiae\tSmikatae\n"
SONIC_HEADER = "group_id\tgroup_size\tsp_in_grp\tseed_ortholog_cnt\tScerevisiae\tSmikatae\n"


def test_root_table_preserves_missing_genes_without_singleton_completion(tmp_path):
    path = tmp_path / "N0.tsv"
    path.write_text(OF_HEADER + "N0.HOG0001\tOG0\tn0\ta\tb\n")
    assert read_species_table(path, "orthofinder_root_hogs", OWNERS, COLUMNS) == {"N0.HOG0001": ["a", "b"]}


@pytest.mark.parametrize("payload", [
    "h\tOG0\tn0\tb\ta\n", "h\tOG0\tn0\ta\n", "h\tOG0\tn0\ta\tb\textra\n",
    "h\tOG0\tn0\ta, a\tb\n", "h\tOG0\tn0\ta,,c\tb\n", "h\tOG0\tn0\talien\tb\n",
    "h\tOG0\tn0\ta\tb\nh2\tOG0\tn1\tc\tb\n", "h\tOG0\tn0\t\t\n",
    "h\tOG0\tn0\ta\t\nh\tOG1\tn1\tc\tb\n",
])
def test_invalid_native_hog_tables_rejected(tmp_path, payload):
    path = tmp_path / "N0.tsv"
    path.write_text(OF_HEADER + payload)
    with pytest.raises(ValueError):
        read_species_table(path, "orthofinder_root_hogs", OWNERS, COLUMNS)


def test_sonic_species_and_counts_checked(tmp_path):
    path = tmp_path / "ortholog_groups.tsv"
    path.write_text(SONIC_HEADER + "g0\t2\t2\t2\ta\tb\ng1\t1\t1\t0\tc\t*\n")
    assert read_species_table(path, "sonicparanoid", OWNERS, COLUMNS) == {"g0": ["a", "b"], "g1": ["c"]}
    path.write_text(SONIC_HEADER + "g0\t3\t2\t2\ta\tb\n")
    with pytest.raises(ValueError, match="counts"):
        read_species_table(path, "sonicparanoid", OWNERS, COLUMNS)


def test_no_silent_species_header_normalization(tmp_path):
    path = tmp_path / "N0.tsv"
    path.write_text(OF_HEADER.replace("Smikatae", "Smikatae.fasta") + "h\tOG0\tn0\ta\tb\n")
    with pytest.raises(ValueError, match="columns"):
        read_species_table(path, "orthofinder_root_hogs", OWNERS, COLUMNS)
    mapping = {"Scerevisiae": "Scerevisiae", "Smikatae.fasta": "Smikatae"}
    assert read_species_table(path, "orthofinder_root_hogs", OWNERS, mapping)["h"] == ["a", "b"]


def test_orthohmm_native_group_and_root_formats(tmp_path):
    path = tmp_path / "groups.txt"
    path.write_text("OG0: a b\nOG1: c\n")
    assert len(read_orthohmm(path, "named_groups", OWNERS)) == 2
    path.write_text("root_hog\tsource_family\tgenes\nh0\tf0\ta,b\nh1\tf0\tc\n")
    assert read_orthohmm(path, "root_hogs", OWNERS)["h0"] == ["a", "b"]


def test_duplicate_species_header_rejected(tmp_path):
    path = tmp_path / "N0.tsv"
    path.write_text(OF_HEADER.replace("Smikatae", "Scerevisiae") + "h\tOG0\tn0\ta\tc\n")
    with pytest.raises(ValueError, match="columns"):
        read_species_table(path, "orthofinder_root_hogs", OWNERS, COLUMNS)
