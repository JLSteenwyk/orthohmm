import pytest

from benchmark_tools.audit_ob_matrix_provenance import SCHEMAS, native_groups, sonic_metadata


def write_table(path, method, rows, species=("a.fa", "b.fa")):
    path.write_text("\t".join((*SCHEMAS[method], *species)) + "\n" + "\n".join(rows) + "\n")


@pytest.mark.parametrize("method,row", [("sonicparanoid", "1\t2\t2\t2\ta\tb"), ("proteinortho", "2\t2\t0.3\ta\tb")])
def test_native_conversion_adds_only_missing_singletons(tmp_path, method, row):
    path = tmp_path / "native.tsv"
    write_table(path, method, [row])
    assert native_groups(path, method, dict(a="a.fa", b="b.fa", c="b.fa")) == ([("a", "b"), ("c",)], 1, 2, 1)


@pytest.mark.parametrize("row", ["1\t2\t2\t2\ta\tx", "1\t2\t2\t2\tb\ta",
                                 "1\t3\t2\t2\ta\tb", "1\t2\t1\t2\ta\tb",
                                 "1\t3\t2\t2\ta,a\tb", "1\t2\t2\t2\ta",
                                 "1\t2\t2\t2\ta\tb\textra"])
def test_bad_native_rows_rejected(tmp_path, row):
    path = tmp_path / "native.tsv"
    write_table(path, "sonicparanoid", [row])
    with pytest.raises(ValueError):
        native_groups(path, "sonicparanoid", dict(a="a.fa", b="b.fa"))


def test_species_inventory_and_duplicate_membership(tmp_path):
    path = tmp_path / "native.tsv"
    write_table(path, "proteinortho", ["2\t2\t0.3\ta\tb"], species=("a.fa", "c.fa"))
    with pytest.raises(ValueError, match="inventory"):
        native_groups(path, "proteinortho", dict(a="a.fa", b="b.fa"))
    write_table(path, "proteinortho", ["2\t2\t0.3\ta\tb", "2\t2\t0.3\ta\tb"])
    with pytest.raises(ValueError, match="Duplicate native membership"):
        native_groups(path, "proteinortho", dict(a="a.fa", b="b.fa"))


def test_metadata_preserves_settings_and_rejects_duplicates():
    text = "SonicParanoid 2.0.9\nInput proteomes:\t12\nThreads:\t32\n"
    assert sonic_metadata(text)["Threads:"] == "32"
    with pytest.raises(ValueError, match="Duplicate"):
        sonic_metadata(text + "Threads:\t16\n")
    with pytest.raises(ValueError, match="version/input"):
        sonic_metadata(text.replace("2.0.9", "2.0.8"))
