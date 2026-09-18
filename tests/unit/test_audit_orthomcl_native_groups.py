import json

import pytest

from benchmark_tools import audit_orthomcl_native_groups as module


def fixture(tmp_path):
    files = [tmp_path / name for name in ("groups", "mcl", "index", "gg")]
    texts = [
        "ORTHOMCL0(3 genes,2 taxa): A(x) B(x) C(y)\n"
        "ORTHOMCL2(2 genes,1 taxa): E(y) F(y)\n",
        "(mclheader\nmcltype matrix\ndimensions 6x3\n)\n(mclmatrix\nbegin\n"
        "0 0 1\n 2 $\n1 3 $\n2 4 5 $\n)\n",
        "0\tA\n1\tB\n2\tC\n3\tD\n4\tE\n5\tF\n",
        "x: A B D\ny: C E F G\n",
    ]
    for path, text in zip(files, texts):
        path.write_text(text)
    return files


def test_partition_singletons_and_unindexed_proteins(tmp_path):
    assert module.validate(*fixture(tmp_path)) == {
        "input_proteins": 7, "input_species": 2, "indexed_proteins": 6,
        "mcl_clusters": 3, "mcl_singleton_clusters": 1, "final_groups": 2,
        "grouped_proteins": 5, "ungrouped_input_proteins": 2,
        "input_proteins_absent_from_index": 1, "single_species_final_groups": 1,
        "cross_species_clique_pairs": 2,
    }


@pytest.mark.parametrize("file,old,new", [
    (0, "3 genes", "2 genes"), (0, "2 taxa", "1 taxa"),
    (0, "A(x)", "A(y)"), (0, "A(x)", "B(x)"),
    (0, "A(x)", "D(x)"), (0, "A(x)", "G(y)"),
    (0, "A(x)", "Z(x)"), (0, "A(x)", "A"),
    (0, "ORTHOMCL2", "ORTHOMCL0"), (0, "ORTHOMCL2", "ORTHOMCL1"),
    (0, "ORTHOMCL2", "ORTHOMCL3"), (0, "ORTHOMCL0", "ORTHOMCL00"),
    (0, "ORTHOMCL2(2 genes,1 taxa): E(y) F(y)\n", ""),
    (1, "6x3", "5x3"), (1, "6x3", "6x4"), (1, "6x3", "6x0"),
    (1, "6x3", "6x7"), (1, "mcltype matrix", "mcltype graph"),
    (1, "begin", ""), (1, "2 4 5 $", "2 4 5"),
    (1, "2 4 5 $", "2 4 $"), (1, "2 4 5 $", "2 4 6 $"),
    (1, "2 4 5 $", "2 4 4 $"), (1, "2 4 5 $", "2 4 0 $"),
    (1, "2 4 5 $", "3 4 5 $"), (1, "2 4 5 $", "2 4 5:1 $"),
    (1, "2 4 5 $", "2 4 -1 $"), (1, "2 4 5 $", "2 4 +5 $"),
    (1, "2 4 5 $", "2 $"), (1, "2 4 5 $", "2 4 5 $ junk"),
    (2, "1\tB", "2\tB"), (2, "1\tB", "1\tA"),
    (2, "1\tB", "1\tZ"), (2, "1\tB", "1 B"),
    (2, "0\tA", "00\tA"), (3, "x: A B D", "x: A A D"),
])
def test_reject_corruption(tmp_path, file, old, new):
    files = fixture(tmp_path)
    files[file].write_text(files[file].read_text().replace(old, new))
    with pytest.raises(ValueError):
        module.validate(*files)


@pytest.mark.parametrize("file", range(4))
def test_empty_input_rejected(tmp_path, file):
    files = fixture(tmp_path)
    files[file].write_text("")
    with pytest.raises(ValueError):
        module.validate(*files)


def test_order_is_not_membership(tmp_path):
    files = fixture(tmp_path)
    lines = files[0].read_text().replace("A(x) B(x) C(y)", "C(y) B(x) A(x)").splitlines()
    files[0].write_text("\n".join(reversed(lines)) + "\n")
    assert module.validate(*files)["cross_species_clique_pairs"] == 2


@pytest.mark.parametrize("suffix", ["\n0 0 $\n", "\n(mclmatrix\nbegin\n)\n"])
def test_partition_trailing_content_rejected(tmp_path, suffix):
    files = fixture(tmp_path)
    files[1].write_text(files[1].read_text() + suffix)
    with pytest.raises(ValueError, match="Trailing"):
        module.validate(*files)


def test_unterminated_partition_rejected(tmp_path):
    files = fixture(tmp_path)
    files[1].write_text(files[1].read_text()[:-2])
    with pytest.raises(ValueError, match="Incomplete"):
        module.validate(*files)


def test_singleton_must_not_be_added(tmp_path):
    files = fixture(tmp_path)
    files[0].write_text(files[0].read_text() + "ORTHOMCL1(1 genes,1 taxa): D(x)\n")
    with pytest.raises(ValueError, match="singleton"):
        module.validate(*files)


def test_all_singletons_allow_empty_final_output(tmp_path):
    files = fixture(tmp_path)
    files[0].write_text("")
    files[1].write_text("(mclheader\nmcltype matrix\ndimensions 6x6\n)\n(mclmatrix\nbegin\n"
                        + "".join(f"{i} {i} $\n" for i in range(6)) + ")\n")
    result = module.validate(*files)
    assert result["final_groups"] == result["grouped_proteins"] == result["cross_species_clique_pairs"] == 0
    assert result["ungrouped_input_proteins"] == 7


def test_raw_identifiers_are_not_normalized(tmp_path):
    files = fixture(tmp_path)
    for path in files:
        path.write_text(path.read_text().replace("A", "sp|ACC|gene"))
    assert module.validate(*files)["grouped_proteins"] == 5


def test_report_is_provenance_bound_not_accuracy_admission(tmp_path):
    output = tmp_path / "report.json"
    args = fixture(tmp_path)
    result = module.audit(*args, output)
    assert json.loads(output.read_text()) == result
    assert result["status"] == "native_final_groups_match_mcl_partition"
    assert result["accuracy_admitted"] is result["publication_ready"] is False
    assert len(result["checked_records"]) == 7
    for item in result["checked_records"]:
        module.check(item)
    with pytest.raises(FileExistsError):
        module.audit(*args, output)


def test_no_report_on_invalid_data(tmp_path):
    args = fixture(tmp_path)
    args[0].write_text("")
    output = tmp_path / "report.json"
    with pytest.raises(ValueError):
        module.audit(*args, output)
    assert not output.exists()


def test_dangling_output_symlink_rejected(tmp_path):
    args = fixture(tmp_path)
    output = tmp_path / "report.json"
    output.symlink_to(tmp_path / "missing")
    with pytest.raises(FileExistsError):
        module.audit(*args, output)
