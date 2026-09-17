import pytest

from benchmark_tools.admit_wgd_comparators import check_root_tree, check_sonic_log


@pytest.mark.parametrize("tree", ["((a,b)N1,c)N0;", "(a,(b,c)N1)N0;"])
def test_root_scope_accepts_different_inferred_topologies(tmp_path, tree):
    path = tmp_path / "tree.txt"
    path.write_text(tree)
    check_root_tree(path, ["a", "b", "c"])


@pytest.mark.parametrize("tree", ["((a,b)N1,c)N2;", "(a,b)N0;", "(a,a,c)N0;", "(a,b,d)N0;"])
def test_wrong_root_or_species_rejected(tmp_path, tree):
    path = tmp_path / "tree.txt"
    path.write_text(tree)
    with pytest.raises(ValueError, match="scope"):
        check_root_tree(path, ["a", "b", "c"])


def sonic_fixture():
    selected = {"copy_inputs_to": "/input", "argv": ["/sonic", "-i", "/input", "-o", "/output", "-t", "32"]}
    text = "\n".join(["SonicParanoid 2.0.9 will be executed with the following parameters:",
                      "Input directory: /input", "Main output directory: /output", "Input proteomes:\t4",
                      "Threads:\t32", "Alignment tool:\tdiamond", "Run mode:\tdefault (Diamond [--very-sensitive])",
                      "MCL inflation:\t1.50", "Perfom only graph-based orthology:\tFalse", "Total elapsed time (seconds):\t12.1"])
    return selected, text


def test_sonic_native_configuration_and_completion():
    selected, text = sonic_fixture()
    check_sonic_log(text, selected, 4)


@pytest.mark.parametrize("old,new", [("\t32", "\t16"), ("2.0.9", "2.0.8"),
                                     ("12.1", "nan"), ("Total elapsed time", "not finished"),
                                     ("\tFalse", "\tTrue"), ("/input", "/different")])
def test_sonic_native_configuration_or_completion_mismatch_rejected(old, new):
    selected, text = sonic_fixture()
    with pytest.raises(ValueError):
        check_sonic_log(text.replace(old, new), selected, 4)


def test_multiple_sonic_runs_in_same_log_rejected():
    selected, text = sonic_fixture()
    with pytest.raises(ValueError):
        check_sonic_log(text + "\n" + text, selected, 4)


def test_repeated_consistent_stage_output_directory_is_valid():
    selected, text = sonic_fixture()
    check_sonic_log(text + "\nMain output directory: /output", selected, 4)
    with pytest.raises(ValueError, match="directory"):
        check_sonic_log(text + "\nMain output directory: /other", selected, 4)
