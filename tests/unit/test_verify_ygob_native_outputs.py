from pathlib import Path

import pytest

from benchmark_tools.verify_ygob_native_outputs import CHECKPOINT_PATTERN, PYTHON, check_groups, expected_command, unique_path, verify_command


@pytest.mark.parametrize("method", ["high_sensitivity", "satellite_v2"])
def test_exact_native_and_harness_commands(method):
    root = Path("/root")
    args = expected_command(root, method)
    metrics = {"command": [PYTHON, str(root / "benchmarks/work/publication_method_native_v2/orthohmm/__main__.py"), *args],
               "harness": {"command": [PYTHON, "-m", "orthohmm", *args]}}
    verify_command(metrics, root, method)
    metrics["harness"]["command"].append("--unfrozen-option")
    with pytest.raises(ValueError, match="command differs"):
        verify_command(metrics, root, method)


def test_unknown_method_rejected():
    with pytest.raises(ValueError):
        expected_command(Path("/root"), "other")


def test_missing_members_are_reported_not_imputed():
    assert check_groups({"g": ["a", "b"]}, {"a", "b", "c"}) == {
        "groups": 1, "genes": 2, "missing_input_genes": 1, "singletons": 0}


@pytest.mark.parametrize("groups", [{"g": ["foreign"]}, {"a": ["a"], "b": ["a"]}])
def test_foreign_and_duplicate_genes_rejected(groups):
    with pytest.raises(ValueError):
        check_groups(groups, {"a"})


def test_ambiguous_or_missing_native_file_rejected(tmp_path):
    with pytest.raises(ValueError):
        unique_path(tmp_path, "*.txt")
    (tmp_path / "a.txt").touch()
    assert unique_path(tmp_path, "*.txt").name == "a.txt"
    (tmp_path / "b.txt").touch()
    with pytest.raises(ValueError):
        unique_path(tmp_path, "*.txt")


def test_checkpoint_selects_species_sequence_ids_not_global_numeric_ids(tmp_path):
    work = tmp_path / "WorkingDirectory"
    work.mkdir()
    (work / "clusters_OrthoFinder_I1.2.txt").touch()
    expected = work / "clusters_OrthoFinder_I1.2.txt_id_pairs.txt"
    expected.touch()
    assert unique_path(tmp_path, CHECKPOINT_PATTERN) == expected
