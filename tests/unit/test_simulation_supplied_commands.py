from copy import deepcopy

import pytest

from benchmark_tools.prepare_simulation_methods import commands
from benchmark_tools.simulation_supplied_commands import fresh_supplied_method


@pytest.fixture
def setup(tmp_path):
    tree = tmp_path / "tree.nwk"
    tree.write_text("((a,b),(c,d));\n")
    baseline = commands({"input": str(tmp_path / "inputs")}, tmp_path / "old",
                        tmp_path / "frozen", tmp_path / "python", tmp_path / "orthofinder")
    return baseline, tree, tmp_path / "new"


@pytest.mark.parametrize("method", ["orthohmm_satellite_v2", "orthofinder_full"])
def test_only_tree_and_output_locations_change(setup, method):
    baseline, tree, output = setup
    original = deepcopy(baseline[method])
    result = fresh_supplied_method(method, original, tree, output)
    restored = list(result["argv"][:-2])
    assert result["argv"][-2:] == (["-s", str(tree)] if method == "orthofinder_full"
                                    else ["--species-tree", str(tree)])
    if method == "orthofinder_full":
        restored[restored.index("-f") + 1] = original["copy_inputs_to"]
        assert result["copy_inputs_from"] == original["copy_inputs_from"]
    else:
        restored[3:5] = original["argv"][3:5]
        restored[restored.index("--species-tree-mode") + 1] = "infer"
    assert restored == original["argv"]
    assert original == baseline[method]
    assert not output.exists()


@pytest.mark.parametrize("method,extra", [
    ("orthohmm_satellite_v2", ["--species-tree=x"]),
    ("orthohmm_satellite_v2", ["--checkpoint-source", "old"]),
    ("orthohmm_satellite_v2", ["--species-tree-mode", "infer"]),
    ("orthofinder_full", ["-ft", "old"]),
    ("orthofinder_full", ["-s", "tree"]),
    ("orthofinder_full", ["-f", "other"]),
])
def test_reject_changed_baseline(setup, method, extra):
    baseline, tree, output = setup
    baseline[method]["argv"].extend(extra)
    with pytest.raises(ValueError):
        fresh_supplied_method(method, baseline[method], tree, output)


@pytest.mark.parametrize("location", ["original", "parent", "child", "input", "tree_parent", "existing"])
def test_preserve_existing_artifacts(setup, location):
    baseline, tree, output = setup
    method = baseline["orthohmm_satellite_v2"]
    old = output.parent / "old/orthohmm_satellite_v2"
    destinations = {"original": old, "parent": old.parent, "child": old / "child",
                    "input": output.parent / "inputs", "tree_parent": tree.parent, "existing": output}
    if location == "existing":
        output.mkdir()
    with pytest.raises((ValueError, FileExistsError)):
        fresh_supplied_method("orthohmm_satellite_v2", method, tree, destinations[location])


def test_missing_tree(setup):
    baseline, tree, output = setup
    tree.unlink()
    with pytest.raises(FileNotFoundError):
        fresh_supplied_method("orthofinder_full", baseline["orthofinder_full"], tree, output)
