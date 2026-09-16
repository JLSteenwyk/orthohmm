from copy import deepcopy
from pathlib import Path

import pytest

from benchmark_tools.run_species_tree_control import control_cell


def inputs():
    return ({"label": "p1_c1_r1", "argv": ["python", "replay.py", "--species-tree-mode", "infer",
             "--output-directory", "/original", "--json", "/original.json", "--membership-constraints", "/constraints",
             "--root-rule", "species_overlap", "--pair-rule", "positive_paralogy", "--cpu", "32"]},
            {"label": "supplied_control", "rooted_rf_clade_distance": 0, "tree": {"path": "/control.nwk"}})


def test_only_mode_tree_cache_and_destinations_change():
    original, tree = inputs()
    unchanged = deepcopy(original)
    result = control_cell(original, tree, Path("/new"))
    assert original == unchanged
    args = result["argv"]
    assert args[args.index("--species-tree-mode") + 1] == "supplied"
    assert args[args.index("--checkpoint-source") + 1] == "/original"
    assert args[args.index("--species-tree") + 1] == "/control.nwk"
    assert args[args.index("--membership-constraints") + 1] == "/constraints"
    assert args[args.index("--cpu") + 1] == "32"
    assert args[args.index("--root-rule") + 1] == "species_overlap"
    assert args[args.index("--pair-rule") + 1] == "positive_paralogy"
    assert result["prediction"].startswith("/new/output/")


@pytest.mark.parametrize("change", ["cell", "tree", "distance", "mode", "checkpoint", "constraints"])
def test_unplanned_control_inputs_rejected(change):
    original, tree = inputs()
    if change == "cell":
        original["label"] = "p0_c0_r1"
    elif change == "tree":
        tree["label"] = "nni1_0"
    elif change == "distance":
        tree["rooted_rf_clade_distance"] = 2
    elif change == "mode":
        original["argv"][3] = "supplied"
    elif change == "checkpoint":
        original["argv"].extend(["--checkpoint-source", "/other"])
    else:
        original["argv"].extend(["--membership-constraints", "/other"])
    with pytest.raises(ValueError):
        control_cell(original, tree, Path("/new"))
