from pathlib import Path

import pytest

from benchmark_tools.run_ob_candidate_neighborhood import VARIANTS, variant_cell


def baseline():
    return {"label": "p1_c1_r1", "argv": ["python", "replay.py", "--candidate-clusters", "old_groups",
        "--membership-constraints", "old_constraints", "--output-directory", "baseline_output",
        "--json", "baseline_metrics", "--species-tree-mode", "infer", "--cpu", "32",
        "--pair-rule", "positive_paralogy"]}


@pytest.mark.parametrize("label", VARIANTS)
def test_changes_only_inputs_outputs_and_checkpoint(label):
    original = baseline()
    arm = {"label": label, "partition": {"path": "new_groups"}, "constraints": {"path": "new_constraints"}}
    result = variant_cell(original, arm, Path("/out"))
    expected = baseline()["argv"]
    for old, new in (("old_groups", "new_groups"), ("old_constraints", "new_constraints"),
                     ("baseline_output", "/out/output"), ("baseline_metrics", "/out/metrics.json")):
        expected[expected.index(old)] = new
    assert result["argv"] == expected + ["--checkpoint-source", "baseline_output"]
    assert original == baseline()
    assert result["candidate_partition"] == "new_groups"


@pytest.mark.parametrize("change", ["baseline", "arm", "supplied", "tree", "checkpoint", "duplicate", "missing"])
def test_rejects_changed_design(change):
    original = baseline()
    arm = {"label": "norm_low", "partition": {"path": "groups"}, "constraints": {"path": "constraints"}}
    if change == "baseline":
        original["label"] = "p0_c1_r1"
    elif change == "arm":
        arm["label"] = "control"
    elif change == "supplied":
        original["argv"][original["argv"].index("infer")] = "supplied"
    elif change in ("tree", "checkpoint"):
        original["argv"] += ["--species-tree" if change == "tree" else "--checkpoint-source", "x"]
    elif change == "duplicate":
        original["argv"] += ["--membership-constraints", "x"]
    else:
        original["argv"].remove("--membership-constraints")
    with pytest.raises(ValueError):
        variant_cell(original, arm, Path("/out"))
