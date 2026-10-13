import copy

import pytest

from benchmark_tools.reconstruct_component_recovery import compare_reconstruction


def fixture():
    old = {"ownership": {"classes": {"0": 4}}, "component_counts": {"families": 1}, "structural_controls": {}}
    for label in ("control", "components500"):
        old["structural_controls"][label] = {"seed": {"sha256": "seed"}, "copy_split_trace": {"sha256": "trace"},
            "candidate": {key: {"sha256": key} for key in ("seed_partition", "candidate_partition", "membership_constraints")}}
    rows = copy.deepcopy(old["structural_controls"])
    diagnostics = {"ownership": {"classes": {0: 4}}, "component_counts": {"families": 1}}
    return old, rows, diagnostics


def test_reconstruction_verifies_every_artifact_and_canonicalizes_json_ownership_only():
    old, rows, diagnostics = fixture()
    before = copy.deepcopy((old, rows, diagnostics))
    calls = []
    checked = compare_reconstruction(old, rows, diagnostics, calls.append)
    assert len(checked) == len(calls) == 20
    assert (old, rows, diagnostics) == before


@pytest.mark.parametrize("key", ["seed", "copy_split_trace", "seed_partition", "candidate_partition", "membership_constraints"])
@pytest.mark.parametrize("label", ["control", "components500"])
def test_any_changed_seed_trace_candidate_or_constraint_is_rejected(label, key):
    old, rows, diagnostics = fixture()
    target = rows[label] if key in ("seed", "copy_split_trace") else rows[label]["candidate"]
    target[key]["sha256"] = "wrong"
    with pytest.raises(ValueError, match="changed"):
        compare_reconstruction(old, rows, diagnostics, lambda pin: None)


@pytest.mark.parametrize("mutation", ["inventory", "ownership", "counts"])
def test_arm_inventory_ownership_and_activation_mismatches_are_rejected(mutation):
    old, rows, diagnostics = fixture()
    if mutation == "inventory":
        del rows["control"]
    elif mutation == "ownership":
        diagnostics["ownership"]["classes"][0] = 5
    else:
        diagnostics["component_counts"]["families"] = 2
    with pytest.raises(ValueError):
        compare_reconstruction(old, rows, diagnostics, lambda pin: None)


def test_invalid_artifact_pin_is_not_accepted_just_because_hashes_agree():
    old, rows, diagnostics = fixture()
    def check(pin):
        raise ValueError("Changed file")
    with pytest.raises(ValueError, match="Changed file"):
        compare_reconstruction(old, rows, diagnostics, check)
