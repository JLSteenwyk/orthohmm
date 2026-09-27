import json

import pytest

from benchmark_tools.run_qfo_fresh_phylogeny import config_arguments, record, validate


def make_plan(tmp_path):
    payload = tmp_path / "input"
    payload.write_text("bound\n")
    plan = dict(arm="retained_order", attempts=1, checkpoint_reuse=False,
                accuracy_evaluated=False, checked_records=[record(payload)])
    path = tmp_path / "plan.json"
    path.write_text(json.dumps(plan))
    return path, plan, payload


def test_bound_plan(tmp_path):
    path, plan, _ = make_plan(tmp_path)
    assert validate(path, record(path)["sha256"]) == plan


def test_changed_hash(tmp_path):
    path, _, _ = make_plan(tmp_path)
    with pytest.raises(ValueError, match="Changed phylogeny plan"):
        validate(path, "wrong")


def test_changed_input(tmp_path):
    path, _, payload = make_plan(tmp_path)
    payload.write_text("altered\n")
    with pytest.raises(ValueError, match="Changed pinned artifact"):
        validate(path, record(path)["sha256"])


@pytest.mark.parametrize("key,value", [("arm", "canonical_order"), ("attempts", 2),
                                      ("checkpoint_reuse", True), ("accuracy_evaluated", True)])
def test_scope_rejected(tmp_path, key, value):
    path, plan, _ = make_plan(tmp_path)
    plan[key] = value
    path.write_text(json.dumps(plan))
    with pytest.raises(ValueError, match="Wrong fresh-arm scope"):
        validate(path, record(path)["sha256"])


def test_exact_frozen_config():
    assert config_arguments(dict(aligner="/mafft", tree_builder="/FastTree")) == dict(
        mode="reconcile", species_tree_mode="infer", species_tree=None,
        aligner="/mafft", tree_builder="/FastTree", root_duplication_rule="species_overlap",
        pair_orthology_rule="positive_paralogy", species_tree_rooting="min_variance")
