import json

import pytest

from benchmark_tools.run_qfo_canonical_phylogeny import record, validate


def fixture(tmp_path):
    input_path = tmp_path / "input"
    input_path.write_text("bound\n")
    plan = dict(arm="canonical_order", attempts=1, reuse_policy="admitted_input_identical_raw_trees",
        species_tree_cache_reuse=False, accuracy_evaluated=False, checked_records=[record(input_path)])
    path = tmp_path / "plan.json"
    path.write_text(json.dumps(plan))
    return path, input_path, plan


def test_canonical_scope_accepts_bound_plan(tmp_path):
    path, _, plan = fixture(tmp_path)
    assert validate(path, record(path)["sha256"]) == plan


@pytest.mark.parametrize("key,value", [("arm","retained_order"), ("attempts",2),
    ("reuse_policy","all_checkpoints"), ("species_tree_cache_reuse",True), ("accuracy_evaluated",True)])
def test_canonical_scope_rejects_mutation(tmp_path, key, value):
    path, _, plan = fixture(tmp_path)
    plan[key] = value
    path.write_text(json.dumps(plan))
    with pytest.raises(ValueError):
        validate(path, record(path)["sha256"])


def test_canonical_input_mutation(tmp_path):
    path, input_path, _ = fixture(tmp_path)
    input_path.write_text("altered")
    with pytest.raises(ValueError):
        validate(path, record(path)["sha256"])


def test_canonical_plan_hash_mutation(tmp_path):
    path, _, _ = fixture(tmp_path)
    with pytest.raises(ValueError):
        validate(path, "wrong")
