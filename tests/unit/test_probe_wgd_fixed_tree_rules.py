"""Focused guards for the fixed-tree diagnostic, not accuracy validation."""

import hashlib

import pytest

from benchmark_tools import probe_wgd_fixed_tree_rules as probe


def test_partition_is_order_independent():
    assert probe.partition([["b", "a"], ["c"]]) == [("a", "b"), ("c",)]


@pytest.mark.parametrize("groups", [[["a", "a"]], [["a"], ["a"]], [[]]])
def test_partition_rejects_invalid_membership(groups):
    with pytest.raises(ValueError):
        probe.partition(groups)


def test_biological_row_ignores_only_display_ids():
    original = {"anchor_groups": ["display"], "coverage_numerator": 3,
                "anchor_group_sizes": [1, 4], "supported_separation_rate": 0}
    result = probe.biological_row(original)
    assert result == {k: v for k, v in original.items() if k != "anchor_groups"}
    assert "anchor_groups" in original
    assert probe.biological_row({**original, "coverage_numerator": 4}) != result


def test_frozen_module_rejects_changed_source_before_execution(tmp_path):
    path = tmp_path / "phylogeny.py"
    path.write_text("raise AssertionError('must not execute')\n")
    with pytest.raises(ValueError, match="source changed"):
        probe.frozen_module(path)


def test_frozen_module_supports_dataclass_loading(tmp_path, monkeypatch):
    path = tmp_path / "phylogeny.py"
    data = b"from dataclasses import dataclass\n@dataclass\nclass Value:\n    count: int\n"
    path.write_bytes(data)
    monkeypatch.setattr(probe, "SOURCE_SHA", hashlib.sha256(data).hexdigest())
    assert probe.frozen_module(path).Value(3).count == 3


def test_run_rejects_changed_protocol_before_loading_inputs(tmp_path):
    results = tmp_path / "benchmark_tools/results"
    results.mkdir(parents=True)
    (results / "WGD_FIXED_TREE_RULE_PROTOCOL_20260928.md").write_text("changed")
    with pytest.raises(ValueError, match="protocol changed"):
        probe.run(tmp_path)
