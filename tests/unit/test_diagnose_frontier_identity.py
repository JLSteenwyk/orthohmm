import copy

import pytest

from benchmark_tools.diagnose_frontier_identity import changes, diagnose


def point(extra=None):
    ids = {"/job": [1, 2], **(extra or {})}
    scopes = sorted(ids)
    counter = lambda scope, n: dict(scope=scope, started_ns=n, finished_ns=n,
                                    raw="usage_usec 0\n")
    return dict(boot_id="boot", target="/job",
                inventory_before=dict(identities=ids),
                inventory_after=dict(identities=copy.deepcopy(ids)),
                root=[counter("/", 0), counter("/", len(scopes)+1)],
                rows=[counter(scope, i+1) for i, scope in enumerate(scopes)])


def test_stable():
    assert not any(changes(point(), point()).values())


def test_add_remove_replace():
    left = point({"/old": [1, 3], "/service": [1, 4]})
    right = point({"/new": [1, 5], "/service": [1, 6]})
    result = changes(left, right)
    assert result["added"] == {"/new": [1, 5]}
    assert result["removed"] == {"/old": [1, 3]}
    assert result["replaced"] == {"/service": dict(before=[1, 4], after=[1, 6])}
    assert not result["boot_changed"] and not result["target_changed"]


def test_boot_change():
    right = point()
    right["boot_id"] = "other"
    assert changes(point(), right)["boot_changed"]


def test_invalid_snapshot_rejected():
    right = point()
    right["inventory_after"]["identities"]["/job"] = [1, 9]
    with pytest.raises(ValueError, match="identity changed"):
        changes(point(), right)


def test_terminal_gate_precedes_native_reads(tmp_path):
    (tmp_path / "accounting.txt").write_text("21889_0|RUNNING|0:0|00:01:00|20|96G|spark-7ff0\n")
    with pytest.raises(ValueError, match="nonterminal"):
        diagnose(tmp_path)
