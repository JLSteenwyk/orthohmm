from copy import deepcopy
import hashlib
import json
from pathlib import Path

import pytest

from benchmark_tools.prepare_scaling_inputs import METHODS, SIZES
from benchmark_tools.prepare_threadripper_overhead import ARMS, arm_order, build, main
from benchmark_tools.prepare_threadripper_run import paths


def parent():
    path = Path(__file__).resolve().parents[2] / "benchmark_tools/results/threadripper_private_commands_20260928.json"
    return json.loads(path.read_text())


def roots(tmp_path):
    return tmp_path / "outputs", Path("/dev/shm") / f"overhead-test-{tmp_path.name}"


def test_complete_paired_plan(tmp_path):
    original = parent()
    before = deepcopy(original)
    out, inp = roots(tmp_path)
    result = build(original, out, inp)
    assert original == before
    assert len(result["runs"]) == 54
    assert result["execution_authorized"] is False
    assert result["scientific_timings_admitted"] is False
    assert not out.exists() and not inp.exists()
    assert result["resources"] == original["resources"]
    assert [r["index"] for r in result["runs"]] == list(range(54))
    assert len({r["run"]["configuration"]["output"] for r in result["runs"]}) == 54
    assert len({r["run"]["prepared_input_directory"] for r in result["runs"]}) == 54
    for pair, old in enumerate(original["runs"]):
        selected = result["runs"][2 * pair:2 * pair + 2]
        assert {r["arm"] for r in selected} == set(ARMS)
        for task in selected:
            assert task["pair"] == pair
            run = task["run"]
            assert run["dataset"] == old["dataset"]
            assert run["cwd"] == old["cwd"]
            assert run["native_argv"][0] == old["native_argv"][0]
            assert run["original_native_argv"] == old["original_native_argv"]
            assert run["input_creation_order"] == old["input_creation_order"]
            paths(run)


@pytest.mark.parametrize("method", METHODS)
@pytest.mark.parametrize("size", SIZES)
def test_precommitted_orders_balance_first_two_repeats(method, size):
    assert arm_order(method, size, 0) == list(reversed(arm_order(method, size, 1)))
    assert arm_order(method, size, 2) == arm_order(method, size, 2)
    assert set(arm_order(method, size, 2)) == set(ARMS)


@pytest.mark.parametrize("change", ["order", "missing", "duplicate", "bool_index", "host", "cpus",
    "memory", "timeout", "authorized", "method", "wrong_native", "input_bytes", "input_order",
    "expected_order", "input_count", "old_output", "old_input", "existing_output", "existing_input",
    "disk_input", "tmpfs_output", "traversal"])
def test_reject_drift_or_reused_storage(tmp_path, change):
    source = parent()
    out, inp = roots(tmp_path)
    if change == "order":
        source["runs"].reverse()
    elif change == "missing":
        source["runs"].pop()
    elif change == "duplicate":
        source["runs"][1] = deepcopy(source["runs"][0])
    elif change == "bool_index":
        source["runs"][0]["index"] = False
    elif change in {"host", "cpus", "memory", "timeout"}:
        key = {"host": "host", "cpus": "native_workers", "memory": "memory_bytes", "timeout": "native_timeout_s"}[change]
        source["resources"][key] = "changed"
    elif change == "authorized":
        source["scientific_execution_authorized"] = True
    elif change == "method":
        source["runs"][0]["method"] = "different"
    elif change == "wrong_native":
        source["runs"][0]["native_method"] = "orthofinder_full"
    elif change == "input_bytes":
        source["runs"][9]["dataset"]["inputs"][0]["bytes"] += 1
    elif change == "input_order":
        source["runs"][9]["input_creation_order"].reverse()
    elif change == "expected_order":
        source["runs"][0]["expected_native_order"].reverse()
    elif change == "input_count":
        source["runs"][0]["dataset"]["inputs"].pop()
    elif change == "old_output":
        out = Path(source["runs"][0]["measurement_directory"]).parent.parent / "nested"
    elif change == "old_input":
        inp = Path(source["runs"][0]["prepared_input_directory"]).parent.parent / "nested"
    elif change == "existing_output":
        out.mkdir()
    elif change == "existing_input":
        inp = Path("/dev/shm")
    elif change == "disk_input":
        inp = tmp_path / "input"
    elif change == "tmpfs_output":
        out = inp / "output"
    else:
        out = tmp_path / "escape" / ".." / "output"
    with pytest.raises(ValueError):
        build(source, out, inp)


def test_reject_symlink_root(tmp_path):
    out, inp = roots(tmp_path)
    destination = tmp_path / "real"
    destination.mkdir()
    out.symlink_to(destination, target_is_directory=True)
    with pytest.raises(ValueError):
        build(parent(), out, inp)


def test_command_values_are_not_rewritten_by_substring(tmp_path):
    source = parent()
    text = "label-containing-" + source["runs"][0]["configuration"]["output"]
    source["runs"][0]["configuration"]["label"] = text
    result = build(source, *roots(tmp_path))
    assert result["runs"][0]["run"]["configuration"]["label"] == text


@pytest.mark.parametrize("failure", [None, "parent_pin", "protocol_pin", "existing", "symlink"])
def test_cli_pins_and_no_overwrite(tmp_path, monkeypatch, failure):
    source = Path(__file__).resolve().parents[2] / "benchmark_tools/results/threadripper_private_commands_20260928.json"
    protocol = tmp_path / "protocol.md"
    protocol.write_text("Synthetic prospective protocol.\n")
    out, inp = roots(tmp_path)
    manifest = tmp_path / "plan.json"
    if failure == "existing":
        manifest.write_text("preserve")
    elif failure == "symlink":
        manifest.symlink_to(tmp_path / "missing")
    parent_sha = hashlib.sha256(source.read_bytes()).hexdigest()
    protocol_sha = hashlib.sha256(protocol.read_bytes()).hexdigest()
    args = ["prepare", "--parent", str(source), "--parent-sha256", parent_sha,
        "--protocol", str(protocol), "--protocol-sha256", protocol_sha,
        "--output-root", str(out), "--input-root", str(inp), "--manifest", str(manifest)]
    if failure == "parent_pin":
        args[4] = "0" * 64
    elif failure == "protocol_pin":
        args[8] = "0" * 64
    monkeypatch.setattr("sys.argv", args)
    if failure is not None:
        with pytest.raises((ValueError, FileExistsError)):
            main()
        if failure == "existing":
            assert manifest.read_text() == "preserve"
        elif failure != "symlink":
            assert not manifest.exists()
    else:
        main()
        value = json.loads(manifest.read_text())
        assert len(value["runs"]) == 54
        assert value["sources"][0]["sha256"] == parent_sha
        assert value["sources"][1]["sha256"] == protocol_sha
        assert len(value["helpers"]) == 6
    assert not out.exists() and not inp.exists()
