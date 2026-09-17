from copy import deepcopy
from pathlib import Path

import pytest

from benchmark_tools.prepare_native_scaling_run import check_copies, prepare, run_directory
from benchmark_tools.snapshot_orthohmm_input_order import record


@pytest.fixture
def fixture(tmp_path):
    inputs = tmp_path / "input"
    inputs.mkdir()
    for name in ("z.fa", "a.fa"):
        (inputs / name).write_text(">" + name + "\nMALW\n")
    records = [record(inputs / name) for name in ("z.fa", "a.fa")]
    root = tmp_path / "run"
    root.mkdir()
    run = {"measurement_directory": str(root / "measurement"),
           "native_method": "orthohmm_high_sensitivity", "native_argv": ["/python", "-m", "orthohmm"],
           "configuration": {"output": str(root / "native"), "metrics": str(root / "metrics.json")},
           "dataset": {"input_directory": str(inputs), "proteomes": 2, "inputs": records}}
    order = {"input_directory": str(inputs), "proteomes": 2, "native_order": ["z.fa", "a.fa"],
             "inputs_in_native_order": records}
    return run, order


def test_hmm_empty_directory_and_exact_wrapper(fixture):
    run, order = fixture
    original = deepcopy(run)
    result = prepare(run, order)
    assert run == original
    assert Path(run["configuration"]["output"]).is_dir()
    assert result["measured_argv"][-3:] == run["native_argv"]
    assert result["expected_native_basename_order"] == ["z.fa", "a.fa"]
    assert not Path(run["measurement_directory"]).exists()
    with pytest.raises(ValueError, match="fresh empty"):
        prepare(run, order)


def orthofinder(run):
    run["native_method"] = "orthofinder_full"
    run["configuration"].pop("metrics")
    run["configuration"].update(copy_inputs_from=run["dataset"]["input_directory"],
                                copy_inputs_to=run["configuration"]["output"] + "/input")


def test_orthofinder_copies_preserve_bytes_and_names(fixture):
    run, order = fixture
    orthofinder(run)
    result = prepare(run, order)
    assert result["expected_native_basename_order"] == ["a.fa", "z.fa"]
    assert len(result["copied_inputs"]) == 2
    assert check_copies(run) == result["copied_inputs"]
    Path(run["configuration"]["copy_inputs_to"], "z.fa").write_text("changed")
    with pytest.raises(ValueError, match="Copied native"):
        check_copies(run)


def test_changed_original_never_prepares_output(fixture):
    run, order = fixture
    Path(run["dataset"]["inputs"][0]["path"]).write_text("changed")
    with pytest.raises(ValueError, match="Original input"):
        prepare(run, order)
    assert not Path(run["configuration"]["output"]).exists()


@pytest.mark.parametrize("key", ["output", "metrics", "copy_inputs_to"])
def test_paths_cannot_escape(fixture, key):
    run, order = fixture
    run["configuration"][key] = "/tmp/unrelated"
    with pytest.raises(ValueError, match="escapes"):
        prepare(run, order)


def test_reject_collector_overlap_and_symlinks(fixture):
    run, order = fixture
    run["configuration"]["output"] = run["measurement_directory"] + "/native"
    with pytest.raises(ValueError, match="overlaps"):
        run_directory(run)
    root = Path(run["measurement_directory"]).parent
    link = root.parent / "link"
    link.symlink_to(root, target_is_directory=True)
    run["measurement_directory"] = str(link / "measurement")
    with pytest.raises(ValueError, match="symlinks"):
        run_directory(run)


def test_wrong_order_identity_rejected(fixture):
    run, order = fixture
    order["native_order"].reverse()
    with pytest.raises(ValueError, match="identities"):
        prepare(run, order)
