from pathlib import Path
import tempfile
import uuid

import pytest

from benchmark_tools import prepare_threadripper_run as module


@pytest.fixture
def setup(tmp_path, monkeypatch):
    if not Path("/dev/shm").is_dir():
        pytest.skip("Native /dev/shm input-copy fixture required")
    tmp_path = tmp_path.resolve()
    raw = tmp_path / "raw"
    raw.mkdir()
    (raw / "a.fa").write_text(">a\nACDEFG\n")
    root = tmp_path / "run_00"
    root.mkdir()
    with tempfile.TemporaryDirectory(prefix="orthohmm_prepare_test_", dir="/dev/shm") as temporary:
        target = Path(temporary) / "run_00/input"
        run = dict(proteomes=1, native_method="orthohmm_high_sensitivity",
            dataset=dict(proteomes=1, input_directory=str(raw), inputs=[module.record(raw / "a.fa")]),
            measurement_directory=str(root / "measurement"), prepared_input_directory=str(target),
            input_creation_order=["a.fa"], expected_native_order=["a.fa"],
            configuration=dict(output=str(root / "output"), copy_inputs_from=str(raw), copy_inputs_to=str(target)),
            native_argv=["/python", "-m", "orthohmm", str(target), "-o", str(root / "output")])
        monkeypatch.setattr(module, "snapshot", lambda *args: {"datasets": [{"native_order": ["a.fa"]}]})
        yield run, dict(core_root="unused", core_sources=[]), root, target


@pytest.mark.parametrize("method", ["orthohmm_high_sensitivity", "orthohmm_satellite_v2", "orthofinder_full"])
def test_all_methods_receive_fresh_copies(setup, method):
    run, baseline, root, target = setup
    run["native_method"] = method
    if method == "orthofinder_full":
        run["native_argv"] = ["/orthofinder", "-f", str(target), "-o", str(root / "output")]
    result = module.prepare(run, baseline)
    assert result["inference_started"] is False
    assert result["input_bytes"] == (target / "a.fa").stat().st_size
    assert (root / "preparation.json").is_file()
    assert (root / "output").exists() == (method != "orthofinder_full")
    module.check_prepared(run, baseline)
    with pytest.raises(ValueError):
        module.prepare(run, baseline)


@pytest.mark.parametrize("change", ["bytes", "extra", "symlink"])
def test_after_run_recheck_rejects_changes(setup, change):
    run, baseline, root, target = setup
    module.prepare(run, baseline)
    if change == "bytes":
        (target / "a.fa").write_text("changed")
    elif change == "extra":
        (target / "extra").mkdir()
    else:
        (target / "a.fa").unlink()
        (target / "a.fa").symlink_to(Path(run["dataset"]["input_directory"]) / "a.fa")
    with pytest.raises(ValueError):
        module.check_prepared(run, baseline)


def test_preserves_failed_copy_attempt(setup, monkeypatch):
    run, baseline, root, target = setup
    monkeypatch.setattr(module, "snapshot", lambda *args: {"datasets": [{"native_order": ["wrong.fa"]}]})
    with pytest.raises(ValueError, match="enumeration"):
        module.prepare(run, baseline)
    assert (target / "a.fa").is_file()
    assert (root / "preparation_failed.json").is_file()


def test_reject_wrong_native_output(setup):
    run, baseline, root, target = setup
    run["native_argv"][-1] = "/tmp/unexpected"
    with pytest.raises(ValueError):
        module.prepare(run, baseline)
    assert not target.exists()


@pytest.fixture
def layout(tmp_path):
    root = tmp_path.resolve() / "run_00"
    target = Path("/dev/shm") / ("orthohmm_layout_" + uuid.uuid4().hex) / root.name / "input"
    config = dict(output=str(root / "output"), copy_inputs_from=str(tmp_path / "raw"),
                  copy_inputs_to=str(target))
    run = dict(measurement_directory=str(root / "measurement"), prepared_input_directory=str(target),
               configuration=config, dataset=dict(input_directory=config["copy_inputs_from"]),
               native_method="orthohmm_high_sensitivity",
               native_argv=["/python", "-m", "orthohmm", str(target), "-o", config["output"]])
    return run, root, target


@pytest.mark.parametrize("method", ["orthohmm_high_sensitivity", "orthohmm_satellite_v2", "orthofinder_full"])
def test_layout_validation_does_not_prepare_inputs_or_claim_tmpfs_availability(layout, method):
    run, root, target = layout
    run["native_method"] = method
    if method == "orthofinder_full":
        run["native_argv"] = ["/orthofinder", "-f", str(target), "-o", run["configuration"]["output"]]
    assert module.paths(run) == (root, target)
    assert not root.exists() and not target.exists()


@pytest.mark.parametrize("change", ["outside_tmpfs", "relative_input", "wrong_input_name", "wrong_run_name",
                                   "traversal", "copy_destination", "copy_source", "unknown_method",
                                   "native_input", "native_output", "persistent_output_on_tmpfs", "prephylogeny"])
def test_layout_rejection_does_not_need_or_create_tmpfs(layout, change):
    run, root, original_target = layout
    target = original_target
    if change == "outside_tmpfs":
        target = root / "input"
    elif change == "relative_input":
        target = Path("run_00/input")
    elif change == "wrong_input_name":
        target = target.with_name("wrong")
    elif change == "wrong_run_name":
        target = target.parent.with_name("other_run") / "input"
    elif change == "traversal":
        target = target.parent / ".." / target.parent.name / "input"
    if target != original_target:
        run["prepared_input_directory"] = run["configuration"]["copy_inputs_to"] = str(target)
        run["native_argv"][3] = str(target)
    if change == "copy_destination":
        run["configuration"]["copy_inputs_to"] = str(root / "different")
    elif change == "copy_source":
        run["configuration"]["copy_inputs_from"] = str(root / "different")
    elif change == "unknown_method":
        run["native_method"] = "unknown"
    elif change == "native_input":
        run["native_argv"][3] = str(root / "different")
    elif change == "native_output":
        run["native_argv"][-1] = str(root / "different")
    elif change == "persistent_output_on_tmpfs":
        forbidden_root = original_target.parent
        run["measurement_directory"] = str(forbidden_root / "measurement")
        run["configuration"]["output"] = str(forbidden_root / "output")
        run["native_argv"][-1] = run["configuration"]["output"]
    elif change == "prephylogeny":
        run["native_method"] = "orthofinder_full"
        run["native_argv"] = ["/orthofinder", "-f", str(target), "-o", run["configuration"]["output"], "-op"]
    with pytest.raises(ValueError):
        module.paths(run)
    assert not root.exists() and not original_target.exists()
