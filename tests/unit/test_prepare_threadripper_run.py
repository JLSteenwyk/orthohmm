from pathlib import Path
import tempfile

import pytest

from benchmark_tools import prepare_threadripper_run as module


@pytest.fixture
def setup(tmp_path, monkeypatch):
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
