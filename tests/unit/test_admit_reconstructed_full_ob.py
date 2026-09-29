import pytest

from benchmark_tools import admit_reconstructed_full_ob as module


def native_files(tmp_path):
    old, new = tmp_path / "old", tmp_path / "new"
    old.mkdir()
    new.mkdir()
    pins = {}
    for name in module.NATIVE_FILES:
        for directory in (old, new):
            (directory / name).write_text("header\noriginal\n")
        pin = module.record(old / name)
        pins[pin["path"]] = pin
    return old, new, pins


def test_all_native_files_compared(tmp_path):
    old, new, pins = native_files(tmp_path)
    result = module.native_comparison(old, new, pins)
    assert len(result) == 4
    assert all(row["byte_equal"] for row in result.values())


@pytest.mark.parametrize("filename", module.NATIVE_FILES)
def test_changed_current_output_retained(tmp_path, filename):
    old, new, pins = native_files(tmp_path)
    (new / filename).write_text("header\nchanged\n")
    result = module.native_comparison(old, new, pins)
    assert not result[filename]["byte_equal"]
    assert sum(row["byte_equal"] for row in result.values()) == 3
    assert result[filename]["current"] == module.record(new / filename)


def test_unpinned_baseline_rejected(tmp_path):
    old, new, pins = native_files(tmp_path)
    pins.pop(next(iter(pins)))
    with pytest.raises(ValueError, match="not bound"):
        module.native_comparison(old, new, pins)


def test_modified_historical_output_rejected(tmp_path):
    old, new, pins = native_files(tmp_path)
    (old / module.NATIVE_FILES[0]).write_text("tampered")
    with pytest.raises(ValueError, match="not bound"):
        module.native_comparison(old, new, pins)


def test_missing_current_file_not_equal(tmp_path):
    old, new, pins = native_files(tmp_path)
    (new / module.NATIVE_FILES[0]).unlink()
    with pytest.raises(FileNotFoundError):
        module.native_comparison(old, new, pins)


def test_existing_audit_not_overwritten(tmp_path):
    with pytest.raises(FileExistsError):
        module.audit(tmp_path / "unused", tmp_path)


def test_scheduler_rejection_precedes_output_creation(tmp_path, monkeypatch):
    def reject(directory):
        raise ValueError("Not completed")
    monkeypatch.setattr(module, "verify", reject)
    output = tmp_path / "audit"
    with pytest.raises(ValueError, match="Not completed"):
        module.audit(tmp_path, output)
    assert not output.exists()


def test_package_paths_come_from_plan(tmp_path, monkeypatch):
    assets, readers = tmp_path / "old-assets", tmp_path / "old-readers"
    lock = assets / "benchmark_tools/results/publication_recovery_requirements_20260926.txt"
    lock.parent.mkdir(parents=True)
    lock.write_text("inference-lock")
    reader_lock = tmp_path / "reader-lock.txt"
    reader_lock.write_text("reader-lock")
    root = tmp_path / "new-run"
    root.mkdir()
    for name in ("inference", "reader"):
        module.save(root / (name + "_install.json"), dict(name=name))
    seen = []
    def wheels(report, directory):
        seen.append((report["name"], directory))
        return []
    locks = []
    monkeypatch.setattr(module, "local_install_wheels", wheels)
    monkeypatch.setattr(module, "verify_lock", lambda text, rows: locks.append(text))
    monkeypatch.setattr(module.subprocess, "check_output", lambda *a, **k: "[]")
    monkeypatch.setattr(module, "verify_inventory", lambda *a: {"verified": True})
    result = module.package_audit(root, dict(command=[
        "--assets", str(assets), "--reader-wheels", str(readers), "--reader-lock", str(reader_lock)]))
    assert seen == [("inference", assets / "wheels"), ("reader", readers)]
    assert locks == ["inference-lock", "reader-lock"]
    assert set(result) == {"inference", "reader"}


def test_changed_baseline_receipt_rejected(tmp_path, monkeypatch):
    monkeypatch.setattr(module, "RESULTS", tmp_path)
    (tmp_path / "integrated_full_ob_result_22337.json").write_text("{}")
    with pytest.raises(ValueError, match="Changed admitted baseline"):
        module.baseline()


def test_all_four_readers_required():
    with pytest.raises(ValueError, match="all four"):
        module.report_pins(dict(reports={}))


def test_conflicting_reader_pins_rejected(tmp_path):
    reports = {}
    for index, name in enumerate(("structure", "sequences", "events", "hierarchy")):
        path = tmp_path / (name + ".json")
        module.save(path, dict(checked_records=[dict(path="same", bytes=index, sha256="a")]))
        reports[name] = module.record(path)
    with pytest.raises(ValueError, match="Conflicting"):
        module.report_pins(dict(reports=reports))
