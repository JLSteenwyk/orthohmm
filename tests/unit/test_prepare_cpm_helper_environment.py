import json
from pathlib import Path

import pytest

from benchmark_tools import prepare_cpm_helper_environment as helper
from benchmark_tools import probe_cpm_private_runtime as control


def wheel_fixture(tmp_path, monkeypatch):
    assets, readers = tmp_path / "assets", tmp_path / "readers"
    (assets / "wheels").mkdir(parents=True)
    readers.mkdir()
    selected = {"numpy": "2.2.6", **{f"example{i}": "1" for i in range(10)}}
    rows, metadata = [], {}
    for directory, packages in ((assets / "wheels", selected),
                                 (readers, {"biopython": "1.87", "numpy": "2.2.6"})):
        for name, version in packages.items():
            path = directory / (name + ".whl")
            path.write_bytes(name.encode())
            rows.append(control.record(path))
            metadata[path] = (name, version)
    monkeypatch.setattr(helper, "wheel_metadata", lambda p: metadata[p])
    plan = dict(command=["python", "--assets", str(assets), "--reader-wheels", str(readers)],
                checked_records=rows)
    audit = dict(inference=dict(inventory=dict(package=[dict(name=k, version=v)
                                                       for k, v in selected.items()])))
    return plan, audit, metadata


def test_only_inference_wheels_and_reader_biopython_selected(tmp_path, monkeypatch):
    plan, audit, _ = wheel_fixture(tmp_path, monkeypatch)
    selected, wheels = helper.select_wheels(plan, audit)
    assert len(selected) == len(wheels) == 12
    assert selected["biopython"] == "1.87"
    assert sum(row["name"] == "numpy" for row in wheels) == 1


@pytest.mark.parametrize("bad", ["missing", "duplicate", "version", "changed", "inventory"])
def test_selection_failures_rejected(tmp_path, monkeypatch, bad):
    plan, audit, metadata = wheel_fixture(tmp_path, monkeypatch)
    if bad == "missing":
        plan["checked_records"].pop(0)
    elif bad == "duplicate":
        plan["checked_records"].append(plan["checked_records"][0])
    elif bad == "version":
        metadata[Path(plan["checked_records"][0]["path"])] = ("numpy", "wrong")
    elif bad == "changed":
        Path(plan["checked_records"][0]["path"]).write_bytes(b"changed")
    else:
        audit["inference"]["inventory"]["package"].pop()
    with pytest.raises(ValueError):
        helper.select_wheels(plan, audit)


def imports(tmp_path):
    prefix, launcher, base = [tmp_path / name for name in ("venv", "source", "base")]
    modules = {name: str((prefix if name == "numpy" else launcher) / (name + ".py"))
               for name in helper.IMPORTS}
    observed = dict(packages={"NumPy": "2.2.6"}, base_prefix=str(base),
                    modules=modules, gc_enabled=True, gc_thresholds=[700, 10, 10])
    return observed, dict(numpy="2.2.6"), prefix, launcher, base


def test_exact_worker_import_roots_and_package_inventory(tmp_path):
    helper.validate_imports(*imports(tmp_path))


@pytest.mark.parametrize("bad", ["extra_package", "wrong_base", "gc_disabled", "gc_thresholds",
                               "missing_import", "outside_module", "nonfrozen_module"])
def test_bad_import_closure_rejected(tmp_path, bad):
    observed, selected, prefix, launcher, base = imports(tmp_path)
    if bad == "extra_package": observed["packages"]["Other"] = "1"
    elif bad == "wrong_base": observed["base_prefix"] = "/other"
    elif bad == "gc_disabled": observed["gc_enabled"] = False
    elif bad == "gc_thresholds": observed["gc_thresholds"] = [1, 2, 3]
    elif bad == "missing_import": observed["modules"].pop("numpy")
    elif bad == "outside_module": observed["modules"]["numpy"] = "/shared/numpy.py"
    else: observed["modules"]["orthohmm.accuracy"] = str(prefix / "accuracy.py")
    with pytest.raises(ValueError):
        helper.validate_imports(observed, selected, prefix, launcher, base)


def prepared_fixture(tmp_path):
    prefix = tmp_path / "venv"
    (prefix / "bin").mkdir(parents=True)
    python = prefix / "bin/python"
    python.write_bytes(b"interpreter")
    log = tmp_path / "imports"
    log.write_bytes(b"{}\n")
    report = dict(status="helper_complete_private_environment_prepared", refinement_attempts=0,
                  seed_admitted=False, accuracy_evaluated=False, publication_ready=False,
                  prefix=str(prefix), interpreter=control.record(python),
                  checked_records=[control.record(python)], import_report=control.record(log),
                  base_python=control.record(python))
    receipt = tmp_path / "result.json"
    receipt.write_text(json.dumps(report))
    return receipt, report


def test_prepared_receipt_is_hash_bound_and_rechecked(tmp_path):
    receipt, report = prepared_fixture(tmp_path)
    expected = control.record(receipt)
    assert control.prepared_runtime(receipt, expected["sha256"]) == (report, expected)
    Path(report["interpreter"]["path"]).write_bytes(b"changed")
    with pytest.raises(ValueError, match="interpreter"):
        control.prepared_runtime(receipt, expected["sha256"])


@pytest.mark.parametrize("bad", ["hash", "status", "attempts", "admitted", "log"])
def test_bad_prepared_receipt_rejected(tmp_path, bad):
    receipt, report = prepared_fixture(tmp_path)
    if bad == "status": report["status"] = "helper_environment_preparation_failed"
    elif bad == "attempts": report["refinement_attempts"] = 1
    elif bad == "admitted": report["seed_admitted"] = True
    elif bad == "log": Path(report["import_report"]["path"]).write_bytes(b"changed")
    receipt.write_text(json.dumps(report))
    sha = "bad" if bad == "hash" else control.record(receipt)["sha256"]
    with pytest.raises(ValueError):
        control.prepared_runtime(receipt, sha)


def test_prepared_receipt_and_hash_must_be_supplied_together(tmp_path, monkeypatch):
    monkeypatch.setattr(control.os, "sched_setaffinity", lambda *_: None)
    with pytest.raises(ValueError, match="both"):
        control.run(tmp_path, tmp_path / "output", "protocol", tmp_path / "receipt")
    assert json.loads((tmp_path / "output/report.json").read_bytes())["refinement_attempts"] == 0


def test_helper_output_cannot_be_overwritten(tmp_path):
    with pytest.raises(FileExistsError):
        helper.prepare(tmp_path, tmp_path, "unused")
