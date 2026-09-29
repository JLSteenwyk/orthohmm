import json
from pathlib import Path
from types import SimpleNamespace

import pytest

from benchmark_tools import check_threadripper_runtime as module


@pytest.fixture
def setup(tmp_path):
    def save(name, value):
        path = tmp_path / name
        path.write_text(json.dumps(value))
        return module.record(path)
    report = dict(python="version", executable="/python", cwd="/core", requested=["core"],
        paths=[], modules={}, mapped_files=[], meta_path=[], path_hooks=[], editable={}, files=[],
        dont_write_bytecode=True, coverage=dict(missing=[], changed=[], all_covered=True),
        scientific_origin=dict(root="/core"))
    expected = save("expected.json", report)
    baseline = save("baseline.json", dict(tool_entrypoints=dict(orthohmm_python=dict(absolute_path="/python"))))
    specs = [["/runtime", "hash"]]
    binding = save("binding.json", dict(runtime_specs=specs))
    source = module.record(Path(module.__file__).with_name("inspect_native_python_lookup.py"))
    receipt = save("receipt.json", dict(status="native_lookup_repeated_identity_match", baseline=baseline,
        binding=binding, source=source, interpreters={n: dict(reports=[expected]) for n in ("orthohmm", "orthofinder")}))
    calls = []
    def runner(command, **kwargs):
        calls.append(command)
        directory = Path(command[-1])
        directory.mkdir()
        for name in ("orthohmm", "orthofinder"):
            (directory / (name + ".json")).write_text(json.dumps(report))
        return SimpleNamespace(returncode=0)
    checker = module.RuntimeChecker(receipt["path"], receipt["sha256"], tmp_path / "checks",
                                    tree_checker=lambda s: {"checked": s}, runner=runner)
    return checker, specs, report, calls


def test_successive_before_after_checks(setup):
    checker, specs, _, calls = setup
    assert checker(specs)["status"] == "runtime_and_lookup_checked"
    assert checker(specs)["status"] == "runtime_and_lookup_checked"
    assert len(calls) == 2 and calls[0][-1] != calls[1][-1]
    assert (checker.output / "checked_02.json").is_file()
    assert calls[0][0] == "/python"


def use_controller(checker, tmp_path):
    controller = tmp_path / "controller-python"
    controller.write_bytes(b"private controller")
    receipt = json.loads(checker.receipt_path.read_text())
    binding_path = Path(receipt["binding"]["path"])
    binding = json.loads(binding_path.read_text())
    binding["controller_python"] = module.record(controller)
    binding_path.write_text(json.dumps(binding))
    receipt["binding"] = module.record(binding_path)
    checker.receipt_path.write_text(json.dumps(receipt))
    checker.receipt_sha256 = module.record(checker.receipt_path)["sha256"]
    return controller


def test_private_controller_selected(setup, tmp_path):
    checker, specs, _, calls = setup
    controller = use_controller(checker, tmp_path)
    checker(specs)
    assert calls[0][0] == str(controller)


def test_controller_drift_prevents_launch(setup, tmp_path):
    checker, specs, _, calls = setup
    controller = use_controller(checker, tmp_path)
    controller.write_bytes(b"changed")
    with pytest.raises(ValueError):
        checker(specs)
    assert calls == []


def test_controller_drift_during_inspection_rejected(setup, tmp_path):
    checker, specs, _, _ = setup
    controller = use_controller(checker, tmp_path)
    previous = checker.runner
    def runner(*args, **kwargs):
        result = previous(*args, **kwargs)
        controller.write_bytes(b"changed")
        return result
    checker.runner = runner
    with pytest.raises(ValueError):
        checker(specs)
    assert not (checker.output / "checked_01.json").exists()


def test_wrong_specs_prevent_inspector(setup):
    checker, _, _, calls = setup
    with pytest.raises(ValueError, match="specifications"):
        checker([["/other", "hash"]])
    assert not calls


def test_lookup_drift_preserves_evidence(setup):
    checker, specs, report, calls = setup
    checker(specs)
    report["modules"] = {"unexpected": "/elsewhere"}
    with pytest.raises(ValueError, match="lookup"):
        checker(specs)
    assert len(calls) == 2
    assert (checker.output / "check_02/orthohmm.json").is_file()
    assert not (checker.output / "checked_02.json").exists()


def test_tree_failure_prevents_inspector(setup):
    checker, specs, _, calls = setup
    def fail(_):
        raise ValueError("runtime changed")
    checker.tree_checker = fail
    with pytest.raises(ValueError, match="runtime changed"):
        checker(specs)
    assert not calls


def test_inspector_failure_retained(setup):
    checker, specs, _, _ = setup
    checker.runner = lambda *a, **k: SimpleNamespace(returncode=7)
    with pytest.raises(RuntimeError):
        checker(specs)
    assert json.loads((checker.output / "process_01.json").read_text())["exit_code"] == 7


@pytest.mark.parametrize("phase", ["before", "after"])
def test_composed_wrapper_blocks_or_preserves_measurement(setup, tmp_path, phase):
    from benchmark_tools.measure_native_scaling_run import run_checked
    checker, specs, report, _ = setup
    launches = []
    if phase == "before":
        report["modules"] = {"unexpected": "/elsewhere"}
    def native(_):
        launches.append(True)
        report["modules"] = {"unexpected": "/elsewhere"}
        return dict(status="command_exited_zero", native=dict(exit_code=0))
    result = run_checked(specs, tmp_path / "run", native, checker)
    if phase == "before":
        assert not launches and "measurement" not in result
        assert result["status"] == "verified_wrapper_failed"
    else:
        assert launches == [True] and result["measurement"]["native"]["exit_code"] == 0
        assert result["status"] == "runtime_changed_or_unverifiable"
