import pytest

from benchmark_tools.probe_verified_threadripper_fixture import METHODS, select_fixture
from benchmark_tools import probe_verified_threadripper_fixture as module


@pytest.mark.parametrize("index,method", enumerate(METHODS))
def test_selects_frozen_method_without_changes(index, method):
    rows = [dict(index=i, native_method=m, native_argv=[m, "frozen-flag"]) for i, m in enumerate(METHODS)]
    assert select_fixture(dict(runs=rows), method) is rows[index]


def test_wrong_method_rejected():
    with pytest.raises(ValueError):
        select_fixture(dict(runs=[]), "unknown")


def test_reordered_block_rejected():
    rows = [dict(index=i, native_method=m) for i, m in enumerate(reversed(METHODS))]
    with pytest.raises(ValueError):
        select_fixture(dict(runs=rows), METHODS[0])


@pytest.mark.parametrize("option", ["--lookup", "--lookup-sha256", "--plan", "--plan-sha256"])
def test_partial_deployment_options_rejected(monkeypatch, option):
    monkeypatch.setattr("sys.argv", ["fixture", "--input", "/unused", "--output", "/unused",
                                   "--tmpfs", "/unused", "--scheduler-command", "none", option, "value"])
    with pytest.raises(SystemExit) as error:
        module.main()
    assert error.value.code == 2


def test_wrong_bound_plan_rejected_before_output(monkeypatch, tmp_path):
    output = tmp_path / "absent"
    monkeypatch.setattr("sys.argv", ["fixture", "--input", "/unused", "--output", str(output),
        "--tmpfs", "/unused", "--scheduler-command", "none", "--lookup", "/lookup",
        "--lookup-sha256", "lookuphash", "--plan", "/plan", "--plan-sha256", "planhash"])
    values = iter([dict(baseline=dict(path="/baseline", sha256="b"), binding=dict(path="/binding", sha256="c")),
                   {}, dict(command_plan={"path": "/different"}), {}])
    monkeypatch.setattr(module, "read_frozen", lambda *a: next(values))
    monkeypatch.setattr(module, "record", lambda *a: {"path": "/plan"})
    with pytest.raises(ValueError, match="Fixture plan"):
        module.main()
    assert not output.exists()
