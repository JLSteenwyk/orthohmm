import json
from pathlib import Path
import subprocess
import sys

import pytest

from benchmark_tools import probe_cpm_parser_imports as module


def prepare(tmp_path, monkeypatch):
    results = tmp_path / "benchmark_tools/results"
    results.mkdir(parents=True)
    prior = results / "qfo_cpm_parser_controls_20260926.json"
    prior.write_text(json.dumps({"checked_records": []}))
    monkeypatch.setattr(module, "PRIOR_SHA", module.helper.record(prior)["sha256"])
    status = tmp_path / "status.json"
    status.write_text(json.dumps(dict(scientific_child_command=[sys.executable],
        runtime_before={"files": []}, checked_records=[])))
    monkeypatch.setattr(module.helper, "STATUS", status.name)
    (tmp_path / module.PROTOCOL).write_text("fixture protocol\n")
    (tmp_path / "benchmarks/work/publication_qfo_replay_native_v1").mkdir(parents=True)
    return prior


@pytest.mark.parametrize("outcome", ["success", "signal", "timeout", "bad_json"])
def test_fixed_arms_and_failure_retention(tmp_path, monkeypatch, outcome):
    prepare(tmp_path, monkeypatch)
    calls = []
    def execute(command, **kwargs):
        calls.append(command[-1])
        assert "-S" not in command
        assert kwargs["env"]["PYTHONMALLOC"] == "debug"
        assert kwargs["env"]["OPENBLAS_NUM_THREADS"] == "1"
        assert kwargs["cwd"].name == "publication_qfo_replay_native_v1"
        assert kwargs["timeout"] == 120
        if outcome == "timeout":
            raise subprocess.TimeoutExpired(command, 120, output=b"partial", stderr=b"before_parser")
        return subprocess.CompletedProcess(command, -11 if outcome == "signal" else 0,
            stdout=b"bad" if outcome == "bad_json" else b"{}", stderr=b"before_parser")
    monkeypatch.setattr(module.subprocess, "run", execute)
    output = tmp_path / "output"
    if outcome == "bad_json":
        with pytest.raises(json.JSONDecodeError):
            module.run(tmp_path, output)
        report = json.loads((output / "report.json").read_text())
        assert report["status"] == "import_controls_failed"
        assert calls == ["site_only"]
    else:
        report = module.run(tmp_path, output)
        assert calls == ["site_only", "frozen_imports"]
        assert all(a["status"] == {"success": "completed", "signal": "failed", "timeout": "timed_out"}[outcome]
                   for a in report["arms"])
    assert report["accuracy_admitted"] is False
    assert all(a["attempts"] == 1 for a in report["arms"])
    assert all(Path(a["stderr"]["path"]).read_bytes() == b"before_parser" for a in report["arms"])
    with pytest.raises(FileExistsError):
        module.run(tmp_path, output)


def test_changed_prior_prevents_execution(tmp_path, monkeypatch):
    prior = prepare(tmp_path, monkeypatch)
    prior.write_text("{}")
    monkeypatch.setattr(module.subprocess, "run", lambda *a, **k: pytest.fail("unexpected launch"))
    with pytest.raises(ValueError, match="Changed parser-only"):
        module.run(tmp_path, tmp_path / "output")
    assert not (tmp_path / "output").exists()


@pytest.mark.parametrize("mode", ["site_only", "frozen_imports"])
def test_worker_controls_import_boundary_and_limits(tmp_path, monkeypatch, capsys, mode):
    limits, imports = [], []
    monkeypatch.setattr(module.resource, "setrlimit", lambda *args: limits.append(args))
    monkeypatch.setattr(module.helper, "pinned_inputs", lambda root: [{"path": str(tmp_path)}] * 3)
    monkeypatch.setattr(module.helper, "check", lambda inputs: None)
    monkeypatch.setattr(module.helper, "parse_only", lambda *args: dict(genes=984137, groups=390845, memberships=984137))
    monkeypatch.setattr(module, "scientific_imports", lambda path: imports.append(path) or [{"path": "frozen"}])
    module.worker(tmp_path, mode)
    output = capsys.readouterr()
    assert output.err.splitlines() == ["before_imports", "before_parser", "after_parser"]
    assert len(imports) == (1 if mode == "frozen_imports" else 0)
    assert json.loads(output.out)["mode"] == mode
    assert limits == [(module.resource.RLIMIT_AS, (4 * 1024**3, 4 * 1024**3)),
                      (module.resource.RLIMIT_CPU, (120, 120)), (module.resource.RLIMIT_CORE, (0, 0))]
