import io
import json
from types import SimpleNamespace

import pytest

from benchmark_tools import prepare_publication_bootstrap as module


@pytest.fixture
def inputs(tmp_path, monkeypatch):
    installer = tmp_path / module.NAME
    installer.write_bytes(b"Synthetic bootstrap, never executed")
    monkeypatch.setattr(module, "EXPECTED", module.identity(installer))
    checksum = b"Synthetic fixed checksum"
    path = tmp_path / (module.NAME + ".sha256")
    path.write_bytes(checksum)
    monkeypatch.setattr(module, "CHECKSUM", module.identity(path))
    monkeypatch.setattr(module.platform, "system", lambda: "Linux")
    monkeypatch.setattr(module.platform, "machine", lambda: "x86_64")
    args = SimpleNamespace(installer=installer, output=tmp_path / "fresh bootstrap", timeout=60,
                           acknowledge_historical_runtime=True)
    return args, checksum


def execution(args, monkeypatch, fail=None):
    calls = []
    def stage(output, name, command, env, timeout):
        calls.append((name, command, env, timeout))
        if fail == name:
            raise RuntimeError("Injected bootstrap stage failure")
        if name == "install":
            prefix = output / "prefix"
            (prefix / "bin").mkdir(parents=True)
            (prefix / "bin/conda").write_text("Synthetic conda, never executed")
            (prefix / "bin/conda").chmod(0o755)
            (prefix / "bin/python-real").write_text("Synthetic Python, never executed")
            (prefix / "bin/python").symlink_to("python-real")
            (prefix / "conda-meta").mkdir()
            (prefix / "conda-meta/conda.json").write_text(json.dumps(dict(name="conda", version="25.3.1", build="fixture")))
        else:
            (output / "conda_version.log").write_text("conda 25.3.1\n")
        return dict(returncode=0, command=command)
    monkeypatch.setattr(module, "stage", stage)
    return calls


def test_mocked_private_bootstrap_install(inputs, monkeypatch):
    args, _ = inputs
    calls = execution(args, monkeypatch)
    result = module.install(args)
    assert result["status"] == "pinned_private_bootstrap_installed"
    assert [c[0] for c in calls] == ["install", "conda_version"]
    assert calls[0][1] == ["/bin/bash", str(args.installer), "-b", "-p", str(args.output / "prefix")]
    assert calls[0][2]["HOME"] == str(args.output / "home") and calls[0][2]["CONDARC"] == "/dev/null"
    assert calls[0][2]["CONDA_EXTRACT_THREADS"] == "1"
    assert "PYTHONPATH" not in calls[0][2] and "LD_PRELOAD" not in calls[0][2]
    assert all(result[k] is False for k in ("retry", "shell_init_requested", "scientific_environment_installed",
        "native_inference_executed", "publication_ready", "security_clearance", "redistribution_clearance"))


@pytest.mark.parametrize("fault", ["ack", "timeout", "name", "hash", "symlink", "host", "existing"])
def test_install_preflight(inputs, tmp_path, monkeypatch, fault):
    args, _ = inputs
    calls = execution(args, monkeypatch)
    if fault == "ack": args.acknowledge_historical_runtime = False
    elif fault == "timeout": args.timeout = True
    elif fault == "name":
        path = tmp_path / "wrong.sh"
        path.write_bytes(args.installer.read_bytes())
        args.installer = path
    elif fault == "hash": args.installer.write_text("changed")
    elif fault == "symlink":
        path = tmp_path / "link.sh"
        path.symlink_to(args.installer)
        args.installer = path
    elif fault == "host": monkeypatch.setattr(module.platform, "machine", lambda: "aarch64")
    else: args.output.mkdir()
    with pytest.raises((ValueError, FileExistsError)):
        module.install(args)
    assert not calls
    if fault != "existing": assert not args.output.exists()


@pytest.mark.parametrize("name", ["install", "conda_version"])
def test_stage_failure_retained_without_retry(inputs, monkeypatch, name):
    args, _ = inputs
    calls = execution(args, monkeypatch, name)
    with pytest.raises(RuntimeError): module.install(args)
    assert calls[-1][0] == name and len({c[0] for c in calls}) == len(calls)
    result = json.loads((args.output / "failed.json").read_bytes())
    assert result["retry"] is False and result["attempts"] == 1
    assert not (args.output / "complete.json").exists()


@pytest.mark.parametrize("url", ["http://github.com/a", "https://evil.example/a", "https://github.com/a?token=x",
    "https://user:password@github.com/a", "https://github.com:444/a", "https://github.com/a#x"])
def test_unsafe_providers(url):
    with pytest.raises(ValueError): module.provider(url)


def test_signed_asset_query_is_supported():
    assert module.provider("https://release-assets.githubusercontent.com/a?temporary=x") == "release-assets.githubusercontent.com"


def mock_downloads(args, checksum, monkeypatch, corrupt=False):
    calls = []
    installer = args.installer.read_bytes()
    class Response(io.BytesIO):
        def __init__(self, payload):
            super().__init__(payload)
            self.headers = {}
        def geturl(self): return "https://release-assets.githubusercontent.com/a?temporary=do-not-record"
    class Opener:
        def open(self, req, timeout):
            calls.append(req.full_url)
            data = checksum if req.full_url.endswith(".sha256") else installer
            return Response(data + b"x" if corrupt else data)
    monkeypatch.setattr(module, "build_opener", lambda *args: Opener())
    return calls


def test_mocked_acquisition_excludes_signed_queries_from_receipt(inputs, monkeypatch):
    args, checksum = inputs
    calls = mock_downloads(args, checksum, monkeypatch)
    result = module.acquire(args)
    assert len(calls) == 2 and result["status"] == "pinned_bootstrap_acquired"
    assert "temporary" not in json.dumps(result)
    assert result["installation_performed"] is False


def test_corrupt_download_retained_without_install(inputs, monkeypatch):
    args, checksum = inputs
    calls = mock_downloads(args, checksum, monkeypatch, True)
    with pytest.raises(ValueError): module.acquire(args)
    assert len(calls) == 1
    assert (args.output / (module.NAME + ".sha256.partial")).exists()
    assert not (args.output / "complete.json").exists()
    assert json.loads((args.output / "failed.json").read_bytes())["retry"] is False
