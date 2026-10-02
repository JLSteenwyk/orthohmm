import io
import json
from pathlib import Path
from types import SimpleNamespace
from urllib.request import Request

import pytest

from benchmark_tools import acquire_publication_base as module


class Response(io.BytesIO):
    def __init__(self, payload, url, headers=None):
        super().__init__(payload)
        self.url = url
        self.headers = headers if headers is not None else {"Content-Length": str(len(payload))}

    def geturl(self):
        return self.url


@pytest.fixture
def acquisition(tmp_path, monkeypatch):
    packages, payloads = [], {}
    for i in range(19):
        filename = f"package{i}-1-main.conda"
        payload = f"Synthetic archive {i}".encode()
        path = tmp_path / filename
        path.write_bytes(payload)
        expected = module.archives.identity(path)
        url = "https://repo.anaconda.com/pkgs/main/linux-64/" + filename
        packages.append(dict(package=dict(name=f"package{i}", fn=filename, subdir="linux-64",
            version="1", build="main", url=url, sha256=expected["sha256"], md5=expected["md5"]),
            archive=dict(bytes=len(payload), sha256=expected["sha256"])))
        payloads[url] = payload
    receipt = tmp_path / "receipt.json"
    receipt.write_text(json.dumps(dict(acquisition=dict(packages=packages))))
    monkeypatch.setattr(module, "RECEIPT_SHA", module.record(receipt)["sha256"])
    wheel = b"Synthetic wheel; never executed"
    monkeypatch.setattr(module, "PIP_BYTES", len(wheel))
    monkeypatch.setattr(module, "PIP_SHA", module.hashlib.sha256(wheel).hexdigest())
    pip_url = "https://files.pythonhosted.org/packages/fixed/" + module.PIP_NAME
    metadata = dict(info=dict(name="pip", version="26.2.1"), urls=[dict(filename=module.PIP_NAME,
        packagetype="bdist_wheel", yanked=False, size=len(wheel), digests=dict(sha256=module.PIP_SHA), url=pip_url)])
    payloads[module.PIP_METADATA] = json.dumps(metadata).encode()
    payloads[pip_url] = wheel
    requests = []

    class Opener:
        def open(self, request, timeout):
            requests.append((request.full_url, timeout, request.get_header("Accept-encoding")))
            return Response(payloads[request.full_url], request.full_url)

    opener = Opener()
    monkeypatch.setattr(module, "network_opener", lambda: opener)
    args = SimpleNamespace(receipt=receipt, output=tmp_path / "fresh acquisition", timeout=60,
                           acknowledge_historical_runtime=True)
    return args, opener, requests, payloads, metadata


def test_acquisition_exact_scope_without_execution(acquisition):
    args, _, requests, _, _ = acquisition
    result = module.acquire(args)
    assert result["status"] == "historical_base_artifacts_acquired"
    assert result["base_packages"] == 19 and len(requests) == len(result["downloads"]) == 21
    assert len(list((args.output / "archives").iterdir())) == 19
    assert all(r[1:] == (60, "identity") for r in requests)
    assert result == json.loads((args.output / "complete.json").read_bytes())
    assert all(result[k] is False for k in ("installation_performed", "native_inference_executed",
        "publication_ready", "security_clearance", "redistribution_clearance", "retry"))
    assert not list(args.output.rglob("*.partial"))


@pytest.mark.parametrize("defect", ["ack", "timeout_zero", "timeout_bool", "changed_receipt", "symlink", "existing", "canonical"])
def test_preflight_before_network_or_creation(acquisition, tmp_path, defect):
    args, _, requests, _, _ = acquisition
    if defect == "ack":
        args.acknowledge_historical_runtime = False
    elif defect == "timeout_zero":
        args.timeout = 0
    elif defect == "timeout_bool":
        args.timeout = True
    elif defect == "changed_receipt":
        args.receipt.write_text("{}")
    elif defect == "symlink":
        alias = tmp_path / "receipt-alias"
        alias.symlink_to(args.receipt)
        args.receipt = alias
    elif defect == "existing":
        args.output.mkdir()
    else:
        alias = tmp_path / "alias"
        alias.symlink_to(tmp_path, target_is_directory=True)
        args.output = alias / "fresh acquisition"
    with pytest.raises((ValueError, FileExistsError)):
        module.acquire(args)
    assert not requests
    if defect != "existing":
        assert not args.output.exists()


@pytest.mark.parametrize("url", ["http://pypi.org/x", "https://evil.example/x", "https://pypi.org:444/x",
    "https://user:password@pypi.org/x", "https://pypi.org/x?token=secret", "https://pypi.org/x#fragment",
    "https://pypi.org/\nx", "https://pypi.org/ x", 1])
def test_provider_url_and_redirect_rejected(url):
    with pytest.raises((ValueError, TypeError)):
        module.provider_url(url)
    with pytest.raises((ValueError, TypeError)):
        module.ProviderRedirect().redirect_request(Request(module.PIP_METADATA), None, 302, "Found", {}, url)


def test_https_redirect_to_allowed_provider():
    result = module.ProviderRedirect().redirect_request(Request(module.PIP_METADATA), None, 302, "Found", {},
        "https://files.pythonhosted.org/packages/x")
    assert result.full_url == "https://files.pythonhosted.org/packages/x"


def test_short_metadata_response_is_not_accepted(acquisition, tmp_path, monkeypatch):
    args, opener, _, _, _ = acquisition
    response = Response(b"{}", module.PIP_METADATA, {"Content-Length": "3"})
    monkeypatch.setattr(opener, "open", lambda *args, **kwargs: response)
    with pytest.raises(ValueError, match="declared size"):
        module.download(opener, module.PIP_METADATA, tmp_path / "metadata.json", 100, args.timeout)
    assert not (tmp_path / "metadata.json").exists()


@pytest.mark.parametrize("defect", ["short", "long", "hash", "network", "metadata", "wheel", "changed_input"])
def test_failures_retained_without_retry(acquisition, monkeypatch, defect):
    args, opener, requests, payloads, _ = acquisition
    first = next(iter(payloads))
    if defect == "short":
        payloads[first] = payloads[first][:-1]
    elif defect == "long":
        payloads[first] += b"x"
    elif defect == "hash":
        payloads[first] = b"x" * len(payloads[first])
    elif defect == "metadata":
        payloads[module.PIP_METADATA] = b"{}"
    elif defect == "wheel":
        url = next(u for u in payloads if u.startswith("https://files.pythonhosted.org/"))
        payloads[url] = b"x" * len(payloads[url])
    elif defect == "network":
        def fail(*args, **kwargs):
            raise OSError("Injected network failure")
        monkeypatch.setattr(opener, "open", fail)
    else:
        download = module.download
        def change_input(*values, **kwargs):
            result = download(*values, **kwargs)
            args.receipt.write_text("changed during acquisition")
            return result
        monkeypatch.setattr(module, "download", change_input)
    with pytest.raises((ValueError, KeyError, OSError)):
        module.acquire(args)
    failed = json.loads((args.output / "failed.json").read_bytes())
    assert failed["retry"] is False and failed["attempts"] == 1
    assert not (args.output / "complete.json").exists()
    assert len({r[0] for r in requests}) == len(requests)
    with pytest.raises(FileExistsError):
        module.acquire(args)


@pytest.mark.parametrize("defect", ["short_without_length", "oversize_without_length", "redirect", "encoding"])
def test_stream_bounds_without_trusting_headers(acquisition, tmp_path, monkeypatch, defect):
    args, opener, _, payloads, _ = acquisition
    url = next(iter(payloads))
    original = payloads[url]
    response = Response(original[:-1] if defect == "short_without_length" else original + b"x",
        "https://evil.example/x" if defect == "redirect" else url,
        {"Content-Encoding": "gzip"} if defect == "encoding" else {})
    monkeypatch.setattr(opener, "open", lambda *args, **kwargs: response)
    expected = dict(bytes=len(original), sha256=module.hashlib.sha256(original).hexdigest())
    with pytest.raises(ValueError):
        module.download(opener, url, tmp_path / "download", len(original), args.timeout, expected)
    assert not (tmp_path / "download").exists()
    assert (tmp_path / "download.partial").exists()


@pytest.mark.parametrize("defect", ["version", "name", "duplicate", "missing", "size", "sha", "yanked", "kind", "host", "filename"])
def test_pip_metadata_must_match_frozen_artifact(acquisition, defect):
    _, _, _, _, metadata = acquisition
    row = metadata["urls"][0]
    if defect in {"version", "name"}:
        metadata["info"][defect] = "different"
    elif defect == "duplicate":
        metadata["urls"].append(dict(row))
    elif defect == "missing":
        metadata["urls"] = []
    elif defect == "size":
        row["size"] += 1
    elif defect == "sha":
        row["digests"]["sha256"] = "0" * 64
    elif defect == "yanked":
        row["yanked"] = True
    elif defect == "kind":
        row["packagetype"] = "sdist"
    elif defect == "host":
        row["url"] = "https://repo.anaconda.com/" + module.PIP_NAME
    else:
        row["url"] = "https://files.pythonhosted.org/packages/different.whl"
    with pytest.raises(ValueError):
        module.pip_selection(metadata)
