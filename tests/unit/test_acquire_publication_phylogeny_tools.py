import io
import json
from pathlib import Path
from types import SimpleNamespace
from urllib.request import Request

import pytest

from benchmark_tools import acquire_publication_phylogeny_tools as module


class Response(io.BytesIO):
    def __init__(self, payload, url, headers=None):
        super().__init__(payload)
        self.url = url
        self.headers = {"Content-Length": str(len(payload))} if headers is None else headers

    def geturl(self):
        return self.url


@pytest.fixture
def acquisition(tmp_path, monkeypatch):
    rows = module.artifacts()
    payloads = {}
    for row in rows:
        payload = ("Synthetic, never executed: " + row["relative"]).encode()
        row.update(bytes=len(payload), sha256=module.base.hashlib.sha256(payload).hexdigest())
        payloads[row["url"]] = payload
    requests = []

    class Opener:
        def open(self, request, timeout):
            requests.append((request.full_url, timeout, request.get_header("Accept-encoding")))
            return Response(payloads[request.full_url], request.full_url)

    opener = Opener()
    monkeypatch.setattr(module, "network_opener", lambda: opener)
    monkeypatch.setattr(module, "artifacts", lambda: rows)
    args = SimpleNamespace(output=tmp_path / "fresh artifacts", timeout=60,
                           acknowledge_historical_runtime=True)
    return args, rows, payloads, opener, requests


def test_exact_frozen_inventory_and_prior_hashes():
    rows = module.artifacts()
    assert len(rows) == 7 and sum(row["bytes"] for row in rows) == 2753786
    assert rows[0]["sha256"] == module.mafft.SHA and rows[0]["bytes"] == 758305
    assert module.fasttree.REVISION == "29c5e62fbcd93230ee325f9c6a17b81f00e3c72a"
    assert {Path(row["relative"]).name: row["sha256"] for row in rows[1:]} == module.fasttree.FILES
    assert {row["role"] for row in rows} == {"mafft_source", "fasttree_binary", "fasttree_source_notice"}
    assert len({row["relative"] for row in rows}) == 7
    for row in rows:
        module.base.provider_url(row["url"], module.PROVIDERS)


def test_complete_without_installed_tools_execution_or_extraction(acquisition):
    args, rows, _, _, requests = acquisition
    result = module.acquire(args)
    assert result["status"] == "frozen_phylogeny_artifacts_acquired"
    assert result["artifact_files"] == len(requests) == 7
    assert result["artifact_bytes"] == sum(row["bytes"] for row in rows)
    assert all(row[1:] == (60, "identity") for row in requests)
    assert result == json.loads((args.output / "complete.json").read_bytes())
    assert result["attempts"] == 1
    for key in ("retry", "installed_tools_required", "extraction_performed", "compilation_performed",
                "installation_performed", "native_code_executed", "scientific_inference_executed",
                "controlled_timing", "publication_ready", "security_clearance", "redistribution_clearance"):
        assert result[key] is False
    for row in result["downloads"]:
        path = Path(row["file"]["path"])
        assert path.stat().st_mode & 0o777 == 0o644
        assert module.base.record(path) == row["file"]
    assert not list(args.output.rglob("*.partial"))


@pytest.mark.parametrize("defect", ["ack", "ack_int", "timeout_zero", "timeout_bool", "exists",
                                     "canonical", "inside_source", "provider"])
def test_preflight_does_not_create_or_download(acquisition, tmp_path, monkeypatch, defect):
    args, rows, _, _, requests = acquisition
    if defect.startswith("ack"):
        args.acknowledge_historical_runtime = 1 if defect == "ack_int" else False
    elif defect.startswith("timeout"):
        args.timeout = True if defect == "timeout_bool" else 0
    elif defect == "exists":
        args.output.mkdir()
    elif defect == "canonical":
        alias = tmp_path / "alias"
        alias.symlink_to(tmp_path, target_is_directory=True)
        args.output = alias / "fresh"
    elif defect == "inside_source":
        monkeypatch.setattr(module, "__file__", str(tmp_path / "source/benchmark_tools/acquire.py"))
        args.output = tmp_path / "source/fresh"
    else:
        rows[0]["url"] = "https://evil.example/archive"
    with pytest.raises((ValueError, FileExistsError)):
        module.acquire(args)
    assert not requests
    if defect != "exists":
        assert not args.output.exists()


@pytest.mark.parametrize("url", ["http://mafft.cbrc.jp/x", "https://evil.example/x",
    "https://raw.githubusercontent.com:444/x", "https://user:pass@mafft.cbrc.jp/x",
    "https://mafft.cbrc.jp/x?token=secret", "https://mafft.cbrc.jp/x#fragment",
    "https://mafft.cbrc.jp/\nx", "https://mafft.cbrc.jp/ x", "https://pypi.org/x"])
def test_tool_redirect_policy_remains_separate_from_base(url):
    with pytest.raises(ValueError):
        module.ToolRedirect().redirect_request(Request(module.mafft.URL), None, 302, "Found", {}, url)


def test_tool_redirect_permits_pinned_hosts_not_base_default():
    for url in (module.mafft.URL, module.fasttree.BASE + "LICENSE"):
        result = module.ToolRedirect().redirect_request(Request(module.mafft.URL), None, 302, "Found", {}, url)
        assert result.full_url == url
        with pytest.raises(ValueError):
            module.base.provider_url(url)
        with pytest.raises(ValueError):
            module.base.ProviderRedirect().redirect_request(Request(module.base.PIP_METADATA), None,
                302, "Found", {}, url)


@pytest.mark.parametrize("defect", ["short", "long", "hash", "network", "final_redirect", "encoded",
                                     "changed_artifact", "changed_source", "executable"])
def test_failed_downloads_preserved_and_never_retried(acquisition, monkeypatch, defect):
    args, rows, payloads, opener, requests = acquisition
    first = rows[0]["url"]
    if defect == "short":
        payloads[first] = payloads[first][:-1]
    elif defect == "long":
        payloads[first] += b"x"
    elif defect == "hash":
        payloads[first] = b"x" * len(payloads[first])
    elif defect == "network":
        def fail(*args, **kwargs):
            raise OSError("Synthetic network failure")
        monkeypatch.setattr(opener, "open", fail)
    elif defect in {"final_redirect", "encoded"}:
        def changed_response(request, timeout):
            requests.append((request.full_url, timeout, request.get_header("Accept-encoding")))
            return Response(payloads[first], "https://evil.example/x" if defect == "final_redirect" else first,
                            {"Content-Encoding": "gzip"} if defect == "encoded" else {})
        monkeypatch.setattr(opener, "open", changed_response)
    elif defect == "changed_source":
        record = module.base.record
        source = str(Path(module.__file__).absolute())
        def changed_record(path):
            result = record(path)
            if requests and str(path) == source:
                result["sha256"] = "0" * 64
            return result
        monkeypatch.setattr(module.base, "record", changed_record)
    else:
        download = module.base.download
        def changed_download(*values, **kwargs):
            result = download(*values, **kwargs)
            path = Path(result["file"]["path"])
            if defect == "executable":
                path.chmod(0o755)
            else:
                path.write_bytes(b"Changed after download")
            return result
        monkeypatch.setattr(module.base, "download", changed_download)
    with pytest.raises((ValueError, OSError)):
        module.acquire(args)
    failed = json.loads((args.output / "failed.json").read_bytes())
    assert failed["attempts"] == 1 and failed["retry"] is False
    assert failed["publication_ready"] is False and failed["native_code_executed"] is False
    assert not (args.output / "complete.json").exists()
    assert len({request[0] for request in requests}) == len(requests)
    with pytest.raises(FileExistsError):
        module.acquire(args)


def test_failure_in_later_download_retains_earlier_artifacts(acquisition, monkeypatch):
    args, _, _, opener, requests = acquisition
    original = opener.open
    def fail_second(request, timeout):
        if requests:
            raise OSError("Second download failure")
        return original(request, timeout)
    monkeypatch.setattr(opener, "open", fail_second)
    with pytest.raises(OSError):
        module.acquire(args)
    failed = json.loads((args.output / "failed.json").read_bytes())
    assert len(failed["downloads"]) == 1
    assert Path(failed["downloads"][0]["file"]["path"]).is_file()


def test_scoped_downloader_checks_final_host_without_broadening_default(tmp_path):
    class Opener:
        def open(self, request, timeout):
            return Response(b"x", module.fasttree.BASE + "LICENSE")
    expected = dict(bytes=1, sha256=module.base.hashlib.sha256(b"x").hexdigest())
    with pytest.raises(ValueError):
        module.base.download(Opener(), module.mafft.URL, tmp_path / "default", 1, 60, expected)
    result = module.base.download(Opener(), module.mafft.URL, tmp_path / "scoped", 1, 60,
                                  expected, hosts=module.PROVIDERS)
    assert result["file"]["sha256"] == expected["sha256"]
