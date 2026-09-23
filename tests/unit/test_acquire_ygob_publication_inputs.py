import hashlib
import io
import json

import pytest

from benchmark_tools import acquire_ygob_publication_inputs as module


@pytest.mark.parametrize("fault", [None, "short", "long", "digest", "network"])
def test_download_is_verified_or_retained_as_failure(tmp_path, monkeypatch, fault):
    payload = b"test input"
    monkeypatch.setattr(module, "FILES", {"AA.fsa": (len(payload), hashlib.sha256(payload).hexdigest())})
    class Response(io.BytesIO):
        status = 200
        def geturl(self): return module.BASE + "AA.fsa"
    def open_url(url, timeout):
        assert timeout == 30 and url == module.BASE + "AA.fsa"
        if fault == "network": raise OSError("fixture")
        data = {"short": payload[:-1], "long": payload+b"x", "digest": b"X"*len(payload)}.get(fault, payload)
        return Response(data)
    monkeypatch.setattr(module, "urlopen", open_url)
    output = tmp_path / "download"
    if fault:
        with pytest.raises((ValueError, OSError)): module.acquire(output)
        assert not (output / "AA.fsa").exists()
        assert json.loads((output / "acquisition.json").read_text())["status"] == "acquisition_failed"
        if fault != "network": assert (output / "AA.fsa.partial").exists()
    else:
        result = module.acquire(output)
        assert result["status"] == "frozen_ygob_inputs_acquired"
        assert result["redistribution_cleared"] is False
        assert (output / "AA.fsa").read_bytes() == payload
    with pytest.raises(FileExistsError): module.acquire(output)


def test_constants_match_retained_candidate_audit():
    from pathlib import Path
    root = Path(__file__).resolve().parents[2]
    report = json.loads((root / "benchmark_tools/results/ygob_overlap_20260916.json").read_text())
    assert module.FILES == {Path(r["path"]).name: (r["bytes"], r["sha256"]) for r in report["candidate_inputs"]}
