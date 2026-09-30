import io
import json
import tarfile
from urllib.error import HTTPError

import pytest

from benchmark_tools.inspect_qfo_darwin_images import cached_layer, context_layers, digest, inspect, inventory, verified
from benchmark_tools import inspect_qfo_darwin_images as module


def test_digest_and_size_both_required():
    data = b"public blob"
    descriptor = {"size": len(data), "digest": digest(data)}
    assert verified(data, descriptor) == data
    with pytest.raises(ValueError):
        verified(data, {**descriptor, "size": len(data) + 1})
    with pytest.raises(ValueError):
        verified(data + b"!", {**descriptor, "size": len(data) + 1})


def test_history_skips_empty_entries_and_maps_copy_only():
    layers = [{"digest": "base"}, {"digest": "context"}, {"digest": "other"}]
    history = [{"created_by": "ADD ubuntu.tar /"},
               {"created_by": "WORKDIR /benchmark", "empty_layer": True},
               {"created_by": "COPY . /benchmark # buildkit"},
               {"created_by": "RUN install"}]
    assert context_layers({"layers": layers}, {"history": history}) == [
        (1, layers[1], "COPY . /benchmark # buildkit")]


def test_bad_history_does_not_imply_no_originals():
    with pytest.raises(ValueError):
        context_layers({"layers": [{}]}, {"history": []})


def test_tar_is_inventory_only_even_with_unsafe_paths_and_links():
    output = io.BytesIO()
    with tarfile.open(fileobj=output, mode="w:gz") as archive:
        for name in ["benchmark/data/treefam/TF0001.nhx.gz",
                     "../treefam2reference.txt", "benchmark/TreeFam-A.json"]:
            member = tarfile.TarInfo(name)
            member.size = 3
            archive.addfile(member, io.BytesIO(b"abc"))
        link = tarfile.TarInfo("benchmark/external")
        link.type = tarfile.SYMTYPE
        link.linkname = "/etc/passwd"
        archive.addfile(link)
    result = inventory(output.getvalue())
    assert result["member_count"] == 4
    assert len(result["original_filename_candidates"]) == 2
    assert len(result["treefam_named_members"]) == 3
    assert result["members"][-1]["linkname"] == "/etc/passwd"


def test_verified_cache_lookup(tmp_path):
    data = b"previous public blob"
    descriptor = {"size": len(data), "digest": digest(data)}
    filename = descriptor["digest"].split(":")[1] + ".tar.gz"
    missing = tmp_path / "missing"
    assert cached_layer(descriptor, [missing]) == (None, None)
    path = tmp_path / filename
    path.write_bytes(data)
    found, ref = cached_layer(descriptor, [missing, tmp_path])
    assert found == data
    assert ref == {"path": str(path), "bytes": len(data), "digest": digest(data)}
    path.write_bytes(b"!" * len(data))
    with pytest.raises(ValueError):
        cached_layer(descriptor, [tmp_path])


def test_cache_rejects_bad_size_and_symlink(tmp_path):
    data = b"blob"
    descriptor = {"size": len(data), "digest": digest(data)}
    path = tmp_path / (descriptor["digest"].split(":")[1] + ".tar.gz")
    path.write_bytes(data + b"!")
    with pytest.raises(ValueError):
        cached_layer(descriptor, [tmp_path])
    path.unlink()
    source = tmp_path / "source"
    source.write_bytes(data)
    path.symlink_to(source)
    with pytest.raises(ValueError):
        cached_layer(descriptor, [tmp_path])


@pytest.mark.parametrize("descriptor", [
    {"digest": "sha256:../../secret", "size": 1},
    {"digest": "sha256:" + "0" * 64, "size": True},
    {"digest": "sha256:" + "0" * 64, "size": -1},
    {"digest": "sha256:" + "0" * 64, "size": 32_000_001},
])
def test_invalid_descriptor_rejected_without_read(descriptor):
    with pytest.raises(ValueError):
        cached_layer(descriptor, [])


@pytest.mark.parametrize("tags", [[], ["2020.1", "2020.1"], ["../outside"], ["tag/other"], [""]])
def test_unsafe_tags_rejected_before_creating_directory(tmp_path, tags):
    output = tmp_path / "output"
    with pytest.raises(ValueError):
        inspect(output, tags)
    assert not output.exists()


class Response(io.BytesIO):
    status = 200


def registry_fixture(monkeypatch, blob, failure=None):
    descriptor = {"size": len(blob), "digest": digest(blob)}
    config = json.dumps({"history": [{"created_by": "COPY . /benchmark"}]}).encode()
    config_ref = {"size": len(config), "digest": digest(config)}
    manifest = json.dumps({"layers": [descriptor], "config": config_ref}).encode()
    hits = []

    class Opener:
        def open(self, request, timeout):
            url = request if isinstance(request, str) else request.full_url
            hits.append(url)
            if url == module.TAGS_URL:
                return Response(json.dumps({"count": 2, "next": None}).encode())
            if url.startswith("https://auth.docker.io/"):
                return Response(b'{"token":"ephemeral-test-token"}')
            if "/manifests/" in url:
                if failure is not None:
                    raise failure
                return Response(manifest)
            if url.endswith(config_ref["digest"]):
                return Response(config)
            if url.endswith(descriptor["digest"]):
                return Response(blob)
            raise AssertionError(url)

    monkeypatch.setattr(module, "build_opener", lambda *_: Opener())
    return descriptor, hits


def test_context_cache_reused_and_reported_without_network_pull(tmp_path, monkeypatch):
    buffer = io.BytesIO()
    with tarfile.open(fileobj=buffer, mode="w:gz") as archive:
        member = tarfile.TarInfo("benchmark/data/TF0001.nhx")
        member.size = 3
        archive.addfile(member, io.BytesIO(b"abc"))
    blob = buffer.getvalue()
    descriptor, hits = registry_fixture(monkeypatch, blob)
    cache = tmp_path / "cache"
    cache.mkdir()
    (cache / (descriptor["digest"].split(":")[1] + ".tar.gz")).write_bytes(blob)
    result = inspect(tmp_path / "output", ["old-a", "old-b"], [cache])
    assert result["downloaded_unique_context_bytes"] == 0
    assert len(result["reused_contexts"]) == 1
    assert all(r["status"] == "selected_context_layers_inspected" for r in result["results"])
    assert all(len(r["layers"][0]["original_filename_candidates"]) == 1 for r in result["results"])
    assert not any(url.endswith(descriptor["digest"]) for url in hits)
    assert "ephemeral-test-token" not in (tmp_path / "output" / "report.json").read_text()


def test_registry_rate_limit_stops_without_retry(tmp_path, monkeypatch):
    _, hits = registry_fixture(monkeypatch, b"unused", HTTPError("url", 429, "Rate limited", {}, None))
    result = inspect(tmp_path / "output", ["old-a", "old-b"])
    assert [r["status"] for r in result["results"]] == ["unresolved", "not_attempted_rate_limit"]
    assert sum("/manifests/" in url for url in hits) == 1


def test_network_timeout_retained_as_unresolved(tmp_path, monkeypatch):
    registry_fixture(monkeypatch, b"unused", TimeoutError("bounded observation timeout"))
    result = inspect(tmp_path / "output", ["old-a"])
    assert result["results"][0]["status"] == "unresolved"
    assert "timeout" in result["results"][0]["error"]
