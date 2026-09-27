import io
import json
import tarfile

import pytest

from benchmark_tools.archive_portable_ob_score import encoded, identity, read_archive, restore, safe_name, write_archive


def payloads():
    data = {"worker.py": b"print('fixture')\n"}
    manifest = dict(schema_version=1, scope="orthobench_scoring_only", publication_ready=False,
                    raw_upstream_inputs_included=False, files={n: identity(p) for n, p in data.items()})
    return {**data, "archive_manifest.json": encoded(manifest)}


def test_deterministic_archive_roundtrip(tmp_path):
    a, b = tmp_path / "a.tar.gz", tmp_path / "b.tar.gz"
    data = payloads()
    write_archive(a, data)
    write_archive(b, data)
    assert a.read_bytes() == b.read_bytes()
    assert read_archive(a, identity(a.read_bytes())["sha256"])[0] == data


@pytest.mark.parametrize("name", ["/absolute", "../escape", "a/../b", "a//b", "a\\b", ""])
def test_unsafe_name(name):
    with pytest.raises(ValueError, match="Unsafe"):
        safe_name(name)


def test_wrong_archive_digest(tmp_path):
    path = tmp_path / "archive.tar.gz"
    write_archive(path, payloads())
    with pytest.raises(ValueError, match="checksum mismatch"):
        read_archive(path, "wrong")


def test_changed_member_digest(tmp_path):
    data = payloads()
    data["worker.py"] = b"changed"
    path = tmp_path / "archive.tar.gz"
    write_archive(path, data)
    with pytest.raises(ValueError, match="member checksum mismatch"):
        read_archive(path, identity(path.read_bytes())["sha256"])


@pytest.mark.parametrize("kind", ["symlink", "duplicate"])
def test_nonregular_or_duplicate_member(tmp_path, kind):
    path = tmp_path / "archive.tar.gz"
    with tarfile.open(path, "w:gz") as stream:
        member = tarfile.TarInfo("worker.py")
        if kind == "symlink":
            member.type = tarfile.SYMTYPE
            member.linkname = "/etc/passwd"
            stream.addfile(member)
        else:
            stream.addfile(member, io.BytesIO(b""))
            stream.addfile(member, io.BytesIO(b""))
    with pytest.raises(ValueError, match="Unexpected archive member"):
        read_archive(path, identity(path.read_bytes())["sha256"])


def test_restore_refuses_existing_directory(tmp_path):
    with pytest.raises(FileExistsError):
        restore(tmp_path / "absent", "wrong", tmp_path, tmp_path / "python", tmp_path)


def test_external_data_mismatch_stops_before_execution(tmp_path):
    acquisition = tmp_path / "acquired"
    source = acquisition / "BENCHMARKS/Input/a.fa"
    source.parent.mkdir(parents=True)
    source.write_bytes(b"changed")
    data = payloads()
    data["inputs.json"] = encoded(dict(files=[dict(relative_path="data/BENCHMARKS/Input/a.fa",
        role="fasta", **identity(b"expected"))]))
    manifest = json.loads(data["archive_manifest.json"])
    manifest["files"]["inputs.json"] = identity(data["inputs.json"])
    data["archive_manifest.json"] = encoded(manifest)
    archive = tmp_path / "archive.tar.gz"
    write_archive(archive, data)
    output = tmp_path / "output"
    with pytest.raises(ValueError, match="Acquired input differs"):
        restore(archive, identity(archive.read_bytes())["sha256"], acquisition, tmp_path / "missing-python", output)
    assert not output.exists()


def test_archive_scope_cannot_claim_publication_readiness(tmp_path):
    data = payloads()
    manifest = json.loads(data["archive_manifest.json"])
    manifest["publication_ready"] = True
    data["archive_manifest.json"] = encoded(manifest)
    path = tmp_path / "archive.tar.gz"
    write_archive(path, data)
    with pytest.raises(ValueError, match="Unexpected archive scope"):
        read_archive(path, identity(path.read_bytes())["sha256"])
