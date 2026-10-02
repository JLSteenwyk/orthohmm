import copy
import gzip
import io
import json
from pathlib import Path
import shutil
import tarfile

import pytest

from benchmark_tools import archive_swiss_raw_sources as module
from benchmark_tools.prepare_ob_candidate_neighborhood import record
from benchmark_tools.relocate_swiss_raw_sources import stage


@pytest.fixture
def panel(tmp_path):
    original = tmp_path / "original"
    original.mkdir()
    a, b = original / "a", original / "b"
    a.write_bytes(b"first raw input")
    b.write_bytes(b"second raw input")
    source = original / "source.json"
    source.write_text(json.dumps({"records": [record(a), record(b), record(a)]}))
    binding = stage(source, record(source)["sha256"], "records", tmp_path / "staged")
    archive = tmp_path / "raw.tar.gz"
    result = module.archive(Path(binding["path"]), binding["sha256"], archive)
    return original, Path(binding["path"]), binding["sha256"], archive, result


def test_deterministic_streamed_archive_and_restore_without_originals(panel, tmp_path):
    original, binding, digest, archive, result = panel
    second = tmp_path / "second.tar.gz"
    assert module.archive(binding, digest, second)["archive"]["sha256"] == result["archive"]["sha256"]
    assert archive.read_bytes() == second.read_bytes()
    assert result["members"] == 3
    assert result["redistribution_authorized"] is False
    document = json.loads(binding.read_text())
    shutil.rmtree(original)
    shutil.rmtree(binding.parent)
    output = tmp_path / "restored"
    restored = module.restore(archive, result["archive"]["sha256"], digest, output)
    assert json.loads((output / "bindings.json").read_text()) == document
    assert len(document["entries"]) == 3
    assert len(list((output / "inputs").iterdir())) == 2
    assert restored["decompressed_bytes"] == 10240
    assert restored["members"] == 3
    assert restored["publication_ready"] is False
    assert restored["native_annotation_admission_rerun"] is False


@pytest.mark.parametrize("operation", ["archive", "restore"])
def test_no_overwrite_empty_directory_or_existing_file(panel, tmp_path, operation):
    _, binding, digest, archive, result = panel
    output = tmp_path / "exists"
    output.mkdir()
    with pytest.raises(FileExistsError):
        if operation == "archive":
            module.archive(binding, digest, output)
        else:
            module.restore(archive, result["archive"]["sha256"], digest, output)


def test_archive_cannot_mutate_input_directory(panel):
    _, binding, digest, _, _ = panel
    with pytest.raises(ValueError, match="mutate"):
        module.archive(binding, digest, binding.parent / "archive.tar.gz")


def test_unexpected_source_member_rejected_before_archive(panel, tmp_path):
    _, binding, digest, _, _ = panel
    (binding.parent / "extra").write_bytes(b"extra")
    output = tmp_path / "rejected.tar.gz"
    with pytest.raises(ValueError, match="filesystem member"):
        module.archive(binding, digest, output)
    assert not output.exists()


@pytest.mark.parametrize("which", ["archive", "binding"])
def test_independent_digests_not_replaced(panel, tmp_path, which):
    _, _, digest, archive, result = panel
    output = tmp_path / "rejected"
    with pytest.raises(ValueError, match="changed"):
        module.restore(archive, "0" * 64 if which == "archive" else result["archive"]["sha256"],
                       "0" * 64 if which == "binding" else digest, output)
    if which == "archive":
        assert not output.exists()


@pytest.mark.parametrize("fault", ["extra", "duplicate", "missing", "traversal", "absolute", "symlink",
    "hardlink", "directory", "mode", "pax", "wrong_size", "changed_bytes", "binding_last"])
def test_rehashed_archive_cannot_bypass_binding_or_member_guards(panel, tmp_path, fault):
    _, _, digest, archive, _ = panel
    with tarfile.open(archive, "r:gz") as reader:
        members = [(copy.copy(m), reader.extractfile(m).read()) for m in reader]
    info, data = members[1]
    if fault == "missing":
        members.pop()
    elif fault == "duplicate":
        members.append((copy.copy(info), data))
    elif fault == "extra":
        extra = copy.copy(info)
        extra.name = "extra"
        members.append((extra, data))
    elif fault == "binding_last":
        members.reverse()
    elif fault in ("traversal", "absolute"):
        info.name = "../escape" if fault == "traversal" else str(tmp_path / "escape")
    elif fault in ("symlink", "hardlink", "directory"):
        info.type = {"symlink": tarfile.SYMTYPE, "hardlink": tarfile.LNKTYPE,
                     "directory": tarfile.DIRTYPE}[fault]
        info.linkname = "../escape"
    elif fault == "mode":
        info.mode = 0o777
    elif fault == "pax":
        info.pax_headers = {"comment": "unexpected"}
    elif fault == "wrong_size":
        info.size += 1
        members[1] = info, data + b"X"
    else:
        members[1] = info, b"X" + data[1:]
    changed = tmp_path / "changed.tar.gz"
    with tarfile.open(changed, "w:gz", format=tarfile.PAX_FORMAT if fault == "pax" else tarfile.USTAR_FORMAT) as writer:
        for item, payload in members:
            writer.addfile(item, io.BytesIO(payload))
    output = tmp_path / "rejected"
    with pytest.raises(ValueError):
        module.restore(changed, record(changed)["sha256"], digest, output)
    assert not (tmp_path / "escape").exists()


@pytest.mark.parametrize("fault", ["truncated", "crc", "extra_padding", "padding_bomb"])
def test_gzip_footer_and_actual_stream_budget_checked(panel, tmp_path, monkeypatch, fault):
    _, _, digest, archive, result = panel
    raw = archive.read_bytes()
    if fault == "truncated":
        raw = raw[:-8]
    elif fault == "crc":
        raw = raw[:-8] + bytes([raw[-8] ^ 1]) + raw[-7:]
    else:
        raw += gzip.compress(b"\0" * (100 if fault == "extra_padding" else 30000))
        monkeypatch.setattr(module, "MAX_BYTES", 20480)
    changed = tmp_path / "changed.tar.gz"
    changed.write_bytes(raw)
    with pytest.raises((ValueError, EOFError, OSError, tarfile.TarError)):
        module.restore(changed, record(changed)["sha256"], digest, tmp_path / "rejected")


def test_compressed_size_guard_before_destination_creation(panel, tmp_path, monkeypatch):
    _, _, digest, archive, result = panel
    monkeypatch.setattr(module, "MAX_BYTES", 1)
    output = tmp_path / "rejected"
    with pytest.raises(ValueError, match="Compressed archive"):
        module.restore(archive, result["archive"]["sha256"], digest, output)
    assert not output.exists()


def test_tar_padding_budget_checked_before_archive_creation(panel, tmp_path, monkeypatch):
    _, binding, digest, _, _ = panel
    monkeypatch.setattr(module, "MAX_BYTES", 9000)
    output = tmp_path / "rejected.tar.gz"
    with pytest.raises(ValueError, match="padding"):
        module.archive(binding, digest, output)
    assert not output.exists()


def test_source_mutation_retains_failed_archive_attempt(panel, tmp_path, monkeypatch):
    _, binding, digest, _, _ = panel
    real = module.inventory
    calls = []
    def mutate_after_check(*args):
        files = real(*args)
        calls.append(True)
        if len(calls) == 1:
            path = next((binding.parent / "inputs").iterdir())
            path.write_bytes(b"X" + path.read_bytes()[1:])
        return files
    monkeypatch.setattr(module, "inventory", mutate_after_check)
    output = tmp_path / "failed.tar.gz"
    with pytest.raises(ValueError, match="Raw payload changed while copying"):
        module.archive(binding, digest, output)
    assert output.exists()


def test_transient_copy_mutation_rejected_even_when_source_checks_match(panel, tmp_path, monkeypatch):
    _, binding, digest, _, _ = panel
    before = module.inventory(binding, digest)
    real_copy = tarfile.copyfileobj
    def transient(src, dst, *args, **kwargs):
        path = Path(getattr(src, "stream", src).name)
        if path.name == "bindings.json":
            return real_copy(src, dst, *args, **kwargs)
        original = path.read_bytes()
        path.write_bytes(b"X" + original[1:])
        try:
            return real_copy(src, dst, *args, **kwargs)
        finally:
            path.write_bytes(original)
    monkeypatch.setattr(module.tarfile, "copyfileobj", transient)
    with pytest.raises(ValueError, match="copying"):
        module.archive(binding, digest, tmp_path / "transient.tar.gz")
    assert module.inventory(binding, digest) == before


@pytest.mark.parametrize("fault", ["rights", "schema", "artifact", "conflict", "budget"])
def test_changed_binding_schema_and_declared_budget_rejected(panel, fault, monkeypatch):
    _, binding, _, _, _ = panel
    doc = json.loads(binding.read_text())
    if fault == "rights":
        doc["redistribution_authorized"] = True
    elif fault == "schema":
        doc["schema_version"] = True
    elif fault == "artifact":
        doc["entries"][0]["artifact"] = "../outside"
    elif fault == "conflict":
        doc["entries"][-1]["original"]["bytes"] += 1
    else:
        monkeypatch.setattr(module, "MAX_RECORDS", 2)
    with pytest.raises(ValueError):
        module.declared_files(doc)
