import hashlib
import zipfile

import pytest

from benchmark_tools import reconstruct_publication_project_wheel as module


@pytest.fixture
def inputs(tmp_path, monkeypatch):
    old, candidate = tmp_path / "old" / module.NAME, tmp_path / "candidate" / module.NAME
    old.parent.mkdir(); candidate.parent.mkdir()
    rows = []
    payloads = {"orthohmm/example.py": b"pass\n", "orthohmm-0.5.0.dist-info/METADATA": b"Name: orthohmm\n"}
    with zipfile.ZipFile(old, "w", compression=zipfile.ZIP_DEFLATED, compresslevel=6) as archive:
        for name, data in payloads.items():
            info = zipfile.ZipInfo(name, (2026, 9, 26, 20, 58, 32))
            info.create_system = 3
            info.external_attr = 0o100644 << 16
            archive.writestr(info, data, compress_type=zipfile.ZIP_DEFLATED, compresslevel=6)
            rows.append([name, list(info.date_time), info.external_attr, len(data), hashlib.sha256(data).hexdigest()])
    with zipfile.ZipFile(candidate, "w") as archive:
        for name, data in payloads.items(): archive.writestr(name, data)
        archive.comment = b"Different candidate wrapper"
    monkeypatch.setattr(module, "MEMBERS", rows)
    monkeypatch.setattr(module, "EXPECTED", {k: module.record(old)[k] for k in ("bytes", "sha256")})
    return old, candidate, tmp_path / "fresh"


def test_wrapper_reconstruction_is_exact_and_executes_no_payload(inputs):
    old, candidate, output = inputs
    before = candidate.read_bytes()
    result = module.run(candidate, output)
    assert (output / module.NAME).read_bytes() == old.read_bytes() and candidate.read_bytes() == before
    assert result["status"] == "historical_project_wheel_byte_reconstructed"
    assert result["matched_payload_members"] == 2 and result["historical_wheel_reproduced"] is True
    assert (output / module.NAME).stat().st_mode & 0o777 == 0o644
    assert all(result[k] is False for k in ("original_wheel_required", "retry", "native_code_executed",
        "historical_admission", "installation_performed", "scientific_inference_executed",
        "controlled_timing", "publication_ready", "security_clearance", "redistribution_clearance"))


@pytest.mark.parametrize("fault", ["missing", "extra", "duplicate", "size", "digest", "name", "symlink", "existing"])
def test_preflight_rejects_changes_before_output(inputs, tmp_path, fault):
    _, candidate, output = inputs
    if fault == "existing": output.mkdir()
    elif fault == "name": candidate = candidate.rename(candidate.with_name("wrong.whl"))
    elif fault == "symlink":
        link = tmp_path / module.NAME; link.symlink_to(candidate); candidate = link
    else:
        rows = module.MEMBERS
        with zipfile.ZipFile(candidate, "w") as archive:
            for i, row in enumerate(rows):
                if fault == "missing" and i == 0: continue
                data = b"pass\n" if i == 0 else b"Name: orthohmm\n"
                if i == 0 and fault == "size": data += b"x"
                if i == 0 and fault == "digest": data = b"fail\n"
                archive.writestr(row[0], data)
            if fault == "extra": archive.writestr("unknown", b"x")
            if fault == "duplicate":
                with pytest.warns(UserWarning): archive.writestr(rows[0][0], b"pass\n")
    with pytest.raises((ValueError, FileExistsError)): module.run(candidate, output)
    if fault != "existing": assert not output.exists()


@pytest.mark.parametrize("fault", ["compression", "input"])
def test_failed_attempt_is_preserved_without_retry(inputs, monkeypatch, fault):
    _, candidate, output = inputs
    original = module.pack
    calls = []
    def pack(path, payloads):
        calls.append(path)
        original(path, payloads)
        if fault == "compression": path.write_bytes(b"Wrong reconstructed wrapper")
        else: candidate.write_bytes(b"Changed supplied candidate")
    monkeypatch.setattr(module, "pack", pack)
    with pytest.raises(ValueError): module.run(candidate, output)
    assert len(calls) == 1 and (output / (module.NAME + ".partial")).exists()
    assert (output / "failed.json").exists() and not (output / "complete.json").exists()
