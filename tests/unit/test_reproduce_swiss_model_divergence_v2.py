"""Versioned archive-header repair; original source and numerical functions preserved."""

import hashlib
import io
import json
from pathlib import Path
import tarfile

import pytest

from benchmark_tools import reproduce_swiss_model_divergence_v2 as replay

ROOT = Path(__file__).resolve().parents[2]


def test_v2_source_is_exactly_the_declared_four_changes():
    receipt = json.loads((ROOT / "benchmark_tools/results/swiss_model_divergence_portable_source_revision_20261007_v2.json").read_text())
    parent = (ROOT / receipt["parent"]).read_bytes()
    assert hashlib.sha256(parent).hexdigest() == receipt["parent_sha256"]
    expected = parent.decode()
    assert len(receipt["edits"]) == 4
    for edit in receipt["edits"]:
        assert expected.count(edit["from"]) == 1
        expected = expected.replace(edit["from"], edit["to"], 1)
    assert expected == (ROOT / receipt["output"]).read_text()
    assert receipt["numerical_functions_changed"] is False
    assert receipt["original_v1_or_evidence_archive_changed"] is False
    assert replay.SCHEMA == "swiss_model_divergence_portable_v2"


def fixture(tmp_path, corruption=None):
    names = [f"invented/{i:03d}.txt" for i in range(183)]
    refs = [dict(path="/old/benchmarks/results/" + name, **replay.identity(name.encode())) for name in names]
    archive = tmp_path / "evidence.tar.gz"
    with tarfile.open(archive, "w:gz") as container:
        header = tarfile.TarInfo("invented")
        header.type, header.mode = tarfile.DIRTYPE, 0o775
        if corruption == "unsafe_directory":
            header.name = "../escape"
        elif corruption == "unknown_directory":
            header.name = "unknown"
        elif corruption == "symlink":
            header.type, header.linkname = tarfile.SYMTYPE, "/outside"
        elif corruption == "hardlink":
            header.type, header.linkname = tarfile.LNKTYPE, "/outside"
        if corruption == "nonzero_directory":
            header.size = 1
            container.addfile(header, io.BytesIO(b"x"))
        else:
            container.addfile(header)
        if corruption == "duplicate_directory":
            container.addfile(header)
        for i, name in enumerate(names):
            if corruption == "missing_file" and i == 182:
                continue
            member = tarfile.TarInfo(names[0] if corruption == "duplicate_file" and i == 182 else name)
            raw = b"changed" if corruption == "checksum" and i == 0 else name.encode()
            member.size = len(raw)
            container.addfile(member, io.BytesIO(raw))
    return archive, dict(checked_inputs=refs)


def test_safe_ancestor_headers_preserve_all_regular_files(tmp_path):
    archive, features = fixture(tmp_path)
    output = tmp_path / "restored"
    assert replay.restore_payloads(archive, features, output) == 183
    assert (output / "invented/000.txt").read_bytes() == b"invented/000.txt"
    assert len(list(output.rglob("*.txt"))) == 183


@pytest.mark.parametrize("corruption", ["unsafe_directory", "unknown_directory", "symlink", "hardlink",
    "nonzero_directory", "duplicate_directory", "missing_file", "duplicate_file", "checksum"])
def test_repair_keeps_safety_and_completeness_gates(tmp_path, corruption):
    archive, features = fixture(tmp_path, corruption)
    output = tmp_path / "restored"
    with pytest.raises(ValueError):
        replay.restore_payloads(archive, features, output)
    assert not output.exists() and not (tmp_path / "escape").exists()


def test_retained_actual_archive_shape_explains_v1_failure_without_replaying_it():
    path = ROOT / "benchmark_tools/results/swiss_model_divergence_evidence_23932_v1.tar.gz"
    assert hashlib.sha256(path.read_bytes()).hexdigest() == "eeed6e6786f5720ee8c19fb6123fad64c38fdb25c33afe57f5e61bcae4a6113d"
    with tarfile.open(path, "r:gz") as archive:
        members = archive.getmembers()
        assert len(members) == 202
        assert sum(m.isfile() for m in members) == 183
        assert sum(m.isdir() for m in members) == 19
        assert all(m.isfile() or m.isdir() for m in members)
    failure = json.loads((ROOT / "benchmark_tools/results/swiss_model_divergence_portable_failed_20261007_v1.json").read_text())
    assert failure["status"] == "portable_replay_failed" and failure["retry"] is False
    assert failure["error"] == "Incorrect archive inventory or types"
