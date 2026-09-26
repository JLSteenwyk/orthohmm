import hashlib
import zipfile

import pytest

from benchmark_tools.audit_frozen_overlay_install import scientific_members, local_install_wheels


@pytest.mark.parametrize("change", [None, "source", "missing", "extra"])
def test_wheel_scientific_bytes_and_inventory(tmp_path, change):
    files = {"orthohmm/__init__.py": {"content": b"# frozen\n", "git_blob": "blob"}}
    wheel = tmp_path / "fixture.whl"
    with zipfile.ZipFile(wheel, "w") as archive:
        if change != "missing":
            archive.writestr("orthohmm/__init__.py", b"changed" if change == "source" else b"# frozen\n")
        for name in ("hmm_viterbi.so", "kmer_prefilter.so", "pair_align.so"):
            archive.writestr("orthohmm/search/csrc/" + name, b"fixture library")
        if change == "extra":
            archive.writestr("orthohmm/extra.py", b"extra")
    if change:
        with pytest.raises(ValueError):
            scientific_members(wheel, files)
    else:
        result = scientific_members(wheel, files)
        assert len(result) == 1 and result[0]["git_blob"] == "blob"


@pytest.mark.parametrize("change", [None, "remote", "wrong_hash", "extra_wheel", "duplicate"])
def test_install_report_binds_local_wheelhouse(tmp_path, change):
    wheel = tmp_path / "package.whl"
    wheel.write_bytes(b"fixture")
    item = dict(metadata={"name": "package", "version": "1"}, download_info={
        "url": wheel.as_uri(), "archive_info": {"hashes": {"sha256": hashlib.sha256(b"fixture").hexdigest()}}})
    report = {"install": [item]}
    if change == "remote":
        item["download_info"]["url"] = "https://example.invalid/package.whl"
    elif change == "wrong_hash":
        item["download_info"]["archive_info"]["hashes"]["sha256"] = "0" * 64
    elif change == "extra_wheel":
        (tmp_path / "extra.whl").write_bytes(b"other")
    elif change == "duplicate":
        report["install"].append(item)
    if change:
        with pytest.raises(ValueError):
            local_install_wheels(report, tmp_path)
    else:
        assert local_install_wheels(report, tmp_path)[0]["wheel"]["path"] == str(wheel)
