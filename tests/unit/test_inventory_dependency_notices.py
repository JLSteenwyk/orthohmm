from copy import deepcopy
import hashlib
import json
import zipfile

import pytest

from benchmark_tools import inventory_dependency_notices as module


def fixture(tmp_path, metadata_extra="", extra=None):
    wheel = tmp_path / "example-1.0-py3-none-any.whl"
    with zipfile.ZipFile(wheel, "w") as archive:
        archive.writestr("example-1.0.dist-info/METADATA",
            "Metadata-Version: 2.4\nName: example\nVersion: 1.0\n" + metadata_extra + "\n")
        for name, content in (extra or {}).items():
            archive.writestr(name, content)
    return dict(download_info=dict(url=wheel.as_uri(), archive_info=dict(
        hashes=dict(sha256=hashlib.sha256(wheel.read_bytes()).hexdigest()))),
        metadata=dict(name="example", version="1.0"))


def test_declared_nested_and_heuristic_notices(tmp_path):
    item = fixture(tmp_path, "License-Expression: MIT\nLicense-File: COPYING\n", {
        "example-1.0.dist-info/licenses/COPYING": "terms", "example/NOTICE.txt": "notice",
        "example/libexample.so.1": b"not executed", "example/COPYRIGHT": "copyright"})
    row = module.inventory(item)
    assert row["license_expression"] == "MIT"
    assert len(row["notice_candidates"]) == 3
    assert row["unresolved_declarations"] == []
    assert row["native_members"] == ["example/libexample.so.1"]
    assert row["redistribution_clearance"] is False


def test_legacy_license_location(tmp_path):
    item = fixture(tmp_path, "License: BSD\nLicense-File: LICENSE.txt\n",
        {"example-1.0.dist-info/LICENSE.txt": "terms"})
    digest = item["download_info"]["archive_info"].pop("hashes")["sha256"]
    item["download_info"]["archive_info"]["hash"] = "sha256=" + digest
    row = module.inventory(item)
    assert row["license_field_short"] == "BSD"
    assert row["unresolved_declarations"] == []


def test_vendored_metadata_does_not_replace_distribution(tmp_path):
    row = module.inventory(fixture(tmp_path, extra={
        "example/_vendor/other-2.0.dist-info/METADATA": "Name: other\nVersion: 2.0\n",
        "example/_vendor/other-2.0.dist-info/licenses/LICENSE": "vendored terms"}))
    assert row["name"] == "example"
    assert len(row["notice_candidates"]) == 1


def test_multiple_top_level_metadata_rejected(tmp_path):
    with pytest.raises(ValueError, match="metadata"):
        module.inventory(fixture(tmp_path, extra={"other-1.0.dist-info/METADATA": "Name: other\n"}))


@pytest.mark.parametrize("extra,count", [({}, 0), ({
    "example-1.0.dist-info/COPYING": "a", "example-1.0.dist-info/licenses/COPYING": "b"}, 2)])
def test_missing_or_ambiguous_is_not_clearance(tmp_path, extra, count):
    row = module.inventory(fixture(tmp_path, "License-File: COPYING\n", extra))
    assert len(row["unresolved_declarations"]) == 1
    assert len(row["unresolved_declarations"][0]["matching_members"]) == count
    assert row["redistribution_clearance"] is False


@pytest.mark.parametrize("value", ["../LICENSE", "/LICENSE", "..\\LICENSE"])
def test_unsafe_declaration_rejected(tmp_path, value):
    with pytest.raises(ValueError):
        module.inventory(fixture(tmp_path, "License-File: " + value + "\n"))


def test_long_license_not_silently_truncated(tmp_path):
    text = "terms " * 100
    row = module.inventory(fixture(tmp_path, "License: " + text + "\n"))
    assert row["license_field_short"] is None
    assert row["license_field_sha256"] == hashlib.sha256(text.encode()).hexdigest()


@pytest.mark.parametrize("kind", ["hash", "legacy_conflict", "version", "name", "remote_url"])
def test_mismatched_source_rejected(tmp_path, kind):
    item = fixture(tmp_path)
    if kind == "hash":
        item["download_info"]["archive_info"]["hashes"]["sha256"] = "0" * 64
    elif kind == "legacy_conflict":
        item["download_info"]["archive_info"]["hash"] = "sha256=" + "0" * 64
    elif kind == "remote_url":
        item["download_info"]["url"] = "https://example.org/example.whl"
    else:
        item["metadata"][kind] = "2.0" if kind == "version" else "other"
    with pytest.raises(ValueError):
        module.inventory(item)


def test_report_pin_and_cli_no_overwrite(tmp_path):
    item = fixture(tmp_path)
    path, output = tmp_path / "install.json", tmp_path / "out.json"
    path.write_text(json.dumps(dict(install=[item])))
    sha = hashlib.sha256(path.read_bytes()).hexdigest()
    argv = ["--install-report", str(path), "--report-sha", sha, "--output", str(output)]
    assert module.main(argv) == 0
    assert json.loads(output.read_text())["publication_ready"] is False
    with pytest.raises(FileExistsError):
        module.main(argv)
    with pytest.raises(ValueError, match="identity"):
        module.audit(path, "0" * 64)


def test_duplicate_distributions_rejected(tmp_path):
    item = fixture(tmp_path)
    path = tmp_path / "install.json"
    path.write_text(json.dumps(dict(install=[item, deepcopy(item)])))
    with pytest.raises(ValueError, match="unique"):
        module.audit(path, hashlib.sha256(path.read_bytes()).hexdigest())
