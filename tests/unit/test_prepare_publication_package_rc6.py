"""rc6 preserves immutable history and separates reporting from scientific admission."""

import copy
import json
from pathlib import Path

import pytest

from benchmark_tools import prepare_publication_package_rc6 as current
from benchmark_tools.bundle_publication_package import archive, build, identity, restore


ROOT = Path(__file__).resolve().parents[2]


def fixture(tmp_path):
    parent = tmp_path / "parent"
    parent.mkdir()
    rows, files = [], []
    for name in ("README.md", "data.json", "PACKAGE_SELECTION.json", "bundle_publication_package.py"):
        path = parent / name
        path.write_text(name)
        rows.append(dict(path=name, mode=0o644, **identity(path)))
        if name in ("README.md", "data.json"):
            original = tmp_path / ("mutable-" + name)
            original.write_text("drifted original bytes")
            files.append(dict(target=name, source=original.name, **identity(path)))
    baseline = dict(schema="publication_package_selection_v1", version="orthohmm-study-2026.10.06-rc5",
        scientific_revision="a" * 40, workflow_revision="b" * 40, files=files,
        publication_ready=False, public_archive_uploaded=False, limitations=["retained limit"])
    index = dict(schema="publication_package_v1", version=baseline["version"],
        scientific_revision=baseline["scientific_revision"], workflow_revision=baseline["workflow_revision"],
        files=rows, publication_ready=False, public_archive_uploaded=False, native_inference_repeated=False)
    (parent / "PACKAGE_INDEX.json").write_text(json.dumps(index))
    (tmp_path / "new-guide.md").write_text("new guide")
    return baseline, index, parent, [("README.md", "new-guide.md", None)]


def test_immutable_inheritance_preserves_every_payload_and_source_revision(tmp_path):
    baseline, index, parent, additions = fixture(tmp_path)
    before = copy.deepcopy(baseline)
    result = current.compose(tmp_path, baseline, index, parent, "c" * 40, additions)
    rows = {row["target"]: row for row in result["files"]}
    assert baseline == before
    assert rows["data.json"]["source"] == "parent/data.json"
    for name in ("README.md", "PACKAGE_SELECTION.json", "bundle_publication_package.py"):
        assert rows["history/rc5/" + name]["source"] == "parent/" + name
        assert rows["history/rc5/" + name]["sha256"] == identity(parent / name)["sha256"]
    assert "history/rc5/PACKAGE_INDEX.json" in rows
    assert result["preserved_parent_selected_payloads"] == 2
    assert result["preserved_parent_indexed_payloads"] == result["preserved_parent_payloads"] == 4
    assert result["version"] == current.VERSION
    assert result["scientific_revision"] == "a" * 40
    assert result["publication_ready"] is result["public_archive_uploaded"] is False
    assert result["limitations"][0] == "retained limit"


@pytest.mark.parametrize("change", [
    "payload", "index", "duplicate_index", "extra_index", "scope", "workflow", "version",
    "mode", "duplicate_target", "escape", "revision", "missing", "symlink", "component",
])
def test_changed_or_unsafe_selection_refused(tmp_path, change):
    baseline, index, parent, additions = fixture(tmp_path)
    revision = "c" * 40
    if change == "payload": (parent / "data.json").write_text("tampered")
    elif change == "index": index["files"][0]["bytes"] += 1
    elif change == "duplicate_index": index["files"].append(index["files"][0])
    elif change == "extra_index": index["files"].append(dict(index["files"][0], path="extra"))
    elif change == "scope": index["public_archive_uploaded"] = True
    elif change == "workflow": index["workflow_revision"] = "d" * 40
    elif change == "version": index["version"] = "orthohmm-study-2026.10.04-rc4"
    elif change == "mode": index["files"][0]["mode"] = 0o755
    elif change == "duplicate_target": additions.append(("data.json", "new-guide.md", None))
    elif change == "escape": additions.append(("../unsafe", "new-guide.md", None))
    elif change == "revision": revision = "main"
    elif change == "missing": additions.append(("other", "absent", None))
    elif change == "symlink":
        link = tmp_path / "linked-guide.md"
        link.symlink_to(tmp_path / "new-guide.md")
        additions.append(("other", link.name, None))
    elif change == "component":
        additions.append(("other", "new-guide.md", dict(bytes=1, sha256="d" * 64)))
    with pytest.raises(ValueError):
        current.compose(tmp_path, baseline, index, parent, revision, additions)


@pytest.mark.parametrize("change", [None, "mode", "duplicate", "index", "scope", "escape"])
def test_direct_inventory_is_preserved_not_augmented(tmp_path, change):
    directory = tmp_path / "component"
    directory.mkdir()
    index = dict(schema="publication_direct_review_v3", publication_ready=False,
        redistribution_clearance=False, transitive_evidence_included=False,
        files=[dict(path="reader.py", mode=0o644, bytes=2, sha256="a" * 64)])
    if change == "mode": index["files"][0]["mode"] = 0o755
    elif change == "duplicate": index["files"].append(index["files"][0])
    elif change == "index": index["files"][0]["path"] = "REVIEW_INDEX.json"
    elif change == "scope": index["transitive_evidence_included"] = True
    elif change == "escape": index["files"][0]["path"] = "../outside"
    if change is not None:
        with pytest.raises(ValueError): current.direct_items(tmp_path, directory, index)
    else:
        items = current.direct_items(tmp_path, directory, index)
        assert items[0][:2] == ("terminal-direct-review/reader.py", "component/reader.py")
        assert items[0][2] == index["files"][0]
        assert items[1][:2] == ("terminal-direct-review/REVIEW_INDEX.json", "component/REVIEW_INDEX.json")


def test_composed_fixture_builds_archives_and_restores_with_unchanged_kernel(tmp_path):
    baseline, index, parent, additions = fixture(tmp_path)
    value = current.compose(tmp_path, baseline, index, parent, "c" * 40, additions)
    selection = tmp_path / "selection.json"
    selection.write_text(json.dumps(value))
    package = tmp_path / "package"
    built = build(tmp_path, selection, identity(selection)["sha256"], package)
    tar = tmp_path / "package.tar.gz"
    archived = archive(package, built["manifest"]["sha256"], tar)
    recovered = tmp_path / "restored"
    restored = restore(tar, archived["archive"]["sha256"], built["manifest"]["sha256"], recovered)
    assert restored == built
    assert (recovered / "history/rc5/README.md").read_bytes() == (parent / "README.md").read_bytes()
    assert restored["native_inference_repeated"] is restored["nested_component_execution_repeated"] is False


def test_actual_indices_and_reader_are_bound():
    for path, sha in ((ROOT / current.PARENT / "PACKAGE_INDEX.json", current.PARENT_INDEX_SHA),
                      (ROOT / current.PARENT / "PACKAGE_SELECTION.json", current.PARENT_SELECTION_SHA),
                      (ROOT / current.DIRECT / "REVIEW_INDEX.json", current.DIRECT_INDEX_SHA),
                      (ROOT / "benchmark_tools/bundle_publication_package.py", current.READER_SHA)):
        assert identity(path)["sha256"] == sha


def test_actual_selection_preserves_parent_and_separate_evidence(tmp_path):
    output = tmp_path / "selection.json"
    result = current.run(ROOT, output, "c" * 40)
    value = json.loads(output.read_text())
    rows = {row["target"]: row for row in value["files"]}
    inherited = [row for row in value["files"] if row["source"].startswith(current.PARENT + "/")]
    assert len(inherited) == 330
    assert value["preserved_parent_selected_payloads"] == 327
    assert value["preserved_parent_indexed_payloads"] == 329
    assert len([name for name in rows if name.startswith("terminal-direct-review/")]) == 147
    assert value["scientific_revision"] == "7f3a9e40dd7e79f842cc2c11fb8b548f9a802806"
    assert result["files"] == len(rows) == 505
    assert "evidence/rc6/results/native_qfo_terminal_failures_20261009_v1_content_readback.json" in rows
    assert "evidence/rc6/results/native_qfo_terminal_failures_20261009_v1_citations.json" in rows
    assert "native-direct-review/REVIEW_INDEX.json" in rows
    with pytest.raises(FileExistsError): current.run(ROOT, output, "c" * 40)
