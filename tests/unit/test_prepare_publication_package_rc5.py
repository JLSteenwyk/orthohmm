"""The new package inherits immutable payloads, never drifted original source paths."""

import copy
import json
from pathlib import Path

import pytest

from benchmark_tools import prepare_publication_package_rc5 as current
from benchmark_tools.bundle_publication_package import identity


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
            original.write_text("changed original worktree bytes")
            files.append(dict(target=name, source=original.name, **identity(path)))
    baseline = dict(schema="publication_package_selection_v1", version="orthohmm-study-2026.10.04-rc4",
        scientific_revision="a" * 40, workflow_revision="b" * 40, files=files,
        publication_ready=False, public_archive_uploaded=False, limitations=["retained limit"])
    index = dict(version=baseline["version"], scientific_revision=baseline["scientific_revision"],
                 files=rows, publication_ready=False, native_inference_repeated=False)
    (parent / "PACKAGE_INDEX.json").write_text(json.dumps(index))
    (tmp_path / "new-guide.md").write_text("new guide")
    return baseline, index, parent, [("README.md", "new-guide.md", None)]


def test_immutable_parent_selected_and_support_bytes_preserved_despite_original_drift(tmp_path):
    baseline, index, parent, additions = fixture(tmp_path)
    before = copy.deepcopy(baseline)
    result = current.compose(tmp_path, baseline, index, parent, "c" * 40, additions)
    rows = {r["target"]: r for r in result["files"]}
    assert baseline == before
    assert rows["data.json"]["source"] == "parent/data.json"
    assert rows["history/rc4/README.md"]["sha256"] == baseline["files"][0]["sha256"]
    assert rows["history/rc4/PACKAGE_SELECTION.json"]["source"] == "parent/PACKAGE_SELECTION.json"
    assert rows["history/rc4/bundle_publication_package.py"]["source"] == "parent/bundle_publication_package.py"
    assert "history/rc4/PACKAGE_INDEX.json" in rows
    assert result["preserved_parent_selected_payloads"] == 2
    assert result["preserved_parent_indexed_payloads"] == result["preserved_parent_payloads"] == 4
    assert result["scientific_revision"] == baseline["scientific_revision"]
    assert result["publication_ready"] is result["public_archive_uploaded"] is False


@pytest.mark.parametrize("change", ["payload", "index", "duplicate", "escape", "revision", "missing", "symlink"])
def test_changed_or_unsafe_inheritance_refused(tmp_path, change):
    baseline, index, parent, additions = fixture(tmp_path)
    revision = "c" * 40
    if change == "payload": (parent / "data.json").write_text("tampered")
    elif change == "index": index["files"][0]["bytes"] += 1
    elif change == "duplicate": additions.append(("data.json", "new-guide.md", None))
    elif change == "escape": additions.append(("../unsafe", "new-guide.md", None))
    elif change == "revision": revision = "main"
    elif change == "missing": additions.append(("other", "absent", None))
    elif change == "symlink":
        link = tmp_path / "linked-guide.md"
        link.symlink_to(tmp_path / "new-guide.md")
        additions.append(("other", link.name, None))
    with pytest.raises(ValueError): current.compose(tmp_path, baseline, index, parent, revision, additions)


def test_direct_component_modes_cannot_be_silently_normalized(tmp_path):
    component = tmp_path / "component"
    component.mkdir()
    index = dict(schema="publication_direct_review_v3", publication_ready=False,
                 redistribution_clearance=False, transitive_evidence_included=False,
                 files=[dict(path="reader.py", mode=0o644, bytes=2, sha256="a" * 64)])
    rows = current.direct_items(tmp_path, component, index)
    assert rows[0][:2] == ("native-direct-review/reader.py", "component/reader.py")
    index["files"][0]["mode"] = 0o755
    with pytest.raises(ValueError): current.direct_items(tmp_path, component, index)


def test_actual_immutable_parent_and_direct_indices_are_bound():
    for path, sha in ((ROOT / current.PARENT / "PACKAGE_INDEX.json", current.PARENT_INDEX_SHA),
                      (ROOT / current.PARENT / "PACKAGE_SELECTION.json", current.PARENT_SELECTION_SHA),
                      (ROOT / current.DIRECT / "REVIEW_INDEX.json", current.DIRECT_INDEX_SHA)):
        assert identity(path)["sha256"] == sha


def test_new_selected_inventory_uses_parent_package_and_preserves_all_counts(tmp_path):
    output = tmp_path / "selection.json"
    result = current.run(ROOT, output, "c" * 40)
    selection = json.loads(output.read_text())
    assert selection["preserved_parent_selected_payloads"] == 181
    assert selection["preserved_parent_indexed_payloads"] == 183
    inherited = [r for r in selection["files"] if r["source"].startswith(current.PARENT + "/")]
    assert len(inherited) == 184
    direct = [r for r in selection["files"] if r["target"].startswith("native-direct-review/")]
    assert len(direct) == 122
    assert result["files"] == len(selection["files"])
    with pytest.raises(FileExistsError): current.run(ROOT, output, "c" * 40)
