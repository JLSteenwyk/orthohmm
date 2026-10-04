"""The new review candidate extends, rather than silently replacing, rc3."""

import json
from pathlib import Path

import pytest

from benchmark_tools import prepare_publication_package_rc4 as prepare
from benchmark_tools.bundle_publication_package import identity


ROOT = Path(__file__).resolve().parents[2]


def fixture(tmp_path):
    files = []
    for target, source in (("README.md", "old-guide.md"), ("data/retained.json", "counts.json")):
        path = tmp_path / source
        path.write_text(source)
        files.append({"target": target, "source": source, **identity(path)})
    (tmp_path / "new-guide.md").write_text("new guide")
    baseline = {"schema": "publication_package_selection_v1", "version": "orthohmm-study-2026.10.04-rc3",
                "scientific_revision": "a" * 40, "workflow_revision": "b" * 40,
                "files": files, "limitations": ["retained limit"], "publication_ready": False,
                "public_archive_uploaded": False}
    return baseline, [("README.md", "new-guide.md")]


def test_parent_scientific_scope_and_payload_bytes_are_preserved(tmp_path):
    baseline, additions = fixture(tmp_path)
    output = prepare.compose(tmp_path, baseline, "c" * 40, additions)
    by_target = {r["target"]: r for r in output["files"]}
    assert by_target["data/retained.json"] == baseline["files"][1]
    assert by_target["evidence/PUBLICATION_PACKAGE_RC3_20261004.md"]["sha256"] == baseline["files"][0]["sha256"]
    assert output["scientific_revision"] == baseline["scientific_revision"]
    assert output["publication_ready"] is False and output["public_archive_uploaded"] is False
    assert output["preserved_parent_payloads"] == output["preserved_parent_selected_payloads"] == 2
    assert "retained limit" in output["limitations"]


@pytest.mark.parametrize("change", ["bytes", "hash", "duplicate", "missing", "revision", "escape", "unsafe_target"])
def test_extension_refuses_changed_inheritance_or_unsafe_payload(tmp_path, change):
    baseline, additions = fixture(tmp_path)
    revision = "c" * 40
    if change == "bytes": baseline["files"][1]["bytes"] += 1
    elif change == "hash": baseline["files"][1]["sha256"] = "d" * 64
    elif change == "duplicate": additions.append(("data/retained.json", "new-guide.md"))
    elif change == "missing": additions.append(("other", "absent"))
    elif change == "revision": revision = "main"
    elif change == "escape":
        outside = tmp_path.parent / "outside-guide"
        outside.write_text("outside")
        additions.append(("other", "../outside-guide"))
    else: additions.append(("../other", "new-guide.md"))
    with pytest.raises(ValueError): prepare.compose(tmp_path, baseline, revision, additions)


def test_actual_rc3_sources_and_reader_still_match_the_prior_payloads():
    path = ROOT / "benchmark_tools/results/publication_package_rc3_selection_20261004.json"
    assert identity(path)["sha256"] == prepare.BASELINE_SHA
    parent = json.loads(path.read_bytes())
    assert len(parent["files"]) == 137
    assert all(identity(ROOT / r["source"]) == {k: r[k] for k in ("bytes", "sha256")} for r in parent["files"])
    index = json.loads((ROOT / "benchmark_tools/results/publication_package_rc3_index_20261004.json").read_bytes())
    reader = next(r for r in index["files"] if r["path"] == "bundle_publication_package.py")
    assert identity(ROOT / "benchmark_tools/bundle_publication_package.py") == {k: reader[k] for k in ("bytes", "sha256")}
    assert reader["sha256"] == prepare.READER_SHA


def test_new_review_and_metadata_are_selected_separately_from_old_artifacts():
    assert "PUBLICATION_MAIN_TEXT_20261004_v3.md" in prepare.EVIDENCE
    assert {p[0] for p in prepare.REPORTS} == {"inventory.json", "diagnostic.json", "linkage.json", "qfo_register.json",
                                            "factorial_resources.json", "all_tool_register.json", "family_evidence.tsv"}
    source = Path(prepare.__file__).read_text()
    assert '"current-review/" + target' in source
    assert '"history/rc3/PACKAGE_SELECTION.json"' in source
    assert '"history/rc3/bundle_publication_package.py"' in source
