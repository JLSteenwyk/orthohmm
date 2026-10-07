"""Check the actual new delivery, immutable inheritance and copied file-access observations."""

import hashlib
import json
import os
from pathlib import Path
import tarfile

import pytest


ROOT = Path(__file__).resolve().parents[2]
BASE = ROOT / "benchmark_tools/results"
INDEX = "publication_package_rc5_index_20261006.json"
SELECTION = "publication_package_rc5_selection_20261006.json"
EXECUTION = "publication_package_rc5_execution_20261006.json"
INDEX_SHA = "66bc02d7821d1097808a870c38198bded9c88c379c24052b1828867af8692922"
ARCHIVE_SHA = "24f27934b6a7f8cb75ace61f1a0448a215f0a6d8f7eb7f6d891ae7f4d49adeb7"


def read(name):
    return json.loads((BASE / name).read_text())


def digest(data):
    return len(data), hashlib.sha256(data).hexdigest()


def test_selection_actual_package_and_execution_match_external_anchors():
    selection, index, execution = read(SELECTION), read(INDEX), read(EXECUTION)
    assert digest((BASE / INDEX).read_bytes()) == (74998, INDEX_SHA)
    assert digest((BASE / SELECTION).read_bytes()) == (110409, execution["selection"]["sha256"])
    assert execution["selection"]["sha256"] == "3825124e5a4397d713911d236766cb9074f29562ad986eb12c59473f6f5cc724"
    assert index["version"] == selection["version"] == execution["version"] == "orthohmm-study-2026.10.06-rc5"
    assert len(selection["files"]) == execution["selected_files"] == 327
    assert len(index["files"]) == execution["indexed_payloads"] == 329
    assert sum(r["bytes"] for r in index["files"]) == execution["payload_bytes"] == 167214767
    assert execution["archive"]["sha256"] == execution["copied_archive"]["sha256"] == ARCHIVE_SHA
    assert execution["archive"]["bytes"] == 138737221
    for stage in execution["executions"].values():
        assert stage["returncode"] == 0 and json.loads(stage["stdout"]) == stage["result"]
    assert execution["executions"]["restore"]["result"] == execution["build"]
    assert execution["executions"]["outer_verify"]["result"] == execution["build"]


def test_every_parent_selected_and_indexed_payload_is_preserved_without_mutable_sources():
    parent = read("publication_package_rc4_index_20261004.json")
    selection = read(SELECTION)
    current = {r["path"]: r for r in read(INDEX)["files"]}
    selected = {r["target"]: r for r in selection["files"]}
    assert len(parent["files"]) == selection["preserved_parent_indexed_payloads"] == 183
    assert selection["preserved_parent_selected_payloads"] == 181
    for row in parent["files"]:
        name = row["path"]
        target = "history/rc4/" + name if name in (
            "README.md", "PACKAGE_SELECTION.json", "bundle_publication_package.py") else name
        assert (current[target]["bytes"], current[target]["sha256"]) == (row["bytes"], row["sha256"])
        assert selected[target]["source"] == "benchmarks/work/orthohmm-study-2026.10.04-rc4/" + name
    assert current["history/rc4/PACKAGE_INDEX.json"]["sha256"] == "550d1b06875a26d5719376d5f8562055f68b87ba5d05619b1273146c4c51cf71"


def test_copied_direct_child_is_verified_after_outer_mode_normalization():
    index, execution = read(INDEX), read(EXECUTION)
    files = {r["path"]: r for r in index["files"]}
    child = read("native_main_review_component_index_20261006_v2.json")
    for row in child["files"]:
        kept = files["native-direct-review/" + row["path"]]
        assert (kept["bytes"], kept["sha256"], kept["mode"]) == (row["bytes"], row["sha256"], row["mode"])
    result = execution["executions"]["child_verify"]["result"]
    assert (result["files"], result["direct_targets"], result["local_html_occurrences"], result["page_count"]) == (121, 99, 102, 19)
    assert files["native-review/document.pdf"]["sha256"] == "e574d508993eb3c8a0a2a02b0a5fd533601b3fe9a3e31e4c88f8109afc7f3901"
    assert execution["presentation"]["assembled_pages"] == 41


def test_trace_absence_covers_every_outer_and_direct_payload_at_new_locations():
    execution, index = read(EXECUTION), read(INDEX)
    directory = Path(execution["restored_directory"])
    assert not directory.is_relative_to(ROOT.parent)
    for ref, scope in zip(execution["file_access_traces"], ("outer", "child")):
        data = (BASE / Path(ref["path"]).name).read_bytes()
        assert digest(data) == (ref["bytes"], ref["sha256"])
        text = data.decode()
        for prefix in execution["forbidden_trace_prefixes"]: assert prefix not in text
        if scope == "outer":
            for row in index["files"]: assert str(directory / row["path"]) in text
        else:
            for row in read("native_main_review_component_index_20261006_v2.json")["files"]:
                assert str(directory / "native-direct-review" / row["path"]) in text
    assert execution["environment"]["PATH"] == "/no-git"
    assert execution["original_paths_observed_in_trace"] == []


def test_full_goal_and_unmet_scientific_requirements_not_redefined_by_delivery():
    execution = read(EXECUTION)
    for flag in ("publication_ready", "public_archive_uploaded", "redistribution_clearance",
                 "hermetic_study_reproduction", "native_inference_repeated", "scoring_or_bootstrap_repeated",
                 "scientific_defaults_changed", "unrelated_workloads_modified"):
        assert execution[flag] is False
    guide = (ROOT / "benchmark_tools/PUBLICATION_PACKAGE_RC5_20261006.md").read_text()
    for phrase in ("Four fresh native cells and twelve contrasts", "independent-family confirmation",
                   "not a hermetic full-study archive", "not proof of biological", "unknown, tool-dependent"):
        assert phrase in " ".join(guide.split())


def test_local_archive_full_inventory_when_explicitly_supplied():
    name = os.environ.get("ORTHOHMM_PUBLICATION_RC5_ARCHIVE")
    if name is None: pytest.skip("Supply the retained local rc5 archive explicitly")
    archive = Path(name)
    assert digest(archive.read_bytes()) == (138737221, ARCHIVE_SHA)
    rows = {r["path"]: r for r in read(INDEX)["files"]}
    with tarfile.open(archive, "r:gz") as stream:
        members = stream.getmembers()
        assert len(members) == 330 and {m.name for m in members} == set(rows) | {"PACKAGE_INDEX.json"}
        for member in members:
            assert member.isfile() and member.mode == 0o644
            data = stream.extractfile(member).read()
            if member.name == "PACKAGE_INDEX.json": assert data == (BASE / INDEX).read_bytes()
            else: assert digest(data) == (rows[member.name]["bytes"], rows[member.name]["sha256"])
