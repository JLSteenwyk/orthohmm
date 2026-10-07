"""Read back the selected committed snapshot and actual copied-verifier observations."""

import hashlib
import json
import os
from pathlib import Path
import subprocess
import tarfile

import pytest


ROOT = Path(__file__).resolve().parents[2]
BASE = ROOT / "benchmark_tools/results"
INDEX = "native_main_review_component_index_20261006_v2.json"
EXECUTION = "native_main_review_component_execution_20261006_v2.json"
INDEX_SHA = "c8cc4899338689b9b6cd033e6ed7017f2bd3c2e152633a55869923cf7bb5d82a"
ARCHIVE_SHA = "889be71851fcfd40dcca6ef9656af286eb5951309303327d5d7513f3917e59bd"


def read(name):
    return json.loads((BASE / name).read_text())


def identity(data):
    return len(data), hashlib.sha256(data).hexdigest()


def test_actual_snapshot_and_restored_execution_are_bound():
    index, execution = read(INDEX), read(EXECUTION)
    assert identity((BASE / INDEX).read_bytes()) == (44163, INDEX_SHA)
    assert execution["published_index"]["sha256"] == INDEX_SHA
    assert execution["archive"]["sha256"] == execution["copied_archive"]["sha256"] == ARCHIVE_SHA
    assert execution["archive"]["bytes"] == execution["copied_archive"]["bytes"] == 3544952
    observed = execution["verification"]
    assert observed["returncode"] == 0 and json.loads(observed["stdout"]) == observed["result"]
    result = observed["result"]
    assert result["files"] == len(index["files"]) == 121
    assert result["payload_bytes"] == sum(r["bytes"] for r in index["files"]) == 7697373
    assert (result["direct_targets"], result["local_html_occurrences"], result["page_count"]) == (99, 102, 19)
    assert execution["archive_member_count"] == 122
    restoration = read("native_main_review_component_restoration_20261006_v2.json")
    assert restoration["payloads"] == 121 and restoration["members"] == 122
    assert restoration["copied_code_executed"] is False


def test_all_payloads_are_the_selected_committed_bytes_and_modes():
    index = read(INDEX)
    for row in index["files"]:
        data = subprocess.check_output(["git", "show", row["git_revision"] + ":" + row["path"]], cwd=ROOT)
        assert identity(data) == (row["bytes"], row["sha256"])
        blob = subprocess.check_output(["git", "rev-parse", row["git_revision"] + ":" + row["path"]],
                                       text=True, cwd=ROOT).strip()
        assert blob == row["git_blob"]
        assert row["mode"] in (0o644, 0o755)
    ledger = next(r for r in index["files"] if r["path"].endswith("/PUBLICATION_PROGRESS.md"))
    assert ledger["git_revision"] == "f5eb601cf30c509e9fc7affe7b72f131598e7fec"


def test_original_path_trace_absence_covers_every_restored_payload():
    execution = read(EXECUTION)
    data = (BASE / "native_main_review_component_verify_20261006_v2.strace").read_bytes()
    assert identity(data) == (execution["published_trace"]["bytes"], execution["published_trace"]["sha256"])
    text = data.decode()
    assert execution["original_paths_observed_in_trace"] == []
    for prefix in execution["forbidden_trace_prefixes"]: assert prefix not in text
    restored = Path(execution["restored_directory"])
    assert not restored.is_relative_to(ROOT.parent)
    for row in read(INDEX)["files"]: assert str(restored / row["path"]) in text
    assert str(restored / "REVIEW_INDEX.json") in text
    assert execution["verification"]["environment"]["PATH"] == "/no-git"
    assert execution["verification"]["command"][6:10] == ["/usr/bin/python3", "-I", "-S", "-B"]


def test_current_native_evidence_and_archive_limitations_not_silently_redefined():
    index, execution = read(INDEX), read(EXECUTION)
    names = {r["path"] for r in index["files"]}
    for name in ("native_qfo_scientific_scores_20261006_v2/report.json",
                 "native_qfo_candidate_swiss_uncertainty_20261006_v1.json",
                 "native_qfo_candidate_swiss_readback_20261006_v1.json",
                 "native_qfo_candidate_vgnc_20261006_v1/report.json",
                 "native_qfo_candidate_vgnc_readback_20261006_v1.json",
                 "native_qfo_candidate_alias_group_20261006_v1/report.json",
                 "native_qfo_candidate_alias_group_20261006_v1/changed_pair_groups.tsv",
                 "native_qfo_candidate_alias_group_readback_20261006_v1.json",
                 "native_qfo_three_cell_figure_20261006_v1/native_qfo_three_cell.pdf"):
        assert "benchmark_tools/results/" + name in names
    assert "benchmark_tools/results/native_main_visual_review_20261006_v2.json" not in names
    assert execution["assembled_review"]["included_in_direct_component"] is False
    for flag in ("publication_ready", "redistribution_clearance", "transitive_evidence_included"):
        assert index[flag] is execution[flag] is execution["verification"]["result"][flag] is False
    assert execution["native_or_scoring_repeated"] is False
    assert execution["new_bootstrap_draws"] == 0 and execution["public_archive_uploaded"] is False


def test_assembled_review_has_seventeen_figures_and_exact_source_preservation_receipt():
    report = read("publication_main_with_figures_20261006_v2/assembly.json")
    assert (report["main_pages"], report["guide_pages"], report["figure_pages"], report["total_pages"]) == (19, 5, 17, 41)
    assert report["original_text_and_figure_pages_pixel_identical"] is True
    assert len(report["preserved_source_pages"]) == 36
    assert report["figures"][-1]["combined_page"] == 41
    assert any(r["combined_page"] == 41 for r in report["redirected_main_figure_links"])
    visual = read("native_assembled_visual_review_20261006_v2.json")
    assert visual["pdf_sha256"] == report["pdf"]["sha256"]
    assert visual["full_document_visual_certification"] is False


def test_actual_local_archive_when_explicitly_supplied():
    name = os.environ.get("ORTHOHMM_NATIVE_REVIEW_ARCHIVE_V2")
    if name is None: pytest.skip("Supply the retained local archive explicitly")
    archive = Path(name)
    assert identity(archive.read_bytes()) == (3544952, ARCHIVE_SHA)
    rows = {r["path"]: r for r in read(INDEX)["files"]}
    with tarfile.open(archive, "r:gz") as stream:
        members = stream.getmembers()
        assert len(members) == 122 and {m.name for m in members} == set(rows) | {"REVIEW_INDEX.json"}
        for member in members:
            assert member.isfile()
            data = stream.extractfile(member).read()
            if member.name == "REVIEW_INDEX.json":
                assert data == (BASE / INDEX).read_bytes() and member.mode == 0o664
            else:
                row = rows[member.name]
                assert identity(data) == (row["bytes"], row["sha256"])
                assert member.mode == row["mode"]
