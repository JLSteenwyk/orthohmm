"""Check the actual committed native-evidence review export and copied execution."""

import hashlib
import json
import os
from pathlib import Path
import subprocess
import tarfile

import pytest


ROOT = Path(__file__).resolve().parents[2]
BASE = ROOT / "benchmark_tools/results"
INDEX = "native_main_review_component_index_20261006.json"
EXECUTION = "native_main_review_component_execution_20261006.json"
INDEX_SHA = "6dffd552fa40a7a3f0f02fae8f291482264d6789773f9c4df096498b2808db14"
ARCHIVE_SHA = "95fed303f76a40f281338c519be7a58304a752b1316c9c0714336d8a9dbe6df3"


def read(name):
    return json.loads((BASE / name).read_text())


def identity(data):
    return len(data), hashlib.sha256(data).hexdigest()


def test_exact_public_index_and_actual_copied_execution_are_bound():
    index = read(INDEX)
    execution = read(EXECUTION)
    assert identity((BASE / INDEX).read_bytes()) == (41307, INDEX_SHA)
    assert execution["published_index"]["sha256"] == execution["index"]["sha256"] == INDEX_SHA
    assert execution["archive"]["sha256"] == execution["copied_archive"]["sha256"] == ARCHIVE_SHA
    assert execution["archive"]["bytes"] == execution["copied_archive"]["bytes"] == 4620877
    assert index["schema"] == "publication_direct_review_v3"
    assert index["main_text"] == "benchmark_tools/results/PUBLICATION_MAIN_TEXT_20261006.md"
    assert execution["source_revision"] == "c2a0ea998c04122f70213116d47d396ae247a048"
    assert execution["ledger_revision"] == "47b7f14ba0c2e8c5b7f070f5b11620dd356a1eb9"
    observed = execution["verification"]
    assert observed["returncode"] == 0 and observed["stderr"] == ""
    assert json.loads(observed["stdout"]) == observed["result"] == execution["build"]["result"]
    assert observed["result"]["files"] == len(index["files"]) == 113
    assert observed["result"]["payload_bytes"] == sum(r["bytes"] for r in index["files"]) == 8143531
    assert observed["result"]["page_count"] == 18
    assert observed["result"]["direct_targets"] == 85
    assert observed["result"]["local_html_occurrences"] == 88
    assert execution["archive_member_count"] == 114
    assert execution["restored_payloads_independently_checked"] == 113


def test_every_selected_payload_is_the_actual_committed_git_blob():
    index = read(INDEX)
    files = {row["path"]: row for row in index["files"]}
    assert len(files) == len(index["files"])
    for name, ref in files.items():
        revision = ref["git_revision"]
        data = subprocess.check_output(["git", "show", revision + ":" + name], cwd=ROOT)
        assert identity(data) == (ref["bytes"], ref["sha256"])
        blob = subprocess.check_output(["git", "rev-parse", revision + ":" + name],
                                       cwd=ROOT, text=True).strip()
        assert blob == ref["git_blob"]
        assert ref["mode"] in (0o644, 0o755)
    ledger = files["benchmark_tools/results/PUBLICATION_PROGRESS.md"]
    assert ledger["git_revision"] == "47b7f14ba0c2e8c5b7f070f5b11620dd356a1eb9"


def test_native_evidence_and_relative_entrypoints_are_not_historical_defaults():
    index = read(INDEX)
    files = {row["path"]: row for row in index["files"]}
    for path in index["entrypoints"].values():
        assert path in files and path.startswith("benchmark_tools/results/")
    for name in (
        "native_factorial_progress_20261005_v6/report.json",
        "native_factorial_uncertainty_binding_20261005_v4.json",
        "native_factorial_uncertainty_readback_20261005_v2.json",
        "native_qfo_scientific_scores_20261006_v1/report.json",
        "native_qfo_swiss_uncertainty_binding_22449_20261006.json",
        "recovered_native_qfo_swiss_readback_22449_20261006.json",
        "native_qfo_p0c0_figure_20261006_v1/native_qfo_p0c0.pdf",
        "native_qfo_functional_pair_sql_readback_20261006.json",
        "NATIVE_QFO_FUNCTIONAL_PAIR_RESULT_20261006.md",
    ):
        assert "benchmark_tools/results/" + name in files
    for role in ("render", "print", "review"):
        assert "20261006_v1" in index["stages"][role]
        assert index["stages"][role] in files
    # Manual review and other transitive documents are not silently included.
    assert "benchmark_tools/results/publication_main_native_visual_review_20261006.json" not in files
    assert index["transitive_evidence_included"] is False


def test_original_pre_extraction_refusal_and_complete_archive_are_preserved():
    refusal = read("native_main_component_initial_restore_refusal_20261006.json")
    execution = read(EXECUTION)
    assert refusal["phase"] == "archive_member_pre_extraction_guard"
    assert refusal["exception"] == "AssertionError"
    assert refusal["archive_extracted"] is refusal["copied_verifier_executed"] is False
    assert refusal["archive"] == execution["archive"]
    assert refusal["observed_member_outside_assumed_modes"] == dict(
        name="REVIEW_INDEX.json", is_regular_file=True, mode_octal="0664")
    assert execution["archive_build_repeated"] is False
    assert execution["initial_refusal_preserved"] is True
    assert execution["archive_construction"]["restore_filter"] == "data"
    assert execution["archive_construction"]["index_mode"] == "0664"


def test_retained_file_trace_checks_actual_copied_payloads_without_checkout():
    execution = read(EXECUTION)
    trace = (BASE / "native_main_review_component_verify_20261006.strace").read_bytes()
    ref = execution["published_trace"]
    assert identity(trace) == (ref["bytes"], ref["sha256"])
    assert ref["sha256"] == execution["file_trace"]["sha256"]
    trace = trace.decode()
    for prefix in execution["forbidden_trace_prefixes"]:
        assert prefix not in trace
    assert execution["original_paths_observed_in_trace"] == []
    observed = execution["verification"]
    assert observed["environment"]["PATH"] == "/no-git"
    command = observed["command"]
    assert command[:4] == ["/usr/bin/strace", "-f", "-e", "trace=%file"]
    assert command[6:10] == ["/usr/bin/python3", "-I", "-S", "-B"]
    restored = Path(execution["restored_directory"])
    assert not restored.is_relative_to(ROOT.parent)
    assert str(restored / "REVIEW_INDEX.json") in trace
    for ref in read(INDEX)["files"]:
        assert str(restored / ref["path"]) in trace


def test_actual_local_archive_matches_all_selected_bytes_when_explicitly_supplied():
    name = os.environ.get("ORTHOHMM_NATIVE_REVIEW_ARCHIVE")
    if name is None:
        pytest.skip("Supply retained local archive explicitly; it is not distributed with raw artifacts")
    archive = Path(name)
    assert identity(archive.read_bytes()) == (4620877, ARCHIVE_SHA)
    refs = {r["path"]: r for r in read(INDEX)["files"]}
    with tarfile.open(archive, "r:gz") as stream:
        members = stream.getmembers()
        assert len(members) == 114 and {m.name for m in members} == set(refs) | {"REVIEW_INDEX.json"}
        for member in members:
            assert member.isfile()
            data = stream.extractfile(member).read()
            if member.name == "REVIEW_INDEX.json":
                assert data == (BASE / INDEX).read_bytes() and member.mode == 0o664
            else:
                ref = refs[member.name]
                assert identity(data) == (ref["bytes"], ref["sha256"])
                assert member.mode == ref["mode"]


def test_archive_scope_does_not_claim_full_study_reproduction_or_readiness():
    index = read(INDEX)
    execution = read(EXECUTION)
    for flag in ("publication_ready", "redistribution_clearance", "transitive_evidence_included"):
        assert index[flag] is execution[flag] is False
    assert execution["native_or_scoring_repeated"] is False
    assert execution["public_archive_uploaded"] is False
    assert execution["new_bootstrap_draws"] == 0
    assert execution["verification"]["result"]["inference_reproduced"] is False
    assert any("not OS containment" in s for s in execution["limitations"])
    assert any("not a direct render target" in s for s in execution["limitations"])
