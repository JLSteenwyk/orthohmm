"""Actual new archive coverage, with no inference or historical-package replay."""

import hashlib
import json
from pathlib import Path, PurePosixPath
import tarfile

ROOT = Path(__file__).resolve().parents[2]
BASE = ROOT / "benchmark_tools/results"
INDEX_SHA = "41d31073e9c563ef0c530187d666a0b7e83d63d894038f19d54f2614877b8d53"
ARCHIVE_SHA = "42c8d5daced8cef058e1f31cd23f80476522be75184b3ef23ea17951d81b3a48"


def test_new_archive_has_exact_anchored_regular_payload_inventory():
    path = BASE / "current_error_strata_direct_review_20261007_v1.tar.gz"
    assert path.stat().st_size == 6354070
    assert hashlib.sha256(path.read_bytes()).hexdigest() == ARCHIVE_SHA
    with tarfile.open(path, "r:gz") as archive:
        members = archive.getmembers()
        names = [m.name for m in members]
        assert len(names) == len(set(names)) == 150
        assert all(m.isfile() and not PurePosixPath(m.name).is_absolute()
                   and ".." not in PurePosixPath(m.name).parts for m in members)
        raw_index = archive.extractfile("REVIEW_INDEX.json").read()
        assert hashlib.sha256(raw_index).hexdigest() == INDEX_SHA
        index = json.loads(raw_index)
        rows = {r["path"]: r for r in index["files"]}
        assert set(names) == set(rows) | {"REVIEW_INDEX.json"}
        assert len(rows) == 149 and sum(r["bytes"] for r in rows.values()) == 9272442
        for member in members:
            if member.name == "REVIEW_INDEX.json":
                continue
            row = rows[member.name]
            assert member.size == row["bytes"] and member.mode == row["mode"]
            assert hashlib.sha256(archive.extractfile(member).read()).hexdigest() == row["sha256"]
        assert index["schema"] == "publication_direct_review_v3"
        assert index["main_text"].endswith("PUBLICATION_MAIN_TEXT_20261007_v3.md")
        for name in (
            "PUBLICATION_MAIN_TEXT_20261007_v3.md",
            "swiss_model_divergence_strata_20261007_v1/TABLE.md",
            "swiss_model_divergence_evidence_23932_v1.tar.gz",
            "native_qfo_three_cell_strata_figure_20261007_v1/fixed_stratum_contrasts.pdf",
            "publication_main_print_20261007_v3/document.pdf",
        ):
            assert "benchmark_tools/results/" + name in rows
        assert "benchmark_tools/results/PUBLICATION_PROGRESS.md" not in rows
        assert "benchmark_tools/results/publication_main_pdf_visual_review_20261007_v3.json" not in rows
        for key in ("publication_ready", "redistribution_clearance", "transitive_evidence_included"):
            assert index[key] is False


def test_new_build_restore_and_actual_copied_verification_scopes_match():
    build, restored, verified = [json.loads((BASE / name).read_text()) for name in (
        "current_error_strata_review_build_20261007_v1.json",
        "current_error_strata_review_restore_20261007_v1.json",
        "current_error_strata_review_verify_20261007_v1.json")]
    assert build["returncode"] == verified["returncode"] == 0
    assert build["result"] == verified["result"]
    assert verified["destination"].startswith("/tmp/")
    assert verified["python_flags"] == ["-I", "-S", "-B"]
    assert restored["status"] == "anchored_payloads_restored"
    assert restored["payloads"] == 149 and restored["members"] == 150
    assert restored["archive"]["sha256"] == ARCHIVE_SHA
    assert restored["index_sha256"] == verified["result"]["manifest"]["sha256"] == INDEX_SHA
    assert restored["copied_code_executed"] is False
    assert verified["status"] == "new_direct_review_copy_verified_outside_checkout"
    assert verified["result"]["inference_reproduced"] is False
    assert build["visual_receipt_in_component"] is False
    assert verified["old_main_or_rc5_rebuilt"] is False
    assert verified["publication_ready"] is False
