from copy import deepcopy
import hashlib

import numpy as np
import pytest

from benchmark_tools.complete_fas_population import (
    ORIGINAL_COMMIT, check_database_state, check_sample, unique_records, validate_prior, verify_git_pin,
)
from benchmark_tools.audit_fas_population import encode
from benchmark_tools.prepare_ob_candidate_neighborhood import record
from tests.unit.test_render_fas_population import fixture


def prior_fixture():
    report, manifest = fixture()
    rows = report["methods"][:6]
    prior = dict(schema="qfo_fas_population_terminal_partial_v1",
        status="timed_out_incomplete_population_recount", job_id=22382,
        source_commit=ORIGINAL_COMMIT, complete_report_produced=False,
        uncertainty_admitted=False, benchmark_scores_changed=False, publication_ready=False,
        completed_partial_rows=[dict(method=row["method"]) for row in rows],
        missing_method_rows=[m["key"] for m in manifest["methods"][6:]])
    return prior, manifest, deepcopy(rows)


def test_reuse_requires_exact_complete_prefix_and_two_missing():
    prior, manifest, rows = prior_fixture()
    validate_prior(prior, manifest, rows)


@pytest.mark.parametrize("change", ["success", "complete", "short", "reorder", "select",
                                  "missing", "score", "sum", "historical"])
def test_rejects_relabeling_or_changed_reused_data(change):
    prior, manifest, rows = prior_fixture()
    if change == "success":
        prior["status"] = "completed"
    elif change == "complete":
        prior["complete_report_produced"] = True
    elif change == "short":
        rows.pop()
    elif change == "reorder":
        rows.reverse()
    elif change == "select":
        rows[5]["method"] = "7"
    elif change == "missing":
        prior["missing_method_rows"].reverse()
    elif change == "score":
        prior["benchmark_scores_changed"] = True
    elif change == "sum":
        rows[0]["precomputed_score_sum"] = True
    else:
        rows[0]["database_historically_hash_bound"] = True
    with pytest.raises(ValueError):
        validate_prior(prior, manifest, rows)


def test_sample_rechecks_unknown_identifiers_without_population_rescan():
    ids = {"A": 0, "B": 1}
    codes, values = np.array([encode(0, 1)], dtype=np.uint64), np.array([.8])
    saved = {("A", "B"): .8, ("C", "D"): .1}
    assert check_sample(saved, 1, ids, codes, values)
    assert ids == {"A": 0, "B": 1, "C": 2, "D": 3}
    assert not check_sample(saved, 0, ids, codes, values)
    assert not check_sample({("A", "B"): .7}, 1, ids, codes, values)


@pytest.mark.parametrize("count", [-1, 2, True, 1.0])
def test_bad_sample_stratum_rejected(count):
    with pytest.raises(ValueError):
        check_sample({("A", "B"): .8}, count, {"A": 0, "B": 1},
                     np.array([encode(0, 1)], dtype=np.uint64), np.array([.8]))


def test_record_deduplication_rejects_identity_conflicts():
    pin = dict(path="/x", bytes=1, sha256="a")
    assert unique_records([pin, dict(pin)]) == [pin]
    with pytest.raises(ValueError):
        unique_records([pin, dict(pin, sha256="b")])


def test_git_proof_metadata_does_not_change_file_pin(tmp_path, monkeypatch):
    path = tmp_path / "source.py"
    path.write_bytes(b"trusted\n")
    pin = dict(path=str(path), bytes=8, sha256=hashlib.sha256(b"trusted\n").hexdigest(),
               git_revision=ORIGINAL_COMMIT, git_blob_matches=True)
    monkeypatch.setattr("benchmark_tools.complete_fas_population.subprocess.check_output",
                        lambda *args, **kwargs: b"trusted\n")
    verify_git_pin(pin, tmp_path)
    monkeypatch.setattr("benchmark_tools.complete_fas_population.subprocess.check_output",
                        lambda *args, **kwargs: b"different\n")
    with pytest.raises(ValueError):
        verify_git_pin(pin, tmp_path)


@pytest.mark.parametrize("suffix", ["-wal", "-shm", "-journal"])
def test_reuse_rejects_sidecars_without_querying_database(tmp_path, suffix):
    path = tmp_path / "database.db"
    path.write_bytes(b"unchanged")
    pin = record(path)
    check_database_state(pin)
    path.with_name(path.name + suffix).write_bytes(b"new state")
    with pytest.raises(ValueError, match="sidecar"):
        check_database_state(pin)
