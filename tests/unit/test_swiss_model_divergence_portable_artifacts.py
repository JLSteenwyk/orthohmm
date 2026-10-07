"""Actual relocated results and serialized component, not another selected replay."""

import hashlib
import json
import math
from pathlib import Path
import tarfile

ROOT = Path(__file__).resolve().parents[2]
BASE = ROOT / "benchmark_tools/results"
INDEX_SHA = "83e3b34b7b99de6f23b6f2c4a4c2cbf28d001c0ac8708b06f427a909346bf3e2"


def load(name):
    return json.loads((BASE / name).read_text())


def test_serialized_new_component_matches_actual_executed_sources_and_inputs():
    path = BASE / "swiss_model_divergence_portable_20261007_v2.tar.gz"
    assert path.stat().st_size == 693861
    assert hashlib.sha256(path.read_bytes()).hexdigest() == "6d995e04580942920301abcbcf3837188c8a8e8480dcd2eb9c13cf30c440fd0c"
    with tarfile.open(path, "r:gz") as archive:
        members = archive.getmembers()
        assert len(members) == len({m.name for m in members}) == 16
        assert all(m.isfile() and m.mode == 0o644 and "/" not in m.name for m in members)
        raw_index = archive.extractfile("REPLAY_INDEX.json").read()
        assert len(raw_index) == 4815 and hashlib.sha256(raw_index).hexdigest() == INDEX_SHA
        index = json.loads(raw_index)
        assert index["schema"] == "swiss_model_divergence_portable_v2"
        assert len(index["files"]) == 15
        assert sum(r["bytes"] for r in index["files"]) == 925023
        assert {m.name for m in members} == {r["path"] for r in index["files"]} | {"REPLAY_INDEX.json"}
        for row in index["files"]:
            data = archive.extractfile(row["path"]).read()
            assert len(data) == row["bytes"]
            assert hashlib.sha256(data).hexdigest() == row["sha256"]
            assert (ROOT / row["git_path"]).read_bytes() == data
            assert row["git_revision"] == "2ea4be041226bdfe0f06bce5a093b6201f18c16b"
        assert hashlib.sha256(archive.extractfile("replay.py").read()).hexdigest() == "a87e81fe6015e3c12ac29ba2572f1494d3fb2e17e90215d170847b51abc9dbc2"


def test_actual_replay_matches_complete_previously_verified_descriptors_and_bins():
    actual = load("swiss_model_divergence_portable_replayed_20261007_v2.json")
    prior = load("swiss_model_divergence_readback_23932_v1.json")
    assert tuple(actual[k] for k in ("families", "proteins", "pairs", "family_rows", "score_rows", "differences")) == (18, 563, 10765, 54, 9, 6)
    assert actual["restored_payloads"] == 183 and actual["component_files"] == 15
    assert actual["bins"] == prior["strata"]
    assert math.isclose(actual["median_family_distance"], prior["median_family_distance"], rel_tol=0, abs_tol=1e-10)
    for family, computed in actual["family_descriptors"].items():
        assert set(computed) == set(prior["features"][family])
        for key, value in computed.items():
            expected = prior["features"][family][key]
            if key in ("members", "pairs", "unit"):
                assert value == expected
            else:
                assert math.isclose(value, expected, rel_tol=0, abs_tol=1e-10)
    for flag in ("native_inference_reproduced", "raw_counts_readmitted", "scientific_timings_admitted",
                 "publication_ready", "independent_confirmation"):
        assert actual[flag] is False
    assert actual["new_bootstrap_draws"] == 0


def test_actual_guarded_execution_has_real_canaries_and_bounded_scope():
    execution = load("swiss_model_divergence_portable_execution_20261007_v2.json")
    actual = load("swiss_model_divergence_portable_replayed_20261007_v2.json")
    assert execution["returncode"] == 0 and execution["failed_v1_preserved"] is True
    assert execution["precheck"]["index_sha256"] == actual["index_sha256"] == INDEX_SHA
    assert execution["result"]["summary"] == {k: v for k, v in actual.items() if k != "family_descriptors"}
    guard = execution["result"]["guard"]
    assert len(guard["canaries"]) == 2 and all(r["event"] == "open" for r in guard["canaries"])
    assert guard["replay_forbidden_events"] == [] and guard["os_containment"] is False
    assert guard["runtime_site_packages_explicitly_allowed"] is True
    assert execution["publication_ready"] is False
    failure = load("swiss_model_divergence_portable_execution_20261007_v1.json")
    assert failure["returncode"] == 1 and failure["automatic_retry"] is False
