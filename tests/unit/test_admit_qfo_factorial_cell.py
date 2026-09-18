from copy import deepcopy
from pathlib import Path

import pytest

from benchmark_tools.admit_qfo_factorial_cell import require_success, check_pairs
from benchmark_tools.validate_factorial_native import check_native_origin


def fixture():
    accounting = "JobID|JobIDRaw|State|ExitCode|Elapsed|NodeList|AllocCPUS\n21671_0|21672|COMPLETED|0:0|00:10:00|bizon|32\n"
    status = {"status": "finished_pending_native_validation", "failed_methods": [],
              "accuracy_evaluated": False, "native_outputs_validated": False,
              "provenance": {"slurm_job_id": "21672", "slurm_array_task_id": "0"}}
    postflight = {"status": "inputs_runtime_sources_reverified", "failed_methods": [],
                  "native_outputs_validated": False, "accuracy_evaluated": False}
    return accounting, status, postflight


def test_terminal_identity_and_postflight():
    accounting, status, postflight = fixture()
    assert require_success(accounting, 0, status, postflight)["JobIDRaw"] == "21672"


@pytest.mark.parametrize("old,new", [("COMPLETED", "RUNNING"), ("0:0", "1:0"), ("bizon", "spark-7ff0"),
                                   ("|32", "|20"), ("21671_0", "21671_1"), ("21672", "99999")])
def test_wrong_scheduler_rejected(old, new):
    accounting, status, postflight = fixture()
    with pytest.raises(ValueError):
        require_success(accounting.replace(old, new), 0, status, postflight)


def test_duplicate_task_and_incomplete_postflight_rejected():
    accounting, status, postflight = fixture()
    with pytest.raises(ValueError):
        require_success(accounting + accounting.splitlines()[1] + "\n", 0, status, postflight)
    for key, value in (("status", "running"), ("failed_methods", ["cell"]), ("accuracy_evaluated", True)):
        changed = deepcopy(postflight)
        changed[key] = value
        with pytest.raises(ValueError):
            require_success(accounting, 0, status, changed)


def pairs(tmp_path, body):
    path = tmp_path / "pairs.tsv"
    path.write_text("gene_a\tspecies_a\tgene_b\tspecies_b\n" + body)
    return path


def test_native_pairs_not_group_cliques(tmp_path):
    path = pairs(tmp_path, "A\tsp1\tB\tsp2\nA\tsp1\tC\tsp3\n")
    owners, candidates = {"A": "sp1", "B": "sp2", "C": "sp3"}, dict.fromkeys("ABC", 0)
    assert check_pairs(path, owners, candidates, 2) == 2
    # All three genes share a candidate, but B-C must not be invented.
    with pytest.raises(ValueError, match="count"):
        check_pairs(path, owners, candidates, 3)


@pytest.mark.parametrize("body", ["B\tsp2\tA\tsp1\n", "A\tsp1\tB\tsp2\nA\tsp1\tB\tsp2\n",
                                  "A\twrong\tB\tsp2\n", "A\tsp1\tD\tsp4\n",
                                  "A\tsp1\tC\tsp1\n"])
def test_invalid_native_pairs_rejected(tmp_path, body):
    with pytest.raises(ValueError):
        check_pairs(pairs(tmp_path, body), {"A": "sp1", "B": "sp2", "C": "sp1"}, dict.fromkeys("ABC", 0), 1)


def test_pairs_cannot_cross_candidates(tmp_path):
    with pytest.raises(ValueError, match="crosses"):
        check_pairs(pairs(tmp_path, "A\tsp1\tB\tsp2\n"), {"A": "sp1", "B": "sp2"}, {"A": 0, "B": 1}, 1)


@pytest.mark.parametrize("revision", ["b66225dc5fc355702575aebe34b644916894f236", "49ab110358c0b4c73806a640de9068494a311f63"])
def test_revision_is_explicit_not_relaxed(revision):
    prepared = {"launcher": {"path": "/frozen/replay.py", "sha256": "source", "bytes": 1}}
    metrics = {"git": {"commit": revision, "dirty": False}, "source": prepared["launcher"], "cwd": "/frozen"}
    check_native_origin(metrics, prepared, Path("/frozen"), revision)
    for key, value in (("git", {"commit": revision, "dirty": True}), ("git", {"commit": "other", "dirty": False}),
                       ("cwd", "/development"), ("source", {"path": "/other"})):
        changed = {**metrics, key: value}
        with pytest.raises(ValueError):
            check_native_origin(changed, prepared, Path("/frozen"), revision)
