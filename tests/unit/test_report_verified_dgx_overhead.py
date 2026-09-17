from copy import deepcopy

import pytest

from benchmark_tools.report_verified_dgx_overhead import PROJECT, RUNTIMES, check_verification, completed_jobs


def accounting():
    return "\n".join(f"21647_{i}|COMPLETED|0:0|00:01:40|20|96G|spark-7ff0" for i in range(6))


def test_complete_scheduler_inventory():
    assert len(completed_jobs(accounting())) == 6


@pytest.mark.parametrize("bad", [lambda x: x.replace("COMPLETED", "RUNNING", 1),
                                 lambda x: x.replace("0:0", "1:0", 1),
                                 lambda x: x.replace("96G", "48G", 1),
                                 lambda x: "\n".join(x.splitlines()[:-1]),
                                 lambda x: x + "\n" + x.splitlines()[0]])
def test_incomplete_failed_changed_or_duplicate_scheduler_rejected(bad):
    with pytest.raises(ValueError):
        completed_jobs(bad(accounting()))


def verification():
    rows = [{"path": PROJECT + "/runtime_inventory_v1/" + name, "sha256": sha,
             "records": count, "scientific_execution_authorized": False,
             "status": "runtime_tree_identity_matches"} for name, sha, count in RUNTIMES]
    return {"status": "command_exited_zero", "before": rows, "after": deepcopy(rows),
            "measurement": {"command": ["load"]}, "before_check_wall_s": 1., "after_check_wall_s": 1.,
            "source_sha256": "36243fcdfe6f49e73dbe1dd1b627bf4a330f19f3cd52dcf94a9618ddf0415cf8"}


def test_exact_before_after_and_embedded_record():
    v = verification()
    check_verification(v, {"command": ["load"]})
    with pytest.raises(ValueError):
        check_verification(v, {"command": ["other"]})


@pytest.mark.parametrize("mutation", [lambda v: v["after"].pop(),
    lambda v: v["after"][0].update(sha256="0" * 64),
    lambda v: v.update(after_check_wall_s=float("nan")),
    lambda v: v.update(source_sha256="0" * 64)])
def test_changed_or_missing_identity_rejected(mutation):
    v = verification()
    mutation(v)
    with pytest.raises(ValueError):
        check_verification(v, v["measurement"])
