from copy import deepcopy
import json
import math

import numpy as np
import pytest

from benchmark_tools.probe_cap400_final_graph import (
    REVISION, gene_table, reconstruct, require_failure, run, verify_rows,
)


def inputs(tmp_path, rows):
    names = tmp_path / "names.txt"
    names.write_text("z\na\nisolate\n")
    abc = tmp_path / "edges.txt"
    abc.write_text(rows)
    return abc, names, tmp_path / "payload"


def test_exact_order_duplicates_orientation_and_isolate(tmp_path):
    weights = [0.0, -0.0, math.nextafter(1.0, 2.0), 5e-324, 1.7976931348623157e308]
    abc, names, payload = inputs(tmp_path, "".join(f"a\tz\t{w}\n" for w in weights))
    result = reconstruct(abc, names, payload, len(weights))
    assert result["genes"] == 3 and result["edges"] == 5
    assert np.load(payload / "sources.npy").tolist() == [1] * 5
    assert np.load(payload / "targets.npy").tolist() == [0] * 5
    assert np.load(payload / "weights.npy").tobytes() == np.array(weights, dtype="float64").tobytes()
    assert verify_rows(abc, names, payload, 5)["status"] == "all_saved_rows_and_float64_bits_verified"
    with pytest.raises(FileExistsError):
        reconstruct(abc, names, payload, 5)


@pytest.mark.parametrize("text", ["a\tz\tnan\n", "a\tz\tinf\n", "a\tz\t1\n",
    "a\tmissing\t1.0\n", "a z 1.0\n", "a\tz\t1.0", "a\tz\t1.0\textra\n"])
def test_reject_invalid_edge_rows(tmp_path, text):
    abc, names, payload = inputs(tmp_path, text)
    with pytest.raises(ValueError):
        reconstruct(abc, names, payload, 1)


@pytest.mark.parametrize("count", [True, 0, -1, 1.0, 2])
def test_reject_wrong_edge_count(tmp_path, count):
    abc, names, payload = inputs(tmp_path, "a\tz\t1.0\n")
    with pytest.raises(ValueError):
        reconstruct(abc, names, payload, count)


def test_reject_excess_rows(tmp_path):
    abc, names, payload = inputs(tmp_path, "a\tz\t1.0\na\tz\t1.0\n")
    with pytest.raises(ValueError):
        reconstruct(abc, names, payload, 1)


@pytest.mark.parametrize("field,value", [("sources", 0), ("targets", 2), ("sources", -1),
                                        ("weights", 2.0), ("weights", float("nan"))])
def test_independent_readback_rejects_mutations(tmp_path, field, value):
    abc, names, payload = inputs(tmp_path, "a\tz\t1.0\n")
    reconstruct(abc, names, payload, 1)
    array = np.load(payload / (field + ".npy"), mmap_mode="r+")
    array[0] = value
    array.flush()
    with pytest.raises(ValueError):
        verify_rows(abc, names, payload, 1)


@pytest.mark.parametrize("text", ["", "a\na\n", "a b\n", "a\n\n"])
def test_invalid_complete_gene_table(tmp_path, text):
    path = tmp_path / "names.txt"
    path.write_text(text)
    with pytest.raises(ValueError):
        gene_table(path)


def test_independent_readback_rejects_signed_zero_change(tmp_path):
    abc, names, payload = inputs(tmp_path, "a\tz\t-0.0\n")
    reconstruct(abc, names, payload, 1)
    weights = np.load(payload / "weights.npy", mmap_mode="r+")
    weights[0] = 0.0
    weights.flush()
    with pytest.raises(ValueError):
        verify_rows(abc, names, payload, 1)


def test_independent_readback_rejects_gene_order_change(tmp_path):
    abc, names, payload = inputs(tmp_path, "a\tz\t1.0\n")
    reconstruct(abc, names, payload, 1)
    (payload / "gene_names.txt").write_text("a\nz\nisolate\n")
    with pytest.raises(ValueError):
        verify_rows(abc, names, payload, 1)


def failure_fixture():
    scheduler = {"JobIDRaw": "24245", "State": "FAILED", "ExitCode": "1:0",
                 "Elapsed": "2-04:30:54", "NodeList": "bizon", "AllocCPUS": "32"}
    plan = {"scientific_revision": REVISION, "case": {"fresh_case": {"command": ["frozen"]}}}
    parent = {"job_id": "24245", "status": "failed_or_unverified_control_continuation",
              "scores": {}, "execution": {"returncode": 1, "command": ["frozen"]}}
    metrics = {"status": "failed", "harness": {"exit_code": 1, "git_commit": REVISION,
               "git_dirty": False, "output_line_counts": {"edges": 33585180}}}
    receipt = {"status": "terminal_cap400_final_leiden_failure_and_retained_checkpoint_verified",
               "optimizer_attempted": False, "scores": {}, "scheduler": deepcopy(scheduler),
               "final_graph_edges_reported_by_failed_harness": 33585180}
    return [plan, parent, metrics, receipt, scheduler]


def test_accept_exact_authenticated_failure():
    require_failure(*failure_fixture())


@pytest.mark.parametrize("index,key,value", [
    (0, "scientific_revision", "other"), (1, "scores", {"F1": 1.0}),
    (1, "status", "running"), (1, "job_id", "24247"),
    (2, "status", "complete"), (3, "optimizer_attempted", True),
    (3, "final_graph_edges_reported_by_failed_harness", 1),
    (4, "State", "RUNNING"), (4, "ExitCode", "0:0"), (4, "AllocCPUS", "1"),
])
def test_reject_changed_failure_gate(index, key, value):
    arguments = failure_fixture()
    arguments[index][key] = value
    with pytest.raises(ValueError):
        require_failure(*arguments)


def test_attempted_namespace_never_reexecutes_optimizer(tmp_path, monkeypatch):
    from pathlib import Path
    import benchmark_tools.probe_cap400_final_graph as probe
    plan = {"source_root": str(Path(probe.__file__).resolve().parent.parent),
            "repo": str(tmp_path), "output": str(tmp_path)}
    path = tmp_path / "execution_plan.json"
    path.write_text(json.dumps(plan))
    (tmp_path / "report.json").write_text("{}")
    def unexpected(*args, **kwargs):
        pytest.fail("Repeated namespace must not invoke a subprocess")
    monkeypatch.setattr(probe.subprocess, "run", unexpected)
    with pytest.raises(ValueError, match="no optimizer retries"):
        run(path)
