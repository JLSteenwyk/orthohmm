import copy
import io
from pathlib import Path

from Bio import Phylo
import pytest

from benchmark_tools import validate_species_tree_perturbations as validator


def accounting():
    return "JobID|JobIDRaw|State|ExitCode|Elapsed\n" + "".join(
        f"21299_{i}|{21300+i}|COMPLETED|0:0|00:03:00\n" for i in range(6))


@pytest.mark.parametrize("change", [None, "running", "failed", "missing", "duplicate"])
def test_all_six_scheduler_tasks_required(change):
    text = accounting()
    if change == "running":
        text = text.replace("COMPLETED", "RUNNING", 1)
    elif change == "failed":
        text = text.replace("0:0", "1:0", 1)
    elif change == "missing":
        text = "\n".join(text.splitlines()[:-1])
    elif change == "duplicate":
        text += text.splitlines()[1] + "\n"
    if change:
        with pytest.raises(ValueError):
            validator.completed_tasks(text)
    else:
        assert set(validator.completed_tasks(text)) == set(range(6))


@pytest.mark.parametrize("change", [None, "index", "tree", "job", "array", "postflight", "scored", "control"])
def test_execution_identity_guards(change):
    cell, tree = {"label": "p1_c1_r1_nni1_0"}, {"label": "nni1_0"}
    preflight = {"job_id": "21300", "array_job_id": "21299", "array_task_id": "0",
                 "perturbation_index": 0, "tree_variant": tree, "cell": cell,
                 "cwd": "/launcher", "control_admission": {"status": "equivalent"}}
    status = {"status": "finished_pending_native_validation", "failed_methods": [],
              "accuracy_evaluated": False, "native_outputs_validated": False,
              "dataset": cell["label"], "methods": {cell["label"]: {}}, "provenance": preflight}
    postflight = {"status": "complete_pending_native_validation", "accuracy_evaluated": False,
                  "baseline_source_unchanged": True, "cell": cell}
    if change == "index":
        preflight["perturbation_index"] = 1
    elif change == "tree":
        preflight["tree_variant"] = {"label": "nni2_0"}
    elif change == "job":
        preflight["job_id"] = "21299"
    elif change == "array":
        preflight["array_task_id"] = "1"
    elif change == "postflight":
        postflight["baseline_source_unchanged"] = False
    elif change == "scored":
        status["accuracy_evaluated"] = True
    elif change == "control":
        preflight["control_admission"]["status"] = "not_equivalent"
    args = status, preflight, postflight, cell, tree, 0, {"JobIDRaw": "21300"}, Path("/launcher")
    if change:
        with pytest.raises(ValueError):
            validator.check_status(*args)
    else:
        validator.check_status(*args)


@pytest.mark.parametrize("change", [None, "distance", "observed", "taxa"])
def test_actual_tree_topology_checked(change):
    source = Phylo.read(io.StringIO("((a,b),(c,d));"), "newick")
    supplied = Phylo.read(io.StringIO("(a,(b,(c,d)));"), "newick")
    observed = copy.deepcopy(supplied)
    distance = 2
    if change == "distance":
        distance = 4
    elif change == "observed":
        observed = source
    elif change == "taxa":
        observed.get_terminals()[0].name = "x"
    if change:
        with pytest.raises(ValueError):
            validator.check_topology(source, supplied, observed, distance)
    else:
        assert validator.check_topology(source, supplied, observed, distance) == 2


def test_live_panel_cannot_reach_input_admission(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(validator.subprocess, "check_output", lambda *a, **k: accounting().replace("COMPLETED", "RUNNING", 1))
    monkeypatch.setattr(validator, "read_frozen", lambda *a: pytest.fail("Passed live scheduler gate"))
    with pytest.raises(ValueError):
        validator.validate(tmp_path)
