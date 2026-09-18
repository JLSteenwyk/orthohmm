from copy import deepcopy

import pytest

from benchmark_tools.admit_qfo_corrected_proteinortho import validate_execution


def fixture():
    def rec(path):
        return {"path": path, "bytes": 10, "sha256": "abc"}
    plan = {"cwd": "/run/input", "output_root": "/run", "native_argv": ["tool"],
            "resources": {"node": "bizon"}, "input_fastas": [rec("/data/s.fasta")],
            "native_outputs": {"pairs": "input/qfo.proteinortho-graph", "groups": "input/qfo.proteinortho.tsv"}}
    plan_record, runtime_record = rec("/plan"), rec("/runtime")
    runtime = {"status": "proteinortho_runtime_frozen_unrun", "plan": plan_record, "source": rec("/runner")}
    execution = {"status": "process_succeeded_pending_native_admission", "exit_code": 0,
                 "plan": plan_record, "runtime": runtime_record, "source": rec("/frozen/runner"),
                 "native_argv": ["tool"], "cwd": "/run/input", "node": "bizon", "job_id": "1",
                 "started_epoch": 1, "finished_epoch": 2, "copied_inputs": [rec("/run/input/s.fasta")],
                 "log": rec("/run/native.log"), "timing": rec("/run/time.txt"),
                 "outputs": [rec("/run/input/s.fasta"), rec("/run/input/qfo.proteinortho-graph"),
                             rec("/run/input/qfo.proteinortho.tsv")]}
    return deepcopy((plan, runtime, execution, plan_record, runtime_record))


def test_valid_binding():
    assert set(validate_execution(*fixture())) == {"pairs", "groups"}


@pytest.mark.parametrize("key,value", [
    ("status", "running"), ("status", "failed"), ("exit_code", 1),
    ("native_argv", ["different"]), ("cwd", "/old/input"), ("node", "other"),
    ("job_id", ""), ("finished_epoch", 0), ("copied_inputs", []), ("outputs", []),
    ("plan", {}), ("runtime", {}), ("source", {"sha256": "different"}),
])
def test_reject_execution_changes(key, value):
    args = fixture()
    args[2][key] = value
    with pytest.raises(ValueError):
        validate_execution(*args)


def test_reject_duplicate_output():
    args = fixture()
    args[2]["outputs"].append(args[2]["outputs"][0])
    with pytest.raises(ValueError, match="inventory"):
        validate_execution(*args)


def test_reject_empty_graph():
    args = fixture()
    args[2]["outputs"][1]["bytes"] = 0
    with pytest.raises(ValueError, match="Missing/empty"):
        validate_execution(*args)


def test_reject_runtime_binding():
    args = fixture()
    args[1]["plan"] = {}
    with pytest.raises(ValueError, match="Runtime plan"):
        validate_execution(*args)
