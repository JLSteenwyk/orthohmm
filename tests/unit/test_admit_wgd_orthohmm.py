from copy import deepcopy
import json
from pathlib import Path

import pytest

from benchmark_tools.admit_wgd_orthohmm import check_receipt, expected_harness
from benchmark_tools.prepare_wgd_commands import commands
from benchmark_tools.snapshot_orthohmm_input_order import record


def fixture(tmp_path):
    manifest = tmp_path / "runtime.json"
    manifest.write_text(json.dumps({"records": [{"path": "/binary"}]}))
    selected = commands(*map(Path, ("/core", "/input", "/out", "/python", "/of", "/sonic")))[0]
    spec = {"environment": {"PYTHONHASHSEED": "0"}, "runtime_manifests": [record(manifest)]}
    checks = [{"records": 1, "scientific_execution_authorized": False, "status": "runtime_tree_identity_matches"}]
    receipt = {"status": "native_exited_zero", "native": {"exit_code": 0, "timed_out": False},
               "method": deepcopy(selected), "hostname": "bizon", "environment": spec["environment"],
               "cwd": str(tmp_path), "executed_argv": ["/usr/bin/time", "-v", "-o", str(tmp_path / "native.time.log"), *selected["argv"]],
               "runtime_before": deepcopy(checks), "runtime_after": deepcopy(checks), "job_id": "999"}
    return receipt, selected, tmp_path, spec, "999|21661_0|COMPLETED|0:0\n", "21661_0"


def test_matching_scheduler_command_and_identity_pass(tmp_path):
    check_receipt(*fixture(tmp_path))


@pytest.mark.parametrize("mutation", ["exit", "timeout", "method", "command", "postflight", "job", "host"])
def test_invalid_receipt_rejected(tmp_path, mutation):
    args = fixture(tmp_path)
    receipt = args[0]
    if mutation == "exit":
        receipt["native"]["exit_code"] = 1
    elif mutation == "timeout":
        receipt["native"]["timed_out"] = True
    elif mutation == "method":
        receipt["method"]["argv"].append("--changed")
    elif mutation == "command":
        receipt["executed_argv"].append("--changed")
    elif mutation == "postflight":
        receipt["runtime_after"] = []
    elif mutation == "job":
        receipt["job_id"] = "998"
    else:
        receipt["hostname"] = "spark-7ff0"
    with pytest.raises(ValueError):
        check_receipt(*args)


def test_harness_phylogeny_flags_are_explicit():
    panel = commands(*map(Path, ("/core", "/input", "/out", "/python", "/of", "/sonic")))
    assert "--phylogeny" not in expected_harness(panel[0])
    command = expected_harness(panel[1])
    assert command[command.index("--phylogeny_candidates") + 1] == "satellite_v2"
    assert command[-2:] == ["--species_tree_rooting", "min_variance"]
