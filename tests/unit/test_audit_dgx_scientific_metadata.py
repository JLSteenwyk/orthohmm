from copy import deepcopy
import json
from pathlib import Path

import pytest

from benchmark_tools.audit_dgx_scientific_metadata import WRAPPER_SHA, runtime_rows, time_command, verify_contract


def fixture(index):
    path = Path(__file__).resolve().parents[2] / "benchmark_tools/results/dgx_scientific_execution_20260917.json"
    spec = json.loads(path.read_text())
    run = spec["runs"][index]
    order = next(r for r in spec["orders"] if r["input_directory"] == run["dataset"]["input_directory"])
    expected = deepcopy(run)
    expected["gnu_time"] = dict(executable="/usr/bin/time", output=str(Path(run["measurement_directory"]).parent / "native.time.tsv"))
    argv = time_command(run["native_argv"], expected["gnu_time"]["output"])
    copies = []
    names = order["native_order"]
    if run["native_method"] == "orthofinder_full":
        copies = sorted([{**r, "path": str(Path(run["configuration"]["copy_inputs_to"]) / Path(r["path"]).name)}
                         for r in run["dataset"]["inputs"]], key=lambda r: r["path"])
        names = sorted(names)
    prepared = dict(run=expected, measured_argv=argv, status="fresh_native_inputs_prepared",
                    inference_started=False, original_inputs=deepcopy(order["inputs_in_native_order"]),
                    copied_inputs=copies, expected_native_basename_order=names, preparation_wall_s=.1)
    measured = dict(command=deepcopy(argv), cwd=run["cwd"], job_id=9000, requested_cpus=20,
                    requested_memory_bytes=96*1024**3, timeout_s=85800, interval_s=1, host_interval_s=30,
                    clock_domain=dict(hostname="spark-7ff0"), status="command_exited_zero", exit_code=0, timed_out=False)
    before = dict(native_order=deepcopy(order["native_order"]), original_inputs=deepcopy(order["inputs_in_native_order"]),
                  runtime=runtime_rows(spec))
    after = deepcopy(before)
    if copies:
        after["copied_inputs"] = deepcopy(copies)
    verified = dict(status="command_exited_zero", measurement=deepcopy(measured), source_sha256=WRAPPER_SHA,
                    before=before, after=after, before_check_wall_s=1, after_check_wall_s=1)
    scheduler = [f"21656_{index}", "9000", "COMPLETED", "0:0", "20", "96G", "spark-7ff0"]
    return spec, run, prepared, verified, measured, scheduler


@pytest.mark.parametrize("index", [0, 1, 2])
def test_each_native_method_contract(index):
    verify_contract(*fixture(index))


@pytest.mark.parametrize("mutation", ["scheduler", "job", "argv", "cwd", "runtime", "order", "input", "cpu", "timeout", "timed_out", "duration", "copies"])
def test_changed_scientific_contract_rejected(mutation):
    spec, run, prepared, verified, measured, scheduler = fixture(2)
    if mutation == "scheduler":
        scheduler[2] = "RUNNING"
    elif mutation == "job":
        scheduler[1] = "9001"
    elif mutation == "argv":
        measured["command"].append("--unexpected")
    elif mutation == "cwd":
        measured["cwd"] = "/tmp"
    elif mutation == "runtime":
        verified["after"]["runtime"][2]["sha256"] = "0"*64
    elif mutation == "order":
        verified["after"]["native_order"].reverse()
    elif mutation == "input":
        verified["before"]["original_inputs"][0]["sha256"] = "0"*64
    elif mutation in ("cpu", "timeout", "timed_out"):
        key, value = {"cpu": ("requested_cpus", 32), "timeout": ("timeout_s", 900), "timed_out": ("timed_out", True)}[mutation]
        measured[key] = value
        verified["measurement"] = deepcopy(measured)
    elif mutation == "duration":
        prepared["preparation_wall_s"] = float("nan")
    else:
        verified["after"]["copied_inputs"].pop()
    with pytest.raises(ValueError):
        verify_contract(spec, run, prepared, verified, measured, scheduler)
