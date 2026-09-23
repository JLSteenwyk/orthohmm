"""Bind prespecified CPM replay settings without altering the baseline plan."""

from copy import deepcopy
import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent))
from prepare_ob_candidate_neighborhood import check, record

PLAN_SHA = "4f25aff958e327a3d133f08326260c600e51e87e1d56c81f114db9993930479b"
PROTOCOL_SHA = "7b66c2f29b0098f26a68cb9f230bd6260e04089db9d40c98f4ce0b9af1924044"
REPLAY_SHA = "e7657732bc94438fb0602f3246680308a075b9fda1c68909ddfac1138b164c0c"
RESOLUTIONS = {"control": .1, "cpm_low": .08, "cpm_high": .12}


def read_frozen(path, digest):
    item = record(path)
    if item["sha256"] != digest:
        raise ValueError("Frozen CPM input changed")
    result = json.loads(path.read_text())
    check(item)
    return result


def derive(root, baseline, arm):
    if arm not in RESOLUTIONS:
        raise ValueError("Unknown prespecified CPM arm")
    launcher = root / "benchmarks/work/publication_qfo_replay_native_v1"
    original_output = root / "benchmarks/results/qfo_corrected_checked_replay_v1"
    directories = {Path(r["path"]).parent for r in baseline["input_fastas"]}
    if len(directories) != 1:
        raise ValueError("Ambiguous corrected FASTA directory")
    checkpoint = baseline["checkpoint_manifest"]
    expected = [sys.executable, str(launcher / "benchmark_tools/replay_high_sensitivity.py"),
        "--accuracy-checkpoint", str(Path(checkpoint["path"]).parent), "--checkpoint-sha256", checkpoint["sha256"],
        "--fasta-directory", str(next(iter(directories))), "--output-directory", str(original_output / "replay"),
        "--json", str(original_output / "replay.json"), "--cpu", "32", "--matrix", "BLOSUM62",
        "--cpm-resolution", "0.1", "--leiden-seed", "4", "--profile-iterations", "1", "--profile-min-species", "1"]
    environment = {"PYTHONPATH": str(launcher), "PYTHONHASHSEED": "0", "OMP_NUM_THREADS": "1",
                   "OPENBLAS_NUM_THREADS": "1", "MKL_NUM_THREADS": "1"}
    if (baseline["native_command"] != expected or baseline["output_root"] != str(original_output)
            or baseline["cwd"] != str(launcher) or baseline["environment_overrides"] != environment
            or baseline["expected_stages"] != ["multipass", "multipass_refined", "strict_profiles", "strict_profiles_refined"]):
        raise ValueError("Baseline command, stages or environment differs")
    # Preserve admitted control and failed v1 variants; replacement variants are fresh.
    version = "v1" if arm == "control" else "v2"
    output = root / f"benchmarks/results/qfo_parameter_cpm_replay_{version}" / arm
    command = list(expected)
    for flag, value in (("--cpm-resolution", str(RESOLUTIONS[arm])),
                        ("--output-directory", str(output / "replay")), ("--json", str(output / "replay.json"))):
        command[command.index(flag) + 1] = value
    return {"arm": arm, "resolution": RESOLUTIONS[arm], "output_root": str(output),
            "native_command": command, "cwd": str(launcher), "environment_overrides": environment,
            "expected_stages": list(baseline["expected_stages"]), "baseline_plan": deepcopy(baseline),
            "metadata": {"cpm_resolution": RESOLUTIONS[arm], "seed": 4, "include_isolates": True,
                         "output_directory": str(output / "replay")}}


def evidence(root, baseline, baseline_record, arm):
    if (baseline_record["sha256"] != REPLAY_SHA
            or baseline_record["path"] != str(root / "benchmarks/work/qfo_corrected_replay_commands_20260918.json")):
        raise ValueError("CPM context requires the pinned corrected replay plan")
    check(baseline_record)
    pinned = read_frozen(Path(baseline_record["path"]), REPLAY_SHA)
    if baseline != pinned:
        raise ValueError("Supplied baseline differs from frozen file")
    results = root / "benchmark_tools/results"
    plan_path = results / "qfo_parameter_neighborhood_plan_20260919.json"
    protocol_path = results / "QFO_PARAMETER_NEIGHBORHOOD_PROTOCOL_20260919.md"
    plan = read_frozen(plan_path, PLAN_SHA)
    protocol = record(protocol_path)
    if protocol["sha256"] != PROTOCOL_SHA:
        raise ValueError("CPM protocol changed")
    expected_arms = [{"label": "control", "delta": {}},
                     {"label": "cpm_low", "delta": {"cpm_resolution": .08}},
                     {"label": "cpm_high", "delta": {"cpm_resolution": .12}}]
    if plan["arms"][:3] != expected_arms or plan["baseline_parameters"]["cpm_resolution"] != .1:
        raise ValueError("CPM neighborhood differs from frozen settings")
    inputs = [baseline_record, record(plan_path), protocol]
    for item in plan["inputs"]:
        actual = record(root / item["path"])
        if actual["sha256"] != item["sha256"]:
            raise ValueError("Parameter baseline input changed")
        inputs.append(actual)
    result = derive(root, baseline, arm)
    result["checked_records"] = inputs
    for item in inputs:
        check(item)
    return result


def settings(manifest, metadata, context):
    if metadata != context["metadata"] or manifest["output_directory"] != context["metadata"]["output_directory"]:
        raise ValueError("CPM payload settings differ from the frozen arm")
