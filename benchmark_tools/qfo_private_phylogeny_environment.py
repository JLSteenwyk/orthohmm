"""Explicit private deployment verification for QfO phylogeny, not timing."""

from copy import deepcopy
import json
import os
from pathlib import Path
import subprocess
import sys

from benchmark_tools.prepare_ob_candidate_neighborhood import check, record

DEPLOYMENT = "benchmark_tools/results/threadripper_private_deployment_20260928.json"
DEPLOYMENT_SHA = "7225658b2152612bc9e0c7e83c7eabb859d99f3e04cec49bca561fdf58c38336"
TREES = "benchmark_tools/results/threadripper_private_trees_20260928.json"
TREES_SHA = "952e45f34ffd6bb679436b592a9b41672ae1006c5966f93cd930722eabe6f973"
LOOKUP = "benchmark_tools/results/threadripper_private_lookup_v2_20260928.json"
LOOKUP_SHA = "5996f36dad39c7a7e38c38f134cbe5e13796c4a23dae0444a77a473f12cfdb54"
BASELINE_SHA = "b39e5f26518191ad53e7a796690f083a0bd5ede2081cbfbcb7b00a86947acdc3"
LEGACY_SHA = "bf677728cf9312edd9755d9ac1ceefbce3abb8c3848d43413631dc00ca00eb2f"
NATIVE_ADMISSION_SHA = "c8a289ac8128711da6c4e93654ce6953e48854b8d21b2adfea074d3feb8c70aa"


def read_pin(path, sha):
    item = record(path)
    if item["sha256"] != sha:
        raise ValueError("Private QfO evidence identity changed: " + str(path))
    return json.loads(Path(item["path"]).read_bytes()), item


def deployment_difference(legacy, ancestral, amended, deployment):
    expected = deepcopy(legacy)
    if expected["environments"]["orthofinder"]["packages"].pop("orthohmm", None) != "0.5.0":
        raise ValueError("Unexpected historical repository-local metadata")
    expected["threadripper_amendment"] = ancestral["threadripper_amendment"]
    if ancestral != expected:
        raise ValueError("Private deployment ancestor changed scientific/native settings")
    expected = deepcopy(ancestral)
    python = deployment["new_interpreter"]
    expected["tool_entrypoints"]["orthohmm_python"] = dict(path="python", absolute_path=python["path"],
        bytes=python["bytes"], sha256=python["sha256"])
    expected["environments"]["orthohmm"] = amended["environments"]["orthohmm"]
    if amended != expected:
        raise ValueError("Private deployment changed more than interpreter/package inventory")


def verify_private_environment(root):
    from benchmark_tools.snapshot_runtime_trees import verify as verify_tree
    from benchmark_tools.run_simulation_methods import verify_environment
    from benchmark_tools.build_private_timing_environment import canonical

    deployment, deployment_pin = read_pin(root / DEPLOYMENT, DEPLOYMENT_SHA)
    private, private_pin = read_pin(Path(deployment["baseline"]["path"]), BASELINE_SHA)
    ancestral, ancestral_pin = read_pin(Path(deployment["previous_baseline"]["path"]),
                                        deployment["previous_baseline"]["sha256"])
    legacy, legacy_pin = read_pin(root / "benchmark_tools/results/publication_variable_native_methods_20260916.json", LEGACY_SHA)
    deployment_difference(legacy, ancestral, private, deployment)
    candidate, candidate_pin = read_pin(Path(deployment["candidate"]["path"]), deployment["candidate"]["sha256"])
    fixture, fixture_pin = read_pin(Path(deployment["fixture"]["path"]), deployment["fixture"]["sha256"])
    lookup, lookup_pin = read_pin(root / LOOKUP, LOOKUP_SHA)
    trees, trees_pin = read_pin(root / TREES, TREES_SHA)
    if (deployment["status"] != "prospective_private_deployment_prepared"
            or fixture["status"] != "both_native_prediction_fixtures_match" or fixture["candidate"] != candidate_pin
            or candidate["status"] != "private_timing_environment_candidate_installed"
            or lookup["status"] != "native_lookup_repeated_identity_match" or lookup["baseline"] != private_pin
            or {canonical(k): v for k, v in private["environments"]["orthohmm"]["packages"].items()}
                != {canonical(k): v for k, v in candidate["selected"].items()}):
        raise ValueError("Private deployment lacks matching installation/fixture/lookup evidence")
    inventory, inventory_pin = read_pin(Path(trees["inventory"]["path"]), trees["inventory"]["sha256"])
    tree_check = verify_tree(inventory)
    records = [deployment_pin, private_pin, ancestral_pin, legacy_pin, candidate_pin, fixture_pin,
               lookup_pin, trees_pin, inventory_pin, deployment["new_interpreter"], candidate["lock"], candidate["import_report"],
               *fixture["checked_records"]]
    for item in records:
        check(item)
    for item in lookup["interpreters"]["orthohmm"]["reports"]:
        check(item)
        records.append(item)
    # The frozen private inventory was captured from core cwd, without repo egg-info.
    cwd = Path.cwd()
    try:
        os.chdir(private["core_root"])
        verify_environment(private)
    finally:
        os.chdir(cwd)
    return private, {"deployment": deployment_pin, "private_baseline": private_pin,
                     "tree_check": tree_check, "checked_records": records}


def verify_baseline(root):
    """Keep historical admission/source checks; replace only deployment verification."""
    from benchmark_tools.run_qfo_corrected_factorial_cell import verify_admission
    from benchmark_tools.run_qfo_factorial_cell import native_command
    from benchmark_tools.validate_simulation_outputs import verify_process

    path = root / "benchmark_tools/results/qfo_corrected_factorial_native_admission_21764.json"
    baseline, baseline_pin = read_pin(path, NATIVE_ADMISSION_SHA)
    if baseline["status"] != "corrected_qfo_native_pair_output_verified" or baseline["cell"] != "p1_c1_r1":
        raise ValueError("Require complete corrected full-pipeline baseline")
    for item in baseline["checked_records"]:
        check(item)
    old = baseline["candidate_admission"]
    _, manifest, original, _, prepared, _ = verify_admission(root, Path(old["path"]), old["sha256"], "21758", 3)
    launcher = Path(manifest["launcher_root"])
    argv, equivalence = native_command(original, launcher, prepared)
    status_pin = baseline["native_group_integrity"]["integrity"]["execution_status"]
    check(status_pin)
    status = json.loads(Path(status_pin["path"]).read_bytes())
    config = {"argv": argv, "output": argv[argv.index("--output-directory") + 1],
              "metrics": argv[argv.index("--json") + 1]}
    verify_process(config, status["methods"][original["label"]])
    environment, runtime = verify_private_environment(root)
    return {"manifest": manifest, "original": original, "launcher": str(launcher), "prepared": str(prepared),
            "environment": environment, "private_runtime": runtime, "original_argv": argv,
            "original_launcher_equivalence": equivalence,
            "checked_records": [baseline_pin, status_pin, *baseline["checked_records"], *runtime["checked_records"]]}


def control_command(original, launcher, prepared, environment, output):
    from benchmark_tools.run_qfo_factorial_cell import native_command

    argv, equivalence = native_command(original, launcher, prepared)
    if (original["label"] != "p1_c1_r1" or "--species-tree" in argv or "--checkpoint-source" in argv
            or "--official-benchmark" in argv or any(argv.count(flag) != 1 for flag in
                ("--output-directory", "--json", "--species-tree-mode", "--cpu"))
            or argv[argv.index("--species-tree-mode") + 1] != "infer"
            or argv[argv.index("--cpu") + 1] != "32"):
        raise ValueError("Unexpected inferred full-pipeline baseline command")
    checkpoint = argv[argv.index("--output-directory") + 1]
    argv[0] = environment["tool_entrypoints"]["orthohmm_python"]["absolute_path"]
    argv[argv.index("--output-directory") + 1] = str(output / "output")
    argv[argv.index("--json") + 1] = str(output / "metrics.json")
    argv.extend(["--checkpoint-source", checkpoint])
    return argv, equivalence


def inspect_launcher(verified, env, output):
    from benchmark_tools.inspect_native_python_lookup import PROBE, scientific_origins

    requested = ["orthohmm.phylogeny_pipeline", "benchmark_tools.replay_phylogeny"]
    python = verified["environment"]["tool_entrypoints"]["orthohmm_python"]["absolute_path"]
    command = [python, "-B", "-c", PROBE, json.dumps(requested)]
    process = subprocess.run(command, cwd=verified["launcher"], env=env, capture_output=True, text=True, timeout=180)
    (output / "lookup_process.json").write_text(json.dumps(dict(command=command, returncode=process.returncode,
        stdout=process.stdout, stderr=process.stderr), indent=2, sort_keys=True) + "\n")
    if process.returncode:
        raise RuntimeError("Private QfO launcher import probe failed")
    report = json.loads(process.stdout)
    origin = scientific_origins(report, "orthohmm", Path(verified["launcher"]) / "orthohmm")
    if report["modules"].get("benchmark_tools.replay_phylogeny") != str(Path(verified["launcher"]) / "benchmark_tools/replay_phylogeny.py"):
        raise ValueError("Wrong frozen QfO replay launcher import")
    paths = sorted(set(report["modules"].values()) | set(report["mapped_files"]) | {report["executable"]})
    if any(path.startswith(("/home/bizon/anaconda3/", "/home/bizon/.local/")) for path in paths):
        raise ValueError("Private QfO launcher uses retired shared/user runtime")
    result = {"origin": origin, "report": report, "checked_records": [record(path) for path in paths],
              "continuous_enforcement": False}
    (output / "lookup.json").write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    return result
