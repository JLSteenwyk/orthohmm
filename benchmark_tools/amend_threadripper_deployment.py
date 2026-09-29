"""Prospective interpreter deployment amendment; never starts timing runs."""

import argparse
from copy import deepcopy
import json
from pathlib import Path
import subprocess

from benchmark_tools.build_private_timing_environment import canonical, record, write
from benchmark_tools.run_simulation_methods import read_frozen


def amend_plan(plan, old_python, new_python):
    amended = deepcopy(plan)
    rows = amended["runs"]
    methods = ("orthohmm_high_sensitivity", "orthohmm_satellite_v2", "orthofinder_full")
    expected = [(methods[(repeat + size_index + offset) % 3], size, repeat)
                for repeat in range(3) for size_index, size in enumerate((4, 8, 12))
                for offset in range(3)]
    if (len(rows) != 27 or [r["index"] for r in rows] != list(range(27))
            or [(r["native_method"], r["proteomes"], r["repeat"]) for r in rows] != expected
            or len({(r["native_method"], r["proteomes"], r["repeat"]) for r in rows}) != 27
            or {r["proteomes"] for r in rows} != {4, 8, 12}
            or {r["repeat"] for r in rows} != {0, 1, 2}):
        raise ValueError("Unexpected frozen panel identities")
    changed = []
    for row in rows:
        if row["native_method"] == "orthofinder_full":
            continue
        for argv in (row["native_argv"], row["configuration"]["argv"]):
            if argv[0] != old_python:
                raise ValueError("Unexpected original OrthoHMM interpreter")
            argv[0] = new_python
        changed.append(row["index"])
    return amended, changed


def run(repo, output):
    if output.exists():
        raise FileExistsError(output)
    directory = repo / "benchmark_tools/results"
    prior_baseline = directory / "threadripper_native_baseline_20260928.json"
    prior_plan = directory / "threadripper_scaling_commands_20260928.json"
    baseline = read_frozen(prior_baseline, "84ec3c10e2075ea564bef049e7cfc5c72aa878c1f050afbff43bff512805d041")
    plan = read_frozen(prior_plan, "c384e27730e3802b39ba14a42f7f50e84da5ce6deb9de9b2c32a74a745aed296")
    candidate_path = directory / "threadripper_patched_runtime_20260928.json"
    parity_path = directory / "threadripper_patched_predictions_20260928.json"
    candidate = json.loads(candidate_path.read_text())
    parity = json.loads(parity_path.read_text())
    if (candidate["status"] != "private_timing_environment_candidate_installed"
            or parity["status"] != "both_native_prediction_fixtures_match"
            or parity["candidate"] != record(candidate_path)
            or candidate["baseline"] != record(prior_baseline)):
        raise ValueError("Candidate lacks matching baseline and fixture evidence")
    for ref in (candidate["lock"], candidate["import_report"]):
        if record(Path(ref["path"])) != ref:
            raise ValueError("Candidate evidence changed")
    python = repo / "benchmarks/work/threadripper_patched_runtime_20260928/venv/bin/python"
    query = ("import importlib.metadata as m,json,sys; "
             "names=sorted({d.metadata['Name'] for d in m.distributions() if d.metadata['Name']}); "
             "print(json.dumps({'python':sys.version,'packages':{n:m.version(n) for n in names}}))")
    command = [str(python), "-I", "-B", "-c", query]
    observed = json.loads(subprocess.check_output(command, text=True, timeout=60))
    if ({canonical(k): v for k, v in observed["packages"].items()}
            != {canonical(k): v for k, v in candidate["selected"].items()}):
        raise ValueError("Candidate distribution drift")
    old_python = baseline["tool_entrypoints"]["orthohmm_python"]["absolute_path"]
    amended, indices = amend_plan(plan, old_python, str(python))
    amended_baseline = deepcopy(baseline)
    item = record(python)
    amended_baseline["tool_entrypoints"]["orthohmm_python"] = dict(
        path="python", absolute_path=str(python), bytes=item["bytes"], sha256=item["sha256"])
    amended_baseline["environments"]["orthohmm"] = observed
    output.mkdir(parents=True)
    write(output / "baseline.json", amended_baseline)
    write(output / "plan.json", amended)
    result = dict(status="prospective_private_deployment_prepared", source=record(Path(__file__)),
        previous_baseline=record(prior_baseline), previous_plan=record(prior_plan),
        candidate=record(candidate_path), fixture=record(parity_path),
        baseline=record(output / "baseline.json"), command_plan=record(output / "plan.json"),
        inventory_command=command, changed_run_indices=indices,
        old_interpreter=old_python, new_interpreter=record(python),
        scientific_execution_authorized=False, production_runs_launched=0,
        limitations=["Only native/configuration argv[0] changes; original_native_argv remains historical provenance.",
                     "Only OrthoHMM interpreter record and environment metadata change in the baseline.",
                     "Runtime-tree binding, repeated import checks, controller isolation and collector v5 validation remain required.",
                     "Fixture parity is not genome-scale equivalence or timing eligibility."])
    write(output / "amendment.json", result)
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    print(run(args.repo.resolve(), args.output.resolve())["status"])
