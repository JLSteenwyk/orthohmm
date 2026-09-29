"""Native prediction parity only; no collector or controlled timing admission."""

import argparse
import json
import os
from pathlib import Path
import shutil
import signal
import subprocess

from benchmark_tools.build_private_timing_environment import record, write
from benchmark_tools.probe_threadripper_native_fixture import fixture_run
from benchmark_tools.run_simulation_methods import read_frozen, execution_environment
from benchmark_tools.validate_scaling_outputs import validate


def check(ref):
    if record(Path(ref["path"])) != ref:
        raise ValueError("Evidence changed: " + ref["path"])


def parity(prior, output):
    rows = []
    for ref in prior["current"]["evidence"]:
        old = Path(ref["path"])
        if old.suffix == ".fa":
            check(ref)
            continue
        check(ref)
        relative = (Path("orthohmm_phylogeny") / old.name
                    if old.parent.name == "orthohmm_phylogeny" else Path(old.name))
        current = record(output / relative)
        identical = (current["sha256"], current["bytes"]) == (ref["sha256"], ref["bytes"])
        rows.append(dict(prior=ref, current=current, identical=identical))
    if not rows or not all(r["identical"] for r in rows):
        raise ValueError("Native prediction bytes differ or no predictions compared")
    return rows


def run(repo, output):
    if output.exists():
        raise FileExistsError(output)
    results = repo / "benchmark_tools/results"
    plan_path = results / "threadripper_scaling_commands_20260928.json"
    plan = read_frozen(plan_path, "c384e27730e3802b39ba14a42f7f50e84da5ce6deb9de9b2c32a74a745aed296")
    candidate_path = results / "threadripper_patched_runtime_20260928.json"
    candidate = json.loads(candidate_path.read_text())
    check(candidate["baseline"])
    check(candidate["import_report"])
    baseline = json.loads(Path(candidate["baseline"]["path"]).read_text())
    source = Path("/tmp/orthohmm-integrated-workflow-20260927-v2/input")
    python = repo / "benchmarks/work/threadripper_patched_runtime_20260928/venv/bin/python"
    env, resolved = execution_environment(baseline)
    env = {k: v for k, v in env.items() if not k.startswith(("PYTHON", "LD_"))}
    env.update(PYTHONPATH=baseline["core_root"], PYTHONDONTWRITEBYTECODE="1", PYTHONHASHSEED="0")
    output.mkdir(parents=True)
    result = dict(status="running", source=record(Path(__file__)), candidate=record(candidate_path),
                  plan=record(plan_path), interpreter=record(python), resolved=resolved, runs=[],
                  scientific_timings_admitted=False, genome_scale_equivalence=False)
    watched = [result["candidate"], result["plan"], result["interpreter"]]
    for ref in baseline["core_sources"]:
        actual = record(Path(ref["absolute_path"]))
        if actual["sha256"] != ref["sha256"] or actual["bytes"] != ref["bytes"]:
            raise ValueError("Frozen core changed")
        watched.append(actual)
    write(output / "started.json", result)
    try:
        for index, label in enumerate(("high", "phylo")):
            prior_path = repo / f"benchmarks/work/threadripper_reporting_{label}_20260928/prior_output_parity.json"
            prior = json.loads(prior_path.read_text())
            for ref in prior["current"]["evidence"]:
                check(ref)
                watched.append(ref)
            watched.append(record(prior_path))
            row = fixture_run(plan["runs"][index], source, output / f"run_{index:02d}", output / "inputs")
            row["native_argv"][0] = str(python)
            target = Path(row["prepared_input_directory"])
            target.mkdir(parents=True)
            for ref in row["dataset"]["inputs"]:
                shutil.copyfile(ref["path"], target / Path(ref["path"]).name)
            root = output / f"run_{index:02d}"
            root.mkdir()
            env.update(NUMBA_CACHE_DIR=str(root / "numba_cache"), HOME=str(root / "home"))
            Path(env["NUMBA_CACHE_DIR"]).mkdir()
            Path(env["HOME"]).mkdir()
            write(root / "started.json", dict(run=row, environment=env, prior=record(prior_path)))
            with (root / "native.log").open("x") as log:
                process = subprocess.Popen(row["native_argv"], cwd=row["cwd"], env=env,
                                           stdout=log, stderr=subprocess.STDOUT, start_new_session=True)
                try:
                    code = process.wait(timeout=900)
                except subprocess.TimeoutExpired:
                    os.killpg(process.pid, signal.SIGKILL)
                    process.wait()
                    raise
            measurement = dict(command=row["native_argv"], cwd=row["cwd"], exit_code=code, timed_out=False)
            write(root / "finished.json", measurement)
            native = validate(row, measurement)
            equal = parity(prior, Path(row["configuration"]["output"]))
            result["runs"].append(dict(method=row["native_method"], native=native, parity=equal,
                                       log=record(root / "native.log")))
        for ref in watched:
            check(ref)
        result.update(status="both_native_prediction_fixtures_match", checked_records=watched)
    except Exception as error:
        result.update(status="failed", error_type=type(error).__name__, error=str(error))
        raise
    finally:
        write(output / "result.json", result)
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    print(run(args.repo.resolve(), args.output.resolve())["status"])
