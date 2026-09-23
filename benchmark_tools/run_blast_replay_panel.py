"""Execute the frozen diagnostic panel once; never admit interrupted results."""

import argparse
import json
import os
from pathlib import Path
import subprocess
import sys
import time

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.prepare_ob_candidate_neighborhood import record, check
from benchmark_tools.run_simulation_methods import read_frozen
from benchmark_tools.run_qfo_corrected_blast import verify, environment

PANEL_SHA = "79a809b43d407d4d5ed705e6840873cede64d7b688eee571869b1607f6688872"


def query_file_records(queries):
    result = []
    for query in queries:
        if (set(query) != {"id", "input_ordinal_0based", "path", "bytes", "sha256"}
                or not isinstance(query["id"], str) or not query["id"]
                or type(query["input_ordinal_0based"]) is not int
                or query["input_ordinal_0based"] < 0):
            raise ValueError("Unexpected diagnostic query metadata")
        result.append({key: query[key] for key in ("path", "bytes", "sha256")})
    return result


def validate_commands(panel, original):
    queries = [panel["combined"], *panel["queries"]]
    names = ["combined", *[f"single_{i:02d}" for i in range(5)]]
    if len(panel["queries"]) != 5 or len(panel["commands"]) != 6:
        raise ValueError("Require exactly six frozen diagnostic commands")
    root = Path(panel["combined"]["path"]).parent
    for name, query, command in zip(names, queries, panel["commands"]):
        expected = list(original)
        expected[expected.index("-i") + 1] = query["path"]
        expected[expected.index("-o") + 1] = str(root / (name + ".blast"))
        if command != {"name": name, "argv": expected}:
            raise ValueError("Changed replay command")
    return root


def run(panel_path):
    if os.environ.get("SLURM_CPUS_PER_TASK") != "180" or not os.environ.get("SLURM_JOB_ID"):
        raise ValueError("Require scheduled 180-CPU allocation")
    panel = read_frozen(panel_path, PANEL_SHA)
    plan_path, runtime_path = (Path(panel[key]["path"]) for key in ("plan", "runtime"))
    plan = verify(plan_path, runtime_path)
    root = validate_commands(panel, plan["search_commands"]["blast"])
    checked = [panel["source"], panel["plan"], panel["runtime"], panel["input"],
               panel["combined"], *query_file_records(panel["queries"]), *panel["database"]]
    for item in checked:
        check(item)
    execution = root / "execution"
    if execution.exists() or any((root / (c["name"] + ".blast")).exists() for c in panel["commands"]):
        raise FileExistsError("No implicit restart or overwrite of diagnostic runs")
    execution.mkdir(exist_ok=False)
    report = {"status": "running", "panel": record(panel_path), "source": record(__file__),
              "job_id": os.environ["SLURM_JOB_ID"], "node": os.uname().nodename,
              "commands": panel["commands"], "environment": environment(), "stages": [],
              "search_admitted": False, "reuse_authorized": False}

    def save():
        temporary = execution / "status.tmp"
        with temporary.open("w") as stream:
            json.dump(report, stream, indent=2, sort_keys=True)
            stream.write("\n")
            stream.flush()
            os.fsync(stream.fileno())
        temporary.replace(execution / "status.json")

    try:
        for command in panel["commands"]:
            stage = {"name": command["name"], "started_epoch": time.time()}
            report["stages"].append(stage)
            save()
            log = execution / (command["name"] + ".log")
            with log.open("xb") as stream:
                result = subprocess.run(command["argv"], cwd=plan["output_root"] + "/work",
                                        env=environment(), stdout=stream, stderr=subprocess.STDOUT)
            stage.update(exit_code=result.returncode, finished_epoch=time.time(), log=record(log))
            output = root / (command["name"] + ".blast")
            if output.exists():
                stage["output"] = record(output)
            save()
        for item in checked:
            check(item)
        verify(plan_path, runtime_path)
        report["status"] = "native_replays_finished_pending_comparison"
        if any(stage["exit_code"] != 0 for stage in report["stages"]):
            raise RuntimeError("At least one native replay exited nonzero; inspect retained diagnostics")
    except BaseException as exc:
        report.update(status="failed_or_interrupted", error=repr(exc))
        raise
    finally:
        save()
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--panel", type=Path, required=True)
    args = parser.parse_args()
    print(json.dumps({"status": run(args.panel.resolve())["status"]}))
