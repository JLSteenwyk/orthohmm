"""Execute one fixed mode-control panel row with all failures retained."""

import argparse
import json
import os
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.prepare_ob_candidate_neighborhood import record, check
from benchmark_tools.run_simulation_methods import read_frozen
from benchmark_tools.run_simulation_tree_mode_control import run, METHOD_SHA, RESULT_SHA
from benchmark_tools.prepare_simulation_mode_panel import panel_rows


def execute(root, manifest_path, digest, index):
    panel = read_frozen(manifest_path, digest)
    if panel["status"] != "mode_panel_frozen_not_executed" or panel["accuracy_evaluated"] is not False:
        raise ValueError("Wrong mode panel state")
    for key in ("source", "method_manifest", "baseline_results", "pilot_admission"):
        check(panel[key])
    rows = panel_rows(read_frozen(Path(panel["method_manifest"]["path"]), METHOD_SHA),
                      read_frozen(Path(panel["baseline_results"]["path"]), RESULT_SHA))
    if rows != panel["rows"] or isinstance(index, bool) or not 0 <= index < len(rows):
        raise ValueError("Panel inventory or array index differs")
    row = rows[index]
    parent = root / "benchmarks/results/simulation_mode_panel_v1"
    status_dir = parent / "tasks" / row["label"]
    status_dir.mkdir(parents=True, exist_ok=False)
    status = {"status": "running", "row": row, "panel": record(manifest_path), "source": record(__file__),
              "job_id": os.environ.get("SLURM_JOB_ID"), "array_job_id": os.environ.get("SLURM_ARRAY_JOB_ID"),
              "array_task_id": os.environ.get("SLURM_ARRAY_TASK_ID"), "accuracy_evaluated": False}
    try:
        if row["reuse_pilot"]:
            admission = json.loads(Path(panel["pilot_admission"]["path"]).read_text())
            check(admission["pilot_report"])
            status.update(status="reused_admitted_pilot", admission=panel["pilot_admission"])
        elif not row["methods"]:
            status["status"] = "no_available_inferred_baseline"
        else:
            result = run(root, row["label"], parent / row["label"], tuple(row["methods"]))
            status.update(status="executed_pending_independent_admission",
                          native_result=record(parent / row["label"] / "results.json"), native_status=result["status"])
        read_frozen(manifest_path, digest)
    except Exception as error:
        status.update(status="failed", error=str(error), error_type=type(error).__name__)
        raise
    finally:
        (status_dir / "status.json").write_text(json.dumps(status, indent=2, sort_keys=True) + "\n")
    return status


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", required=True, type=Path)
    parser.add_argument("--manifest", required=True, type=Path)
    parser.add_argument("--manifest-sha256", required=True)
    parser.add_argument("--index", required=True, type=int)
    args = parser.parse_args()
    execute(args.root.resolve(), args.manifest.resolve(), args.manifest_sha256, args.index)
