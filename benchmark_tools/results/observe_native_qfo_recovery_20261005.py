"""Bounded collector readback after an agent command-launcher outage."""

from collections import Counter
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import time


ROOT = Path(__file__).resolve().parents[2]
RESULTS = ROOT / "benchmark_tools/results"
OUTPUT = RESULTS / "native_qfo_agent_observation_recovery_20261005.json"
LO, HI = 1791227023.2308886, 1791227816.754092


def record(path):
    path = Path(path).resolve()
    with path.open("rb") as stream:
        digest = hashlib.file_digest(stream, "sha256").hexdigest()
    return {"path": str(path), "bytes": path.stat().st_size, "sha256": digest}


def run(argv):
    result = subprocess.run(argv, text=True, capture_output=True, check=True)
    return {"command": argv, "stdout": result.stdout, "stderr": result.stderr,
            "returncode": result.returncode}


assert not OUTPUT.exists()
plan_ref = record(RESULTS / "native_factorial_receipt_amendment_20261004/plan.json")
assert plan_ref["sha256"] == "6c87babcbb5581830e0b9e7b9bf9aaba30a85bde1c4ab465e561017e67e9c89b"
plan = json.loads(Path(plan_ref["path"]).read_text())
assert len(plan["helper_sources"]) == 920
for ref in plan["helper_sources"]:
    assert record(ref["path"]) == ref
request = record(ROOT / "benchmarks/work/native_factorial_launch_20261004/request_07_receipt_amended.json")
assert request["sha256"] == "1078eb4bb39bd2f068e3faae3df1ca0131c2cf848c333aa0bc18424fe82213c7"
attachment = Path("/home/bizon/.codex/attachments/c53784d4-1029-4fd4-9dfb-dea6a2632e17/pasted-text-1.txt")
assert attachment.read_bytes() == (ROOT / "benchmark_tools/PUBLICATION_GOAL_20261003.txt").read_bytes()
paths = sorted((ROOT / "benchmarks/results/native_factorial_cost_v2_20261004/run_07/measurement").glob("point_*.json"))[:-2]
first = max(i for i, path in enumerate(paths) if path.stat().st_mtime <= LO)
last = next(i for i, path in enumerate(paths) if path.stat().st_mtime >= HI)
selected = paths[first:last + 1]
indices = [int(path.stem.split("_")[1]) for path in selected]
assert indices == list(range(indices[0], indices[-1] + 1))
refs, samples, parent_counts = [], [], Counter()
for path in selected:
    ref = record(path)
    point = json.loads(path.read_text())
    assert record(path) == ref and point["schema"] == "native_lineage_v1"
    assert point["native_membership"].startswith("0::/system.slice/slurmstepd.scope/job_22437/step_0/")
    parent = [row for row in point["parent"] if row["scope"] == "/system.slice/slurmstepd.scope/job_22437"]
    child = [row for row in point["children"] if row["scope"] == "/system.slice/slurmstepd.scope/job_22437/step_0"]
    assert len(parent) in (1, 2) and len(child) == 1
    parent_counts[len(parent)] += 1
    usage = int(dict(line.split() for line in parent[0]["raw"].splitlines())["usage_usec"])
    step_usage = int(dict(line.split() for line in child[0]["raw"].splitlines())["usage_usec"])
    affinity = point["thread_affinity"]
    assert affinity["status"] == "observed_within_affinity" and not affinity["errors"]
    refs.append(ref)
    samples.append({"index": indices[len(samples)], "file_mtime_unix": path.stat().st_mtime,
                    "started_monotonic_ns": parent[0]["started_ns"], "job_usage_usec": usage,
                    "native_step_usage_usec": step_usage})
cadence = [(b["started_monotonic_ns"] - a["started_monotonic_ns"]) / 1e9
           for a, b in zip(samples, samples[1:])]
assert cadence and all(t > 0 for t in cadence)
assert all(b["job_usage_usec"] >= a["job_usage_usec"] and b["native_step_usage_usec"] >= a["native_step_usage_usec"]
           for a, b in zip(samples, samples[1:]))
accounting = run(["sacct", "-j", "22437,22440", "--noheader", "--parsable2",
                  "--format=JobIDRaw,State,ExitCode,Elapsed,AllocCPUS,NodeList"])
rows = {row.split("|")[0]: row.split("|") for row in accounting["stdout"].splitlines()}
assert rows["22437"][1] == "RUNNING" and rows["22440"][1] == "PENDING"
controller = run(["scontrol", "show", "job", "22437", "-o"])
assert "Restarts=0" in controller["stdout"] and "Requeue=0" in controller["stdout"]
assert "Comment=" + request["sha256"] in controller["stdout"]
reviewer = run(["scontrol", "show", "job", "22440", "-o"])
assert "Dependency=afterany:22437" in reviewer["stdout"] and "Reason=Dependency" in reviewer["stdout"]
python = str(ROOT / "benchmarks/work/native_factorial_review_py310_20261004/bin/python")
probe = run([python, "-B", "-c", "import json,psutil; p=psutil.Process(265470); print(json.dumps(dict(pid=p.pid,created=p.create_time(),affinity=p.cpu_affinity(),argv=p.cmdline(),workers=[dict(pid=c.pid,created=c.create_time(),cpu=c.cpu_times().user+c.cpu_times().system,status=c.status()) for c in p.children(recursive=True)])))"])
live = json.loads(probe["stdout"])
assert live["created"] == 1791220584.63 and live["affinity"] == list(range(32)) and live["argv"][-1] == "7"
assert [w["pid"] for w in live["workers"]] == list(range(265473, 265481))
result = {"schema": "native_qfo_agent_observation_recovery_v1", "observed_unix_ns": time.time_ns(),
          "source": record(__file__), "observer_command": sys.orig_argv, "plan": plan_ref, "request": request,
          "goal_source": record(attachment), "helper_sources_checked": 920,
          "previous_agent_failure": {"error": "bwrap: loopback: Failed RTM_NEWADDR: Operation not permitted",
                                     "consecutive_goal_turns": 3, "source": "retained conversation tool outputs",
                                     "current_native_job_state_was_unverified": True},
          "earlier_readback_failure": {"exception": "AssertionError", "assumption": "exactly one parent counter row",
                                       "observed": "some points have two bracketing parent rows", "output_written": False},
          "observation_bounds": {"last_pre_gap_process_observation_unix": LO,
                                 "first_post_recovery_process_observation_unix": HI,
                                 "exact_launcher_outage_start_end_known": False},
          "collector_window": {"first": samples[0], "last": samples[-1], "points": len(samples),
                               "contiguous_indices": True, "minimum_cadence_seconds": min(cadence),
                               "maximum_cadence_seconds": max(cadence), "parent_rows_per_point": dict(parent_counts),
                               "cpu_counters_nondecreasing": True, "observed_affinity_errors": 0, "records": refs},
          "accounting": accounting, "native_controller": controller, "reviewer_controller": reviewer,
          "original_process": probe, "new_job_submitted": False, "native_job_restarted": False,
          "full_resource_review": False, "scientific_accuracy_admitted": False, "publication_ready": False,
          "limitations": ["Bounded collector continuity/counter/affinity check, not full raw-resource replay or terminal admission.",
                          "One/two parent rows are preserved as recorded; the first parent reading supplies the cadence/counter observation.",
                          "Current process identities and zero scheduler restarts do not prove every aspect of uninterrupted execution.",
                          "External agent environment change restored observation; no sandbox or host service was repaired by this observer.",
                          "Shared-host timing effects remain unknown and potentially tool-dependent; no corrected time or isolated ranking.",
                          "Original results, failure history, frozen sources and successor gates remain unchanged."]}
with OUTPUT.open("x") as stream:
    json.dump(result, stream, indent=2, sort_keys=True)
    stream.write("\n")
print(json.dumps({"output": record(OUTPUT), "points": len(samples), "cadence_min": min(cadence),
                  "cadence_max": max(cadence), "parent_rows_per_point": dict(parent_counts)}))
