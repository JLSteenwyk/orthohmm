"""Independently validate completed corrected-QfO candidate parameter preparation."""

import argparse
import json
import math
from pathlib import Path
import subprocess
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.prepare_ob_candidate_neighborhood import ARMS, check, record
from benchmark_tools.prepare_qfo_candidate_neighborhood import ENVIRONMENT, PLAN_SHA, PROTOCOL_SHA
from benchmark_tools.run_simulation_methods import read_frozen
from benchmark_tools.verify_ygob_validation import require_completed_job

EXECUTOR_COMMIT = "6405da57618cc703fc1ca0f51e450d3cf6e2a05f"
SOURCE_SHA = "e6e7d889fe0e0cc545ac9c01c84ae146ca80828324bb49e4eb9023c4e3263aeb"


def validate_report(report, baseline, scheduler):
    if (scheduler["State"], scheduler["ExitCode"], scheduler["NodeList"], scheduler["AllocCPUS"]) != (
            "COMPLETED", "0:0", "bizon", "2"):
        raise ValueError("Require successful terminal 2-CPU preparation")
    if (report["status"] != "corrected_qfo_five_candidate_neighborhood_arms_prepared_unscored"
            or report["job_id"] != scheduler["JobIDRaw"] or report["accuracy_evaluated"] is not False
            or report["publication_ready"] is not False or report["environment"] != ENVIRONMENT):
        raise ValueError("Wrong preparation status, job or environment")
    if [row["label"] for row in report["arms"]] != [label for label, _ in ARMS]:
        raise ValueError("Candidate panel is incomplete or reordered")
    parameters = baseline["expansion"]["parameters"]
    for row, (_, delta) in zip(report["arms"], ARMS):
        if (row["status"] != "candidate_prepared_unscored" or row["delta"] != delta
                or row["applied_parameters"] != {**parameters, **delta}
                or type(row["engine_calls"]) is not int or row["engine_calls"] != 1
                or row["engine_fixed_profile_report"]["parameters"] != parameters
                or row["candidate_arm"]["expansion"] != row["engine_fixed_profile_report"]):
            raise ValueError("Incorrect candidate parameters or execution state")
        if (type(row["incremental_seconds"]) not in (int, float)
                or not math.isfinite(row["incremental_seconds"]) or row["incremental_seconds"] < 0):
            raise ValueError("Invalid incremental time")
    control = report["arms"][0]
    if control.get("baseline_byte_equivalent") is not True:
        raise ValueError("Control was not admitted as byte-equivalent")
    for key in ("candidate_partition", "membership_constraints"):
        if any(control["candidate_arm"][key][k] != baseline[key][k] for k in ("sha256", "bytes")):
            raise ValueError("Control differs from pinned baseline")


def admit(root, job, manifest_sha, output):
    if output.exists():
        raise FileExistsError(output)
    accounting = subprocess.check_output(["sacct", "-j", str(job), "--parsable2",
        "--format=JobIDRaw,State,ExitCode,Elapsed,NodeList,AllocCPUS"], text=True)
    scheduler = require_completed_job(accounting, job)
    executor = root / "benchmarks/work/publication_qfo_parameter_candidates_v1"
    if subprocess.check_output(["git", "-C", str(executor), "rev-parse", "HEAD"], text=True).strip() != EXECUTOR_COMMIT:
        raise ValueError("Preparation executor revision changed")
    subprocess.run(["git", "-C", str(executor), "diff", "--exit-code", "HEAD", "--", "benchmark_tools", "orthohmm"], check=True)
    source = record(executor / "benchmark_tools/prepare_qfo_candidate_neighborhood.py")
    if source["sha256"] != SOURCE_SHA:
        raise ValueError("Wrong preparation source")
    results = root / "benchmark_tools/results"
    plan_path = results / "qfo_parameter_neighborhood_plan_20260919.json"
    protocol_path = results / "QFO_PARAMETER_NEIGHBORHOOD_PROTOCOL_20260919.md"
    plan = read_frozen(plan_path, PLAN_SHA)
    if record(protocol_path)["sha256"] != PROTOCOL_SHA:
        raise ValueError("Protocol changed")
    expected = [record(plan_path), record(protocol_path)]
    loaded = []
    for item in plan["inputs"]:
        path = root / item["path"]
        if record(path)["sha256"] != item["sha256"]:
            raise ValueError("Baseline evidence changed")
        expected.append(record(path))
        if path.suffix == ".json":
            loaded.append(json.loads(path.read_text()))
    _replay, baseline_manifest, native, _score = loaded
    baseline = baseline_manifest["candidate_arms"]["p1_c1"]
    runtime_path = results / "publication_native_runtime_20260916.json"
    expected.extend([native["candidate_admission"], baseline_manifest["plan"], record(runtime_path), baseline["seed_partition"]])
    directory = root / "benchmarks/results/qfo_parameter_candidates_v1"
    path = directory / "manifest.json"
    report = read_frozen(path, manifest_sha)
    validate_report(report, baseline, scheduler)
    if (report["source"] != source or report["inputs"][:len(expected)] != expected
            or report["inputs"][-1] != source or len(report["inputs"]) <= len(expected)):
        raise ValueError("Preparation provenance inventory differs")
    for item in report["inputs"][len(expected):]:
        if Path(item["path"]).parent != executor / "benchmark_tools":
            raise ValueError("Preparation helper outside frozen executor")
    records = [record(path), source, *report["inputs"]]
    for item in records:
        check(item)
    from benchmark_tools.run_qfo_corrected_factorial_cell import verify_admission
    admission = native["candidate_admission"]
    _, verified, _, _, _, _ = verify_admission(root, Path(admission["path"]), admission["sha256"], "21758", 3)
    if verified != baseline_manifest:
        raise ValueError("Wrong independently admitted baseline")
    from benchmark_tools.verify_qfo_replay_launcher import verify
    runtime = verify(Path(baseline_manifest["core_root"]), Path(baseline_manifest["launcher_root"]), runtime_path)
    if not runtime == report["runtime_before"] == report["runtime_after"]:
        raise ValueError("Preparation runtime differs")
    from benchmark_tools.audit_accuracy_checkpoint import audit as audit_numeric
    numeric = report["numeric_checkpoint"]
    replay_plan = json.loads(Path(baseline_manifest["plan"]["path"]).read_text())
    checkpoint = replay_plan["checkpoint_manifest"]
    if numeric["manifest"] != checkpoint:
        raise ValueError("Different numeric checkpoint")
    checked_numeric = audit_numeric(Path(checkpoint["path"]).parent, checkpoint["sha256"])
    if not checked_numeric["summary"] == numeric["summary"] == baseline_manifest["numeric_checkpoint"]["summary"]:
        raise ValueError("Numeric checkpoint summary differs")
    universe = set((Path(checkpoint["path"]).parent / "gene_names.txt").read_text().splitlines())
    if len(universe) != 984137:
        raise ValueError("Wrong corrected universe")
    from benchmark_tools.audit_candidate_arm import audit as audit_arm
    from benchmark_tools.trace_ob_families import partition, validate_merge_reconstruction
    helpers = [record(module.__file__) for name, module in sorted(sys.modules.items())
               if name.startswith("benchmark_tools.") and getattr(module, "__file__", None)]
    records.extend(helpers)
    seeds, _ = partition(Path(baseline["seed_partition"]["path"]), "plain", universe)
    verification = []
    for row in report["arms"]:
        arm = row["candidate_arm"]
        observed = audit_arm(arm, baseline["seed_partition"], directory / row["label"], universe, True)
        if observed != arm["content_audit"]:
            raise ValueError("Candidate content audit differs")
        groups, _ = partition(Path(arm["candidate_partition"]["path"]), "plain", universe)
        events = json.loads(Path(arm["membership_constraints"]["path"]).read_text())
        validate_merge_reconstruction(events, seeds, groups, set())
        verification.append({"label": row["label"], "genes": len(universe), "seed_groups": len(seeds),
                             "candidate_groups": len(groups), "reconstructed_merges": len(events)})
        records.extend(arm["output_files"])
    for item in records:
        check(item)
    result = {"status": "corrected_qfo_candidate_neighborhood_admitted_unscored",
              "accuracy_evaluated": False, "publication_ready": False, "source": record(__file__),
              "preparation": record(path), "scheduler": scheduler, "accounting": accounting,
              "arms": report["arms"], "verification": verification, "checked_records": records,
              "numeric_recheck": checked_numeric, "helpers": helpers,
              "limitations": ["Content/merge validation, not independent rescoring of candidate search evidence.",
                  "Four variants require native phylogeny and scoring; two CPM variants remain separate."]}
    with output.open("x") as handle:
        handle.write(json.dumps(result, indent=2, sort_keys=True) + "\n")
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("root", "output"):
        parser.add_argument("--" + name, required=True, type=Path)
    parser.add_argument("--job", required=True)
    parser.add_argument("--manifest-sha256", required=True)
    args = parser.parse_args()
    admit(args.root.resolve(), args.job, args.manifest_sha256, args.output.resolve())
