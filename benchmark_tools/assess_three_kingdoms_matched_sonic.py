"""Validate and score the frozen matched-input SonicParanoid run after completion."""

import argparse
import json
from pathlib import Path
import subprocess
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.prepare_ob_candidate_neighborhood import record, check
from benchmark_tools.run_simulation_methods import read_frozen
from benchmark_tools.run_three_kingdoms_matched_sonic import verify
from benchmark_tools.verify_ygob_validation import require_completed_job
from benchmark_tools.validate_sonicparanoid_groups import validate
from benchmark_tools.audit_three_kingdoms_pair_counts import count, compare
from benchmark_tools.summarize_three_kingdoms_parity import parse_score

PLAN_SHA = "5c3ce4608b54e8e550f0ea855754879211b842f0175b2e8e864cf4b5c2998283"
EXECUTOR = "7c3d0784d1a22da0a0e37db6860977c6649fd2ca"
JOB = 21795


def validate_execution(plan, execution, scheduler, plan_record, runner_record):
    expected_scheduler = {"JobIDRaw": str(JOB), "State": "COMPLETED", "ExitCode": "0:0",
                          "NodeList": "bizon", "AllocCPUS": "32", "ReqMem": "192G"}
    if any(scheduler.get(k) != v for k, v in expected_scheduler.items()):
        raise ValueError("Wrong terminal scheduler state or resources")
    if (execution.get("status") != "process_succeeded_pending_native_admission"
            or execution.get("exit_code") != 0 or execution.get("accuracy_admitted") is not False
            or execution.get("job_id") != str(JOB) or execution.get("node") != "bizon"):
        raise ValueError("Native execution not successful or wrong identity")
    if execution["plan"] != plan_record or execution["source"] != runner_record:
        raise ValueError("Execution source/plan changed")
    if execution["native_argv"] != plan["native_argv"]:
        raise ValueError("Native command changed")
    if not execution["finished_epoch"] >= execution["started_epoch"] > 0:
        raise ValueError("Invalid native timestamps")
    if not execution["runtime_before"] or execution["runtime_before"] != execution["runtime_after"]:
        raise ValueError("Native runtime changed")
    root = Path(plan["output_root"])
    expected = [{**item, "path": str(root / "input" / Path(item["path"]).name)} for item in plan["inputs"]]
    if execution["copied_inputs"] != expected:
        raise ValueError("Input copies differ from frozen plan")
    for key, name in (("native_log", "native.log"), ("timing", "time.txt")):
        if execution[key]["path"] != str(root / name) or execution[key]["bytes"] <= 0:
            raise ValueError("Invalid native log or timing")


def assess(repo, destination):
    if destination.exists():
        raise FileExistsError(destination)
    accounting = subprocess.check_output(["sacct", "-j", str(JOB), "--parsable2",
        "--format=JobIDRaw,State,ExitCode,Elapsed,NodeList,AllocCPUS,ReqMem"], text=True)
    scheduler = require_completed_job(accounting, JOB)
    plan_path = repo / "benchmark_tools/results/three_kingdoms_sonic_matched_commands_20260918.json"
    plan = read_frozen(plan_path, PLAN_SHA)
    root = Path(plan["output_root"])
    executor = repo / "benchmarks/work/publication_three_kingdoms_sonic_matched_v1"
    if subprocess.check_output(["git", "-C", str(executor), "rev-parse", "HEAD"], text=True).strip() != EXECUTOR:
        raise ValueError("Inference executor revision changed")
    subprocess.run(["git", "-C", str(executor), "diff", "--exit-code", "HEAD", "--", "benchmark_tools"], check=True)
    status = record(root / "execution.json")
    execution = json.loads(Path(status["path"]).read_text())
    validate_execution(plan, execution, scheduler, record(plan_path),
                       record(executor / "benchmark_tools/run_three_kingdoms_matched_sonic.py"))
    _, _, env, runtime = verify(plan_path, PLAN_SHA)
    if execution["environment"] != env or execution["runtime_after"] != runtime:
        raise ValueError("Effective runtime/environment changed")
    outputs = [record(p) for p in sorted((root / "output").rglob("*")) if p.is_file()]
    if not outputs or outputs != execution["outputs"]:
        raise ValueError("Native output inventory changed")
    if {str(p) for p in (root / "input").iterdir()} != {r["path"] for r in execution["copied_inputs"]}:
        raise ValueError("Input directory membership changed")
    tables = list((root / "output").rglob("ortholog_groups.tsv"))
    if len(tables) != 1:
        raise ValueError("Require exactly one native ortholog_groups.tsv")
    records = [status, record(plan_path), execution["source"], execution["native_log"], execution["timing"],
               *plan["checked_records"], *execution["copied_inputs"], *outputs]
    for item in records:
        check(item)
    destination.mkdir(parents=True, exist_ok=False)
    result = {"status": "assessment_started", "accuracy_admitted": False, "source": record(__file__),
              "scheduler": scheduler, "accounting": accounting, "checked_records": records}
    try:
        normalized, score_path = destination / "orthogroups.txt", destination / "score.txt"
        python = json.loads(Path(plan["runtime"]["path"]).read_text())["python"]["path"]
        conversion = [python, "-B", "-s", str(repo / "benchmark_tools/normalize_three_kingdoms_orthogroups.py"),
                      "sonicparanoid", str(tables[0]), str(normalized)]
        result["conversion_argv"] = conversion
        with (destination / "conversion.log").open("x") as stream:
            subprocess.run(conversion, env=env, cwd=destination, stdout=stream, stderr=subprocess.STDOUT, check=True)
        result["native_validation"] = validate(tables[0], [Path(r["path"]) for r in execution["copied_inputs"]],
                                               root / "output/snapshot.tsv", normalized)
        reference = repo / "three_kingdoms/busco/reference_orthogroups.txt"
        scoring = [python, "-B", "-s", str(repo / "three_kingdoms/score_against_busco.py"),
                   "--predictions", str(normalized), "--reference", str(reference),
                   "--label", "SonicParanoid matched-input Three Kingdoms"]
        result["scoring_argv"] = scoring
        with score_path.open("x") as stream:
            subprocess.run(scoring, env=env, cwd=destination, stdout=stream, stderr=subprocess.STDOUT, check=True)
        independent = count(reference, normalized)
        compare(independent, parse_score(score_path))
        if (independent["reference_orthogroups"], independent["reference_genes"],
                independent["true_positive_gene_pairs"] + independent["false_negative_gene_pairs"]) != (255, 2035, 7352):
            raise ValueError("Reference universe changed")
        for item in records + result["native_validation"]["evidence"]:
            check(item)
        if {str(p) for p in (root / "output").rglob("*") if p.is_file()} != {r["path"] for r in outputs}:
            raise ValueError("Native output membership changed during assessment")
        result.update(status="matched_three_kingdoms_sonic_score_verified", accuracy_admitted=True,
                      counts=independent, normalized=record(normalized), score=record(score_path),
                      conversion_log=record(destination / "conversion.log"),
                      limitations=["Contemporary matched-input run, not an isolated stop-marker causal test.",
                                   "BUSCO reference-gene co-membership only; no penalty outside reference.",
                                   "Shared-host timing is descriptive, not dedicated scaling."])
    except Exception as exc:
        result.update(status="assessment_failed", error=str(exc), accuracy_admitted=False)
        raise
    finally:
        (destination / "report.json").write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    assess(args.repo.resolve(), args.output.resolve())
