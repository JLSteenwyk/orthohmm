"""Score authenticated completed cap400 inference without rerunning science."""

import argparse
import hashlib
import importlib
import json
import os
from pathlib import Path
import subprocess
import sys
import time

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.audit_cap400_scoring_parent import audit
from benchmark_tools.prepare_ob_candidate_neighborhood import check, record
from benchmark_tools.recover_full_qfo_search import short_id
from benchmark_tools.report_hillclimb import endpoint_vector, readback_scores
from benchmark_tools.run_full_search_hillclimb import runtime
from benchmark_tools.run_membership_hillclimb import score_qfo, write_json
from benchmark_tools.run_search_recovery import host_activity

PROTOCOL = "hillclimb_cap400_scoring_only_recovery_20261016.json"
ADMISSION = "benchmark_tools/results/cap400_scoring_parent_admission_20261016.json"
READY = "completed_original_cap400_native_reverified_scoring_only_pending"
SUCCESS = "original_cap400_all_six_scores_recovered_pending_independent_admission"
SCORING_HELPERS = ("run_membership_hillclimb.py", "prepare_qfo_factorial_pairs.py",
    "qfo_filter_pairs.py", "qfo_summarize_scores.py", "run_qfo_recovered_assessment.py",
    "validate_qfo_native_assessment.py", "report_hillclimb.py")


def require_admission(admission, evidence):
    if (admission.get("status") != READY or evidence.get("status") != READY
            or admission.get("scores") != {} or evidence.get("scores") != {}
            or admission.get("native_inference_rerun") is not False
            or evidence.get("native_inference_rerun") is not False
            or admission.get("producer") != evidence.get("producer")
            or admission.get("scheduler") != evidence.get("scheduler")
            or admission.get("scientific_revision") != evidence.get("scientific_revision")
            or admission.get("coverage") != {k: evidence["native"][k]
                for k in ("genes", "species", "groups", "pairs")}
            or (evidence["native"]["genes"], evidence["native"]["species"]) != (984137, 78)):
        raise ValueError("Require exact unscored admission of completed original cap400 inference")


def read_admission(pin):
    check(pin)
    admission = json.loads(Path(pin["path"]).read_text())
    for key in ("audit", "source", "producer"):
        check(admission[key])
    evidence = json.loads(Path(admission["audit"]["path"]).read_text())
    require_admission(admission, evidence)
    for item in evidence["checked"]:
        check(item)
    return admission, evidence


def verify_helpers(source, parent, revision):
    actual = subprocess.check_output(["git", "-C", str(source), "rev-parse", "HEAD"], text=True).strip()
    if actual != revision:
        raise ValueError("Frozen recovery executor revision changed")
    pins = [record(p) for p in sorted((source / "benchmark_tools").glob("*.py"))]
    pins += [record(p) for p in sorted((source / "orthohmm").rglob("*.py"))]
    pins += [record(source / "qfo_benchmark/og_to_pairwise.py"),
             record(source / "benchmark_tools" / PROTOCOL),
             record(source / "benchmark_tools/recover_cap400_scores.slurm")]
    for pin in pins:
        relative = Path(pin["path"]).relative_to(source)
        blob = subprocess.check_output(["git", "-C", str(source), "show", f"{revision}:{relative}"])
        if hashlib.sha256(blob).hexdigest() != pin["sha256"]:
            raise ValueError("Recovery executor differs from committed content")
    for name in SCORING_HELPERS:
        original = Path(parent["source_root"]) / "benchmark_tools" / name
        if record(original)["sha256"] != record(source / "benchmark_tools" / name)["sha256"]:
            raise ValueError("Recovery changed an original scoring helper")
    converter = source / "qfo_benchmark/og_to_pairwise.py"
    blob = subprocess.check_output(["git", "-C", parent["scientific_source"], "show",
        parent["scientific_revision"] + ":qfo_benchmark/og_to_pairwise.py"])
    if hashlib.sha256(blob).hexdigest() != record(converter)["sha256"]:
        raise ValueError("Converter differs from original scientific revision")
    return pins


def prepare(repo, source, output):
    if source == repo or output.exists():
        raise ValueError("Require new output and a separate frozen executor")
    admission_pin = record(repo / ADMISSION)
    if record(source / ADMISSION)["sha256"] != admission_pin["sha256"]:
        raise ValueError("Frozen admission differs from retained evidence")
    admission, evidence = read_admission(admission_pin)
    parent = json.loads(Path(evidence["execution_plan"]["path"]).read_text())
    if Path(parent["repo"]).resolve() != repo:
        raise ValueError("Parent belongs to another repository")
    revision = subprocess.check_output(["git", "-C", str(source), "rev-parse", "HEAD"], text=True).strip()
    pins = verify_helpers(source, parent, revision)
    identifier = short_id(repo, output, 1)
    work = repo / "qfo_benchmark/w" / ("hc1_" + identifier)
    if work.exists() or runtime() != parent["runtime"]:
        raise ValueError("Scoring work already exists or original runtime changed")
    result = {"status": "prepared_cap400_scoring_only_recovery", "repo": str(repo),
        "source_root": str(source), "source_revision": revision, "output": str(output),
        "admission": admission_pin, "parent_root": parent["output"],
        "parent_plan": evidence["execution_plan"], "parent_report": evidence["producer"],
        "sources": pins, "runtime": parent["runtime"], "scoring_id": identifier,
        "scientific_revision": parent["scientific_revision"], "native_inference_rerun": False,
        "scope": "Scoring-only recovery; original failed fresh and downstream attempts retained; no uninterrupted pipeline timing."}
    output.mkdir(parents=True)
    write_json(output / "execution_plan.json", result)
    return result


def verify_imported_converter(source):
    converter = importlib.import_module("qfo_benchmark.og_to_pairwise")
    helper = importlib.import_module("benchmark_tools.prepare_qfo_factorial_pairs")
    expected = source / "qfo_benchmark/og_to_pairwise.py"
    if (Path(converter.__file__).resolve() != expected
            or Path(helper.__file__).resolve() != source / "benchmark_tools/prepare_qfo_factorial_pairs.py"
            or helper._strip_to_uniprot is not converter._strip_to_uniprot):
        raise ValueError("Scoring imported a different converter or helper")
    return record(expected)


def run(path):
    plan_pin = record(path)
    plan = json.loads(path.read_text())
    repo, source, output = (Path(plan[k]) for k in ("repo", "source_root", "output"))
    if (Path(__file__).resolve() != source / "benchmark_tools/recover_cap400_scores.py"
            or path.resolve() != output / "execution_plan.json" or (output / "results.json").exists()):
        raise ValueError("Wrong frozen executor or attempted recovery namespace")
    if (os.environ.get("SLURM_CPUS_PER_TASK"), os.environ.get("SLURM_MEM_PER_NODE"),
            os.environ.get("SLURMD_NODENAME"), os.environ.get("CUDA_VISIBLE_DEVICES")) != (
            "8", "65536", "bizon", ""):
        raise ValueError("Require CPU-only 8-CPU 64-GiB bizon allocation")
    result = {"status": "validating_original_cap400_native_for_scoring", "execution_plan": plan_pin,
        "job_id": os.environ.get("SLURM_JOB_ID"), "source_revision": plan["source_revision"],
        "scientific_revision": plan["scientific_revision"], "scores": {}, "native_inference_rerun": False,
        "defaults_promoted": False, "goal_completed": False, "host_before": host_activity()}
    write_json(output / "results.json", result)
    checked = [plan_pin, plan["admission"], plan["parent_plan"], plan["parent_report"],
               *plan["sources"], *plan["runtime"]["files"]]
    started = time.monotonic()
    try:
        for pin in checked:
            check(pin)
        admission, prior = read_admission(plan["admission"])
        parent = json.loads(Path(plan["parent_plan"]["path"]).read_text())
        if (prior["execution_plan"] != plan["parent_plan"] or prior["producer"] != plan["parent_report"]
                or parent["output"] != plan["parent_root"] or parent["repo"] != str(repo)
                or plan["scientific_revision"] != parent["scientific_revision"]
                or plan.get("native_inference_rerun") is not False
                or short_id(repo, output, 1) != plan["scoring_id"]
                or runtime() != plan["runtime"] or plan["runtime"] != parent["runtime"]):
            raise ValueError("Recovery changed original identity, runtime or namespace")
        if verify_helpers(source, parent, plan["source_revision"]) != plan["sources"]:
            raise ValueError("Recovery executor inventory differs from preparation")
        result["converter"] = verify_imported_converter(source)
        fresh = audit(repo, Path(plan["parent_root"]))
        require_admission(admission, fresh)
        if fresh["native"] != prior["native"]:
            raise ValueError("Native content differs from independent retained audit")
        checked += [admission["audit"], *fresh["checked"]]
        native = fresh["native"]
        result.update(status="scoring_original_cap400_native", native_audit=admission["audit"],
            native_products=native["products"], native_summary=native["summary"],
            coverage={k: native[k] for k in ("genes", "species", "groups", "pairs")},
            native_verification_wall_s=time.monotonic() - started)
        write_json(output / "results.json", result)
        scoring_started = time.monotonic()
        scored = score_qfo(repo, Path(native["phylogeny"]), output, plan["scoring_id"])
        scores, details, pins = readback_scores(repo, "QfO", scored, native, output,
                                               "ohmm_hillclimb_r1_" + plan["scoring_id"])
        endpoint_vector("QfO", scores)
        checked += pins
        result.update(status=SUCCESS, scores=scores, details=details, scoring=scored,
            deltas={k: value - parent["control"]["scores"][k] for k, value in scores.items()},
            scoring_and_readback_wall_s=time.monotonic() - scoring_started)
    except Exception as error:
        result.update(status="failed_cap400_scoring_only_recovery", error_type=type(error).__name__, error=str(error))
    finally:
        try:
            for pin in checked:
                check(pin)
            if runtime() != plan["runtime"]:
                raise ValueError("Original runtime changed during scoring recovery")
            result["postflight_pins_verified"] = True
        except Exception as error:
            result.update(status="failed_cap400_scoring_only_recovery", postflight_error=str(error))
        result.update(total_recovery_wall_s=time.monotonic() - started, host_after=host_activity())
        write_json(output / "results.json", result)
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="operation", required=True)
    preparation = sub.add_parser("prepare")
    for name in ("repo", "source", "output"):
        preparation.add_argument("--" + name, required=True, type=Path)
    execution = sub.add_parser("run")
    execution.add_argument("--execution-plan", required=True, type=Path)
    args = parser.parse_args()
    result = (prepare(args.repo.resolve(), args.source.resolve(), args.output.resolve())
              if args.operation == "prepare" else run(args.execution_plan.resolve()))
    print(result["status"])
    raise SystemExit(0 if result["status"] in ("prepared_cap400_scoring_only_recovery", SUCCESS) else 1)
