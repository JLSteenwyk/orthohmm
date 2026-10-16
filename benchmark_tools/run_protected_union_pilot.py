"""Run the precommitted complete-proteome OrthoBench protected-union pilot."""

import argparse
import csv
import hashlib
import io
import json
import math
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.fresh_protected_union_search import SUCCESS as FRESH_SUCCESS
from benchmark_tools.prepare_ob_candidate_neighborhood import check, record
from benchmark_tools.report_hillclimb import endpoint_vector, readback_scores, replay_native
from benchmark_tools.run_boundary_clear_pipeline import route_totals
from benchmark_tools.run_full_search_hillclimb import native_command, runtime, validate_native, verify_input_inventory
from benchmark_tools.run_membership_hillclimb import PRODUCTS, argument, score_orthobench, write_json
from benchmark_tools.run_search_recovery import host_activity
from benchmark_tools.run_three_kingdoms_membership_supplement import sources, timed_stage
from benchmark_tools.validate_profile_runtime import require_profile_runtime

PROTOCOL = "hillclimb_protected_union_fresh_pilot_20261016.json"
SLURM = "run_protected_union_pilot.slurm"
FRESH_GATE = "benchmark_tools/results/protected_union_fresh_readback_20261016.json"
QFO_ADMISSION = "benchmark_tools/results/hillclimb_comparison_readback_v26_20261016/shared_qfo_admission_readback.json"
CONTROL_LAUNCH = "benchmark_tools/results/canonical_downstream_submission_20261014.json"
SUCCESS = "protected_union_orthobench_native_and_score_verified_pending_independent_admission"


def committed(repo, source, revision, name):
    pin = record(repo / name)
    blob = subprocess.check_output(["git", "-C", str(source), "show", f"{revision}:{name}"])
    if hashlib.sha256(blob).hexdigest() != pin["sha256"]:
        raise ValueError("Prerequisite differs from committed admission: " + name)
    value = json.loads(Path(pin["path"]).read_text())
    check(pin)
    return value, pin


def terminal_dependencies():
    text = subprocess.check_output(["sacct", "-j", "24240,24360,24362", "-X", "-P", "-o",
        "JobIDRaw,State,ExitCode,AllocCPUS,NodeList"], text=True)
    rows = list(csv.DictReader(io.StringIO(text), delimiter="|"))
    expected = [{"JobIDRaw": job, "State": state, "ExitCode": code, "AllocCPUS": cpu, "NodeList": "bizon"}
        for job, state, code, cpu in (("24240", "FAILED", "1:0", "16"),
                                     ("24360", "COMPLETED", "0:0", "8"),
                                     ("24362", "COMPLETED", "0:0", "8"))]
    if len(rows) != len(expected) or sorted(rows, key=lambda r: r["JobIDRaw"]) != expected:
        raise ValueError("Require terminal original control, admitted QfO recovery and fresh-search gate jobs")
    return expected


def scientific_source(source, fresh):
    previous = Path(fresh["source_root"])
    original = {str(Path(p["path"]).relative_to(previous)): p for p in fresh["sources"]
                if Path(p["path"]).is_relative_to(previous / "orthohmm")}
    paths = [*sorted((source / "orthohmm").rglob("*.py")),
             *sorted((source / "orthohmm/search/csrc").glob("*.c")),
             *sorted((source / "orthohmm/search/csrc").glob("*.so"))]
    actual = {str(path.relative_to(source)): record(path) for path in paths}
    if not original or set(original) != set(actual) or any(
            (p["bytes"], p["sha256"]) != (actual[k]["bytes"], actual[k]["sha256"]) for k, p in original.items()):
        raise ValueError("Scientific code or CPU libraries changed after the fresh gate")
    return list(actual.values())


def control_readback(repo, table, parent):
    methods = [m for m in table["methods"] if m["key"] == "canonical_collection:canonical"]
    if len(methods) != 1:
        raise ValueError("Require a unique admitted canonical control")
    cells = [c for c in methods[0]["cells"] if c["dataset"] == "OrthoBench"
             and c["status"] == "independent_content_readback_verified"]
    if len(cells) != 1:
        raise ValueError("Require a unique verified OrthoBench canonical cell")
    cell = cells[0]
    for pin in cell["checked"]:
        check(pin)
    phylo = Path(cell["products"][PRODUCTS[0]]["path"]).parent
    folder = Path(parent["output"]) / "OrthoBench/canonical"
    preflight_pin = record(folder / "preflight.json")
    if phylo != folder / "native/orthohmm_phylogeny" or preflight_pin not in cell["checked"]:
        raise ValueError("Admitted control belongs to another planned native run")
    preflight = json.loads(Path(preflight_pin["path"]).read_text())
    row = {"outputs": list(cell["products"].values()),
        "coverage": {"all_candidate_genes_preserved": True,
                     "candidate_genes": cell["coverage"]["genes"], "prediction_genes": cell["coverage"]["genes"]},
        "scores": {"OrthoBench F1": cell["scores"]["OrthoBench: OrthoBench F1"]}}
    native = replay_native(row, {"command": preflight["native"]})
    scores, details, pins = readback_scores(repo, "OrthoBench", row, native, folder, "unused_for_orthobench")
    if (endpoint_vector("OrthoBench", scores) != cell["scores"] or details != cell["details"]
            or {k: native[k] for k in ("genes", "species", "groups", "pairs")} != cell["coverage"]
            or native["products"] != cell["products"] or cell["source_revision"] != parent["source_revision"]):
        raise ValueError("Independently revalidated OrthoBench control differs from admitted cell")
    checked = [preflight_pin, *cell["checked"], *native["checked"], *pins]
    for pin in checked:
        check(pin)
    return {"scores": scores, "details": details, "coverage": cell["coverage"],
        "products": native["products"], "species_tree": record(phylo / "species_tree.rooted.nwk"),
        "source_revision": cell["source_revision"], "resources": cell["resources"],
        "checked": list({p["path"]: p for p in checked}.values()),
        "scope": "Retained min4/cap100 canonical accepted-hit replay control; not fresh matched-resource timing."}


def prerequisites(repo, source, revision):
    accounting = terminal_dependencies()
    gate, gate_pin = committed(repo, source, revision, FRESH_GATE)
    admission, admission_pin = committed(repo, source, revision, QFO_ADMISSION)
    launch, launch_pin = committed(repo, source, revision, CONTROL_LAUNCH)
    if (gate.get("status") != "fresh_protected_union_all_144_directions_independently_verified"
            or gate.get("advance_to_orthobench_pilot") is not True or gate.get("directions") != 144
            or gate.get("raw_outputs_independently_recounted") != 288
            or gate.get("observed_cpu_scored_candidates") != 4560971
            or admission.get("status") != "shared_qfo_all_six_scores_admitted_and_independently_crosschecked"
            or admission.get("complete_eight_endpoint_contrast") is not True
            or admission.get("native_inference_rerun") is not False
            or gate.get("defaults_promoted") is not False or admission.get("defaults_promoted") is not False):
        raise ValueError("Fresh equivalence or prior QfO closure is not independently admitted")
    pins = [gate_pin, admission_pin, launch_pin, launch["execution_plan"],
        *[gate[k] for k in ("producer", "execution_plan", "kernel_routes", "time")],
        *[admission[k] for k in ("summary", "score_table", "producer", "normalization", "reporter_tests")]]
    for pin in pins:
        check(pin)
    fresh = json.loads(Path(gate["execution_plan"]["path"]).read_text())
    result = json.loads(Path(gate["producer"]["path"]).read_text())
    if (result.get("status") != FRESH_SUCCESS or result.get("postflight_pins_verified") is not True
            or result["execution_plan"] != gate["execution_plan"]
            or result["source_revision"] != fresh["source_revision"] or result["source_revision"] != gate["source_revision"]
            or len(result.get("directions", [])) != 144 or result.get("actual_kernel_calls") != 288
            or result.get("scored_candidates") != 4560971 or result.get("fresh_initial_search_executed") is not True):
        raise ValueError("Original fresh-search producer changed after admission")
    pins += [*fresh["inputs"], *fresh["sources"], *fresh["runtime"]["files"],
        *[r[label]["raw"] for r in result["directions"] for label in ("control", "protected_union")]]
    parent = json.loads(Path(launch["execution_plan"]["path"]).read_text())
    if (parent["repo"] != str(repo) or parent["source_root"] != launch["frozen_source"]
            or parent["source_revision"] != launch["source_revision"] or runtime() != parent["runtime"]):
        raise ValueError("Matched-control source or scientific runtime changed")
    pins += [*parent["sources"], *parent["runtime"]["files"]]
    table = json.loads(Path(admission["summary"]["path"]).read_text())
    control = control_readback(repo, table, parent)
    if (control["coverage"]["genes"], control["coverage"]["species"]) != (251378, 12):
        raise ValueError("Wrong complete OrthoBench control universe")
    pins += control["checked"]
    for pin in list({p["path"]: p for p in pins}.values()):
        check(pin)
    return parent, fresh, control, list({p["path"]: p for p in pins}.values()), accounting


def command(source, fasta, folder, tools, routes, cap=100, cpu=16):
    argv = native_command(source, fasta, folder, tools, cap, cpu)
    return [*argv, "--search-result-order", "canonical", "--search-candidate-policy", "protected_union",
            "--native-kernel-trace", str(routes)]


def validate_plan(plan, parent):
    cases = [c for c in parent["cases"] if c["dataset"] == "OrthoBench"]
    if len(cases) != 1:
        raise ValueError("Ambiguous original OrthoBench input case")
    old, case = cases[0], plan["case"]
    tools = [record(shutil.which(old[k]) or old[k]) for k in ("aligner", "tree_builder")]
    inputs = {p["path"]: p for p in parent["inputs"]}
    root, source = Path(plan["output"]), Path(plan["source_root"])
    expected = command(source, Path(old["fasta_directory"]), root / "native",
                       [p["path"] for p in tools], root / "native_routes")
    if (plan["cpu"] != 16 or plan["tools"] != tools or any(inputs.get(p["path"]) != p for p in tools)
            or case["dataset"] != "OrthoBench" or case["label"] != "protected_union"
            or case["fasta"] != old["fasta_directory"] or case["inputs"] != old["fastas"]
            or case["output"] != str(root / "native") or plan["native_routes"] != str(root / "native_routes")
            or case["command"] != expected or case["search_candidate_policy"] != "protected_union"
            or case["search_min_total_hits"] != 4 or case["search_max_candidates_per_query"] != 100):
        raise ValueError("Pilot plan differs from original inputs, tools or prespecified intervention")


def pilot_gate(scores, details, control):
    endpoint_vector("OrthoBench", scores)
    endpoint_vector("OrthoBench", control["scores"])
    for axes in (details, control["details"]):
        if any(isinstance(axes[k], bool) or not isinstance(axes[k], (int, float))
               or not math.isfinite(axes[k]) or not 0 <= axes[k] <= 1 for k in ("precision", "recall")):
            raise ValueError("Require normalized finite OrthoBench precision and recall")
    deltas = {"F1": scores["OrthoBench F1"] - control["scores"]["OrthoBench F1"],
              **{k: details[k] - control["details"][k] for k in ("precision", "recall")}}
    passed = (deltas["F1"] > 0 and all(details[k] >= control["details"][k] - .005 for k in ("precision", "recall")))
    return {"passes_development_gate": passed, "deltas": deltas,
        "scope": "Strict F1 gain with neither precision nor recall down more than 0.005; development screen only, not equivalence or default promotion."}


def prepare(repo, source, output):
    if source == repo or output.exists():
        raise ValueError("Require separate frozen source and a fresh full-pilot namespace")
    revision = subprocess.check_output(["git", "-C", str(source), "rev-parse", "HEAD"], text=True).strip()
    parent, fresh, control, pins, accounting = prerequisites(repo, source, revision)
    scientific_source(source, fresh)
    rows = [c for c in parent["cases"] if c["dataset"] == "OrthoBench"]
    if len(rows) != 1:
        raise ValueError("Ambiguous complete OrthoBench input case")
    old = rows[0]
    tools = [record(shutil.which(old[k]) or old[k]) for k in ("aligner", "tree_builder")]
    input_index = {p["path"]: p for p in parent["inputs"]}
    if any(input_index.get(p["path"]) != p for p in (*tools, *old["fastas"])):
        raise ValueError("Input proteomes or native tools differ from matched control")
    folder, routes = output / "native", output / "native_routes"
    case = {"dataset": "OrthoBench", "label": "protected_union", "output": str(folder),
        "fasta": old["fasta_directory"], "inputs": old["fastas"], "search_candidate_policy": "protected_union",
        "search_min_total_hits": 4, "search_max_candidates_per_query": 100,
        "golden": {**control["products"], "species_tree.rooted.nwk": control["species_tree"]},
        "command": command(source, Path(old["fasta_directory"]), folder, [p["path"] for p in tools], routes)}
    verify_input_inventory(case)
    source_pins = [*sources(source), record(source / "benchmark_tools" / PROTOCOL), record(source / "benchmark_tools" / SLURM)]
    for pin in source_pins:
        relative = Path(pin["path"]).relative_to(source)
        if relative.suffix != ".so":
            blob = subprocess.check_output(["git", "-C", str(source), "show", f"{revision}:{relative}"])
            if hashlib.sha256(blob).hexdigest() != pin["sha256"]:
                raise ValueError("Pilot source differs from committed frozen content")
    for pin in (*pins, *source_pins, *tools):
        check(pin)
    plan = {"status": "prepared_full_orthobench_protected_union_pilot", "repo": str(repo), "output": str(output),
        "source_root": str(source), "source_revision": revision, "sources": source_pins,
        "inputs": list({p["path"]: p for p in (*pins, *tools, *old["fastas"])}.values()),
        "runtime": parent["runtime"], "case": case, "tools": tools, "cpu": 16, "control": control,
        "native_routes": str(routes), "dependency_accounting": accounting,
        "profile_runtime": require_profile_runtime(source), "defaults_promoted": False, "goal_completed": False,
        "scope": "Single fresh full-proteome OrthoBench initial-HMM/graph/profile/phylogeny pilot on shared Threadripper; cached accuracy control has unmatched timing scope."}
    validate_plan(plan, parent)
    output.mkdir(parents=True)
    write_json(output / "execution_plan.json", plan)
    return plan


def validate_variant(plan):
    case, source = plan["case"], Path(plan["source_root"])
    native = validate_native(case, plan["cpu"])
    meta, harness = native["metrics"]["metadata"], native["metrics"]["harness"]
    if (meta.get("search_result_order") != "canonical" or meta.get("search_candidate_policy") != "protected_union"
            or meta.get("singleton_tie_policy", "first_hit") != "first_hit"
            or meta.get("gene_tree_support_mode", "legacy") != "legacy"
            or (native["genes"], native["species"]) != (plan["control"]["coverage"]["genes"], plan["control"]["coverage"]["species"])):
        raise ValueError("Actual intervention or complete gene universe differs")
    expected_sources = [*sorted((source / "orthohmm").rglob("*.py")), source / "benchmark_tools/benchmark_production.py"]
    manifest = harness["source_manifest"]
    if (harness.get("git_commit") != plan["source_revision"]
            or len({p["path"] for p in manifest}) != len(manifest)
            or {p["path"]: (p["bytes"], p["sha256"]) for p in manifest}
            != {str(p.relative_to(source)): (record(p)["bytes"], record(p)["sha256"]) for p in expected_sources}
            or argument(harness["command"], "--refinement_profile") != "default"):
        raise ValueError("Actual native scientific source or refinement profile differs")
    expected_inputs = {Path(p["path"]).name: (p["bytes"], p["sha256"]) for p in case["inputs"]}
    if (len(harness["input_manifest"]) != len(expected_inputs)
            or {p["path"]: (p["bytes"], p["sha256"]) for p in harness["input_manifest"]} != expected_inputs):
        raise ValueError("Actual native FASTA manifests differ from trial")
    library = record(source / "orthohmm/search/csrc/hmm_viterbi.so")
    routes = route_totals(Path(plan["native_routes"]), library, native["metrics"]["counts"])
    output = Path(case["output"]) / "output"
    pins = [record(Path(case["output"]) / "metrics.json"), *native["products"].values(), *routes["files"]]
    for item in harness["output_manifest"]:
        pin = record(output / item["path"])
        if (pin["bytes"], pin["sha256"]) != (item["bytes"], item["sha256"]):
            raise ValueError("Completed native output manifest changed")
        pins.append(pin)
    native.update(prediction=str(output / "orthohmm_working_res/orthohmm_edges_clustered.txt"),
                  phylogeny=str(output / "orthohmm_phylogeny"))
    return native, routes, pins


def run(path):
    plan_pin = record(path)
    plan = json.loads(path.read_text())
    source, repo, output = (Path(plan[k]) for k in ("source_root", "repo", "output"))
    if (Path(__file__).resolve() != source / "benchmark_tools/run_protected_union_pilot.py"
            or path.resolve() != output / "execution_plan.json" or (output / "results.json").exists()):
        raise ValueError("Require frozen runner and preserve every pilot attempt")
    if (os.environ.get("SLURM_CPUS_PER_TASK"), os.environ.get("SLURM_MEM_PER_NODE"),
            os.environ.get("SLURMD_NODENAME"), os.environ.get("CUDA_VISIBLE_DEVICES")) != ("16", "131072", "bizon", ""):
        raise ValueError("Require CPU-only 16-CPU 128-GiB local allocation")
    result = {"status": "validating_full_orthobench_pilot", "execution_plan": plan_pin,
        "source_revision": plan["source_revision"], "job_id": os.environ.get("SLURM_JOB_ID"), "scores": {},
        "defaults_promoted": False, "goal_completed": False, "host_before": host_activity(), "scope": plan["scope"]}
    write_json(output / "results.json", result)
    checked = [plan_pin, *plan["inputs"], *plan["sources"], *plan["runtime"]["files"]]
    started = time.monotonic()
    try:
        for pin in checked:
            check(pin)
        if runtime() != plan["runtime"]:
            raise ValueError("Matched scientific runtime changed")
        revision = subprocess.check_output(["git", "-C", str(source), "rev-parse", "HEAD"], text=True).strip()
        if revision != plan["source_revision"] or plan["cpu"] != 16:
            raise ValueError("Frozen executor revision or pilot allocation changed")
        parent, fresh, control, pins, _ = prerequisites(repo, source, revision)
        scientific_source(source, fresh)
        if control != plan["control"]:
            raise ValueError("Matched control changed since pilot preparation")
        validate_plan(plan, parent)
        checked += pins
        verify_input_inventory(plan["case"])
        folder = Path(plan["case"]["output"])
        folder.mkdir(parents=True, exist_ok=False)
        result.update(status="fresh_full_orthobench_native_running", control=control)
        write_json(output / "results.json", result)
        env = {**os.environ, "PYTHONPATH": str(source), "PYTHONHASHSEED": "0", "OMP_NUM_THREADS": "1",
            "NUMBA_NUM_THREADS": "16", "OPENBLAS_NUM_THREADS": "1", "MKL_NUM_THREADS": "1"}
        result["execution"] = timed_stage(plan["case"]["command"], folder, "native", source, env)
        if result["execution"]["returncode"]:
            raise RuntimeError("Fresh full OrthoBench protected-union pipeline failed")
        native, routes, native_pins = validate_variant(plan)
        checked += native_pins
        result.update(status="native_verified_scoring_pending", coverage={k: native[k] for k in ("genes", "species", "groups", "pairs")},
            native_products=native["products"], native_summary=native["summary"], route_evidence=routes)
        write_json(output / "results.json", result)
        scored = score_orthobench(repo, Path(native["prediction"]))
        scores, details, score_pins = readback_scores(repo, "OrthoBench", scored, native, folder, "unused_for_orthobench")
        checked += score_pins
        result.update(status=SUCCESS, scores=scores, details=details, scoring=scored,
                      gate=pilot_gate(scores, details, control))
    except Exception as error:
        result.update(status="failed_full_orthobench_protected_union_pilot", error_type=type(error).__name__, error=str(error))
    finally:
        try:
            for pin in list({p["path"]: p for p in checked}.values()):
                check(pin)
            if runtime() != plan["runtime"]:
                raise ValueError("Runtime changed during full pilot")
            result["postflight_pins_verified"] = True
        except Exception as error:
            result.update(status="failed_full_orthobench_protected_union_pilot", postflight_error=str(error))
        result.update(total_pilot_wall_s=time.monotonic() - started, host_after=host_activity())
        write_json(output / "results.json", result)
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="operation", required=True)
    preparation = sub.add_parser("prepare")
    for name in ("repo", "source", "output"):
        preparation.add_argument("--" + name, type=Path, required=True)
    execution = sub.add_parser("run")
    execution.add_argument("--execution-plan", type=Path, required=True)
    args = parser.parse_args()
    result = (prepare(args.repo.resolve(), args.source.resolve(), args.output.resolve())
              if args.operation == "prepare" else run(args.execution_plan.resolve()))
    print(result["status"])
    raise SystemExit(0 if result["status"] in ("prepared_full_orthobench_protected_union_pilot", SUCCESS) else 1)
