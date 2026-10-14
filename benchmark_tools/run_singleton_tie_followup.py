"""Execute the fixed singleton-policy follow-up with the original scientific modules."""

import argparse
import hashlib
import importlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
from types import SimpleNamespace

REVISION = "1225e3f9b504ee68e90632bf4efdc92e8232703a"
UNIVERSES = {"QfO": (984137, 78), "ThreeKingdoms": (443217, 12)}
POLICIES = ("first_hit", "defer_ambiguous")


def science_api(source):
    source = source.resolve()
    for name, module in list(sys.modules.items()):
        if name.split(".")[0] in ("orthohmm", "benchmark_tools") and getattr(module, "__file__", None):
            if not Path(module.__file__).resolve().is_relative_to(source):
                raise ValueError("Scientific modules already imported from another source")
    sys.path.insert(0, str(source))
    exports = {
        "benchmark_tools.prepare_ob_candidate_neighborhood": ("check", "record"),
        "benchmark_tools.run_membership_hillclimb": ("PRODUCTS", "compare_golden", "score_qfo",
            "score_three_kingdoms", "validate_partition", "write_json"),
        "benchmark_tools.report_hillclimb": ("replay_native", "readback_scores", "endpoint_vector"),
        "benchmark_tools.run_canonical_downstream": ("commands",),
        "benchmark_tools.prepare_orthobench_factorial": ("prepare_partition",),
        "benchmark_tools.replay_phylogeny": ("load_membership_constraints",),
        "benchmark_tools.run_full_search_hillclimb": ("runtime",),
        "benchmark_tools.run_search_recovery": ("host_activity",),
        "benchmark_tools.run_three_kingdoms_membership_supplement": ("timed_stage",),
        "benchmark_tools.run_qfo_recovered_assessment": ("require_darwin_work_path",),
        "benchmark_tools.validate_profile_runtime": ("require_profile_runtime",),
        "orthohmm.accuracy": ("load_accuracy_checkpoint",),
        "orthohmm.orthohmm": ("_expand_phylogeny_candidates",),
    }
    api = {}
    for name, symbols in exports.items():
        module = importlib.import_module(name)
        if not Path(module.__file__).resolve().is_relative_to(source):
            raise ValueError("Scientific module resolved outside the frozen source")
        api.update({symbol: getattr(module, symbol) for symbol in symbols})
    return SimpleNamespace(**api)


def selected_case(prepared, dataset):
    gate = prepared["pilot_gate"]
    if (not 0 <= gate["control_f1"] <= 1 or not 0 <= gate["variant_f1"] <= 1
            or gate["delta"] != gate["variant_f1"] - gate["control_f1"]):
        raise ValueError("Pilot decline must reproduce from valid admitted scores")
    if (prepared["status"] != "prespecified_followup_inputs_prepared_not_submitted"
            or prepared["scientific_revision"] != REVISION or prepared["cpu"] != 16
            or not prepared["pilot_gate"]["eligible"]
            or prepared["pilot_gate"]["delta"] < -.005
            or [case["dataset"] for case in prepared["cases"]] != list(UNIVERSES)):
        raise ValueError("Require the unchanged prespecified follow-up and eligible pilot")
    if dataset not in UNIVERSES:
        raise ValueError("Unsupported singleton follow-up dataset")
    return next(case for case in prepared["cases"] if case["dataset"] == dataset)


def scoring_id(output):
    return "sf_" + hashlib.sha256(str(output).encode()).hexdigest()[:12]


def freeze_reporting_evidence(api, inputs, repo, driver):
    frozen, origins = [], []
    for pin in inputs:
        api.check(pin)
        path = Path(pin["path"])
        if path.parent == repo / "benchmark_tools" and path.suffix == ".py":
            target = driver / "evidence/benchmark_tools" / path.name
            target.parent.mkdir(parents=True, exist_ok=True)
            if not target.exists():
                shutil.copy2(path, target)
            captured = api.record(target)
            if (captured["bytes"], captured["sha256"]) != (pin["bytes"], pin["sha256"]):
                raise ValueError("Frozen reporting evidence differs from the admitted version")
            frozen.append(captured)
            origins.append({"original": pin, "frozen": captured})
        else:
            frozen.append(pin)
    return frozen, origins


def prepare(prepared_path, driver, output, dataset):
    if output.exists():
        raise FileExistsError("Preserve every prior execution namespace")
    prepared = json.loads(prepared_path.read_text())
    case = selected_case(prepared, dataset)
    science = Path(prepared["scientific_source"])
    api = science_api(science)
    api.check(prepared["protocol"])
    protocol = json.loads(Path(prepared["protocol"]["path"]).read_text())
    if protocol["datasets"] != list(UNIVERSES) or protocol["arms"] != list(POLICIES):
        raise ValueError("Wrong fixed follow-up protocol inventory")
    source_pins = [api.record(driver / "benchmark_tools" / name)
                  for name in ("run_singleton_tie_followup.py", "run_singleton_tie_followup.slurm")]
    repo = Path(prepared["repo"])
    revision = subprocess.check_output(["git", "-C", str(repo), "rev-parse", "HEAD"], text=True).strip()
    for pin in source_pins:
        relative = Path(pin["path"]).relative_to(driver)
        blob = subprocess.check_output(["git", "-C", str(repo), "show", revision + ":" + str(relative)])
        if hashlib.sha256(blob).hexdigest() != pin["sha256"]:
            raise ValueError("Execution driver differs from committed source")
    checkpoint = Path(case["checkpoint"])
    if api.record(checkpoint / "manifest.json")["sha256"] != case["checkpoint_sha"]:
        raise ValueError("Canonical accepted-hit checkpoint differs")
    inputs = [api.record(prepared_path), prepared["preparer_source"], *prepared["inputs"],
              *[api.record(path) for path in sorted(checkpoint.iterdir()) if path.is_file()]]
    phylo = Path(case["case"]["baseline_output"]) / "orthohmm_phylogeny"
    for name in ("alignments", "gene_trees", "checkpoints", "species_tree_inference"):
        inputs.extend(api.record(path) for path in sorted((phylo / name).rglob("*")) if path.is_file())
    inputs = list({pin["path"]: pin for pin in inputs}.values())
    inputs, origins = freeze_reporting_evidence(api, inputs, repo, driver)
    for pin in (*inputs, *source_pins):
        api.check(pin)
    if api.runtime() != prepared["runtime"]:
        raise ValueError("Scientific runtime differs from the pilot")
    profile_probe = api.require_profile_runtime(science)
    identifier = scoring_id(output)
    if dataset == "QfO":
        work = repo / "qfo_benchmark/w" / ("hc1_" + identifier)
        api.require_darwin_work_path(work)
        if work.exists():
            raise FileExistsError("Variant QfO scoring namespace already exists")
    plan = {"status": "prepared_singleton_followup_unrun", "repo": str(repo), "output": str(output),
            "driver_root": str(driver), "driver_revision": revision, "driver_sources": source_pins,
            "scientific_source": str(science), "scientific_revision": REVISION,
            "prepared_inputs": api.record(prepared_path), "inputs": inputs,
            "frozen_reporting_evidence": origins,
            "runtime": prepared["runtime"], "tool_path": prepared["tool_path"], "cpu": 16,
            "profile_runtime": profile_probe,
            "dataset": dataset, "case": case, "scoring_id": identifier,
            "scope": "Matched cached downstream singleton-policy contrast; not fresh HMM search or isolated timing.",
            "defaults_promoted": False, "goal_completed": False}
    output.mkdir(parents=True)
    api.write_json(output / "execution_plan.json", plan)
    return plan


def control_scores(api, plan, row, native, folder):
    repo = Path(plan["repo"])
    if plan["dataset"] == "QfO":
        recovery = repo / "benchmarks/results/canonical_qfo_scoring_recovery_r1_v1_20261014"
        path = recovery / "results.json"
        if api.record(path) not in plan["inputs"]:
            raise ValueError("Recovered QfO control assessment not bound to the execution plan")
        report = json.loads(path.read_text())
        recovery_path = recovery / "execution_plan.json"
        if report["execution_plan"] != api.record(recovery_path) or report["execution_plan"] not in plan["inputs"]:
            raise ValueError("Recovered assessment plan differs from its pinned producer")
        recovery_plan = json.loads(recovery_path.read_text())
        if report["status"] != "canonical_qfo_scoring_recovered_independently_verified":
            raise ValueError("QfO control assessment is not independently verified")
        for key in ("scores", "conversion", "metric_files", "references"):
            if key in report:
                row[key] = report[key]
        participant = "ohmm_hillclimb_r1_" + recovery_plan["scoring_id"]
        scores, details, pins = api.readback_scores(repo, "QfO", row, native, recovery, participant)
        row["control_scoring_reuse"] = api.record(path)
    else:
        row.update(api.score_three_kingdoms(repo, Path(native["prediction"]), folder))
        scores, details, pins = api.readback_scores(repo, "ThreeKingdoms", row, native, folder, None)
    if api.endpoint_vector(plan["dataset"], scores) != plan["case"]["control"]["scores"]:
        raise ValueError("Independently read-back control scores differ")
    return scores, details, pins


def run(path):
    plan = json.loads(path.read_text())
    science, driver, output = (Path(plan[key]) for key in ("scientific_source", "driver_root", "output"))
    if (Path(__file__).resolve() != driver / "benchmark_tools/run_singleton_tie_followup.py"
            or (output / "results.json").exists()):
        raise ValueError("Require the frozen standalone driver and preserve previous attempts")
    api = science_api(science)
    pins = [api.record(path), *plan["inputs"], *plan["driver_sources"]]
    for pin in pins:
        api.check(pin)
    api.check(plan["prepared_inputs"])
    prepared = json.loads(Path(plan["prepared_inputs"]["path"]).read_text())
    case = selected_case(prepared, plan["dataset"])
    if (case != plan["case"] or plan["scientific_revision"] != REVISION or plan["cpu"] != 16
            or plan["scoring_id"] != scoring_id(output)
            or plan["scientific_source"] != prepared["scientific_source"]
            or plan["runtime"] != prepared["runtime"] or plan["tool_path"] != prepared["tool_path"]):
        raise ValueError("Execution plan changes the fixed scientific contrast")
    if api.runtime() != plan["runtime"] or os.environ.get("CUDA_VISIBLE_DEVICES") != "":
        raise ValueError("Scientific runtime or CUDA visibility changed")
    checkpoint = Path(case["checkpoint"])
    names, species, q, t, scores = api.load_accuracy_checkpoint(checkpoint, verify=True)
    universe = UNIVERSES[plan["dataset"]]
    if (len(names), len(set(species.tolist()))) != universe:
        raise ValueError("Wrong complete accepted-hit universe")
    result = {"status": "running", "execution_plan": api.record(path), "job_id": os.environ.get("SLURM_JOB_ID"),
              "dataset": plan["dataset"], "source_revision": REVISION, "driver_revision": plan["driver_revision"],
              "arms": [], "host_before": api.host_activity(), "scope": plan["scope"],
              "defaults_promoted": False, "goal_completed": False}
    api.write_json(output / "results.json", result)
    env = {**os.environ, "PATH": plan["tool_path"], "PYTHONPATH": str(science), "PYTHONHASHSEED": "0",
           "OMP_NUM_THREADS": "1", "NUMBA_NUM_THREADS": "16", "OPENBLAS_NUM_THREADS": "1", "MKL_NUM_THREADS": "1"}
    control_ok = False
    for policy in POLICIES:
        row = {"dataset": plan["dataset"], "label": policy, "status": "validating_inputs", "scores": {}}
        result["arms"].append(row)
        api.write_json(output / "results.json", result)
        if policy != POLICIES[0] and not control_ok:
            row["status"] = "not_run_failed_control_gate"
            api.write_json(output / "results.json", result)
            continue
        folder = output / policy
        folder.mkdir(exist_ok=False)
        try:
            for pin in pins:
                api.check(pin)
            if api.runtime() != plan["runtime"]:
                raise ValueError("Runtime changed before arm execution")
            graph, command = api.commands(science, plan, case["case"], folder, checkpoint, case["checkpoint_sha"])
            graph.extend(["--singleton-tie-policy", policy])
            api.write_json(folder / "preflight.json", {"graph": graph, "native": command})
            row["status"] = "graph_running"
            api.write_json(output / "results.json", result)
            row["graph_execution"] = api.timed_stage(graph, folder, "graph", science, env)
            if row["graph_execution"]["returncode"]:
                raise RuntimeError("Singleton-policy graph replay failed")
            metrics_path = folder / "graph_metrics.json"
            parameters = json.loads(metrics_path.read_text())["parameters"]
            expected = {"singleton_tie_policy": policy, "cpm_resolution": .1, "leiden_seed": 4,
                        "profile_iterations": 1, "profile_min_species": 1, "matrix": "BLOSUM62"}
            if any(parameters.get(key) != value for key, value in expected.items()):
                raise ValueError("Actual singleton graph settings differ")
            row["candidate"] = api.prepare_partition(folder / "graph/orthogroups_profiles_refined.txt", folder / "candidate",
                names, species, (q, t, scores), True, api._expand_phylogeny_candidates, api.load_membership_constraints)
            if policy == POLICIES[0]:
                for key in ("candidate_partition", "membership_constraints"):
                    a, b = row["candidate"][key], case["candidate"][key]
                    if (a["bytes"], a["sha256"]) != (b["bytes"], b["sha256"]):
                        raise ValueError("Canonical candidate control does not reproduce")
            row["status"] = "native_running"
            api.write_json(output / "results.json", result)
            row["native_execution"] = api.timed_stage(command, folder, "native", science, env)
            if row["native_execution"]["returncode"]:
                raise RuntimeError("Singleton-policy native replay failed")
            phylo = folder / "native/orthohmm_phylogeny"
            prediction = folder / "native/orthohmm_working_res/orthohmm_edges_clustered.txt"
            row["coverage"] = api.validate_partition(Path(row["candidate"]["candidate_partition"]["path"]), prediction)
            row["outputs"] = [api.record(phylo / name) for name in api.PRODUCTS]
            native = api.replay_native(row, {"command": command})
            if (native["genes"], native["species"]) != universe:
                raise ValueError("Native complete gene universe differs")
            row["status"] = "native_verified_scoring_pending"
            api.write_json(output / "results.json", result)
            if policy == POLICIES[0]:
                row["golden_matches"] = api.compare_golden(phylo, case["golden"])
                actual, details, score_pins = control_scores(api, plan, row, native, folder)
            else:
                if plan["dataset"] == "QfO":
                    row.update(api.score_qfo(Path(plan["repo"]), phylo, folder, plan["scoring_id"]))
                    participant = "ohmm_hillclimb_r1_" + plan["scoring_id"]
                else:
                    row.update(api.score_three_kingdoms(Path(plan["repo"]), prediction, folder))
                    participant = None
                actual, details, score_pins = api.readback_scores(Path(plan["repo"]), plan["dataset"], row, native, folder, participant)
            checked = [api.record(metrics_path), *native["checked"], *score_pins]
            for pin in (*pins, *checked):
                api.check(pin)
            if api.runtime() != plan["runtime"]:
                raise ValueError("Runtime changed during the arm")
            row.update(status="native_and_score_independently_verified", details=details,
                       native_summary=native["summary"], checked=checked,
                       metric_vector=api.endpoint_vector(plan["dataset"], actual))
            if policy == POLICIES[0]:
                control_ok = True
        except Exception as error:
            row.update(status="failed", error_type=type(error).__name__, error=str(error))
        api.write_json(output / "results.json", result)
    result["status"] = ("both_singleton_policies_scored" if all(row["status"] == "native_and_score_independently_verified"
                        for row in result["arms"]) else "finished_with_explicit_gaps")
    if result["status"] == "both_singleton_policies_scored":
        result["deltas"] = {key: value - result["arms"][0]["metric_vector"][key]
                            for key, value in result["arms"][1]["metric_vector"].items()}
    result["host_after"] = api.host_activity()
    api.write_json(output / "results.json", result)
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    prep = sub.add_parser("prepare")
    for name in ("prepared-inputs", "driver-root", "output"):
        prep.add_argument("--" + name, type=Path, required=True)
    prep.add_argument("--dataset", choices=list(UNIVERSES), required=True)
    worker = sub.add_parser("run")
    worker.add_argument("--execution-plan", type=Path, required=True)
    args = parser.parse_args()
    result = (prepare(args.prepared_inputs.resolve(), args.driver_root.resolve(), args.output.resolve(), args.dataset)
              if args.command == "prepare" else run(args.execution_plan.resolve()))
    print(result["status"])
    raise SystemExit(0 if result["status"] in ("prepared_singleton_followup_unrun", "both_singleton_policies_scored") else 1)
