"""Transfer admitted projection arms without importing newer scientific modules.

Execute this file directly. Authentication runs in the frozen planner process;
only then may the scientific namespace bind to the original recovery checkout.
"""

import argparse
from copy import deepcopy
import hashlib
import importlib
import importlib.util
import json
import os
from pathlib import Path
import resource
import subprocess
import sys
import time


SCIENCE = "a6a56ac50ef7d89a2b70428b28914f8bcdce876d"
PREPARER = "98e05ee6d4ce11fad138eb59f1fcd8fcce2a994e"
PROTOCOL = "benchmark_tools/hmm_projection_followup_20261019.json"
PREPARED = "benchmarks/results/projection_followup_inputs_v1_20261019/prepared_inputs.json"
PLANNER = "benchmarks/work/projection_followup_plan_v1_20261019"
UNIVERSES = {"QfO": (984137, 78), "ThreeKingdoms": (443217, 12)}
ARMS = ("legacy_off", "projection_off", "projection_builtin_reselect",
        "projection_hmmer_score_only", "projection_hmmer_reselect")
SUCCESS = "native_and_score_independently_verified"


def helper(driver):
    spec = importlib.util.spec_from_file_location("_projection_weighting_bridge",
        driver / "benchmark_tools/run_profile_weighting_followup.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def science_api(source, driver):
    api = helper(driver).science_api(source)
    for name, symbols in {
        "benchmark_tools.run_profile_weighting_pilot": ("recovery_configuration",),
        "benchmark_tools.native_pair_turnover": ("native_pair_turnover",),
    }.items():
        module = importlib.import_module(name)
        if not Path(module.__file__).resolve().is_relative_to(source):
            raise ValueError("Projection scientific helper resolved outside original checkout")
        for symbol in symbols: setattr(api, symbol, getattr(module, symbol))
    return api


def authenticate_prepared(repo, prepared):
    if (prepared.get("status") != "projection_transfer_inputs_prepared_execution_driver_required"
            or prepared.get("scientific_revision") != SCIENCE or prepared.get("preparer_revision") != PREPARER
            or prepared.get("cpu") != 16 or prepared.get("repo") != str(repo)
            or prepared.get("output") != str((repo / PREPARED).parent)
            or prepared.get("defaults_promoted") is not False or prepared.get("goal_completed") is not False):
        raise ValueError("Require original prepared transfer, not a partial, failed or no-followup receipt")
    planner = repo / PLANNER
    helper(Path(__file__).resolve().parents[1]).successful_job("24517", 1, "16G")
    if subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=planner, text=True).strip() != PREPARER:
        raise ValueError("Frozen original preparer revision changed")
    # Keep current orchestration imports out of the original scientific process.
    code = """import hashlib,json,sys
from pathlib import Path
from types import SimpleNamespace
from benchmark_tools import plan_projection_followup as original
repo=Path(sys.argv[1]); driver=Path.cwd()
r,p,s,pins,provenance=original.authenticate(repo)
api=SimpleNamespace(record=original.record,check=original.check,verified_cell=original.verified_cell,PRODUCTS=original.PRODUCTS)
cases=[]
for dataset in original.DATASETS if s['execution_arms'] else ():
    case,extra=original.selected_case(api,repo,dataset); cases.append(case); pins+=extra
    pins += [original.record(path) for path in sorted(Path(case['checkpoint']).iterdir()) if path.is_file()]
    phylo=Path(case['case']['baseline_output'])/'orthohmm_phylogeny'
    for name in ('alignments','gene_trees','checkpoints','species_tree_inference'):
        pins += [original.record(path) for path in sorted((phylo/name).rglob('*')) if path.is_file()]
protocol=original.record(driver/'benchmark_tools'/original.PROTOCOL)
files=original.sources(driver)+[protocol]
pins=original.campaign.unique_pins(pins+files)
evidence=dict(inputs=pins,cases=cases,preparer_files=files,protocol=protocol)
digest=hashlib.sha256(json.dumps(evidence,sort_keys=True,separators=(',',':')).encode()).hexdigest()
print(json.dumps(dict(selection=s,pilot_readback=provenance['reader']['reader'],pilot_plan=r['execution_plan'],scientific_source=p['source_root'],runtime=p['runtime'],tool_path=p['tool_path'],scientific_sources=p['sources'],hmmer_tools=p['hmmer_tools'],prepared_evidence_digest=digest)))
"""
    env = {**os.environ, "PYTHONPATH": str(planner), "CUDA_VISIBLE_DEVICES": "",
        "OMP_NUM_THREADS": "1", "OPENBLAS_NUM_THREADS": "1", "MKL_NUM_THREADS": "1", "NUMBA_NUM_THREADS": "1"}
    audit = json.loads(subprocess.check_output([sys.executable, "-c", code, str(repo)], cwd=planner, env=env, text=True))
    evidence = {key: prepared[key] for key in ("inputs", "cases", "preparer_files", "protocol")}
    digest = hashlib.sha256(json.dumps(evidence, sort_keys=True, separators=(",", ":")).encode()).hexdigest()
    if digest != audit.pop("prepared_evidence_digest"):
        raise ValueError("Prepared source/data/cache inventory differs from complete original inputs")
    for key in ("selection", "pilot_readback", "pilot_plan", "scientific_source", "runtime", "tool_path"):
        if prepared[key] != audit[key]: raise ValueError("Prepared transfer differs from reauthenticated " + key)
    if not audit["selection"]["execution_arms"]:
        raise ValueError("No original candidate or applicability follow-up requires inference")
    return audit


def scoring_id(output, label):
    return "pr_" + hashlib.sha256(str(output / label).encode()).hexdigest()[:12]


def validate_plan(plan):
    chosen, dataset, output = plan["selection"], plan["dataset"], Path(plan["output"])
    arms = chosen["execution_arms"]
    if (plan["status"] != "prepared_projection_followup_unrun" or plan["scientific_revision"] != SCIENCE
            or plan["cpu"] != 16 or dataset not in UNIVERSES or not arms
            or arms != [label for label in ARMS if label in arms] or arms[:2] != list(ARMS[:2])
            or plan["arms"] != arms or set(chosen["roles"]) != set(arms)
            or set(chosen["control_dependencies"]) != set(arms)
            or plan.get("defaults_promoted") is not False or plan.get("goal_completed") is not False):
        raise ValueError("Changed fixed projection transfer source, universe, roles, order or budget")
    expected = {"legacy_off": [], "projection_off": ["legacy_off"],
        "projection_builtin_reselect": ["legacy_off", "projection_off"],
        "projection_hmmer_score_only": ["projection_off"],
        "projection_hmmer_reselect": ["legacy_off", "projection_off", "projection_builtin_reselect", "projection_hmmer_score_only"]}
    if chosen["control_dependencies"] != {label: expected[label] for label in arms}:
        raise ValueError("Changed matched-control dependencies")
    if any(dependency not in arms[:arms.index(label)] for label in arms for dependency in expected[label]):
        raise ValueError("Missing or out-of-order matched control")
    identifiers = {label: scoring_id(output, label) for label in arms if label not in ("legacy_off", "projection_hmmer_score_only")}
    if plan["scoring_ids"] != identifiers or len(set(identifiers.values())) != len(identifiers):
        raise ValueError("Changed or colliding assessment namespaces")
    case = plan["case"]
    if (case["dataset"] != dataset or (case["control"]["coverage"]["genes"], case["control"]["coverage"]["species"]) != UNIVERSES[dataset]
            or len(case["case"]["fastas"]) != UNIVERSES[dataset][1]
            or any(Path(pin["path"]).parent != Path(case["case"]["fasta_directory"]) for pin in case["case"]["fastas"])
            or any(Path(pin["path"]).parent != Path(case["case"]["baseline_output"]) / "orthohmm_phylogeny" for pin in case["golden"].values())
            or Path(case["checkpoint_manifest"]["path"]) != Path(case["checkpoint"]) / "manifest.json"
            or case["checkpoint_manifest"]["sha256"] != case["checkpoint_sha"]):
        raise ValueError("Changed complete transfer case or checkpoint")
    needed = [plan["prepared_inputs"], plan["pilot_readback"], plan["pilot_plan"], case["checkpoint_manifest"],
        *case["case"]["fastas"], *case["golden"].values(), case["candidate"]["candidate_partition"],
        case["candidate"]["membership_constraints"], *plan["scientific_sources"], *plan["runtime"]["files"], *plan["hmmer_tools"].values()]
    if any(pin not in plan["inputs"] for pin in needed):
        raise ValueError("Projection transfer is missing authenticated source/data/tool inputs")


def prepare(repo, driver, output, dataset):
    if output.exists(): raise FileExistsError("Preserve every projection execution namespace")
    if driver == repo or Path(__file__).resolve() != driver / "benchmark_tools/run_projection_followup.py":
        raise ValueError("Require frozen standalone projection execution driver")
    prepared_path = repo / PREPARED
    prepared = json.loads(prepared_path.read_text())
    audit = authenticate_prepared(repo, prepared)
    science = Path(audit["scientific_source"])
    api = science_api(science, driver)
    if dataset not in UNIVERSES or [case["dataset"] for case in prepared["cases"]] != list(UNIVERSES):
        raise ValueError("Require fixed QfO then Three Kingdoms transfer cases")
    case, _ = helper(driver).selected_case(api, repo, dataset)
    if case != next(c for c in prepared["cases"] if c["dataset"] == dataset):
        raise ValueError("Prepared case differs from original admitted canonical case")
    revision = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=driver, text=True).strip()
    driver_files = [api.record(driver / path) for path in ("benchmark_tools/run_projection_followup.py",
        "benchmark_tools/run_projection_followup.slurm", "benchmark_tools/run_profile_weighting_followup.py", PROTOCOL)]
    for pin in driver_files:
        blob = subprocess.check_output(["git", "show", revision + ":" + str(Path(pin["path"]).relative_to(driver))], cwd=driver)
        if hashlib.sha256(blob).hexdigest() != pin["sha256"]:
            raise ValueError("Projection driver differs from committed source")
    if api.record(driver / PROTOCOL)["sha256"] != prepared["protocol"]["sha256"]:
        raise ValueError("Transfer protocol changed after pilot selection")
    plan = {"status": "prepared_projection_followup_unrun", "repo": str(repo), "output": str(output),
        "driver_root": str(driver), "driver_revision": revision, "driver_sources": driver_files,
        "scientific_revision": SCIENCE, "prepared_inputs": api.record(prepared_path), "cpu": 16,
        **audit, "case": case, "dataset": dataset, "arms": prepared["selection"]["execution_arms"],
        "scoring_ids": {label: scoring_id(output, label) for label in prepared["selection"]["execution_arms"]
                        if label not in ("legacy_off", "projection_hmmer_score_only")},
        "inputs": [api.record(prepared_path), *prepared["inputs"]],
        "defaults_promoted": False, "goal_completed": False,
        "scope": "Fixed source/selection/controls; complete-proteome cached-search and compatible native-cache replay with inferred phylogeny. Shared-host incremental timings, not isolated fresh-pipeline costs. Assessment reuse is labeled and not an independent replicate."}
    validate_plan(plan)
    for pin in [*plan["inputs"], *driver_files]: api.check(pin)
    if api.runtime() != plan["runtime"]: raise ValueError("Scientific runtime changed")
    if dataset == "QfO":
        for identifier in plan["scoring_ids"].values():
            work = repo / "qfo_benchmark/w" / ("hc1_" + identifier)
            api.require_darwin_work_path(work)
            if work.exists(): raise FileExistsError("Preserve existing QfO assessment namespace")
    plan["profile_runtime"] = api.require_profile_runtime(science)
    output.mkdir(parents=True)
    api.write_json(output / "execution_plan.json", plan)
    return plan


def arm_commands(api, plan, folder, label):
    case, science = plan["case"], Path(plan["scientific_source"])
    graph, native = api.commands(science, plan, case["case"], folder, Path(case["checkpoint"]), case["checkpoint_sha"])
    policy, mode = api.recovery_configuration(label)
    if policy != "legacy": graph += ["--profile-projection-policy", policy]
    if mode != "off": graph += ["--selective-profile-mode", mode]
    return graph, native


def candidate_equivalent(row, reference):
    for key in ("candidate_partition", "membership_constraints"):
        if any(row["candidate"][key][field] != reference["candidate"][key][field] for field in ("bytes", "sha256")):
            raise ValueError("Exact control candidate or constraints differ")


def score_arm(api, plan, row, native, folder, projection):
    repo, dataset, label = Path(plan["repo"]), plan["dataset"], row["label"]
    if label == "legacy_off": return api.control_scores(api, plan, row, native, folder)
    if label == "projection_hmmer_score_only":
        if projection is None: raise ValueError("Measurement scoring requires verified repaired control")
        for key in ("scores", "conversion", "metric_files", "references"):
            if key in projection: row[key] = deepcopy(projection[key])
        original = Path(plan["output"]) / "projection_off"
        row["measurement_scoring_reuse"] = api.record(original / "assessment_receipt.json")
        participant = "ohmm_hillclimb_r1_" + plan["scoring_ids"]["projection_off"] if dataset == "QfO" else None
        actual, details, pins = api.readback_scores(repo, dataset, row, native, original, participant)
        if (actual != projection["scores"] or details != projection["details"]
                or api.endpoint_vector(dataset, actual) != projection["metric_vector"]):
            raise ValueError("Measurement assessment differs from repaired control")
        return actual, details, [row["measurement_scoring_reuse"], *pins]
    if dataset == "QfO":
        identifier = plan["scoring_ids"][label]
        row.update(api.score_qfo(repo, Path(native["phylogeny"]), folder, identifier))
        participant = "ohmm_hillclimb_r1_" + identifier
    else:
        row.update(api.score_three_kingdoms(repo, Path(native["prediction"]), folder))
        participant = None
    return api.readback_scores(repo, dataset, row, native, folder, participant)


def execute_arms(api, plan, path, names, species, hits, ownership):
    output, source = Path(plan["output"]), Path(plan["scientific_source"])
    pins = [api.record(path), *plan["inputs"], *plan["driver_sources"]]
    result = {"status": "running", "execution_plan": api.record(path), "job_id": os.environ.get("SLURM_JOB_ID"),
        "source_revision": SCIENCE, "driver_revision": plan["driver_revision"], "dataset": plan["dataset"],
        "selection": plan["selection"], "ownership": ownership, "arms": [], "host_before": api.host_activity(),
        "defaults_promoted": False, "goal_completed": False, "scope": plan["scope"]}
    env = {**os.environ, "PATH": plan["tool_path"], "PYTHONPATH": str(source), "PYTHONHASHSEED": "0",
        "CUDA_VISIBLE_DEVICES": "", "OMP_NUM_THREADS": "1", "NUMBA_NUM_THREADS": "16", "OPENBLAS_NUM_THREADS": "1", "MKL_NUM_THREADS": "1"}
    completed, projection, digest = {}, None, None
    for label in plan["arms"]:
        row = {"label": label, "dataset": plan["dataset"], "role": plan["selection"]["roles"][label], "status": "validating_inputs", "scores": {}}
        result["arms"].append(row)
        missing = [name for name in plan["selection"]["control_dependencies"][label] if name not in completed]
        if missing:
            row.update(status="not_run_failed_required_control", missing_controls=missing)
            api.write_json(output / "results.json", result)
            continue
        folder = output / label
        folder.mkdir(exist_ok=False)
        started = time.monotonic()
        api.write_json(output / "results.json", result)
        try:
            for pin in pins: api.check(pin)
            if api.runtime() != plan["runtime"]: raise ValueError("Runtime changed before arm")
            graph, native_command = arm_commands(api, plan, folder, label)
            api.write_json(folder / "preflight.json", {"graph": graph, "native": native_command})
            row["status"] = "graph_running"
            api.write_json(output / "results.json", result)
            row["graph_execution"] = api.timed_stage(graph, folder, "graph", source, env)
            if row["graph_execution"]["returncode"]: raise RuntimeError("Projection graph replay failed")
            metrics = folder / "graph_metrics.json"
            graph_report = json.loads(metrics.read_text())
            row["activation"] = api.validate_activation(graph_report, label, "projection_recovery",
                folder=folder, gene_names=names, hmmer_tools=plan["hmmer_tools"])
            signature = [(stage["label"], stage["output"]["bytes"], stage["output"]["sha256"]) for stage in graph_report["stages"]]
            if label not in ("legacy_off", "projection_off"):
                candidate_digest = row["activation"]["candidate_digest"]
                if digest is not None and digest != candidate_digest: raise ValueError("Active candidate generation differs")
                digest = candidate_digest
            if label == "projection_hmmer_score_only" and signature != projection["graph_stage_signature"]:
                raise ValueError("Measurement-only graph differs from repaired control")
            row.update(graph_report=api.record(metrics), graph_stage_signature=signature)
            for stage in graph_report["stages"]:
                api.check(stage["output"])
                api.read_partition(Path(stage["output"]["path"]), set(names))
            before, candidate_started = resource.getrusage(resource.RUSAGE_SELF), time.monotonic()
            row["candidate"] = api.prepare_partition(folder / "graph/orthogroups_profiles_refined.txt", folder / "candidate",
                names, species, hits, True, api._expand_phylogeny_candidates, api.load_membership_constraints)
            after = resource.getrusage(resource.RUSAGE_SELF)
            row["candidate_resources"] = {"wall_s": time.monotonic() - candidate_started,
                "cpu_s": after.ru_utime + after.ru_stime - before.ru_utime - before.ru_stime,
                "peak_parent_process_rss_GiB": after.ru_maxrss / 1024 ** 2,
                "memory_scope": "Parent lifetime maximum, not simultaneous process-tree sum"}
            reference = plan["case"] if label == "legacy_off" else projection if label == "projection_hmmer_score_only" else None
            if reference is not None: candidate_equivalent(row, reference)
            row["status"] = "native_running"
            api.write_json(output / "results.json", result)
            row["native_execution"] = api.timed_stage(native_command, folder, "native", source, env)
            if row["native_execution"]["returncode"]: raise RuntimeError("Projection native replay failed")
            phylo = folder / "native/orthohmm_phylogeny"
            row["coverage"] = api.validate_partition(Path(row["candidate"]["candidate_partition"]["path"]),
                folder / "native/orthohmm_working_res/orthohmm_edges_clustered.txt")
            row["outputs"] = [api.record(phylo / name) for name in api.PRODUCTS]
            row["rooted_species_tree"] = api.record(phylo / "species_tree.rooted.nwk")
            native = api.replay_native(row, {"command": native_command})
            if (native["genes"], native["species"]) != UNIVERSES[plan["dataset"]]:
                raise ValueError("Native complete universe differs")
            if label == "legacy_off": row["golden_matches"] = api.compare_golden(phylo, plan["case"]["golden"])
            elif label == "projection_hmmer_score_only":
                reference = {"candidate": projection["candidate"], "golden": {Path(p["path"]).name: p for p in [*projection["outputs"], projection["rooted_species_tree"]]}}
                row["golden_matches"] = api.compare_golden(phylo, reference["golden"])
            row["native_pair_turnover"] = api.native_pair_turnover(
                Path(plan["case"]["case"]["baseline_output"]) / "orthohmm_phylogeny/orthohmm_pairwise_orthologs_confidence.tsv",
                phylo / "orthohmm_pairwise_orthologs_confidence.tsv")
            row["status"] = "native_verified_scoring_pending"
            api.write_json(output / "results.json", result)
            actual, details, score_pins = score_arm(api, plan, row, native, folder, projection)
            if actual != row["scores"]: raise ValueError("Independent assessment scores differ")
            row.update(details=details, metric_vector=api.endpoint_vector(plan["dataset"], actual), native_summary=native["summary"],
                checked=[api.record(metrics), *native["checked"], *score_pins, *row["native_pair_turnover"]["files"]])
            for pin in [*pins, *row["checked"]]: api.check(pin)
            if api.runtime() != plan["runtime"]: raise ValueError("Runtime changed during arm")
            api.write_json(folder / "assessment_receipt.json", {key: row[key] for key in ("scores", "details", "metric_vector", "checked")})
            row["checked"].append(api.record(folder / "assessment_receipt.json"))
            row["status"] = SUCCESS
            completed[label] = row
            if label == "projection_off": projection = row
        except Exception as error:
            row.update(status="failed", error_type=type(error).__name__, error=str(error))
        finally:
            row["total_incremental_wall_s"] = time.monotonic() - started
            api.write_json(output / "results.json", result)
    result["status"] = "all_projection_followup_arms_scored" if all(row["status"] == SUCCESS for row in result["arms"]) else "finished_with_explicit_gaps"
    result["host_after"] = api.host_activity()
    api.write_json(output / "results.json", result)
    return result


def run(path):
    plan = json.loads(path.read_text())
    validate_plan(plan)
    source, driver, output, repo = (Path(plan[key]) for key in ("scientific_source", "driver_root", "output", "repo"))
    if Path(__file__).resolve() != driver / "benchmark_tools/run_projection_followup.py" or (output / "results.json").exists():
        raise ValueError("Require frozen driver and preserve prior attempts")
    prepared_path = Path(plan["prepared_inputs"]["path"])
    if prepared_path != repo / PREPARED: raise ValueError("Unexpected original preparation path")
    audit = authenticate_prepared(repo, json.loads(prepared_path.read_text()))
    if any(plan[key] != value for key, value in audit.items()): raise ValueError("Execution plan changes original authenticated transfer")
    api = science_api(source, driver)
    for pin in [api.record(path), *plan["inputs"], *plan["driver_sources"]]: api.check(pin)
    if (api.runtime() != plan["runtime"] or os.environ.get("CUDA_VISIBLE_DEVICES") != ""
            or subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=source, text=True).strip() != SCIENCE
            or subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=driver, text=True).strip() != plan["driver_revision"]
            or api.sources(source) + [api.record(source / "benchmark_tools/run_profile_weighting_pilot.slurm")] != plan["scientific_sources"]):
        raise ValueError("Original science, driver or runtime changed")
    case, _ = helper(driver).selected_case(api, repo, plan["dataset"])
    if case != plan["case"]: raise ValueError("Execution case changes original admitted canonical case")
    names, species, q, t, scores = api.load_accuracy_checkpoint(Path(case["checkpoint"]), verify=True)
    owners = {}
    for pin in case["case"]["fastas"]:
        for sequence in api.SeqIO.parse(pin["path"], "fasta"):
            if sequence.id in owners or not sequence.seq: raise ValueError("Duplicate or empty complete input gene")
            owners[sequence.id] = Path(pin["path"]).name
    ownership = api.validate_ownership(names, species, owners)
    if (len(names), len(ownership["classes"])) != UNIVERSES[plan["dataset"]]:
        raise ValueError("Canonical complete checkpoint universe differs")
    return execute_arms(api, plan, path, names, species, (q, t, scores), ownership)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    prep = commands.add_parser("prepare")
    for name in ("repo", "driver", "output"): prep.add_argument("--" + name, required=True, type=Path)
    prep.add_argument("--dataset", required=True, choices=UNIVERSES)
    worker = commands.add_parser("run")
    worker.add_argument("--execution-plan", required=True, type=Path)
    args = parser.parse_args()
    result = prepare(args.repo.resolve(), args.driver.resolve(), args.output.resolve(), args.dataset) if args.command == "prepare" else run(args.execution_plan.resolve())
    print(result["status"])
    raise SystemExit(0 if result["status"] in ("prepared_projection_followup_unrun", "all_projection_followup_arms_scored") else 1)
