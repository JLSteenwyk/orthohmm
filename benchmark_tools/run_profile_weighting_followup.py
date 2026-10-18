"""Run unchanged admitted weighting arms with the original pilot scientific modules.

Execute this file directly, not with -m: the scientific namespace must be
bound to the old immutable checkout before importing any benchmark package.
"""

import argparse
import csv
import hashlib
import importlib
import io
import json
import os
from pathlib import Path
import resource
import subprocess
import sys
import time


REVISION = "620bfdf91d3d1b8c3a5d054df36b0f199c11de55"
ADMISSION_REVISION = "d8c6c281a4c1cfc4cfa85cda65ed66276e0d4aee"
ADMISSION = "benchmark_tools/results/profile_weighting_admission_20261018.json"
PROTOCOL = "benchmark_tools/hmm_profile_weighting_followup_20261018.json"
CONTEXT = "benchmark_tools/results/hillclimb_comparison_readback_v27_20261016/summary.json"
CANONICAL = "benchmarks/results/canonical_downstream_r1_v1_20261014"
UNIVERSES = {"QfO": (984137, 78), "ThreeKingdoms": (443217, 12)}
ARMS = ("none", "henikoff", "henikoff_cluster")


def science_api(source):
    source = source.resolve()
    for name, module in list(sys.modules.items()):
        if name.split(".")[0] in ("orthohmm", "benchmark_tools") and getattr(module, "__file__", None):
            if not Path(module.__file__).resolve().is_relative_to(source):
                raise ValueError("Scientific modules already imported from another source")
    sys.path.insert(0, str(source))
    bridge = importlib.import_module("benchmark_tools.run_singleton_tie_followup")
    if not Path(bridge.__file__).resolve().is_relative_to(source):
        raise ValueError("Wrong scientific API bridge")
    api = bridge.science_api(source)
    additions = {
        "benchmark_tools.run_profile_weighting_pilot": ("validate_activation", "advancement"),
        "benchmark_tools.plan_singleton_tie_followup": ("verified_cell",),
        "benchmark_tools.run_three_kingdoms_membership_supplement": ("validate_ownership", "sources"),
        "benchmark_tools.audit_historical_profile_ablation": ("read_partition",),
    }
    for name, symbols in additions.items():
        module = importlib.import_module(name)
        if not Path(module.__file__).resolve().is_relative_to(source):
            raise ValueError("Scientific helper resolved outside the pilot checkout")
        for symbol in symbols:
            setattr(api, symbol, getattr(module, symbol))
    api.control_scores = bridge.control_scores
    api.SeqIO = importlib.import_module("Bio.SeqIO")
    return api


def committed(api, repo, relative, revision=ADMISSION_REVISION):
    path = repo / relative
    blob = subprocess.check_output(["git", "show", revision + ":" + relative], cwd=repo)
    pin = api.record(path)
    if hashlib.sha256(blob).hexdigest() != pin["sha256"]:
        raise ValueError("Changed committed prerequisite: " + relative)
    return json.loads(blob), pin


def successful_job(job_id, cpu, memory):
    command = ["sacct", "-j", str(job_id), "-X", "-P", "-o", "JobIDRaw,State,ExitCode,AllocCPUS,ReqMem,NodeList"]
    rows = list(csv.DictReader(io.StringIO(subprocess.check_output(command, text=True)), delimiter="|"))
    expected = {"JobIDRaw": str(job_id), "State": "COMPLETED", "ExitCode": "0:0", "AllocCPUS": str(cpu),
                "ReqMem": memory, "NodeList": "bizon"}
    if rows != [expected]:
        raise ValueError("Follow-up prerequisite job is not uniquely terminal successful")
    return expected


def pilot_inputs(api, repo):
    admission, admission_pin = committed(api, repo, ADMISSION)
    api.check(admission["readback"])
    readback = json.loads(Path(admission["readback"]["path"]).read_text())
    if (admission["status"] != "profile_weighting_orthobench_independently_admitted"
            or admission["scientific_revision"] != REVISION or readback["source_revision"] != REVISION
            or readback["status"] != "profile_pilot_independently_read_back" or readback["accuracy_result_admitted"] is not True
            or [cell["label"] for cell in readback["cells"]] != ["none", "henikoff", "henikoff_cluster", "species_henikoff_cluster"]):
        raise ValueError("Require the complete separately admitted weighting pilot")
    cells = readback["cells"]
    for cell in cells:
        if (cell["status"] != "independent_content_readback_verified" or cell["source_revision"] != REVISION
                or (cell["coverage"]["genes"], cell["coverage"]["species"]) != (251378, 12)):
            raise ValueError("Wrong pilot source or complete gene universe")
    def normalized(cell):
        return {"f1": cell["scores"]["OrthoBench: OrthoBench F1"], "precision": cell["details"]["precision"],
                "recall": cell["details"]["recall"]}
    control, advancing = normalized(cells[0]), []
    for compact, cell in zip(admission["arms"], cells):
        actual = normalized(cell)
        if compact["label"] != cell["label"] or any(compact[key] != value for key, value in actual.items()):
            raise ValueError("Compact admission and full pilot scores differ")
        eligible = False
        if cell["label"] != "none":
            gate = api.advancement(control, actual)
            if gate != cell["advancement"]:
                raise ValueError("Original pilot gate does not reproduce")
            eligible = gate["eligible"]
        if compact["advances"] != eligible:
            raise ValueError("Compact admission changes original eligibility")
        if eligible:
            advancing.append(cell["label"])
    if (len(admission["arms"]) != len(cells) or advancing != readback["advancing_labels"]
            or advancing != admission["advancing_labels"] or ("none", *advancing) != ARMS):
        raise ValueError("Follow-up must include every and only prespecified passing arm")
    accounting = [successful_job("24477", 16, "128G"), successful_job("24484", 1, "16G")]
    api.check(readback["execution_plan"])
    pilot = json.loads(Path(readback["execution_plan"]["path"]).read_text())
    if pilot["source_revision"] != REVISION:
        raise ValueError("Pilot scientific source differs")
    pins = [admission_pin, admission["readback"], admission["score_table"], readback["producer"],
            readback["execution_plan"], *readback["checked"]]
    return admission, readback, pilot, pins, accounting


def selected_case(api, repo, dataset):
    if dataset not in UNIVERSES:
        raise ValueError("Unsupported follow-up dataset")
    summary, context_pin = committed(api, repo, CONTEXT)
    control = api.verified_cell(summary, "canonical_collection:canonical", dataset)
    if (control["coverage"]["genes"], control["coverage"]["species"]) != UNIVERSES[dataset]:
        raise ValueError("Canonical control universe differs")
    root = repo / CANONICAL
    plan_path, report_path = root / "execution_plan.json", root / "results.json"
    canonical, report = json.loads(plan_path.read_text()), json.loads(report_path.read_text())
    api.check(report["execution_plan"])
    if report["execution_plan"] != api.record(plan_path):
        raise ValueError("Canonical report is not bound to its original plan")
    rows = [row for row in report["arms"] if (row["dataset"], row["label"]) == (dataset, "canonical")]
    configs = [case for case in canonical["cases"] if case["dataset"] == dataset]
    if len(rows) != 1 or len(configs) != 1:
        raise ValueError("Require a unique canonical follow-up case")
    row, config = rows[0], configs[0]
    native = root / dataset / "canonical/native"
    golden = {name: api.record(native / "orthohmm_phylogeny" / name) for name in (*api.PRODUCTS, "species_tree.rooted.nwk")}
    for name in api.PRODUCTS:
        if any(golden[name][key] != control["products"][name][key] for key in ("bytes", "sha256")):
            raise ValueError("Canonical native control differs from its admission")
    if golden["species_tree.rooted.nwk"] not in control["checked"]:
        raise ValueError("Canonical species tree is not the admitted control")
    manifest = row["ordering"]["manifest"]
    case = {"dataset": dataset, "case": {**config, "baseline_output": str(native)},
        "checkpoint": str(Path(manifest["path"]).parent), "checkpoint_sha": manifest["sha256"],
        "checkpoint_manifest": manifest, "control": control, "golden": golden, "candidate": row["candidate"],
        "original_producer_status": row["status"]}
    pins = [context_pin, summary["score_table"], api.record(plan_path), api.record(report_path), manifest,
        row["candidate"]["candidate_partition"], row["candidate"]["membership_constraints"],
        *config["fastas"], *config["scorers"], *golden.values(), *control["checked"]]
    if dataset == "QfO":
        recovery = repo / "benchmarks/results/canonical_qfo_scoring_recovery_r1_v1_20261014"
        pins += [api.record(recovery / name) for name in ("results.json", "execution_plan.json")]
    return case, pins


def scoring_id(output, mode):
    return "pw_" + hashlib.sha256(str(output / mode).encode()).hexdigest()[:12]


def prepare(repo, driver, output, dataset):
    if output.exists():
        raise FileExistsError("Preserve every follow-up execution namespace")
    if driver == repo or Path(__file__).resolve() != driver / "benchmark_tools/run_profile_weighting_followup.py":
        raise ValueError("Require the frozen standalone follow-up driver")
    science = repo / "benchmarks/work/profile_weighting_pilot_v1_20261018"
    api = science_api(science)
    admission, readback, pilot, pins, accounting = pilot_inputs(api, repo)
    protocol_path = driver / PROTOCOL
    protocol = json.loads(protocol_path.read_text())
    if (protocol["scientific_revision"] != REVISION or protocol["arms"] != list(ARMS)
            or protocol["datasets"] != list(UNIVERSES) or protocol["resources"]["cpus_per_run"] != 16):
        raise ValueError("Changed fixed weighting follow-up protocol")
    if (Path(pilot["source_root"]) != science or api.runtime() != pilot["runtime"]
            or subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=science, text=True).strip() != REVISION):
        raise ValueError("Original pilot source or scientific runtime changed")
    if api.sources(science) + [api.record(science / "benchmark_tools/run_profile_weighting_pilot.slurm")] != pilot["sources"]:
        raise ValueError("Original pilot scientific source inventory changed")
    case, case_pins = selected_case(api, repo, dataset)
    revision = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=driver, text=True).strip()
    driver_sources = [api.record(driver / name) for name in ("benchmark_tools/run_profile_weighting_followup.py",
        "benchmark_tools/run_profile_weighting_followup.slurm", PROTOCOL)]
    for pin in driver_sources:
        blob = subprocess.check_output(["git", "show", revision + ":" + str(Path(pin["path"]).relative_to(driver))], cwd=driver)
        if hashlib.sha256(blob).hexdigest() != pin["sha256"]:
            raise ValueError("Execution driver differs from committed source")
    pins += [*case_pins, *pilot["sources"], *pilot["runtime"]["files"]]
    pins += [api.record(path) for path in sorted(Path(case["checkpoint"]).iterdir()) if path.is_file()]
    phylo = Path(case["case"]["baseline_output"]) / "orthohmm_phylogeny"
    for name in ("alignments", "gene_trees", "checkpoints", "species_tree_inference"):
        pins += [api.record(path) for path in sorted((phylo / name).rglob("*")) if path.is_file()]
    for pin in [*pins, *driver_sources]:
        api.check(pin)
    pins = list({pin["path"]: pin for pin in pins}.values())
    identifiers = {mode: scoring_id(output, mode) for mode in ARMS[1:]}
    if len(set(identifiers.values())) != len(identifiers):
        raise ValueError("Scoring namespace collision")
    if dataset == "QfO":
        for identifier in identifiers.values():
            work = repo / "qfo_benchmark/w" / ("hc1_" + identifier)
            api.require_darwin_work_path(work)
            if work.exists():
                raise FileExistsError("QfO scoring namespace already exists")
    plan = {"status": "prepared_weighting_followup_unrun", "repo": str(repo), "output": str(output),
        "driver_root": str(driver), "driver_revision": revision, "driver_sources": driver_sources,
        "scientific_source": str(science), "scientific_revision": REVISION, "pilot_admission": api.record(repo / ADMISSION),
        "pilot_readback": admission["readback"], "pilot_plan": readback["execution_plan"], "prerequisite_accounting": accounting,
        "inputs": pins, "scientific_sources": pilot["sources"], "runtime": pilot["runtime"], "tool_path": pilot["tool_path"],
        "profile_runtime": api.require_profile_runtime(science), "cpu": 16, "dataset": dataset, "case": case,
        "arms": list(ARMS), "scoring_ids": identifiers, "defaults_promoted": False, "goal_completed": False,
        "scope": "Unchanged weighted-profile transfer from authenticated canonical initial hits, with exact-compatible native caches and inferred phylogeny. Shared-host incremental wall/CPU/max-process-RSS, not fresh initial search or isolated timing. Exact QfO control assessment may be reused only after native equivalence."}
    validate_plan(plan)
    output.mkdir(parents=True)
    api.write_json(output / "execution_plan.json", plan)
    return plan


def validate_plan(plan):
    if (plan["status"] != "prepared_weighting_followup_unrun" or plan["scientific_revision"] != REVISION
            or plan["cpu"] != 16 or tuple(plan["arms"]) != ARMS or plan["dataset"] not in UNIVERSES
            or set(plan["scoring_ids"]) != set(ARMS[1:])
            or any(plan["scoring_ids"][mode] != scoring_id(Path(plan["output"]), mode) for mode in ARMS[1:])):
        raise ValueError("Unexpected follow-up source, arms, budget, dataset or scoring IDs")
    case = plan["case"]
    if case["dataset"] != plan["dataset"] or (case["control"]["coverage"]["genes"], case["control"]["coverage"]["species"]) != UNIVERSES[plan["dataset"]]:
        raise ValueError("Wrong complete canonical follow-up universe")
    needed = [plan["pilot_admission"], plan["pilot_readback"], plan["pilot_plan"], case["checkpoint_manifest"],
        *case["case"]["fastas"], *case["golden"].values(), case["candidate"]["candidate_partition"], case["candidate"]["membership_constraints"]]
    if any(pin not in plan["inputs"] for pin in needed):
        raise ValueError("Follow-up case is not bound to authenticated inputs")
    if (Path(case["checkpoint_manifest"]["path"]) != Path(case["checkpoint"]) / "manifest.json"
            or case["checkpoint_manifest"]["sha256"] != case["checkpoint_sha"]
            or len(case["case"]["fastas"]) != UNIVERSES[plan["dataset"]][1]
            or any(Path(pin["path"]).parent != Path(case["case"]["fasta_directory"]) for pin in case["case"]["fastas"])
            or any(Path(pin["path"]).parent != Path(case["case"]["baseline_output"]) / "orthohmm_phylogeny"
                   for pin in case["golden"].values()) or len(set(plan["scoring_ids"].values())) != len(ARMS) - 1):
        raise ValueError("Follow-up paths, checkpoint or scoring namespaces differ")


def run(path):
    plan = json.loads(path.read_text())
    validate_plan(plan)
    science, driver, output = (Path(plan[key]) for key in ("scientific_source", "driver_root", "output"))
    if Path(__file__).resolve() != driver / "benchmark_tools/run_profile_weighting_followup.py" or (output / "results.json").exists():
        raise ValueError("Require frozen driver and preserve all prior attempts")
    api = science_api(science)
    pins = [api.record(path), *plan["inputs"], *plan["driver_sources"]]
    for pin in pins:
        api.check(pin)
    if (api.runtime() != plan["runtime"] or os.environ.get("CUDA_VISIBLE_DEVICES") != ""
            or subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=science, text=True).strip() != REVISION):
        raise ValueError("Original scientific source/runtime or CUDA visibility changed")
    admission, _, pilot, _, accounting = pilot_inputs(api, Path(plan["repo"]))
    if (accounting != plan["prerequisite_accounting"] or pilot["sources"] != plan["scientific_sources"]
            or admission["readback"] != plan["pilot_readback"] or plan["runtime"] != pilot["runtime"]
            or plan["tool_path"] != pilot["tool_path"]
            or api.sources(science) + [api.record(science / "benchmark_tools/run_profile_weighting_pilot.slurm")] != pilot["sources"]):
        raise ValueError("Admitted fixed pilot prerequisites changed")
    case = plan["case"]
    expected_case, _ = selected_case(api, Path(plan["repo"]), plan["dataset"])
    if expected_case != case:
        raise ValueError("Planned dataset differs from the admitted canonical case")
    names, species, q, t, scores = api.load_accuracy_checkpoint(Path(case["checkpoint"]), verify=True)
    owners = {}
    for pin in case["case"]["fastas"]:
        for sequence in api.SeqIO.parse(pin["path"], "fasta"):
            if sequence.id in owners or not sequence.seq:
                raise ValueError("Duplicate or empty complete input gene")
            owners[sequence.id] = Path(pin["path"]).name
    ownership = api.validate_ownership(names, species, owners)
    if (len(names), len(ownership["classes"])) != UNIVERSES[plan["dataset"]]:
        raise ValueError("Canonical initial-hit universe differs")
    report = {"status": "running", "execution_plan": api.record(path), "job_id": os.environ.get("SLURM_JOB_ID"),
        "dataset": plan["dataset"], "source_revision": REVISION, "driver_revision": plan["driver_revision"],
        "arms": [], "ownership": ownership, "host_before": api.host_activity(), "scope": plan["scope"],
        "defaults_promoted": False, "goal_completed": False}
    api.write_json(output / "results.json", report)
    env = {**os.environ, "PATH": plan["tool_path"], "PYTHONPATH": str(science), "PYTHONHASHSEED": "0",
        "CUDA_VISIBLE_DEVICES": "", "OMP_NUM_THREADS": "1", "NUMBA_NUM_THREADS": "16", "OPENBLAS_NUM_THREADS": "1", "MKL_NUM_THREADS": "1"}
    control_ok = False
    for mode in ARMS:
        row = {"dataset": plan["dataset"], "label": mode, "status": "validating_inputs", "scores": {}}
        report["arms"].append(row)
        api.write_json(output / "results.json", report)
        if mode != "none" and not control_ok:
            row["status"] = "not_run_failed_control_gate"
            api.write_json(output / "results.json", report)
            continue
        folder = output / mode
        folder.mkdir(exist_ok=False)
        started = time.monotonic()
        try:
            for pin in pins:
                api.check(pin)
            graph, command = api.commands(science, plan, case["case"], folder, Path(case["checkpoint"]), case["checkpoint_sha"])
            if mode != "none":
                graph.extend(["--profile-sequence-weighting", mode])
            api.write_json(folder / "preflight.json", {"graph": graph, "native": command})
            row["status"] = "graph_running"
            api.write_json(output / "results.json", report)
            row["graph_execution"] = api.timed_stage(graph, folder, "graph", science, env)
            if row["graph_execution"]["returncode"]:
                raise RuntimeError("Weighted-profile graph replay failed")
            metrics = folder / "graph_metrics.json"
            graph_report = json.loads(metrics.read_text())
            row["activation"] = api.validate_activation(graph_report, mode)
            row["graph_report"] = api.record(metrics)
            for stage in graph_report["stages"]:
                api.check(stage["output"])
                api.read_partition(Path(stage["output"]["path"]), set(names))
            before, candidate_started = resource.getrusage(resource.RUSAGE_SELF), time.monotonic()
            row["candidate"] = api.prepare_partition(folder / "graph/orthogroups_profiles_refined.txt", folder / "candidate",
                names, species, (q, t, scores), True, api._expand_phylogeny_candidates, api.load_membership_constraints)
            after = resource.getrusage(resource.RUSAGE_SELF)
            row["candidate_resources"] = {"wall_s": time.monotonic() - candidate_started,
                "cpu_s": after.ru_utime + after.ru_stime - before.ru_utime - before.ru_stime,
                "peak_parent_process_rss_GiB": after.ru_maxrss / 1024 ** 2,
                "memory_scope": "Parent lifetime maximum, not simultaneous process-tree memory"}
            if mode == "none":
                for key in ("candidate_partition", "membership_constraints"):
                    if any(row["candidate"][key][field] != case["candidate"][key][field] for field in ("bytes", "sha256")):
                        raise ValueError("Unchanged canonical candidate control differs")
            row["status"] = "native_running"
            api.write_json(output / "results.json", report)
            row["native_execution"] = api.timed_stage(command, folder, "native", science, env)
            if row["native_execution"]["returncode"]:
                raise RuntimeError("Weighted-profile native replay failed")
            phylo = folder / "native/orthohmm_phylogeny"
            prediction = folder / "native/orthohmm_working_res/orthohmm_edges_clustered.txt"
            row["coverage"] = api.validate_partition(Path(row["candidate"]["candidate_partition"]["path"]), prediction)
            row["outputs"] = [api.record(phylo / name) for name in api.PRODUCTS]
            native = api.replay_native(row, {"command": command})
            if (native["genes"], native["species"]) != UNIVERSES[plan["dataset"]]:
                raise ValueError("Native complete universe differs")
            row["status"] = "native_verified_scoring_pending"
            api.write_json(output / "results.json", report)
            if mode == "none":
                row["golden_matches"] = api.compare_golden(phylo, case["golden"])
                actual, details, score_pins = api.control_scores(api, plan, row, native, folder)
            else:
                if plan["dataset"] == "QfO":
                    identifier = plan["scoring_ids"][mode]
                    row.update(api.score_qfo(Path(plan["repo"]), phylo, folder, identifier))
                    participant = "ohmm_hillclimb_r1_" + identifier
                else:
                    row.update(api.score_three_kingdoms(Path(plan["repo"]), prediction, folder))
                    participant = None
                actual, details, score_pins = api.readback_scores(Path(plan["repo"]), plan["dataset"], row, native, folder, participant)
            row["checked"] = [api.record(metrics), *native["checked"], *score_pins]
            for pin in [*pins, *row["checked"]]:
                api.check(pin)
            if api.runtime() != plan["runtime"]:
                raise ValueError("Scientific runtime changed during arm")
            row.update(status="native_and_score_independently_verified", details=details, native_summary=native["summary"],
                       metric_vector=api.endpoint_vector(plan["dataset"], actual))
            if mode == "none":
                control_ok = True
        except Exception as error:
            row.update(status="failed", error_type=type(error).__name__, error=str(error))
        finally:
            row["total_incremental_wall_s"] = time.monotonic() - started
            api.write_json(output / "results.json", report)
    report["status"] = ("all_weighting_followup_arms_scored" if all(row["status"] == "native_and_score_independently_verified"
        for row in report["arms"]) else "finished_with_explicit_gaps")
    report["deltas"] = {row["label"]: {key: value - report["arms"][0]["metric_vector"][key]
        for key, value in row["metric_vector"].items()} for row in report["arms"][1:]
        if row["status"] == "native_and_score_independently_verified"}
    report["host_after"] = api.host_activity()
    api.write_json(output / "results.json", report)
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    subcommands = parser.add_subparsers(dest="command", required=True)
    prep = subcommands.add_parser("prepare")
    for name in ("repo", "driver-root", "output"):
        prep.add_argument("--" + name, required=True, type=Path)
    prep.add_argument("--dataset", choices=UNIVERSES, required=True)
    worker = subcommands.add_parser("run")
    worker.add_argument("--execution-plan", required=True, type=Path)
    args = parser.parse_args()
    result = (prepare(args.repo.resolve(), args.driver_root.resolve(), args.output.resolve(), args.dataset)
              if args.command == "prepare" else run(args.execution_plan.resolve()))
    print(result["status"])
    raise SystemExit(0 if result["status"] in ("prepared_weighting_followup_unrun", "all_weighting_followup_arms_scored") else 1)
