"""Read back the two fixed dependency replays without scoring or inference."""

import argparse
import json
from pathlib import Path
import shutil
import subprocess
import tempfile

import numpy as np

from benchmark_tools.prepare_ob_candidate_neighborhood import record, check
from benchmark_tools.probe_installed_ob_clustering import validate_graph, write_json
from benchmark_tools.compare_installed_ob_search import partition
from benchmark_tools.audit_installed_orthobench import compare_partitions
from benchmark_tools.audit_candidate_arm import audit as audit_candidates

PLAN_SHA = "950a472f1ca767851c69446412edc0b26ea5fc5646160c1d385f6318327b87ba"
STAGES = ("multipass", "multipass_refined", "strict_profiles", "strict_profiles_refined")
ENV = dict(PATH="/usr/bin:/bin", HOME="/home/bizon", LANG="C.UTF-8",
           OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1", PYTHONHASHSEED="0")


def validate_execution(execution, arm, plan):
    command = [arm["python"], "-I", plan["source"]["path"], "--driver", arm["directory"],
               "--repo", plan["repo"], "--checkpoint", plan["checkpoint"], "--inputs", plan["inputs"],
               "--cpu", str(plan["cpu"])]
    if (execution["returncode"] != 0 or execution["timed_out"] is not False
            or execution["attempts"] != 1 or execution["job_id"] != "22320"
            or execution["command"] != command or execution["env"] != ENV):
        raise ValueError("Execution does not match the single planned successful arm")


def compare_graphs(left, right, names):
    arrays = []
    for path in (left, right):
        with np.load(path, allow_pickle=False) as data:
            if set(data.files) != {"sources", "targets", "weights"}:
                raise ValueError("Unexpected graph fields")
            values = [data[k] for k in ("sources", "targets", "weights")]
        validate_graph(names, *values)
        arrays.append(values)
    a, b = arrays
    same_endpoints = all(np.array_equal(x, y) for x, y in zip(a[:2], b[:2]))
    return dict(left_edges=len(a[0]), right_edges=len(b[0]), endpoints_equal=same_endpoints,
                arrays_equal=same_endpoints and np.array_equal(a[2], b[2]))


def candidate_readback(directory, stage, seed_record, universe):
    working = directory / "replay/orthohmm_working_res"
    names = ("orthohmm_edges_clustered.txt", "phylogeny_candidate_superfamilies.txt",
             "phylogeny_candidate_seeds.tsv", "phylogeny_candidate_merges.json")
    originals = [record(working / name) for name in names]
    details = stage["candidates"]
    for key, name in zip(("candidate_checkpoint", "seed_sidecar", "merge_trace_sidecar"), names[1:]):
        if details[key] != str(working / name):
            raise ValueError("Candidate path differs")
    # Adapt only paths to the existing strict candidate-content auditor.
    with tempfile.TemporaryDirectory(prefix="ob-dependency-candidate-") as temporary:
        target = Path(temporary)
        dest = target / "orthohmm_working_res"
        dest.mkdir()
        for name in names:
            shutil.copyfile(working / name, dest / name)
        expansion = dict(details)
        for key, name in zip(("candidate_checkpoint", "seed_sidecar", "merge_trace_sidecar"), names[1:]):
            expansion[key] = str(dest / name)
        arm = dict(seed_partition=seed_record, candidate_expansion=True, expansion=expansion,
                   candidate_partition=record(dest / names[0]), membership_constraints=record(dest / names[3]),
                   output_files=[record(p) for p in sorted(target.rglob("*")) if p.is_file()])
        result = audit_candidates(arm, seed_record, target, universe, True)
    result.pop("checked_records")
    result["original_records"] = originals
    for item in originals:
        check(item)
    return result


def audit(repo, directory):
    records = []
    def read(path, sha=None):
        item = record(path)
        if sha is not None and item["sha256"] != sha:
            raise ValueError("Changed pinned evidence")
        records.append(item)
        return json.loads(Path(path).read_text())

    plan = read(directory / "plan.json", PLAN_SHA)
    terminal = read(directory / "execution.json")
    if terminal != dict(status="native_complete_pending_readback", plan=record(directory / "plan.json")):
        raise ValueError("Missing matching terminal receipt")
    command = ["sacct", "-j", "22320", "--format=JobID,State,ExitCode", "-n", "-P"]
    scheduler = subprocess.check_output(command, text=True, timeout=60)
    if [line for line in scheduler.splitlines() if line.startswith("22320|")] != ["22320|COMPLETED|0:0"]:
        raise ValueError("Scheduler does not confirm success")
    if [a["label"] for a in plan["arms"]] != ["leiden012", "leiden011"]:
        raise ValueError("Wrong arm inventory")
    records.extend(plan["checked_records"])
    for item in records:
        check(item)
    pinned = {r["path"]: r for r in records}
    clean = read(repo / "benchmark_tools/results/installed_ob_clustering_probe_20260926.json",
                 "2fd69ad44dfe5ff5d1bb6f4d82b31523d6cce06d2a2e8173ad036ef68f527f13")
    overlay = read(repo / "benchmark_tools/results/ob_leiden_overlay_probe_20260926.json",
                   "3cc5b1795d1d0b361a02e4d18fc98e0a1bc5307ab6c2c086fe28e7bc0e0c4e4a")
    native_pins = {r["path"]: r for r in clean["runtime"]["files"]}
    native_pins.update({r["copy"]["path"]: r["copy"] for r in overlay["copied"]})
    names_path = Path(plan["checkpoint"]) / "gene_names.txt"
    names = names_path.read_text().splitlines()
    if len(names) != 251378 or len(set(names)) != len(names):
        raise ValueError("Wrong checkpoint universe")
    universe = set(names)
    groups, graphs, arms = {}, {}, []
    for arm in plan["arms"]:
        label, root = arm["label"], Path(arm["directory"])
        if root != directory / label:
            raise ValueError("Arm directory differs")
        records.append(arm["pth"])
        execution = read(root / "execution.json")
        validate_execution(execution, arm, plan)
        records.append(execution["log"])
        runtime = read(root / "runtime.json")
        stage, replay = read(root / "stage_report.json"), read(root / "replay.json")
        expected_version = "0.12.0" if label == "leiden012" else "0.11.0"
        if (runtime != stage["runtime"] or runtime["executable"] != arm["python"]
                or runtime["versions"] != dict(numpy="2.2.6", igraph="1.0.0", leidenalg=expected_version)
                or stage["status"] != "prephylogeny_dependency_replay_complete"
                or stage["accuracy_evaluated"] is not False or stage["phylogeny_run"] is not False):
            raise ValueError("Runtime or diagnostic scope differs")
        if len(runtime["sources"]) != 29 or len(runtime["native"]) != 2:
            raise ValueError("Unexpected runtime inventory")
        for item in runtime["sources"]:
            if pinned.get(item["path"]) != item:
                raise ValueError("Scientific source differs from plan")
        for item in runtime["native"]:
            if native_pins.get(item["path"]) != item:
                raise ValueError("Native extension differs from previous distribution evidence")
        externals = next(r["path"] for r in runtime["sources"] if r["path"].endswith("/externals.py"))
        if runtime["plain_child_import_probe"] != [expected_version, externals, runtime["native"][1]["path"]]:
            raise ValueError("Child import differs")
        records.extend(runtime["sources"] + runtime["native"])
        if (replay["parameters"] != dict(accuracy_profile="high_sensitivity", cpm_resolution=0.1,
                jackknife_profile_thresholds=False, jackknife_single_copy_profiles=False, leiden_seed=4,
                matrix="BLOSUM62", profile_expansion=True, profile_iterations=1, profile_min_species=1)
                or replay["input"]["manifest"] != pinned[str(Path(plan["checkpoint"]) / "manifest.json")]
                or replay["source"] != pinned[str(repo / "benchmark_tools/replay_high_sensitivity.py")]):
            raise ValueError("Replay settings or source differ")
        if [s["index"] for s in stage["snapshots"]] != list(range(4)):
            raise ValueError("Expected exactly four clustering calls")
        groups[label], graphs[label] = {}, {}
        for snapshot in stage["snapshots"]:
            i = snapshot["index"]
            before = read(root / f"clustering_{i}/before.json")
            if before != dict(graph=snapshot["graph"], resolution=0.1, include_isolates=True, seed=4):
                raise ValueError("Clustering settings differ")
            for key, filename in (("graph", "graph.npz"), ("partition", "partition.txt")):
                if snapshot[key] != record(root / f"clustering_{i}" / filename):
                    raise ValueError("Snapshot path or content differs")
                records.append(snapshot[key])
            groups[label][f"clustering_{i}"] = partition(Path(snapshot["partition"]["path"]), universe)
            graphs[label][i] = snapshot["graph"]["path"]
        if tuple(s["label"] for s in replay["stages"]) != STAGES:
            raise ValueError("Unexpected retained refinement stages")
        for entry in replay["stages"]:
            records.append(entry["output"])
            value = partition(Path(entry["output"]["path"]), universe)
            if len(value) != entry["clusters"]:
                raise ValueError("Stage count differs")
            groups[label][entry["label"]] = value
        for i, name in ((1, "multipass"), (3, "strict_profiles")):
            if groups[label][f"clustering_{i}"] != groups[label][name]:
                raise ValueError("Replay stage does not match clustering snapshot")
        candidates = stage["candidate_partition"]
        if candidates != record(root / "replay/orthohmm_working_res/phylogeny_candidate_superfamilies.txt"):
            raise ValueError("Candidate checkpoint differs")
        records.append(candidates)
        groups[label]["candidates"] = partition(Path(candidates["path"]), universe)
        candidate_check = candidate_readback(root, stage, replay["stages"][-1]["output"], universe)
        records.extend(candidate_check["original_records"])
        arms.append(dict(label=label, groups={k: len(v) for k, v in groups[label].items()},
                         candidate_readback=candidate_check, wall_seconds_descriptive=execution["wall_seconds"]))
    comparisons = {k: compare_partitions(groups["leiden012"][k], groups["leiden011"][k])
                   for k in groups["leiden012"]}
    graph_comparisons = {str(i): compare_graphs(graphs["leiden012"][i], graphs["leiden011"][i], names)
                         for i in range(4)}
    if not graph_comparisons["0"]["arrays_equal"]:
        raise ValueError("Initial graphs differ; distribution contrast is confounded")
    previous = read(repo / "benchmark_tools/results/installed_ob_search_comparison_20260926.json",
                    "247133956257e25cfbf5b2b79a19544ea382de07f9a7393f275d707ee8a658c5")
    candidate_records = {r["sha256"]: r for r in previous["checked_records"]
                         if r["path"].endswith("/phylogeny_candidate_superfamilies.txt")}
    baseline_comparisons = {}
    for label, sha in (("historical", "44f201dbc7e4eccf06d9401e5ce1b20fdf6ad523e15ee6007fe1c58f3941842d"),
                       ("fresh_installed", "117519224fb5c720de25412ebc989a99fc8d119e77cbb221aedcc03530738c63")):
        item = candidate_records[sha]
        records.append(item)
        reference = partition(Path(item["path"]), universe)
        baseline_comparisons[label] = {arm: compare_partitions(reference, values["candidates"])
                                       for arm, values in groups.items()}
    for name in ("prepare_ob_candidate_neighborhood.py", "probe_installed_ob_clustering.py",
                 "compare_installed_ob_search.py", "audit_installed_orthobench.py",
                 "audit_candidate_arm.py", "audit_historical_profile_ablation.py"):
        records.append(record(repo / "benchmark_tools" / name))
    for item in records:
        check(item)
    return dict(status="dependency_replay_content_read_back", source=record(__file__), checked_records=records,
                scheduler=dict(command=command, output=scheduler), arms=arms, comparisons=comparisons,
                graph_comparisons=graph_comparisons, candidate_baseline_comparisons=baseline_comparisons,
                accuracy_evaluated=False, historical_scores_replaced=False,
                limitations=["Pre-phylogeny stage comparison, not attribution of final F1 changes.",
                             "One execution per distribution; no universal determinism claim.",
                             "Candidate merge consistency checked, not independently recomputed support.",
                             "Shared-host durations are descriptive, not comparative efficiency evidence."])


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", type=Path, required=True)
    parser.add_argument("--directory", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    write_json(args.output, audit(args.repo.resolve(), args.directory.resolve()))
