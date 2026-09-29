"""One explicit-GC stage diagnostic; never admits high-CPM scientific results."""

import argparse
import gc
import importlib.util
import json
import os
from pathlib import Path
import resource
import subprocess
import sys
import time

spec = importlib.util.spec_from_file_location("checkpoint_setup", Path(__file__).with_name("probe_cpm_checkpoint_setup.py"))
setup = importlib.util.module_from_spec(spec)
spec.loader.exec_module(setup)
helper = setup.helper
PRIOR = "benchmark_tools/results/qfo_cpm_checkpoint_setup_20260926.json"
PRIOR_SHA = "794bca9664f0e862ea6bcaac2ff4131af4b54cb77afab3e871c032ffa8109701"
PROTOCOL = "benchmark_tools/results/QFO_CPM_REFINEMENT_GC_PROTOCOL_20260928.md"


def collect(stage):
    setup.observe("before_gc_" + stage)
    unreachable = gc.collect(2)
    setup.observe("after_gc_" + stage, unreachable=unreachable)


def stages(replay, reader, np, checkpoint, native, output, collect_fn=collect):
    names, species, queries, targets, scores, numeric = replay.load_replay_input(
        checkpoint=Path(checkpoint["path"]).parent, checkpoint_sha256=checkpoint["sha256"])
    universe = set(names)
    if len(universe) != len(names) or (native / "payload/gene_names.txt").read_text().splitlines() != list(names):
        raise ValueError("Changed graph gene universe")
    reader(native / "orthogroups_profiles.txt", universe)
    clusters = replay.read_index_clusters(str(native / "orthogroups_profiles.txt"), names)
    arrays = [np.load(native / "payload" / name, mmap_mode="r", allow_pickle=False)
              for name in ("sources.npy", "targets.npy", "weights.npy")]
    hits = replay.production_refinement_hits(queries, targets, scores, species)
    collect_fn("before_refinement")
    setup.observe("before_refinement")
    refined = replay.refine_cluster_indices(clusters, *hits, arrays[0], arrays[1], species,
                                           rbnh_scores=arrays[2])
    setup.observe("after_refinement", groups=len(refined))
    collect_fn("after_refinement")
    destination = output / "refined.txt"
    if destination.exists() or destination.is_symlink():
        raise FileExistsError(destination)
    setup.observe("before_write")
    replay.write_clusters(destination, refined, names)
    setup.observe("after_write")
    collect_fn("after_write")
    setup.observe("before_readback")
    groups = reader(destination, universe)
    setup.observe("after_readback", groups=len(groups))
    return dict(genes=len(names), species=len(np.unique(species)), groups=len(groups),
                refinement_directed_hits=len(hits[2]), numeric_checkpoint=numeric,
                output=helper.record(destination))


def worker(root, output):
    resource.setrlimit(resource.RLIMIT_AS, (64 * 1024**3, 64 * 1024**3))
    resource.setrlimit(resource.RLIMIT_CPU, (300, 300))
    resource.setrlimit(resource.RLIMIT_CORE, (0, 0))
    os.sched_setaffinity(0, {min(os.sched_getaffinity(0))})
    if not gc.isenabled() or gc.get_threshold() != (700, 10, 10):
        raise ValueError("Require unchanged enabled GC defaults")
    launcher = root / "benchmarks/work/publication_qfo_replay_native_v1"
    modules = setup.context.scientific_imports(launcher)
    from benchmark_tools import replay_high_sensitivity as replay
    from benchmark_tools.audit_historical_profile_ablation import read_partition
    import numpy as np
    checkpoint = json.loads((root / setup.PLAN).read_text())["checkpoint_manifest"]
    result = stages(replay, read_partition, np, checkpoint, root / setup.NATIVE, output)
    helper.check(modules)
    result.update(scientific_sources=modules, affinity=sorted(os.sched_getaffinity(0)),
                  gc_enabled=gc.isenabled(), gc_thresholds=list(gc.get_threshold()),
                  accuracy_admitted=False)
    print(json.dumps(result, sort_keys=True))


def run(root, output, protocol_sha):
    if output.exists() or output.is_symlink():
        raise FileExistsError(output)
    prior_ref, protocol_ref = helper.record(root / PRIOR), helper.record(root / PROTOCOL)
    if prior_ref["sha256"] != PRIOR_SHA or protocol_ref["sha256"] != protocol_sha:
        raise ValueError("Prior setup or protocol changed")
    prior = json.loads((root / PRIOR).read_text())
    records = [prior_ref, protocol_ref, helper.record(__file__), *prior["checked_records"],
               *prior["result"]["scientific_sources"]]
    expected = helper.record(root / setup.NATIVE / "orthogroups_profiles_refined.txt")
    if expected["sha256"] != "f4c6f1973bc9636828baf1e6d9be3f416a18fc8ada5302fb502c082495fc1811":
        raise ValueError("Reference partition changed")
    records.append(expected)
    helper.check(records)
    env = {k: v for k, v in os.environ.items() if k not in ("PYTHONHOME", "LD_PRELOAD", "LD_LIBRARY_PATH")}
    env.update(prior["environment_overrides"])
    command = [prior["command"][0], "-B", str(Path(__file__).resolve()), "--root", str(root),
               "--output", str(output), "--worker"]
    output.mkdir(parents=True)
    report = dict(status="started", attempts=1, command=command, cwd=prior["cwd"],
        checked_records=records, environment_overrides=prior["environment_overrides"],
        accuracy_admitted=False, publication_ready=False,
        limits=dict(cpu_seconds=300, wall_seconds=360, address_space_bytes=64 * 1024**3, affinity_cpus=1),
        limitations=["Forced collections and stage logging alter allocation history.",
                     "One diagnostic only; no optimizer, default change, retry or scientific admission.",
                     "Shared-host elapsed time is descriptive, not comparative timing."])
    started = time.monotonic()
    try:
        try:
            done = subprocess.run(command, cwd=prior["cwd"], env=env, capture_output=True, timeout=360)
            stdout, stderr = done.stdout, done.stderr
            report.update(returncode=done.returncode, status="completed" if done.returncode == 0 else "failed")
        except subprocess.TimeoutExpired as error:
            stdout, stderr = error.stdout or b"", error.stderr or b""
            report.update(returncode=None, status="timed_out")
        for label, content in (("stdout", stdout), ("stderr", stderr)):
            path = output / label
            path.write_bytes(content)
            report[label] = helper.record(path)
        if report["status"] == "completed":
            result = json.loads(stdout)
            if ((result["genes"], result["species"], result["groups"]) != setup.EXPECTED_COUNTS
                    or result["refinement_directed_hits"] != 0
                    or result["accuracy_admitted"] is not False
                    or len(result["affinity"]) != 1
                    or not result["gc_enabled"] or result["gc_thresholds"] != [700, 10, 10]):
                raise ValueError("Unexpected diagnostic completion")
            actual = helper.record(output / "refined.txt")
            if result["output"] != actual or any(actual[k] != expected[k] for k in ("bytes", "sha256")):
                raise ValueError("Diagnostic partition bytes differ")
            helper.check(result["scientific_sources"])
            report["result"] = result
        helper.check(records)
    except Exception as error:
        report.update(status="validation_failed", error_type=type(error).__name__, error=str(error))
        raise
    finally:
        report["wall_seconds_descriptive_only"] = time.monotonic() - started
        with (output / "report.json").open("x") as handle:
            json.dump(report, handle, indent=2, sort_keys=True)
            handle.write("\n")
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--protocol-sha256")
    parser.add_argument("--worker", action="store_true")
    args = parser.parse_args()
    if args.worker:
        worker(args.root.resolve(), args.output.absolute())
    else:
        if not args.protocol_sha256:
            parser.error("--protocol-sha256 is required")
        result = run(args.root.resolve(), args.output.absolute(), args.protocol_sha256)
        raise SystemExit(0 if result["status"] == "completed" else 1)
