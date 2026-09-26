"""One checkpoint/setup diagnostic, stopped before high-CPM refinement."""

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

spec = importlib.util.spec_from_file_location("parser_import_control", Path(__file__).with_name("probe_cpm_parser_imports.py"))
context = importlib.util.module_from_spec(spec)
spec.loader.exec_module(context)
helper = context.helper
PRIOR = "benchmark_tools/results/qfo_cpm_parser_imports_20260926.json"
PRIOR_SHA = "80ca73e6a9557a00c51940c3778ae0765bf9f8db19891f09175de34f5f2788b9"
PLAN = "benchmarks/work/qfo_corrected_replay_commands_20260918.json"
PLAN_SHA = "e7657732bc94438fb0602f3246680308a075b9fda1c68909ddfac1138b164c0c"
PROTOCOL = "benchmark_tools/results/QFO_CPM_CHECKPOINT_SETUP_PROTOCOL_20260926.md"
NATIVE = "benchmarks/results/qfo_cpm_checkpoint_recovery_v1"
EXPECTED_COUNTS = (984137, 78, 390845)
SETUP_PINS = {
    "orthogroups_profiles.txt": "6ad67311d9caea560cec0ff429a1185e0babb1c9aaa883486b9980ee72ee260a",
    "payload/sources.npy": "b7ba1120342b9d2becb8281155619d3ae356a8197001decb198061e10c01ba5c",
    "payload/targets.npy": "707c679c0c9360282cb0c9a6b1d635c60c284302f5d9916f009d82e0e42f5e42",
    "payload/weights.npy": "c78363f7b2a856d8d947e7e2be47616230323839ee7ad43d16914bb33ade3d96",
}


def observe(stage, **values):
    print(json.dumps(dict(stage=stage, gc_enabled=gc.isenabled(), gc_stats=gc.get_stats(),
                          **values), sort_keys=True), file=sys.stderr, flush=True)


def setup(replay, read_partition, np, checkpoint, native, final):
    observe("before_checkpoint")
    names, species, queries, targets, scores, numeric = replay.load_replay_input(
        checkpoint=Path(checkpoint["path"]).parent, checkpoint_sha256=checkpoint["sha256"])
    observe("after_checkpoint", genes=len(names), hits=len(scores))
    if (len(names), len(np.unique(species))) != EXPECTED_COUNTS[:2]:
        raise ValueError("Wrong checkpoint universe")
    universe = set(names)
    observe("before_checkpoint_parser")
    first = read_partition(final, universe)
    observe("after_checkpoint_parser", groups=len(first))
    del first
    if (native / "payload/gene_names.txt").read_text().splitlines() != list(names):
        raise ValueError("Checkpoint and graph gene order differ")
    read_partition(native / "orthogroups_profiles.txt", universe)
    clusters = replay.read_index_clusters(str(native / "orthogroups_profiles.txt"), names)
    arrays = [np.load(native / "payload" / name, mmap_mode="r", allow_pickle=False)
              for name in ("sources.npy", "targets.npy", "weights.npy")]
    hits = replay.production_refinement_hits(queries, targets, scores, species)
    observe("after_setup", seed_groups=len(clusters), graph_edges=len(arrays[0]), refinement_hits=len(hits[2]))
    groups = read_partition(final, universe)
    observe("after_setup_parser", groups=len(groups))
    if len(groups) != EXPECTED_COUNTS[2] or len(hits[2]) != 0:
        raise ValueError("Unexpected setup/parser outcome")
    # Keep checkpoint, graph and cluster objects alive through the final readback.
    result = dict(genes=len(names), groups=len(groups), seed_groups=len(clusters), numeric_checkpoint=numeric,
        refinement_directed_hits=len(hits[2]), arrays=[dict(shape=list(a.shape), dtype=str(a.dtype),
            nbytes=int(a.nbytes), memory_mapped=isinstance(a, np.memmap)) for a in (*arrays, species, queries, targets, scores)])
    observe("stopped_before_refinement")
    return result


def worker(root):
    resource.setrlimit(resource.RLIMIT_AS, (8 * 1024**3, 8 * 1024**3))
    resource.setrlimit(resource.RLIMIT_CPU, (300, 300))
    resource.setrlimit(resource.RLIMIT_CORE, (0, 0))
    launcher = root / "benchmarks/work/publication_qfo_replay_native_v1"
    modules = context.scientific_imports(launcher)
    from benchmark_tools import replay_high_sensitivity as replay
    from benchmark_tools.audit_historical_profile_ablation import read_partition
    import numpy as np
    checkpoint = json.loads((root / PLAN).read_text())["checkpoint_manifest"]
    result = setup(replay, read_partition, np, checkpoint, root / NATIVE,
                   root / NATIVE / "orthogroups_profiles_refined.txt")
    helper.check(modules)
    result.update(scientific_sources=modules, accuracy_admitted=False, publication_ready=False)
    print(json.dumps(result, sort_keys=True))


def run(root, output):
    if output.exists() or output.is_symlink():
        raise FileExistsError(output)
    prior_record, plan_record = helper.record(root / PRIOR), helper.record(root / PLAN)
    if (prior_record["sha256"], plan_record["sha256"]) != (PRIOR_SHA, PLAN_SHA):
        raise ValueError("Changed import controls or checkpoint plan")
    prior = json.loads((root / PRIOR).read_text())
    plan = json.loads((root / PLAN).read_text())
    checkpoint = plan["checkpoint_manifest"]
    helper.check([checkpoint])
    manifest = json.loads(Path(checkpoint["path"]).read_text())
    inputs = [prior_record, plan_record, *prior["checked_records"], checkpoint,
              *prior["arms"][1]["result"]["scientific_sources"], helper.record(__file__), helper.record(root / PROTOCOL)]
    inputs.extend(dict(path=str(Path(checkpoint["path"]).parent / name), **item)
                  for name, item in manifest["files"].items())
    for name, sha in SETUP_PINS.items():
        item = helper.record(root / NATIVE / name)
        if item["sha256"] != sha:
            raise ValueError("Changed saved setup input")
        inputs.append(item)
    helper.check(inputs)
    original = json.loads((root / helper.STATUS).read_text())
    launcher = root / "benchmarks/work/publication_qfo_replay_native_v1"
    env = {k: v for k, v in os.environ.items() if k not in ("PYTHONHOME", "LD_PRELOAD", "LD_LIBRARY_PATH")}
    overrides = prior["arms"][1]["environment_overrides"]
    env.update(overrides)
    command = [original["scientific_child_command"][0], "-B", str(Path(__file__).resolve()),
               "--root", str(root), "--worker"]
    output.mkdir(parents=True)
    report = dict(status="checkpoint_setup_running", checked_records=inputs, command=command, cwd=str(launcher),
        environment_overrides=overrides, attempts=1, accuracy_admitted=False, publication_ready=False,
        limits=dict(cpu_seconds=300, address_space_bytes=8 * 1024**3, wall_seconds=300),
        limitations=["Additional parser probes alter allocation history; this is not full refinement replay.",
                     "No refinement, optimizer, score admission or comparative timing; one attempt only."])
    started = time.monotonic()
    try:
        try:
            done = subprocess.run(command, cwd=launcher, env=env, capture_output=True, timeout=300)
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
            report["result"] = json.loads(stdout)
        helper.check(inputs)
    except Exception as error:
        report.update(status="checkpoint_setup_failed", error_type=type(error).__name__, error=str(error))
        raise
    finally:
        report["wall_seconds_descriptive_only"] = time.monotonic() - started
        with (output / "report.json").open("x") as stream:
            json.dump(report, stream, indent=2, sort_keys=True)
            stream.write("\n")
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--worker", action="store_true")
    args = parser.parse_args()
    if args.worker:
        worker(args.root.resolve())
    else:
        if args.output is None:
            parser.error("--output is required")
        run(args.root.resolve(), args.output.absolute())
