"""Reconstruct original high-CPM multipass refinement, without clustering."""

import argparse
import json
import os
from pathlib import Path
import subprocess
import sys

PREFIX_SHA = "f038a84f8dc2429f276a647962435a3b35f4f8b10d5a5a3ed2662c6fb7ef63a4"


def reconstruct(replay, read_partition, names, species, queries, targets, scores,
                sources, destinations, weights, partition, expected, output):
    if output.exists() or output.is_symlink():
        raise FileExistsError(output)
    universe = set(names)
    if len(universe) != len(names):
        raise ValueError("Duplicate gene universe")
    read_partition(partition, universe)
    reference = {frozenset(group) for group in read_partition(expected, universe)}
    clusters = replay.read_index_clusters(str(partition), names)
    refinement_hits = replay.production_refinement_hits(queries, targets, scores, species)
    refined = replay.refine_cluster_indices(clusters, *refinement_hits, sources, destinations,
                                            species, rbnh_scores=weights)
    replay.write_clusters(output, refined, names)
    observed = {frozenset(group) for group in read_partition(output, universe)}
    if observed != reference:
        raise ValueError("Reconstructed refinement memberships differ")
    return {"genes": len(names), "groups": len(observed), "partition_equal": True,
            "refinement_directed_hits": len(refinement_hits[2])}


def worker(root, output):
    launcher = root / "benchmarks/work/publication_qfo_replay_native_v1"
    if Path.cwd() != launcher:
        raise ValueError("Wrong frozen reconstruction working directory")
    sys.path.insert(0, str(launcher))
    from benchmark_tools import replay_high_sensitivity as replay
    from benchmark_tools.audit_historical_profile_ablation import read_partition
    from benchmark_tools.orthobench_stage_diagnostics import file_provenance
    import orthohmm.accuracy
    import orthohmm.refinement
    import numpy as np

    modules = [replay, orthohmm.accuracy, orthohmm.refinement,
               sys.modules[read_partition.__module__], sys.modules[replay.audit_numeric_checkpoint.__module__]]
    if any(not Path(module.__file__).resolve().is_relative_to(launcher) for module in modules):
        raise ValueError("Nonfrozen scientific reconstruction import")
    plan = json.loads((root / "benchmarks/work/qfo_corrected_replay_commands_20260918.json").read_text())
    checkpoint = plan["checkpoint_manifest"]
    names, species, queries, targets, scores, numeric = replay.load_replay_input(
        checkpoint=Path(checkpoint["path"]).parent, checkpoint_sha256=checkpoint["sha256"])
    if len(names) != 984137 or len(np.unique(species)) != 78:
        raise ValueError("Wrong reconstruction universe")
    directory = root / "benchmarks/results/qfo_parameter_cpm_replay_v3/cpm_high"
    payload = directory / "clustering/cluster_1_multipass/payload"
    if (payload / "gene_names.txt").read_text().splitlines() != list(names):
        raise ValueError("Graph and numeric gene order differ")
    sources, destinations, weights = [np.load(payload / name, mmap_mode="r", allow_pickle=False)
                                     for name in ("sources.npy", "targets.npy", "weights.npy")]
    result = reconstruct(replay, read_partition, names, species, queries, targets, scores,
        sources, destinations, weights, directory / "replay/orthogroups_multipass.txt",
        directory / "replay/orthogroups_multipass_refined.txt", output / "reconstructed.txt")
    result.update(numeric_checkpoint=numeric, modules=[file_provenance(Path(m.__file__)) for m in modules],
                  output=file_provenance(output / "reconstructed.txt"), accuracy_evaluated=False)
    with (output / "worker.json").open("x") as stream:
        json.dump(result, stream, indent=2, sort_keys=True)
        stream.write("\n")


def run(root, output):
    sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
    from benchmark_tools.prepare_ob_candidate_neighborhood import check, record
    from benchmark_tools.run_simulation_methods import read_frozen
    from benchmark_tools.run_blast_recovery_batch import save_status
    from benchmark_tools.verify_qfo_replay_launcher import verify
    from benchmark_tools.checked_replay_payload_worker import corrected_evidence
    from benchmark_tools.cpm_replay_context import REPLAY_SHA

    if output.exists() or output.is_symlink():
        raise FileExistsError(output)
    prefix_path = root / "benchmark_tools/results/qfo_cpm_recovery_predecessors_20260923.json"
    prefix = read_frozen(prefix_path, PREFIX_SHA)
    if (prefix["status"] != "cpm_recovery_predecessors_verified_not_admitted"
            or prefix["recovery_authorized"] is not False or prefix["accuracy_evaluated"] is not False):
        raise ValueError("Wrong predecessor evidence")
    plan, plan_record, admission_record, names_record = corrected_evidence(
        root / "benchmarks/work/qfo_corrected_replay_commands_20260918.json", REPLAY_SHA)
    launcher = root / "benchmarks/work/publication_qfo_replay_native_v1"
    core = root / "benchmarks/work/publication_method_native_v2"
    runtime_path = root / "benchmark_tools/results/publication_native_runtime_20260916.json"
    runtime = verify(core, launcher, runtime_path)
    if runtime != plan["runtime"]:
        raise ValueError("Changed reconstruction runtime")
    directory = root / "benchmarks/results/qfo_parameter_cpm_replay_v3/cpm_high"
    expected = record(directory / "replay/orthogroups_multipass_refined.txt")
    records = [record(__file__), record(prefix_path), plan_record, admission_record, names_record,
               expected, *prefix["checked_records"]]
    for item in records:
        check(item)
    output.mkdir()
    env = os.environ.copy()
    env.update(PYTHONPATH=str(launcher), PYTHONHASHSEED="0", OMP_NUM_THREADS="1",
               OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1", PYTHONNOUSERSITE="1")
    for key in ("PYTHONHOME", "LD_PRELOAD", "LD_LIBRARY_PATH"):
        env.pop(key, None)
    command = [sys.executable, "-B", str(Path(__file__).resolve()), "--root", str(root),
               "--output", str(output), "--worker"]
    report = dict(status="reconstructing_original_refinement", checked_records=records,
                  source=records[0], command=command, runtime_before=runtime, expected=expected,
                  job_id=os.environ.get("SLURM_JOB_ID"), recovery_authorized=False,
                  accuracy_evaluated=False, publication_ready=False)
    save_status(output / "status.json", report)
    try:
        with (output / "worker.log").open("x") as log:
            done = subprocess.run(command, cwd=launcher, env=env, stdout=log, stderr=subprocess.STDOUT)
        report["returncode"] = done.returncode
        if done.returncode:
            raise RuntimeError(f"Refinement reconstruction exited {done.returncode}")
        result = json.loads((output / "worker.json").read_text())
        if (result["partition_equal"] is not True or result["genes"] != 984137
                or result["accuracy_evaluated"] is not False or result["refinement_directed_hits"] != 0
                or result["output"] != record(output / "reconstructed.txt")):
            raise ValueError("Invalid refinement reconstruction result")
        report["runtime_after"] = verify(core, launcher, runtime_path)
        if report["runtime_after"] != runtime:
            raise ValueError("Runtime changed during reconstruction")
        for item in [*records, *result["modules"], result["output"]]:
            check(item)
        report.update(status="original_cpm_refinement_reproduced_unscored", result=result,
                      worker_report=record(output / "worker.json"))
    except BaseException as error:
        report.update(status="refinement_reconstruction_failed", error_type=type(error).__name__, error=str(error))
        raise
    finally:
        report["worker_log"] = record(output / "worker.log")
        save_status(output / "status.json", report)
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--worker", action="store_true")
    args = parser.parse_args()
    (worker if args.worker else run)(args.root.resolve(), args.output.absolute())
