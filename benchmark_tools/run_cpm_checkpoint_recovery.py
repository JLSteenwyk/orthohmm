"""One fixed-setting high-CPM continuation; never admits downstream scoring."""

import argparse
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time


def refinement_worker(root, output, repeat=False):
    launcher = root / "benchmarks/work/publication_qfo_replay_native_v1"
    if Path.cwd() != launcher:
        raise ValueError("Wrong frozen refinement working directory")
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
        raise ValueError("Nonfrozen recovery refinement import")
    plan = json.loads((root / "benchmarks/work/qfo_corrected_replay_commands_20260918.json").read_text())
    checkpoint = plan["checkpoint_manifest"]
    names, species, queries, targets, scores, numeric = replay.load_replay_input(
        checkpoint=Path(checkpoint["path"]).parent, checkpoint_sha256=checkpoint["sha256"])
    if len(names) != 984137 or len(np.unique(species)) != 78:
        raise ValueError("Wrong recovery refinement universe")
    payload = output / "payload"
    if (payload / "gene_names.txt").read_text().splitlines() != list(names):
        raise ValueError("Recovery graph and numeric gene order differ")
    partition = output / "orthogroups_profiles.txt"
    read_partition(partition, set(names))
    clusters = replay.read_index_clusters(str(partition), names)
    sources, targets_graph, weights = [np.load(payload / name, mmap_mode="r", allow_pickle=False)
                                     for name in ("sources.npy", "targets.npy", "weights.npy")]
    hits = replay.production_refinement_hits(queries, targets, scores, species)
    refined = replay.refine_cluster_indices(clusters, *hits, sources, targets_graph, species, rbnh_scores=weights)
    label = "refinement_repeat" if repeat else "refinement"
    destination = output / ("refinement_repeat.txt" if repeat else "orthogroups_profiles_refined.txt")
    if destination.exists() or destination.is_symlink():
        raise FileExistsError(destination)
    replay.write_clusters(destination, refined, names)
    groups = read_partition(destination, set(names))
    result = dict(genes=len(names), groups=len(groups), refinement_directed_hits=len(hits[2]),
        numeric_checkpoint=numeric, modules=[file_provenance(Path(m.__file__)) for m in modules],
        output=file_provenance(destination), accuracy_evaluated=False)
    with (output / f"{label}.json").open("x") as stream:
        json.dump(result, stream, indent=2, sort_keys=True)
        stream.write("\n")


def optimizer_evidence(root, output, preflight, environment):
    from benchmark_tools.prepare_ob_candidate_neighborhood import check, record
    from benchmark_tools.repeat_qfo_saved_graph import check_worker
    from benchmark_tools.probe_leiden_boundary import saved_fingerprint
    from benchmark_tools.audit_historical_profile_ablation import read_partition

    payload = output / "payload"
    paths = [payload / name for name in ("worker_before.json", "native_boundary.json",
             "constructor_adapter.json", "constructor_parity.json")]
    records = [record(path) for path in paths]
    snapshot, boundary, adapter, parity = [json.loads(path.read_text()) for path in paths]
    launcher = Path(preflight["context"]["cwd"])
    overrides = {key: environment[key] for key in (
        "PYTHONPATH", "PYTHONHASHSEED", "OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS")}
    check_worker(snapshot, launcher, payload, overrides, expected_cpm_resolution=.12)
    expected_inputs = [record(payload / name) for name in ("gene_names.txt", "sources.npy", "targets.npy", "weights.npy")]
    if (len(snapshot["cpu_affinity"]) != 1 or snapshot["inputs"] != expected_inputs
            or snapshot["observer"] != record(Path(__file__).with_name("repeat_qfo_saved_graph.py"))
            or snapshot["python"] != record(sys.executable)):
        raise ValueError("Recovery worker affinity, inputs or observer differs")
    saved = saved_fingerprint(payload)
    if saved != preflight["saved_graph"]:
        raise ValueError("Recovery saved graph changed")
    arguments = dict(initial_membership=None, weights="weight", n_iterations=2, max_comm_size=0,
        seed=4, kwargs={"resolution_parameter": .12}, partition_type="leidenalg.VertexPartition.CPMVertexPartition")
    if boundary != {"accuracy_evaluated": False, "calls": [dict(arguments=arguments, before=saved,
            saved=saved, after=saved, status="optimizer_returned")]}:
        raise ValueError("Missing or changed recovery optimizer observation")
    expected_constructor = dict(status="constructor_returned", shape=[saved["edges"], 2], dtype="int32",
        ordered_input_bytes_sha256="93b36aad4916b85195c4dcc7dccb37d944cb911f7b051cf989ad1479729fded4",
        n=saved["vertices"], directed=False)
    if adapter != dict(format="python_pairs", accuracy_evaluated=False, calls=[expected_constructor]):
        raise ValueError("Recovery constructor observation differs")
    if (parity["status"] != "constructor_parity_verified_before_optimizer"
            or parity["accuracy_evaluated"] is not False or parity["fingerprint"] != saved or parity["saved"] != saved
            or any(parity["differences"][name]["different_edges"] != 0 for name in (
                "native_vs_saved", "constructor_vs_saved", "native_vs_constructor"))):
        raise ValueError("Recovery exhaustive endpoint comparison failed")
    names = (payload / "gene_names.txt").read_text().splitlines()
    if len(names) != 984137 or len(set(names)) != len(names):
        raise ValueError("Wrong recovered gene universe")
    partition = output / "orthohmm_working_res/orthohmm_edges_clustered.txt"
    groups = read_partition(partition, set(names))
    records.extend([record(partition), *snapshot["modules"].values(), *snapshot["native_libraries"], snapshot["python"]])
    for item in records:
        check(item)
    return dict(partition=record(partition), genes=len(names), groups=len(groups), saved_graph=saved,
                checked_records=records, accuracy_evaluated=False)


def run(root, commit):
    sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
    from benchmark_tools.prepare_cpm_checkpoint_recovery import prepare
    from benchmark_tools.prepare_ob_candidate_neighborhood import check, record
    from benchmark_tools.run_blast_recovery_batch import save_status
    from benchmark_tools.validate_cpm_refinement_reconstruction import compare_partitions
    from benchmark_tools.verify_qfo_replay_launcher import verify

    if (os.environ.get("SLURM_CPUS_PER_TASK") != "1" or os.environ.get("SLURM_MEM_PER_NODE") != "65536"
            or os.environ.get("SLURMD_NODENAME") != "bizon" or not os.environ.get("SLURM_JOB_ID")):
        raise ValueError("Require one-CPU 64-GiB bizon allocation")
    executor = Path(__file__).resolve().parent.parent
    revision = subprocess.check_output(["git", "-C", str(executor), "rev-parse", "HEAD"], text=True).strip()
    if revision != commit:
        raise ValueError("Recovery executor revision differs")
    subprocess.run(["git", "-C", str(executor), "diff", "--exit-code", "HEAD", "--", "benchmark_tools", "orthohmm"], check=True)
    output = root / "benchmarks/results/qfo_cpm_checkpoint_recovery_v1"
    if output.exists() or output.is_symlink():
        raise FileExistsError(output)
    preflight_directory = root / "benchmarks/results/qfo_cpm_checkpoint_recovery_preflight_v1"
    preflight = prepare(root, preflight_directory)
    if preflight["status"] != "cpm_checkpoint_preflight_verified_unscored" or preflight["preflight_passed"] is not True:
        raise ValueError("Recovery preflight incomplete")
    records = [record(preflight_directory / "status.json"), *preflight["checked_records"],
               *[record(path) for path in sorted((executor / "benchmark_tools").glob("*.py"))]]
    reference_record = record(root / "benchmarks/work/qfo_cpm_refinement_check_22153/worker.json")
    if reference_record not in preflight["checked_records"]:
        raise ValueError("Missing admitted reconstruction reference")
    reference_refinement = json.loads(Path(reference_record["path"]).read_text())
    for item in records:
        check(item)
    output.mkdir()
    payload = output / "payload"
    payload.mkdir()
    for name in ("gene_names.txt", "sources.npy", "targets.npy", "weights.npy"):
        (payload / name).symlink_to(Path(preflight["saved_payload"]) / name)
    save_status(payload / "metadata.json", dict(cpm_resolution=.12, seed=4, include_isolates=True,
                                               output_directory=str(output)))
    records.append(record(payload / "metadata.json"))
    launcher = Path(preflight["context"]["cwd"])
    env = os.environ.copy()
    env.update(PYTHONPATH=str(launcher), PYTHONHASHSEED="0", PYTHONNOUSERSITE="1",
               OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1")
    for key in ("PYTHONHOME", "LD_PRELOAD", "LD_LIBRARY_PATH"):
        env.pop(key, None)
    report = dict(status="checkpoint_recovery_running", source=record(__file__), executor_commit=commit,
        job_id=os.environ["SLURM_JOB_ID"], preflight=records[0], checked_records=records, phases=[],
        optimizer_attempted=False, accuracy_evaluated=False, downstream_admitted=False, publication_ready=False,
        missing_original_statistics=["profile_counters", "successful_stage_timings", "full_run_time"])
    save_status(output / "status.json", report)

    def execute(mode):
        command = [sys.executable, "-B", str(Path(__file__).resolve()), "--root", str(root),
                   "--output", str(output), "--mode", mode]
        phase = dict(mode=mode, command=command, status="running")
        report["phases"].append(phase)
        if mode == "optimize":
            report["optimizer_attempted"] = True
        save_status(output / "status.json", report)
        started = time.monotonic()
        with (output / f"{mode}.log").open("x") as log:
            done = subprocess.run(command, cwd=launcher, env=env, stdout=log, stderr=subprocess.STDOUT)
        phase.update(returncode=done.returncode, wall_s=time.monotonic() - started,
                     log=record(output / f"{mode}.log"), status="completed" if done.returncode == 0 else "failed")
        save_status(output / "status.json", report)
        if done.returncode:
            raise RuntimeError(f"Recovery {mode} exited {done.returncode}; no retry")

    try:
        execute("optimize")
        report["optimizer"] = optimizer_evidence(root, output, preflight, env)
        shutil.copyfile(report["optimizer"]["partition"]["path"], output / "orthogroups_profiles.txt")
        execute("refine")
        execute("repeat-refinement")
        reports = [json.loads((output / name).read_text()) for name in ("refinement.json", "refinement_repeat.json")]
        for result, filename in zip(reports, ("orthogroups_profiles_refined.txt", "refinement_repeat.txt")):
            if result["genes"] != 984137 or result["refinement_directed_hits"] != 0 or result["accuracy_evaluated"] is not False:
                raise ValueError("Wrong refinement result")
            if (result["output"] != record(output / filename)
                    or result["numeric_checkpoint"] != reference_refinement["numeric_checkpoint"]
                    or result["modules"] != reference_refinement["modules"]):
                raise ValueError("Refinement source, checkpoint or output differs")
            for item in [result["output"], *result["modules"]]:
                check(item)
        if any(reports[0][key] != reports[1][key] for key in ("genes", "groups", "numeric_checkpoint", "modules")):
            raise ValueError("Independent refinement evidence differs")
        names = (payload / "gene_names.txt").read_text().splitlines()
        report["refinement_comparison"] = compare_partitions(output / "orthogroups_profiles_refined.txt",
            output / "refinement_repeat.txt", names, reports[0]["groups"])
        original = Path(preflight["context"]["output_root"]) / "replay"
        report["stages"] = [dict(label=label, origin=origin, output=record(path)) for label, origin, path in (
            ("multipass", "reused", original / "orthogroups_multipass.txt"),
            ("multipass_refined", "reused", original / "orthogroups_multipass_refined.txt"),
            ("strict_profiles", "recovered", output / "orthogroups_profiles.txt"),
            ("strict_profiles_refined", "recovered", output / "orthogroups_profiles_refined.txt"))]
        report["runtime_after"] = verify(root / "benchmarks/work/publication_method_native_v2", launcher,
                                        root / "benchmark_tools/results/publication_native_runtime_20260916.json")
        if report["runtime_after"] != preflight["runtime"]:
            raise ValueError("Runtime changed during recovery")
        for item in [*records, *report["optimizer"]["checked_records"], *[r["output"] for r in report["stages"]]]:
            check(item)
        report.update(status="cpm_checkpoint_recovered_pending_independent_admission",
            refinement_reports=[record(output / name) for name in ("refinement.json", "refinement_repeat.json")],
            limitations=["Original job remains failed; independent recovery admission and all downstream stages remain required.",
                         "Shared-host incremental phase timings are not controlled end-to-end timing."])
    except BaseException as error:
        report.update(status="checkpoint_recovery_failed", error_type=type(error).__name__, error=str(error))
        raise
    finally:
        save_status(output / "status.json", report)
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--commit")
    parser.add_argument("--mode", choices=("optimize", "refine", "repeat-refinement"))
    args = parser.parse_args()
    root = args.root.resolve()
    if args.mode:
        if args.output is None or args.commit is not None:
            parser.error("Child mode requires output and no commit")
        if args.mode == "optimize":
            sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
            from benchmark_tools.guard_cpm_optimizer import run_frozen_worker
            run_frozen_worker(root / "benchmarks/work/publication_qfo_replay_native_v1", args.output.absolute() / "payload")
        else:
            refinement_worker(root, args.output.absolute(), args.mode == "repeat-refinement")
    else:
        if args.commit is None or args.output is not None:
            parser.error("Parent requires commit and uses fixed output paths")
        run(root, args.commit)
