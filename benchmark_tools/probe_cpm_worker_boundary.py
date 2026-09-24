"""One frozen-worker diagnostic terminating before high-CPM optimization."""

import argparse
import faulthandler
import json
import os
from pathlib import Path
import subprocess
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
sys.path.insert(0, str(Path(__file__).resolve().parent))
from benchmark_tools.prepare_ob_candidate_neighborhood import check, record
from benchmark_tools.probe_cpm_high_constructor import AUDIT_SHA, HELPERS, FILES, require_audit, require_result

PRIOR_SHA = "2e1aa3998e29c14c22babbd859dab45c5692d4f4e4aa1ec094982c3154c628c5"
PROTOCOL_SHA = "69425260bfa456036189ea9ae6b272c4025c91889510ff90dec8ee1d7cf2f930"
EXTRA_HELPERS = {"checked_python_pair_worker.py": "443f2cb6b1f6f1fcc4d69a4410f0eeddaaa3cfee83fd0013f475221f7a39684f",
                 "stop_before_leiden.py": "52ad4ac4bd21b9be8389640b47afc58cd73197b10817127f8a6d59e4bc04f633"}


def require_stop(result, adapter, expected):
    if (result["status"] != "stopped_before_optimizer" or result["optimizer_called"] is not False
            or result["accuracy_evaluated"] is not False or result["fingerprint"] != expected
            or result["saved"] != expected):
        raise ValueError("Invalid or changed worker-boundary observation")
    if (adapter["format"] != "python_pairs" or len(adapter["calls"]) != 1
            or adapter["calls"][0]["status"] != "constructor_returned"):
        raise ValueError("Missing completed Python-pair constructor")
    for name in ("native_vs_saved", "constructor_vs_saved", "native_vs_constructor"):
        if result["differences"][name]["different_edges"] != 0:
            raise ValueError("Worker-boundary endpoint mismatch")


def worker(root, payload):
    faulthandler.enable(all_threads=True)
    from repeat_qfo_saved_graph import worker as frozen_worker, set_worker_affinity
    from checked_python_pair_worker import python_pair_constructor
    from benchmark_tools.stop_before_leiden import stop_before_partition, DiagnosticStop
    set_worker_affinity([min(os.sched_getaffinity(0))])
    import igraph
    import leidenalg
    launcher = root / "benchmarks/work/publication_qfo_replay_native_v1"
    if Path.cwd() != launcher:
        raise ValueError("Wrong frozen-worker working directory")
    try:
        with python_pair_constructor(igraph, payload / "constructor_adapter.json"):
            with stop_before_partition(leidenalg, payload):
                frozen_worker(launcher, payload, native_boundary=False, expected_cpm_resolution=.12)
    except DiagnosticStop:
        if not (payload / "preoptimizer_stop.json").is_file():
            raise ValueError("Missing explicit diagnostic-stop observation")
        return
    raise RuntimeError("Frozen worker returned without expected diagnostic stop")


def run(root, output):
    from benchmark_tools.run_blast_recovery_batch import save_status
    from benchmark_tools.run_simulation_methods import read_frozen
    from benchmark_tools.verify_blast_recovery_panel import unique_records
    if output.exists() or output.is_symlink():
        raise FileExistsError(output)
    results = root / "benchmark_tools/results"
    audit_path = results / "qfo_cpm_high_failed_payload_audit_20260923.json"
    prior_path = results / "qfo_cpm_high_frozen_constructor_22121.json"
    protocol = results / "QFO_CPM_WORKER_BOUNDARY_PROTOCOL_20260923.md"
    audit = read_frozen(audit_path, AUDIT_SHA)
    prior = read_frozen(prior_path, PRIOR_SHA)
    require_audit(audit)
    require_result(prior["result"], "frozen_imports")
    if record(protocol)["sha256"] != PROTOCOL_SHA:
        raise ValueError("Changed diagnostic protocol")
    checked = [record(p) for p in (audit_path, prior_path, protocol, Path(__file__))]
    checked.extend([*audit["checked_records"], *prior["checked_records"], *prior["observations"], prior["worker_log"]])
    for name, digest in {**HELPERS, **EXTRA_HELPERS}.items():
        item = record(Path(__file__).with_name(name))
        if item["sha256"] != digest:
            raise ValueError("Changed diagnostic helper")
        checked.append(item)
    source = root / "benchmarks/results/qfo_parameter_cpm_replay_v3/cpm_high/clustering/cluster_3_profile_expanded/payload"
    inputs = [record(source / name) for name in FILES]
    if any(item not in audit["checked_records"] for item in inputs):
        raise ValueError("Graph input not covered by failed-payload audit")
    metadata_path = source / "metadata.json"
    checked.append(record(metadata_path))
    metadata = json.loads(metadata_path.read_text())
    from repeat_qfo_saved_graph import validate_clustering_metadata
    validate_clustering_metadata(metadata, .12)
    checked = unique_records(checked)
    for item in checked:
        check(item)
    output.mkdir()
    payload = output / "payload"
    payload.mkdir()
    for name in FILES:
        (payload / name).symlink_to(source / name)
    metadata["output_directory"] = str(output)
    save_status(payload / "metadata.json", metadata)
    checked.append(record(payload / "metadata.json"))
    launcher = root / "benchmarks/work/publication_qfo_replay_native_v1"
    env = os.environ.copy()
    env.update(PYTHONPATH=str(launcher), PYTHONHASHSEED="0", PYTHONFAULTHANDLER="1",
               OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1")
    for name in ("PYTHONHOME", "LD_PRELOAD", "LD_LIBRARY_PATH"):
        env.pop(name, None)
    command = [sys.executable, "-B", "-X", "faulthandler", str(Path(__file__).resolve()),
               "--root", str(root), "--output", str(output), "--worker"]
    report = dict(status="worker_boundary_diagnostic_running", checked_records=checked, command=command,
        job_id=os.environ.get("SLURM_JOB_ID"), optimizer_called=False, accuracy_evaluated=False,
        rerun_authorized=False, publication_ready=False)
    save_status(output / "status.json", report)
    try:
        with (output / "worker.log").open("x") as log:
            done = subprocess.run(command, cwd=launcher, env=env, stdout=log, stderr=subprocess.STDOUT)
        report["returncode"] = done.returncode
        if done.returncode:
            raise RuntimeError(f"Frozen-worker diagnostic exited {done.returncode}")
        result = json.loads((payload / "preoptimizer_stop.json").read_text())
        adapter = json.loads((payload / "constructor_adapter.json").read_text())
        require_stop(result, adapter, prior["result"]["saved"])
        from repeat_qfo_saved_graph import check_worker
        snapshot = json.loads((payload / "worker_before.json").read_text())
        check_worker(snapshot, launcher, payload, {key: env[key] for key in (
            "PYTHONPATH", "PYTHONHASHSEED", "OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS")},
            expected_cpm_resolution=.12)
        observed_records = [*snapshot["modules"].values(), *snapshot["native_libraries"], snapshot["python"]]
        checked = unique_records([*checked, *observed_records])
        report["checked_records"] = checked
        if set(p.name for p in output.iterdir()) != {"payload", "worker.log", "status.json"}:
            raise ValueError("Unexpected diagnostic output, possible clustering artifact")
        for item in checked:
            check(item)
        report.update(status="frozen_worker_stopped_before_optimizer_unscored", result=result,
            observations=[record(payload / name) for name in (
                "worker_before.json", "constructor_adapter.json", "preoptimizer_stop.json")])
    except BaseException as error:
        report.update(status="worker_boundary_diagnostic_failed", error_type=type(error).__name__, error=str(error))
        raise
    finally:
        if (output / "worker.log").is_file():
            report["worker_log"] = record(output / "worker.log")
        save_status(output / "status.json", report)
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--worker", action="store_true")
    args = parser.parse_args()
    if args.worker:
        worker(args.root.resolve(), args.output.resolve() / "payload")
    else:
        run(args.root.resolve(), args.output.absolute())
