"""One isolated constructor observation on the failed high-CPM graph, without optimization."""

import argparse
import faulthandler
import json
import os
from pathlib import Path
import subprocess
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent))
from prepare_ob_candidate_neighborhood import check, record

AUDIT_SHA = "836c39dd63fe491d62e80c382ded24443b0ae6923965116cf0a092f540746b6b"
HELPERS = {
    "probe_qfo_direct_graph.py": "6d8e2fe3799211f05528a6fe5714c50fa164c33fc562bfe9a84114b571004d4b",
    "probe_leiden_boundary.py": "c9dc0dd592e9a21c621b302500740de55fe14204228961766d4ac36fdc380f14",
    "repeat_qfo_saved_graph.py": "6d80d4e92eccf3f32a62072df39050ccd4573b3bd326e32977302830cbf1f0f4",
}
FILES = ("gene_names.txt", "sources.npy", "targets.npy", "weights.npy")


def require_audit(report):
    observed = report["observed"]
    if (report["status"] != "failed_high_cpm_payload_checked_without_execution"
            or report["rerun_authorized"] is not False
            or report["accuracy_evaluated"] is not False
            or observed["vertices"] != 984137 or observed["edges"] != 25501180
            or any(observed[k] != 0 for k in ("negative_weights", "zero_weights", "self_edges"))
            or observed["constructor_int32_bytes_sha256"] !=
            "93b36aad4916b85195c4dcc7dccb37d944cb911f7b051cf989ad1479729fded4"):
        raise ValueError("Wrong failed-graph audit")


def require_result(result):
    if (result["status"] != "direct_construction_observed"
            or result["mode"] != "minimal_imports" or result["edge_format"] != "python_pairs"
            or result["optimizer_called"] is not False or result["accuracy_evaluated"] is not False
            or result["after_weights"]["fingerprint"] != result["saved"]):
        raise ValueError("Constructor observation does not preserve the saved graph")
    for phase in ("before_weights", "after_weights"):
        differences = result[phase]["differences"]
        for comparison in ("native_vs_saved", "constructor_vs_saved", "native_vs_constructor"):
            if differences[comparison]["different_edges"] != 0:
                raise ValueError("Native constructor endpoint mismatch")


def run(root, output):
    if output.exists() or output.is_symlink():
        raise FileExistsError(output)
    audit_path = root / "benchmark_tools/results/qfo_cpm_high_failed_payload_audit_20260923.json"
    audit_record = record(audit_path)
    if audit_record["sha256"] != AUDIT_SHA:
        raise ValueError("Changed read-only audit")
    audit = json.loads(audit_path.read_text())
    require_audit(audit)
    helpers = [record(Path(__file__).with_name(name)) for name in HELPERS]
    if [r["sha256"] for r in helpers] != list(HELPERS.values()):
        raise ValueError("Changed constructor helpers")
    records = [audit_record, record(__file__), *helpers, *audit["checked_records"]]
    for item in records:
        check(item)
    source = root / "benchmarks/results/qfo_parameter_cpm_replay_v3/cpm_high/clustering/cluster_3_profile_expanded/payload"
    inputs = [record(source / name) for name in FILES]
    if any(item not in audit["checked_records"] for item in inputs):
        raise ValueError("Input not covered by failed-graph audit")
    output.mkdir(parents=True)
    payload = output / "payload"
    payload.mkdir()
    for name in FILES:
        (payload / name).symlink_to(source / name)
    launcher = root / "benchmarks/work/publication_qfo_replay_native_v1"
    env = os.environ.copy()
    env.update(PYTHONPATH=str(launcher), PYTHONHASHSEED="0", PYTHONFAULTHANDLER="1",
               OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1")
    for key in ("PYTHONHOME", "LD_PRELOAD", "LD_LIBRARY_PATH"):
        env.pop(key, None)
    command = [sys.executable, "-B", "-X", "faulthandler", str(Path(__file__).resolve()),
               "--root", str(root), "--output", str(output), "--worker-payload", str(payload)]
    report = dict(status="constructor_diagnostic_running", inputs=inputs, checked_records=records,
                  command=command, job_id=os.environ.get("SLURM_JOB_ID"),
                  accuracy_evaluated=False, optimizer_called=False, rerun_authorized=False,
                  publication_ready=False, limitations=[
                      "One fresh minimal-import constructor; not an exact replay of the failed worker's allocation history.",
                      "Success cannot rule out intermittent crashes or establish the original cause.",
                      "No optimizer, groups, predictions, default changes or authorization of full replay."])
    def save():
        with (output / "status.json").open("w") as stream:
            json.dump(report, stream, indent=2, sort_keys=True)
            stream.write("\n")
            stream.flush()
            os.fsync(stream.fileno())
    save()
    try:
        with (output / "worker.log").open("x") as log:
            completed = subprocess.run(command, cwd=launcher, env=env, stdout=log, stderr=subprocess.STDOUT)
        report["returncode"] = completed.returncode
        report["worker_log"] = record(output / "worker.log")
        for item in records:
            check(item)
        if completed.returncode:
            raise RuntimeError(f"Constructor subprocess exited {completed.returncode}; preserve its traceback")
        result = json.loads((payload / "result.json").read_text())
        require_result(result)
        report.update(status="constructor_diagnostic_complete_unscored", result=result,
                      observations=[record(payload / name) for name in
                                    ("worker_before.json", "before_weights.json", "result.json")])
    except BaseException as error:
        report.update(status="failed", error_type=type(error).__name__, error=str(error))
        raise
    finally:
        save()
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("root", "output"):
        parser.add_argument("--" + name, type=Path, required=True)
    parser.add_argument("--worker-payload", type=Path)
    args = parser.parse_args()
    if args.worker_payload is not None:
        faulthandler.enable(all_threads=True)
        from probe_qfo_direct_graph import worker
        worker(args.root.resolve(), args.worker_payload.resolve(), "minimal_imports", "python_pairs")
    else:
        run(args.root.resolve(), args.output.absolute())
