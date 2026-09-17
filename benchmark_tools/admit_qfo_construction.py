"""Verify completed construction observations without rerunning native construction."""

import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import numpy as np
from benchmark_tools.admit_qfo_affinity import IDENTITY_KEYS
from benchmark_tools.prepare_ob_candidate_neighborhood import check, record
from benchmark_tools.probe_leiden_boundary import saved_fingerprint
from benchmark_tools.probe_qfo_construction import PREVIOUS_SHA
from benchmark_tools.repeat_qfo_saved_graph import check_worker
from benchmark_tools.run_simulation_methods import read_frozen
from benchmark_tools.verify_qfo_replay_launcher import verify
from benchmark_tools.verify_ygob_validation import require_completed_job

REPORT_SHA = "7f33724b03afe50049225c854f4da93e05c430cee2e9a8cddc8130d100914254"


def check_observation(row, saved, sources, targets):
    mode, result = row["mode"], row["result"]
    if (mode not in ("original_int32", "explicit_int64") or result["mode"] != mode
            or result["status"] != "construction_observed_without_optimization"
            or result["accuracy_evaluated"] is not False or result["optimizer_called"] is not False
            or result["saved"] != saved):
        raise ValueError("Invalid construction observation")
    dtype = "int32" if mode == "original_int32" else "int64"
    original, converted = result["differences"], result["converted_differences"]
    if result["converted_dtype"] != dtype:
        raise ValueError("Changed conversion dtype")
    for differences, expected_dtype in ((original, "int32"), (converted, dtype)):
        if (differences["constructor_dtype"] != expected_dtype or differences["constructor_c_contiguous"] is not True
                or differences["constructor_vs_saved"] != {"different_edges": 0, "examples": []}
                or differences["native_vs_constructor"] != original["native_vs_saved"]
                or differences["native_vs_saved"] != original["native_vs_saved"]):
            raise ValueError("Input integrity or native difference reports disagree")
    differences = original["native_vs_saved"]
    examples = differences["examples"]
    if differences["different_edges"] != len(examples) or len(examples) > 20:
        raise ValueError("Require complete bounded mismatch witnesses")
    patches = {}
    for example in examples:
        index = example["edge_index"]
        if type(index) is not int or not 0 <= index < len(sources) or index in patches:
            raise ValueError("Invalid or repeated mismatch index")
        expected = sorted([int(sources[index]), int(targets[index])])
        observed = example["left"]
        if (example["right"] != expected or observed == expected or len(observed) != 2
                or observed != sorted(observed)
                or any(type(v) is not int or not 0 <= v < saved["vertices"] for v in observed)):
            raise ValueError("Invalid mismatch endpoints")
        patches[index] = observed
    digest = hashlib.sha256()
    for start in range(0, len(sources), 100000):
        end = min(start + 100000, len(sources))
        endpoints = np.column_stack((sources[start:end], targets[start:end])).astype("<i8")
        endpoints.sort(axis=1)
        for index, pair in patches.items():
            if start <= index < end:
                endpoints[index - start] = pair
        digest.update(endpoints.tobytes(order="C"))
    if result["native"] != {**saved, "ordered_endpoints_sha256": digest.hexdigest()}:
        raise ValueError("Witnesses do not reconstruct full native graph fingerprint")
    witnesses = result["witnesses"]
    if len(witnesses) != len(examples):
        raise ValueError("Missing independent-access witnesses")
    for witness, example in zip(witnesses, examples):
        if (witness["edge_index"] != example["edge_index"] or sorted(witness["tuple"]) != example["left"]
                or witness["source_target"] != witness["tuple"] or witness["expected_pair_edge_id"] != -1):
            raise ValueError("Native edge-access witnesses disagree")
    return {"mode": mode, "index": row["index"], "mismatched_edges": len(examples),
            "complete_endpoint_hash_reconstructed": True}


def admit(root, output):
    if output.exists():
        raise FileExistsError(output)
    path = root / "benchmarks/results/qfo_construction_v1/results.json"
    report = read_frozen(path, REPORT_SHA)
    if (report["status"] != "six_construction_workers_complete_unscored" or report["job_id"] != "21323"
            or report["optimizer_called"] is not False or report["accuracy_evaluated"] is not False
            or [(row["index"], row["mode"]) for row in report["workers"]] !=
            [(i, mode) for i in range(3) for mode in ("original_int32", "explicit_int64")]):
        raise ValueError("Wrong or incomplete construction panel")
    accounting = subprocess.check_output(["sacct", "-j", "21323", "--parsable2",
        "--format=JobIDRaw,State,ExitCode,Elapsed"], text=True)
    scheduler = require_completed_job(accounting, 21323)
    previous = read_frozen(Path(report["previous"]["path"]), PREVIOUS_SHA)
    if report["graph_inputs"] != previous["graph_inputs"]:
        raise ValueError("Changed saved graph inputs")
    executor = root / "benchmarks/work/publication_qfo_construction_v1"
    revision = subprocess.check_output(["git", "-C", str(executor), "rev-parse", "HEAD"], text=True).strip()
    expected = subprocess.check_output(["git", "-C", str(root), "rev-parse", "f2827a6^{commit}"], text=True).strip()
    if revision != expected or report["source"] != record(executor / "benchmark_tools/probe_qfo_construction.py"):
        raise ValueError("Changed frozen executor")
    subprocess.run(["git", "-C", str(executor), "diff", "--exit-code", "HEAD", "--", "benchmark_tools"], check=True)
    launcher = root / "benchmarks/work/publication_qfo_replay_native_v1"
    runtime = verify(root / "benchmarks/work/publication_method_native_v2", launcher,
                     root / "benchmark_tools/results/publication_native_runtime_20260916.json")
    if runtime != report["runtime"]:
        raise ValueError("Frozen runtime differs")
    records = [record(path), report["source"], report["previous"], *report["graph_inputs"]]
    overrides = {"PYTHONPATH": str(launcher), "PYTHONHASHSEED": "0", "OMP_NUM_THREADS": "1",
                 "OPENBLAS_NUM_THREADS": "1", "MKL_NUM_THREADS": "1"}
    summaries, first = [], None
    for row in report["workers"]:
        payload = path.parent / f"{row['mode']}_{row['index']}" / "payload"
        if row["observation"] != record(payload / "construction.json") or row["snapshot"] != record(payload / "worker_before.json"):
            raise ValueError("Changed native observation files")
        if json.loads((payload / "construction.json").read_text()) != row["result"]:
            raise ValueError("Report differs from native observation")
        worker = json.loads((payload / "worker_before.json").read_text())
        check_worker(worker, launcher, payload, overrides)
        if worker["inputs"] != report["graph_inputs"] or len(worker["cpu_affinity"]) != 1:
            raise ValueError("Changed worker inputs or affinity")
        if first is not None and any(worker[key] != first[key] for key in (*IDENTITY_KEYS, "cpu_affinity")):
            raise ValueError("Worker software or CPU identity differs")
        first = worker
        saved = saved_fingerprint(payload)
        summaries.append(check_observation(row, saved, np.load(payload / "sources.npy", mmap_mode="r"),
                                           np.load(payload / "targets.npy", mmap_mode="r")))
        records.extend([row["observation"], row["snapshot"], *worker["modules"].values(),
                        *worker["native_libraries"], worker["python"], worker["observer"]])
    unique = {}
    for item in records:
        if item["path"] in unique and unique[item["path"]] != item:
            raise ValueError("Conflicting file provenance")
        unique[item["path"]] = item
        check(item)
    result = {"status": "construction_observations_verified", "accuracy_evaluated": False,
              "optimizer_called": False, "publication_ready": False, "scheduler": scheduler,
              "source_report": record(path), "source": record(__file__), "workers": summaries,
              "provenance_checked": list(unique.values()), "native_report": report,
              "limitations": ["Validates preserved observations and complete mismatch/hash consistency, not independent live native graph inspection.",
                  "One explicit-int64 worker still mismatched; integer-width conversion alone is not a sufficient fix in this experiment.",
                  "Five matching workers do not establish general correctness; this panel does not isolate a library defect or hardware cause."]}
    output.mkdir(parents=True)
    (output / "results.json").write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    admit(args.root.resolve(), args.output.resolve())
