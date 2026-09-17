"""Admit the fixed six-worker direct-construction import panel."""

import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import numpy as np
from benchmark_tools.prepare_ob_candidate_neighborhood import check, record
from benchmark_tools.probe_leiden_boundary import saved_fingerprint
from benchmark_tools.run_simulation_methods import read_frozen
from benchmark_tools.verify_qfo_replay_launcher import verify
from benchmark_tools.verify_ygob_validation import require_completed_job

ADMISSION_SHA = "499f90e06b6da311356c80f435149ba0cf132f75b13e2eacf720e02ea8f2382e"
MODES = ("minimal_imports", "frozen_imports")


def reconstructed_hash(differences, saved, sources, targets):
    delta = differences["native_vs_saved"]
    if (differences["constructor_dtype"] != "int32" or differences["constructor_c_contiguous"] is not True
            or differences["constructor_vs_saved"] != {"different_edges": 0, "examples": []}
            or differences["native_vs_constructor"] != delta
            or type(delta["different_edges"]) is not int
            or delta["different_edges"] != len(delta["examples"]) or len(delta["examples"]) > 20):
        raise ValueError("Invalid or incomplete stage difference evidence")
    patches = {}
    for example in delta["examples"]:
        index, observed = example["edge_index"], example["left"]
        if type(index) is not int or not 0 <= index < len(sources) or index in patches:
            raise ValueError("Invalid or repeated edge index")
        expected = sorted([int(sources[index]), int(targets[index])])
        if (example["right"] != expected or observed == expected or len(observed) != 2
                or observed != sorted(observed)
                or any(type(v) is not int or not 0 <= v < saved["vertices"] for v in observed)):
            raise ValueError("Invalid endpoint witnesses")
        patches[index] = observed
    digest = hashlib.sha256()
    for start in range(0, len(sources), 100000):
        end = min(start + 100000, len(sources))
        pairs = np.column_stack((sources[start:end], targets[start:end])).astype("<i8")
        pairs.sort(axis=1)
        for index, observed in patches.items():
            if start <= index < end:
                pairs[index - start] = observed
        digest.update(pairs.tobytes())
    return digest.hexdigest()


def check_stages(result, mode, saved, sources, targets):
    if (result["status"] != "direct_construction_observed" or result["mode"] != mode
            or result["accuracy_evaluated"] is not False or result["optimizer_called"] is not False
            or result["saved"] != saved):
        raise ValueError("Invalid direct observation")
    before, after = result["before_weights"], result["after_weights"]
    if any(before[key] != saved[key] for key in ("vertices", "edges", "directed")):
        raise ValueError("Pre-weight graph shape differs")
    before_hash = reconstructed_hash(before["differences"], saved, sources, targets)
    after_hash = reconstructed_hash(after["differences"], saved, sources, targets)
    if after["fingerprint"] != {**saved, "ordered_endpoints_sha256": after_hash}:
        raise ValueError("Post-weight native fingerprint not reconstructed")
    return {"before_mismatched_edges": before["differences"]["native_vs_saved"]["different_edges"],
            "after_mismatched_edges": after["differences"]["native_vs_saved"]["different_edges"],
            "before_witness_implied_hash": before_hash, "after_native_hash_reconstructed": after_hash,
            "stage_witnesses_equal": before["differences"] == after["differences"]}


def check_panel(report):
    if (report["status"] != "six_direct_workers_complete_unscored" or report["job_id"] != "21326"
            or report["optimizer_called"] is not False or report["accuracy_evaluated"] is not False
            or [(row["index"], row["mode"]) for row in report["workers"]] !=
            [(i, mode) for i in range(3) for mode in MODES]
            or any(row["exit_code"] != 0 for row in report["workers"])):
        raise ValueError("Incomplete or wrong direct-stage panel")


def admit(root, output, report_sha256):
    if output.exists():
        raise FileExistsError(output)
    path = root / "benchmarks/results/qfo_direct_graph_v1/results.json"
    report = read_frozen(path, report_sha256)
    check_panel(report)
    accounting = subprocess.check_output(["sacct", "-j", "21326", "--parsable2",
        "--format=JobIDRaw,State,ExitCode,Elapsed"], text=True)
    scheduler = require_completed_job(accounting, 21326)
    admission_path = root / "benchmark_tools/results/qfo_construction_verified_20260916.json"
    prior = read_frozen(admission_path, ADMISSION_SHA)
    if report["admission"] != record(admission_path) or report["graph_inputs"] != prior["native_report"]["graph_inputs"]:
        raise ValueError("Changed admitted graph")
    executor = root / "benchmarks/work/publication_qfo_direct_graph_v1"
    revision = subprocess.check_output(["git", "-C", str(executor), "rev-parse", "HEAD"], text=True).strip()
    expected = subprocess.check_output(["git", "-C", str(root), "rev-parse", "ab364e4^{commit}"], text=True).strip()
    if revision != expected or report["source"] != record(executor / "benchmark_tools/probe_qfo_direct_graph.py"):
        raise ValueError("Changed executor")
    subprocess.run(["git", "-C", str(executor), "diff", "--exit-code", "HEAD", "--", "benchmark_tools"], check=True)
    launcher = root / "benchmarks/work/publication_qfo_replay_native_v1"
    if verify(root / "benchmarks/work/publication_method_native_v2", launcher,
              root / "benchmark_tools/results/publication_native_runtime_20260916.json") != report["runtime"]:
        raise ValueError("Changed frozen runtime")
    records = [record(path), report["source"], report["admission"], *prior["provenance_checked"]]
    first, summaries = {}, []
    common = ("cpu_affinity", "cwd", "host", "platform", "python", "environment", "versions")
    overrides = {"PYTHONPATH": str(launcher), "PYTHONHASHSEED": "0", "OMP_NUM_THREADS": "1",
                 "OPENBLAS_NUM_THREADS": "1", "MKL_NUM_THREADS": "1"}
    for row in report["workers"]:
        mode = row["mode"]
        payload = path.parent / f"{mode}_{row['index']}" / "payload"
        result = json.loads((payload / "result.json").read_text())
        worker = json.loads((payload / "worker_before.json").read_text())
        if (row["observation"] != record(payload / "result.json") or row["result"] != result
                or result["snapshot"] != record(payload / "worker_before.json")
                or json.loads((payload / "before_weights.json").read_text()) != result["before_weights"]):
            raise ValueError("Preserved native files disagree")
        if (worker["mode"] != mode or worker["source"] != report["source"] or worker["inputs"] != report["graph_inputs"]
                or worker["accuracy_evaluated"] is not False or worker["optimizer_called"] is not False
                or len(worker["cpu_affinity"]) != 1 or worker["cwd"] != str(launcher)
                or any(worker["environment"][key] != value for key, value in overrides.items())):
            raise ValueError("Worker execution identity differs")
        scientific = [name for name in worker["modules"] if name.split(".")[0] in {"orthohmm", "leidenalg"}]
        if mode == "minimal_imports" and scientific:
            raise ValueError("Minimal worker imported scientific modules")
        if mode == "frozen_imports" and Path(worker["modules"]["orthohmm.leiden_worker"]["path"]) != launcher / "orthohmm/leiden_worker.py":
            raise ValueError("Wrong frozen import")
        if first and any(worker[key] != next(iter(first.values()))[key] for key in common):
            raise ValueError("Common worker identity changed")
        if mode in first and any(worker[key] != first[mode][key] for key in ("modules", "native_libraries")):
            raise ValueError("Same-mode runtime identity changed")
        first.setdefault(mode, worker)
        saved = saved_fingerprint(payload)
        summary = check_stages(result, mode, saved, np.load(payload / "sources.npy", mmap_mode="r"),
                               np.load(payload / "targets.npy", mmap_mode="r"))
        summaries.append({"index": row["index"], "mode": mode, **summary})
        records.extend([row["observation"], result["snapshot"], record(payload / "before_weights.json"),
                        *worker["modules"].values(), *worker["native_libraries"], worker["python"]])
    unique = {}
    for item in records:
        if item["path"] in unique and unique[item["path"]] != item:
            raise ValueError("Conflicting file provenance")
        unique[item["path"]] = item
        check(item)
    result = {"status": "direct_stage_observations_verified", "publication_ready": False,
              "accuracy_evaluated": False, "optimizer_called": False, "scheduler": scheduler,
              "source_report": record(path), "source": record(__file__), "workers": summaries,
              "provenance_checked": list(unique.values()), "native_report": report,
              "limitations": ["Preserved observation audit, not an independent historical live-graph inspection.",
                  "Before-weight hashes are implied by complete bounded witnesses; only post-weight native hashes were recorded directly.",
                  "Direct setup and observation change allocation/timing; this does not isolate a library defect or hardware cause."]}
    output.mkdir(parents=True)
    (output / "results.json").write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--report-sha256", required=True)
    args = parser.parse_args()
    admit(args.root.resolve(), args.output.resolve(), args.report_sha256)
