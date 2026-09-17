"""Independently admit the fixed checked QfO cached replay, retaining historical disagreement."""

import argparse
import json
from pathlib import Path
import subprocess
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.admit_qfo_checked_repeats import reconstruct
from benchmark_tools.prepare_ob_candidate_neighborhood import check, record
from benchmark_tools.run_qfo_publication_replay import command_for, compare_partition
from benchmark_tools.run_simulation_methods import read_frozen
from benchmark_tools.verify_qfo_replay_launcher import verify
from benchmark_tools.verify_ygob_validation import require_completed_job

ADMISSION_SHA = "78f51f5ce703caf39307c5e518ad737acdf655dbe0926a34f2c20bbdec6d03f1"
CALLS = ("initial", "multipass", "profile_base", "profile_expanded")
STAGES = ("multipass", "multipass_refined", "profiles", "profiles_refined")
RUNS = {"v1": (21329, "96333fd"), "v2": (21333, "c93eb2c8cf5671b9e99b3e534f68d11fac6282d7")}


def check_inventory(parent, worker, replay, version="v1"):
    if version not in RUNS:
        raise ValueError("Unknown fixed replay version")
    if (parent["status"] != "full_checked_replay_complete_unscored" or parent["job_id"] != str(RUNS[version][0])
            or parent["exit_code"] != 0 or parent["accuracy_evaluated"] is not False
            or worker["status"] != "checked_full_replay_returned" or worker["accuracy_evaluated"] is not False
            or [(r["index"], r["stage"]) for r in worker["calls"]] != list(enumerate(CALLS))
            or any(r["status"] != "checked" or r["exit_code"] != 0 or r["accuracy_evaluated"] is not False for r in worker["calls"])
            or [r["label"] for r in replay["stages"]] != list(STAGES)
            or replay["counts"]["genes"] != 976504 or replay["counts"].get("profiles_built", 0) <= 0):
        raise ValueError("Require complete fixed checked full replay inventory")
    expected = {"accuracy_profile": "high_sensitivity", "cpm_resolution": .1, "profile_expansion": True,
                "profile_iterations": 1, "jackknife_profile_thresholds": False,
                "jackknife_single_copy_profiles": False, "profile_min_species": 1, "matrix": "BLOSUM62", "leiden_seed": 4}
    if replay["parameters"] != expected or any("official_orthobench" in r for r in replay["stages"]):
        raise ValueError("Replay settings or scoring scope changed")
    if version == "v2":
        for row in worker["calls"]:
            overrides = {k: "1" for k in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS")}
            inherited = {**overrides, "OMP_NUM_THREADS": "32" if row["index"] >= 2 else "1"}
            if row.get("thread_environment") != {"inherited": inherited, "child_overrides": overrides}:
                raise ValueError("Explicit clustering child thread isolation differs")


def admit(root, output, digest, version="v1"):
    if output.exists():
        raise FileExistsError(output)
    if version not in RUNS:
        raise ValueError("Unknown fixed replay version")
    job_id, revision_pin = RUNS[version]
    directory = root / f"benchmarks/results/qfo_checked_full_replay_{version}"
    path = directory / "results.json"
    report = read_frozen(path, digest)
    worker = json.loads((directory / "checked_worker.json").read_text())
    replay = json.loads((directory / "replay.json").read_text())
    check_inventory(report, worker, replay, version)
    accounting = subprocess.check_output(["sacct", "-j", str(job_id), "--parsable2", "--format=JobIDRaw,State,ExitCode,Elapsed,MaxRSS"], text=True)
    scheduler = require_completed_job(accounting, job_id)
    executor = root / f"benchmarks/work/publication_qfo_checked_full_replay_{version}"
    revision = subprocess.check_output(["git", "-C", str(executor), "rev-parse", "HEAD"], text=True).strip()
    expected_revision = subprocess.check_output(["git", "-C", str(root), "rev-parse", revision_pin + "^{commit}"], text=True).strip()
    if revision != expected_revision or report["executor_commit"] != revision:
        raise ValueError("Frozen executor revision changed")
    subprocess.run(["git", "-C", str(executor), "diff", "--exit-code", "HEAD", "--", "benchmark_tools"], check=True)
    launcher = root / "benchmarks/work/publication_qfo_replay_native_v1"
    core = root / "benchmarks/work/publication_method_native_v2"
    runtime = root / "benchmark_tools/results/publication_native_runtime_20260916.json"
    if verify(core, launcher, runtime) != report["runtime"]:
        raise ValueError("Frozen native runtime changed")
    source = record(executor / "benchmark_tools/run_qfo_checked_full_replay.py")
    if report["source"] != source or worker["source"] != source:
        raise ValueError("Wrong replay executor source")
    expected_sources = [record(executor / "benchmark_tools" / name) for name in
        ("run_qfo_checked_full_replay.py", "checked_replay_interceptor.py", "validate_checked_replay_payload.py",
         "checked_replay_payload_worker.py", "checked_python_pair_worker.py", "repeat_qfo_saved_graph.py", "probe_leiden_boundary.py")]
    if report["executor_sources"] != expected_sources or worker["helpers"] != expected_sources[1:4]:
        raise ValueError("Incomplete or changed executor source inventory")
    expected_worker_command = [sys.executable, source["path"], "--root", str(root), "--output", str(directory), "--replay-worker"]
    if report["worker_command"] != expected_worker_command:
        raise ValueError("Parent worker command differs")
    command = command_for(root, launcher, directory)
    if report["replay_command"] != command or worker["command"] != command or replay["command"] != command:
        raise ValueError("Replay command differs")
    if worker["replay_source"] != record(launcher / "benchmark_tools/replay_high_sensitivity.py") or replay["source"] != worker["replay_source"]:
        raise ValueError("Wrong replay scientific source")
    if report["worker"] != record(directory / "checked_worker.json") or report["replay"] != record(directory / "replay.json"):
        raise ValueError("Changed parent-linked outputs")
    admission_path = root / "benchmark_tools/results/qfo_checked_repeats_verified_20260917.json"
    prior = read_frozen(admission_path, ADMISSION_SHA)
    if report["admission"] != record(admission_path):
        raise ValueError("Changed initial-repeat admission")
    records = [record(path), report["worker"], report["replay"], report["admission"], *report["executor_sources"],
               *worker["helpers"], *prior["provenance_checked"], record(directory / "replay_command.json")]
    summaries = []
    for row in worker["calls"]:
        stage_dir = directory / "clustering" / f"cluster_{row['index']}_{row['stage']}"
        payload = stage_dir / "payload"
        if json.loads((stage_dir / "execution.json").read_text()) != row or row["manifest"] != record(stage_dir / "payload_manifest.json"):
            raise ValueError("Preserved stage execution disagrees with worker report")
        manifest = json.loads((stage_dir / "payload_manifest.json").read_text())
        if manifest["index"] != row["index"] or manifest["stage"] != row["stage"] or manifest["output_directory"] != str(directory / "replay"):
            raise ValueError("Stage manifest identity differs")
        actual_inputs = [record(payload / name) for name in ("gene_names.txt", "sources.npy", "targets.npy", "weights.npy", "metadata.json")]
        if manifest["inputs"] != actual_inputs or len(manifest["original_payload"]) != 5:
            raise ValueError("Stage payload inventory differs")
        if any(actual_inputs[0][k] != prior["native_report"]["graph_inputs"][0][k] for k in ("bytes", "sha256")):
            raise ValueError("Stage gene order differs from admitted QfO universe")
        for original, copied in zip(manifest["original_payload"], actual_inputs):
            if any(original[k] != copied[k] for k in ("bytes", "sha256")):
                raise ValueError("Original and preserved payload bytes differ")
        command = [replay["command"][0], str(executor / "benchmark_tools/checked_replay_payload_worker.py"), "--root", str(root),
                   "--payload", str(payload), "--manifest", str(stage_dir / "payload_manifest.json"),
                   "--manifest-sha256", row["manifest"]["sha256"]]
        if row["command"] != command:
            raise ValueError("Adapted native command differs")
        # The replay's mutable output now holds the last call; validate retained native
        # observations independently instead of asking the in-run validator to reread it.
        saved, oriented, universe = reconstruct(payload)
        if saved["vertices"] != 976504:
            raise ValueError("Incomplete QfO graph universe")
        from benchmark_tools.admit_qfo_checked_repeats import check_native
        boundary = json.loads((payload / "native_boundary.json").read_text())
        adapter = json.loads((payload / "constructor_adapter.json").read_text())
        check_native(boundary, adapter, saved, oriented)
        validation = row["validation"]
        observed = json.loads((payload / "worker_before.json").read_text())
        if validation["status"] != "payload_checked" or validation["saved_graph"] != saved or validation["worker"] != observed:
            raise ValueError("Preserved validation disagrees with native records")
        metadata = json.loads((payload / "metadata.json").read_text())
        expected_metadata = {"cpm_resolution": .1, "seed": 4, "include_isolates": True, "output_directory": str(directory / "replay")}
        overrides = {"PYTHONPATH": str(launcher), "PYTHONHASHSEED": "0", "OMP_NUM_THREADS": "1", "OPENBLAS_NUM_THREADS": "1", "MKL_NUM_THREADS": "1"}
        if (metadata != expected_metadata or observed["metadata"] != metadata
                or observed["status"] != "before_native_clustering" or observed["accuracy_evaluated"] is not False
                or any(observed["environment"][k] != v for k, v in overrides.items()) or not observed["native_libraries"]):
            raise ValueError("Worker settings or environment differ")
        provenance = json.loads((payload / "checked_payload_provenance.json").read_text())
        helpers = [expected_sources[i] for i in (5, 4, 6)]
        if (provenance["source"] != expected_sources[3] or provenance["helpers"] != helpers
                or provenance["manifest"] != row["manifest"] or provenance["admission"] != report["admission"]
                or provenance["inputs"] != actual_inputs or provenance["stage"] != row["stage"]
                or provenance["accuracy_evaluated"] is not False or observed["observer"] != helpers[0]):
            raise ValueError("Native helper provenance differs")
        for name in ("orthohmm.leiden_worker", "orthohmm.externals", "orthohmm.helpers"):
            if observed["modules"][name] != record(launcher / (name.replace(".", "/") + ".py")):
                raise ValueError("Native module differs")
        if observed["inputs"] != manifest["inputs"][:4] or len(observed["cpu_affinity"]) != 1 or observed["cwd"] != str(launcher):
            raise ValueError("Native worker input, affinity or working directory differs")
        partition = stage_dir / "partition.txt"
        if row["partition"] != record(partition):
            raise ValueError("Retained native partition changed")
        coverage = compare_partition(partition, partition, universe)
        if coverage["observed_groups"] != validation["groups"] or validation["genes"] != len(universe):
            raise ValueError("Native partition counts differ")
        if row["index"] == 0:
            for a, b in zip(manifest["inputs"][:4], prior["native_report"]["graph_inputs"]):
                if any(a[k] != b[k] for k in ("bytes", "sha256")):
                    raise ValueError("Initial payload differs from admitted graph")
        records.extend([row["manifest"], row["partition"], record(stage_dir / "execution.json"), record(stage_dir / "worker.log"),
                        *actual_inputs, *observed["modules"].values(), *observed["native_libraries"], observed["python"],
                        *[record(payload / name) for name in ("worker_before.json", "native_boundary.json", "constructor_adapter.json", "checked_payload_provenance.json")],
                        *validation["provenance_checked"]])
        summaries.append({"index": row["index"], "stage": row["stage"], "saved_graph": saved, "coverage": coverage})
    if (replay["counts"]["rbnh_edges"] != summaries[0]["saved_graph"]["edges"]
            or replay["counts"]["multipass_edges"] != summaries[1]["saved_graph"]["edges"]):
        raise ValueError("Replay edge counts disagree with preserved graphs")
    coverage = []
    for stage in replay["stages"]:
        check(stage["output"])
        p = Path(stage["output"]["path"])
        coverage.append({"stage": stage["label"], **compare_partition(p, p, universe)})
        records.append(stage["output"])
    if coverage != report["coverage"]:
        raise ValueError("Final stage coverage comparisons differ")
    by_label = {r["label"]: Path(r["output"]["path"]) for r in replay["stages"]}
    for index, label in ((1, "multipass"), (3, "profiles")):
        if not compare_partition(Path(worker["calls"][index]["partition"]["path"]), by_label[label], universe)["byte_equal"]:
            raise ValueError("Frozen replay stage not copied from corresponding checked partition")
    before = json.loads((directory / "inputs_before.json").read_text())
    after = json.loads((directory / "inputs_after.json").read_text())
    if before != after:
        raise ValueError("Original input audits disagree")
    for key in ("inputs_before", "inputs_after"):
        if report[key] != record(directory / (key + ".json")):
            raise ValueError("Input audit file changed")
        records.append(report[key])
    output.mkdir(parents=True)
    fresh = output / "inputs_rechecked.json"
    subprocess.run([sys.executable, str(executor / "benchmark_tools/audit_qfo_replay_inputs.py"), "--root", str(root), "--output", str(fresh)], check=True)
    if json.loads(fresh.read_text()) != before:
        raise ValueError("Fresh independent input audit differs")
    historical = compare_partition(Path(before["target_partition"]["path"]), by_label["profiles_refined"], universe)
    initial = compare_partition(Path(prior["native_report"]["repeats"][0]["partition"]["path"]), Path(worker["calls"][0]["partition"]["path"]), universe)
    if initial != report["initial_versus_checked_repeat"]:
        raise ValueError("Initial partition comparison differs")
    records.extend([record(fresh), before["target_partition"]])
    unique = {}
    for item in records:
        if item["path"] in unique and unique[item["path"]] != item:
            raise ValueError("Conflicting provenance records")
        unique[item["path"]] = item
        check(item)
    if verify(core, launcher, runtime) != report["runtime"]:
        raise ValueError("Native runtime changed during admission")
    result = {"status": "checked_full_replay_verified", "publication_ready": False, "accuracy_evaluated": False,
              "run_version": version,
              "source": record(__file__), "source_report": record(path), "scheduler": scheduler, "scheduler_accounting": accounting,
              "native_report": report, "native_replay": replay, "clustering": summaries, "coverage": coverage,
              "initial_versus_checked_repeat": initial, "final_versus_historical": historical,
              "provenance_checked": list(unique.values()),
              "limitations": ["One cached replay, not general determinism or a controlled end-to-end resource measurement.",
                  "Historical disagreement is retained, not a reason for retry or accuracy-based selection.",
                  "Native hashes validate preserved observations, not retrospective live-memory inspection."]}
    (output / "results.json").write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--report-sha256", required=True)
    parser.add_argument("--run-version", choices=tuple(RUNS), default="v1")
    args = parser.parse_args()
    admit(args.root.resolve(), args.output.resolve(), args.report_sha256, args.run_version)
