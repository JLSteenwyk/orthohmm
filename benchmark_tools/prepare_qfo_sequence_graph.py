"""Freeze both corrected sequence-control graph commands after resource review."""

import argparse
import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.admit_qfo_sequence_numeric import completed
from benchmark_tools.compare_qfo_search_coverage import frozen
from benchmark_tools.prepare_ob_candidate_neighborhood import record, check
from benchmark_tools.run_sequence_graph_control import graph_command
from benchmark_tools.run_simulation_methods import read_frozen
from benchmark_tools.verify_qfo_replay_launcher import verify

PLANNER = "928eba0430f2348647d30819d5a6fa23873cb19d"
VARIANTS = ("all_hits", "top100")


def validate_payload(report, source, scheduler, memory_gib):
    if (report["status"] != "admitted_qfo_graph_payload_estimated" or report["source"] != source
            or report["job_id"] != scheduler["JobIDRaw"]
            or set(report["variants"]) != set(VARIANTS)
            or any(report[key] is not False for key in (
                "graph_launched", "graph_feasibility_admitted", "accuracy_evaluated", "publication_ready"))):
        raise ValueError("Wrong admitted payload report")
    if set(memory_gib) != set(VARIANTS):
        raise ValueError("Require explicit resource review for both variants")
    for label in VARIANTS:
        size = memory_gib[label]
        result = report["variants"][label]
        estimate = result["estimate"]
        if (type(size) is not int or size <= 0
                or result["status"] != "frozen_rbnh_named_array_payload_bound"
                or estimate["genes"] != 984137 or estimate["species_slot_extent"] != 78
                or estimate["graph_feasibility_admitted"] is not False
                or estimate["total_peak_ram_bound_available"] is not False):
            raise ValueError("Wrong checkpoint dimensions or resource review")
        if size * 1024**3 <= estimate["named_array_payload_lower_bytes"]:
            raise ValueError("Allocation cannot contain even the modeled snapshot")
    # Passing this necessary check does not establish that total peak memory fits.


def prepare(root, payload_path, payload_sha, job, memory_gib, output, destination):
    if output.exists() or destination.exists():
        raise FileExistsError("Require fresh output root and plan")
    scheduler = completed(job, 2, "64G")
    planner = frozen(root, "publication_qfo_graph_payload_v1", PLANNER)
    payload = read_frozen(payload_path, payload_sha)
    validate_payload(payload, record(planner / "benchmark_tools/plan_qfo_sequence_graph_memory.py"),
                     scheduler, memory_gib)
    if completed(payload["admission_scheduler"]["JobIDRaw"], 2, "192G") != payload["admission_scheduler"]:
        raise ValueError("Numeric admission accounting changed")
    conversion_path = root / "benchmarks/results/qfo_sequence_numeric_v1/manifest.json"
    conversion_record = record(conversion_path)
    if conversion_record not in payload["checked_records"]:
        raise ValueError("Conversion not bound by payload estimate")
    conversion = read_frozen(conversion_path, conversion_record["sha256"])
    checked = [record(payload_path), *payload["checked_records"], payload["source"], *payload["helpers"]]
    core = root / "benchmarks/work/publication_method_native_v2"
    launcher = root / "benchmarks/work/publication_qfo_replay_native_v1"
    runtime_path = root / "benchmark_tools/results/publication_native_runtime_20260916.json"
    for item in checked:
        check(item)
    runtime = verify(core, launcher, runtime_path)
    variants = {}
    for label in VARIANTS:
        estimated = payload["variants"][label]
        manifest = estimated["checkpoint_audit"]["manifest"]
        expected = conversion["variants"][label]
        if manifest != expected["manifest"] or expected["cap"] != (None if label == "all_hits" else 100):
            raise ValueError("Variant checkpoint identity differs")
        checkpoint = Path(expected["checkpoint"])
        if manifest["path"] != str(checkpoint / "manifest.json"):
            raise ValueError("Checkpoint location differs")
        checked.extend([manifest, estimated["core"], estimated["source"], *estimated["checkpoint_files"]])
        names = [r for r in estimated["checkpoint_files"] if r["path"] == str(checkpoint / "gene_names.txt")]
        if len(names) != 1:
            raise ValueError("Require exact gene-order evidence")
        variants[label] = dict(checkpoint_manifest=manifest, gene_names=names[0], cap=expected["cap"],
            output_root=str(output / label), cpu=32, requested_memory_gib=memory_gib[label],
            native_command=graph_command(launcher, checkpoint, manifest["sha256"], output / label),
            expected_clustering_calls=["initial", "multipass"],
            expected_stages=["multipass", "multipass_refined"],
            expected_genes=984137, expected_species=78)
    for item in checked:
        check(item)
    report = dict(status="corrected_sequence_graph_commands_frozen_unrun", variants=variants,
        source=record(__file__), helpers=[record(Path(__file__).with_name(name)) for name in (
            "run_sequence_graph_control.py", "verify_qfo_replay_launcher.py",
            "admit_qfo_sequence_numeric.py", "compare_qfo_search_coverage.py")],
        payload_report=record(payload_path), payload_scheduler=scheduler,
        checked_records=checked, runtime=runtime, python=record(sys.executable), cwd=str(launcher),
        environment_overrides=dict(PYTHONPATH=str(launcher), PYTHONHASHSEED="0", OMP_NUM_THREADS="1",
                                   OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1"),
        execution_authorized=False, accuracy_evaluated=False, graph_feasibility_admitted=False,
        remaining_gates=["Use a separately frozen checked-clustering executor; never execute native_command alone.",
                         "Require two checked native boundaries and complete gene coverage for each variant.",
                         "Review allocations beyond the named-array snapshot; retain resource failures.",
                         "Independent output admission and benchmark scoring remain required.",
                         "No profile refinement, candidate expansion or reconciliation in this P0C0R0 control.",
                         "Incremental shared-host replay is not dedicated end-to-end timing."])
    with destination.open("x") as stream:
        json.dump(report, stream, indent=2, sort_keys=True)
        stream.write("\n")
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("root", "payload", "output-root", "manifest"):
        parser.add_argument("--" + name, type=Path, required=True)
    parser.add_argument("--payload-sha256", required=True)
    parser.add_argument("--payload-job", required=True)
    parser.add_argument("--all-hits-memory-gib", type=int, required=True)
    parser.add_argument("--top100-memory-gib", type=int, required=True)
    args = parser.parse_args()
    prepare(args.root.resolve(), args.payload.resolve(), args.payload_sha256, args.payload_job,
            dict(all_hits=args.all_hits_memory_gib, top100=args.top100_memory_gib),
            args.output_root.resolve(), args.manifest.resolve())
