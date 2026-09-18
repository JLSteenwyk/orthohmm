"""Check native inference using copied, prevalidated fixture indexes without rebuilding."""

import argparse
import json
from pathlib import Path
import subprocess
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.prepare_ob_candidate_neighborhood import record, check
from benchmark_tools.prepare_qfo_corrected_orthomcl_native import write_sources, TOOL, RUNTIME_SHA, HELPERS_SHA
from benchmark_tools.probe_orthomcl_bpo_parity import expected_indexes
from benchmark_tools.probe_orthomcl_native_inference import partition
from benchmark_tools.stage_orthomcl_native_inputs import stage, NAMES
from benchmark_tools.run_qfo_corrected_blast import environment
from benchmark_tools.run_simulation_methods import read_frozen
from benchmark_tools.snapshot_runtime_trees import verify
from benchmark_tools.verify_orthomcl_python_runtime import verify_runtime
from benchmark_tools.orthomcl_matrix_to_pairwise import load_species


def run(root, output):
    if output.exists():
        raise FileExistsError(output)
    runtime_before = verify_runtime(root)
    results = root / "benchmark_tools/results"
    manifests = [read_frozen(results / name, sha) for name, sha in (
        ("qfo_corrected_orthomcl_perl_runtime_20260918.json", RUNTIME_SHA),
        ("qfo_corrected_orthomcl_system_helpers_20260918.json", HELPERS_SHA))]
    for manifest in manifests:
        verify(manifest)
    old_path = results / "orthomcl_native_inference_probe_20260918.json"
    old = read_frozen(old_path, "838eb833fe0b18b429db1c9bb1ad759313638dbbd2dcbcde7eb1c8270f21dc4e")
    previous = old["arms"]["native_serial"]
    originals = {Path(r["path"]).name: r for r in previous["outputs"] if Path(r["path"]).parent.name == "data"}
    inputs = {key: originals[name] for key, name in NAMES.items()}
    checked = [record(old_path), *inputs.values(), record(__file__),
               record(Path(__file__).with_name("stage_orthomcl_native_inputs.py"))]
    for item in checked:
        check(item)
    bpo, gg = Path(inputs["bpo"]["path"]), Path(inputs["species"]["path"])
    offsets, ranges = expected_indexes(bpo)
    index_summary = {"status": "native_bpo_indexes_verified", "records": len(offsets) - 1,
                     "queries": len(ranges), "offset_entries_including_eof": len(offsets), "bpo_bytes": bpo.stat().st_size}
    owners = load_species(gg)
    output.mkdir(parents=True, exist_ok=False)
    report = {"status": "running", "checked_records": checked, "runtime_before": runtime_before,
              "accuracy_admitted": False, "publication_ready": False}
    try:
        staged = stage(inputs, output / "inputs", len(owners), len(set(owners.values())), index_summary)
        report["staging"] = staged
        mtimes = {key: Path(item["path"]).stat().st_mtime_ns for key, item in staged["staged"].items()}
        report["configured_sources"] = write_sources(output / "tool", output / "inputs", 2, True)
        argv = [str(TOOL / "venv_orthomcl/bin/perl"), "-I" + str(output / "tool"),
                str(Path(__file__).with_name("run_orthomcl_perl_script.pl")), str(output / "tool/orthomcl.pl"),
                "--mode", "4", "--bpo_file", str(output / "inputs/all.bpo"), "--gg_file", str(output / "inputs/all.gg")]
        report.update(command=argv, cwd=str(output), environment={**environment(), "ORTHOMCL_PAIR_WORKERS": "2"})
        with (output / "native.log").open("xb") as log:
            done = subprocess.run(argv, cwd=output, env=report["environment"], stdout=log, stderr=subprocess.STDOUT, timeout=120)
        report["exit_code"] = done.returncode
        if done.returncode:
            raise ValueError("Staged native inference failed")
        native = list((output / "tool").glob("*/all_orthomcl.out"))
        baseline = [r for r in previous["outputs"] if Path(r["path"]).name == "all_orthomcl.out"]
        if len(native) != 1 or len(baseline) != 1:
            raise ValueError("Require unique fixture partitions")
        check(baseline[0])
        groups = partition(native[0], set(owners))
        if groups != partition(Path(baseline[0]["path"]), set(owners)):
            raise ValueError("Staged native partition differs from serial baseline")
        for item in staged["staged"].values():
            check(item)
        if mtimes != {key: Path(item["path"]).stat().st_mtime_ns for key, item in staged["staged"].items()}:
            raise ValueError("Native execution rewrote a staged input or index")
        for manifest in manifests:
            verify(manifest)
        for item in checked:
            check(item)
        report.update(status="staged_native_fixture_partition_and_input_preservation_verified",
                      groups=len(groups), grouped_proteins=sum(map(len, groups)), input_mtimes_ns=mtimes,
                      runtime_after=verify_runtime(root), outputs=[record(p) for p in sorted(output.rglob("*")) if p.is_file()],
                      limitations=["Small fixture only, not corrected production inference or 64-worker schedule validation.",
                                   "Baseline itself differs from the bundled historical example; that discrepancy remains retained."])
    except BaseException as error:
        report.update(status="failed", error_type=type(error).__name__, error=str(error))
        raise
    finally:
        with (output / "report.json").open("x") as stream:
            json.dump(report, stream, indent=2, sort_keys=True)
            stream.write("\n")
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("root", "output"):
        parser.add_argument("--" + name, type=Path, required=True)
    args = parser.parse_args()
    result = run(args.root.resolve(), args.output.resolve())
    print(json.dumps({"status": result["status"], "groups": result["groups"]}))
