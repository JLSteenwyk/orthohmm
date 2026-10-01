"""Validate private QfO phylogeny deployment against the entire fixed baseline."""

import argparse
import json
import os
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.prepare_ob_candidate_neighborhood import check, record
from benchmark_tools.qfo_private_phylogeny_environment import verify_baseline, control_command, inspect_launcher

PROTOCOL = "benchmark_tools/results/QFO_PRIVATE_PHYLOGENY_CONTROL_PROTOCOL_20261001.md"
OUTPUT = "benchmarks/results/qfo_private_phylogeny_control_v1"
NATIVE_FILES = ("orthohmm_root_hogs.tsv", "orthohmm_pairwise_orthologs.tsv",
    "orthohmm_pairwise_orthologs_confidence.tsv", "orthohmm_reconciliation_nodes.tsv",
    "orthohmm_hierarchical_orthogroups.tsv", "species_tree.rooted.nwk")


def compare_outputs(old, new):
    result = {}
    for name in NATIVE_FILES:
        historical, current = record(old / name), record(new / name)
        result[name] = {"historical": historical, "current": current,
                       "byte_equal": (historical["bytes"], historical["sha256"]) == (current["bytes"], current["sha256"])}
    return result


def run(root, protocol_sha):
    if (not os.environ.get("SLURM_JOB_ID") or os.environ.get("SLURM_CPUS_PER_TASK") != "32"
            or os.environ.get("SLURM_MEM_PER_NODE") != "196608" or os.environ.get("SLURM_JOB_NODELIST") != "bizon"
            or os.environ.get("SLURM_ARRAY_TASK_ID") or sys.executable != "/home/bizon/anaconda3/bin/python"):
        raise ValueError("Require standalone 32-CPU/192-GiB bizon control with historical verifier interpreter")
    output = root / OUTPUT
    if output.exists() or output.is_symlink():
        raise FileExistsError(output)
    protocol = record(root / PROTOCOL)
    if protocol["sha256"] != protocol_sha:
        raise ValueError("Unreviewed full private phylogeny control protocol")
    verified = verify_baseline(root)
    from benchmark_tools.run_simulation_methods import execution_environment, execute

    launcher, prepared = Path(verified["launcher"]), Path(verified["prepared"])
    argv, equivalence = control_command(verified["original"], launcher, prepared, verified["environment"], output)
    env, resolved = execution_environment(verified["environment"])
    env.update(verified["manifest"]["environment_overrides"], PYTHONPATH=str(launcher),
        PYTHONDONTWRITEBYTECODE="1", PYTHONNOUSERSITE="1", PYTHONPYCACHEPREFIX=str(output / "bytecode_cache"),
        NUMBA_CACHE_DIR=str(output / "numba_cache"), NUMBA_CACHE_LOCATOR_CLASSES="UserProvidedCacheLocator")
    for key in ("PYTHONHOME", "LD_PRELOAD", "LD_LIBRARY_PATH", "LD_AUDIT"):
        env.pop(key, None)
    helpers = [record(path) for path in sorted(Path(__file__).parent.glob("*.py"))]
    provenance = {"source": record(__file__), "helpers": helpers, "protocol": protocol,
        "verified": verified, "executed_argv": argv, "launcher_equivalence": equivalence, "resolved_tools": resolved,
        "cwd": str(launcher), "job_id": os.environ["SLURM_JOB_ID"],
        "scope": "Entire frozen QfO p1_c1_r1 private deployment parity; validated checkpoint reuse; unscored incremental shared-host control"}
    output.mkdir(parents=True, exist_ok=False)
    (output / "numba_cache").mkdir()
    (output / "preflight.json").write_text(json.dumps(provenance, indent=2, sort_keys=True) + "\n")
    postflight = {"status": "running", "accuracy_evaluated": False, "publication_ready": False,
                  "recovered_cpm_inference_authorized": False, "controlled_timing": False}
    cwd = Path.cwd()
    try:
        lookup = inspect_launcher(verified, env, output)
        method = {"argv": argv, "output": str(output / "output"), "metrics": str(output / "metrics.json")}
        inputs = {"status": "ready", "inputs": [{**r, "absolute_path": r["path"]} for r in verified["manifest"]["input_fastas"]]}
        label = "private_phylogeny_control"
        os.chdir(launcher)
        execution = execute({"label": label, "methods": {label: method}}, [label], env,
                            output / "execution", inputs, provenance)
        os.chdir(cwd)
        if execution["failed_methods"]:
            raise RuntimeError("Private QfO control failed; preserve without retry")
        if verify_baseline(root) != verified:
            raise ValueError("Private QfO control baseline/runtime changed during execution")
        for item in [provenance["source"], protocol, *helpers, *verified["checked_records"], *lookup["checked_records"]]:
            check(item)
        for pair in equivalence:
            check(pair["prepared"])
            check(pair["executed"])
        old = Path(argv[argv.index("--checkpoint-source") + 1]) / "orthohmm_phylogeny"
        comparisons = compare_outputs(old, output / "output/orthohmm_phylogeny")
        postflight["native_comparison"] = comparisons
        if not all(item["byte_equal"] for item in comparisons.values()):
            raise ValueError("Private QfO baseline native outputs differ; preserve negative control")
        for item in comparisons.values():
            check(item["historical"])
            check(item["current"])
        postflight.update(status="private_qfo_baseline_parity_complete_pending_admission")
    except BaseException as error:
        postflight.update(status="failed", error_type=type(error).__name__, error=str(error))
        raise
    finally:
        os.chdir(cwd)
        (output / "postflight.json").write_text(json.dumps(postflight, indent=2, sort_keys=True) + "\n")
    return postflight


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--protocol-sha256", required=True)
    args = parser.parse_args()
    run(args.root.resolve(), args.protocol_sha256)
