"""Run the fixed recovered high-CPM arm in the admitted private native runtime."""

import argparse
import json
import os
from pathlib import Path
import subprocess
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from benchmark_tools.prepare_ob_candidate_neighborhood import check, record

ADMISSION = "benchmarks/work/qfo_private_phylogeny_control_admission_v2_20261001.json"
ADMISSION_SHA = "6865539281f9cb8c89c04b8ca16e2277056b4cc2b3293c330a12a3115879af97"
READBACK = "benchmark_tools/results/qfo_private_phylogeny_admission_readback_20261001.json"
READBACK_SHA = "713718ada456ad731df02b18f3c22219d8a76a1d7c2c45365d63f7df4ead18fc"
READBACK_SOURCE_SHA = "982d6de892e790a5c568e2f3078d3f89e4df86ba15f4db7ff31e9b306e675717"
ADMISSION_COMMIT = "9f0fcb68b849aa8af6180c73c9ba82d3e4f4b8b1"
EXECUTOR = "benchmarks/work/qfo_private_phylogeny_admission_v2_executor_20261001"
PROTOCOL = "benchmark_tools/results/QFO_PRIVATE_CPM_PHYLOGENY_PROTOCOL_20261001.md"
OUTPUT = "benchmarks/results/qfo_cpm_private_phylogeny_v1/cpm_high"
PYTHON = "/home/bizon/anaconda3/bin/python"


def pinned(root, relative, sha):
    pin = record(root / relative)
    if pin["sha256"] != sha:
        raise ValueError("Private recovered-arm evidence changed: " + relative)
    return json.loads(Path(pin["path"]).read_bytes()), pin


def verify_private_control(root):
    source = record(root / "benchmark_tools/readback_qfo_private_phylogeny_control.py")
    if source["sha256"] != READBACK_SOURCE_SHA:
        raise ValueError("Private admission readback verifier changed")
    from benchmark_tools.readback_qfo_private_phylogeny_control import completed

    accounting = subprocess.check_output(["sacct", "-j", "22387,22389", "-P",
        "--format=JobID,JobIDRaw,State,ExitCode,Elapsed,AllocCPUS,ReqMem,NodeList"], text=True)
    scheduler = completed(accounting)
    report, admission_pin = pinned(root, ADMISSION, ADMISSION_SHA)
    readback, readback_pin = pinned(root, READBACK, READBACK_SHA)
    if (report["status"] != "private_qfo_phylogeny_deployment_admitted_unscored"
            or report["job_id"] != "22387" or report["scheduler"] != scheduler[0]
            or report["native_group_integrity"]["native_outputs_validated"] is not True
            or readback["status"] != "private_qfo_phylogeny_admission_independently_read_back_unscored"
            or readback["admission"] != admission_pin or readback["source"] != source
            or readback["source_revision"] != ADMISSION_COMMIT or readback["scheduler"] != scheduler
            or readback["native_pair_count"] != report["native_pair_count"]
            or readback["native_comparison"] != report["native_comparison"]
            or readback["cache_use"] != report["cache_use"]
            or any(item["recovered_cpm_inference_authorized"] is not True for item in (report, readback))
            or any(item[k] is not False for item in (report, readback) for k in
                ("accuracy_evaluated", "scoring_admitted", "controlled_timing", "publication_ready"))):
        raise ValueError("Private QfO deployment admission/readback disagrees")
    executor = root / EXECUTOR
    if subprocess.check_output(["git", "-C", str(executor), "rev-parse", "HEAD"], text=True).strip() != ADMISSION_COMMIT:
        raise ValueError("Private admission executor revision changed")
    subprocess.run(["git", "-C", str(executor), "diff", "--exit-code", "HEAD", "--", "benchmark_tools", "orthohmm"],
                   check=True, capture_output=True)
    submission = json.loads(Path(readback["submission"]["path"]).read_bytes())
    if (submission["job_id"] != "22389" or submission["executor"] != str(executor)
            or submission["executor_commit"] != ADMISSION_COMMIT
            or report["source"] != submission["source_records"][0]
            or readback["source_git_bindings"] != [{**item, "git_revision": ADMISSION_COMMIT}
                                                  for item in submission["source_records"]]):
        raise ValueError("Private admission source/Git bindings disagree")
    for item in readback["source_git_bindings"]:
        path = Path(item["path"])
        if record(path) != {k: item[k] for k in ("path", "bytes", "sha256")} or path.read_bytes() != subprocess.check_output(
                ["git", "-C", str(executor), "show", f"{ADMISSION_COMMIT}:{path.relative_to(executor)}"]):
            raise ValueError("Private admission source differs from frozen Git blob")
    records = [admission_pin, readback_pin, source, readback["submission"],
        submission["parent_submission"], submission["preserved_failed_attempt"],
        *submission["source_records"], *report["checked_records"]]
    unique = {}
    for item in records:
        if item["path"] in unique and unique[item["path"]] != item:
            raise ValueError("Conflicting private deployment identity")
        unique[item["path"]] = item
        check(item)
    return dict(admission=admission_pin, readback=readback_pin, scheduler=scheduler,
                accounting=accounting, checked_records=list(unique.values()))


def verify_sources(root):
    authority = verify_private_control(root)
    from benchmark_tools.run_helper_cpm_phylogeny import verify_candidates
    from benchmark_tools.qfo_private_phylogeny_environment import verify_baseline

    candidates = verify_candidates(root)
    baseline = verify_baseline(root)
    return dict(private_control=authority, candidates=candidates, baseline=baseline,
        checked_records=[*authority["checked_records"], *candidates["checked_records"], *baseline["checked_records"]])


def native_command(verified, output):
    from benchmark_tools.run_ob_candidate_neighborhood import variant_cell
    from benchmark_tools.run_qfo_factorial_cell import native_command as frozen_command

    baseline = verified["baseline"]
    planned = variant_cell(baseline["original"], verified["candidates"]["arm"], output)
    argv, equivalence = frozen_command(planned, Path(baseline["launcher"]), Path(baseline["prepared"]))
    if (argv.count("--cpu") != 1 or argv[argv.index("--cpu") + 1] != "32"
            or argv.count("--species-tree-mode") != 1 or argv[argv.index("--species-tree-mode") + 1] != "infer"
            or "--species-tree" in argv or "--official-benchmark" in argv
            or argv.count("--checkpoint-source") != 1
            or argv[argv.index("--checkpoint-source") + 1] != planned["checkpoint_source"]):
        raise ValueError("Unexpected recovered inferred-phylogeny command")
    argv[0] = baseline["environment"]["tool_entrypoints"]["orthohmm_python"]["absolute_path"]
    return {**planned, "argv": argv}, planned, argv, equivalence


def run(root, protocol_sha):
    if (not os.environ.get("SLURM_JOB_ID") or os.environ.get("SLURM_CPUS_PER_TASK") != "32"
            or os.environ.get("SLURM_MEM_PER_NODE") != "196608" or os.environ.get("SLURM_JOB_NODELIST") != "bizon"
            or os.environ.get("SLURM_ARRAY_TASK_ID") or sys.executable != PYTHON):
        raise ValueError("Require standalone 32-CPU/192-GiB bizon task with historical verification controller")
    output = root / OUTPUT
    if output.exists() or output.is_symlink():
        raise FileExistsError(output)
    protocol = record(root / PROTOCOL)
    if protocol["sha256"] != protocol_sha:
        raise ValueError("Unreviewed private recovered-phylogeny protocol")
    verified = verify_sources(root)
    baseline, candidates = verified["baseline"], verified["candidates"]
    cell, planned, argv, equivalence = native_command(verified, output)
    from benchmark_tools.run_simulation_methods import read_frozen, execution_environment, execute
    from benchmark_tools.qfo_private_phylogeny_environment import inspect_launcher

    env, resolved = execution_environment(baseline["environment"])
    env.update(baseline["manifest"]["environment_overrides"], PYTHONPATH=baseline["launcher"],
        PYTHONDONTWRITEBYTECODE="1", PYTHONNOUSERSITE="1", PYTHONPYCACHEPREFIX=str(output / "bytecode_cache"),
        NUMBA_CACHE_DIR=str(output / "numba_cache"), NUMBA_CACHE_LOCATOR_CLASSES="UserProvidedCacheLocator")
    for key in ("PYTHONHOME", "LD_PRELOAD", "LD_LIBRARY_PATH", "LD_AUDIT"):
        env.pop(key, None)
    helpers = [record(path) for path in sorted(Path(__file__).parent.glob("*.py"))]
    provenance = dict(source=record(__file__), protocol=protocol, helpers=helpers, verified=verified,
        cell=cell, planned_cell=planned, executed_argv=argv, launcher_source_equivalence=equivalence,
        resolved_tools=resolved, cwd=baseline["launcher"], job_id=os.environ["SLURM_JOB_ID"],
        seed_handoff="explicit_helper_runtime_seed_amendment",
        native_handoff="explicit_admitted_private_phylogeny_deployment",
        scope="Unscored fixed recovered high-CPM arm; inferred phylogeny with validated raw-tree checkpoint reuse; incremental shared-host execution")
    output.mkdir(parents=True, exist_ok=False)
    (output / "numba_cache").mkdir()
    (output / "preflight.json").write_text(json.dumps(provenance, indent=2, sort_keys=True) + "\n")
    postflight = dict(status="running", accuracy_evaluated=False, native_outputs_validated=False, publication_ready=False)
    cwd = Path.cwd()
    try:
        lookup = inspect_launcher(baseline, env, output)
        fresh = output / "fresh_candidate_admission.json"
        admission = candidates["admission"]
        command = [sys.executable, "-B", str(Path(candidates["admission_executor"]) / "benchmark_tools/admit_helper_cpm_candidates.py"),
            "--root", str(root), "--preparation-sha256", admission["preparation"]["sha256"],
            "--protocol-sha256", admission["protocol"]["sha256"], "--output", str(fresh)]
        postflight["admission_command"] = command
        with (output / "admission.log").open("x") as log:
            subprocess.run(command, check=True, stdout=log, stderr=subprocess.STDOUT)
        fresh_pin = record(fresh)
        if read_frozen(fresh, fresh_pin["sha256"]) != admission:
            raise ValueError("Fresh recovered candidate admission disagrees")
        postflight["fresh_candidate_admission"] = fresh_pin
        method = dict(argv=argv, output=str(output / "output"), metrics=str(output / "metrics.json"))
        inputs = dict(status="ready", inputs=[{**r, "absolute_path": r["path"]} for r in baseline["manifest"]["input_fastas"]])
        os.chdir(baseline["launcher"])
        execution = execute(dict(label=cell["label"], methods={cell["label"]: method}), [cell["label"]],
            env, output / "execution", inputs, provenance)
        os.chdir(cwd)
        if execution["failed_methods"]:
            raise RuntimeError("Private recovered CPM phylogeny failed; preserve without retry")
        if verify_sources(root) != verified:
            raise ValueError("Private recovered phylogeny evidence/runtime changed during execution")
        for item in [provenance["source"], protocol, fresh_pin, *helpers,
                     *verified["checked_records"], *lookup["checked_records"]]:
            check(item)
        for pair in equivalence:
            check(pair["prepared"])
            check(pair["executed"])
        postflight.update(status="complete_pending_native_validation", cell=cell)
    except BaseException as error:
        postflight.update(status="failed", error_type=type(error).__name__, error=str(error))
        raise
    finally:
        os.chdir(cwd)
        (output / "postflight.json").write_text(json.dumps(postflight, indent=2, sort_keys=True) + "\n")
    return postflight


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", required=True, type=Path)
    parser.add_argument("--protocol-sha256", required=True)
    args = parser.parse_args()
    run(args.root.resolve(), args.protocol_sha256)
